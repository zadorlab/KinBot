from copy import deepcopy
import json

import pytest
from ase.units import Hartree, invcm

from kinbot.anl import dispatch, live_control, vpt2_refinement
from kinbot.anl.results import parse_gaussian_harmonic
from kinbot.execution import atomic_json, file_sha256
from tests.test_vpt2_common_recovery import graph, fail
from tests.test_anl_dispatch import _mock_execution


NATIVE = ('Full point group CS NOp 2\n'
          'Framework definition based on atomic masses: C2V\n'
          'The system is set in Eckart orientation for the anharmonic treatment.\n'
          'ERROR: Inconsistency found in framework group definition:\n'
          'New: C2V - Old: CS\nError termination via Lnk1e in l717.exe\n')
TRIED = ['NONE', 'PG=CS', 'PG=C2V', 'PG=C2V,LOOSE', 'PG=C1', 'NOSYMM', 'PG=CS,LOOSE,CENTEROFMASS']


def exhausted(tmp_path, monkeypatch):
    root = graph(tmp_path, symmetry='Symmetry=(PG=CS,Loose,CenterOfMass)')
    _, spec, state = dispatch._load(root)
    spec['intent'] = {'stationary_point': 'minimum'}
    task = next(t for t in spec['tasks'] if t['id'] == 'custom_anharmonic')
    task['input_template'] = task['input_template'].replace('SCF=XQC', 'SCF=(XQC,Conver=8) Integral=UltraFine')
    task['recovery'] = dict(attempted_symmetry_settings=TRIED)
    live_control.commit(root, spec, state, {'custom_anharmonic'}, {'custom_anharmonic': 'test prior ladder'})
    dispatch.advance(root)
    fail(root, NATIVE)
    _mock_execution(root, 'done')
    state = dispatch._load(root)[2]
    state['tasks']['live'].update(status='submitted', job_id='43')
    atomic_json(root/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda job: job == '43')
    return root


def harmonic_output(natoms=5):
    frequencies = [1000.]*(3*natoms-6)
    zpe = sum(frequencies)*.5*invcm/Hartree
    return (' #p B2PLYP/cc-pVTZ Freq EmpiricalDispersion=GD3BJ Integral=UltraFine\n\n'
            ' Frequencies -- '+' '.join(map(str, frequencies))+'\n'
            f' Zero-point correction= {zpe:.9f}\nNormal termination of Gaussian\n')


def complete_harmonic(root, ident):
    _mock_execution(root, ident)
    d = root/'tasks'/ident
    (d/'vpt2.log').write_text(harmonic_output())
    p = d/'execution.json'; e = json.loads(p.read_text())
    e['artifacts']['vpt2.log'] = file_sha256(d/'vpt2.log')
    e['details']['parsed_result'] = parse_gaussian_harmonic(harmonic_output(), method='B2PLYP', basis='cc-pVTZ', dispersion='GD3BJ', natoms=5)
    atomic_json(p, e)


def test_refines_once_requires_hessian_and_restarts_ladder_without_touching_l3(tmp_path, monkeypatch):
    root = exhausted(tmp_path, monkeypatch)
    siblings = {p:file_sha256(p) for n in ('l2_geometry', 'live', 'done')
                for p in (root/'tasks'/n).iterdir() if p.is_file()}
    state = dispatch.advance(root)
    geom, harmonic = 'custom_anharmonic_refine_geometry', 'custom_anharmonic_refine_harmonic'
    assert state['tasks'][geom]['status'] == 'staged'
    assert state['tasks'].get('custom_anharmonic', {}).get('status') in (None, 'blocked')
    _, spec, _ = dispatch._load(root)
    tasks = {t['id']:t for t in spec['tasks']}
    assert tasks[geom]['optimizer']['fmax'] == 1.e-5
    assert 'Conver=10' in tasks[geom]['profile']['calculator_kwargs']['scf']
    assert tasks['live']['geometry_from'] == 'initial'
    assert 'recovery' not in tasks['custom_anharmonic']
    assert 'Symm' not in tasks['custom_anharmonic']['input_template']
    _mock_execution(root, geom, geometry='final.xyz')
    state = dispatch.advance(root)
    assert state['tasks'][harmonic]['status'] == 'staged'
    assert state['tasks'].get('custom_anharmonic', {}).get('status') in (None, 'blocked')
    complete_harmonic(root, harmonic)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'staged'
    from kinbot.anl.workflow import task_component
    assert task_component(root, harmonic, key='harmonic_zpe', state_id='test').value_hartree > 0
    # Byte-identical optimized coordinates still constitute a fresh generation.
    fail(root, NATIVE)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'staged'
    assert 'Symmetry=None' in (root/'tasks/custom_anharmonic/vpt2.com').read_text()
    _, spec, state = dispatch._load(root)
    task = next(t for t in spec['tasks'] if t['id'] == 'custom_anharmonic')
    task['recovery']['attempted_symmetry_settings'] = TRIED
    live_control.commit(root, spec, state, {'custom_anharmonic'}, {'custom_anharmonic': 'remaining fixture attempts'})
    dispatch.advance(root); fail(root, NATIVE)
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['status'] == 'unavailable'
    assert state['tasks']['custom_anharmonic']['vpt2_refinement_exhausted']
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'unavailable'
    assert state['tasks']['live']['job_id'] == '43'
    assert siblings == {p:file_sha256(p) for p in siblings}
    assert len([t for t in spec['tasks'] if t['id'].endswith('_refine_geometry')]) == 1


@pytest.mark.parametrize('damage', ['saddle', 'active', 'custom'])
def test_refinement_applicability_fences(tmp_path, monkeypatch, damage):
    root = exhausted(tmp_path, monkeypatch)
    _, spec, state = dispatch._load(root)
    if damage == 'saddle': spec['intent']['stationary_point'] = 'transition_state'
    if damage == 'active': monkeypatch.setattr(dispatch, '_job_active', lambda job: True)
    if damage == 'custom':
        task = next(t for t in spec['tasks'] if t['id'] == 'custom_anharmonic')
        task['input_template'] = task['input_template'].replace('Integral=UltraFine', 'Integral=UltraFine SCRF=Water')
        with pytest.raises(ValueError, match='Custom VPT2'):
            vpt2_refinement._tasks(spec, task)
        return
    _, after = vpt2_refinement.recover(root, spec, state)
    assert not any(n.endswith('_refine_geometry') for n in after['tasks'])


@pytest.mark.parametrize('change', [lambda t:'', lambda t:t.replace('1000.0','-1000.0',1),
    lambda t:t.replace('Freq ', 'Freq=Anharmonic '), lambda t:t.replace('GD3BJ','GD3'),
    lambda t:t.replace('Normal termination','Error termination')])
def test_refined_hessian_must_be_a_valid_minimum(change):
    with pytest.raises(ValueError):
        parse_gaussian_harmonic(change(harmonic_output()), method='B2PLYP', basis='cc-pVTZ', dispersion='GD3BJ', natoms=5)


def test_final_unavailability_releases_cbh_owner_despite_general_failure_hold(tmp_path, monkeypatch):
    from kinbot.anl import automatic, refinement
    from kinbot import live_control as control
    spec = {'tasks':[], 'intent':{}}
    state = {'tasks':{'gaussian_vpt2':{'status':'unavailable', 'vpt2_refinement_exhausted':True}}}
    monkeypatch.setattr(automatic, '_load', lambda p:(p, spec, state))
    monkeypatch.setattr(automatic, 'advance', lambda *a,**k:state)
    monkeypatch.setattr(automatic, '_retry_l3_execution_failures', lambda *a:state)
    monkeypatch.setattr(automatic, '_retry_vpt2_framework_failures', lambda *a:state)
    monkeypatch.setattr(live_control, 'process', lambda p:None)
    monkeypatch.setattr(refinement, 'window', lambda *a,**k:True)
    monkeypatch.setattr(control, 'holding', lambda *a,**k:True)
    monkeypatch.setattr(automatic.time, 'sleep', lambda *a:pytest.fail('Terminal VPT2 must not leave CBH waiting'))
    assert automatic._drive_graphs_locked({'case':tmp_path}, 1, .01) == {'case':state}


def test_refined_l2_harmonic_fallback_keeps_existing_l3_sp_geometry(tmp_path, monkeypatch):
    from dataclasses import replace
    from kinbot.anl import automatic
    molecule = dict(symbols=['C','H','H','H','H'], charge=0, multiplicity=1,
                    source=dict(geometry_sha256='a'*64, record_sha256='b'*64))
    source = dict(molecule=molecule, smiles='C', energy_hartree=-40., harmonic_zpe_hartree=.03,
                  l2_theory=dict(composite_method='B2PLYP-D3BJ', basis='cc-pVTZ', label='B2PLYP-D3BJ/cc-pVTZ'))
    receipt = dict(geometry_task='tight_geometry', harmonic_task='tight_hessian')
    spec = dict(name='anl-automatic-post-l2', molecule=molecule, tasks=[dict(id='gaussian_vpt2',
                vpt2_geometry_refinement=receipt)])
    spec['tasks'].extend(dict(id=n) for n in ('tight_geometry','tight_hessian','l3_sp_l2'))
    state = dict(tasks={n:dict(status=s) for n,s in dict(gaussian_vpt2='unavailable',
                tight_geometry='complete', tight_hessian='complete', l3_sp_l2='complete').items()})
    monkeypatch.setattr(automatic, '_load', lambda p:(tmp_path,spec,state))
    base = automatic._l2_source_component(source, key='electronic', quantity='electronic',
            value=-40., method='B2PLYP-D3BJ', basis='cc-pVTZ', state_id='test')
    def component(path, task, *, key, state_id):
        if task == 'l3_sp_l2': return replace(base, method='CCSD(T)', value_hartree=-40.1)
        assert task == 'tight_hessian'
        return replace(base, key=key, quantity='zpe', value_hartree=.029, geometry_sha256='c'*64)
    monkeypatch.setattr(automatic, 'task_component', component)
    monkeypatch.setattr(automatic, '_ase_energy_component', lambda *a,**k:replace(base, value_hartree=-40.001, geometry_sha256='c'*64))
    profiles = automatic.assemble_available_profiles(tmp_path,source,tmp_path/'profiles',state_id='test')
    records = [json.loads(__import__('pathlib').Path(p).read_text()) for p in profiles['thermochemistry'].values()]
    l2 = next(r for r in records if r['recipe'].startswith('L2:'))
    l3 = next(r for r in records if r['recipe'].startswith('L3:'))
    assert 'vpt2_correction' not in l2['components']
    assert 'vpt2_correction' not in l3['components']
    assert l2['zero_k_hartree'] == pytest.approx(-40.001+.029)
    assert l3['zero_k_hartree'] == pytest.approx(-40.1+.03)
    assert {c['geometry_sha256'] for c in l2['components'].values()} == {'c'*64}
    assert {c['geometry_sha256'] for c in l3['components'].values()} == {'a'*64}
