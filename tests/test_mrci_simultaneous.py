"""Root-protection input/selection checks, not claims of native convergence."""
import tempfile
import json
from pathlib import Path
import re

from ase import Atoms
import pytest

from kinbot.molpro import _vrc_multireference_method
from kinbot import molpro_reconnaissance as recon
from kinbot.orbital_monitor import _native_progress
from test_molpro_confirmation_recovery import fixture
from test_molpro_confirmation_package import PACKAGE, module


def model():
    return {'state_average': [{'multiplicity': 2, 'roots': 3},
                              {'multiplicity': 4, 'roots': 3}],
            'selected_roots': {'2': 3, '4': 2},
            'correlation_options': {'reference_policy': 'state_average',
                                    'state_policy': 'simultaneous',
                                    'max_external_vectors': 20}}


@pytest.mark.parametrize('spin,slot', [(1, 3), (3, 2)])
def test_simultaneous_production_selects_requested_slot_and_keeps_relaxed_davidson(spin, slot):
    text, variable = _vrc_multireference_method('mrci-f12+q(5,4)', 25,
                                               spin=spin, state_model=model())
    assert 'noswap' in text and 'state,3,1,2,3;' in text
    assert 'refstate,3,1,2,3;option,nstati=3;' in text
    assert f'KBVRCSEL=energd({slot})' in text
    assert f'KBVRCVAR=energy({slot})' in text
    assert 'option,cluster=0;' in text and variable == 'KBVRCSEL'
    assert 'maxdav=20' in text


def test_simultaneous_policy_does_not_silently_change_caspt2_or_reference():
    with pytest.raises(ValueError, match='Simultaneous correlated states require MRCI'):
        _vrc_multireference_method('caspt2-f12(5,4)', 25, spin=1, state_model=model())
    invalid = model()
    invalid['correlation_options']['reference_policy'] = 'selected_root'
    with pytest.raises(ValueError, match='reference_policy'):
        _vrc_multireference_method('mrci-f12+q(5,4)', 25, spin=1, state_model=invalid)


def test_all_resolution_decks_group_roots_and_survive_native_transform(tmp_path, monkeypatch):
    from types import SimpleNamespace
    adapter, source, target, _, review = fixture(tmp_path, monkeypatch)
    record = adapter.prepare(source, target, PACKAGE/'multistate_resolution_plan.json', review)
    assert sum(s['jobs'] for s in record['systems']) == 7
    assert record['execution']['max_jobs'] == 3
    assert record['execution']['memory_mb'] == 104000
    monkeypatch.setattr(recon, 'qc_runtime_environment', lambda *a, **k: ({}, {'scratch_directory': tempfile.mkdtemp(dir=tmp_path)}))
    monkeypatch.setattr(recon, 'cleanup_qc_runtime', lambda *a: None)
    inspected = []
    def launch(command, **kwargs):
        text = (kwargs['cwd']/command[-1]).read_text()
        blocks = re.findall(r'\{mrci-f12[^{}]*\}', text)
        assert len(blocks) == 2
        expected = 6 if 'ch3ss_ss' in str(job['directory']) else 3
        assert len(re.findall(r'^reconc\d+=', text, re.M)) == 2*expected
        for block in blocks:
            assert f'state,{expected},'+','.join(map(str,range(1,expected+1)))+';' in block
            assert f'refstate,{expected},'+','.join(map(str,range(1,expected+1)))+';' in block
            assert f'nstati={expected}' in block and 'maxdav=20' in block
            assert 'noswap' in block and 'cluster=0' in block
            assert 'thrvar=1.d-10' in block and 'maxvi=60' in block
            cabs = 1 if 'simultaneous_control' in str(job['directory']) else 3
            assert f'cabs_singles={cabs}' in block
        for slot in range(1,expected+1):
            assert text.count(f'KBVRCSEL=energd({slot})') == 2
        inspected.append(command)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', launch)
    for system in record['systems']:
        m=json.loads(Path(system['manifest']).read_text())
        for job in m['jobs']:
            job['retry_solver_options']=m['solver_options']
            if job.get('orbital_predecessor'):
                (Path(job['directory'])/(job['name']+'.wfu')).write_bytes(b'synthetic checkpoint')
            recon._execute(job)
    assert len(inspected) == 7


def test_grouped_monitor_tracks_second_spin_not_second_root(tmp_path):
    job={'name':'point','state_model':model(), 'correlation_checks':[{'method':'mrci-f12+q(5,4)'}],
         'root_states':[{'multiplicity':m,'root':r} for m in (2,4) for r in (1,2,3)]}
    (tmp_path/'point.out').write_text('PROGRAM * MULTI\nPROGRAM * MRCI-F12\n'
        'RECONC1 = -20.3 AU\nRECONC2 = -20.2 AU\nRECONC3 = -20.1 AU\nPROGRAM * MRCI-F12\n')
    progress=_native_progress(tmp_path,job)
    assert progress['printed']==3 and progress['calls_started']==2
    assert progress['current']['multiplicity']==4
    assert progress['current']['root']=='1,2,3 (simultaneous)'


def test_coupled_cabs_review_does_not_add_a_separate_correction():
    review=module('review_export')
    text='''PROGRAM * MRCI-F12
Cluster corrected energies -20.5 (Davidson, relaxed reference)
SETTING KBVRCSEL = -20.5 AU
'''
    checks=review.mrci_energy_checks(text,'{mrci-f12,cabs_singles=3;}')
    assert checks[0]['rounding_residual_hartree']==0
    assert checks[0]['separate_cabs_correction'] is None
    with pytest.raises(ValueError,match='lacks'):
        review.mrci_energy_checks(text,'{mrci-f12,cabs_singles=1;}')
    ambiguous=text.replace('SETTING KBVRCSEL',
        'Cluster corrected energies -20.6 (Davidson, relaxed reference)\nSETTING KBVRCSEL')
    with pytest.raises(ValueError,match='Ambiguous'):
        review.mrci_energy_checks(ambiguous,'{mrci-f12,cabs_singles=3;state,1,1;}')


def test_simultaneous_numerical_gradient_uses_same_state_as_energy():
    from kinbot.ase_modules.calculators.molpro import render_input
    # A displaced geometry must differentiate the requested state, not element 1.
    text=render_input(Atoms('CH3S'), charge=0, mult=2, basis='avtz-f12',
                     geometry_name='point.xyz', method='mrci-f12+q(5,4)',state_model=model())
    assert 'KBVRCSEL=energd(3)' in text
    assert 'forces,numerical,variable=energd(3)' in text
