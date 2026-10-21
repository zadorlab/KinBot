"""Production boundaries: portable checkpoints, same-level gaps and resource limits."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from kinbot.molpro import _vrc_multireference_method, stage_vrc_state_model
from kinbot.vrc_correction import resolved_state_model, spin_splitting_settings, assemble_spin_splitting_reference
from kinbot.vrc_tst_scan import vrc_dispatch_resources

ROOT = Path(__file__).parents[1]
PACKAGE = ROOT / 'examples/anl/tools'


@pytest.fixture(scope='module')
def installed(tmp_path_factory):
    """Copy real, tracked seeds using only the single per-PES input."""
    import shutil
    from kinbot.execution import file_sha256
    directory = tmp_path_factory.mktemp('portable_orbitals')
    (directory/'orbitals').mkdir()
    count = 0
    for case in ('ch3s', 'ch3ss'):
        example = ROOT/f'examples/anl/{case}_anl1_f12'
        data = json.loads((example/f'{case}.json').read_text())
        for channel in data['vrc_tst_surface_model']['channels']:
            for stage in channel['bottlenecks']:
                for model in stage['parameters']['vrc_tst_state_model']['basis_models'].values():
                    selection = model['orbital_selection']
                    source = example/selection['path']
                    assert source.is_file()
                    assert file_sha256(source) == selection['sha256']
                    shutil.copyfile(source, directory/'orbitals'/source.name)
                    count += 1
    assert count == 6
    return directory



def inline_selection(case, installed):
    molecule = case.split('_')[0]
    data = json.loads((ROOT/f'examples/anl/{molecule}_anl1_f12/{molecule}.json').read_text())
    index = 1 if case == 'ch3ss_ss' else 0
    payload = data['vrc_tst_surface_model']['channels'][index]['bottlenecks'][0]['parameters']['vrc_tst_state_model']
    for model in payload['basis_models'].values():
        selection = model['orbital_selection']
        selection['path'] = str(installed/'orbitals'/Path(selection['path']).name)
    return payload


@pytest.mark.parametrize('case,nel,norb,closed,occupied,roots', [
    ('ch3s_cs', 5, 4, 10, 14, 3), ('ch3ss_cs', 7, 5, 17, 22, 1),
    ('ch3ss_ss', 7, 5, 17, 22, 6)])
@pytest.mark.parametrize('basis', ['avdz-f12', 'avtz-f12'])
def test_actual_six_reviewed_seeds_reach_production_renderer(installed, tmp_path, case, nel, norb,
                                                           closed, occupied, roots, basis):
    payload = inline_selection(case, installed)
    parameters = {'vrc_tst_active_space': payload['active_space'],
                  'vrc_tst_state_model': payload}
    settings = {'target_multiplicity': 2, 'reference_multiplicity': 4}
    symbols = payload['molecule']['symbols']
    model = resolved_state_model(parameters, symbols, 0, settings, basis=basis)
    staged = stage_vrc_state_model(model, tmp_path, 'test')
    staged.pop('energy_selection', None)  # Explicit-root override path.
    for mult in (2, 4):
        staged['selected_roots'][mult] = roots
        deck, _ = _vrc_multireference_method(f'caspt2-f12({nel},{norb})',
            25 if case=='ch3s_cs' else 41, mult-1, state_model=staged)
        assert 'start,2140.2,type=natural' in deck
        assert f'occ,{occupied}\nclosed,{closed}\n' in deck
        assert 'diab,' not in deck and 'rotate,' not in deck
        assert deck.count('canonical,2141.2') == deck.count('natorb,2140.2') == 1
        assert 'shift=0.1,ipea=0.25' in deck
        assert f'state,1,{roots};refstate,{roots},' in deck
        assert ';orbital,2141.2,type=canonical;' in deck
        assert 'maxit,999;' in deck and 'maxit=2000,maxiti=2000,maxvi=60' in deck
        assert 'thrvar=1.d-10' in deck
        assert 'thrden=1.d-7' in deck
        assert 'accuracy,energy=1.d-8,gradient=1.d-5' in deck
    with pytest.raises(ValueError, match='atom order'):
        resolved_state_model(parameters, list(reversed(symbols)), 0, settings, basis=basis)


def test_ground_reference_evaluates_the_whole_manifold_with_unchanged_h0(installed, tmp_path):
    payload = inline_selection('ch3s_cs', installed)
    par = {'vrc_tst_active_space': payload['active_space'], 'vrc_tst_state_model': payload}
    model = resolved_state_model(par, payload['molecule']['symbols'], 0,
        {'target_multiplicity': 2, 'reference_multiplicity': 4}, basis='avtz-f12')
    model = stage_vrc_state_model(model, tmp_path, 'minimum')
    deck, selected = _vrc_multireference_method('caspt2-f12(5,4)', 25, 3, state_model=model)
    assert selected == 'KBVRCSEL'
    assert deck.count('{rs2-f12,shift=0.1,ipea=0.25,') == 3
    assert deck.count('refstate,3,1,2,3;weight,') == 3
    for root in (1, 2, 3):
        assert f'state,1,{root};refstate,3,1,2,3;' in deck
        assert f'KBGAPS({root})=energy(1)' in deck
    assert 'KBGAPS=0\n' in deck and 'KBVRCSEL=min(KBGAPS)' in deck
    target, _ = _vrc_multireference_method('caspt2-f12(5,4)', 25, 1, state_model=model)
    assert target.count('{rs2-f12,') == 1 and 'min(KBGAPS)' not in target
    from kinbot.ase_modules.calculators.molpro import _mr_gradient_variable
    from kinbot.molpro import vrc_multireference_energy_contract
    with pytest.raises(ValueError, match='energies only'):
        _mr_gradient_variable(vrc_multireference_energy_contract('caspt2-f12(5,4)'), model, 4)
    assert _mr_gradient_variable(vrc_multireference_energy_contract('caspt2-f12(5,4)'), model, 2) == 'ENERGY(1)'
    with pytest.raises(ValueError, match='requires CASPT2'):
        _vrc_multireference_method('mrci-f12+q(5,4)', 25, 3, state_model=model)


def test_cbs_gap_commutes_with_subtraction_and_keeps_tz_difference():
    settings = spin_splitting_settings({})
    energies = {'reference:0': [-49.8, -48.8],
        'reference:1': [-50., -49.],
        'gap:caspt2_f12:target:0': [-30.2, -29.1],
        'gap:caspt2_f12:reference:0': [-30., -29.],
        'gap:caspt2_f12:target:1': [-31.23, -30.12],
        'gap:caspt2_f12:reference:1': [-31., -30.]}
    curve, proof = assemble_spin_splitting_reference([3., 30.], energies, settings)
    dz, tz = np.array([-.2, -.1]), np.array([-.23, -.12])
    gap = tz + (tz-dz)*2**3.7/(3**3.7-2**3.7)
    reference = np.array([-50., -49.]) - .2*2**3.7/(3**3.7-2**3.7)
    assert curve == pytest.approx(reference + gap)
    assert proof['gap_cbs']['caspt2_f12']['cbs_minus_upper_gap_hartree'] == pytest.approx(gap-tz)
    assert 'no uncertainty' in proof['method_sensitivity']['interpretation']


def test_memory_budget_fits_documented_node_capacity():
    from kinbot.anl.dispatch import _molpro_stack_mw
    for name in ('ch3s','ch3ss'):
        par = json.loads((ROOT/f'examples/anl/{name}_anl1_f12/{name}.json').read_text())
        resources = vrc_dispatch_resources({**par, 'vrc_tst_walltime': '24:00:00'})
        assert resources['memory_mb'] == 104000 < 126000
        stack = _molpro_stack_mw({**resources, 'cores': 8})
        assert stack >= par['vrc_tst_min_stack_mw']
        assert (stack+200)*8*8 <= .85*104000
        assert par['queue_job_limit'] == par['rotdpy_max_jobs'] == par['vrc_tst_max_nodes'] == 10


@pytest.mark.parametrize('name', ['ch3s', 'ch3ss'])
def test_production_model_routing_and_pes_parameter_round_trip(name, tmp_path, installed):
    from kinbot.parameters import Parameters
    from kinbot.vrc_multistate import load_surface_model, state_parameters
    from kinbot.vrc_correction import correction_components
    path = ROOT/f'examples/anl/{name}_anl1_f12/{name}.json'
    par = Parameters(str(path)).par
    saved = tmp_path/'worker.json'
    saved.write_text(json.dumps(par))
    worker = Parameters(str(saved)).par
    assert worker['epsilon'] == par['epsilon']
    assert worker['epsilon_unit'] == 'cm-1'
    assert worker['simultaneous_kinbot'] == 1
    assert worker['barrier_threshold_L2'] == 400
    assert worker['barrier_threshold'] == 410
    assert worker['TemperatureList'] == list(range(300, 5001, 100))
    model, _ = load_surface_model(worker)
    assert len(model['channels']) == (1 if name == 'ch3s' else 2)
    for channel in model['channels']:
        assert channel['combination'] == 'single' and 'equilibration' not in channel
        stage, = channel['bottlenecks']
        # The checkout carries each seed directly; this test also exercises relocation.
        for basis_model in stage['parameters']['vrc_tst_state_model']['basis_models'].values():
            selection = basis_model['orbital_selection']
            selection['path'] = str(installed/'orbitals'/Path(selection['path']).name)
        for state in stage['states']:
            scalar = state_parameters(worker, stage, state, [2, 3])
            symbols = ['C', 'S'] + (['S'] if name == 'ch3ss' else []) + ['H']*3
            settings = spin_splitting_settings(scalar, parent_multiplicity=2,
                                                fragment_multiplicities=[2, 3])
            components = correction_components(scalar, symbols=symbols, charge=0,
                                               settings=settings)
            assert len(components) == 7  # frozen sample, two CC references, D/Q in each basis
            for component in components:
                if component['method'].startswith('caspt2'):
                    deck, _ = _vrc_multireference_method(component['method'],
                        25 if name == 'ch3s' else 41, component['multiplicity']-1,
                        state_model=component['state_model'],
                        shift=component.get('shift'), ipea=component.get('ipea'))
                    assert 'shift=0.1,ipea=0.25' in deck
                    assert 'refstate,' in deck and ';weight,' in deck
                    assert 'thrvar=1.d-10' in deck and 'maxit,999;' in deck


@pytest.mark.parametrize('case', ['ch3s', 'ch3ss'])
def test_each_example_runs_from_one_json_and_its_own_seed_data(case, tmp_path):
    """No other example, inventory, review JSON or installer supplies settings."""
    import shutil
    from kinbot.parameters import Parameters
    from kinbot.vrc_multistate import load_surface_model
    from kinbot.execution import file_sha256
    source = ROOT/f'examples/anl/{case}_anl1_f12'
    assert sorted(p.name for p in source.iterdir() if p.is_file()) == [case+'.json']
    destination = tmp_path/'standalone'
    destination.mkdir()
    shutil.copyfile(source/f'{case}.json', destination/f'{case}.json')
    shutil.copytree(source/'orbital_reviews/orbitals', destination/'orbital_reviews/orbitals')
    par = Parameters(str(destination/f'{case}.json')).par
    surface, _ = load_surface_model(par)
    for channel in surface['channels']:
        for stage in channel['bottlenecks']:
            envelope = stage['parameters']['vrc_tst_state_model']
            for basis, model in envelope['basis_models'].items():
                selection = model['orbital_selection']
                seed = Path(selection['path'])
                assert seed.is_relative_to(destination)
                assert file_sha256(seed) == selection['sha256']
                resolved = resolved_state_model(stage['parameters'], envelope['molecule']['symbols'], 0,
                    {'target_multiplicity': 2, 'reference_multiplicity': 4}, basis=basis)
                space = envelope['active_space']
                method = f"caspt2-f12({space['active_electrons']},{space['active_orbitals']})"
                deck, _ = _vrc_multireference_method(method, 25 if case == 'ch3s' else 41, 1,
                    state_model=stage_vrc_state_model(resolved, tmp_path/basis/('_'.join(map(str, channel['bond']))), 'standalone'))
                assert 'start,2140.2,type=natural' in deck
                assert 'shift=0.1,ipea=0.25' in deck


@pytest.mark.parametrize('case', ['ch3s', 'ch3ss'])
def test_retained_native_evidence_hashes_and_seed_sources(case):
    import tarfile
    import hashlib
    from kinbot.execution import file_sha256
    review = ROOT/f'examples/anl/{case}_anl1_f12/orbital_reviews'
    inventory = json.loads((review/'native_inventory.json').read_text())
    archive = review/inventory['archive']
    assert file_sha256(archive) == inventory['archive_sha256']
    expected = {entry['member']: entry for entry in inventory['files']}
    observed = set()
    with tarfile.open(archive) as tar:
        for member in tar:
            data = tar.extractfile(member).read()
            assert hashlib.sha256(data).hexdigest() == expected[member.name]['sha256']
            assert len(data) == expected[member.name]['bytes']
            observed.add(member.name)
    assert observed == set(expected)
    seed_inventory = json.loads((review/'orbital_inventory.json').read_text())
    native_hashes = {e['sha256'] for e in expected.values() if e['member'].endswith('.out')}
    for seed in seed_inventory['seeds']:
        assert file_sha256(review/seed['file']) == seed['sha256']
        assert seed['source_output_sha256'] in native_hashes
