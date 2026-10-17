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
PACKAGE = ROOT / 'examples/anl/caspt2_production'


@pytest.fixture(scope='module')
def installed(tmp_path_factory):
    import shutil
    directory = tmp_path_factory.mktemp('portable_orbitals')
    for name in ('orbital_inventory.json', 'reviewed_orbitals.tar.gz',
                 'ch3s_cs_state_model.json', 'ch3ss_cs_state_model.json', 'ch3ss_ss_state_model.json'):
        shutil.copyfile(PACKAGE/name, directory/name)
    spec = importlib.util.spec_from_file_location('install_seeds', PACKAGE/'install_orbitals.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.install(directory)['installed'] == 6
    assert mod.install(directory)['installed'] == 6  # immutable/idempotent
    return directory


@pytest.mark.parametrize('case,nel,norb,closed,occupied,roots', [
    ('ch3s_cs', 5, 4, 10, 14, 3), ('ch3ss_cs', 7, 5, 17, 22, 1),
    ('ch3ss_ss', 7, 5, 17, 22, 6)])
@pytest.mark.parametrize('basis', ['avdz-f12', 'avtz-f12'])
def test_actual_six_reviewed_seeds_reach_production_renderer(installed, tmp_path, case, nel, norb,
                                                           closed, occupied, roots, basis):
    payload = json.loads((installed/(case+'_state_model.json')).read_text())
    parameters = {'vrc_tst_active_space': payload['active_space'],
                  'vrc_tst_state_model_file': str(installed/(case+'_state_model.json'))}
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
    with pytest.raises(ValueError, match='atom order'):
        resolved_state_model(parameters, list(reversed(symbols)), 0, settings, basis=basis)


def test_ground_reference_evaluates_the_whole_manifold_with_unchanged_h0(installed, tmp_path):
    path = installed/'ch3s_cs_state_model.json'
    payload = json.loads(path.read_text())
    par = {'vrc_tst_active_space': payload['active_space'], 'vrc_tst_state_model_file': str(path)}
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
    energies = {'reference': [-50., -49.],
        'gap:caspt2_f12:target:0': [-30.2, -29.1],
        'gap:caspt2_f12:reference:0': [-30., -29.],
        'gap:caspt2_f12:target:1': [-31.23, -30.12],
        'gap:caspt2_f12:reference:1': [-31., -30.]}
    curve, proof = assemble_spin_splitting_reference([3., 30.], energies, settings)
    dz, tz = np.array([-.2, -.1]), np.array([-.23, -.12])
    gap = tz + (tz-dz)*2**3.7/(3**3.7-2**3.7)
    assert curve == pytest.approx(np.array([-50., -49.]) + gap)
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
    assert worker['simultaneous_kinbot'] == 1 and worker['barrier_threshold'] == 1000
    model, _ = load_surface_model(worker)
    assert len(model['channels']) == (1 if name == 'ch3s' else 2)
    for channel in model['channels']:
        assert channel['combination'] == 'single' and 'equilibration' not in channel
        stage, = channel['bottlenecks']
        # A fresh checkout carries the archive, not extracted checkpoints.
        model_name = Path(stage['parameters']['vrc_tst_state_model_file']).name
        assert (installed/model_name).is_file()
        stage['parameters']['vrc_tst_state_model_file'] = str(installed/model_name)
        for state in stage['states']:
            scalar = state_parameters(worker, stage, state, [2, 3])
            symbols = ['C', 'S'] + (['S'] if name == 'ch3ss' else []) + ['H']*3
            settings = spin_splitting_settings(scalar, parent_multiplicity=2,
                                                fragment_multiplicities=[2, 3])
            components = correction_components(scalar, symbols=symbols, charge=0,
                                               settings=settings)
            assert len(components) == 6  # frozen sample, CC, D/Q in each basis
            for component in components:
                if component['method'].startswith('caspt2'):
                    deck, _ = _vrc_multireference_method(component['method'],
                        25 if name == 'ch3s' else 41, component['multiplicity']-1,
                        state_model=component['state_model'],
                        shift=component.get('shift'), ipea=component.get('ipea'))
                    assert 'shift=0.1,ipea=0.25' in deck
                    assert 'refstate,' in deck and ';weight,' in deck
                    assert 'thrvar=1.d-10' in deck and 'maxit,999;' in deck
