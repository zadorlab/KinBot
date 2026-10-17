"""Automatic CBH/ANL planning starts from the accepted KinBot network."""

import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule

from kinbot.anl.automatic import (DEFAULT_SIZE_POLICY, _model_record,
                                  _component_family_omissions,
                                  _drive_graphs, _ordered_methods,
                                  _reference_input,
                                  _run_reference_l2,
                                  _reference_thermochemistry_models,
                                  _vpt2_zpe_components,
                                  automatic_requested,
                                  build_network_request,
                                  finalize_automatic_run, highest_tier,
                                  validate_automatic_environment)
from kinbot.anl.validation import (automatic_composite_spec,
                                   molecule_from_smiles)
from kinbot.anl.dispatch import (_geometry_hash, SchedulerUnavailableError,
                                 validate_spec)


def _mock_scientific_runtime(monkeypatch, tmp_path, found):
    """Provide the licensed-program/cache contract for environment unit tests."""
    found.update({
        'g16': '/opt/gaussian/g16/g16',
        'molpro': '/opt/molpro/bin/molpro',
        'xcfour': '/opt/cfour/bin/xcfour',
        'dmrcc': '/opt/mrcc/bin/dmrcc',
        'formchk': '/opt/gaussian/g16/formchk',
        'cubegen': '/opt/gaussian/g16/cubegen',
    })
    genbas = tmp_path / 'GENBAS'
    genbas.write_text('test basis library\n')
    checkpoint = tmp_path / 'uma-s-1p2.pt'
    checkpoint.write_bytes(b'test checkpoint')
    monkeypatch.setenv('CFOUR_GENBAS', str(genbas))
    monkeypatch.setitem(sys.modules, 'fairchem', SimpleNamespace())
    monkeypatch.setattr(
        'kinbot.fairchem_utils.resolve_model_checkpoint',
        lambda model: checkpoint)


def _species(name, smiles, atoms, multiplicity=1, energy=-1., zpe=.01):
    return SimpleNamespace(
        name=name, chemid=name, wellorts=1, smiles=smiles,
        atom=list(atoms), geom=[[float(i), 0., 0.] for i in range(len(atoms))],
        charge=0, mult=multiplicity, energy=energy, zpe=zpe,
        freq=[500., 1000.], source_job=name + '_high', source_row_id=7)


def test_explicit_spin_orbit_correction_survives_l2_fallback_expression(tmp_path, monkeypatch):
    from kinbot.anl.automatic import assemble_available_profiles
    from kinbot.anl.cbh import zero_k_from_composite_record
    molecule = {'symbols': ['S'], 'positions': [[0., 0., 0.]],
                'charge': 0, 'multiplicity': 3, 'smiles': '[S]',
                'source': {'geometry_sha256': 'a' * 64, 'record_sha256': 'b' * 64}}
    source = {'molecule': molecule, 'smiles': '[S]', 'energy_hartree': -398.,
              'harmonic_zpe_hartree': 0.,
              'l2_theory': {'composite_method': 'B2PLYP-D3BJ', 'basis': 'cc-pVTZ',
                            'label': 'B2PLYP-D3BJ/cc-pVTZ'},
              'spin_orbit_correction': {'value_hartree': -.001,
                                         'source': 'test complete-term centroid'}}
    spec = {'name': 'anl-automatic-post-l2', 'molecule': molecule, 'tasks': []}
    monkeypatch.setattr('kinbot.anl.automatic._load', lambda path: (tmp_path, spec, {'tasks': {}}))
    profiles = assemble_available_profiles(tmp_path, source, tmp_path / 'profiles', state_id='sulfur')
    for role in ('thermochemistry', 'barrier'):
        for path in profiles[role].values():
            record = json.loads(Path(path).read_text())
            assert record['zero_k_hartree'] == pytest.approx(-398.001)
            assert record['components']['spin_orbit']['source'] == 'test complete-term centroid'
            assert zero_k_from_composite_record('[S]', path).hartree == pytest.approx(-398.001)


def test_network_request_generates_recursive_cbh_states_without_user_species():
    acetaldehyde = _species(
        'acetaldehyde', 'CC=O', ['C', 'C', 'O', 'H', 'H', 'H', 'H'],
        energy=-153.)
    methyl = _species('methyl', '[CH3]', ['C', 'H', 'H', 'H'], 2,
                      energy=-39.)
    formyl = _species('formyl', '[CH]=O', ['C', 'O', 'H'], 2,
                      energy=-113.)
    reaction = SimpleNamespace(
        instance_name='acetaldehyde_hom_sci_0_1', family_name='hom_sci',
        products=[methyl, formyl])
    acetaldehyde.reac_obj = [reaction]
    acetaldehyde.reac_ts_done = [-1]
    acetaldehyde.reac_type = ['hom_sci']

    request = build_network_request(acetaldehyde, {
        'composite_method': 'ANL0-F12',
        'anl_size_policy': DEFAULT_SIZE_POLICY,
        'queue_name': 'day-long-cpu', 'ppn': 8, 'queue_job_limit': 8,
        'atct_version': '1.220',
    })
    assert request['resources']['global_max_jobs'] == 8
    assert len(request['network']['stable_species']) == 3
    assert request['network']['reactions'][0]['barrierless'] is True
    assert request['network']['transition_states'] == {}
    assert set(request['cbh']['targets']) == {
        'acetaldehyde', 'methyl', 'formyl'}
    acetaldehyde_target = request['cbh']['targets']['acetaldehyde']
    assert request['requested_cbh_rung'] == 'auto'
    assert acetaldehyde_target['rung_policy'] == 'auto'
    assert acetaldehyde_target['selected_structural_rung'] == 1
    assert acetaldehyde_target['available_structural_rungs'] == [1, 0]
    assert acetaldehyde_target['reaction']['stoichiometry'] == {
        'C': -1, 'C=O': 1, 'CC': 1, 'CC=O': -1}
    calculated = {entry['smiles'] for entry in
                  request['cbh']['calculated_states']}
    assert {'CC=O', 'C', 'C=O', 'CC'} <= calculated
    assert all(entry['fallback_tiers'][-1] == 'L2'
               for entry in request['cbh']['calculated_states'])


def test_public_reaction_identity_controls_barrierless_classification():
    acetaldehyde = _species(
        'acetaldehyde', 'CC=O', ['C', 'C', 'O', 'H', 'H', 'H', 'H'],
        energy=-153.)
    methyl = _species('methyl', '[CH3]', ['C', 'H', 'H', 'H'], 2,
                      energy=-39.)
    formyl = _species('formyl', '[CH]=O', ['C', 'O', 'H'], 2,
                      energy=-113.)
    reaction = SimpleNamespace(
        instance_name='acetaldehyde_hom_sci_1_2', family_name='hs',
        instance=[0, 1], products=[methyl, formyl])
    acetaldehyde.reac_obj = [reaction]
    acetaldehyde.reac_inst = [[0, 1]]
    acetaldehyde.reac_ts_done = [-1]
    acetaldehyde.reac_type = ['hom_sci']

    request = build_network_request(acetaldehyde, {
        'composite_method': 'ANL0-F12',
        'selected_reactions': {'acetaldehyde': {'hom_sci': [[0, 1]]}},
        'anl_size_policy': DEFAULT_SIZE_POLICY,
        'queue_name': 'day-long-cpu', 'ppn': 8, 'queue_job_limit': 8,
        'atct_version': '1.220',
    })

    record = request['network']['reactions'][0]
    assert record['family'] == 'hom_sci'
    assert record['barrierless'] is True
    assert record['transition_state'] is None


def test_every_generated_acetaldehyde_composite_graph_validates():
    root = _species(
        'acetaldehyde', 'CC=O', ['C', 'C', 'O', 'H', 'H', 'H', 'H'],
        energy=-153.)
    methyl = _species('methyl', '[CH3]', ['C', 'H', 'H', 'H'], 2,
                      energy=-39.)
    formyl = _species('formyl', '[CH]=O', ['C', 'O', 'H'], 2,
                      energy=-113.)
    methane = _species('methane', 'C', ['C', 'H', 'H', 'H', 'H'],
                       energy=-40.)
    carbon_monoxide = _species('carbon_monoxide', '[C-]#[O+]', ['C', 'O'],
                               energy=-113.)
    ketene = _species('ketene', 'C=C=O',
                      ['C', 'C', 'O', 'H', 'H'], energy=-152.)
    hydrogen = _species('hydrogen', '[H][H]', ['H', 'H'], energy=-1.)
    insertion = SimpleNamespace(
        instance_name='acetaldehyde_r12_insertion_R_7_2_1',
        family_name='r12_insertion_R', products=[methane, carbon_monoxide],
        ts=_species('insertion_ts', '', root.atom, energy=-152.8))
    elimination = SimpleNamespace(
        instance_name='acetaldehyde_h2_elim_4_7', family_name='h2_elim',
        products=[ketene, hydrogen],
        ts=_species('elimination_ts', '', root.atom, energy=-152.7))
    scission = SimpleNamespace(
        instance_name='acetaldehyde_hom_sci_1_2', family_name='hom_sci',
        products=[methyl, formyl])
    root.reac_obj = [insertion, elimination, scission]
    root.reac_ts_done = [-1, -1, -1]
    root.reac_type = ['r12_insertion_R', 'h2_elim', 'hom_sci']
    plan = build_network_request(root, {
        'composite_method': 'ANL0-F12',
        'anl_size_policy': DEFAULT_SIZE_POLICY,
        'queue_name': 'day-long-cpu', 'ppn': 8, 'queue_job_limit': 8,
        'atct_version': '1.220',
    })

    def imported_molecule(entry):
        record = entry.get('network_record')
        if record is not None:
            return record['molecule']
        value = molecule_from_smiles(
            entry['smiles'], charge=entry['charge'],
            multiplicity=entry['multiplicity'])
        value['source'] = {
            'geometry_sha256': _geometry_hash(Atoms(
                value['symbols'], positions=value['positions'])),
            'record_sha256': hashlib.sha256(
                entry['smiles'].encode()).hexdigest(),
            'energy_hartree': -1.,
            'l2_theory': dict(plan['l2_theory']),
        }
        return value

    generated = []
    for entry in plan['cbh']['calculated_states']:
        generated.append(automatic_composite_spec(
            imported_molecule(entry), highest_tier=entry['highest_tier'],
            max_nodes=4, partition='day-long-cpu',
            vpt2_method=plan['l2_theory']['method'],
            vpt2_basis=plan['l2_theory']['basis'],
            vpt2_dispersion=plan['l2_theory']['dispersion']))
    for entry in plan['network']['transition_states'].values():
        generated.append(automatic_composite_spec(
            entry['molecule'], highest_tier=entry['highest_tier'],
            max_nodes=4, partition='day-long-cpu',
            stationary_point='transition_state',
            vpt2_method=plan['l2_theory']['method'],
            vpt2_basis=plan['l2_theory']['basis'],
            vpt2_dispersion=plan['l2_theory']['dispersion']))
    assert generated
    for specification in generated:
        for task in specification['tasks']:
            task['resources'].update(
                cores=1, memory_mb=64000, partition='test')
        validate_spec(specification)


def test_tight_transition_state_uses_explicit_reactant_energy_identity():
    acetaldehyde = _species(
        'acetaldehyde', 'CC=O', ['C', 'C', 'O', 'H', 'H', 'H', 'H'],
        energy=-153.)
    methane = _species('methane', 'C', ['C', 'H', 'H', 'H', 'H'],
                       energy=-40.)
    formaldehyde = _species('formaldehyde', 'C=O', ['C', 'O', 'H', 'H'],
                          energy=-113.)
    transition_state = _species(
        'acetaldehyde_tight_ts', '',
        ['C', 'C', 'O', 'H', 'H', 'H', 'H'], energy=-152.9)
    reaction = SimpleNamespace(
        instance_name='acetaldehyde_r12_insertion_R_7_2_1',
        family_name='r12_insertion_R', products=[methane, formaldehyde],
        ts=transition_state)
    acetaldehyde.reac_obj = [reaction]
    acetaldehyde.reac_ts_done = [-1]
    acetaldehyde.reac_type = ['r12_insertion_R']

    request = build_network_request(acetaldehyde, {
        'composite_method': 'ANL0-F12', 'cbh_rung': 1,
        'anl_size_policy': DEFAULT_SIZE_POLICY,
        'queue_name': 'day-long-cpu', 'ppn': 8, 'queue_job_limit': 8,
        'atct_version': '1.220',
    })
    record = request['network']['transition_states'][
        reaction.instance_name]
    assert record['smiles'] == 'CC=O'
    assert record['routing_name'] == reaction.instance_name
    assert record['identity_kind'] == 'reactant_connectivity_anchor'


def test_selected_reaction_failure_blocks_automatic_network_planning():
    acetaldehyde = _species(
        'acetaldehyde', 'CC=O', ['C', 'C', 'O', 'H', 'H', 'H', 'H'],
        energy=-153.)
    acetaldehyde.short_name = 'w1'
    reaction = SimpleNamespace(
        instance_name='acetaldehyde_h2_elim_4_7',
        family_name='h2elim', instance=[3, 6], products=[])
    acetaldehyde.reac_obj = [reaction]
    acetaldehyde.reac_inst = [[3, 0, 1, 6]]
    acetaldehyde.reac_ts_done = [-999]
    acetaldehyde.reac_type = ['h2_elim']
    parameters = {
        'composite_method': 'ANL0-F12',
        'selected_reactions': {
            'initial': {'h2_elim': [[3, 0, 1, 6]]}},
        'anl_size_policy': DEFAULT_SIZE_POLICY,
        'queue_name': 'day-long-cpu', 'ppn': 8, 'queue_job_limit': 8,
        'atct_version': '1.220',
    }
    with pytest.raises(RuntimeError, match='explicitly selected.*h2_elim'):
        build_network_request(acetaldehyde, parameters)

    acetaldehyde.reac_ts_done = [-1]
    reaction.ts = _species(
        reaction.instance_name, '',
        ['C', 'C', 'O', 'H', 'H', 'H', 'H'], energy=-152.8)
    request = build_network_request(acetaldehyde, parameters)
    assert [item['family'] for item in request['network']['reactions']] == [
        'h2_elim']


def test_size_policy_downgrades_monotonically():
    tier, metrics = highest_tier('CC=O', 0, 1, DEFAULT_SIZE_POLICY)
    assert tier == 'ANL0-F12'
    assert metrics['heavy_atoms'] == 3
    restricted = {
        'full_anl0_f12': {
            'max_atoms': 1, 'max_heavy_atoms': 1, 'max_heavy_bonds': 1},
        'l3_geometry': {
            'max_atoms': 2, 'max_heavy_atoms': 2, 'max_heavy_bonds': 2},
        'l3_single_point': {
            'max_atoms': 3, 'max_heavy_atoms': 3, 'max_heavy_bonds': 3},
    }
    assert highest_tier('CC=O', 0, 1, restricted)[0] == 'L2'
    assert highest_tier(
        'CC=O', 0, 1, DEFAULT_SIZE_POLICY,
        requested_method='ANL0')[0] == 'ANL0'


def test_automatic_anl0_and_f12_graphs_have_distinct_reference_fallbacks():
    atoms = molecule('CH4')
    source = {
        'geometry_sha256': '1' * 64, 'record_sha256': '2' * 64,
        'energy_hartree': -40., 'harmonic_zpe_hartree': .04,
    }
    mol = {
        'symbols': atoms.get_chemical_symbols(),
        'positions': atoms.positions.tolist(), 'charge': 0, 'multiplicity': 1,
        'source': source,
    }
    f12 = automatic_composite_spec(mol, highest_tier='ANL0-F12')
    conventional = automatic_composite_spec(mol, highest_tier='ANL0')
    f12_tasks = {task['id']: task for task in f12['tasks']}
    assert 'Symmetry=None' in f12_tasks['gaussian_vpt2']['input_template']
    conventional_ids = {task['id'] for task in conventional['tasks']}
    assert {'f12_tz', 'f12_qz', 'anl0_qz', 'anl0_5z'} <= set(f12_tasks)
    assert {'f12_5z', 'anl1_5z', 'anl1_6z'}.isdisjoint(f12_tasks)
    assert {'anl0_qz', 'anl0_5z'} <= conventional_ids
    assert not {'f12_tz', 'f12_qz'} & conventional_ids
    assert f12_tasks['anl0_qz']['failure_skips'] == ['anl0_5z']
    assert f12_tasks['anl0_5z']['failure_skips'] == ['anl0_qz']
    assert 'basis=aug-cc-pVQZ' in f12_tasks['anl0_qz']['input_template']


def test_automatic_f12_augmentation_updates_native_tasks_and_intent():
    atoms = molecule('CH4')
    mol = {
        'symbols': atoms.get_chemical_symbols(),
        'positions': atoms.positions.tolist(), 'charge': 0, 'multiplicity': 1,
        'source': {
            'geometry_sha256': '1' * 64, 'record_sha256': '2' * 64,
            'energy_hartree': -40., 'harmonic_zpe_hartree': .04,
        },
    }
    spec = automatic_composite_spec(
        mol, highest_tier='ANL0-F12', f12_augmented=True)
    tasks = {task['id']: task for task in spec['tasks']}
    assert tasks['f12_tz']['result_parser']['basis'] == 'aug-cc-pVTZ-F12'
    assert tasks['f12_qz']['result_parser']['basis'] == 'aug-cc-pVQZ-F12'
    assert 'basis=aug-cc-pVTZ-F12' in tasks['f12_tz']['input_template']
    assert spec['intent']['f12_augmented'] is True
    assert tasks['l3_geometry']['optimizer']['fmax'] == 1e-4
    tight = automatic_composite_spec(mol, highest_tier='ANL1-F12', geometry_fmax=5e-5)
    assert all(t['optimizer']['fmax'] == 5e-5 for t in tight['tasks'] if t['kind'] == 'ase_optimize')
    assert all(t['optimizer']['minimum_recovery']['attempts'] == 2
               for t in tight['tasks'] if t['kind'] == 'ase_optimize')
    disabled = automatic_composite_spec(mol, highest_tier='ANL1-F12',
                                       minimum_recovery={'attempts': 0})
    assert all(t['optimizer']['minimum_recovery']['attempts'] == 0
               for t in disabled['tasks'] if t['kind'] == 'ase_optimize')
    saddle = automatic_composite_spec(mol, highest_tier='ANL1-F12',
                                     stationary_point='transition_state')
    assert all('minimum_recovery' not in t['optimizer']
               for t in saddle['tasks'] if t['kind'] == 'ase_optimize')


def test_automatic_environment_checks_final_solvers_before_qc(
        monkeypatch, tmp_path):
    found = {
        'mess': '/opt/mess/bin/mess', 'messpf': '/opt/mess/bin/messpf',
        'sbatch': '/usr/bin/sbatch', 'squeue': '/usr/bin/squeue'}
    _mock_scientific_runtime(monkeypatch, tmp_path, found)
    monkeypatch.setattr(
        'kinbot.anl.automatic.shutil.which', lambda value: found.get(value))
    checked = validate_automatic_environment({
        'cbh_rung': 1, 'composite_method': 'ANL0-F12',
        'me': 1, 'run_me': 1, 'run_messpf': 1,
        'mess_command': 'mess', 'messpf_command': 'messpf',
        'queuing': 'slurm',
    })
    assert checked['mess_command'] == '/opt/mess/bin/mess'
    assert checked['messpf_command'] == '/opt/mess/bin/messpf'
    assert checked['sbatch'] == '/usr/bin/sbatch'
    gaussian = {key: value for key, value in checked.items()
                if key.endswith('_gaussian')}
    assert gaussian
    assert set(gaussian.values()) == {'/opt/gaussian/g16/g16'}
    assert 'molpro' not in checked
    assert 'xcfour' not in checked
    assert 'dmrcc' not in checked


def test_automatic_environment_honors_solver_environment(monkeypatch, tmp_path):
    found = {
        '/site/mess/bin/mess': '/site/mess/bin/mess',
        '/site/mess/bin/messpf': '/site/mess/bin/messpf'}
    _mock_scientific_runtime(monkeypatch, tmp_path, found)
    monkeypatch.setenv('KINBOT_MESS_COMMAND', '/site/mess/bin/mess')
    monkeypatch.setenv('KINBOT_MESSPF_COMMAND', '/site/mess/bin/messpf')
    monkeypatch.setattr(
        'kinbot.anl.automatic.shutil.which', lambda value: found.get(value))
    checked = validate_automatic_environment({
        'composite_method': 'ANL0-F12', 'me': 1, 'run_me': 1,
        'run_messpf': 1, 'mess_command': 'missing-mess',
        'messpf_command': 'missing-messpf', 'queuing': 'local'})
    assert checked['mess_command'] == '/site/mess/bin/mess'
    assert checked['messpf_command'] == '/site/mess/bin/messpf'
    assert checked['atct_client'] == 'available'


def test_automatic_environment_finds_unambiguous_user_prefix(
        monkeypatch, tmp_path):
    prefix = tmp_path / '.local' / 'kinbot-mess' / 'bin'
    prefix.mkdir(parents=True)
    for name in ('mess', 'messpf'):
        executable = prefix / name
        executable.write_text('#!/bin/sh\n')
        executable.chmod(0o755)
    found = {}
    _mock_scientific_runtime(monkeypatch, tmp_path, found)
    monkeypatch.setenv('HOME', str(tmp_path))
    monkeypatch.delenv('KINBOT_MESS_COMMAND', raising=False)
    monkeypatch.delenv('KINBOT_MESSPF_COMMAND', raising=False)
    monkeypatch.setattr(
        'kinbot.anl.automatic.shutil.which', lambda value: found.get(value))
    parameters = {
        'composite_method': 'ANL0-F12', 'me': 1, 'run_me': 1,
        'run_messpf': 1, 'mess_command': 'mess',
        'messpf_command': 'messpf', 'queuing': 'local'}
    checked = validate_automatic_environment(parameters)
    assert checked['mess_command'] == str((prefix / 'mess').resolve())
    assert checked['messpf_command'] == str((prefix / 'messpf').resolve())
    assert parameters['mess_command'] == str((prefix / 'mess').resolve())
    assert parameters['messpf_command'] == str((prefix / 'messpf').resolve())


def test_automatic_environment_checks_rotdpy_before_qc(monkeypatch, tmp_path):
    provenance = {'git_revision': 'pinned'}
    found = {}
    _mock_scientific_runtime(monkeypatch, tmp_path, found)
    monkeypatch.setattr(
        'kinbot.anl.automatic.shutil.which', lambda value: found.get(value))
    monkeypatch.setattr(
        'kinbot.rotdpy.ensure_available', lambda: provenance)
    checked = validate_automatic_environment({
        'cbh_rung': 1, 'composite_method': 'ANL0-F12',
        'rotdpy_run': 1, 'me': 0, 'queuing': 'local',
    })
    assert checked['rotdpy'] == provenance


def test_cbh0_is_an_active_automatic_request(monkeypatch):
    assert automatic_requested({
        'cbh_rung': 0, 'composite_method': 'ANL0-F12',
        'anl_handoff_manifest': ''})
    assert not automatic_requested({
        'cbh_rung': 0, 'composite_method': 'ANL0-F12',
        'anl_handoff_manifest': 'complete.json'})
    assert not automatic_requested({
        'cbh_rung': 0, 'composite_method': 'ANL0-F12',
        'anl_handoff_manifest': '',
        'thermochemistry_model_output': 'reference-model.json'})
    assert automatic_requested({
        'cbh_rung': 0, 'composite_method': 'ANL1-F12',
        'anl_handoff_manifest': ''})
    monkeypatch.setattr('kinbot.anl.automatic.shutil.which', lambda value: None)
    with pytest.raises(FileNotFoundError, match='messpf_command'):
        validate_automatic_environment({
            'cbh_rung': 0, 'composite_method': 'ANL0',
            'anl_handoff_manifest': '', 'run_messpf': 1,
            'messpf_command': 'messpf', 'me': 0, 'queuing': 'local',
        })


def test_automatic_environment_rejects_missing_messpf_before_qc(monkeypatch):
    monkeypatch.setattr('kinbot.anl.automatic.shutil.which', lambda value: None)
    with pytest.raises(FileNotFoundError, match='messpf_command'):
        validate_automatic_environment({
            'cbh_rung': 1, 'composite_method': 'ANL0-F12',
            'me': 2, 'run_me': 0, 'run_messpf': 1,
            'messpf_command': 'messpf', 'queuing': 'local',
        })


def test_reference_input_is_an_ordinary_one_slot_kinbot_run(tmp_path):
    source = tmp_path / 'parent.json'
    source.write_text(json.dumps({
        'smiles': 'CC=O', 'charge': 0, 'mult': 1,
        'reaction_search': 1, 'me': 1, 'cbh_rung': 1,
        'composite_method': 'ANL0-F12',
        'fc_model_path': '/cache/uma.pt'}))
    directory = tmp_path / 'reference'
    directory.mkdir()
    output = _reference_input(
        source, {'smiles': 'C=O', 'charge': 0, 'multiplicity': 1},
        directory, directory / 'model.json')
    generated = json.loads(output.read_text())
    assert generated['smiles'] == 'C=O'
    assert generated['reaction_search'] == generated['me'] == 0
    assert generated['cbh_rung'] == 0
    assert generated['conformer_search'] == generated['high_level'] == 1
    assert generated['rotor_scan'] == 1
    assert generated['queue_job_limit'] == 1


def test_reference_input_uses_the_shared_automatic_branch_limit(tmp_path):
    source = tmp_path / 'parent.json'
    source.write_text(json.dumps({
        'smiles': 'CC=O', 'charge': 0, 'mult': 1,
        'reaction_search': 1, 'me': 1, 'cbh_rung': 1,
        'composite_method': 'ANL0-F12'}))
    directory = tmp_path / 'reference'
    directory.mkdir()
    output = _reference_input(
        source, {'smiles': 'C=O', 'charge': 0, 'multiplicity': 1},
        directory, directory / 'model.json', queue_job_limit=4)
    assert json.loads(output.read_text())['queue_job_limit'] == 4


def test_inherited_reference_controller_is_replaced_once_if_it_disappears(
        monkeypatch, tmp_path):
    source = tmp_path / 'parent.json'
    source.write_text(json.dumps({
        'smiles': 'CC=O', 'charge': 0, 'mult': 1,
        'reaction_search': 1, 'me': 0, 'cbh_rung': 1,
        'composite_method': 'ANL0-F12'}))
    root = tmp_path / 'references'
    directory = root / ('state_' + hashlib.sha256(
        b'C|0|1').hexdigest()[:16])
    directory.mkdir(parents=True)
    (directory / 'driver.pid').write_text('321\n')

    matches = iter((True, False))
    monkeypatch.setattr(
        'kinbot.anl.automatic._pid_matches_reference',
        lambda *args: next(matches))
    monkeypatch.setattr('kinbot.anl.automatic.os.kill', lambda *args: None)
    clock = iter((0., 61.))
    monkeypatch.setattr(
        'kinbot.anl.automatic.time.monotonic', lambda: next(clock))
    monkeypatch.setattr('kinbot.anl.automatic.time.sleep', lambda value: None)

    launched = []

    class FinishedProcess:
        pid = 654

        def poll(self):
            return 0

    def launch(*args, **kwargs):
        launched.append((args, kwargs))
        (directory / 'thermochemistry_model.json').write_text('{}')
        return FinishedProcess()

    monkeypatch.setattr('kinbot.anl.automatic.subprocess.Popen', launch)
    monkeypatch.setattr(
        'kinbot.anl.automatic._model_record',
        lambda model, entry, theory: {'record': str(model)})

    records = _run_reference_l2(
        [{'smiles': 'C', 'charge': 0, 'multiplicity': 1}],
        source, root, max_jobs=1, l2_theory={})

    assert len(launched) == 1
    assert (directory / 'driver.pid').read_text() == '654\n'
    assert records['C|0|1'] == {
        'record': str(directory / 'thermochemistry_model.json')}


def test_exported_reference_model_becomes_a_hash_bound_l2_source(tmp_path):
    atoms = molecule('CH4')
    from kinbot.anl.dispatch import _geometry_hash
    model = {
        'name': 'methane', 'smiles': 'C', 'charge': 0, 'multiplicity': 1,
        'symbols': atoms.get_chemical_symbols(),
        'positions': atoms.positions.tolist(),
        'frequencies_cm_inverse': [1000.] * 9,
        'source': {'electronic_energy_hartree': -40.,
                   'harmonic_zpe_hartree': .04,
                   'geometry_sha256': _geometry_hash(atoms)},
    }
    path = tmp_path / 'model.json'
    path.write_text(json.dumps({'schema': 1, 'species': {'C': model}}))
    record = _model_record(
        path, {'smiles': 'C', 'charge': 0, 'multiplicity': 1},
        {'calculator': 'gaussian', 'method': 'B2PLYP',
         'basis': 'cc-pVTZ', 'dispersion': 'GD3BJ',
         'composite_method': 'B2PLYP-D3BJ',
         'label': 'B2PLYP-D3(BJ)/cc-pVTZ'})
    assert record['molecule']['source']['provider'] == \
        'automatic-kinbot-reference-l2'
    assert len(record['record_sha256']) == 64
    assert record['energy_hartree'] == -40.
    assert record['molecule']['source']['geometry_hash_format'] == \
        'symbols-and-positions-json'

    # G16's accepted HSS reference was exported with a different localized
    # bond/radical assignment. Keep the requested reference and provenance.
    hss = Atoms('SSH', positions=[[0, 0, 0], [1.98, 0, 0], [2.3, 1.3, 0]])
    alternative = dict(model, smiles='S=[SH]', multiplicity=2,
                       symbols=list(hss.symbols), positions=hss.positions.tolist(),
                       source=dict(model['source'], geometry_sha256=_geometry_hash(hss)))
    other = tmp_path / 'hss.json'
    other.write_text(json.dumps({'schema': 1, 'species': {'S=[SH]': alternative}}))
    resolved = _model_record(other, {'smiles': '[S]S', 'charge': 0, 'multiplicity': 2},
                             record['l2_theory'])
    assert resolved['smiles'] == resolved['thermochemistry_model']['smiles'] == '[S]S'
    assert resolved['thermochemistry_model']['identity_reconciliation']['perceived_smiles'] == 'S=[SH]'
    assert resolved['thermochemistry_model']['source'] == alternative['source']
    with pytest.raises(ValueError, match='expected one model'):
        _model_record(other, {'smiles': '[S][S]', 'charge': 0, 'multiplicity': 3}, record['l2_theory'])
    alternative['smiles'] = '[S][S]'
    other.write_text(json.dumps({'schema': 1, 'species': {'S=[SH]': alternative}}))
    with pytest.raises(ValueError, match='model identity differs'):
        _model_record(other, {'smiles': '[S]S', 'charge': 0, 'multiplicity': 2}, record['l2_theory'])

    model['source']['geometry_sha256'] = hashlib.sha256(np.asarray(
        atoms.positions, dtype='<f8').tobytes()).hexdigest()
    path.write_text(json.dumps({'schema': 1, 'species': {'C': model}}))
    legacy = _model_record(
        path, {'smiles': 'C', 'charge': 0, 'multiplicity': 1},
        record['l2_theory'])
    assert legacy['molecule']['source']['geometry_hash_format'] == \
        'legacy-position-bytes'

    model['source']['geometry_sha256'] = '0' * 64
    path.write_text(json.dumps({'schema': 1, 'species': {'C': model}}))
    with pytest.raises(ValueError, match='geometry hash differs'):
        _model_record(
            path, {'smiles': 'C', 'charge': 0, 'multiplicity': 1},
            record['l2_theory'])


def test_terminal_atct_reference_needs_l2_model_without_composite_graph():
    entry = {'smiles': '[H][H]', 'charge': 0, 'multiplicity': 1}
    key = '[H][H]|0|1'
    source = {
        'thermochemistry_model': {
            'smiles': '[H][H]', 'charge': 0, 'multiplicity': 1,
            'symbols': ['H', 'H'],
            'positions': [[0., 0., 0.], [0., 0., .74]],
            'frequencies_cm_inverse': [4400.],
        },
    }
    models = _reference_thermochemistry_models(
        [entry], {key: source}, {})
    assert models['[H][H]']['frequencies_cm_inverse'] == [4400.]


def test_unaccepted_vpt2_result_falls_back_to_harmonic(monkeypatch):
    monkeypatch.setattr(
        'kinbot.anl.automatic._verified_task_result',
        lambda *args, **kwargs: (
            None, None, None, None, None, {'review_required': True}))
    assert _vpt2_zpe_components('run', state_id='state') is None


def test_automatic_final_report_is_self_contained_and_idempotent(tmp_path):
    automatic_dir = tmp_path / 'anl_automatic'
    automatic_dir.mkdir()
    me = tmp_path / 'me'
    me.mkdir()
    plan = automatic_dir / 'network_request.json'
    plan.write_text(json.dumps({
        'schema': 1, 'requested_method': 'ANL0-F12',
        'requested_cbh_rung': 1,
        'network': {'initial': 'root', 'reactions': [],
                    'stable_species': [{'routing_name': 'root'}],
                    'transition_states': {}},
    }))
    handoff = automatic_dir / 'handoff.json'
    handoff.write_text(json.dumps({'schema': 1, 'mode': 'cbh-anl'}))
    from hashlib import sha256
    handoff_hash = sha256(handoff.read_bytes()).hexdigest()
    (me / 'anl_handoff.json').write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'manifest_sha256': handoff_hash,
    }))
    execution = automatic_dir / 'execution.json'
    execution.write_text(json.dumps({
        'schema': 1, 'status': 'complete', 'plan': str(plan),
        'handoff_manifest': str(handoff), 'graphs': {}, 'formations': {},
    }))
    automatic = {
        'status': 'complete', 'plan': str(plan),
        'execution_record': str(execution),
    }
    first = finalize_automatic_run(
        automatic, {'me': 0, 'run_messpf': 0, 'rotdpy_run': 0},
        run_dir=tmp_path)
    second = finalize_automatic_run(
        automatic, {'me': 0, 'run_messpf': 0, 'rotdpy_run': 0},
        run_dir=tmp_path)
    assert first == second
    report = json.loads(Path(first).read_text())
    assert report['status'] == 'complete'
    assert report['handoff']['sha256'] == handoff_hash
    assert report['network']['stable_species'] == ['root']


def test_automatic_audit_handles_multistate_stages_and_only_selected_reactions(tmp_path):
    from kinbot.anl.automatic import _audit_rotd_results
    from kinbot.vrc_multistate import stage_channel
    from test_vrc_multistate import _completed_channel
    rotd = tmp_path / 'rotdPy'
    rotd.mkdir()
    manifest = _completed_channel(rotd)
    me = tmp_path / 'me'
    provenance = stage_channel(manifest, me, required_profile='interface')
    (me / 'rotdpy_sources.json').write_text(json.dumps({'reactions': {'r': provenance}}))
    plan = {'network': {'reactions': [{'name': 'r', 'barrierless': True},
                                      {'name': 'unselected', 'barrierless': True}]}}
    parameters = {'me': 1, 'rotdpy_mess_mode': 'required',
                  'vrc_tst_scan': {'parent': ['r']}}
    report = _audit_rotd_results(plan, parameters, tmp_path)
    assert set(report) == {'r'}
    assert report['r']['kind'] == 'multistate'
    assert {s['name'] for s in report['r']['bottlenecks']} == {'inner', 'outer'}
    assert sum(len(s['states']) for s in report['r']['bottlenecks']) == 4
    assert not (rotd / 'r.py').exists()  # only scalar physical-state leaves exist
    table = me / provenance['staged_bottlenecks'][0]['file']
    table.write_text(table.read_text() + '\n')
    with pytest.raises(ValueError, match='flux table changed'):
        _audit_rotd_results(plan, parameters, tmp_path)


def test_automatic_report_resolves_and_verifies_mess_paths_from_run_root(
        tmp_path):
    automatic_dir = tmp_path / 'anl_automatic'
    automatic_dir.mkdir()
    me = tmp_path / 'me'
    me.mkdir()
    plan = automatic_dir / 'network_request.json'
    plan.write_text(json.dumps({
        'schema': 1, 'requested_method': 'ANL0-F12',
        'requested_cbh_rung': 1,
        'network': {'initial': 'root', 'reactions': [],
                    'stable_species': [{'routing_name': 'root'}],
                    'transition_states': {}},
    }))
    handoff = automatic_dir / 'handoff.json'
    handoff.write_text(json.dumps({'schema': 1, 'mode': 'cbh-anl'}))
    handoff_hash = hashlib.sha256(handoff.read_bytes()).hexdigest()
    (me / 'anl_handoff.json').write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'manifest_sha256': handoff_hash,
    }))
    mess_input = me / 'mess_0000.inp'
    mess_output = me / 'mess_0000.out'
    mess_input.write_text('MESS input\n')
    mess_output.write_text('MESS output\n')
    mess_record = me / 'mess_execution.json'
    mess_record.write_text(json.dumps({
        'schema': 1, 'status': 'complete', 'queue': 'slurm',
        'command': ['mess'],
        'jobs': {'mess_0000': {
            'input': 'me/mess_0000.inp',
            'input_sha256': hashlib.sha256(
                mess_input.read_bytes()).hexdigest(),
            'output': 'me/mess_0000.out',
            'output_sha256': hashlib.sha256(
                mess_output.read_bytes()).hexdigest(),
            'output_bytes': mess_output.stat().st_size,
        }},
    }))
    execution = automatic_dir / 'execution.json'
    execution.write_text(json.dumps({
        'schema': 1, 'status': 'complete', 'plan': str(plan),
        'handoff_manifest': str(handoff), 'graphs': {}, 'formations': {},
    }))
    automatic = {
        'status': 'complete', 'plan': str(plan),
        'execution_record': str(execution),
    }
    result = finalize_automatic_run(
        automatic, {'me': 1, 'run_me': 1, 'run_messpf': 0,
                    'rotdpy_run': 0}, run_dir=tmp_path)
    report = json.loads(Path(result).read_text())
    assert report['mess']['status'] == 'complete'

    mess_input.write_text('changed\n')
    with pytest.raises(ValueError, match='MESS native input hash differs'):
        finalize_automatic_run(
            automatic, {'me': 1, 'run_me': 1, 'run_messpf': 0,
                        'rotdpy_run': 0}, run_dir=tmp_path)


def test_automatic_report_accepts_verified_network_with_no_reactions(
        tmp_path):
    automatic_dir = tmp_path / 'anl_automatic'
    automatic_dir.mkdir()
    me = tmp_path / 'me'
    me.mkdir()
    plan = automatic_dir / 'network_request.json'
    plan.write_text(json.dumps({
        'schema': 1, 'requested_method': 'ANL0-F12',
        'requested_cbh_rung': 1,
        'network': {'initial': 'root', 'reactions': [],
                    'stable_species': [{'routing_name': 'root'}],
                    'transition_states': {}},
    }))
    handoff = automatic_dir / 'handoff.json'
    handoff.write_text(json.dumps({'schema': 1, 'mode': 'cbh-anl'}))
    (me / 'anl_handoff.json').write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'manifest_sha256': hashlib.sha256(
            handoff.read_bytes()).hexdigest(),
    }))
    (me / 'mess_networks.json').write_text(json.dumps([{
        'stem': 'mess_0000', 'status': 'no_reactions',
    }]))
    (me / 'mess_execution.json').write_text(json.dumps({
        'schema': 1, 'status': 'complete', 'queue': 'slurm',
        'command': ['mess'], 'jobs': {},
    }))
    execution_record = automatic_dir / 'execution.json'
    execution_record.write_text(json.dumps({
        'schema': 1, 'status': 'complete', 'plan': str(plan),
        'handoff_manifest': str(handoff), 'graphs': {}, 'formations': {},
    }))
    result = finalize_automatic_run({
        'status': 'complete', 'plan': str(plan),
        'execution_record': str(execution_record),
    }, {'me': 1, 'run_me': 1, 'run_messpf': 0, 'rotdpy_run': 0},
       run_dir=tmp_path)
    report = json.loads(Path(result).read_text())
    assert report['mess']['status'] == 'no_reactions'


def test_exact_profile_order_preserves_harmonic_fallback_without_ambiguity():
    labels = {
        'L2:B2PLYP-D3BJ/cc-pVTZ:harmonic-ZPE': None,
        'L3:CCSD(T)/cc-pVTZ//B2PLYP-D3BJ/cc-pVTZ:harmonic-ZPE': None,
        'L3:CCSD(T)/cc-pVTZ//CCSD(T)/cc-pVTZ:harmonic-ZPE': None,
        'profiled:ANL0-F12:scaled-triples:harmonic-ZPE': None,
        'profiled:ANL0-F12:scaled-triples:B2PLYP-D3BJ': None,
        'profiled:ANL0:B2PLYP-D3BJ': None,
    }
    ordered = _ordered_methods({'species': labels})
    assert ordered[0].endswith('B2PLYP-D3BJ')
    assert ordered[1].endswith('harmonic-ZPE')
    assert ordered[2].startswith('profiled:ANL0:')
    assert '//CCSD(T)/cc-pVTZ' in ordered[3]
    assert ordered[-1].startswith('L2:')


def test_reduced_profiles_precede_lower_tiers_and_are_common_subset_complete():
    full = 'profiled:ANL0-F12:scaled-triples:B2PLYP-D3BJ'
    reduced_hoe = full + ':reduced[without=higher_order]'
    reduced_dboc = full + ':reduced[without=dboc]'
    reduced_two = full + ':reduced[without=higher_order+dboc]'
    conventional = 'profiled:ANL0:B2PLYP-D3BJ'
    ordered = _ordered_methods({'species': {
        conventional: None, reduced_two: None, reduced_dboc: None,
        reduced_hoe: None, full: None}})
    assert ordered == (
        full, reduced_hoe, reduced_dboc, reduced_two, conventional)

    omissions = _component_family_omissions({
        'higher_order': False,
        'dboc': True,
        'scalar_relativistic': True,
        'core_valence': False,
    })
    assert len(omissions) == 4
    assert omissions[0] == ('higher_order', 'core_valence')
    assert all('higher_order' in item and 'core_valence' in item
               for item in omissions)
    assert ('higher_order', 'dboc', 'scalar_relativistic',
            'core_valence') in omissions


def test_multi_graph_driver_waits_for_scheduler_recovery(monkeypatch):
    responses = iter([
        SchedulerUnavailableError('slurmctld unavailable'),
        {'tasks': {'energy': {'status': 'complete'}}},
    ])

    def advance(*args, **kwargs):
        response = next(responses)
        if isinstance(response, Exception):
            raise response
        return response

    monkeypatch.setattr('kinbot.anl.automatic.advance', advance)
    sleeps = []
    monkeypatch.setattr(
        'kinbot.anl.automatic.time.sleep', lambda value: sleeps.append(value))
    result = _drive_graphs({'state': Path('/graph')}, 1, interval=7.)
    assert result == {'state': {'tasks': {
        'energy': {'status': 'complete'}}}}
    assert sleeps == [7.]
