"""Automatic CBH/ANL planning starts from the accepted KinBot network."""

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase.build import molecule

from kinbot.anl.automatic import (DEFAULT_SIZE_POLICY, _model_record,
                                  _component_family_omissions,
                                  _ordered_methods, _reference_input,
                                  _reference_thermochemistry_models,
                                  _vpt2_zpe_components,
                                  automatic_requested,
                                  build_network_request,
                                  finalize_automatic_run, highest_tier,
                                  validate_automatic_environment)
from kinbot.anl.validation import automatic_composite_spec


def _species(name, smiles, atoms, multiplicity=1, energy=-1., zpe=.01):
    return SimpleNamespace(
        name=name, chemid=name, wellorts=1, smiles=smiles,
        atom=list(atoms), geom=[[float(i), 0., 0.] for i in range(len(atoms))],
        charge=0, mult=multiplicity, energy=energy, zpe=zpe,
        freq=[500., 1000.], source_job=name + '_high', source_row_id=7)


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
    conventional_ids = {task['id'] for task in conventional['tasks']}
    assert {'f12_tz', 'f12_qz', 'anl0_qz', 'anl0_5z'} <= set(f12_tasks)
    assert {'anl0_qz', 'anl0_5z'} <= conventional_ids
    assert not {'f12_tz', 'f12_qz'} & conventional_ids
    assert f12_tasks['anl0_qz']['failure_skips'] == ['anl0_5z']
    assert f12_tasks['anl0_5z']['failure_skips'] == ['anl0_qz']
    assert 'basis=aug-cc-pVQZ,h=cc-pVQZ' in \
        f12_tasks['anl0_qz']['input_template']


def test_automatic_environment_checks_final_solvers_before_qc(monkeypatch):
    found = {
        'mess': '/opt/mess/bin/mess', 'messpf': '/opt/mess/bin/messpf',
        'sbatch': '/usr/bin/sbatch', 'squeue': '/usr/bin/squeue'}
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


def test_automatic_environment_honors_solver_environment(monkeypatch):
    found = {
        '/site/mess/bin/mess': '/site/mess/bin/mess',
        '/site/mess/bin/messpf': '/site/mess/bin/messpf'}
    monkeypatch.setenv('KINBOT_MESS_COMMAND', '/site/mess/bin/mess')
    monkeypatch.setenv('KINBOT_MESSPF_COMMAND', '/site/mess/bin/messpf')
    monkeypatch.setattr(
        'kinbot.anl.automatic.shutil.which', lambda value: found.get(value))
    checked = validate_automatic_environment({
        'composite_method': 'ANL0-F12', 'me': 1, 'run_me': 1,
        'run_messpf': 1, 'mess_command': 'missing-mess',
        'messpf_command': 'missing-messpf', 'queuing': 'local'})
    assert checked == {
        'mess_command': '/site/mess/bin/mess',
        'messpf_command': '/site/mess/bin/messpf',
        'atct_client': 'available'}


def test_automatic_environment_checks_rotdpy_before_qc(monkeypatch):
    provenance = {'git_revision': 'pinned'}
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
    assert not automatic_requested({
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
