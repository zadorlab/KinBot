"""Molpro orbital reconnaissance is explicit, reviewable, and hash bound."""

import hashlib
import json
from pathlib import Path

import pytest

from kinbot.molpro_reconnaissance import audit, prepare, select_record


def _configuration(path):
    path.write_text(json.dumps({
        'schema': 1,
        'title': 'test_radical',
        'symbols': ['C', 'S', 'H'],
        'geometry_angstrom': [[0., 0., 0.], [1.8, 0., 0.], [0., 1., 0.]],
        'charge': 0,
        'multiplicity': 2,
        'bases': [{'name': 'avtz', 'basis': 'aug-cc-pVTZ'}],
        'anchors': [
            {'name': 'equilibrium'},
            {'name': 'dissociated', 'bond': [1, 2],
             'moving_atoms': [2], 'distance_angstrom': 5.0},
        ],
        'partition': 'day-long-cpu',
        'ppn': 8,
    }))


def test_prepare_renders_c1_rohf_jobs_and_displaced_anchor(tmp_path):
    source = tmp_path / 'input.json'
    _configuration(source)
    manifest_path = prepare(source, tmp_path / 'recon')
    manifest = json.loads(manifest_path.read_text())
    assert manifest['status'] == 'prepared'
    assert len(manifest['jobs']) == 2
    equilibrium = (tmp_path / 'recon' / 'avtz__equilibrium' /
                   'avtz__equilibrium.inp').read_text()
    displaced = (tmp_path / 'recon' / 'avtz__dissociated' /
                 'avtz__dissociated.inp').read_text()
    for deck in (equilibrium, displaced):
        assert 'symmetry,nosym' in deck
        assert 'orient,noorient' in deck
        assert 'basis=aug-cc-pVTZ' in deck
        assert '{rhf;wf,23,1,1;orbital,2100.2}' in deck
        assert 'gprint,orbitals,orben' in deck
        assert 'put,molden,' in deck
    sulfur_x = [line for line in displaced.splitlines()
                if line.startswith('S ')][0].split()[1]
    assert float(sulfur_x) == 5.0
    setup = (tmp_path / 'recon' / 'site_setup.sh').read_text()
    assert 'Generated from the environment visible during prepare' in setup
    assert 'KINBOT_BACKEND' in setup


def test_audit_requires_native_output_wfu_and_molden(tmp_path):
    source = tmp_path / 'input.json'
    _configuration(source)
    manifest_path = prepare(source, tmp_path / 'recon')
    for directory in (tmp_path / 'recon').glob('avtz__*'):
        stem = directory.name
        (directory / f'{stem}.out').write_text(
            'SETTING RECONDONE = 1\nMolpro calculation terminated\n')
        (directory / f'{stem}.wfu').write_bytes(b'wfu ' + stem.encode())
        (directory / f'{stem}.molden').write_text('[Molden Format]\n')
    result = audit(manifest_path)
    assert result['status'] == 'complete_review_required'
    assert all(job['wfu_sha256'] for job in result['jobs'])


def test_select_hash_binds_a_reviewed_active_orbital_record(tmp_path):
    wfu = tmp_path / 'reviewed.wfu'
    wfu.write_bytes(b'reviewed active orbitals')
    output = select_record(
        wfu, tmp_path / 'state.json', record='2140.2',
        active_electrons=5, active_orbitals=4,
        state_average=[
            {'multiplicity': 2, 'roots': 3},
            {'multiplicity': 4, 'roots': 3},
        ], selected_roots={2: 1, 4: 1},
        review_note='C-S sigma/sigma-star and sulfur 3p manifold reviewed')
    payload = json.loads(output.read_text())
    assert payload['active_space']['active_electrons'] == 5
    assert payload['state_model']['orbital_selection']['sha256'] == \
        hashlib.sha256(wfu.read_bytes()).hexdigest()
    assert payload['state_model']['selected_roots'] == {'2': 1, '4': 1}


def test_select_accumulates_basis_specific_reviewed_records(tmp_path):
    output = tmp_path / 'state.json'
    for basis in ('avtz-f12', 'avqz-f12'):
        wfu = tmp_path / f'{basis}.wfu'
        wfu.write_bytes(f'reviewed {basis}'.encode())
        select_record(
            wfu, output, record='2140.2', basis=basis,
            active_electrons=9, active_orbitals=7,
            state_average=[
                {'multiplicity': 2, 'roots': 3},
                {'multiplicity': 4, 'roots': 3},
            ], selected_roots={2: 1, 4: 1},
            review_note='same reviewed C-S/S-S manifold')
    payload = json.loads(output.read_text())
    assert payload['schema'] == 2
    assert set(payload['basis_models']) == {'avtz-f12', 'avqz-f12'}
    assert all(
        model['orbital_selection']['sha256']
        for model in payload['basis_models'].values())


@pytest.mark.parametrize('name,channels', [('ch3s', 1), ('ch3ss', 2)])
def test_sulfur_recon_covers_all_eight_bases_and_derives_channel_states(tmp_path, name, channels):
    source = Path(__file__).parents[1] / 'examples' / 'anl' / f'{name}_anl1_f12' / 'orbital_reconnaissance.json'
    manifest = json.loads(prepare(source, tmp_path / name).read_text())
    jobs = manifest['jobs']
    assert len({job['basis'] for job in jobs}) == 8
    assert len({job['channel'] for job in jobs if job['channel']}) == channels
    assert max(job['distance_angstrom'] for job in jobs) >= 30.
    for job in jobs:
        assert job['spin_context']['multiplicities'] == [2, 4]
        roots = job['spin_context']['roots']
        channel = job['channel'] or job['orientation']['channel']
        count = (1 if channel == 'cs' else 6) if name == 'ch3ss' else 3
        assert roots == {'2': count, '4': count}
        deck = (Path(job['directory']) / (job['name'] + '.inp')).read_text()
        assert '{multi,cas-ci;' in deck
        assert 'orbital,2140.2,type=natural' in deck
        assert f'basis={job["basis"]}\n' in deck
        assert '(rev2)' not in deck
        from kinbot.molpro import validate_input_blocks
        validate_input_blocks(deck, job['name'])
    job = jobs[0]
    batch = (Path(job['directory']) / (job['name'] + '.sbatch')).read_text()
    assert batch.index('#SBATCH --mem=') < batch.index('set -euo pipefail')
    assert batch.index('export KINBOT_BACKEND=molpro') < batch.index('source ')
    assert 'molpro_reconnaissance run' in (tmp_path / name / 'submit.sh').read_text()


def test_generated_rohf_block_is_balanced_and_malformed_deck_is_rejected(tmp_path, monkeypatch):
    from kinbot import molpro_reconnaissance as recon
    from kinbot.molpro import validate_input_blocks
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    data = recon._load(source)
    deck = recon._deck(data, data['bases'][0]['basis'], data['anchors'][0], 'probe')
    validate_input_blocks(deck)
    with pytest.raises(ValueError, match='unmatched closing brace'):
        validate_input_blocks(deck.replace('orbital,2100.2}', 'orbital,2100.2}}'))
    original = recon._deck
    monkeypatch.setattr(recon, '_deck', lambda *args, **kwargs: original(*args, **kwargs) + '}')
    with pytest.raises(ValueError, match='unmatched closing brace'):
        recon.prepare(source, tmp_path / 'invalid')


def test_execute_selects_its_own_wavefunction_repository(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    path = recon.prepare(source, tmp_path / 'recon')
    manifest = json.loads(path.read_text())
    job = manifest['jobs'][0]
    monkeypatch.setattr(recon, 'qc_runtime_environment', lambda *args, **kwargs: ({}, {'scratch_directory': str(tmp_path)}))
    monkeypatch.setattr(recon, 'cleanup_qc_runtime', lambda *args: None)
    def execute(command, **kwargs):
        assert command[command.index('-W') + 1] == job['directory']
        assert kwargs['cwd'] == Path(job['directory'])
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', execute)
    recon.execute(path, job['name'])
    assert json.loads((Path(job['directory']) / 'execution.json').read_text())['status'] == 'complete'


def test_first_native_probe_gates_other_submissions_and_reports_its_failure(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _configuration(source)
    path = recon.prepare(source, tmp_path / 'recon')
    submitted = []
    def scheduler(command, **kwargs):
        if command[0] == 'squeue':
            return SimpleNamespace(stdout='', returncode=0)
        submitted.append(command)
        return SimpleNamespace(stdout='42\n', stderr='', returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    report = recon.submit(path, max_jobs=10)
    assert len(submitted) == 1
    first = report['jobs'][0]
    (Path(first['directory']) / 'execution.json').write_text(json.dumps({
        'status': 'failed', 'error': 'Molpro: UNMATCHED }'}))
    with pytest.raises(RuntimeError, match='UNMATCHED'):
        recon.submit(path, max_jobs=10)
    assert len(submitted) == 1


def _coordinate_configuration(source):
    _configuration(source)
    data = json.loads(source.read_text())
    data.pop('anchors')
    data.update(probe_active_space=True, coordinates=[{
        'name': 'cs', 'bond': [1, 2], 'distances_angstrom': [1.8, 2., 3., 30.],
    }])
    source.write_text(json.dumps(data))


def _complete_probes(manifest_path):
    manifest = json.loads(manifest_path.read_text())
    prior_hashes = {}
    for index, job in enumerate(manifest['jobs']):
        directory, stem = Path(job['directory']), job['name']
        if job['orbital_predecessor']:
            digest = prior_hashes[job['orbital_predecessor']]
            inp = directory / f'{stem}.inp'
            inp.write_text(inp.read_text().replace('{{PREVIOUS_WFU_SHA256}}', digest))
            job['initial_wfu_sha256'] = digest
            job['input_sha256'] = hashlib.sha256(inp.read_bytes()).hexdigest()
        energy = -10. if index == 0 else -9.9
        (directory / f'{stem}.out').write_text(
            f'SETTING RECONENERGY = {energy} AU\n'
            + ''.join(f'SETTING RECONR{root} = {energy + .1 * root} AU\n'
                      for root in range(1, sum(job['spin_context']['roots'].values()) + 1))
            + 'SETTING RECONDONE = 1\nMolpro calculation terminated\n')
        (directory / f'{stem}.wfu').write_bytes(stem.encode())
        prior_hashes[stem] = hashlib.sha256(stem.encode()).hexdigest()
        (directory / f'{stem}.molden').write_text('[Molden Format]\n')
        (directory / f'{stem}.natural.molden').write_text(
            '[MO]\n' + 'Occup=2\n' * 11 + 'Occup=1\n')
    manifest_path.write_text(json.dumps(manifest))
    return manifest


def test_refinement_preserves_completed_checkpoint_chain_and_exposes_candidates(tmp_path):
    from kinbot.molpro_reconnaissance import refine
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    manifest_path = prepare(source, tmp_path / 'recon')
    before = _complete_probes(manifest_path)
    result = audit(manifest_path)
    assert result['status'] == 'coverage_refinement_required'
    coverage = next(iter(result['coverage'].values()))
    assert coverage['asymptote_plateau']
    assert coverage['occupation_candidates_by_geometry'][0]['orbitals_one_based'] == [12]
    refine(manifest_path)
    refined = json.loads(manifest_path.read_text())
    by_name = {job['name']: job for job in refined['jobs']}
    for job in before['jobs']:
        assert by_name[job['name']]['orbital_predecessor'] == job['orbital_predecessor']
        assert by_name[job['name']]['input_sha256'] == job['input_sha256']
    assert len(refined['jobs']) > len(before['jobs'])


def test_changed_input_is_not_reblessed_by_prepare_and_failed_audit_clears_evidence(tmp_path):
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    path = prepare(source, tmp_path / 'recon')
    manifest = _complete_probes(path)
    audit(path)
    first = manifest['jobs'][0]
    inp = Path(first['directory']) / (first['name'] + '.inp')
    inp.write_text(inp.read_text() + '\nchanged\n')
    with pytest.raises(ValueError, match='input changed'):
        prepare(source, tmp_path / 'recon')
    (inp.parent / (first['name'] + '.out')).unlink()
    failed = audit(path)['jobs'][0]
    assert failed['status'] == 'incomplete'
    assert 'energy_hartree' not in failed
    assert 'wfu_sha256' not in failed


def test_light_atom_dissociation_includes_its_antibonding_projection(tmp_path):
    from kinbot.molpro_reconnaissance import _load
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    data = json.loads(source.read_text())
    data['coordinates'][0]['bond'] = [1, 3]
    source.write_text(json.dumps(data))
    space = _load(source)['_probe_space']
    assert any(center['atom'] == 3 for center in space['orbital_selection']['centers'])
    assert space['active_orbitals'] == 9


def test_second_controller_cannot_rewrite_submission_journal(tmp_path):
    from kinbot.molpro_reconnaissance import _manifest_lock, submit
    source = tmp_path / 'input.json'
    _configuration(source)
    manifest = prepare(source, tmp_path / 'recon')
    with _manifest_lock(manifest.parent):
        with pytest.raises(RuntimeError, match='Another reconnaissance controller'):
            submit(manifest)


def test_parent_spin_is_first_in_probe_energy_array(tmp_path):
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    data = json.loads(source.read_text())
    data['multiplicity'] = 4
    source.write_text(json.dumps(data))
    manifest = json.loads(prepare(source, tmp_path / 'recon').read_text())
    job = manifest['jobs'][0]
    deck = (Path(job['directory']) / (job['name'] + '.inp')).read_text()
    multi = deck.split('{multi,cas-ci;', 1)[1]
    assert multi.index('wf,23,1,3') < multi.index('wf,23,1,1')
    assert 'reconr1=energy(1)' in deck
    assert job['root_states'] == [{'multiplicity': 4, 'root': 1},
                                  {'multiplicity': 2, 'root': 1}]


def test_asymptote_checks_every_root_even_when_parent_energy_is_flat(tmp_path):
    from kinbot.molpro_reconnaissance import _coverage
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    path = prepare(source, tmp_path / 'recon')
    _complete_probes(path)
    report = audit(path)
    coverage = next(iter(report['coverage'].values()))
    assert coverage['asymptote_plateau']
    states = coverage['state_energy_profiles'][-1]['states']
    assert states[0]['multiplicity'] == 2
    assert states[1]['multiplicity'] == 4
    # Only a higher root is still drifting. Ground-state-only plateau
    # detection previously reported this asymptote as settled.
    report['jobs'][-1]['root_energies_hartree'][-1] += .001
    coverage = next(iter(_coverage(report).values()))
    assert not coverage['asymptote_plateau']
    assert coverage['suggested_far_distance_angstrom'] > 30.


def test_orientation_probes_preserve_bond_and_fragment_geometry(tmp_path):
    import numpy as np
    from kinbot.molpro_reconnaissance import _load, _anchor_geometry
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    data = json.loads(source.read_text())
    data['orientation_probes'] = {'distance_scales': [1., 2.], 'angles_degrees': [60.],
                                  'axes': ['perpendicular_1', 'perpendicular_2'],
                                  'include_asymptote': True}
    source.write_text(json.dumps(data))
    loaded = _load(source)
    probes = [anchor for anchor in loaded['anchors'] if anchor.get('_orientation')]
    assert probes and max(anchor['distance_angstrom'] for anchor in probes) == 30.
    for anchor in probes:
        geometry = _anchor_geometry(loaded, anchor)
        assert np.linalg.norm(geometry[0] - geometry[1]) == pytest.approx(anchor['distance_angstrom'])
        assert np.linalg.norm(geometry[0] - geometry[2]) == pytest.approx(1.)
    manifest = json.loads(prepare(source, tmp_path / 'recon').read_text())
    assert any(job['orientation'] for job in manifest['jobs'])


def test_ionic_dissociation_uses_explicit_charge_partition_not_neutral_guess(tmp_path):
    from kinbot.molpro_reconnaissance import _load
    source = tmp_path / 'ionic.json'
    data = {'schema': 1, 'title': 'methane_cation', 'symbols': ['C', 'H', 'H', 'H', 'H'],
            'geometry_angstrom': [[0, 0, 0], [0, 0, 1.09], [1.03, 0, -.36],
                                 [-.51, .89, -.36], [-.51, -.89, -.36]],
            'charge': 1, 'multiplicity': 2,
            'bases': [{'name': 'test', 'basis': 'aug-cc-pVDZ'}],
            'coordinates': [{'name': 'ch', 'bond': [1, 2],
                             'distances_angstrom': [1.09, 2., 4., 30.]}]}
    source.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='explicit fragment'):
        _load(source)
    data['coordinates'][0]['fragment_states'] = [
        {'atom': 1, 'charge': 0, 'multiplicity': 2, 'smiles': '[CH3]'},
        {'atom': 2, 'charge': 1, 'multiplicity': 1, 'smiles': '[H+]'},
    ]
    source.write_text(json.dumps(data))
    context = _load(source)['anchors'][0]['_context']
    assert context['multiplicities'] == [2]
    assert sorted(fragment['charge'] for fragment in context['fragment_states']) == [0, 1]


def test_completed_controller_waits_for_allocations_before_next_system(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    path = tmp_path / 'reconnaissance.json'
    result = {'status': 'complete_review_required', 'jobs': [
        {'job_id': '42', 'name': 'probe', 'status': 'complete'}]}
    monkeypatch.setattr(recon, '_submit', lambda *args, **kwargs: result)
    queue = iter(['42\n', ''])
    monkeypatch.setattr(recon.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout=next(queue)))
    sleeps = []
    monkeypatch.setattr(recon.time, 'sleep', lambda delay: sleeps.append(delay))
    assert recon._run(path) is result
    assert sleeps == [20.]
