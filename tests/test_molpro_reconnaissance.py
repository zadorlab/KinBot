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


@pytest.mark.parametrize('ending', [
    'NO CONVERGENCE! MAXIMAL NUMBER OF ITERATIONS REACHED!',
    'No convergence in MULTI', 'MCSCF not converged',
    'No orbitals optimized', ''])
def test_casscf_requires_positive_native_convergence_even_after_normal_exit(ending):
    from kinbot.molpro_reconnaissance import convergence_diagnostics
    text = ('PROGRAM * RHF-SCF\nCONVERGENCE REACHED!\nPROGRAM * MULTI\n'
            + ending + '\nSETTING RECONDONE = 1\nMolpro calculation terminated\n')
    assert convergence_diagnostics(text, {'probe_method': 'casscf'})['errors']


@pytest.mark.parametrize('failure', ['?No convergence in rhfpr',
    'NO CONVERGENCE! MAXIMAL NUMBER OF ITERATIONS REACHED!'])
def test_rhf_failure_is_only_a_warning_after_converged_casscf(failure):
    from kinbot.molpro_reconnaissance import convergence_diagnostics
    text = ('Program * Restricted Hartree-Fock\n' + failure
            + '\nPROGRAM * MULTI\nCONVERGENCE REACHED!\n')
    job = {'probe_method': 'casscf', 'natural_orbital_probe': True, 'avas_reference': 'rhf'}
    result = convergence_diagnostics(text, job)
    assert not result['errors'] and result['warnings']
    job['probe_method'] = 'casci'
    assert convergence_diagnostics(text, job)['errors']
    assert convergence_diagnostics(text, {'natural_orbital_probe': False})['errors']


def test_rhf_so_failure_is_not_hidden_by_successful_frozen_casci():
    from kinbot.molpro_reconnaissance import convergence_diagnostics
    text = ('Program * Restricted Hartree-Fock\n'
            'NO CONVERGENCE! MAXIMAL NUMBER OF ITERATIONS REACHED!\n'
            '!RHF STATE 1.1 Energy -834.421740430818\n'
            'PROGRAM * MULTI (Direct Multiconfiguration SCF)\nNo orbital optimization\n'
            'SETTING RECONDONE = 1\nMolpro calculation terminated\n')
    result = convergence_diagnostics(text, {'probe_method': 'casci',
        'natural_orbital_probe': True, 'avas_reference': 'rhf'})
    assert result['errors'] == ['RHF starting calculation did not converge.']


def test_so_sci_failure_reports_its_actual_threshold_after_converged_rhf():
    from kinbot.molpro_reconnaissance import convergence_diagnostics
    text = ('Program * Restricted Hartree-Fock\nCONVERGENCE REACHED!\n'
            'PROGRAM * MULTI (Direct Multiconfiguration SCF)\n'
            'Convergence thresholds  0.10E-04 (gradient)  0.10E-11 (energy)\n'
            'First-order MCSCF: SO-SCI\n'
            'NO CONVERGENCE! MAXIMAL NUMBER OF ITERATIONS REACHED!\n'
            'SETTING RECONDONE = 1\nMolpro calculation terminated\n')
    result = convergence_diagnostics(text, {'probe_method': 'casscf',
        'natural_orbital_probe': True, 'avas_reference': 'rhf'})
    assert len(result['errors']) == 1 and not result['warnings']
    assert '0.10E-11 (energy)' in result['errors'][0]


@pytest.mark.parametrize('solver', ['second_order', 'so_sci'])
@pytest.mark.parametrize('custom', [False, True])
def test_casscf_sets_local_energy_accuracy_independently_of_solver(tmp_path, solver, custom):
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    data = json.loads(source.read_text())
    data.update(probe_method='casscf', scf_convergence=1.e-10,
                casscf_options={'solver': solver})
    if custom:
        data['casscf_options'].update(energy_convergence=1.e-9, gradient_convergence=1.e-6)
    source.write_text(json.dumps(data))
    data = recon._load(source)
    for previous in (None, 'previous'):
        deck = recon._deck(data, 'avtz-f12', data['anchors'][0], 'probe', previous=previous)
        before, multi = deck.split('{multi,', 1)
        assert 'gthresh,energy=1.d-10,coeff=1.d-10' in before
        assert 'accuracy,' not in before
        requested = 'accuracy,energy=1.d-9,gradient=1.d-6' if custom else 'accuracy,energy=1.d-10'
        assert requested + '\n' in multi
        assert multi.count('accuracy,') == 1
        recon.validate_input_blocks(deck)


@pytest.mark.parametrize('key', ['energy_convergence', 'gradient_convergence'])
@pytest.mark.parametrize('value', [True, 0., -1.e-8, float('nan'), float('inf'), '1e-8', 3.e-8])
def test_casscf_accuracy_rejects_invalid_or_unrepresentable_targets(tmp_path, key, value):
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    data = json.loads(source.read_text())
    data.update(probe_method='casscf', casscf_options={key: value})
    source.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='casscf_options.' + key):
        recon._load(source)


@pytest.mark.parametrize('verbosity', [True, -1, 5, 1.5, '4'])
def test_casscf_rejects_invalid_diagnostic_verbosity(tmp_path, verbosity):
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    data = json.loads(source.read_text())
    data.update(probe_method='casscf', casscf_options={'verbosity': verbosity})
    source.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='casscf_options.verbosity'):
        recon._load(source)


@pytest.mark.parametrize('execution_format', ['none', 'legacy', 'wrong_hash',
    'missing_hash', 'wrong_command', 'nonzero_returncode', 'partial_modern',
    'null_native_path', 'missing_native_hash', 'wrong_native_hash'])
def test_real_blodgett_molden_pilot_execution_evidence(tmp_path, monkeypatch, execution_format):
    source = (Path(__file__).parents[1] / 'examples/anl/ch3s_anl1_f12/orbital_reconnaissance.json')
    path = prepare(source, tmp_path / 'recon')
    job = json.loads(path.read_text())['jobs'][0]
    directory, name = Path(job['directory']), job['name']
    energies = [-436.30070077, -436.30070077, -436.14287849,
                -435.90611238, -435.76240012, -435.76240011]
    (directory / (name + '.out')).write_text(
        ' SETTING RECONENERGY    =      -436.30070077  HARTREE\n'
        + ''.join(f' SETTING RECONR{index}        =      {energy:.8f}  HARTREE\n'
                  for index, energy in enumerate(energies, 1))
        + ' SETTING RECONDONE      =         1.00000000\n'
        ' Molpro calculation terminated\n')
    (directory / (name + '.wfu')).write_bytes(b'completed native WFU')
    (directory / (name + '.molden')).write_text('[Molden Format]\n')
    occupations = [2.] * 6 + [1.99814, 1.99213, 1.99109, 1.99109,
                              1.50808, 1.50808, 1.49864, .51273]
    (directory / (name + '.natural.molden')).write_text(
        '[MO]\n' + ''.join(f'Occup= {value}\n' for value in occupations))
    if execution_format != 'none':
        # Schema emitted before native XML input copies were introduced.
        execution = {'schema': 1, 'status': 'complete', 'returncode': 0,
                     'command': ['molpro', '-d', str(tmp_path / 'scratch'),
                                 '-W', str(directory), '-t', '1', '-n', '8', name + '.inp'],
                     'input_sha256': job['input_sha256'], 'runtime': {}}
        if execution_format == 'wrong_hash':
            execution['input_sha256'] = '0' * 64
        elif execution_format == 'missing_hash':
            execution.pop('input_sha256')
        elif execution_format == 'wrong_command':
            execution['command'][-1] = 'other.inp'
        elif execution_format == 'nonzero_returncode':
            execution['returncode'] = 1
        elif execution_format == 'partial_modern':
            execution['input_transform'] = 'native_xml_compact_title_v2'
        elif execution_format == 'null_native_path':
            execution['native_input'] = None
        elif execution_format in ('missing_native_hash', 'wrong_native_hash'):
            execution['native_input'] = str(directory / (name + '.inp'))
            if execution_format == 'wrong_native_hash':
                execution['native_input_sha256'] = '0' * 64
        (directory / 'execution.json').write_text(json.dumps(execution))
    result = audit(path)['jobs'][0]
    if execution_format not in ('none', 'legacy'):
        assert result['status'] == 'incomplete'
        assert result['validation_errors']
        assert 'wfu_sha256' not in result
        return
    assert result['status'] == 'complete'
    assert result['root_energies_hartree'] == energies
    assert result['state_energies_hartree'][3]['multiplicity'] == 4
    if execution_format == 'legacy':
        from types import SimpleNamespace
        from kinbot import molpro_reconnaissance as recon
        # Reproduce the stopped pilot manifest after the erroneous rejection.
        # A restart must return completion without submitting or rewriting it.
        manifest = json.loads(path.read_text())
        manifest['jobs'] = [dict(result, status='incomplete', job_id='52907000', missing_since=0.)]
        manifest.pop('scan_policy', None)
        manifest['status'] = 'incomplete'
        path.write_text(json.dumps(manifest))
        evidence = {file: file.read_bytes() for file in directory.iterdir()
                    if file.suffix in ('.inp', '.wfu') or file.name == 'execution.json'}
        def scheduler(command, **kwargs):
            assert command[0] == 'squeue', 'Completed legacy pilot must never be resubmitted'
            return SimpleNamespace(stdout='', returncode=0)
        monkeypatch.setattr(recon.subprocess, 'run', scheduler)
        resumed = recon.run(path)
        assert resumed['status'] == 'complete_review_required'
        assert 'missing_since' not in resumed['jobs'][0]
        assert all(file.read_bytes() == original for file, original in evidence.items())


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
def test_sulfur_recon_covers_production_basis_and_derives_channel_states(tmp_path, name, channels):
    source = Path(__file__).parents[1] / 'examples' / 'anl' / f'{name}_anl1_f12' / 'orbital_reconnaissance.json'
    manifest = json.loads(prepare(source, tmp_path / name).read_text())
    jobs = manifest['jobs']
    assert {job['basis'] for job in jobs} == {'avtz-f12'}
    assert len(jobs) == (31 if name == 'ch3s' else 75)
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


def test_launcher_stderr_and_host_are_reported_when_no_native_output_exists(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    path = recon.prepare(source, tmp_path / 'recon')
    job = json.loads(path.read_text())['jobs'][0]
    directory = Path(job['directory'])
    monkeypatch.setattr(recon, 'qc_runtime_environment', lambda *args, **kw: ({}, {'scratch_directory': str(tmp_path)}))
    monkeypatch.setattr(recon, 'cleanup_qc_runtime', lambda *args: None)
    message = 'mktemp: failed to create directory: Input/output error'
    calls = []
    def failed(command, **kwargs):
        calls.append(command)
        kwargs['stderr'].write(message + '\n')
        kwargs['stdout'].write('Molpro launcher\n')
        return SimpleNamespace(returncode=1)
    monkeypatch.setattr(recon.subprocess, 'run', failed)
    with pytest.raises(RuntimeError, match='Input/output error'):
        recon.execute(path, job['name'])
    execution = json.loads((directory / 'execution.json').read_text())
    assert execution['host']
    assert execution['status'] == 'failed' and execution['returncode'] == 1
    assert 'No native output' in execution['error']
    assert 'Launcher stderr' in execution['error'] and 'Molpro launcher' in execution['error']
    assert Path(execution['launcher_stderr']).read_text() == message + '\n'
    assert len(calls) == 1  # A launcher error does not silently retry a calculation.


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


def test_uncertain_later_submission_blocks_every_new_probe(tmp_path, monkeypatch):
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _configuration(source)
    path = prepare(source, tmp_path / 'recon')
    manifest = json.loads(path.read_text())
    manifest['jobs'][-1]['submission_status'] = 'submitting'
    path.write_text(json.dumps(manifest))
    from types import SimpleNamespace
    submitted = []
    def scheduler(command, **kwargs):
        if command[0] == 'sbatch':
            submitted.append(command)
        return SimpleNamespace(stdout='42' if command[0] == 'sbatch' else '',
                               returncode=0, stderr='')
    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    with pytest.raises(RuntimeError, match='reconciliation'):
        recon.submit(path, max_jobs=10)
    assert submitted == []


def test_submission_timeout_keeps_uncertain_journal_and_prevents_retry(tmp_path, monkeypatch):
    import subprocess
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _configuration(source)
    path = prepare(source, tmp_path / 'recon')
    submitted = []
    def scheduler(command, **kwargs):
        assert kwargs['timeout'] == 30
        if command[0] == 'squeue':
            return SimpleNamespace(stdout='', returncode=0)
        submitted.append(command)
        raise subprocess.TimeoutExpired(command, 30)
    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    with pytest.raises(RuntimeError, match='uncertain'):
        recon.submit(path)
    assert json.loads(path.read_text())['jobs'][0]['submission_status'] == 'submitting'
    with pytest.raises(RuntimeError, match='reconciliation'):
        recon.submit(path)
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


def test_live_predecessor_is_not_consumed_or_hashed_during_retry(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    path = prepare(source, tmp_path / 'recon')
    manifest = _complete_probes(path)
    first = manifest['jobs'][0]
    first['job_id'] = '42'
    first['submission_status'] = 'submitted'
    # Valid-looking artifacts can belong to an attempt still owned by the
    # guard. The following checkpoint-dependent job must wait for queue exit.
    for job in manifest['jobs'][1:]:
        directory = Path(job['directory'])
        (directory / (job['name'] + '.out')).unlink()
        (directory / (job['name'] + '.wfu')).unlink()
    path.write_text(json.dumps(manifest))
    submitted = []
    def scheduler(command, **kwargs):
        if command[0] == 'squeue':
            return SimpleNamespace(stdout='42\n', returncode=0)
        submitted.append(command)
        return SimpleNamespace(stdout='43\n', stderr='', returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    report = recon.submit(path, max_jobs=10)
    assert report['jobs'][0]['status'] == 'incomplete'
    assert 'wfu_sha256' not in report['jobs'][0]
    assert submitted == []


@pytest.mark.parametrize('status', ['running', 'failed'])
def test_legacy_exports_do_not_override_unfinished_native_execution(tmp_path, status):
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    path = prepare(source, tmp_path / 'recon')
    manifest = _complete_probes(path)
    directory = Path(manifest['jobs'][0]['directory'])
    (directory / 'execution.json').write_text(json.dumps({'status': status}))
    report = audit(path)
    assert report['jobs'][0]['status'] == 'incomplete'
    assert 'wfu_sha256' not in report['jobs'][0]


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


@pytest.mark.parametrize('limit,startup_checked,publication_lag,chain', [
    (20, True, False, False), (3, True, False, False), (1, True, False, False),
    (3, True, True, False), (20, False, False, False), (20, True, False, True)])
def test_shared_campaign_capacity_reuses_live_jobs_and_keeps_dependencies(
        tmp_path, monkeypatch, limit, startup_checked, publication_lag, chain):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    paths = []
    for index in range(2):
        source = tmp_path / f'input{index}.json'
        (_coordinate_configuration if chain else _configuration)(source)
        paths.append(prepare(source, tmp_path / f'recon{index}'))
    existing = json.loads(paths[0].read_text())
    existing['jobs'][0].update(job_id='7', submission_status='submitted', status='incomplete')
    paths[0].write_text(json.dumps(existing))
    live, outstanding = {'7'}, {'7'}
    submissions, peaks = [], []
    jobs_count = sum(len(json.loads(p.read_text())['jobs']) for p in paths)

    def scheduler(command, **kwargs):
        if command[0] == 'squeue':
            return SimpleNamespace(stdout='\n'.join(sorted(live)), returncode=0)
        assert command[0] == 'sbatch'
        directory = Path(kwargs['cwd'])
        manifest = json.loads((directory.parent / 'reconnaissance.json').read_text())
        job = next(j for j in manifest['jobs'] if j['directory'] == str(directory))
        if job['orbital_predecessor']:
            prior = next(j for j in manifest['jobs'] if j['name'] == job['orbital_predecessor'])
            assert prior['status'] == 'complete' and prior['job_id'] not in outstanding
        identifier = str(100 + len(submissions))
        submissions.append((directory.parent.name, job['name']))
        outstanding.add(identifier)
        live.add(identifier)
        assert len(outstanding) <= limit
        return SimpleNamespace(stdout=identifier, stderr='', returncode=0)

    def advance(_):
        peaks.append(len(outstanding))
        assert len(peaks) < 20, 'scheduler failed to make progress'
        if publication_lag and len(peaks) == 1:
            live.clear()  # Slurm no longer shows jobs; final outputs are not published yet.
            return
        if publication_lag and len(peaks) == 2:
            assert len(submissions) == limit - 1  # Missing allocations still reserve slots.
        for path in paths:
            for job in json.loads(path.read_text())['jobs']:
                if job.get('job_id') not in outstanding:
                    continue
                directory, name = Path(job['directory']), job['name']
                (directory / (name + '.out')).write_text(
                    'SETTING RECONENERGY = -10.0 AU\n'
                    + ''.join(f'SETTING RECONR{index} = -10.0 AU\n'
                              for index, _ in enumerate(job['root_states'], 1))
                    + 'SETTING RECONDONE = 1\nMolpro calculation terminated\n')
                (directory / (name + '.wfu')).write_bytes(name.encode())
                (directory / (name + '.molden')).write_text('[Molden Format]\n')
                if job['natural_orbital_probe']:
                    (directory / (name + '.natural.molden')).write_text(
                        '[MO]\n' + 'Occup=2\n' * 11 + 'Occup=1\n')
        outstanding.clear()
        live.clear()

    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    monkeypatch.setattr(recon.time, 'sleep', advance)
    reports = recon.run_many(paths, max_jobs=limit, startup_checked=startup_checked)
    assert all(j['status'] == 'complete' for m in reports for j in m['jobs'])
    assert len(submissions) == jobs_count - 1  # Adopt, never resubmit, the live job.
    assert reports[0]['jobs'][0]['job_id'] == '7'
    assert peaks[0] == (min(limit, jobs_count) if startup_checked and not chain else 2)


@pytest.mark.parametrize('failure', ['uncertain', 'scheduler'])
def test_shared_scheduler_fails_before_any_submission(tmp_path, monkeypatch, failure):
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _configuration(source)
    paths = [prepare(source, tmp_path / f'recon{i}') for i in range(2)]
    if failure == 'uncertain':
        manifest = json.loads(paths[1].read_text())
        manifest['jobs'][0]['submission_status'] = 'submitting'
        paths[1].write_text(json.dumps(manifest))
    calls = []
    def scheduler(command, **kwargs):
        calls.append(command[0])
        if command[0] != 'squeue':
            pytest.fail('submitted without knowing campaign capacity')
        raise RuntimeError('scheduler unavailable')
    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    with pytest.raises(RuntimeError, match='reconciliation|scheduler unavailable'):
        recon.run_many(paths, startup_checked=True)
    assert calls == ([] if failure == 'uncertain' else ['squeue'])


def test_shared_capacity_includes_allocations_from_another_phase(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _configuration(source)
    paths = [prepare(source, tmp_path / f'recon{i}') for i in range(3)]
    outside = json.loads(paths[2].read_text())
    outside['jobs'][0].update(job_id='77', status='incomplete', submission_status='submitted')
    paths[2].write_text(json.dumps(outside))
    submitted = []
    def scheduler(command, **kwargs):
        if command[0] == 'squeue':
            return SimpleNamespace(stdout='77\n')
        submitted.append(kwargs['cwd'])
        return SimpleNamespace(stdout='88', stderr='', returncode=0)
    def stop_after_tick(_):
        raise RuntimeError('end of first polling cycle')
    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    monkeypatch.setattr(recon.time, 'sleep', stop_after_tick)
    with pytest.raises(RuntimeError, match='end of first polling cycle'):
        recon.run_many(paths[:2], capacity_manifests=paths, max_jobs=2, startup_checked=True)
    assert len(submitted) == 1  # One old allocation plus one new allocation.
    assert json.loads(paths[2].read_text())['jobs'][0]['job_id'] == '77'


@pytest.mark.parametrize('initial_failure', ['validation', 'execution'])
@pytest.mark.parametrize('remaining_result', ['valid', 'invalid', 'missing'])
def test_shared_failure_drains_other_phases_and_collects_all_results(
        tmp_path, monkeypatch, initial_failure, remaining_result):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _configuration(source)
    paths = [prepare(source, tmp_path / f'recon{i}') for i in range(2)]
    for i, path in enumerate(paths):
        manifest = json.loads(path.read_text())
        manifest['jobs'][0].update(job_id=str(i + 1), submission_status='submitted')
        if i == 0:
            manifest['jobs'][0]['missing_since'] = 0.
            if initial_failure == 'execution':
                (Path(manifest['jobs'][0]['directory']) / 'execution.json').write_text(
                    json.dumps({'status': 'failed', 'error': 'Native runtime failed.'}))
        path.write_text(json.dumps(manifest))
    protected = {p: p.read_bytes() for p in tmp_path.rglob('*.inp')}
    tick = [0]
    def scheduler(command, **kwargs):
        assert command[0] == 'squeue', 'No submissions or cancellations after probe failure'
        # Include an unrelated allocation: it must not delay campaign exit.
        return SimpleNamespace(stdout='2\n999\n' if tick[0] == 0 else '999\n')
    def advance(_):
        tick[0] += 1
        assert tick[0] < 8, 'controller waited for unsubmitted work or an unrelated job'
        if tick[0] != 2 or remaining_result == 'missing':
            return
        # Queue exit preceded file publication by one polling cycle.
        job = json.loads(paths[1].read_text())['jobs'][0]
        directory, name = Path(job['directory']), job['name']
        (directory / (name + '.out')).write_text(
            ('No convergence in rhfpr\n' if remaining_result == 'invalid' else '')
            + 'SETTING RECONDONE = 1\nMolpro calculation terminated\n')
        (directory / (name + '.wfu')).write_bytes(b'completed checkpoint')
        (directory / (name + '.molden')).write_text('[Molden Format]\n')
    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    monkeypatch.setattr(recon.time, 'sleep', advance)
    monkeypatch.setattr(recon.time, 'time', lambda: 100. + 20. * tick[0])
    # Simulate restarting at the pilot phase while full-panel jobs still run.
    with pytest.raises(recon.ProbeFailure) as error:
        recon.run_many(paths[:1], capacity_manifests=paths, startup_checked=True)
    first, second = [json.loads(p.read_text()) for p in paths]
    assert first['status'] == 'failed'
    assert first['jobs'][0]['terminal_failure']
    assert all('job_id' not in m['jobs'][1] for m in (first, second))
    assert all(p.read_bytes() == content for p, content in protected.items())
    job = second['jobs'][0]
    if remaining_result == 'valid':
        assert tick[0] == 2
        assert job['status'] == 'complete' and job['wfu_sha256']
        assert 'terminal_failure' not in job
        assert str(paths[1].parent) not in str(error.value)
    else:
        assert tick[0] == 5  # Publication grace retained for the second probe.
        assert job['status'] == 'incomplete' and job['terminal_failure']
        assert 'wfu_sha256' not in job
        assert str(paths[0].parent) in str(error.value)
        assert str(paths[1].parent) in str(error.value)


def test_retry_archives_known_failure_and_preserves_completed_probe(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    path = prepare(source, tmp_path / 'recon')
    _complete_probes(path)
    report = audit(path)
    first, failed = report['jobs'][:2]
    completed_wfu = Path(first['wfu']).read_bytes()
    directory = Path(failed['directory'])
    (directory / (failed['name'] + '.out')).write_text('I/O error. Perhaps full disk?\n')
    (directory / 'execution.json').write_text(json.dumps({'status': 'failed', 'error': 'scratch I/O'}))
    report['jobs'][1].update(job_id='42', submission_status='submitted')
    path.write_text(json.dumps(report))
    monkeypatch.setattr(recon.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout='42\n'))
    with pytest.raises(RuntimeError, match='allocations'):
        recon.retry_failed(path)
    monkeypatch.setattr(recon.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout=''))
    result = recon.retry_failed(path)
    assert result['retry_ready'] == [failed['name']]
    after = json.loads(path.read_text())
    assert after['jobs'][0]['status'] == 'complete'
    assert Path(first['wfu']).read_bytes() == completed_wfu
    assert 'job_id' not in after['jobs'][1]
    assert after['jobs'][1]['input_sha256'] == failed['input_sha256']
    assert after['jobs'][1]['orbital_predecessor'] == failed['orbital_predecessor']
    assert (directory / 'attempts/1/execution.json').is_file()
    assert (directory / 'attempts/1' / (failed['name'] + '.inp')).is_file()
    assert recon.retry_failed(path)['retry_ready'] == []


@pytest.mark.parametrize('long_names', [False, True])
def test_xml_orbitals_support_native_high_angular_momentum_evidence(tmp_path, monkeypatch, long_names):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    data = json.loads(source.read_text())
    data['title'] = 'ch3ss_orbital_reconnaissance_pilot'
    if long_names:
        data['coordinates'][0]['name'] = 'long_sulfur_orientation_checkpoint'
    source.write_text(json.dumps(data))
    path = prepare(source, tmp_path / 'recon')
    job = json.loads(path.read_text())['jobs'][0]
    directory, name = Path(job['directory']), job['name']
    original = (directory / (name + '.inp')).read_bytes()
    monkeypatch.setattr(recon, 'qc_runtime_environment', lambda *args, **kwargs: ({}, {'scratch_directory': str(tmp_path)}))
    monkeypatch.setattr(recon, 'cleanup_qc_runtime', lambda *args: None)
    def native(command, **kwargs):
        deck = (directory / command[-1]).read_text()
        assert 'put,molden' not in deck.lower()
        assert 'put,xml,' in deck.lower()
        assert len(deck.splitlines()[0].split(',', 1)[1]) <= 32
        assert 'ch3ss_orbital_reconnaissance_pilot orbital reconnaissance:' in original.decode()
        assert command[command.index('-o') + 1] == name + '.out'
        (directory / (name + '.out')).write_text(
            'NATURAL ORBITALS (STATE AVERAGED)\n===\n'
            + ''.join(f'{index}.1 2.00000 -1.0\n' for index in range(1, 12))
            + '12.1 1.00000 -0.1\nNatural orbital dump\n'
            'SETTING RECONENERGY = -436.1 HARTREE\n'
            + ''.join(f'SETTING RECONR{index} = -436.1 HARTREE\n'
                      for index in range(1, len(job['root_states']) + 1))
            + 'SETTING RECONDONE = 1\nMolpro calculation terminated\n')
        checkpoint = deck.split('file,2,', 1)[1].splitlines()[0]
        assert len(checkpoint) <= 31 and '/' not in checkpoint
        assert command[command.index('-W') + 1] == str(directory)
        (directory / checkpoint).write_bytes(b'h functions retained in native WFU')
        for suffix in ('.orbitals.xml', '.natural.orbitals.xml'):
            (directory / (name + suffix)).write_text('<molpro><orbital occupation="1"/></molpro>')
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', native)
    recon.execute(path, name)
    audited = audit(path)['jobs'][0]
    assert audited['status'] == 'complete'
    assert audited['occupation_source'] == 'native_state_averaged_table'
    assert len(audited['orbital_xml']) == 2
    assert (directory / (name + '.inp')).read_bytes() == original
    execution = json.loads((directory / 'execution.json').read_text())
    assert execution['input_transform'] == 'native_xml_short_wfu_v3'
    assert execution['input_sha256'] == hashlib.sha256(original).hexdigest()
    (directory / (name + '.natural.orbitals.xml')).write_text('incomplete XML')
    assert audit(path)['jobs'][0]['status'] == 'incomplete'


@pytest.mark.parametrize('checkpoint_form', ['truncated31', 'truncated32', 'short_native'])
def test_completed_long_name_probe_recovers_exact_native_checkpoint(tmp_path, checkpoint_form):
    """Regression for Blodgett job 52907323: all evidence but FILE,2 passed."""
    from kinbot.molpro import molpro_wfu_name
    source = Path(__file__).parents[1] / 'examples/anl/ch3s_anl1_f12/orbital_reconnaissance.json'
    path = prepare(source, tmp_path / 'recon')
    manifest = json.loads(path.read_text())
    job = next(job for job in manifest['jobs'] if job['orientation'])
    directory, name = Path(job['directory']), job['name']
    requested = name + '.wfu'
    native = directory / (name + '.native.inp')
    native.write_text((directory / (name + '.inp')).read_text())
    if checkpoint_form == 'short_native':
        checkpoint = directory / molpro_wfu_name(name)
        native.write_text(native.read_text().replace(requested, checkpoint.name))
    else:
        checkpoint = directory / requested[:int(checkpoint_form[-2:])]
    checkpoint.write_bytes(b'final state averaged native checkpoint')
    occupations = [2.] * 6 + [1.99808, 1.99228, 1.98755, 1.98755,
                              1.51137, 1.51137, 1.4945, .51729]
    energies = [-436.35736908, -436.35736908, -436.22996844,
                -436.04457692, -435.93856107, -435.93856107]
    (directory / (name + '.out')).write_text(
        'NATURAL ORBITALS (STATE AVERAGED)\n'
        + ''.join(f'{index}.1 {value:.5f} -1.0\n' for index, value in enumerate(occupations, 1))
        + 'Natural orbital dump\nSETTING RECONENERGY = -436.35736908 HARTREE\n'
        + ''.join(f'SETTING RECONR{index} = {value:.8f} HARTREE\n'
                  for index, value in enumerate(energies, 1))
        + 'SETTING RECONDONE = 1.00000000\nMolpro calculation terminated\n')
    exports = {}
    for kind in ('scf', 'natural'):
        export = directory / (kind + '.xml')
        export.write_text('<molpro/>')
        exports[kind] = str(export)
    execution = {'status': 'complete', 'returncode': 0,
                 'input_sha256': job['input_sha256'], 'native_input': str(native),
                 'native_input_sha256': hashlib.sha256(native.read_bytes()).hexdigest(),
                 'orbital_exports': exports}
    record = directory / 'execution.json'
    record.write_text(json.dumps(execution))
    # A staged predecessor with the expected long name must not mask the
    # actual final checkpoint produced by Molpro's truncated FILE name.
    (directory / requested).write_bytes(b'staged starting orbitals')
    job['initial_wfu_sha256'] = hashlib.sha256(b'staged starting orbitals').hexdigest()
    path.write_text(json.dumps(manifest))
    audited = next(item for item in audit(path)['jobs'] if item['name'] == name)
    assert audited['status'] == 'complete'
    assert audited['wfu'] == str(checkpoint)
    assert audited['natural_occupations'] == occupations
    assert audited['root_energies_hartree'] == energies
    # Missing final WFU cannot be replaced by the staged predecessor.
    checkpoint.unlink()
    rejected = next(item for item in audit(path)['jobs'] if item['name'] == name)
    assert rejected['status'] == 'incomplete'
    assert any('checkpoint' in error for error in rejected['validation_errors'])


def test_native_checkpoint_recovery_rejects_ambiguous_or_unbound_files(tmp_path):
    from kinbot.molpro_reconnaissance import _checkpoint_evidence
    name = 'avd__cs_orient_f1_bond_a60_r1p776838825'
    native = tmp_path / 'probe.native.inp'
    native.write_text(f'file,2,{name}.wfu\n')
    job = {'name': name, 'input_sha256': 'prepared'}
    execution = {'input_sha256': 'prepared', 'native_input': str(native),
                 'native_input_sha256': hashlib.sha256(native.read_bytes()).hexdigest()}
    for length in (31, 32):
        (tmp_path / name[:length]).write_bytes(b'checkpoint')
    assert 'Ambiguous' in _checkpoint_evidence(tmp_path, job, execution)[1]
    native.write_text('changed native input')
    assert 'SHA-256 mismatch' in _checkpoint_evidence(tmp_path, job, execution)[1]


def test_run_reports_specific_validation_failure(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    manifest = {'status': 'incomplete', 'jobs': [{'name': 'probe', 'job_id': '42',
        'status': 'incomplete', 'directory': str(tmp_path), 'missing_since': 0.,
        'validation_errors': ['Missing or empty FILE,2 checkpoint: probe.wfu.']}]}
    monkeypatch.setattr(recon, '_submit', lambda *args, **kwargs: manifest)
    monkeypatch.setattr(recon.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout=''))
    with pytest.raises(RuntimeError, match='Missing or empty FILE,2 checkpoint'):
        recon._run(tmp_path / 'reconnaissance.json')


@pytest.mark.parametrize('failure', [
    '? You have requested Molden output, but the basis set contains h or higher basis functions\n',
    'Dump information in style XML to probe.orbitals.xml\n'
    '? String overflow inserting 16 into ch3ss_orbital_reconnaissance_pilot orbital reconnaissance\n'
    '? The problem occurs in strings::string_insert\n',
])
def test_export_only_failure_recovers_the_converged_scf_record(tmp_path, monkeypatch, failure):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    path = prepare(source, tmp_path / 'recon')
    manifest = json.loads(path.read_text())
    job = manifest['jobs'][0]
    job.update(job_id='42', submission_status='submitted')
    path.write_text(json.dumps(manifest))
    directory, name = Path(job['directory']), job['name']
    (directory / (name + '.out')).write_text(
        '!RHF STATE 1.1 Energy -437.09556216\n'
        'PROGRAMS * TOTAL RHF-SCF INT FILE\n' + failure)
    (directory / (name + '.wfu')).write_bytes(b'converged RHF canonical record')
    (directory / 'execution.json').write_text(json.dumps({
        'status': 'failed', 'input_sha256': job['input_sha256'], 'error': 'orbital export'}))
    monkeypatch.setattr(recon.subprocess, 'run', lambda *args, **kw: SimpleNamespace(stdout=''))
    assert recon.retry_failed(path)['retry_ready'] == [name]
    saved = json.loads(path.read_text())['jobs'][0]
    assert saved['recovered_scf']['record'] == '2100.2'
    assert saved['recovered_scf']['energy_hartree'] == -437.09556216
    archive = Path(saved['recovered_scf']['path'])
    (directory / (name + '.wfu')).write_bytes(archive.read_bytes())
    monkeypatch.setattr(recon, 'qc_runtime_environment', lambda *args, **kw: ({}, {'scratch_directory': str(tmp_path)}))
    monkeypatch.setattr(recon, 'cleanup_qc_runtime', lambda *args: None)
    def native(command, **kwargs):
        deck = (directory / command[-1]).read_text()
        assert '{rhf;' not in deck
        assert 'orbital,2100.2,type=canonical' in deck
        assert 'reconenergy=-437.09556216' in deck
        assert '{multi,cas-ci;' in deck and 'wf,23,1,1' in deck
        assert len(deck.splitlines()[0].split(',', 1)[1]) <= 32
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', native)
    recon.execute(path, name)


@pytest.mark.parametrize('output', [
    # A generic string overflow is not evidence of an XML-only failure.
    '!RHF STATE 1.1 Energy -437.09556216\nPROGRAMS * TOTAL RHF-SCF INT FILE\n'
    '? String overflow inserting 16\n? The problem occurs in strings::string_insert\n',
    # An XML error without a converged RHF energy cannot recover a checkpoint.
    'PROGRAMS * TOTAL RHF-SCF INT FILE\nDump information in style XML\n'
    '? String overflow inserting 16\n? The problem occurs in strings::string_insert\n',
])
def test_retry_does_not_recover_unproven_export_only_scf(tmp_path, monkeypatch, output):
    from types import SimpleNamespace
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    path = prepare(source, tmp_path / 'recon')
    manifest = json.loads(path.read_text())
    job = manifest['jobs'][0]
    job.update(job_id='42', submission_status='submitted')
    path.write_text(json.dumps(manifest))
    directory, name = Path(job['directory']), job['name']
    (directory / (name + '.out')).write_text(output)
    (directory / (name + '.wfu')).write_bytes(b'unverified checkpoint')
    (directory / 'execution.json').write_text(json.dumps({
        'status': 'failed', 'input_sha256': job['input_sha256']}))
    monkeypatch.setattr(recon.subprocess, 'run', lambda *args, **kw: SimpleNamespace(stdout=''))
    assert recon.retry_failed(path)['retry_ready'] == [name]
    assert 'recovered_scf' not in json.loads(path.read_text())['jobs'][0]


def test_runtime_setup_failure_is_recorded_before_molpro_launch(tmp_path, monkeypatch):
    from kinbot import molpro_reconnaissance as recon
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    path = prepare(source, tmp_path / 'recon')
    job = json.loads(path.read_text())['jobs'][0]
    def unavailable(*args, **kwargs):
        raise RuntimeError('scratch capacity requirement not met')
    monkeypatch.setattr(recon, 'qc_runtime_environment', unavailable)
    with pytest.raises(RuntimeError, match='scratch capacity'):
        recon.execute(path, job['name'])
    record = json.loads((Path(job['directory']) / 'execution.json').read_text())
    assert record['status'] == 'failed'
    assert record['error'] == 'scratch capacity requirement not met'
    assert record['input_sha256'] == job['input_sha256']
