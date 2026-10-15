import json
import os
from pathlib import Path
import tarfile

import pytest

from kinbot import orbital_campaign as campaign
from test_molpro_reconnaissance import _coordinate_configuration, _complete_probes


def test_preflight_and_batch_activate_the_molpro_backend(tmp_path, monkeypatch):
    import subprocess
    import shlex
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    commands = tmp_path / 'commands'
    commands.mkdir()
    for name, body in {
        'molpro': 'printf "%s" "$KINBOT_BACKEND" > native-backend\n',
        'sbatch': 'exit 0\n', 'squeue': 'exit 0\n',
    }.items():
        command = commands / name
        command.write_text('#!/bin/bash\n' + body)
        command.chmod(0o700)
    setup = ('case "${KINBOT_BACKEND:-}" in\n'
             f'molpro) export PATH={shlex.quote(str(commands))}:"$PATH" ;;\nesac\n')
    monkeypatch.setattr(campaign.recon, 'render_site_setup', lambda _: setup)
    monkeypatch.setenv('PATH', str(commands) + ':' + os.environ['PATH'])
    monkeypatch.setenv('KINBOT_MOLPRO_SCRATCH', str(tmp_path / 'scratch'))
    root = tmp_path / 'campaign'
    saved = json.loads(campaign.prepare(root, [source]).read_text())
    result = campaign.preflight(root)
    assert result['status'] == 'ready_for_native_pilots'
    assert all(item['molpro'] == str(commands / 'molpro') for item in result['setups'])
    manifest = json.loads(Path(saved['systems'][0]['pilot_manifest']).read_text())
    job = manifest['jobs'][0]
    subprocess.run(['bash', job['name'] + '.sbatch'], cwd=job['directory'], check=True)
    assert (Path(job['directory']) / 'native-backend').read_text() == 'molpro'
    (commands / 'molpro').unlink()
    with pytest.raises(RuntimeError, match='site setup failed'):
        campaign.preflight(root)


def test_campaign_snapshots_context_and_prepares_a_native_pilot(tmp_path):
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    context = tmp_path / 'context.json'
    context.write_text('{}')
    data = json.loads(source.read_text())
    data['kinbot_input'] = context.name
    source.write_text(json.dumps(data))
    root = tmp_path / 'campaign'
    saved = json.loads(campaign.prepare(root, [source]).read_text())
    record = saved['systems'][0]
    snap = json.loads(Path(record['configuration']).read_text())
    assert Path(snap['kinbot_input']).is_relative_to(root)
    pilot = json.loads(Path(record['pilot_manifest']).read_text())
    assert len(pilot['jobs']) == 1
    assert pilot['jobs'][0]['root_states']
    assert (root / '.kinbot_run.json').is_file()
    assert campaign.progress(root)['systems'][0]['manifest']['complete'] == 0
    with pytest.raises(FileExistsError):
        campaign.prepare(root, [source])


def test_campaign_controller_has_one_global_cap_and_is_run_scoped(tmp_path, monkeypatch):
    monkeypatch.setattr(campaign, 'preflight', lambda root: None)
    configurations = []
    for number in range(2):
        source = tmp_path / f'input{number}.json'
        _coordinate_configuration(source)
        configurations.append(source)
    root = tmp_path / 'campaign'
    campaign.prepare(root, configurations)
    saved = json.loads((root / 'campaign.json').read_text())
    saved['error'] = 'previous interruption'
    (root / 'campaign.json').write_text(json.dumps(saved))
    calls = []
    original = Path.cwd()
    def run(path, *, max_jobs):
        assert Path.cwd() == root
        assert 'error' not in json.loads((root / 'campaign.json').read_text())
        calls.append((Path(path).parent.name, max_jobs))
        return {'status': 'complete_review_required', 'jobs': [{'status': 'complete'}]}
    monkeypatch.setattr(campaign.recon, 'run', run)
    result = campaign.run(root, max_jobs=10)
    assert calls == [('01_pilot', 10), ('02_pilot', 10), ('01_recon', 10), ('02_recon', 10)]
    assert result['status'] == 'complete_review_required'
    assert Path.cwd() == original
    with campaign._lock(root):
        with pytest.raises(RuntimeError, match='Another controller'):
            campaign.run(root)


def test_failed_pilot_stops_before_full_jobs_and_records_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(campaign, 'preflight', lambda root: None)
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    root = tmp_path / 'campaign'
    campaign.prepare(root, [source])
    calls = []
    def fail(path, **kwargs):
        calls.append(path)
        raise RuntimeError('native pilot failed')
    monkeypatch.setattr(campaign.recon, 'run', fail)
    with pytest.raises(RuntimeError, match='native pilot failed'):
        campaign.run(root)
    assert len(calls) == 1
    assert json.loads((root / 'campaign.json').read_text())['status'] == 'failed'


def test_campaign_rejects_changed_context_before_submission(tmp_path, monkeypatch):
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    root = tmp_path / 'campaign'
    saved = json.loads(campaign.prepare(root, [source]).read_text())
    Path(saved['systems'][0]['configuration']).write_text('{}')
    monkeypatch.setattr(campaign.recon, 'run', lambda *args, **kw: pytest.fail('changed input submitted'))
    with pytest.raises(ValueError, match='configuration changed'):
        campaign.run(root)


def test_excluded_directories_are_archived_idempotently_and_collected(tmp_path, monkeypatch):
    from types import SimpleNamespace
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    config = json.loads(source.read_text())
    config['bases'].append({'name': 'large', 'basis': 'avqz-f12'})
    source.write_text(json.dumps(config))
    root = tmp_path / 'campaign'
    campaign.prepare(root, [source])
    monkeypatch.setattr(campaign.subprocess, 'run', lambda *a, **k:
                        SimpleNamespace(stdout='', returncode=0))
    campaign.exclude_bases(root, ['avqz-f12'])
    saved = json.loads((root / 'campaign.json').read_text())
    manifest_path = Path(saved['systems'][0]['manifest'])
    manifest = json.loads(manifest_path.read_text())
    excluded = manifest['excluded_jobs'][0]
    evidence = Path(excluded['directory']) / (excluded['name'] + '.out')
    evidence.write_text('native diagnostic to preserve')
    active = Path(manifest['jobs'][0]['directory'])
    result = campaign.archive_excluded(root)
    assert result['files_deleted'] == 0
    assert active.is_dir() and not evidence.parent.exists()
    assert campaign.archive_excluded(root) == result
    moved = json.loads(manifest_path.read_text())['excluded_jobs'][0]
    assert moved['directory'] == excluded['directory']
    assert (Path(moved['archived_directory']) / evidence.name).read_text() == 'native diagnostic to preserve'
    archive = campaign.collect(root, tmp_path / 'review.tar.gz')
    with tarfile.open(archive) as bundle:
        assert any(name.endswith(evidence.name) and name.startswith('archived_excluded/')
                   for name in bundle.getnames())


def test_archive_checks_live_excluded_jobs_before_moving_anything(tmp_path, monkeypatch):
    from types import SimpleNamespace
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    config = json.loads(source.read_text())
    config['bases'].append({'name': 'large', 'basis': 'avqz-f12'})
    source.write_text(json.dumps(config))
    root = tmp_path / 'campaign'
    saved = json.loads(campaign.prepare(root, [source]).read_text())
    monkeypatch.setattr(campaign.subprocess, 'run', lambda *a, **k:
                        SimpleNamespace(stdout='', returncode=0))
    campaign.exclude_bases(root, ['avqz-f12'])
    path = Path(saved['systems'][0]['manifest'])
    manifest = json.loads(path.read_text())
    manifest['excluded_jobs'][0]['job_id'] = '99'
    path.write_text(json.dumps(manifest))
    monkeypatch.setattr(campaign.subprocess, 'run', lambda *a, **k:
                        SimpleNamespace(stdout='99\n', returncode=0))
    with pytest.raises(RuntimeError, match='allocations'):
        campaign.archive_excluded(root)
    assert not (root / 'archived_excluded').exists()


@pytest.mark.parametrize('execution_status', [None, 'running', 'complete'])
def test_explicit_interrupted_retry_preserves_evidence_and_rejects_parser_only_failure(
        tmp_path, monkeypatch, execution_status):
    from types import SimpleNamespace
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    root = tmp_path / 'campaign'
    saved = json.loads(campaign.prepare(root, [source]).read_text())
    path = Path(saved['systems'][0]['manifest'])
    manifest = json.loads(path.read_text())
    job = manifest['jobs'][0]
    job.update(job_id='99', submission_status='submitted')
    path.write_text(json.dumps(manifest))
    directory = Path(job['directory'])
    (directory / (job['name'] + '.out')).write_text('partial native result')
    if execution_status:
        (directory / 'execution.json').write_text(json.dumps({'status': execution_status}))
    monkeypatch.setattr(campaign.subprocess, 'run', lambda *a, **k:
                        SimpleNamespace(stdout='', returncode=0))
    campaign.retry_failed(root)
    assert not (directory / 'attempts').exists()
    campaign.retry_failed(root, include_interrupted=True)
    updated = json.loads(path.read_text())['jobs'][0]
    if execution_status == 'complete':
        assert updated['job_id'] == '99'
        assert not (directory / 'attempts').exists()
    else:
        assert 'job_id' not in updated
        assert updated['submission_status'] == 'retry_ready'
        assert (directory / 'attempts/1' / (job['name'] + '.out')).read_text() == 'partial native result'


def test_collection_keeps_binary_wavefunctions_on_cluster(tmp_path):
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    root = tmp_path / 'campaign'
    saved = json.loads(campaign.prepare(root, [source]).read_text())
    manifest = json.loads(Path(saved['systems'][0]['manifest']).read_text())
    first = manifest['jobs'][0]
    directory = Path(first['directory'])
    (directory / (first['name'] + '.wfu')).write_bytes(b'large binary')
    (directory / (first['name'] + '.out')).write_text('incomplete native output')
    archive = Path(campaign.collect(root, tmp_path / 'review.tar.gz'))
    with tarfile.open(archive) as stream:
        names = stream.getnames()
        assert 'campaign.json' in names
        assert any(name.endswith('.out') for name in names)
        assert not any(name.endswith('.wfu') for name in names)


def test_stopped_campaign_basis_selection_preserves_results_and_survives_refinement(tmp_path, monkeypatch):
    from types import SimpleNamespace
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    config = json.loads(source.read_text())
    config['bases'].append({'name': 'large', 'basis': 'av5z-f12'})
    source.write_text(json.dumps(config))
    root = tmp_path / 'campaign'
    prepared = json.loads(campaign.prepare(root, [source]).read_text())
    record = prepared['systems'][0]
    for key in ('pilot_manifest', 'manifest'):
        _complete_probes(Path(record[key]))
    snapshots = {path: path.read_bytes() for path in root.rglob('*') if path.is_file()
                 and path.suffix in ('.inp', '.wfu', '.out')}
    fingerprint = Path(record['configuration']).read_bytes()
    monkeypatch.setattr(campaign.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout=''))
    selected = campaign.exclude_bases(root, ['av5z-f12'])
    assert [entry['retained'] for entry in selected] == [1, 4]
    assert all(path.read_bytes() == contents for path, contents in snapshots.items())
    assert Path(record['configuration']).read_bytes() == fingerprint
    progress = campaign.progress(root)['systems'][0]
    assert progress['pilot_manifest']['complete'] == 1
    assert progress['manifest']['excluded_jobs'] == 4
    assert all(entry['newly_excluded'] == 0 for entry in campaign.exclude_bases(root, ['av5z-f12']))
    campaign.recon.refine(record['manifest'])
    refined = json.loads(Path(record['manifest']).read_text())
    assert len(refined['jobs']) > 4
    assert {job['basis'] for job in refined['jobs']} == {'aug-cc-pVTZ'}
    assert len(refined['excluded_jobs']) == 4
    archive = campaign.collect(root, tmp_path / 'review.tar.gz')
    with tarfile.open(archive) as stream:
        assert any('/large__' in name and name.endswith('.out') for name in stream.getnames())
        assert not any(name.endswith('.wfu') for name in stream.getnames())


def test_basis_selection_checks_all_live_allocations_before_mutation(tmp_path, monkeypatch):
    from types import SimpleNamespace
    source = tmp_path / 'input.json'
    _coordinate_configuration(source)
    data = json.loads(source.read_text())
    data['bases'].append({'name': 'large', 'basis': 'av5z-f12'})
    source.write_text(json.dumps(data))
    root = tmp_path / 'campaign'
    saved = json.loads(campaign.prepare(root, [source]).read_text())
    last = Path(saved['systems'][0]['manifest'])
    manifest = json.loads(last.read_text())
    manifest['jobs'][-1].update(job_id='77', submission_status='submitted')
    last.write_text(json.dumps(manifest))
    original = {path: path.read_bytes() for path in root.rglob('reconnaissance.json')}
    monkeypatch.setattr(campaign.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout='77\n'))
    with pytest.raises(RuntimeError, match='all campaign allocations'):
        campaign.exclude_bases(root, ['av5z-f12'])
    assert all(path.read_bytes() == contents for path, contents in original.items())
    monkeypatch.setattr(campaign.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(stdout=''))
    with pytest.raises(ValueError, match='Unknown reconnaissance'):
        campaign.exclude_bases(root, ['not_a_basis'])
    with pytest.raises(ValueError, match='manifest empty'):
        campaign.exclude_bases(root, ['aug-cc-pVTZ', 'av5z-f12'])
