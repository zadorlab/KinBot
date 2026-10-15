import json
import os
from pathlib import Path
import tarfile

import pytest

from kinbot import orbital_campaign as campaign
from test_molpro_reconnaissance import _coordinate_configuration


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
