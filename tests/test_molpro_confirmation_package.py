"""The example transfer recipe must render and resume through shared execution.

Synthetic WFU bytes exercise provenance only; no native result is fabricated.
"""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tarfile
from types import SimpleNamespace

import pytest

from kinbot.execution import file_sha256

PACKAGE = Path(__file__).parents[1] / 'examples/anl/orbital_reviews/molpro_confirmation'


def module(name):
    spec = importlib.util.spec_from_file_location(name, PACKAGE / (name + '.py'))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def build_fixture(tmp_path):
    data = json.loads((PACKAGE / 'geometry_and_seed_catalog.json').read_text())
    for case in data['cases']:
        p = tmp_path / (case['id'] + '.wfu')
        p.write_bytes(b'Synthetic checkpoint for input/provenance tests, not Molpro evidence.')
        case['seed'].update(path=str(p), sha256=file_sha256(p))
    catalog = tmp_path / 'catalog.json'
    catalog.write_text(json.dumps(data))
    return module('prepare_confirmation').build(catalog, tmp_path / 'run')


def test_confirmation_covers_both_bases_and_sparse_reference_without_repeated_swaps(tmp_path):
    path = build_fixture(tmp_path)
    record = json.loads(path.read_text())
    assert len(record['systems']) == 9
    assert sum(s['jobs'] for s in record['systems']) == 39
    for system in record['systems']:
        cfg = json.loads(Path(system['configuration']).read_text())
        manifest = json.loads(Path(system['manifest']).read_text())
        stage = cfg['confirmation_provenance']['stage']
        assert len(manifest['jobs']) == (5 if stage == 'caspt2' else 3)
        if stage == 'mrci':
            assert cfg['bases'][0]['basis'] == 'avtz-f12'
        for i, job in enumerate(manifest['jobs']):
            inp = (Path(job['directory']) / (job['name'] + '.inp')).read_text()
            assert ('type=natural,ovl' in inp) == (i == 0)
            assert (job['initial_record'] is not None) == (i == 0)
            if 'ch3ss_cs' in system['title']:
                assert ('rotate,17.1,19.1,0' in inp) == (i == 0)
            assert 'maxit,999' in inp
            assert 'thrvar=1.d-10' in inp
            assert 'maxiti=2000' in inp
            if stage == 'mrci':
                assert 'mrci-f12,scalf12=1,cabs_singles=1' in inp
            else:
                assert 'cabs_singles=' not in inp
        if stage == 'caspt2':
            assert manifest['jobs'][1]['geometry_angstrom'] == manifest['jobs'][3]['geometry_angstrom']


def test_confirmation_run_preserves_files_and_rejects_changed_inputs(tmp_path, monkeypatch):
    path = build_fixture(tmp_path)
    run = module('run_confirmation')
    monkeypatch.setenv('SBATCH_MEM_PER_NODE', '200000')
    monkeypatch.setenv('SBATCH_MEM_PER_CPU', '32000')
    def execute(paths, **kwargs):
        assert kwargs['max_jobs'] == 9 and kwargs['allocation_retries'] == 2
        assert kwargs['keep_going'] is True
        assert os.environ['SBATCH_MEM_PER_NODE'] == '104000'
        assert 'SBATCH_MEM_PER_CPU' not in os.environ
        return [{'jobs': [{'status': 'complete'}]} for _ in paths]
    monkeypatch.setattr(run.recon, 'run_many', execute)
    before = {p: p.read_bytes() for p in path.parent.rglob('*.inp')}
    result = run.run(path.parent)
    assert result['status'] == 'complete_review_required' and not result['production_approved']
    assert 'controller_pid' not in result
    assert os.environ['SBATCH_MEM_PER_NODE'] == '200000'
    assert os.environ['SBATCH_MEM_PER_CPU'] == '32000'
    assert all(p.read_bytes() == data for p, data in before.items())
    cfg = Path(result['systems'][0]['configuration'])
    cfg.write_text(cfg.read_text() + '\n')
    with pytest.raises(ValueError, match='Immutable'):
        run.run(path.parent)


def test_confirmation_catalog_controls_reach_native_inputs(tmp_path):
    data = json.loads((PACKAGE / 'geometry_and_seed_catalog.json').read_text())
    data['settings']['solver_options']['scf'].update(
        energy_convergence=1.e-7, max_iterations=351)
    data['settings']['reference_policies']['mrci'] = 'state_average'
    data['settings']['mrci_points'] = ['far', 'tilted']
    prepare = module('prepare_confirmation')
    cfg = next(c for c in prepare.configurations(data)
               if c['title'] == 'ch3s_cs_atz_mrci')
    assert cfg['scf_convergence'] == 1.e-7
    assert [p['name'] for p in cfg['anchors']] == ['far', 'tilted']
    # Exercise real input generation, not just adapter dictionary values.
    seed = tmp_path / 'seed.wfu'
    seed.write_bytes(b'Synthetic input-rendering fixture only')
    cfg['initial_records']['atz__far'].update(path=str(seed), sha256=file_sha256(seed))
    path = tmp_path / 'config.json'
    path.write_text(json.dumps(cfg))
    manifest = json.loads(prepare.prepare(path, tmp_path / 'run').read_text())
    job = manifest['jobs'][0]
    text = (Path(job['directory']) / (job['name'] + '.inp')).read_text()
    assert 'maxit=351' in text
    assert 'gthresh,energy=1.d-7' in text
    assert 'refstate,3,1,2,3;' in text


def test_confirmation_export_preserves_inputs_and_hashes_seeds(tmp_path, monkeypatch):
    path = build_fixture(tmp_path)
    run = module('run_confirmation')
    monkeypatch.setattr(run.subprocess, 'run', lambda *a, **k: SimpleNamespace(stdout=''))
    # No native calculation was run: all 39 points must remain incomplete.
    before = {p: p.read_bytes() for p in path.parent.rglob('*.inp')}
    target = tmp_path / 'evidence.tar.gz'
    result = run.collect(path.parent, target)
    assert result['complete_jobs'] == 0 and result['total_jobs'] == 39
    assert result['production_approved'] is False
    assert result['sha256'] == file_sha256(target)
    with tarfile.open(target) as archive:
        index = json.load(archive.extractfile('export_inventory.json'))
        seeds = [f for f in index['files'] if f['archive_path'].startswith('reviewed_seeds/')]
        assert len(seeds) == 3
        import hashlib
        for artifact in index['files']:
            assert hashlib.sha256(archive.extractfile(artifact['archive_path']).read()).hexdigest() == artifact['sha256']
    assert all(p.read_bytes() == original for p, original in before.items())
    with pytest.raises(FileExistsError):
        run.collect(path.parent, target)


def test_confirmation_export_refuses_live_or_unavailable_scheduler(tmp_path, monkeypatch):
    path = build_fixture(tmp_path)
    run = module('run_confirmation')
    record = json.loads(path.read_text())
    manifest_path = Path(record['systems'][0]['manifest'])
    manifest = json.loads(manifest_path.read_text())
    manifest['jobs'][0]['job_id'] = '123'
    manifest_path.write_text(json.dumps(manifest))
    monkeypatch.setattr(run.subprocess, 'run', lambda *a, **k: SimpleNamespace(stdout='123\n'))
    with pytest.raises(RuntimeError, match='remain live'):
        run.collect(path.parent, tmp_path / 'live.tar.gz')
    def unavailable(*a, **k):
        raise subprocess.CalledProcessError(1, a[0])
    monkeypatch.setattr(run.subprocess, 'run', unavailable)
    with pytest.raises(subprocess.CalledProcessError):
        run.collect(path.parent, tmp_path / 'unknown.tar.gz')


def retirement_fixture(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(PACKAGE))
    adapter = module('retire_previous')
    path = build_fixture(tmp_path)
    catalog = path.parent / 'catalog.json'
    root = tmp_path / 'old'
    root.mkdir()
    (root / '.kinbot_run.json').write_text('{}')
    directory = root / '01_pilot' / 'old_job'
    directory.mkdir(parents=True)
    output = directory / 'old_job.out'
    output.write_text('RECONC1 = -100.0 AU\nRetained native fixture\n')
    manifest = directory.parent / 'reconnaissance.json'
    job = {'name': 'old_job', 'directory': str(directory), 'status': 'complete'}
    manifest.write_text(json.dumps({'jobs': [job]}))
    (root / 'campaign.json').write_text(json.dumps({'systems': [
        {'pilot_manifest': str(manifest), 'manifest': str(manifest)}]}))
    events = []
    def inspect(root):
        events.append('inspect')
        return {'jobs': [{'job_id': '123'}], 'processes': []}
    monkeypatch.setattr(adapter.run_cleanup, 'inspect', inspect)
    monkeypatch.setattr(adapter, 'stop_controller', lambda root: events.append('controller'))
    def stop(root, *, archive, timeout):
        assert archive is False and timeout == 300.
        events.append('stop')
        return {'jobs': [], 'status': 'stopped'}
    monkeypatch.setattr(adapter.run_cleanup, 'stop', stop)
    monkeypatch.setattr(adapter.recon, 'audit', lambda p: json.loads(p.read_text()))
    return adapter, root, catalog, manifest, output, events


def test_retire_preserves_native_and_metadata_and_scopes_shutdown(tmp_path, monkeypatch):
    adapter, root, catalog, manifest, output, events = retirement_fixture(tmp_path, monkeypatch)
    before = output.read_bytes()
    report = tmp_path / 'preserved'
    result = adapter.retire(root, report, catalog)
    assert events == ['inspect', 'controller', 'stop']
    assert result['status'] == 'retired_ready_for_confirmation'
    assert result['counts'] == {'complete': 1}
    assert output.read_bytes() == before
    assert (report / 'metadata_before' / manifest.relative_to(root)).read_bytes() == manifest.read_bytes()
    assert len(result['reviewed_seeds_reused']) == 3


def test_retire_refuses_scheduler_failure_before_cancellation(tmp_path, monkeypatch):
    adapter, root, catalog, _, _, events = retirement_fixture(tmp_path, monkeypatch)
    def unavailable(root):
        raise RuntimeError('scheduler unavailable')
    monkeypatch.setattr(adapter.run_cleanup, 'inspect', unavailable)
    with pytest.raises(RuntimeError, match='scheduler unavailable'):
        adapter.retire(root, tmp_path / 'preserved', catalog)
    assert events == []
    assert not (tmp_path / 'preserved').exists()


def test_retire_requires_matching_seed_before_cancellation(tmp_path, monkeypatch):
    adapter, root, catalog, _, _, events = retirement_fixture(tmp_path, monkeypatch)
    data = json.loads(catalog.read_text())
    Path(data['cases'][0]['seed']['path']).write_bytes(b'changed')
    with pytest.raises(ValueError, match='seed missing or changed'):
        adapter.retire(root, tmp_path / 'preserved', catalog)
    assert events == []


def test_retire_flags_duplicate_geometry_and_space_for_reuse_review(tmp_path, monkeypatch):
    adapter, root, catalog, manifest, _, events = retirement_fixture(tmp_path, monkeypatch)
    config = adapter.configurations(json.loads(catalog.read_text()))[0]
    data = json.loads(manifest.read_text())
    data['jobs'][0].update(symbols=config['symbols'], basis=config['bases'][0]['basis'],
        geometry_angstrom=config['anchors'][0]['geometry_angstrom'],
        probe_active_space=config['probe_definition'])
    manifest.write_text(json.dumps(data))
    result = adapter.retire(root, tmp_path / 'preserved', catalog)
    assert result['status'] == 'reuse_review_required'
    assert result['possible_duplicate_points']
    assert events == ['inspect', 'controller', 'stop']
    data['jobs'][0]['geometry_angstrom'] = [[x] for x,y,z in config['anchors'][0]['geometry_angstrom']]
    assert adapter.matching_points(data['jobs'][0], [config]) == []
