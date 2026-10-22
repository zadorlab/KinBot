"""Evidence transport must not mutate QC records or silently drop native files."""
import importlib.util
import json
from pathlib import Path
import tarfile

import pytest

SOURCE = Path(__file__).resolve().parents[1] / 'examples/anl/tools/export_orbital_evidence.py'
spec = importlib.util.spec_from_file_location('export_orbital_evidence', SOURCE)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)


@pytest.fixture
def evidence(tmp_path):
    root = tmp_path / 'runs'
    job = root / 'campaign/native/job'
    job.mkdir(parents=True)
    (job / 'probe.out').write_text('Native result\n')
    (job / 'probe.wfu').write_bytes(b'\x00\xfforbital record')
    (job / 'execution.json').write_text('{"status":"complete"}\n')
    plan = tmp_path / 'plan.json'
    plan.write_text(json.dumps(dict(sources=[dict(root='runs', path='campaign')])))
    return plan, {'runs': root}, tmp_path / 'export/pack.tar.gz', job


def test_archive_preserves_all_native_bytes_and_source_records(evidence):
    plan, roots, output, job = evidence
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in job.iterdir()}
    report = helper.export(plan, roots, output)
    assert not report['issues']
    assert helper.verify(output) == report
    with tarfile.open(output) as archive:
        for p, (data, mtime) in before.items():
            name = 'runs/' + str(p.relative_to(roots['runs']))
            assert archive.extractfile(name).read() == data
            assert p.read_bytes() == data and p.stat().st_mtime_ns == mtime
    assert output.with_name(output.name + '.sha256').is_file()
    with pytest.raises(FileExistsError):
        helper.export(plan, roots, output)
    with pytest.raises(ValueError, match='outside'):
        helper.export(plan, roots, roots['runs'] / 'campaign/recursive.tar.gz')


def test_missing_sources_and_symlinks_are_explicit_without_traversal(evidence):
    plan, roots, output, job = evidence
    outside = output.parent.parent / 'outside'
    outside.mkdir()
    (outside / 'unrelated.wfu').write_bytes(b'not in scope')
    (job / 'scratch').symlink_to(outside, target_is_directory=True)
    data = json.loads(plan.read_text())
    data['sources'].extend([dict(root='runs', path='missing'),
                            dict(root='runs', path='optional_missing', optional=True)])
    plan.write_text(json.dumps(data))
    report = helper.export(plan, roots, output)
    assert {v['reason'] for v in report['issues']} == {'symlink_not_followed', 'required_source_missing'}
    assert len(report['issues']) == 2
    assert not any('unrelated' in f['member'] or 'scratch' in f['member'] for f in report['files'])


def test_replaced_source_is_flagged_and_archived_bytes_still_verify(evidence, monkeypatch):
    plan, roots, output, job = evidence
    original = helper.HashReader.read

    # fdopen streams have an integer .name; associate the native source using
    # the observed content, not a file-system path from the stream wrapper.
    def mutate_read(self, size=-1):
        result = original(self, size)
        if result == b'Native result\n':
            replacement = job / 'replacement'
            replacement.write_text('new generation\n')
            replacement.replace(job / 'probe.out')
        return result

    monkeypatch.setattr(helper.HashReader, 'read', mutate_read)
    report = helper.export(plan, roots, output)
    assert report['issues'] == [dict(source=str(job / 'probe.out'), reason='source_changed_during_copy')]
    assert helper.verify(output) == report


def test_unsafe_plan_path_is_refused(evidence):
    plan, roots, output, job = evidence
    plan.write_text(json.dumps(dict(sources=[dict(root='runs', path='../outside')])))
    with pytest.raises(ValueError, match='Invalid evidence selection'):
        helper.export(plan, roots, output)
    assert not output.exists()
