"""Read-only report and bounded native export of this validation panel."""
from __future__ import annotations
import argparse
from collections import deque
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import tarfile
import tempfile

SUFFIXES = {'.json', '.com', '.log', '.inp', '.out', '.xml', '.xyz', '.traj',
            '.py', '.slurm', '.sh', '.txt', '.err', '.stdout', '.stderr', '.exitcode', '.db', '.molden'}
DIAGNOSTIC = re.compile(
    r'Normal termination|Error termination|SCF Done|EUMP2|E\(B2PLYP\)|'
    r'Frequencies\s+--|Zero-point correction|S\*\*2|'
    r'path calculation complete|maximum number of|convergence|converged|'
    r'T1 diagnostic|D1 diagnostic|RHF-UCC|Final.*energy|Imaginary|'
    r'Number of steps|NET REACTION COORDINATE|Point Number|Step Number', re.I)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):h.update(block)
    return h.hexdigest()


def selected(root):
    """Known panel directories only; never traverse QC checkpoint/scratch stores."""
    root = root.resolve(); files = set()
    folders = [root, root/'remaining_l2', root/'native']
    manifest = root/'remaining_l2/campaign.json'
    if manifest.is_file():
        for case in json.loads(manifest.read_text())['cases']:
            folders.append(root/'remaining_l2'/case['name'])
    else:
        # A partially staged panel has no complete campaign receipt yet.
        folders.extend(p for p in (root/'remaining_l2').glob('*') if p.is_dir() and not p.is_symlink())
    receipt = root/'native/stage.json'
    if receipt.is_file():
        names = list(json.loads(receipt.read_text())['workflows'])
        for filename in ('reconciliation_v1.json', 'reconciliation_repairs_v2.json', 'reconciliation_electronic_v3.json', 'reconciliation_departure_v3.json'):
            decision = root/'native'/filename
            # Only explicitly declared supplemental graphs; never recursively
            # discover attempt archives or checkpoint stores.
            if decision.is_file(): names.extend(json.loads(decision.read_text()).get('curvature',{}))
        for name in names:
            if not isinstance(name,str) or Path(name).name != name:
                raise ValueError('Invalid graph name in validation manifest.')
            graph = root/'native'/name
            folders.append(graph)
            spec = graph/'workflow.json'
            if spec.is_file():
                folders.extend(graph/'tasks'/task['id'] for task in json.loads(spec.read_text())['tasks'])
    for folder in folders:
        if not folder.is_dir() or folder.is_symlink():continue
        for path in folder.iterdir():
            if (path.is_file() and not path.is_symlink() and path.suffix in SUFFIXES
                    and path.resolve().is_relative_to(root)):
                files.add(path)
    return sorted(files)


def report(root, stream):
    root = root.resolve()
    print('VALIDATION', root, file=stream)
    print('Execution evidence only; no automatic scientific approval or MESS update.', file=stream)
    for path in selected(root):
        if path.name in ('progress.json', 'completion.json', 'results.json'):
            print('\nRECORD', path.relative_to(root), file=stream)
            print(path.read_text(), file=stream)
        if path.name == 'execution.json':
            record = json.loads(path.read_text())
            print('TASK', path.relative_to(root), json.dumps({k: record.get(k) for k in
                ('status', 'error', 'details')}), file=stream)
        if path.suffix in ('.log', '.out'):
            print('\nNATIVE', path.relative_to(root), 'SHA256', digest(path), file=stream)
            # Exact diagnostic lines, not a guessed convergence verdict. Full
            # paths, geometries and mode vectors are in the companion archive.
            with path.open(errors='replace') as handle:
                lines = deque(maxlen=40); count = 0
                for number, line in enumerate(handle, 1):
                    if DIAGNOSTIC.search(line):
                        lines.append((number, line.rstrip())); count += 1
                print('DIAGNOSTIC_LINES', count, 'LAST', len(lines), file=stream)
                for number, line in lines:print(number, line, file=stream)


def collect(root, output_parent):
    root = root.resolve(); output_parent.mkdir(parents=True, exist_ok=True)
    out = Path(tempfile.mkdtemp(prefix='r8-validation-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-', dir=output_parent))
    files = selected(root)
    print('Collecting', len(files), 'explicit native/metadata files; checkpoint stores excluded.', flush=True)
    manifest = []; issues = []
    archive = out/'validation-evidence.tar.gz'
    with tempfile.TemporaryDirectory(dir=out) as temporary:
        with tarfile.open(archive, 'w:gz', compresslevel=1) as tar:
            for index, source in enumerate(files):
                relative = str(source.relative_to(root)); staged = Path(temporary)/'snapshot'
                # Consistent read-only database snapshot, even if QC is finishing.
                if source.suffix == '.db':
                    if staged.exists():staged.unlink()
                    with sqlite3.connect(source.as_uri()+'?mode=ro', uri=True) as src:
                        with sqlite3.connect(staged) as dest:src.backup(dest)
                    sha = digest(staged); size = staged.stat().st_size
                    tar.add(staged, arcname=relative)
                else:
                    sha = digest(source); size = source.stat().st_size
                    tar.add(source, arcname=relative)
                    if digest(source) != sha:issues.append(dict(path=relative, issue='changed during capture'))
                manifest.append(dict(path=relative, bytes=size, sha256=sha, sqlite_snapshot=source.suffix=='.db'))
                if index and index % 100 == 0:print('CAPTURED', index, '/', len(files), flush=True)
    # Verify the actual exported bytes rather than equating archive creation
    # with a successful export. Scientific review is a separate subsequent step.
    expected = {entry['path']: entry['sha256'] for entry in manifest}
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            h = hashlib.sha256()
            with tar.extractfile(member) as handle:
                for block in iter(lambda: handle.read(1024*1024), b''):h.update(block)
            if h.hexdigest() != expected.pop(member.name):issues.append(dict(path=member.name, issue='archive hash mismatch'))
    if expected:issues.append(dict(issue='missing archive members', paths=list(expected)))
    (out/'manifest.json').write_text(json.dumps(dict(files=manifest, issues=issues), indent=2)+'\n')
    with (out/'review.txt').open('w') as stream:report(root, stream)
    print(json.dumps(dict(directory=str(out), archive=str(archive), sha256=digest(archive),
        report=str(out/'review.txt'), files=len(manifest), issues=issues, scientific_approval=False), indent=2))
    return out


if __name__ == '__main__':
    import sys
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('report', 'collect'))
    parser.add_argument('root', type=Path)
    parser.add_argument('--output-parent', type=Path, default=Path.home()/'kinbot-evidence')
    args = parser.parse_args()
    if args.action == 'report':report(args.root, sys.stdout)
    else:collect(args.root, args.output_parent)
