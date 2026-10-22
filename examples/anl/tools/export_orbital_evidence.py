#!/usr/bin/env python3
"""Archive an explicit evidence plan without importing or changing KinBot runs.

All source opens are read-only. No scheduler, manifest audit, native program,
controller lock, or live-refresh consumer is called. WFU files are included.
Production metadata is a collection of individual snapshots, not a transaction.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import stat
import tarfile


def digest(stream):
    h = hashlib.sha256()
    for block in iter(lambda: stream.read(1024 * 1024), b''):
        h.update(block)
    return h.hexdigest()


class HashReader:
    def __init__(self, stream):
        self.stream = stream
        self.hash = hashlib.sha256()

    def read(self, size=-1):
        block = self.stream.read(size)
        self.hash.update(block)
        return block


def stamp(s):
    return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)


def collect(plan, roots):
    files, issues, selections = {}, [], []

    def visit(path, base, label):
        member = str(Path(label) / path.relative_to(base))
        try:
            info = path.lstat()
            if stat.S_ISLNK(info.st_mode):
                issues.append(dict(source=str(path), reason='symlink_not_followed',
                                   target=os.readlink(path)))
            elif stat.S_ISDIR(info.st_mode):
                # Never chase a directory symlink into scratch or another run.
                for child in sorted(path.iterdir()):
                    visit(child, base, label)
            elif stat.S_ISREG(info.st_mode):
                files[member] = path
            else:
                issues.append(dict(source=str(path), reason='non_regular_file'))
        except OSError as error:
            issues.append(dict(source=str(path), reason=str(error)))

    for selection in plan['sources']:
        label, pattern = selection['root'], selection['path']
        if label not in roots or Path(pattern).is_absolute() or '..' in Path(pattern).parts:
            raise ValueError(f'Invalid evidence selection: {selection}')
        base = roots[label]
        matches = sorted(base.glob(pattern))
        selections.append(dict(selection, matches=[str(p) for p in matches]))
        if not matches and not selection.get('optional', False):
            issues.append(dict(source=str(base / pattern), reason='required_source_missing'))
        for path in matches:
            # A glob can pass through symlinked parents even if the leaf is regular.
            if any(p.is_symlink() for p in (path, *path.parents) if p != base and base in p.parents):
                issues.append(dict(source=str(path), reason='symlink_not_followed'))
                continue
            visit(path, base, label)
    return files, issues, selections


def verify(path):
    """Verify every archived byte against the manifest; never extract paths."""
    hashes = {}
    with tarfile.open(path, 'r|gz') as archive:
        manifest = None
        for member in archive:
            if not member.isfile() or member.name in hashes:
                raise ValueError(f'Unexpected archive member: {member.name}')
            stream = archive.extractfile(member)
            if member.name == 'MANIFEST.json':
                if manifest is not None:
                    raise ValueError('Repeated manifest')
                manifest = json.load(stream)
            else:
                hashes[member.name] = (member.size, digest(stream))
    if manifest is None:
        raise ValueError('Archive has no manifest')
    expected = {f['member']: (f['bytes'], f['sha256']) for f in manifest['files']}
    if len(expected) != len(manifest['files']) or hashes != expected:
        raise ValueError('Archive content does not match MANIFEST.json')
    return manifest


def export(plan_path, roots, output):
    plan_path = Path(plan_path).expanduser().resolve()
    roots = {k: Path(v).expanduser().resolve() for k, v in roots.items()}
    output = Path(output).expanduser().resolve()
    if output.exists():
        raise FileExistsError(output)
    # Prevent archive recursion even when a selected source does not exist yet.
    if any(output.is_relative_to(base) for base in roots.values()):
        raise ValueError('Write the export outside all source roots (for example /tmp).')
    plan = json.loads(plan_path.read_text())
    files, issues, selections = collect(plan, roots)
    files['export/plan.json'] = plan_path
    files['export/export_orbital_evidence.py'] = Path(__file__).resolve()
    estimated = sum(p.stat().st_size for p in files.values())
    print(f'Collecting {len(files)} files, {estimated / 1024**3:.2f} GiB before compression; '
          f'{len(issues)} collection issue(s).', flush=True)
    report = dict(schema=1, captured_utc=datetime.now(timezone.utc).isoformat(),
                  scope=plan.get('scope'), roots={k: str(v) for k, v in roots.items()},
                  selections=selections, issues=issues, files=[],
                  note='Read-only per-file snapshots; not a cross-file transaction or QC validation.')
    output.parent.mkdir(parents=True, exist_ok=True)
    # Reserve the whole destination as a new directory in main; library calls
    # also refuse replacement. A failed write remains explicitly .partial.
    partial = output.with_name(output.name + '.partial')
    with partial.open('xb') as raw, tarfile.open(fileobj=raw, mode='w:gz', compresslevel=1) as archive:
        for member, path in sorted(files.items()):
            fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
            with os.fdopen(fd, 'rb') as stream:
                before = os.fstat(stream.fileno())
                if not stat.S_ISREG(before.st_mode):
                    raise ValueError(f'Source is no longer a regular file: {path}')
                info = tarfile.TarInfo(member)
                info.size, info.mtime, info.mode = before.st_size, before.st_mtime, 0o600
                reader = HashReader(stream)
                archive.addfile(info, reader)
                try:
                    unchanged = stamp(before) == stamp(os.fstat(stream.fileno())) == stamp(path.lstat())
                except OSError:
                    unchanged = False
                if not unchanged:
                    issues.append(dict(source=str(path), reason='source_changed_during_copy'))
                report['files'].append(dict(member=member, source=str(path), bytes=info.size,
                                            sha256=reader.hash.hexdigest(), source_unchanged=unchanged))
        data = json.dumps(report, indent=2).encode() + b'\n'
        info = tarfile.TarInfo('MANIFEST.json')
        info.size, info.mode = len(data), 0o600
        archive.addfile(info, io.BytesIO(data))
    print('Verifying archive contents on local export disk...', flush=True)
    verify(partial)
    # link() refuses to overwrite a concurrently created destination.
    os.link(partial, output)
    partial.unlink()
    with output.open('rb') as stream:
        checksum = digest(stream)
    output.with_name(output.name + '.sha256').write_text(f'{checksum}  {output.name}\n')
    output.with_name(output.name + '.manifest.json').write_bytes(data)
    print(f'ARCHIVE VERIFIED: {output}\nSHA256: {checksum}', flush=True)
    print(f'COLLECTION ISSUES: {len(issues)} (see .manifest.json; verification is not scientific approval)')
    for issue in issues[:12]:
        print(issue['reason'] + ': ' + issue['source'])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--production', type=Path, required=True)
    parser.add_argument('--output-directory', type=Path, required=True,
                        help='New, nonexistent directory outside the source roots.')
    args = parser.parse_args()
    destination = args.output_directory.expanduser().resolve()
    roots = {k: getattr(args, k).expanduser().resolve() for k in ('repo', 'runs', 'production')}
    if any(destination.is_relative_to(base) for base in roots.values()):
        parser.error('The output directory must be outside the source roots.')
    destination.mkdir(parents=True, exist_ok=False)
    result = export(args.plan, roots, destination / 'orbital-evidence.tar.gz')
    return 2 if result['issues'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
