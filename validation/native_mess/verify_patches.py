"""Reconstruct the recorded experimental MESS sources in a temporary directory.

No installed MESS source or executable is changed. This checks source identity
and patch application, not a build or numerical validation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile


def verify(archive, report_path):
    directory = Path(__file__).resolve().parent
    manifest = json.loads((directory/'experimental_source_manifest.json').read_text())
    actual = hashlib.sha256(archive.read_bytes()).hexdigest()
    if actual != manifest['source_archive_sha256']:
        raise ValueError('Source archive differs from the recorded native version')
    results = []
    with tempfile.TemporaryDirectory(prefix='kinbot-mess-patches-') as tmp:
        source = Path(tmp)/'source'
        source.mkdir()
        with tarfile.open(archive) as members:
            for member in members:
                if not member.isfile() or not member.name.startswith('mess-source/src/'):
                    continue
                relative = Path(member.name[len('mess-source/src/'):])
                if relative.is_absolute() or '..' in relative.parts:
                    raise ValueError('Unsafe source archive member')
                destination = source/relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(members.extractfile(member).read())
        layers = [
            ('compiler_compatibility', None),
            ('local_accelerate_lp64', 'mess_mp_lp64'),
            ('finite_pressure_experimental', 'mess_finite'),
            ('direct_hybrid_experimental', 'mess_direct_hybrid'),
            ('direct_conservative_experimental', 'mess_direct_conservative'),
            ('flux_cutoff_experimental', 'mess_flux_cutoff'),
            ('flux_labels_diagnostic', 'mess_flux_labels'),
        ]
        for name, build in layers:
            if name == 'flux_labels_diagnostic':
                source = Path(tmp)/'finite'
            patch = directory/(name+'.patch')
            if hashlib.sha256(patch.read_bytes()).hexdigest() != manifest['patch_sha256'][patch.name]:
                raise ValueError(f'Patch identity mismatch: {name}')
            result = subprocess.run(['patch', '-p1', '--batch', '--forward', '-i', str(patch)],
                                    cwd=source, capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(f'{name}: {result.stdout}\n{result.stderr}')
            if build:
                expected = manifest['builds'][build]['source_files']
            else:
                baseline = manifest['builds']['mess_mp_lp64']['source_files']
                expected = {key: baseline[key] for key in ('libmess/mpack.hh', 'libmess/mess.cc')}
            for filename, sha in expected.items():
                if hashlib.sha256((source/filename).read_bytes()).hexdigest() != sha:
                    raise ValueError(f'Reconstructed file differs: {name}/{filename}')
            results.append(dict(patch=name, applied=True, verified_source_files=len(expected),
                                after_hashes_matched=True))
            if name == 'finite_pressure_experimental':
                shutil.copytree(source, Path(tmp)/'finite')
    report = dict(scope='Temporary source reconstruction, patch identity, and exact source hashes; not build or rate validation.',
                  source_archive_sha256=actual, results=results)
    report_path.write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_archive', type=Path)
    parser.add_argument('--report', required=True, type=Path)
    args = parser.parse_args()
    result = verify(args.source_archive, args.report)
    print('Verified source layers:', len(result['results']))
