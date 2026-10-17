"""Install the six reviewed, basis-specific checkpoints without overwriting.

Run from any directory. No jobs are submitted or cancelled. Models retain
relative paths so the complete directory can be copied to another machine.
"""
import json
from pathlib import Path
import tarfile

from kinbot.execution import file_sha256


def install(root=None):
    root = Path(root or Path(__file__).parent).resolve()
    inventory = json.loads((root / 'orbital_inventory.json').read_text())
    archive = root / 'reviewed_orbitals.tar.gz'
    if file_sha256(archive) != inventory['archive_sha256']:
        raise ValueError('Orbital archive hash differs from the reviewed inventory.')
    expected = {entry['file']: entry for entry in inventory['seeds']}
    with tarfile.open(archive, 'r:gz') as tar:
        members = tar.getmembers()
        if len(members) != len(expected) or {item.name for item in members} != set(expected):
            raise ValueError('Unexpected orbital archive members.')
        for item in members:
            entry = expected[item.name]
            target = root / item.name
            if (not item.isfile() or target.parent != root / 'orbitals'
                    or item.size != entry['bytes']):
                raise ValueError('Unsafe or inconsistent orbital archive member.')
            if target.exists():
                if file_sha256(target) != entry['sha256']:
                    raise ValueError('Existing checkpoint differs: ' + str(target))
                continue
            target.parent.mkdir(exist_ok=True)
            data = tar.extractfile(item).read()
            import hashlib
            if hashlib.sha256(data).hexdigest() != entry['sha256']:
                raise ValueError('Orbital member hash mismatch: ' + item.name)
            with target.open('xb') as stream:
                stream.write(data)
    return {'installed': len(expected), 'archive_sha256': inventory['archive_sha256']}


if __name__ == '__main__':
    print(json.dumps(install(), indent=2))
