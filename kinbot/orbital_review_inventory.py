"""Read-only request, model and endpoint evidence for a complete VRC review list."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3

from ase.db import connect
from kinbot.execution import file_sha256


def inventory(root):
    root = Path(root).expanduser().resolve(strict=True)
    result = dict(schema=1, root=str(root), requests=[], issues=[])
    for case in sorted(root.iterdir()):
        if not case.is_dir() or not (case/'input.json').is_file():
            continue
        config = json.loads((case/'input.json').read_text())
        models = config.get('vrc_tst_surface_models', {})
        for path in sorted(case.glob('[0-9]*/vrctst/orbital_review/*/request.json')):
            request = json.loads(path.read_text()); worker = path.parents[3]
            name = request['reaction']
            record = dict(case=case.name, request=request, path=str(path),
                          sha256=file_sha256(path), endpoints=[], network=[])
            model = models.get(request['parent'], {})
            record['installed_channels'] = [c for c in model.get('channels', [])
                if sorted(c['bond']) == sorted(request['bond'])
                and c.get('reaction', name) == name]
            network = worker/'accepted_l2_network.json'
            if network.exists():
                record['network'] = [r for r in json.loads(network.read_text()).get(
                    'network', {}).get('reactions', []) if r.get('name') == name]
                record['network_sha256'] = file_sha256(network)
            database = worker/'kinbot.db'
            # mode=ro and explicit ASE connection avoid creating a database,
            # journals or .lock files in an archived calculation tree.
            if database.is_file():
                db = connect(str(database), use_lock_file=False)
                db.connection = sqlite3.connect(database.resolve().as_uri()+'?mode=ro', uri=True)
                db.change_count = 0
                try:
                    for side in ('F', 'R'):
                        for level in ('', '_high'):
                            job = name+'_IRC_'+side+'_prod'+level
                            rows = list(db.select(name=job, sort='-id', limit=1))
                            if not rows:
                                continue
                            row = rows[0]; atoms = row.toatoms()
                            record['endpoints'].append(dict(name=job, row_id=row.id,
                                status=row.data.get('status'), symbols=atoms.get_chemical_symbols(),
                                positions=atoms.positions.tolist()))
                finally:
                    db.connection.close()
            elif request.get('status') == 'needs_review':
                result['issues'].append(dict(path=str(database), error='Missing endpoint database.'))
            result['requests'].append(record)
    result['counts'] = dict(Counter(r['request'].get('status', 'unknown') for r in result['requests']))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    print(json.dumps(inventory(args.root), indent=2))


if __name__ == '__main__':
    main()
