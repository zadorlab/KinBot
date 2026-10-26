"""Read current requests/endpoints and bounded failure evidence; submit nothing.

Print JSON to stdout. SQLite connections are URI mode=ro with query_only.
No manifest audit, controller refresh, scheduler mutation or pickle loading.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3

from ase.db import connect


def capture(production):
    root = Path(production).expanduser().resolve(strict=True)
    result = dict(captured_utc=datetime.now(timezone.utc).isoformat(),
                  root=str(root), requests=[], files=[], errors=[],
                  note='Read-only per-file snapshots, not an atomic whole-run snapshot.')
    def read_json(path):
        try:
            return json.loads(path.read_text())
        except (OSError, ValueError) as error:
            result['errors'].append(dict(path=str(path), error=str(error)))
            return {}
    def read(path, limit=131072):
        try:
            with path.open('rb') as stream:
                before = path.stat()
                stream.seek(max(0, before.st_size-limit)); content = stream.read(limit)
                after = path.stat()
            result['files'].append(dict(path=str(path.relative_to(root)), bytes=before.st_size,
                tail_bytes=len(content), tail_sha256=hashlib.sha256(content).hexdigest(),
                mtime_ns=after.st_mtime_ns,
                changed=(before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns),
                text=content.decode(errors='replace')))
        except OSError as error:
            result['errors'].append(dict(path=str(path), error=str(error)))
    for case in ('ch3s', 'ch3ss'):
        d = root/case
        for w in sorted(d.glob('[0-9]*')):
            if not w.is_dir(): continue
            pending=[]
            for p in sorted((w/'vrctst/orbital_review').glob('*/request.json')):
                q=read_json(p)
                if q.get('status')=='needs_review': pending.append((p,q))
            if not pending: continue
            network_path=w/'accepted_l2_network.json'
            network=read_json(network_path) if network_path.exists() else {}
            entries=[]
            for path,q in pending:
                entry=dict(path=str(path.relative_to(root)),request=q,geometries=[],
                    network=[r for r in network.get('network',{}).get('reactions',[])
                             if r.get('name')==q.get('reaction')])
                entries.append(entry)
            result['requests'].extend(entries)
            db_path=w/'kinbot.db'; db=None
            try:
                db=connect(str(db_path),use_lock_file=False)
                db.connection=sqlite3.connect(db_path.as_uri()+'?mode=ro',uri=True)
                db.connection.execute('PRAGMA query_only=ON'); db.connection.execute('BEGIN')
                db.change_count=0
                for entry in entries:
                    q=entry['request']
                    for end in ('_IRC_F_prod_high','_IRC_R_prod_high'):
                        for row in db.select(name=q['reaction']+end,sort='-id',limit=1):
                            atoms=row.toatoms()
                            entry['geometries'].append(dict(name=q['reaction']+end,id=row.id,
                                status=row.data.get('status'),symbols=atoms.get_chemical_symbols(),
                                positions=atoms.positions.tolist()))
            except (OSError,ValueError,sqlite3.Error) as error:
                result['errors'].append(dict(path=str(db_path),error=str(error)))
            finally:
                if db is not None and db.connection is not None: db.connection.close()
    patterns=[
        '.kinbot_resource_policy.json', '.kinbot_live_control.json',
        '.kinbot_controls/requests/*.json', '.kinbot_controls/scheduler/*.json',
        '*/pes_exploration.json', '*/.kinbot_input_refresh.json',
        '*/[0-9]*/kinbot.log', '*/rotdPy/*.execution.json',
        '*/rotdPy/*.campaign.json', '*/rotdPy/*.rotdpy.stderr',
        '*/rotdPy/kb_*/rotdpy.log',
        '*/rotdPy/kb_*/Surface_11[345]/jobs/*.err',
        '*/rotdPy/kb_*/Surface_11[345]/jobs/*.stderr',
        '*/rotdPy/kb_*/Surface_11[345]/jobs/*.stdout',
        '*/rotdPy/kb_*/Surface_11[345]/jobs/*.sh',
        '*/rotdPy/kb_*/Surface_11[345]/jobs/.kinbot_slurm_io_retries/**/execution.json',
        '*/rotdPy/kb_*/Surface_11[345]/**/recovery.json',
    ]
    # Native failed searches identified in the 08:30 monitor, kept outside
    # installed code. Read their attempts; never convert them to accepted TSs.
    failures=[
        ('791671560700020000002-sdc19d1516ced03fc','h2_elim_4_5'),
        ('791671560700020000002-sedd8f7729b69b19b','h2_elim_4_5'),
        ('791671580700000000002-s8fc26f5f496eeb4e','r12_insertion_R_1_2_6'),
        ('791871560500020000002-s254020c87a151164','r12_insertion_R_1_2_3'),
        ('791871560500020000002-sb6e34ded5baaf03b','h2_elim_5_6_high'),
    ]
    for parent,suffix in failures:
        for directory in ('','perm/'):
            for extension in ('.log','.out','.err','.stderr','.result.json','.exitcode','.com','_sella.log'):
                patterns.append(f'ch3ss/{parent}/{directory}{parent}_{suffix}{extension}')
    files={p for pattern in patterns for p in root.glob(pattern) if p.is_file() and not p.is_symlink()}
    for p in sorted(files): read(p)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('production',type=Path)
    args=parser.parse_args()
    print(json.dumps(capture(args.production),indent=2))
