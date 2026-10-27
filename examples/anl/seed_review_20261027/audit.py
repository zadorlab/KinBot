"""Read a verified extracted seed audit; write only an evidence report to stdout.

Example-specific interpretation lives beside this script, outside installed code.
Usage: python -B audit.py /path/to/extracted-seed-evidence
"""
from collections import Counter
import json
from pathlib import Path
import re
import sqlite3
import sys

import numpy as np
from ase.io import read
from kinbot.execution import file_sha256
from kinbot.anl.dispatch import _geometry_hash
from kinbot.search_state import contract

ARCHIVE_SHA256='e1d601d46211fc12106df565bc6fb6676dba1ec30c6ec30bed08f3665ad74732'
MANIFEST_SHA256='b8ae51dd99807505ee125b8cb0090157392cfa51d9f2aec563df31f48dbcf9b7'


def audit(base):
    assert file_sha256(base/'manifest.json')==MANIFEST_SHA256
    manifest=json.loads((base/'manifest.json').read_text())
    assert manifest['status']=='captured' and not manifest['issues']
    for item in manifest['files']:
        p=base/item['path']
        assert p.resolve().is_relative_to(base.resolve()) and not p.is_symlink()
        assert file_sha256(p)==item['sha256'],p
    root=base/'run'
    result=dict(schema=1,archive_sha256=ARCHIVE_SHA256,verified_files=len(manifest['files']),
                source_root=manifest['root'],cases={},geometry_candidates=[],mrcc_candidates=[],
                reuse_approved=False,native_restart_verified=False)
    def proof(path):return dict(path=str(path.relative_to(base)),sha256=file_sha256(path))
    for case in ('ch3s','ch3ss'):
        workers=[];total=Counter();components=Counter();db_rows=0
        for path in sorted((root/case).glob('*/.kinbot_search_state.json')):
            data=json.loads(path.read_text());inp=path.parent/(path.parent.name+'.json')
            counts=Counter(str(x['status']) for x in data['reactions'].values());total.update(counts)
            db=path.parent/'kinbot.db'
            # The snapshot is quiescent; source DB and sidecars were hash-verified.
            with sqlite3.connect(db.as_uri()+'?mode=ro',uri=True) as con:
                check=con.execute('PRAGMA quick_check').fetchall();assert check==[('ok',)]
                n=con.execute('SELECT count(*) FROM systems').fetchone()[0];db_rows+=n
            workers.append(dict(worker=path.parent.name,search_statuses=dict(counts),
                contract_matches=contract(inp)==data['contract'],database_rows=n,
                database_check='ok',ledger=proof(path),database=proof(db),
                pending_submissions={k:v['pending_submission'] for k,v in data['reactions'].items()
                                     if v.get('pending_submission')}))
        for path in sorted((root/case).rglob('state.json')):
            components.update(x['status'] for x in json.loads(path.read_text()).get('tasks',{}).values())
        result['cases'][case]=dict(workers=workers,search_statuses=dict(total),database_rows=db_rows,
                                   recorded_component_statuses=dict(components))
    for path in sorted(root.rglob('optimization.traj')):
        if not path.stat().st_size:continue
        log=path.with_suffix('.log')
        rows=[s.split() for s in log.read_text().splitlines() if s.startswith('Sella ')]
        assert rows and all(len(row)==9 for row in rows)
        row=rows[-1];index=int(row[-1]);step=int(row[1]);frames=read(path,':');point=frames[index]
        taskdir=next(p for p in path.parents if (p/'task.json').is_file())
        record=json.loads((taskdir/'task.json').read_text());original=read(taskdir/'geometry.xyz')
        assert np.allclose(frames[0].positions,original.positions,atol=1e-6,rtol=0)
        assert point.get_chemical_symbols()==record['molecule']['symbols']
        assert _geometry_hash(original)==record['geometry_sha256']
        assert abs(point.get_potential_energy()-float(row[3]))<2e-6
        result['geometry_candidates'].append(dict(task=str(taskdir.relative_to(root)),
            evidence={k:proof(p) for k,p in dict(trajectory=path,log=log,task=taskdir/'task.json',
                                               starting_geometry=taskdir/'geometry.xyz').items()},
            evaluations=len(frames),last_logged_step=step,accepted_frame=index,
            energy_ev=point.get_potential_energy(),
            force_max_ev_angstrom=float(np.linalg.norm(point.get_forces(),axis=1).max()),
            target_fmax=record['task']['optimizer']['fmax'],positions=point.positions.tolist(),
            symbols=point.get_chemical_symbols(),
            status='geometry_seed_candidate_recompute_force_and_hessian'))
    for path in sorted((base/'evidence/checkpoints/nodes').glob('*/preservation.json')):
        receipt=json.loads(path.read_text())
        for f in receipt['files']:
            if not f['source'].endswith('/fort.16') or not f['runtime_sources']:continue
            native_input=next(x for x in receipt['files'] if
                Path(x['source']).parent==Path(f['source']).parent and Path(x['source']).name=='MINP')
            entry=next(x for x in manifest['files'] if x['source']==native_input['target'])
            inp=base/entry['path'];out=next(inp.parent.glob('*.out'));text=out.read_text(errors='replace')
            for runtime in f['runtime_sources']:
                rp=root/runtime
                if not rp.is_file():continue
                runtime_data=json.loads(rp.read_text())
                assert Path(f['source']).parent==Path(runtime_data['native_work_directory'])
                result['mrcc_candidates'].append(dict(
                    runtime=proof(rp),native_input=proof(inp),output=proof(out),receipt=proof(path),
                    preserved_checkpoint=f,
                    last_iteration=re.findall(r'Iteration\s+\d+\s+CC energy:[^\n]+',text)[-1:],
                    converged_ccsdt_energy=re.findall(r'Total CCSDT energy\s*\[au\]:\s*([-\d.]+)',text),
                    spin_cases_started=re.findall(r'Spin case\s+(\d+)',text),
                    input_text=inp.read_text(),status='attributed_candidate_native_restart_not_verified'))
    return result


if __name__=='__main__':
    print(json.dumps(audit(Path(sys.argv[1]).resolve()),indent=2,sort_keys=True))
