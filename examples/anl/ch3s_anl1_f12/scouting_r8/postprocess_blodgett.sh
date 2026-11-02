#!/usr/bin/env bash
# Rebuild completed R8 MESS handoffs only, retaining original inputs and QC.
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
KB_PY="$PWD/.venv/bin/python"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8"
KB_EXPORT=$("$KB_PY" -B - "$KB_RUN" <<'PY'
from pathlib import Path
import json, shutil, sys, tempfile
from kinbot.execution import atomic_json
from kinbot.run_cleanup import _processes, _jobs
r=Path(sys.argv[1]).resolve(strict=True)
assert not _processes(r), 'A local owner is still active; leave it running.'
assert not _jobs(r), 'Slurm work is still active; leave it running.'
cases=[]
for case in ('ch3s','ch3ss'):
 d=r/case; source=d/'input.json'; par=json.loads(source.read_text())
 assert par.get('me')==1 and not par.get('run_me')
 assert not par.get('L3_calc') and not par.get('rotdpy_run')
 assert not par.get('vrc_tst_scan') and not par.get('vrc_tst_noscan')
 assert par.get('rotdpy_mess_mode')=='off'
 workers=(d/'chemids').read_text().split()
 assert workers, d
 for w in workers:
  assert json.loads((d/w/'kinbot_status.json').read_text())['status']=='complete', w
 cases.append((d,par))
b=Path.home()/'kinbot-evidence';b.mkdir(exist_ok=True)
out=Path(tempfile.mkdtemp(prefix='r8-mess-refresh.',dir=b))
for d,par in cases:
 archive=out/d.name;archive.mkdir()
 for name in ('input.json','mess_scout.json','pes_status.json','pes.log'):
  if (d/name).is_file():shutil.copy2(d/name,archive/name)
 if (d/'me').is_dir():shutil.copytree(d/'me',archive/'previous_me')
 # Input-only rendering; no scheduler template discovery is needed.
 atomic_json(d/'mess_scout.json',dict(par,correct_submerged=1,run_me=0,run_messpf=0,queuing='local'))
print(out)
PY
)
printf 'EVIDENCE %s\n' "$KB_EXPORT"
for KB_CASE in ch3s ch3ss; do
  (
    cd "$KB_RUN/$KB_CASE"
    "$KB_PY" -B -m kinbot.pes mess_scout.json no-kinbot \
      > "$KB_EXPORT/$KB_CASE/postprocess.log" 2>&1
  ) || { printf 'POSTPROCESS_FAILED %s LOG %s/%s/postprocess.log\n' "$KB_CASE" "$KB_EXPORT" "$KB_CASE"; exit 1; }
done
"$KB_PY" -B - "$KB_RUN" <<'PY'
from pathlib import Path
import hashlib,json,sys
for case in ('ch3s','ch3ss'):
 d=Path(sys.argv[1])/case
 assert json.loads((d/'pes_status.json').read_text())['status']=='complete'
 for n in json.loads((d/'me/mess_networks.json').read_text()):
  p=d/'me'/(n['stem']+'.inp')
  print('MESS_READY',case,p,'WELLS',len(n['wells']),'EXITS',len(n['products']),
        'BARRIERS',len(n['barriers']),'SHA256',hashlib.sha256(p.read_bytes()).hexdigest())
  for line in p.read_text().splitlines():
   if 'correct_submerged=1:' in line:print('MODEL_ADJUSTMENT',line)
PY
