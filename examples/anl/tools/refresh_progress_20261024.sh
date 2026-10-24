#!/usr/bin/env bash
# Example-specific staging; installed recovery/control modules remain general.
set -euo pipefail
if [[ $# != 1 ]]; then
    echo 'Usage: refresh_progress_20261024.sh PRODUCTION_ROOT' >&2
    exit 2
fi
KB_CHECKOUT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
KB_RUN="$(cd "$1" && pwd)"
KB_REVIEW="$(dirname "$KB_RUN")/r7_progress_review_20261024"
KB_PY="$KB_CHECKOUT/.venv/bin/python"
KB_REVISION="$(git -C "$KB_CHECKOUT" rev-parse HEAD)"
export KB_CHECKOUT KB_RUN KB_REVIEW KB_REVISION
export PYTHONPATH="$KB_CHECKOUT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$KB_REVIEW"
exec 8>"$KB_REVIEW/.handoff.lock"
flock -n 8 || { echo 'This progress handoff already has an owner; nothing changed.' >&2; exit 1; }

# Validate all requested identities/decks and stage finite diagnostics before
# changing controller ownership. Existing native screens are verified/reused.
"$KB_PY" -B -u - <<'PY'
import json,os,sys
from pathlib import Path
checkout=Path(os.environ['KB_CHECKOUT']);root=Path(os.environ['KB_RUN']);output=Path(os.environ['KB_REVIEW'])
sys.path.insert(0,str(checkout/'examples/anl/tools'))
import requested_orbital_recon as screens
import vrc_relaxation_recon as relaxation
plan=checkout/'examples/anl/orbital_review_batches/r7_progress_20261024'
done=output/'refresh_complete.json'
if done.exists() and json.loads(done.read_text())['revision']!=os.environ['KB_REVISION']:
    raise RuntimeError('This diagnostic handoff belongs to another checkout revision; inspect its receipt.')
print('RELAXATION PANEL',json.dumps(relaxation.prepare(root,plan/'relaxations.json',output/'relaxations')),flush=True)
if (output/'screens/review_launch.json').exists():
    saved=screens._saved(output/'screens')
    if json.loads((output/'screens/plan.json').read_text())!=json.loads((plan/'screens.json').read_text()):
        raise RuntimeError('Existing screen plan changed; inspect before resuming.')
    print('ORBITAL SCREENS retained',len(saved['systems']),flush=True)
else:
    print('ORBITAL SCREENS',json.dumps(screens.batch(root,plan/'screens.json',output/'screens')),flush=True)
PY

# Reassert the previously requested C+SH3 exclusion through the owning driver.
export KB_DECISIONS_FILE="$KB_CHECKOUT/examples/anl/orbital_review_batches/r7_progress_20261024/decisions.json"
if [[ ! -f "$KB_REVIEW/refresh_complete.json" ]]; then
    bash "$KB_CHECKOUT/examples/anl/tools/refresh_retained_run.sh" "$KB_RUN" 22
    "$KB_PY" -B - <<'PY'
import os,json
from pathlib import Path
from kinbot.execution import atomic_json
root=Path(os.environ['KB_RUN'])
record=json.loads((root/'controller_refresh.json').read_text())
atomic_json(Path(os.environ['KB_REVIEW'])/'refresh_complete.json',dict(revision=os.environ['KB_REVISION'],code=record['code']))
PY
else
    echo 'Completed controller handoff retained; resuming only its diagnostic owners.'
fi
KB_PIN="$("$KB_PY" -B -c 'import json,sys;print(json.load(open(sys.argv[1]))["code"])' "$KB_REVIEW/refresh_complete.json")"
export PYTHONPATH="$KB_PIN"
cd "$KB_REVIEW"
"$KB_PY" -B -c 'import sys;from kinbot.rotdpy_worker_code import validate;validate(sys.argv[1])' "$KB_PIN"

# Share the three reserved review slots by running these finite campaigns in
# sequence. Production continues on its existing admission pool throughout.
echo 'Running three orbital screens (24 fixed-geometry jobs), then eight relaxed-geometry probes.'
"$KB_PY" -B -u "$KB_CHECKOUT/examples/anl/tools/requested_orbital_recon.py" run "$KB_REVIEW/screens" --report-failures
"$KB_PY" -B -u "$KB_CHECKOUT/examples/anl/tools/requested_orbital_recon.py" evidence "$KB_REVIEW/screens" > "$KB_REVIEW/orbital_evidence.txt"
"$KB_PY" -B -u "$KB_CHECKOUT/examples/anl/tools/vrc_relaxation_recon.py" run "$KB_REVIEW/relaxations"
echo "Diagnostic campaigns finished. Review $KB_REVIEW/orbital_evidence.txt and $KB_REVIEW/relaxations/report.json before publishing new surface models."
