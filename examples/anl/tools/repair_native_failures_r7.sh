#!/usr/bin/env bash
# r7 handoff: preserve QC, load parser/recovery fixes, retry only the diagnosed MPI attempt.
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 || ( $# -eq 2 && "$2" != --apply ) ]]; then
  echo "Usage: bash repair_native_failures_r7.sh RUN_ROOT [--apply]" >&2
  exit 2
fi
export KINBOT_NATIVE_RECOVERY_ROOT="$1"
export KINBOT_NATIVE_RECOVERY_APPLY="${2:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$HOME/KinBot/.venv/bin/python" - <<'PY'
from pathlib import Path
import json, os
from kinbot import live_control
from kinbot.anl.dispatch import _load, _verify_stage_files, _verify_execution

root = Path(os.environ['KINBOT_NATIVE_RECOVERY_ROOT']).resolve(strict=True)
target = root / ('ch3s/470831050000000000002/vrctst/molpro_scan/'
    '470831050000000000002_hom_sci_1_2__single_m2_r2_vts_pt14_8605fa2c4a1a0ec2/tasks/relaxed_scan')
graph, spec, state = _load(target.parent.parent)
entry = state['tasks']['relaxed_scan']
reason = 'One serial diagnostic retry of r7 job 52918353 MPI abort; unchanged scientific model.'
previous = [json.loads(p.read_text()) for p in (root/'.kinbot_controls/requests').glob('*.json')]
previous = [r for r in previous if r.get('target') == str(target) and r.get('reason') == reason]
if previous:
    print('Serial retry receipt:', [(r['id'], r['status']) for r in previous], flush=True)
    if any(r['status'] == 'blocked' for r in previous):
        raise SystemExit('Prior retry request is blocked; inspect run_control show. No repeat request published.')
elif entry.get('status') not in ('failed', 'unavailable') or str(entry.get('job_id')) != '52918353' or entry.get('attempt') != 2:
    print('Original failed generation has changed; retained current work:', entry.get('status'), entry.get('job_id'), flush=True)
else:
    task = next(t for t in spec['tasks'] if t['id'] == 'relaxed_scan')
    execution = json.loads((target/'execution.json').read_text())
    _verify_stage_files(graph, task, entry)
    _verify_execution(graph, task, entry, execution)
    if (execution.get('status') != 'failed' or str(execution.get('slurm_job_id')) != '52918353'
            or execution.get('failure_phase') != 'execution'
            or 'relaxed_scan_step_0025 exited with status 255' not in execution.get('error', '')):
        raise SystemExit('Native failure changed; no diagnostic retry published.')
    fallback = {'schema': 1, 'resources': {'cores': 1},
                'policy': {'optimization_step_retries': 0}}
    print('Serial diagnostic fallback:', fallback, flush=True)
    if os.environ['KINBOT_NATIVE_RECOVERY_APPLY']:
        result = live_control.request(root, target, reason=reason,
            expected={k: entry.get(k) for k in ('attempt', 'job_id', 'task_sha256', 'status')},
            fallback=fallback)
        print('Queued:', result['id'], '(running allocations cannot be cancelled)', flush=True)
PY
"$HOME/KinBot/.venv/bin/python" -m kinbot.run_refresh "$1" ${2:+--apply}
if [[ "${2:-}" == --apply ]]; then
  "$HOME/KinBot/.venv/bin/python" -m kinbot.run_control reconcile "$1" \
    --select 'ch3s/*/anl_automatic/states/state_e3f670553d74c260/composite' \
    --select 'ch3s/*/anl_automatic/states/state_ef9e765099fc5f69/composite' --apply
  "$HOME/KinBot/.venv/bin/python" -m kinbot.run_control show "$1"
  "$HOME/KinBot/.venv/bin/python" -m kinbot.run_monitor "$1"
fi
