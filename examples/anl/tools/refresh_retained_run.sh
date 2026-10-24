#!/usr/bin/env bash
# Controller-only recovery of an existing prepared run. Never reopen rejected searches.
# Usage: bash refresh_retained_run.sh ROOT MAX_JOBS [OLD_HANDOFF_PID OLD_HANDOFF_LOG]
set -euo pipefail
if [[ $# != 2 && $# != 4 ]]; then
    echo 'Usage: refresh_retained_run.sh ROOT MAX_JOBS [OLD_HANDOFF_PID OLD_HANDOFF_LOG]' >&2
    exit 2
fi
KB_CHECKOUT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$KB_CHECKOUT"
KB_RUN="$(cd "$1" && pwd)"
export KB_RUN KB_CAP="$2" KB_WAIT_PID="${3:-}" KB_WAIT_LOG="${4:-}"
KB_PY="$KB_CHECKOUT/.venv/bin/python"
export PYTHONPATH="$KB_CHECKOUT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec 9>"$KB_RUN/.controller_handoff.lock"
flock -n 9 || { echo 'Another complete handoff is active; nothing changed.' >&2; exit 1; }

"$KB_PY" -B -u - <<'PY'
from pathlib import Path
import os, time
from kinbot.pes_phase import process_identity
root = Path(os.environ['KB_RUN']).resolve(strict=True)
if not (root/'prepared.json').is_file() or int(os.environ['KB_CAP']) < 0:
    raise ValueError('Expected a prepared run and nonnegative job cap.')
if os.environ.get('KB_DECISIONS_FILE'):
    from kinbot.reaction_decisions import load_batch
    load_batch(root, os.environ['KB_DECISIONS_FILE'])
from kinbot.rotdpy_worker_code import preflight
preflight(root)
print('Handoff input preflight passed.', flush=True)
pid = int(os.environ['KB_WAIT_PID'] or 0)
identity = process_identity(pid) if pid else None
if identity:
    try:
        output = os.readlink(f'/proc/{pid}/fd/1')
    except FileNotFoundError:
        output = None
    if output == os.environ['KB_WAIT_LOG']:
        while process_identity(pid) == identity:
            print(f'WAIT: previous handoff {pid} is still active; no processes signalled.', flush=True)
            time.sleep(30)
    else:
        print('Old PID no longer identifies the recorded handoff; no signal sent.', flush=True)
# Also wait for an orphaned reconciliation writer if its parent exited.
while True:
    writers = []
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():
            continue
        try:
            if p.stat().st_uid != os.getuid():
                continue
            args = (p/'cmdline').read_bytes().decode().strip('\0').split('\0')
            if '-m' not in args or '--apply' not in args:
                continue
            tail = args[args.index('-m')+1:]
            if (tail[:2] == ['kinbot.run_control', 'reconcile'] and len(tail) > 2
                    and Path(tail[2]).resolve() == root and process_identity(int(p.name))):
                writers.append(p.name)
        except (FileNotFoundError, ProcessLookupError, PermissionError, UnicodeError):
            continue
    if not writers:
        break
    print('WAIT: existing reconciliation writer(s): '+', '.join(writers), flush=True)
    time.sleep(30)
PY

"$KB_PY" -B -u - <<'PY'
from pathlib import Path
import os, re, subprocess
from kinbot import run_policy
root = Path(os.environ['KB_RUN'])
text = subprocess.run(['scontrol', 'show', 'nodes', '-o'], check=True,
                      capture_output=True, text=True, timeout=30).stdout
rows = [dict(re.findall(r'(\w+)=(\S+)', line)) for line in text.splitlines() if line.strip()]
if not rows or any(not r.get('NodeName') or not r.get('State') for r in rows):
    raise RuntimeError('Cannot verify node states; nothing refreshed.')
bad = sorted({r['NodeName'] for r in rows if
              re.search(r'DOWN|UNKNOWN|DRAIN|FAIL|NOT_RESPOND|NO_RESPOND', r['State'])})
policy = run_policy.update(root, max_jobs=int(os.environ['KB_CAP']), exclude=bad, manage_pending=True)
print('POOL', policy['max_jobs'], 'EXCLUDED', ','.join(policy['excluded_nodes']), flush=True)
PY

"$KB_PY" -B -u -m kinbot.run_refresh "$KB_RUN" --apply --wait-idle 1800
KB_PIN="$("$KB_PY" -B -c 'import json,sys; from pathlib import Path; print(json.loads((Path(sys.argv[1])/"controller_refresh.json").read_text())["code"])' "$KB_RUN")"
export PYTHONPATH="$KB_PIN"
cd "$KB_RUN"  # Keep the checkout cwd from taking precedence over the code snapshot.

# Retained sampling owners cache their scheduler template, but each allocation
# sources this generated setup. Route future workers through the verified pin.
"$KB_PY" -B -u -m kinbot.rotdpy_worker_code "$KB_RUN" "$KB_PIN"

# Optional operator decisions are consumed by the replacement owners. Never
# publish new control semantics to old pinned drivers before replacement.
"$KB_PY" -B -u - <<'PY'
from pathlib import Path
import json, os, time
from kinbot import live_control as control
from kinbot.reaction_decisions import load_batch, read
from kinbot.pes_phase import active_worker
source = os.environ.get('KB_DECISIONS_FILE')
if source:
    root = Path(os.environ['KB_RUN']).resolve()
    decisions = load_batch(root, source)
    requests = []
    for item in decisions:
        worker = (root/item['worker']).resolve(strict=True)
        choices = {item['reaction']:item['decision']}
        if read(worker)['reactions'].get(item['reaction'], {}).get('decision') == item['decision']:
            print('REACTION DECISION RETAINED', item['reaction'], item['decision'], flush=True)
            continue
        pending = [(p,r) for p,r in control.requests(worker, kinds={'network'})
                   if r['target'] == str(worker) and r['options'].get('reaction_decisions') == choices]
        if pending:
            requests.append(pending[0][0])
            continue
        deadline = time.monotonic()+1800
        while not active_worker(worker):
            status = worker.parent/'pes_status.json'
            owner = json.loads(status.read_text()) if status.exists() else {}
            if owner.get('status') == 'failed' or time.monotonic() > deadline:
                raise RuntimeError('Decision needs a live replacement worker: '+str(worker))
            print('WAIT: replacement worker registration', worker.name, flush=True)
            time.sleep(10)
        record = control.request(root, worker, kind='network', reason=item['reason'],
            options={'reaction_decisions':choices}, cancel_running=item['cancel_running'])
        requests.append(root/'.kinbot_controls/requests'/(record['id']+'.json'))
        print('REACTION DECISION QUEUED', record['id'], item['reaction'], item['decision'], flush=True)
    deadline = time.monotonic()+1800
    while requests:
        for path in requests[:]:
            record = json.loads(path.read_text())
            if record['status'] == 'applied':
                print('REACTION DECISION APPLIED', record['id'], flush=True)
                requests.remove(path)
            elif record['status'] not in ('queued','applying'):
                raise RuntimeError('Reaction decision needs attention: '+json.dumps(record))
        if not requests:
            break
        if time.monotonic() > deadline:
            raise RuntimeError('Decision remains owned/queued; inspect run_control show before another refresh.')
        print('WAIT: replacement drivers are applying reaction decisions.', flush=True)
        time.sleep(15)
PY

# Resume independent campaigns before a potentially slow, whole-run result sweep.
"$KB_PY" -B -u - <<'PY'
from pathlib import Path
import json, os
from kinbot import rotdpy_campaign as campaigns
from kinbot.execution import atomic_json
root = Path(os.environ['KB_RUN']).resolve()
rows = campaigns.resume_retained(root)
atomic_json(root/'sampling_resume.json', dict(schema=1, campaigns=rows))
labels = dict(excluded='CHANNEL EXCLUDED', owner_retained='OWNER RETAINED',
              owner_resumed='OWNER RESUMED', complete_retained='COMPLETE RETAINED',
              held='HELD: SCAN NOT READY')
for row in rows:
    print('ROTD', labels[row['status']], Path(row['input']).name,
          json.dumps({k:v for k,v in row.items() if k not in ('input','status')}), flush=True)
PY

check_owners() {
    "$KB_PY" -B -u - <<'PY'
from pathlib import Path
import json, os, time
from kinbot import pes_phase
root = Path(os.environ['KB_RUN']).resolve()
inputs = [Path(p).resolve() for p in json.loads((root/'prepared.json').read_text())['inputs']]
time.sleep(15)
for inp in inputs:
    if not inp.is_relative_to(root):
        raise RuntimeError(f'Input is outside run: {inp}')
    p = inp.parent/'pes_status.json'
    status = json.loads(p.read_text()) if p.exists() else {}
    p = inp.parent/'pes_exploration.json'
    phase = json.loads(p.read_text()) if p.exists() else {}
    alive = pes_phase.owner_alive(phase)
    print('PES HEALTH', inp.parent.name, 'status', status.get('status'), 'owner_alive', alive,
          'phase', phase.get('phase'), 'error', status.get('error'))
    if status.get('status') == 'failed' or (not alive and status.get('status') != 'complete'):
        raise RuntimeError(f'{inp.parent}: replacement owner is not healthy; inspect driver.log.')
PY
}
check_owners

KB_REPORT="$KB_RUN/outage_reconcile.$(date -u +%Y%m%dT%H%M%SZ).$$.json"
echo "Reconciling saved results and automatic retries; report: $KB_REPORT"
"$KB_PY" -B -u -m kinbot.run_control reconcile "$KB_RUN" --apply > "$KB_REPORT"
"$KB_PY" -B -u - "$KB_REPORT" <<'PY'
import json, sys
report = json.load(open(sys.argv[1]))
for row in report['graphs']:
    if row.get('changes') or row['status'] == 'blocked':
        print('RECONCILE', json.dumps(row))
for warning in report.get('warnings', []):
    print('RECONCILE WARNING', warning)
PY
check_owners
"$KB_PY" -B -m kinbot.run_monitor "$KB_RUN"
echo 'Handoff finished. Review any held/blocked tasks in the report and monitor; live jobs were retained.'
