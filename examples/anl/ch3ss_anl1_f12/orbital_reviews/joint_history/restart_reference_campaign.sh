#!/bin/bash
# Run on Blodgett after logging in yourself. Scope: this example campaign only.
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
RUN=/home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015

"$KB_PY" -m kinbot.orbital_campaign stop-controller "$RUN"
"$KB_PY" -m kinbot.run_cleanup "$RUN" --stop --timeout 300

"$KB_PY" - "$RUN" <<'PY'
import json, shutil, sys
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from kinbot.execution import atomic_json, file_sha256
from kinbot.orbital_campaign import load_settings, retry_failed, retry_startup_failures
from kinbot.molpro_solver import apply_solver_options, validate_solver_options
from kinbot.molpro_reconnaissance import validate_input_blocks

root = Path(sys.argv[1])
path = root / 'campaign_settings.json'
settings = json.loads(path.read_text()) if path.exists() else {'schema': 1}
settings.update(max_jobs=20, pilot_only=True, probe_method='casscf', max_refinements=0,
                keep_going=True, allocation_retries=2, startup_retries=2, io_retries=2,
                memory_mb=104000)
settings.setdefault('scratch_root', '/home/hram/.cache/kinbot/molpro')
settings['excluded_nodes'] = sorted(set(settings.get('excluded_nodes', [])) |
    {'node1','node9','node10','node13','node14','node16','node21','node23','node24','node30','node33'})
settings.pop('reference_solver_options', None)
settings['solver_options'] = {
    'scf': {'energy_convergence': 1e-8, 'max_iterations': 2000},
    'casscf': {'energy_convergence': 1e-8, 'gradient_convergence': 1e-5, 'max_iterations': 999},
    'correlation': {'energy_convergence': 1e-8, 'coefficient_convergence': 1e-5,
                    'max_iterations': 2000, 'max_reference_iterations': 2000,
                    'max_reference_vectors': 60}}
validate_solver_options(settings['solver_options'])
campaign = json.loads((root / 'campaign.json').read_text())
paths = sorted({Path(s[k]) for s in campaign['systems'] for k in ('pilot_manifest','manifest')})
for p in paths:
    for j in json.loads(p.read_text())['jobs']:
        inp = Path(j['directory']) / (j['name'] + '.inp')
        if file_sha256(inp) != j['input_sha256']:
            raise RuntimeError(f'Immutable input changed: {inp}')
        validate_input_blocks(apply_solver_options(inp.read_text(), settings['solver_options']))
stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
backup = root / ('restart_metadata.' + stamp)
backup.mkdir()
for p in [root/'campaign.json', path, *paths]:
    if p.exists():
        target = backup / p.relative_to(root)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, target)
atomic_json(path, settings)
load_settings(root)
Path(settings['scratch_root']).mkdir(parents=True, exist_ok=True)
print(json.dumps(retry_startup_failures(root), indent=2))
ids = {str(j['job_id']) for p in paths for j in json.loads(p.read_text())['jobs'] if j.get('job_id')}
if ids:
    print(json.dumps(retry_failed(root, include_interrupted=True, job_ids=ids), indent=2))
for p in paths:
    m = json.loads(p.read_text())
    changed = False
    for j in m['jobs']:
        if j.get('job_id') or j.get('status') == 'complete':
            continue
        if j.get('retry_solver_options') != settings['solver_options']:
            j.setdefault('solver_settings_history', []).append({
                'previous': j.get('retry_solver_options'), 'time_utc': stamp,
                'reason': 'Explicit restart with reviewed numerical thresholds'})
            j['retry_solver_options'] = deepcopy(settings['solver_options'])
            changed = True
        for key in ('max_reference_iterations', 'max_reference_vectors'):
            if 'retry_' + key in j:
                j['retry_' + key] = settings['solver_options']['correlation'][key]
                changed = True
    if changed:
        atomic_json(p, m)
print('Metadata backup:', backup)
print('Effective settings:', json.dumps(load_settings(root), indent=2))
PY

"$KB_PY" -m kinbot.orbital_campaign refresh-exclusions "$RUN"
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$RUN" \
    >> "$RUN/controller.log" 2>&1 < /dev/null &
KB_CONTROLLER_PID=$!
printf '%s\n' "$KB_CONTROLLER_PID" > "$RUN/controller.pid"
sleep 3
if ! kill -0 "$KB_CONTROLLER_PID" 2>/dev/null; then
    tail -n 50 "$RUN/controller.log"
    exit 1
fi
"$KB_PY" -m kinbot.orbital_campaign monitor "$RUN" --details
