#!/usr/bin/env bash
# Independent diagnostic pool. The retired source is read only.
set -euo pipefail
cd "$HOME/KinBot"
test "$(git branch --show-current)" = composite
export PYTHONPATH="$PWD"
KB_SOURCE="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_POOL="$HOME/kinbot-runs/ch3s_ch3ss_orbital_reviews_20261027"
KB_TOOLS="$PWD/examples/anl/orbital_review_batches/r7_state_support_20261026"
KB_PANEL="$KB_POOL/state_support"
"$PWD/.venv/bin/python" -B - "$KB_SOURCE" "$KB_POOL" <<'PY'
from pathlib import Path
import sys
from kinbot import run_policy
source, pool = map(Path, sys.argv[1:])
policy = run_policy.read(source)
assert policy, 'Missing source resource policy'
assert not pool.resolve().is_relative_to(source.resolve())
pool.mkdir(parents=True, exist_ok=True)
if not (pool/run_policy.NAME).exists():
    run_policy.update(pool, max_jobs=3, exclude=policy['excluded_nodes'])
else:
    current = run_policy.read(pool)
    assert current['max_jobs'] <= 3 and set(policy['excluded_nodes']) <= set(current['excluded_nodes'])
PY
KB_PROOF="$KB_POOL/requests.$(date -u +%Y%m%dT%H%M%SZ).json"
"$PWD/.venv/bin/python" -B -m kinbot.orbital_review_inventory "$KB_SOURCE" >"$KB_PROOF"
"$PWD/.venv/bin/python" -B "$KB_TOOLS/check_new_requests.py" "$KB_SOURCE"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_SOURCE" "$KB_TOOLS/screens.json" "$KB_PANEL" \
  --shared-pool --resource-root "$KB_POOL" --preflight --resume
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_SOURCE" "$KB_TOOLS/screens.json" "$KB_PANEL" \
  --shared-pool --resource-root "$KB_POOL" --resume
KB_LOG="$KB_POOL/run.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.requested_orbital_recon \
  run "$KB_PANEL" --report-failures >"$KB_LOG" 2>&1 </dev/null &
printf 'Review PID: %s\nLog: %s\nComplete request inventory: %s\n' "$!" "$KB_LOG" "$KB_PROOF"
