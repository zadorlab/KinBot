#!/usr/bin/env bash
# New diagnostics only. The retired production source remains read-only.
set -euo pipefail
cd "$HOME/KinBot"
test "$(git branch --show-current)" = composite
export PYTHONPATH="$PWD"
KB_SOURCE="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_POOL="$HOME/kinbot-runs/ch3s_ch3ss_orbital_reviews_20261027"
KB_TOOLS="$PWD/examples/anl/orbital_review_batches/r7_closure_20261027"
KB_PANEL="$KB_POOL/closure"
"$PWD/.venv/bin/python" -B - "$KB_POOL" <<'PY'
from pathlib import Path
import sys
from kinbot import run_policy
pool=Path(sys.argv[1]);policy=run_policy.read(pool)
assert policy and 0 < policy['max_jobs'] <= 3, 'Use the existing three-slot diagnostic pool.'
PY
"$PWD/.venv/bin/python" -B examples/anl/orbital_review_batches/r7_state_support_20261026/check_new_requests.py "$KB_SOURCE"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
 batch "$KB_SOURCE" "$KB_TOOLS/screens.json" "$KB_PANEL" \
 --shared-pool --resource-root "$KB_POOL" --preflight --resume
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
 batch "$KB_SOURCE" "$KB_TOOLS/screens.json" "$KB_PANEL" \
 --shared-pool --resource-root "$KB_POOL" --resume
KB_LOG="$KB_POOL/closure.run.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.requested_orbital_recon \
 run "$KB_PANEL" --report-failures >"$KB_LOG" 2>&1 </dev/null &
printf 'Review PID: %s\nLog: %s\nPanel: %s\n' "$!" "$KB_LOG" "$KB_PANEL"
