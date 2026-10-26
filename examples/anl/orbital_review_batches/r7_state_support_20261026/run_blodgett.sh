#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/KinBot"
test "$(git branch --show-current)" = composite
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_TOOLS="$PWD/examples/anl/orbital_review_batches/r7_state_support_20261026"
KB_PANEL="$KB_RUN/orbital_recon/state_support_20261026"
KB_PROOF="$HOME/state_support_requests.$(date -u +%Y%m%dT%H%M%SZ).json"
"$PWD/.venv/bin/python" -B "$KB_TOOLS/check_new_requests.py" "$KB_RUN" >"$KB_PROOF"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_RUN" "$KB_TOOLS/screens.json" "$KB_PANEL" --shared-pool --preflight --resume
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_RUN" "$KB_TOOLS/screens.json" "$KB_PANEL" --shared-pool --resume
cp "$KB_PROOF" "$KB_PANEL/new_requests_preflight.json"
"$PWD/.venv/bin/python" -B "$KB_TOOLS/capture_status.py" "$KB_RUN" >"$KB_PANEL/production_status_before.json"
KB_LOG="$KB_PANEL/run.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.requested_orbital_recon \
  run "$KB_PANEL" --report-failures >"$KB_LOG" 2>&1 </dev/null &
printf 'Review PID: %s\nLog: %s\nRead-only status capture: %s\n' "$!" "$KB_LOG" "$KB_PANEL/production_status_before.json"
