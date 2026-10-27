#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_SOURCE="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_POOL="$HOME/kinbot-runs/ch3s_ch3ss_orbital_reviews_20261027"
KB_TOOLS="$PWD/examples/anl/orbital_review_batches/r7_additional_requests_20261027"
KB_PANEL="$KB_POOL/additional_requests"
test -f "$KB_POOL/.kinbot_resource_policy.json"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
 batch "$KB_SOURCE" "$KB_TOOLS/screens.json" "$KB_PANEL" \
 --shared-pool --resource-root "$KB_POOL" --preflight --resume
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
 batch "$KB_SOURCE" "$KB_TOOLS/screens.json" "$KB_PANEL" \
 --shared-pool --resource-root "$KB_POOL" --resume
KB_LOG="$KB_POOL/additional.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.requested_orbital_recon \
 run "$KB_PANEL" --report-failures >"$KB_LOG" 2>&1 </dev/null &
printf 'Additional review PID: %s\nLog: %s\nBoth panels share the same three-job admission cap.\n' "$!" "$KB_LOG"
