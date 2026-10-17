#!/usr/bin/env bash
# User runs this on Blodgett; no SSH or transfer is performed here.
set -euo pipefail
KB_PACKAGE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
KB_REPO=$(cd -- "$KB_PACKAGE/../../../.." && pwd)
cd "$KB_REPO"
export PYTHONPATH="$KB_REPO${PYTHONPATH:+:$PYTHONPATH}"
KB_PY=${KB_PY:-"$KB_REPO/.venv/bin/python"}
KB_OLD=${1:-/home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015}
KB_CONFIRM=${2:-/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016}
KB_PRESERVE="${KB_OLD}.handoff.$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -e "$KB_CONFIRM" ]]; then
    printf 'Run directory already exists: %s\nUse the documented resume command; no jobs stopped.\n' "$KB_CONFIRM" >&2
    exit 1
fi
# Complete preparation and seed/input validation before cancelling old work.
"$KB_PY" "$KB_PACKAGE/prepare_confirmation.py" "$KB_CONFIRM"
"$KB_PY" "$KB_PACKAGE/retire_previous.py" "$KB_OLD" "$KB_PRESERVE" \
    --catalog "$KB_CONFIRM/catalog.json"

nohup "$KB_PY" -u "$KB_PACKAGE/run_confirmation.py" run "$KB_CONFIRM" \
    >> "$KB_CONFIRM/controller.log" 2>&1 < /dev/null &
KB_CONTROLLER_PID=$!
sleep 3
if ! kill -0 "$KB_CONTROLLER_PID" 2>/dev/null; then
    tail -n 60 "$KB_CONFIRM/controller.log"
    exit 1
fi
printf 'Preserved old evidence: %s\nNew confirmation: %s\nController PID: %s\n' \
    "$KB_PRESERVE" "$KB_CONFIRM" "$KB_CONTROLLER_PID"
"$KB_PY" "$KB_PACKAGE/run_confirmation.py" monitor "$KB_CONFIRM"
