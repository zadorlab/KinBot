#!/usr/bin/env bash
# Run by the user on Blodgett. Preserve the failed confirmation in place.
set -euo pipefail
KB_PACKAGE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
KB_REPO=$(cd -- "$KB_PACKAGE/../../../.." && pwd)
cd "$KB_REPO"
export PYTHONPATH="$KB_REPO${PYTHONPATH:+:$PYTHONPATH}"
KB_PY=${KB_PY:-"$KB_REPO/.venv/bin/python"}
KB_OLD=${1:-/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016}
KB_CONFIRM=${2:-/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2}

"$KB_PY" "$KB_PACKAGE/restart_confirmation.py" "$KB_OLD" "$KB_CONFIRM"
nohup "$KB_PY" -u "$KB_PACKAGE/run_confirmation.py" run "$KB_CONFIRM" \
    >> "$KB_CONFIRM/controller.log" 2>&1 < /dev/null &
KB_CONTROLLER_PID=$!
sleep 3
if ! kill -0 "$KB_CONTROLLER_PID" 2>/dev/null; then
    tail -n 60 "$KB_CONFIRM/controller.log"
    exit 1
fi
printf 'Old evidence retained: %s\nReplacement: %s\nController PID: %s\n' \
    "$KB_OLD" "$KB_CONFIRM" "$KB_CONTROLLER_PID"
"$KB_PY" "$KB_PACKAGE/run_confirmation.py" monitor "$KB_CONFIRM"
