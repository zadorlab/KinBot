#!/usr/bin/env bash
# User-run; preserve every old campaign, including the independent MRCI retry.
set -euo pipefail
KB_PACKAGE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
KB_REPO=$(cd -- "$KB_PACKAGE/../../../.." && pwd)
cd "$KB_REPO"
export PYTHONPATH="$KB_REPO${PYTHONPATH:+:$PYTHONPATH}"
KB_PY=${KB_PY:-"$KB_REPO/.venv/bin/python"}
KB_SOURCE=${1:-/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_recovery}
KB_RUN=${2:-/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_shift01}
"$KB_PY" "$KB_PACKAGE/prepare_recovery.py" "$KB_SOURCE" "$KB_RUN" \
    --plan "$KB_PACKAGE/caspt2_shift01_plan.json" \
    --review "$KB_PACKAGE/native_recovery_review.json"
nohup "$KB_PY" -u "$KB_PACKAGE/run_confirmation.py" run "$KB_RUN" \
    >> "$KB_RUN/controller.log" 2>&1 < /dev/null &
KB_CONTROLLER_PID=$!
sleep 3
if ! kill -0 "$KB_CONTROLLER_PID" 2>/dev/null; then
    tail -n 60 "$KB_RUN/controller.log"
    exit 1
fi
printf 'All previous campaigns preserved. Matched shift check: %s\nController PID: %s\n' "$KB_RUN" "$KB_CONTROLLER_PID"
"$KB_PY" "$KB_PACKAGE/run_confirmation.py" monitor "$KB_RUN" --all --roots
