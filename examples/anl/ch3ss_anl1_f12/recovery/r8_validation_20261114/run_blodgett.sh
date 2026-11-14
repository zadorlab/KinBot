#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_SOURCE="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8"
KB_PANEL="examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114"
KB_VALIDATION="$KB_SOURCE/validation_20261114"
test -f "$KB_SOURCE/.kinbot_resource_policy.json"
mkdir -p "$KB_VALIDATION"
"$PWD/.venv/bin/python" -B -m kinbot.saddle_recheck stage \
  "$KB_PANEL/remaining_l2.json" "$KB_VALIDATION/remaining_l2" --source-root "$KB_SOURCE"
KB_STAMP=$(date -u +%Y%m%dT%H%M%SZ)
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.saddle_recheck run \
  "$KB_VALIDATION/remaining_l2" > "$KB_VALIDATION/l2.$KB_STAMP.log" 2>&1 < /dev/null &
printf 'L2_CONTROLLER %s\n' "$!"
nohup "$PWD/.venv/bin/python" -B -u "$KB_PANEL/campaign.py" run \
  "$KB_VALIDATION/native" > "$KB_VALIDATION/native.$KB_STAMP.log" 2>&1 < /dev/null &
printf 'NATIVE_CONTROLLER %s\nVALIDATION %s\n' "$!" "$KB_VALIDATION"
printf 'Both controllers use the existing shared node pool. No PES controller refresh or MESS execution.\n'
