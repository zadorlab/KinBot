#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_SOURCE="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8"
KB_CHECK="$KB_SOURCE/direct_l2_checks_20261113"
KB_PLAN="examples/anl/ch3ss_anl1_f12/recovery/r8_direct_l2_20261113/campaign.json"
test -f "$KB_SOURCE/.kinbot_resource_policy.json"
"$PWD/.venv/bin/python" -B -m kinbot.saddle_recheck stage \
  "$KB_PLAN" "$KB_CHECK" --source-root "$KB_SOURCE"
KB_LOG="$KB_CHECK/controller.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.saddle_recheck run "$KB_CHECK" \
  > "$KB_LOG" 2>&1 < /dev/null &
KB_PID=$!
printf 'CONTROLLER_PID %s\nCAMPAIGN %s\nLOG %s\n' "$KB_PID" "$KB_CHECK" "$KB_LOG"
