#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_RECIPE="$PWD/examples/anl/ch3ss_anl1_f12/orbital_reviews/r9_completion_20261112"
KB_POOL="$HOME/kinbot-runs/ch3s_ch3ss_orbital_reviews_20261027"
KB_PANEL="$KB_POOL/r9_completion_20261112"
"$PWD/.venv/bin/python" -B -m kinbot.run_policy set "$KB_POOL" --max-jobs 20
"$PWD/.venv/bin/python" -B "$KB_RECIPE/campaign.py" stage "$KB_PANEL" --pool "$KB_POOL"
KB_LOG="$KB_POOL/r9_completion.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u "$KB_RECIPE/campaign.py" run "$KB_PANEL" > "$KB_LOG" 2>&1 < /dev/null &
printf 'Diagnostic PID %s\nLog %s\nCampaign %s\n' "$!" "$KB_LOG" "$KB_PANEL"
