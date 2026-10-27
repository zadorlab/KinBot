#!/usr/bin/env bash
# Explicit diagnostic launch only. Never starts production or publishes models.
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_RECIPE="$PWD/examples/anl/ch3ss_anl1_f12/orbital_reviews/verification_20261027"
KB_POOL="$HOME/kinbot-runs/ch3s_ch3ss_orbital_reviews_20261027"
KB_PANEL="$KB_POOL/verification"
"$PWD/.venv/bin/python" -B -m kinbot.run_policy set "$KB_POOL" --max-jobs 23
"$PWD/.venv/bin/python" -B "$KB_RECIPE/panel.py" stage "$KB_PANEL" --pool "$KB_POOL"
KB_LOG="$KB_POOL/verification.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u "$KB_RECIPE/panel.py" run "$KB_PANEL" >"$KB_LOG" 2>&1 </dev/null &
printf 'Diagnostic PID: %s\nLog: %s\nPanel: %s\n' "$!" "$KB_LOG" "$KB_PANEL"
