#!/usr/bin/env bash
# Run from the updated KinBot checkout on Blodgett. Only the panel owner is refreshed.
set -euo pipefail
KB_REPO=$(git rev-parse --show-toplevel)
cd "$KB_REPO"
test "$(git branch --show-current)" = composite
export PYTHONPATH="$KB_REPO"
KB_NATIVE=${1:-"$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114/native"}
KB_PANEL="$KB_REPO/examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114"
"$KB_REPO/.venv/bin/python" -B -m kinbot.run_control enable "$KB_NATIVE"
"$KB_REPO/.venv/bin/python" -B "$KB_PANEL/reconcile.py" prepare "$KB_NATIVE"
"$KB_REPO/.venv/bin/python" -B "$KB_PANEL/resume.py" refresh "$KB_NATIVE"
