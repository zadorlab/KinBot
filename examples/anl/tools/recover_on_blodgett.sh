#!/usr/bin/env bash
# Run after git pull. Repeating this command never duplicates an r7 driver.
set -euo pipefail
printf 'Preparing detached r7 handoff; r6 will be stopped and archived after preflight.\n'
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
KB_TOOL=/home/hram/KinBot/examples/anl/tools/recover_launch.py
KB_RUNS=/home/hram/kinbot-runs
KB_R7=$KB_RUNS/ch3s_ch3ss_caspt2_production_20261017_r7_dft_scf8
export PATH=/home/hram/KinBot/.venv/bin:$PATH
export PYTHONPATH=/home/hram/KinBot${PYTHONPATH:+:$PYTHONPATH}
# Each worker discovers /scratch/<its OS account> on its compute node.
# Do not select/probe compute scratch from this login-node script.
unset KINBOT_MOLPRO_SCRATCH KINBOT_MOLPRO_SCRATCH_FALLBACKS
mkdir -p "$KB_RUNS"
KB_LOG=$(mktemp "$KB_RUNS/r7-restart.XXXXXX.log")
nohup "$KB_PY" -u "$KB_TOOL" "$KB_R7" \
    --retire "$KB_RUNS/ch3s_ch3ss_caspt2_production_20261016_r5_dft" \
    --retire "$KB_RUNS/ch3s_ch3ss_caspt2_production_20261017_r6_dft_scf8" \
    >"$KB_LOG" 2>&1 </dev/null &
KB_PID=$!
printf 'Detached handoff PID: %s\nProgress log: %s\n' "$KB_PID" "$KB_LOG"
printf 'Follow progress: tail -n 35 -F %q\n' "$KB_LOG"
printf 'This confirms only that the handoff was started, not that the PES drivers are running.\n'
