#!/usr/bin/env bash
# Resume the reviewed r5 downstream failures using accepted calculation caches.
# Both runs must be inactive. No jobs are cancelled by this command.
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
KB_TOOL=/home/hram/KinBot/examples/anl/tools/run_examples.py
KB_RUN=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016_r5_dft
export PATH=/home/hram/KinBot/.venv/bin:$PATH
export PYTHONPATH=/home/hram/KinBot${PYTHONPATH:+:$PYTHONPATH}
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
mkdir -p "$KINBOT_MOLPRO_SCRATCH"
"$KB_PY" -u "$KB_TOOL" resume "$KB_RUN"
"$KB_PY" "$KB_TOOL" monitor "$KB_RUN"
