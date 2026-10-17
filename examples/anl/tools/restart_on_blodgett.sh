#!/usr/bin/env bash
# Optional clean replacement of r4 after the atomic-SCF fix. Run after git pull.
# A pull alone is sufficient for new QC workers in the continuing r4 runs.
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
KB_TOOL=/home/hram/KinBot/examples/anl/tools/run_examples.py
KB_PREVIOUS=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016_r4_dft
KB_REPLACEMENT=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016_r5_dft
export PATH=/home/hram/KinBot/.venv/bin:$PATH
export PYTHONPATH=/home/hram/KinBot${PYTHONPATH:+:$PYTHONPATH}
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
mkdir -p "$KINBOT_MOLPRO_SCRATCH"

# Preparation and all preflight checks precede any cancellation/archive.
# Stops only the two previous PES trees, waits for their jobs/processes to
# disappear, then moves their common parent to a dated .stopped.* sibling.
# Replacement controllers detach on this login host; only QC jobs use Slurm.
"$KB_PY" -u "$KB_TOOL" launch "$KB_REPLACEMENT" --replace-production "$KB_PREVIOUS"
"$KB_PY" "$KB_TOOL" monitor "$KB_REPLACEMENT"
