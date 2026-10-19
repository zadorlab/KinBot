#!/usr/bin/env bash
# Clean replacement of r5 with SCF 1e-8, CC 1e-7 and tight geometries.
# Run after git pull. Never delete the previous native evidence.
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
KB_TOOL=/home/hram/KinBot/examples/anl/tools/run_examples.py
KB_PREVIOUS=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016_r5_dft
KB_REPLACEMENT=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261017_r6_dft_scf8
export PATH=/home/hram/KinBot/.venv/bin:$PATH
export PYTHONPATH=/home/hram/KinBot${PYTHONPATH:+:$PYTHONPATH}
# Each worker discovers /scratch/<its OS account> on its compute node.
# Do not select/probe compute scratch from this login-node script.
unset KINBOT_MOLPRO_SCRATCH KINBOT_MOLPRO_SCRATCH_FALLBACKS

# Preparation and all preflight checks precede any cancellation/archive.
# Stops only the two previous PES trees, waits for their jobs/processes to
# disappear, then moves their common parent to a dated .stopped.* sibling.
# Replacement controllers detach on this login host; only QC jobs use Slurm.
"$KB_PY" -u "$KB_TOOL" launch "$KB_REPLACEMENT" --replace-production "$KB_PREVIOUS"
"$KB_PY" "$KB_TOOL" monitor "$KB_REPLACEMENT"
