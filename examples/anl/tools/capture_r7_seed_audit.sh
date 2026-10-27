#!/usr/bin/env bash
# Read-only source capture: no QC, drivers, submissions, cancellations or edits.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_EVIDENCE="$HOME/kinbot-evidence/r7-outage.OMot4WR0"
test -f "$KB_EVIDENCE/retirement.json"
"$PWD/.venv/bin/python" -B -m kinbot.run_seed_audit "$KB_RUN" \
  --output-parent "$HOME/kinbot-evidence" --evidence "$KB_EVIDENCE" \
  --task ch3s/470831050000000000002/anl_automatic/states/ts_72ecf8afb235bd9e/composite/tasks/l3_geometry \
  --task ch3s/471030810040000000002/anl_automatic/states/ts_0450141d9a891e9c/composite/tasks/l3_geometry \
  --task ch3s/470831050000000000002/anl_automatic/states/state_684eb53c4637f96a/composite/tasks/ccsdtq_tz
