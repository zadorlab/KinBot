#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_SOURCE="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_POOL="$HOME/kinbot-runs/ch3s_ch3ss_orbital_reviews_20261027"
for KB_NAME in state_support additional_requests; do
  if test -f "$KB_POOL/$KB_NAME/review_launch.json"; then
    "$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
      monitor "$KB_POOL/$KB_NAME" >"$KB_POOL/$KB_NAME/selection_report.txt"
    "$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
      character "$KB_POOL/$KB_NAME" >>"$KB_POOL/$KB_NAME/selection_report.txt"
  fi
done
KB_EXPORT="$HOME/kinbot-evidence/detached-orbitals-$(date -u +%Y%m%dT%H%M%SZ)"
"$PWD/.venv/bin/python" -B examples/anl/tools/export_orbital_evidence.py \
 examples/anl/orbital_review_batches/r7_state_support_20261026/export_detached.json \
 --repo "$PWD" --runs "$HOME/kinbot-runs" --production "$KB_SOURCE" \
 --output-directory "$KB_EXPORT"
printf 'Evidence directory: %s\n' "$KB_EXPORT"
