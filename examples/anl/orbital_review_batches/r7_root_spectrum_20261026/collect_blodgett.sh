#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_TOOLS="$PWD/examples/anl/orbital_review_batches/r7_root_spectrum_20261026"
KB_PANEL="$KB_RUN/orbital_recon/root_spectrum_20261026"
KB_REPORT="$KB_PANEL/selection_report.txt"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon monitor "$KB_PANEL" >"$KB_REPORT"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon character "$KB_PANEL" >>"$KB_REPORT"
"$PWD/.venv/bin/python" -B "$KB_TOOLS/extract_results.py" "$KB_RUN" >>"$KB_REPORT"
"$PWD/.venv/bin/python" -B "$KB_TOOLS/capture_status.py" "$KB_RUN" >"$KB_PANEL/production_status_after.json"
KB_EXPORT="$HOME/kinbot-evidence/root-spectrum-$(date -u +%Y%m%dT%H%M%SZ)"
"$PWD/.venv/bin/python" -B examples/anl/tools/export_orbital_evidence.py \
  "$KB_TOOLS/export.json" --repo "$PWD" --runs "$HOME/kinbot-runs" \
  --production "$KB_RUN" --output-directory "$KB_EXPORT"
"$PWD/.venv/bin/python" -B -m kinbot.run_monitor "$KB_RUN" --all --details >"$KB_EXPORT/production_monitor.txt"
printf 'Report: %s\nVerified export directory: %s\n' "$KB_REPORT" "$KB_EXPORT"
