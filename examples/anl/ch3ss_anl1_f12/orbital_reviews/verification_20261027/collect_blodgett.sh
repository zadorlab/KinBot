#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_RECIPE="$PWD/examples/anl/ch3ss_anl1_f12/orbital_reviews/verification_20261027"
KB_PANEL="$HOME/kinbot-runs/ch3s_ch3ss_orbital_reviews_20261027/verification"
mkdir -p "$HOME/kinbot-evidence"
KB_DEST=$(mktemp -d "$HOME/kinbot-evidence/vrc-verification.XXXXXXXX")
"$PWD/.venv/bin/python" -B "$KB_RECIPE/analyze.py" "$KB_PANEL" --report "$KB_DEST/analysis.json.gz"
KB_EXPORTER="$PWD/examples/anl/ch3s_anl1_f12/tools/export_orbital_evidence.py"
if ! test -f "$KB_EXPORTER"; then
 KB_EXPORTER="$PWD/examples/anl/tools/export_orbital_evidence.py"
fi
"$PWD/.venv/bin/python" -B "$KB_EXPORTER" \
 "$KB_RECIPE/export.json" --repo "$PWD" --runs "$HOME/kinbot-runs" \
 --production "$KB_PANEL" --output-directory "$KB_DEST/native"
printf 'Analysis and verified native evidence: %s\n' "$KB_DEST"
