# Lean CH3SS R9 review stage

The requested native archive was reviewed on 2026-11-12. The current decisions,
single follow-up campaign and text-report command are in
`../orbital_reviews/r9_completion_20261112/REVIEW.md` and its `README.md`.
The collection instructions below describe the earlier archive request.

This package stages five populations, four product sets and nine connections.
`scope.json` contains the exact identities, including both encounter complexes
and their connecting reverse abstraction saddle.

The shared seven-channel CAS readiness table, reuse inventory description,
settings and submission requirements are in
`../../ch3s_anl1_f12/production_r9/README.md`.

`input.draft.json` carries the existing CH3 + S2 CAS(7,5) and the previously
installed HSCH2S to CH2S + SH CAS(7,6). `prior_reviews/` preserves the existing
candidate records for the remaining outer channels. `new_h_loss_reviews.json`
records the exact CH3SS to CH2SS + H review still needed.

The already-submitted verification cases needed by this lean scope are 7, 8,
13 and 25. `existing_reviews_export.json` selects their fixed, relaxed and
post-relaxation native evidence plus campaign metadata. Collection is read-only
for the calculation directories. The exporter records missing or changing files
and verifies the resulting archive. Calculation completion and scientific
selection are recorded separately in `readiness.json`.

To collect those existing results on Blodgett, run this after updating `composite`.
The detached process writes its log and archive under `~/kinbot-evidence`.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
KB_PY="$PWD/.venv/bin/python"
KB_PLAN="$PWD/examples/anl/ch3ss_anl1_f12/production_r9/existing_reviews_export.json"
KB_EXPORTER="$PWD/examples/anl/ch3s_anl1_f12/tools/export_orbital_evidence.py"
if ! test -f "$KB_EXPORTER"; then
  KB_EXPORTER="$PWD/examples/anl/tools/export_orbital_evidence.py"
fi
mkdir -p "$HOME/kinbot-evidence"
KB_DEST=$(mktemp -d "$HOME/kinbot-evidence/r9-selected-reviews.XXXXXXXX")
nohup "$KB_PY" -B "$KB_EXPORTER" "$KB_PLAN" \
  --repo "$PWD" --runs "$HOME/kinbot-runs" \
  --production "$HOME/kinbot-runs/ch3s_ch3ss_orbital_reviews_20261027/verification" \
  --output-directory "$KB_DEST/native" \
  > "$KB_DEST/collect.log" 2>&1 < /dev/null &
printf 'Collector PID %s\nLog %s\nArchive when verified %s\n' \
  "$!" "$KB_DEST/collect.log" "$KB_DEST/native/orbital-evidence.tar.gz"
BASH
```
