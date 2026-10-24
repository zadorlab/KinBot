# Retain the native selection evidence

Run these commands yourself on Blodgett after the standalone reviews have been
stopped and the publication receipt exists. The exporter includes every file in
the six named selection/reference campaigns: native `.wfu`, `.inp`, `.out`,
audits, manifests and retained failed attempts. It also includes the initial
curated checkpoints, the published model/scan plan, and current production
scan-progress/review/correction metadata. It does not launch QC or change a run.

```bash
cd "$HOME/KinBot"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_EXPORT="$HOME/kinbot-evidence/vrc-selection-$(date -u +%Y%m%dT%H%M%SZ)"
"$PWD/.venv/bin/python" -B examples/anl/tools/export_orbital_evidence.py \
  examples/anl/orbital_review_batches/r7_progress_20261024/export.json \
  --repo "$PWD" --runs "$HOME/kinbot-runs" --production "$KB_RUN" \
  --output-directory "$KB_EXPORT"
printf 'Evidence directory: %s\n' "$KB_EXPORT"
```

The exporter verifies the archive bytes against its manifest. Exit status 2
means missing, changed or uncollected evidence; inspect its issue list and repeat
into a new directory once that cause is resolved. A controller can update a live
metadata file during export; completed native campaign files should be immutable.
The archive is a collection of individually verified snapshots, not an atomic
snapshot of the running calculation.

Use your usual SSH transfer from your Mac to copy the three generated files
(`orbital-evidence.tar.gz`, `.sha256`, `.manifest.json`) into a new local directory
under `examples/anl/orbital_review_batches/r7_progress_20261024/native_exports/`.
For example, substitute the printed remote directory and your SSH alias:

```bash
cd /Users/hrishikeshram/KinBot
KB_LOCAL="examples/anl/orbital_review_batches/r7_progress_20261024/native_exports/selection-20261024"
mkdir -p "$KB_LOCAL"
scp 'YOUR_SSH_ALIAS:/home/hram/kinbot-evidence/PRINTED_DIRECTORY/orbital-evidence.*' "$KB_LOCAL/"
.venv/bin/python -B - "$KB_LOCAL/orbital-evidence.tar.gz" <<'PY'
from pathlib import Path
import sys
sys.path.insert(0,str(Path('examples/anl/tools').resolve()))
from export_orbital_evidence import verify
m=verify(Path(sys.argv[1]))
print('Verified native archive:',len(m['files']),'files')
PY
```

Do not describe new remote files as locally archived until this transfer and
verification succeed. The already retained initial and 21 October archives remain
in their existing documentation directories; this export adds later evidence.

After production scans finish, repeat this metadata export and separately retain
the native scan directories named in the correction/progress records. Those
production WFUs can be large and are deliberately not included in this selection
archive. Recheck atom maps, anchor distances versus COM sampling distances,
fragment connectivity, orbital subspaces and root correspondence in both bases,
asymptotic plateaus, and the location of the variational minimum. The saved
selection audit alone does not certify those later results.
