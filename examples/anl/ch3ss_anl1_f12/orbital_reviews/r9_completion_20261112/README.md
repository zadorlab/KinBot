# R9 orbital completion campaign

The native review and scope decisions are in `REVIEW.md`. This campaign adds
the two new H-loss reviews, sixteen necessary angular-state checks for case 7,
and eight exact endpoint-binding checks for case 13. It reuses the original
hash-bound seeds. The unfinished diagnostic optimizations receive zero retries.

The dependent screening and full-path checks run under one controller. The
controller retains native completion records across repeated invocation and
uses the existing shared resource policy capped at twenty eight-core exclusive
jobs. `completion.json` reports native completion and outstanding trial issues.
Scientific selection remains an explicit review of the resulting evidence.

Local validation passed 171 focused tests and five discovery subtests. A clean
published-source snapshot with this patch passed 37 additional tests, including
complete simulated orchestration, interruption/reentry, completed-work adoption,
failed prerequisites, equivalent-space deduplication and the read-only report.
All 92 initial native inputs were rendered and restaged successfully. These
tests cover execution logic and input construction. The new native outcomes
will come from Blodgett.

## Launch on Blodgett

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
test "$(git branch --show-current)" = composite
git pull --ff-only origin composite
bash examples/anl/ch3ss_anl1_f12/orbital_reviews/r9_completion_20261112/run_blodgett.sh
BASH
```

The command prints the controller PID, log and campaign directory. The source
seed files must match the exact hashes in `plan.json`. An existing campaign
owner prevents another owner from staging over it.

## Read-only text handoff

After the campaign completes, the following detached reader writes a text
report under `~/kinbot-evidence`. It checks native input/output/checkpoint
hashes, prints root-resolved energies and native orbital/CI evidence for the
selected diagnostic paths, and computes same-geometry return subspace overlaps.
The report also includes every screening energy and failed or pending status.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_RECIPE="$PWD/examples/anl/ch3ss_anl1_f12/orbital_reviews/r9_completion_20261112"
KB_PANEL="$HOME/kinbot-runs/ch3s_ch3ss_orbital_reviews_20261027/r9_completion_20261112"
mkdir -p "$HOME/kinbot-evidence"
KB_REPORT=$(mktemp "$HOME/kinbot-evidence/r9-orbital-review.XXXXXXXX.txt")
nohup "$PWD/.venv/bin/python" -B -u "$KB_RECIPE/report.py" "$KB_PANEL" \
  > "$KB_REPORT" 2> "$KB_REPORT.log" < /dev/null &
printf 'Reader PID %s\nText report %s\nReader log %s\n' \
  "$!" "$KB_REPORT" "$KB_REPORT.log"
BASH
```

The final `SUMMARY` line marks the end of the report. Transfer that text file
to local `examples/anl/tmp`. Its native orbital coefficients, CI tables,
root energies, source hashes and overlap diagnostics support the review while
the binary WFU records stay available on Blodgett for model installation.
