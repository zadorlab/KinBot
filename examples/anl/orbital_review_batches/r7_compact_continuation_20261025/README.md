# Same-space continuation of two unresolved orbital reviews

**Native archive now received:** see the [current audit and handoff](../r7_compact_union_20261025/README.md).

**Completed:** the supplied monitor at 16:20:12 PDT on 25 October reports
14/14 jobs complete, with no recorded retries. The
[completed-panel audit](continuation_audit.md) compares both initialization
routes and explains the remaining scientific questions. The launch commands
below are retained as provenance; rerunning the completed panel is unnecessary.

The [completed-panel audit](../r7_compact_resolution_20261025/resolution_audit.md)
records why the two near solutions still need examination. This finite plan
retains their existing compact active spaces and state averages. It changes
only the initialization route and adds explicit intermediate separations.

| Request | Unchanged model | Source and new inward points, in each basis |
| --- | --- | --- |
| 791671450790040000002 / hom_sci_1_2 | CAS(3,3), one Q and one D, total weights 0.5 each | Completed 30 Å natural orbitals 20–22 → 6 → 3.5 → 2.5 → original 2 Å target; eight jobs |
| 791871710370000000002 / R_Addition_CSm_R_5_1_2 | CAS(7,5), two D roots, weights 0.5 each | Completed 30 Å natural orbitals 18–22 → 6 → 3.5 → original L2 product complex; six jobs |

Distances are the existing mapped interfragment anchor distances, not COM
dividing-surface distances. Internal fragment geometry, atom order and
orientation are preserved. The last outer-exit point is the actual product
complex, not an imposed 2 Å geometry. This does not alter production grids.

The four source WFU hashes are bound to the supplied completed report. Record
2140.2 is used without SCF-index repacking; closure and discard counts are zero.
Each subsequent point consumes the immediately preceding completed checkpoint
in its own basis. The existing executor journals and checks that WFU hash.
No basis projection, root/weight change, orbital freezing, full-valence space,
CASPT2 calculation, or automatic scientific acceptance is introduced.

The existing recon generator uses `DIAB` for propagated checkpoints. That
operation aligns orbitals within the active span; it does not constrain
active/inactive exchange or prove root continuity. This follows the
[Molpro MULTI documentation](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#diabatic_orbitals).
Assess the native AO character, CI vectors and multiplicity-labeled energies
at every point. Agreement with the earlier near solution would make a direct
30 Å SCF-to-near initialization explanation less likely; disagreement would
identify initialization dependence. Neither outcome alone validates a rate.

## Blodgett commands, through the user's SSH session

Pulling loads the small general continuation-interface extension into these
new review workers. Existing production controllers need no restart for this
diagnostic. The new campaign has a three-job cap within the existing shared
production policy (requested total cap 23), with all existing exclusions.
It neither rewrites that policy nor cancels allocations. Preflight submits
nothing; the second command stages the finite campaign; `run` submits it.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
test "$(git branch --show-current)" = composite
git pull --ff-only origin composite
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_PLAN="$PWD/examples/anl/orbital_review_batches/r7_compact_continuation_20261025/screens.json"
KB_PANEL="$KB_RUN/orbital_recon/compact_continuation_20261025"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_RUN" "$KB_PLAN" "$KB_PANEL" --shared-pool --preflight --resume
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_RUN" "$KB_PLAN" "$KB_PANEL" --shared-pool --resume
KB_LOG="$KB_PANEL/run.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.requested_orbital_recon \
  run "$KB_PANEL" --report-failures >"$KB_LOG" 2>&1 </dev/null &
printf 'Review PID: %s\nLog: %s\n' "$!" "$KB_LOG"
BASH
```

The four recent near jobs took about 16–22 seconds in DZ and 155–192 seconds
in TZ, including their native SCF step. Continuation iteration counts may
differ; these timings do not include scheduler waiting and are not a completion
forecast. There is no reason to restart completed parts of the 60-job panel.

Read-only extraction after completion:

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_PANEL="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/orbital_recon/compact_continuation_20261025"
KB_REPORT="$HOME/compact_continuation_character.$(date -u +%Y%m%dT%H%M%SZ).txt"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon monitor "$KB_PANEL" >"$KB_REPORT"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon character "$KB_PANEL" >>"$KB_REPORT"
printf 'Report: %s\n' "$KB_REPORT"
cat "$KB_REPORT"
BASH
```

Retain the full native report
`$HOME/compact_resolution_native.20261025T224531Z.txt` for transfer with the
next evidence export. The locally saved 60-job character report includes WFU
and XML hashes, not the WFU/XML contents themselves.

The supplied `export.json` collects the new reviews and their source campaigns,
including full inputs/outputs, WFU/XML files and retained attempts. Run after
the continuation finishes, then transfer the printed directory through your
usual SSH/SCP setup into a fresh local evidence directory:

```bash
cd "$HOME/KinBot"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_EXPORT="$HOME/kinbot-evidence/compact-continuation-$(date -u +%Y%m%dT%H%M%SZ)"
"$PWD/.venv/bin/python" -B examples/anl/tools/export_orbital_evidence.py \
  examples/anl/orbital_review_batches/r7_compact_continuation_20261025/export.json \
  --repo "$PWD" --runs "$HOME/kinbot-runs" --production "$KB_RUN" \
  --output-directory "$KB_EXPORT"
"$PWD/.venv/bin/python" -B -m kinbot.run_monitor "$KB_RUN" --all --details \
  >"$KB_EXPORT/production_monitor.txt"
printf 'Evidence directory: %s\n' "$KB_EXPORT"
```

The exporter reads source files without changing jobs, verifies archived bytes
against its manifest, and writes `.sha256` and `.manifest.json` alongside the
archive. Collection issues give exit status 2; inspect them before treating the
export as complete. This supplements the earlier verified selection archive.
The separate `production_monitor.txt` is a fresh read-only production snapshot;
it is saved alongside the archive and is not covered by the archive manifest.
Transfer the entire printed directory, including that snapshot and sidecars.

## Local validation

Regression tests exercise the real plan with synthetic checkpoints: whole-panel
preflight, unchanged counts/roots/weights and near geometries, inward dependency
order, explicit source hashes, rejection of broken chains, preservation of prior
files, and repeat staging. A simulated scheduler checks that every successor
waits for its completed predecessor and receives its exact WFU, retains 22 live
production jobs under cap23, releases admissions, and submits no duplicates on
repeat. These are software/orchestration checks, not native Molpro executions.
