# Four natural-orbital initialization checks, 25 October 2026

Eleven of the twelve new channel selections are prepared. The final
C(SH2) + SH product exit remains held. The ten-job panel is complete and does
not need repeating. This panel contains **four jobs only**, with **CAS(9,6),
two doublets and weights 0.5/0.5 unchanged** from the last comparison.

## What the new evidence establishes

The full same-geometry XML comparison supplied by the user shows that the
previous far-canonical-17 initialization did not yield the intended six-orbital
span. New active 17 is mostly previously inactive occupied sigma character.
Against the identified old near natural 17–22 union, the minimum overlap
singular value is **0.121844 DZ / 0.122409 TZ**. Replacing new active 17 with
its inactive SH-p orbital 16 raises this to **0.990473 / 0.990240**. New 18–22
are otherwise nearly identical to old 18–22. Reconstruction error is below
5.6e-15 and orthogonality error below 9.4e-13. The dominant CI strings retain
the added orbital doubly occupied. These results support correcting the seed
identity; they do not establish the correctness of either optimized branch.

The [numeric evidence and hashes](source_evidence.json) are an explicitly
labeled transcription of the supplied comparison. Original source hashes were
checked against the already verified compact-continuation archive. The last
canonical-trial XML itself still awaits local export. Its native completion
is distinct from scientific acceptance.

## Exact initialization

Use the completed **independent near CAS(7,5)** natural checkpoint from
`compact_resolution_20261025/csh2_sh_exit_resolution/`
`csh2_sh_exit_resolution_union/{adz,atz}__exit_l2`, record **2140.2**.
This source has 17 closed orbitals and active 18–22. Its inactive **17 is the
identified SH-p direction**. Explicitly open that occupied orbital, retaining
natural **17,18,19,20,21,22**: 41 total electrons minus 32 inactive electrons
gives CAS(9,6). The deck reads `start,2140.2,type=natural`, `closed,16`, `occ,22`.
The source is not the canonical checkpoint used in the preceding trial.

| Source checkpoint | SHA-256 |
| --- | --- |
| DZ near WFU | `20cb03e26c6500c89dfdb767297d1c93f094cc2c884842e65ba7c8e79c7e5491` |
| TZ near WFU | `6916ea49d341f6a60ccfa703557203a1de1d06c51f1a05395a5acf09e46198f5` |

Each basis has two ordered jobs: the **unchanged L2 product-complex geometry**
(S2–S3 separation 2.4857243464 Å), followed by the **unchanged 30 Å endpoint**
using the first job's completed checkpoint. Fragment atoms are C1/S2/H4/H6 and
S3/H5. Only the latter fragment translates. There is no internal-geometry
optimization in these orbital probes. Source promotion happens at the first
point only; the second consumes the target CAS. Solver, tolerances, roots and
weights remain unchanged; no CASPT2 job or orbital freeze is added.

The shared `orbital_partition.open_highest_closed` option implements the
explicit promotion. It checks integer counts, available closed orbitals and
source/target consistency, and cannot be combined with `close_lowest > 0`.
Checkpoint hashes, native completion and model counts are checked before
staging/submission. This option performs no automatic active-space expansion
or retry. Molpro defines `CLOSED` orbitals as doubly occupied but optimized;
lowering that boundary makes the specified directions active without freezing
the rest. See the [MULTI orbital-subspace documentation](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#defining_the_orbital_subspaces).

After completion, inspect active **and inactive** native orbitals, near/far
state character, CI vectors and the same-geometry overlap with the intended
source union. Successful termination or a lower energy alone is insufficient.
If the direction leaves the active space again, these results document that
failure; this recipe does not expand the CAS or accept it automatically.

## Launch from the user's Blodgett SSH session

Pull the accompanying commit first. These commands prepare and run only this
four-job diagnostic campaign. The shared cap, node exclusions and production
allocations remain in force; up to three diagnostic jobs are allowed, with
the two basis chains also imposing their predecessor dependencies. No PES
controller restart, production publication or cancellation occurs here.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
test "$(git branch --show-current)" = composite
git pull --ff-only origin composite
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_PLAN="$PWD/examples/anl/orbital_review_batches/r7_natural_union_20261025/screens.json"
KB_PANEL="$KB_RUN/orbital_recon/natural_union_20261025"
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

The campaign lock and `--resume` retain existing prepared/running/completed
work if invoked again. Preflight checks the whole plan before staging it.

## Read-only progress and final extraction

```bash
"$HOME/KinBot/.venv/bin/python" -B -m kinbot.requested_orbital_recon monitor \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/orbital_recon/natural_union_20261025"
```

After completion, capture the character summary, full native MULTI sections,
CI vectors, hashes, same-geometry overlap and latest production monitor.
This reads existing calculations and writes only the report. Transfer that text
file rather than truncating the native sections in a terminal paste.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_PANEL="$KB_RUN/orbital_recon/natural_union_20261025"
KB_REPORT="$HOME/natural_union_character.$(date -u +%Y%m%dT%H%M%SZ).txt"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon monitor "$KB_PANEL" >"$KB_REPORT"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon character "$KB_PANEL" >>"$KB_REPORT"
"$PWD/.venv/bin/python" -B examples/anl/orbital_review_batches/r7_natural_union_20261025/extract_results.py \
  "$KB_RUN" >>"$KB_REPORT"
"$PWD/.venv/bin/python" -B -m kinbot.run_monitor "$KB_RUN" >>"$KB_REPORT"
printf 'Report: %s\n' "$KB_REPORT"
BASH
```

Export both the ten preceding jobs and these four with their WFU/XML/native
files and current metadata. Earlier source campaigns are already archived:

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_EXPORT="$HOME/kinbot-evidence/natural-union-$(date -u +%Y%m%dT%H%M%SZ)"
"$PWD/.venv/bin/python" -B examples/anl/tools/export_orbital_evidence.py \
  examples/anl/orbital_review_batches/r7_natural_union_20261025/export.json \
  --repo "$PWD" --runs "$HOME/kinbot-runs" --production "$KB_RUN" \
  --output-directory "$KB_EXPORT"
"$PWD/.venv/bin/python" -B -m kinbot.run_monitor "$KB_RUN" --all --details \
  >"$KB_EXPORT/production_monitor.txt"
printf 'Evidence directory: %s\n' "$KB_EXPORT"
BASH
```

Transfer the verified export directory into `examples/anl/tools/` locally.
Archive verification establishes file integrity, not scientific acceptance.

## Combined production handoff still pending

The [eleven prepared selections](../r7_compact_union_20261025/selections) and
[combined publication plan](../r7_compact_union_20261025/publication_plan.json)
remain held until this last channel is resolved. Candidate-aware publication
now uses the explicitly named model, refuses ambiguity and changed reviewed
hashes, and prevents held plans from publishing partially.

The final handoff will publish all twelve new selections, reconcile the five
previously published channels, and restart/refresh the owning controllers to
load the endpoint-discovery and publication fixes. It must retain live Slurm
allocations, resource policy, operator decisions and RotdPy campaigns. The
endpoint fix affects VRC preparation; **no earlier PES/IRC/L2/ANL/CBH rerun is
indicated** by the audited defect. Actual post-refresh scan dispatch and
sampling adoption still require remote verification.

Local validation: 209 tests passed for reconnaissance, explicit candidate
selection, combined publication, the actual four-job plan and predecessor
transport; 174 passed and one Linux-only test was skipped for controller
refresh, model inheritance, RotdPy adoption and vdW dispatch. These include
simulated staging/submission races, live-work preservation, decision handling
and repeated invocation. No native result or remote refresh is claimed by
these software tests.
