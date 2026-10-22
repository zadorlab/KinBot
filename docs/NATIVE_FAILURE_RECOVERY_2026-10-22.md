# Native failure recovery, r7, 2026-10-22

The carbon job requested CCSDTQ(P), but MRCC printed converged CCSDTQ and
normal termination. Its native counts are six total electrons, two frozen,
three correlated alpha and one correlated beta. Quintuple excitations cannot
exist in that four-electron space. The saved -37.765267479485 Eh result is
therefore valid for the requested frozen-core CCSDTQ(P) component; it needs
reparsing, not another QC calculation. The CCSDTQ minus CCSDT(Q) increment
is still evaluated from the two actual energies, not set to zero.

## General rank handling

Both CCSDT(Q) and CCSDTQ(P) now accept native lower-rank iterative results
only with unambiguous total/core/alpha/beta counts, a sufficient excitation
rank, explicit convergence, normal termination and the required native
settings. Missing or contradictory evidence remains a failure. Conventional
CCSD(T) is not substituted for iterative CCSDT with three electrons.

New tasks can declare the total molecular electron count: one maps to HF,
two to CCSD, three to CCSDT, and four maps CCSDTQ(P) to CCSDTQ while leaving
CCSDT(Q) intact. Total count is a conservative bound; the planner does not
guess how MRCC selects a frozen core. Reductions due to frozen electrons
require native evidence. Bare HF acceptance additionally requires one total
electron and zero frozen electrons. Requested/effective methods and count
evidence remain in result provenance and component assembly. The existing
one/two-electron ANL shortcut still omits all exactly zero higher-order
increment jobs.

## Gaussian VPT2

The CH2 triplet retry failed in Link 717 with `New: C2V - Old: CS` after the
earlier empty-group failure had selected `PG=CS`. Recovery now permits the
new native group, with a persistent history, at most two distinct named-group
attempts, and no repeat/cycle. The native old group must match the previous
automatic assignment. Method, geometry and all other input options are
retained; each failed attempt is archived. These rules apply to every Gaussian
VPT2 dispatcher task. They do not classify a generic segmentation message as
a symmetry failure.

## Molpro diagnostic retry and handoff

Job 52918353 aborted at numerical-force evaluation 25 after rank 5 called
MPI_Abort. Other ranks' signal-9 messages do not establish OOM or a bad node.
Slurm accounting is disabled. The handoff queues one serial retry of that
exact failed generation, with the original 100 optimization steps and no step
extensions, same initial seed, method, root, active space and thresholds.
Its existing partition-maximum walltime is retained. Attempt 3 is beyond the
automatic two-attempt transient budget; a repeat failure remains for review.
This is a diagnostic, not a proven fix for the native abort or the observed
force irregularity. New accepted/running generations are left alone.

Run on Blodgett:

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
bash examples/anl/tools/repair_native_failures_r7.sh \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7" --apply
BASH
```

Omitting `--apply` previews the retry and refresh without submitting or
cancelling jobs (the refresh may create an immutable code snapshot). The
serial request is queued before the controller replacement so the new scan
owner sees it on replay. A changed generation is retained; request ownership
and artifact hashes are checked again by the consuming driver. Repeating the
handoff does not publish another request for the diagnosed attempt.

The refresh replaces only controllers to load the new Python code. It retains
running/allocated QC, completed outputs, the existing direct CASPT2-F12 CBS
inputs and the installed CH + H2S model. No new scientific migration flags or
orbital proposal are needed. For future runs with this version, the MRCC and
VPT2 recoveries are automatic; explicit operator retries use `run_control`.

The monitor should show the carbon component recovered, the CH2 VPT2 retry
staged/submitted/complete, and the serial VRC request applied (or current work
explicitly retained). `owned_by_controller` during reconciliation means its
live driver owns recovery. Staging alone is not proof of native success.
The C + SH3 scan connectivity hold remains and must be reviewed before rates.

## Follow-up: reviewed capacity, duplicate graphs, and archived abort evidence

The effective cap of 20 was the configured 23 minus three orbital-review
slots. Installing the CH + H2S model did not acknowledge the old request until
reaction replay reached VRC. Workers now validate installed, already reviewed
channels before replay, and at supported live model-update boundaries. They
check parent identity, fragment terms, all required basis checkpoints and
state/root contracts. Spin-splitting models additionally require their native
reference review. A valid installed model changes the orbital request to
`model_accepted`; the shared policy then releases the reservation automatically.
Invalid models retain the hold. This does not accept a scan or waive its
connectivity review. The monitor displays configured capacity, the reservation,
and the requests responsible when any slots remain reserved.

Deduplication alone now preserves L1 reuse-only settings and existing bounded
search-retry receipts. It also preserves and uses configured comparison
tolerances instead of resetting them to defaults; inconsistent tolerances
within one PES stop the preflight. Exact graph reuse and outer-PES lowest-L2
pathway selection remain separate. Raw reaction discoveries remain in history.
Distinct saddles or stereochemical pathways are not made equivalent simply
because their products match.

The one serial diagnostic is now attempt 3, Slurm job 52919565. It started
at 12:29:50 with the original scientific model. The initial force is still
large (2.0869 eV/angstrom), so serialization alone does not establish smooth
state-specific numerical derivatives. The saved attempt-2 files contain no
more specific native error than the rank-5 MPI abort. The initiating cause is
unresolved; this update does not claim a native Molpro repair. It retains the
running diagnostic and never queues another copy.

The next supplied read shows serial steps 0 through 4, ending at 12:47:24,
without a failure receipt; convergence is still unproven. That read also shows
capacity 23, no outstanding orbital-review requests and zero reserved slots:
the old worker eventually acknowledged CH + H2S. The new acknowledgement code
prevents the delay on future updates; it did not cause this already observed
remote transition.

New failure receipts retain error contexts from a bounded 16 MiB window per
file, in addition to the existing 64 KiB tails. Aggregate context reads are
bounded, and both dotted and underscore-numbered native logs are recognized.
Truncation is recorded. A read-only report searches archived attempts and
separates their failed job IDs from the current retry:

```bash
python -m kinbot.native_failure_report /path/to/graph/tasks/relaxed_scan \
  --label relaxed_scan_step_0025
```

For r7, the following full controller refresh retains VRC work (including the
serial diagnostic), completed QC and current direct CASPT2-F12 CBS models.
Only allocations belonging to verified duplicate ANL graphs can be cancelled;
their files are archived. The preflight is repeated with controllers quiesced
before changing ownership. No rejected-channel revisit is requested.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_PY="$HOME/KinBot/.venv/bin/python"
"$KB_PY" -u -m kinbot.run_refresh "$KB_RUN" --deduplicate-anl --apply
"$KB_PY" -m kinbot.run_monitor "$KB_RUN"
BASH
```

The handoff prints verified duplicate allocation IDs and retained IDs. Once
the CH + H2S worker acknowledges its model, expect capacity 23 unless another
current channel needs orbital review. A launched handoff is not evidence that
the serial diagnostic or a VPT2 retry has converged; inspect their current
execution receipts separately.

## C + SH3 and serial follow-up at 13:25 PDT

The new geometry audit shows all three S-H bonds retained at 2.0-2.6 A:
the shortest C-H distance grows from 2.4290 to 2.8746 A. At 2.8, 3.0 and
3.2 A, both the primary and fragment-seeded attempts instead form C-H near
1.1098-1.1099 A and stretch the corresponding S-H to 2.9422-3.2559 A.
The running 3.4 A trajectory has C-H 1.1099 A and S-H 3.4108 A already.
These are rejected H-transfer structures, not usable C + SH3 scan points.
No grid, dividing surface, orbital model or review approval is changed.
The earlier intact geometries are not thereby certified as constrained
minima or as a complete physical dissociation path.

After a gap the primary already starts from validated asymptotic fragments.
Previously a connectivity rejection could trigger another optimization from
that same geometry and orbital/state seed. The scan now suppresses that
redundant fallback, records `fragment_retry_skipped`, and continues subsequent
points while retaining the gap and whole-path review. Geometry comparison
uses native input precision, with the complete state model required to match.
A different seed still gets its bounded retry. Existing retry directories
remain eligible for normal reattachment (including interrupted preparations,
allocated jobs and explicit live retries); explicit primary-task fallbacks
are not deduplicated from their original logical seed.

Serial job 52919565 remains active with no new failure receipt in this read.
Its maximum force decreased from 11.2188 to 0.0080 eV/A by step 12, with
CH3 intact, but the requested threshold is 0.0005 eV/A. The last three
completed central energy evaluations identify RSPT2 STATE 2.1 and agree with
the saved high-precision energy. This does not yet establish convergence,
displaced-state continuity, or the cause of the earlier MPI abort.

### Separate fragment-branch diagnostic

`python -m kinbot.vrc_branch_probe run SCAN_REVIEW_JSON --distance 6`
starts one independent constrained relaxation from the hash-verified, completed
asymptotic geometry and final orbital checkpoint. It preserves the native
method, basis, state average, selected root, convergence thresholds, and internal
fragment geometry. It does not freeze additional bonds or publish production
scan geometries, corrections, or approvals.

The default budget is 50 steps and one hour in one allocation under the existing
production shared pool. The existing bounded fragment-retry dispatcher policy
prevents automatic task/walltime extensions. Explicit shared I/O recovery, if
triggered inside that allocation, remains subject to its walltime. Repeating the
same command reattaches to the same graph, including a failed graph. Source
checkpoint, geometry, task and context hashes are checked before staging; an
explicitly modified source task is refused. Diagnostic graphs live separately
under `diagnostics/vrc_branch_probes`.

Use `python -m kinbot.vrc_branch_probe status PRODUCTION_ROOT` for a read-only
summary of task IDs, recent optimizer steps and, on completion, all pair
distances. A converged intact result establishes an outer constrained solution
at that orientation/model; it does not establish a barrierless channel or
active-space completeness. A transfer or budget failure likewise does not prove
that no intact solution exists. Keep the production path under review until
the branch and electronic-state continuity are established.

### Serial result and monitoring after reattachment

Job 52919565 subsequently completed and its original graph was marked
`complete`. The saved final projected force was 0.0002543938466342233 eV/A,
below the requested 0.0005, with constraint residual 4.688294197308227e-10.
This verifies the serial workaround at this point; it does not identify the
earlier MPI abort's root cause. Sella's printed final force can lag the final
convergence check, so the saved final force/constraint test is authoritative.

Dispatcher completion and VRC acceptance are separate. The scan must still
validate connectivity and its orbital checkpoint, publish the geometry, and
record the point as `accepted`. The r1 log is not evidence of acceptance or
continuation on r2. Resource-only retries retain their graph identity and
do not change the default resources for subsequent points.

The standard monitor now includes `diagnostics` alongside the prepared PES
cases, lists allocations elsewhere under the monitored root, and displays
current-pass versus saved scan progress. Replayed acceptance can reuse QC;
it is not a count of new calculations. Diagnostic convergence remains
separate from production-model acceptance. Updating this read-only monitor
requires a code pull, not replacement of the running controllers.

### Bounded retries and unattended continuation

The 14:34 PDT snapshot confirms r2 point 14 (4.8 A) and point 15 (5.0 A)
are accepted, with point 16 (5.2 A) running. The serial result has therefore
been consumed by the production scan. The separate C + SH3 6 A probe ended
with `ASE geometry optimization did not converge` at its 50-step budget.
It remains a failed diagnostic, not an approved outer point. Its optimizer
tail alone does not establish final connectivity or state continuity.

Review of the failure path exposed a hold interaction: after an explicitly
bounded fragment retry failed, `_drive_vrc_dispatch` could honor the global
`hold_task_failures` setting forever instead of returning the failure to
`ScanProgress`. Graphs marked `vrc_fragment_retry` now return their terminal
failure immediately, without another submission. The scan records the gap
as `needs_review`, retains the failed attempt, and continues from the verified
fragment seed at later distances. The complete-grid requirement still blocks
production use of the incomplete path. Ordinary VRC tasks retain their failure
hold and existing transient-retry behavior. This controller change requires
replacement of running controllers; `run_refresh ROOT --apply` alone retains
existing QC allocations, task results and scientific inputs.

The same snapshot has 21 allocations under a cap of 23 and no reserved
orbital-review slots. Two completed receipts are intentionally deferred:
CH3 `anl1_6z` job 52920059 awaits outer-PES L2 selection, while CH3SS scan
point 5 job 52919120 awaits its replacement worker returning from reaction
search to VRC. Neither is evidence for rerunning the native calculation.
The C2V CH2 and C1 VPT2 recoveries are staged behind the ANL gate; their
automatic framework-recovery path remains enabled. Historical rejected TS
searches and the C + SH3 connectivity gaps remain explicit scientific
limitations, rather than being relabeled as successful for the handoff.

Validation: 82 tests passed across VRC dispatch, scan recovery, parallel
channels and branch probes. The regressions cover an exhausted native
fragment retry under a global hold, continuation to the next point without
accepting the gap, and preservation of ordinary task holds.
