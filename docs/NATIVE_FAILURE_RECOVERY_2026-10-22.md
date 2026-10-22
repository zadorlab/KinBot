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

### Independent sampling and live recovery

Managed PES workers now launch each complete VRC surface from the case's
canonical `rotdpy_directory`, including an individual physical root whose
sibling remains under review. The complete correction and its scan gates are
required before launch. Every sampling owner uses the existing shared pool.
The final PES aggregator generates the identical input and adopts the same
campaign. Surface numbering is local to each input, independent of completion
order. Changed inputs or correction hashes cannot overwrite existing samples.

The monitor displays each campaign, its executor status and its correction
review status. For a generated production input, native curve/state review
can occur while sampling runs. Approval records the exact input/correction
contract and does not stop or restart sampling:

```bash
python -m kinbot.rotdpy_campaign accept /path/to/case/rotdPy/SURFACE.py \
  --note 'Describe the completed native correction and state-continuity review'
```

This approval is required before the sampled number-of-states table enters
MESS. It is not a substitute for successful QC, resolved scan points, correct
connectivity, the production profile, or valid native sampling statistics.
If review changes the scientific model, the old samples are retained and a
new sampling directory is required. Approval of unchanged evidence does not
interrupt the running executor.

An exhausted or failed sampling owner can be resumed without replacing PES
or KB controllers:

```bash
python -m kinbot.rotdpy_campaign status /path/to/case/rotdPy/SURFACE.py
python -m kinbot.rotdpy_campaign retry /path/to/case/rotdPy/SURFACE.py
```

The existing bounded native sample retries and shared I/O retries remain
automatic. `retry` above resumes the same hash-bound database after those
budgets are exhausted, never starts another live owner, and never loops on
unclassified failures. The final aggregator waits for live recovery when
`hold_task_failures` or `hold_failures` is enabled. It always waits for required
correction review. Existing `run_control retry` task/QC/reaction/conformer
requests continue to be consumed by their live owners.

Dispatcher-managed Molpro force jobs now have one automatic serial resource
fallback for a recorded communication-error/MPI-abort combination, with the
same method, orbital model, geometry and optimizer settings. Stage/execution
provenance and inactive allocation are checked, and the failed attempt is
archived. Generic aborts, explicit I/O/OOM/convergence failures, used fallback
budgets and deliberately bounded fragment probes do not trigger this rule.
The already accepted r7 serial result is retained. VPT2 framework recovery,
MRCC exact-rank parsing, timeout recovery and optimizer continuation retain
their existing evidence checks and finite budgets.

### Reconciliation of search rejections supplied after 15:04 PDT

CH3S's 35 rejected searches consist of 15 wrong-IRC-endpoint searches,
8 returning to the initial well, 7 invalid path geometries, 2 exhausted
750-step optimizations, 2 high-level TS identity mismatches and 1 failed
stationary-point frequency check. CH3SS's initial well has 13 rejected
searches: 3 wrong-endpoint, 9 initial-well returns and 1 exhausted optimizer.
Other CH3SS wells were still exploring; these are not final PES totals.

The seven previously highlighted unsuccessful production attempts overlap
these counts: five belong to CH3S's 35 (one frequency rejection, two exhausted
optimizations, two invalid H2-elimination path guesses), one belongs to
CH3SS's 13, and one is a failed TS conformer rather than a whole-reaction
rejection. Successful sibling conformers must still pass normal validation.
These search outcomes do not block unrelated calculations. An inconclusive
search is not evidence that its proposed physical channel cannot exist.

The ad-hoc shell extraction compared each exit receipt against the literal
`0`, incorrectly labeling successful `JOB_ID 0` receipts as failures. The
maintained `run_pes_audit` now uses the authoritative generation-aware worker
receipt parser, includes conformer failures, and reports explicit rejection
categories. Its JSON retains the individual log evidence. Reports remain
read-only with respect to the run and are written outside it.

The C+SH3 well's two tight H-transfer discoveries have the same path class,
near-identical energies and imaginary modes around -747 cm^-1. Outer PES L2
selection retains the reverse discovery
`471230610040000000002_intra_H_migration_2_5`; the tight path is represented.
That path is distinct from dissociation to C + SH3. The separate 6 A probe's
last frame retains S-H distances 1.32892, 1.47851 and 1.47848 A, with all C-H
distances above 6 A. It exhausted its 50-step optimization budget. This is
evidence of intact geometry at that saved iterate, not a converged outer
branch, nor evidence that C + SH3 is absent. The production channel remains
held while independent surfaces may sample.

### Read-only CH3SS comparison with the paper

`python -m kinbot.paper_pes_check ROOT` compares the current CH3SS discovery
records with Fig. 6 of [DOI 10.1039/d4cp04873c](https://pubs.rsc.org/en/content/articlehtml/2026/cp/d4cp04873c).
Add `--json` for the complete per-search evidence. It prints to stdout, reads
existing files and ASE databases in read-only mode, and does not call Slurm,
dispatch, recovery, or a QC program. It can be run before or after a controller
refresh, and does not need controllers to be restarted.

The comparison covers W1-W4, P1-P4, TS1-TS6, entrance association and W4
scission. It matches element/H connectivity and unordered endpoints, with
separate tight-TS versus homolytic-scission classifications. The paper's W2
has C-S-S connectivity and two hydrogens on carbon; W3 has S-C-S connectivity.
Bond-order, spin and conformer identity are not inferred from this match.
Paper energies are reference labels, not numerical acceptance thresholds.

Accepted searches, pending candidates, rejected candidates and outer-PES
selection are reported separately. Current monitor acceptance takes precedence
over older plans. Product identities missing from manifests can be reconstructed
from existing native well geometries; the latest failed cache row never revives
an older success. W1 requires a recorded complex, and its particular paper
orientation remains unverified. Separated CH2S + SH alone does not establish W1,
and matching net endpoints does not prove that the same paper TS was found.

Additional accepted species/connections and unfinished extra wells/candidates
are listed separately. Multiple searches are retained as evidence, not counted
as unique saddles. Missing or changing files and unmapped endpoints remain
visible; a snapshot with such gaps cannot establish that a paper feature is
absent. Normal search rejections do not become new MESS approval requirements.

Within the configured L1/L2 discovery protocol, rejected candidates and failed
individual conformers are terminal search outcomes. They do not each require
repair before assembling the retained network. Required ANL components,
selected unresolved VRC channels, correction review and valid sampling results
remain explicit completion requirements for that network.

Validation: the controller/recovery regression suite passed 420 tests with one
skip. Ten additional paper-comparison tests cover connectivity, reverse paths,
complex versus separated fragments, path type, stale selection, missing records,
duplicate searches, unfinished candidates, and read-only native cache access.
