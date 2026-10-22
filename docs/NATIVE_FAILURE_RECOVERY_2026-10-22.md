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
