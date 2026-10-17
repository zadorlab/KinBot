# S–S simultaneous MRCI: first-iteration DSYGV failure

## Evidence and current status

**Closed after the user's 17:32 PDT snapshot:** the Molpro23 build probe
also fails with `dsygv failed` in `cif12solve` after the first iteration,
with zero of twelve correlated markers (job 52908248, node20). The reported
reference energies and first-iteration diagnostics match the Molpro24 retry.
This provides no working alternate build and does not isolate the exact
internal cause; both attempts ran on node20. The user chose CASPT2-only
production. No further MRCI retry is required or scheduled, and MRCI is no
longer a production prerequisite. Preserve this evidence and the earlier
successful benchmarks. The new result is from the pasted monitor; no full
Molpro23 native archive has been independently reviewed.

**Update at the user's 16:45 PDT snapshot:** the 0.1 Eh denominator-shift
trial also failed. Job **52908243**, node20, in
`_r2_multistate_solver_retry/ch3ss_ss_atz_mrci_preconditioned/atz__center`
again reports `dsygv failed` in `cif12solve` immediately after iteration 1 of
the six-root doublet call, with 0/12 correlated markers. The last printed
root-six denominator is now -0.49172383 Eh (previously -0.47453845), but the
failure location is unchanged. The trial has **one failed, one waiting**;
tilted has not run because center supplies its checkpoint. This is a failed
preconditioning hypothesis, not a completed repair. No completed result has
been invalidated or deleted. This update is from the pasted monitor; the
full failed native input/output still needs export.

The user supplied the native tail and solver context of job **52908238**,
`_r2_multistate/ch3ss_ss_atz_mrci_simultaneous_cabs3/atz__center`, on node20.
Molpro reports **2024.1, linked January 31, 2024, 64-bit MPP**. The previous
far point completed according to the controller; its full native export still
needs independent local review. Four CH3S points already passed that review.
Thus the seven-point supplement has **five controller-complete, one failed,
one waiting**. Tilted needs center's completed orbital checkpoint. It has not
failed an electronic-structure calculation.

The center input and output agree on six doublet states and references:
`NSTATE=NSTATI=NSTATR=6`, CABS=3, SFIX, NOSWAP, THRDEN=1e-8,
THRVAR=1e-10, MAXVI=60. Input MAXDAV=20 is printed as 24 by this executable.
The native calculation has 26,112,200 contracted configurations. It completes
reference rotation and F12 integral setup, prints iteration 1 for all six
states (zero correlation energies, VAR(S)=0.039–0.043, VAR(P)≈0.14), then
aborts with `dsygv failed` in `cif12solve`. No correlated selection marker
is produced. This is earlier than the former post-CI `mrf12_singles.F90`
failure. MPI abort/killed-rank messages follow the native error.

These observations come from the pasted native report, not an independently
hashed export of this point. Do not invent its file hash or claim the new
S–S far spectrum has been scientifically accepted.

## What is known about the cause

LAPACK DSYGV solves a real symmetric generalized eigenproblem. Its numeric
INFO distinguishes invalid arguments, eigensolver nonconvergence, and failure
of the positive-definite factorization of the metric. Molpro did not print
INFO here. An ill-conditioned or non-positive metric is a plausible mechanism,
but the message alone does not establish it or identify which internal matrix
was passed. [LAPACK DSYGV source and INFO contract](https://www.netlib.org/lapack/lapack-3.1.1/html/dsygv.f.html).

There is no evidence here for an iteration-cap failure, I/O/permission failure,
out-of-memory kill, or a missing reference state. Raising MAXIT, excluding
node20, repeating the old IPEA series, or accepting the printed reference
energies would not address the reported failure. Nor does this solver abort
establish that CAS(7,5) is chemically inadequate. The same root/method setup
passed the far point, but that does not guarantee numerical robustness at an
interacting geometry. The exact implementation-level cause remains unresolved.

## Historical preconditioning trial and rationale

The only new electronic-input options are:

```text
shiftp=0.1,shifts=0.1,shifti=0.1
```

Molpro documents these as shifts in Davidson update denominators for pairs,
singles and internal configurations, respectively. They act on the iterative
solver, rather than the CASPT2 level shift. The manual states that ordinary
MRCI options also apply to MRCI-F12.
[Molpro MRCI denominator shifts](https://www.molpro.net/manual/doku.php?id=the_mrci_program#denominator_shifts).

**Scope of the imaginary-shift advice.** The manual also advises zero real
SHIFTS/SHIFTP and imaginary SHIFTSI/SHIFTPI if needed under its *Perturbative
F12 correction* section. That section describes Brueckner/orbital-optimized
CCD/DCD followed by a perturbative correction (the example uses
`DCSD-F12B,PERTCALC=2`). It is not the present variational MRCI-F12 calculation.
Do not transplant those imaginary-shift keywords into MRCI without documented
support. [The exact F12 section](https://www.molpro.net/manual/doku.php?id=explicitly_correlated_methods#perturbative_f12_correction).

Preconditioning is an established way to alter Davidson convergence; Morgan
and Scott analyze this connection in *SIAM J. Sci. Stat. Comput.* **7**, 817–825
(1986). The accessible publisher abstract supports that general mechanism,
not a repair for this Molpro error.
[Morgan and Scott, DOI 10.1137/0907054](https://epubs.siam.org/doi/10.1137/0907054).

**Our inference/trial choice:** changing the update denominators changes the
trial vectors entering later projected solves and may avoid this numerical
breakdown. The value 0.1 Eh is the user-requested modest first trial, not a universal
standard, a published recommendation for this case, or a guaranteed cure. A converged result targets
the same MRCI eigenproblem, so the successful unshifted far energy remains the
appropriate reference. No new full far calculation is scheduled.

The trial retains CAS(7,5), six roots per spin calculated simultaneously,
CABS=3, SFIX, reference policy, basis, geometries, tolerances, and
`OPTION,CLUSTER=0`/relaxed-Davidson indexing. It does not revive independent
excited-root calls that collapsed. It does not change pair-deletion thresholds
without evidence that they govern the failing matrix.

The [plan](../molpro_confirmation/mrci_preconditioner_plan.json)
prepares only **center → tilted**, in a new sibling directory. Before doing
so, the adapter checks that the source is stopped, revalidates its completed
records on temporary manifest copies, checks input/output/checkpoint hashes,
and verifies that the center still has the diagnosed error. It preserves the
original manifests and native files. Its retained `source_completion_snapshot.json`
is explicitly a completion/checkpoint check, not a scientific acceptance report.
The already-selected far WFU is copied without applying the old full-valence
orbital permutation again. Memory remains 104000 MB/node, eight ranks and
1000 MW/rank; concurrency is one. Other campaigns are unaffected.

The general JSON control is
`correlation_options.denominator_shifts={"pairs":0.1,"singles":0.1,"internal":0.1}`.
Only supplied keys are rendered; defaults remain unchanged. Values must be
finite nonnegative numbers. This control is rejected for CASPT2 and is never
enabled automatically on an I/O retry.

## Completion criterion and next branch

A successful trial must finish both spin blocks, yield all twelve energies,
pass native convergence/termination checks, and preserve the active subspace.
Review independent state character/overlaps, relaxed-Davidson selectors and
far-referenced interactions together with the completed far export. Passing
iteration 1 alone is not acceptance. Full production remains pending.

If the same first-iteration DSYGV error persists, stop this branch rather than
cycling arbitrary thresholds. The missing diagnostic is LAPACK INFO and the
projected matrix conditioning; obtaining it requires a diagnostic Molpro build
or developer support, or a controlled comparison with another available build.
No such build or native local validation is claimed here. Do not bypass the
failure with NOCHECK/IGNORE_ERROR or promote partial energies.

## Next discriminating check: one point with the installed 2023 build

The user reports modules `molpro19` through `molpro24`. The newest installed
module is not proof that it contains every subsequent upstream fix. The
observed failed executable identifies itself as **2024.1 linked 2024-01-31**.
The upstream change page lists later releases, but the public documentation
search did **not** establish a fix for this particular DSYGV failure.
Do not call this a known Molpro bug or promise a version change will cure it.
[Molpro release information](https://www.molpro.net/manual/doku.php?id=recent_changes).

The available 2023 build permits a bounded comparison without guessing more
electronic-structure parameters. The
[build-probe plan](../molpro_confirmation/mrci_build_probe_plan.json)
prepares **center only**, with the same two six-root spin blocks, copied from
the failed 0.1-shift trial. The launcher verifies identical prepared-input and
original-seed hashes, checks the generated batch setup's version, and records
the executable/module information. It removes the old explicitly prepended
executable directory before loading the new module. Local transformation
tests also verify byte-identical effective native decks. Both runs keep the
same MPI rank count, basis, SFIX, CABS=3, roots, tolerances and shifts. No
Molpro option is added as a purported DSYGV fix, and no tilted job is queued.

This compares complete builds, potentially including their math libraries;
it does not isolate a source-code regression. A checkpoint-reading or input
compatibility error in 2023 is a separate result, not reproduction of DSYGV.
A successful center must reproduce the intended CASSCF subspace/reference
model and finish all twelve correlated energies. It would justify reviewing
that executable as a workaround, followed by matched-build far/tilted checks;
it does not authorize mixing versions in a production correction potential.
If the same first-iteration error recurs, obtain INFO/matrix diagnostics from
a diagnostic build or Molpro support before another settings sweep.

The documentation identifies THRDLP with pair-overlap deletion and NPUPD
with the pair-update representation; neither identifies the matrix in this
abort. SCALF12=0 changes the ansatz, and lowering the root count changes the
modeled manifold. Those are not substituted automatically for the current
model. The multi-reference F12 formulation itself explicitly supports mixed
reference functions; near-degeneracy alone is not a reason to abandon the
required states. [Shiozaki and Werner, JCP 134, 184104 (2011)](https://doi.org/10.1063/1.3587632).

The [handoff](../../../../../../MOLPRO_CONFIRMATION_HANDOFF.md) contains the comparison,
monitor and export commands. Local tests cover unchanged default decks,
MRCI-only validation, both transformed trial inputs, copied-seed provenance,
preservation and refusal of changed/live source evidence; they are not native
Molpro convergence tests.
