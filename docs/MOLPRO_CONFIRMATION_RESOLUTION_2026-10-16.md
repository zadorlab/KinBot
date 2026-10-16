# Molpro confirmation: evidence, resolved cause and targeted recovery

**CABS=3 update at 11:42:40 PDT:** the controller validates the CH3S six-root
far point; the S–S twelve-root far point is running with output age 18 s and
no reported error. Its last line does not identify which root is active.
The optional `monitor --roots` view now scans emitted energy markers and
program headers to show that context. Full recovery and CABS=3 native exports
are not yet present locally; independent energy review remains pending.

**Latest native-tail update, October 16:** both CABS=2 asymptotes have failed
after converged MRCI. CH3S completes three doublet calls and fails on quartet
root 1; S–S fails on doublet root 1. Both report `strange..` in
`mrf12_singles.F90`, with matching NSTATI/NSTATR. CABS=2 is rejected for
production in these tested conditions. The failed files remain evidence;
their pre-CABS energies are not final results. The CASPT2 recovery continues
independently; the later monitor accounts for MAXDAV20 below. Two isolated CABS=3 asymptotes are
prepared below; neither launches a dependent scan nor establishes a proven fix.

**Subsequent monitor at 10:12:32 PDT:** 19 validated, two failed, four waiting,
one running. Only S–S aTZ shared-reference CASPT2 `return` is running. Its
`tilted` successor accounts for one waiting point; the three CABS dependents
account for the others. With the known 26-point plan, this implies the
MAXDAV20 point completed and 18 CASPT2 points completed. These are inferences
from the supplied monitor/plan, not a new native energy audit. Retain the last
CASPT2 chain until it finishes, then collect the recovery export. The
[10:12 snapshot](../examples/anl/orbital_reviews/molpro_confirmation/recovery_status_20261016_101232.json)
keeps this progress separate from the two error reports.

The earlier 09:45:56 PDT monitor reported ten validated supplemental points,
five running, ten waiting and one failure. No newer complete count is inferred
from the two-job report itself. The historical
[snapshot record](../examples/anl/orbital_reviews/molpro_confirmation/recovery_status_20261016_094556.json)
and the subsequent [native-tail report](../examples/anl/orbital_reviews/molpro_confirmation/cabs2_failure_report_20261016.txt)
are retained separately. The [structured review](../examples/anl/orbital_reviews/molpro_confirmation/cabs2_failure_review_20261016.json)
records their provenance and the comparisons with the original native export.

The selected minimum spaces remain **CH3S C–S (5,4), CH3SS C–S (7,5),
CH3SS S–S (7,5)**. The original `_r2` run has 35 validated points, two failed
points and two unsubmitted dependents. Numerical completion is not approval
of the old CASPT2 reference policy. The policy responsible for the large
CH3S asymptotic splitting has a demonstrated remedy in an existing control.
The supplemental run tests that remedy on the current geometries/bases and
completes the MRCI checks. No old calculation is overwritten or relabelled.

## Evidence and reproducibility

Uploaded archive:
`examples/anl/orbital_reviews/minimum_active_space_molpro_confirmation_20261016_r2.evidence.tar.gz`.
SHA-256: `4c80871bfb7ddb8f20d930a076bfabe8906bdba0df10b8d76c168f47ac2625fb`.
All 833 inventory files passed size/hash checks. The export was captured at
2026-10-16T16:02:59.076892+00:00. Its original absolute paths are mapped through
`export_inventory.json`; the native evidence itself is not rewritten.

The extracted native evidence is retained in
`local_runs/molpro_confirmation_r2_review/`. The full archive stays outside
the Git payload. Reproduce the independent review from the repository root:

```bash
.venv/bin/python examples/anl/orbital_reviews/molpro_confirmation/review_export.py \
  local_runs/molpro_confirmation_r2_review \
  --output examples/anl/orbital_reviews/molpro_confirmation/native_r2_review.json
```

[The machine-readable review](../examples/anl/orbital_reviews/molpro_confirmation/native_r2_review.json)
contains native energy arrays, active orbital tables, return comparisons,
output/WFU hashes and error sentences. All 35 completed outputs have normal
termination, completed execution records, completion markers, the expected
root counts, matching manifest energies, and no detected native warnings or
convergence errors. Inherited variables printed by `PROGRAM * RESTART` are
excluded from energy extraction. These checks do not prove electronic-state
identity by CI-vector overlap, nor certify every possible production geometry.

## The CASPT2 splitting has a demonstrated remedy

The report supplied by the user comes from the older reference campaign,
`18_pilot/avt_f12__cs_r_30`, marked complete. It is a separate campaign from
the uploaded 39-point confirmation export. Its extracted values are retained
in [prior_reference_controls.json](../examples/anl/orbital_reviews/molpro_confirmation/prior_reference_controls.json),
explicitly labelled user-supplied output, not independently checksummed native
files. No new IPEA calculations are needed to establish the following result.

| IPEA / real shift, Eh | Correlated reference density | Doublet spatial span, kcal/mol | Quartet spatial span, kcal/mol |
|---|---|---:|---:|
| 0.25 / 0.30 | Selected root | 1.174804 | 1.174735 |
| 0 / 0.30 | Selected root | 0.000151 | 0.000207 |
| 0 / 0.10 | Selected root | 0.000151 | 0.000220 |
| 0.25 / 0.30 | Three equal-weight same-spin references | 0.000144 | 0.000144 |
| 0 / 0.30 | Three equal-weight same-spin references | 0.000144 | 0.000151 |

This directly contradicts the earlier statement that the large CH3S splitting
survived removal of IPEA in this comparison. That statement is withdrawn.
The ~1.175 kcal/mol split between equivalent sulfur components is not the
physical CH3S fragment excitation relevant to the different S–S channel.
The control changes the correlated reference treatment without enlarging the
CAS and restores the expected degeneracy while retaining IPEA 0.25.

**Chosen remedy:** keep IPEA 0.25 and real shift 0.3; use the declared complete,
equal-weight same-spin reference manifold in the CASPT2 call. For CH3S this is
`REFSTATE,3,1,2,3; WEIGHT,1/3,1/3,1/3` (the generated input uses decimals).
For S–S it is six roots per spin with weight 1/6; both CH3S fragment components
remain in that manifold. CH3SS C–S has one root per spin, so this change is
identical to its existing treatment. The CASSCF averages and orbital spaces
are unchanged. This is not a switch to MS/XMS-CASPT2 or a different energy
selector, and it does not force excited fragment states to the same zero.

Why the remedy is physically appropriate: Molpro constructs H0 from the
selected state's density unless an averaged reference density is requested.
Its IPEA implementation is not exactly orbital invariant. Equal averaging over
the complete sulfur 3p manifold gives an isotropic density in that subspace;
using its common H0 avoids assigning different corrections to equivalent
components. The controlled results support this explanation; they do not
constitute a source-code diagnosis of every intermediate matrix.
[Molpro CASPT2 documentation](https://www.molpro.net/manual/doku.php?id=multireference_rayleigh_schroedinger_perturbation_theory).

Size consistency is a different question: whether separated-fragment energies
are additive. A common MRCI size-consistency error can leave spatial components
degenerate. It does not imply that MRCI must split them more than CASPT2.
Also, conventional CASPT2 is not unconditionally exactly size consistent for
all open-shell separations; the manual discusses modified H0 treatments for
that problem. Neither observation justifies the present nonphysical spatial
splitting. The relevant direct evidence is the controlled table above.

The peer-reviewed XMS work explains why invariance to rotations of reference
states matters near degeneracies and why an averaged Fock operator and balanced
reference treatment are useful. It supports the conceptual distinction; it
does not prove that ordinary shared-reference RS2-F12 is equivalent to XMS or
that an untested XMS-F12 switch is supported here.
[Shiozaki et al., JCP 135, 081106 (2011)](https://doi.org/10.1063/1.3633329),
[Shiozaki, Woywod and Werner, PCCP 15, 262–269 (2013)](https://doi.org/10.1039/C2CP43381H).

In the current `_r2` output, CH3S aTZ CASSCF and MRCI+Q within-spin spans are
approximately 0.00003–0.00004 kcal/mol, whereas selected-root CASPT2 gives
about 1.175. In S–S, the two triplet centroids are separated by about 1.287
kcal/mol at aTZ (about 1.316 in the local ORCA check), while each atomic-sulfur
triplet has a much smaller ~0.011 kcal/mol residual spread. The centroid gap
represents different CH3S fragment components; the within-triplet spread is
what the shared-reference check must remove. Do not conflate those quantities.

## What the successful orbitals establish

All labels below are **one-based Molpro C1 labels**. `2140.2` is the natural
orbital record and `2141.2` the pseudo-canonical record. Assignments refer to
subspaces; rotations inside a degenerate pair/triplet are allowed.
Basis-specific file paths, hashes and assignments are in
[native_orbital_orderings.json](../examples/anl/orbital_reviews/molpro_confirmation/native_orbital_orderings.json).

| Channel | CLOSED / OCC | Natural active order at 30 Å | Pseudo-canonical active order at 30 Å |
|---|---|---|---|
| CH3S C–S (5,4) | 10 / 14 | 11–13: departing S 3p; 14: methyl radical | Same groups, internal S rotations/permutations |
| CH3SS C–S (7,5) | 17 / 22 | 18–19: S2 pi; 20–21: S2 pi*; 22: methyl radical | Same groups, rotations within each pi pair |
| CH3SS S–S (7,5) | 17 / 22 | 18–19: CH3S nonbonding/radical pair; 20–22: departing S 3p | **18–20: departing S 3p; 21–22: CH3S pair** |

These groups occur in both aDZ and aTZ outputs. The S–S reordering between
natural and canonical records remains entirely inside the five-orbital CAS.
It is not loss of an orbital from the active subspace. For C–S in CH3SS, the
original seed's 17.1↔19.1 swap has already been applied in the completed WFU;
do not apply it again when starting from that WFU. Orbital indices have no
portable meaning without the corresponding basis, geometry, file and record.

The six center→inner→return checks have minimum active-subspace singular value
0.999999999992. This uses the full native MO expansions at identical geometry
and AO basis, with rank, reconstruction and orthogonality checks. It assumes
Molpro's MO normalization; it is not an independently integrated AO-overlap
calculation and is not a cross-geometry CI-state overlap.

| Channel | aDZ maximum CASPT2 return discrepancy | aTZ maximum CASPT2 return discrepancy |
|---|---:|---:|
| CH3S C–S | 0.007474 kcal/mol | 0.002887 kcal/mol |
| CH3SS C–S | 0 at printed precision | 0.0000063 kcal/mol |
| CH3SS S–S | 0.0000126 kcal/mol | 0.0000126 kcal/mol |

There is no detected active-subspace loss or core intrusion in these tested
points. The new common-H0 checks repeat the energy return tests because changing
the correlated treatment requires verifying its resulting surface, even when
the underlying CASSCF orbital controls remain unchanged.

## MRCI failure 1: external CI stagnation, not an I/O failure

Native file: `ch3s_cs_atz_mrci/atz__tilted/atz__tilted.out`, job 52908112.
The solve reaches 2000 iterations with energy near -437.41629999 Eh,
VAR(S) 0.94e-9 and VAR(P) 1.0e-9, above its 1e-10 squared-coefficient target.
The small energy changes do not satisfy the independent vector criterion.
Native MAXDAV is only 4; MAXVI=60 enlarges a different, internal space.
Near-degenerate reference energies and sustained cycling are consistent with
an inadequate external iteration subspace. They do not alone prove root flipping.

**Remedy:** `max_external_vectors: 20` → `OPTION,MAXDAV=20`; retain CABS=1,
2000 iterations, MAXVI=60, reference roots and existing convergence targets.
Start from the completed center checkpoint, not from the failed final WFU.
This enlarges the solver search space without changing the intended Hamiltonian.

The Davidson algorithm builds an iterative subspace for large eigenproblems;
additional retained directions can resolve slow/competing roots. The original
paper supports the algorithmic rationale, not a universal optimal value of 20.
[Davidson, J. Comput. Phys. 17, 87–94 (1975)](https://doi.org/10.1016/0021-9991(75)90065-0).
The particular MAXDAV=20 trial also follows a Molpro developer's advice for
an MRCI convergence/degeneracy case; that advice is primary developer guidance,
not a peer-reviewed guarantee for this system.
[Werner, Molpro user discussion, May 2018](https://www.molpro.net/pipermail/molpro-user/2018-May/007247.html).

## MRCI failure 2: the post-MRCI singles-CI correction

Native file: `ch3ss_ss_atz_mrci/atz__far/atz__far.out`, job 52908084.
The MRCI solve converges in eight iterations, giving -835.141055331546 Eh,
and prints the pre-CABS relaxed-Davidson value -835.21982764 Eh. It then
aborts with an orthogonalization error in `mrf12_singles.F90` before a singles
correction and final selected-energy marker are produced.

The earlier AO/CABS construction succeeds: 348 CABS vectors are retained,
with reported minimum CABS metric eigenvalue 2.69e-7. The reference and external
CI phases also finish. Thus the failure is localized to the later singles-CI
correction, rather than the initial auxiliary basis construction, Slurm, an
iteration limit, or demonstrated memory exhaustion. Linear dependence or
near-redundancy in that later space is a plausible mechanism, but the log
does not print its failed matrix/eigenvalues. The exact internal cause cannot
be asserted from this output or the accessible documentation.

**Original trial, now unsuccessful:** use `cabs_singles: 2` for the complete S–S far→center→tilted
chain. This retains a CABS correction using the documented perturbational
form instead of the failing singles-CI form. It is an approximation change,
not just a numerical retry. Two already successful CH3S far/center geometries
are computed with the same alternative to quantify its effect on interaction
energies and spin gaps. No mixed CABS=1/2 interaction curve is constructed.

Molpro distinguishes external MAXDAV from internal MAXVI; CABS choices 1, 2,
and 3 mean singles CI, perturbational, and coupled treatments. Default 1 is
not being redefined globally. `CLUSTER=0; ENERGD(1)` remains the relaxed-Davidson
selection after CABS processing.
[Molpro MRCI documentation](https://www.molpro.net/manual/doku.php?id=the_mrci_program).

The multireference F12 literature develops both the correlated ansatz and
singles corrections, including balanced treatment of multiple references.
It establishes these as method choices; it does not diagnose this specific
2024.1 error or guarantee that choosing 2 avoids every internal orthogonalization.
[Shiozaki, Knizia and Werner, JCP 134, 034113 (2011)](https://doi.org/10.1063/1.3528720),
[Shiozaki and Werner, JCP 134, 184104 (2011)](https://doi.org/10.1063/1.3587632),
[Shiozaki and Werner, Mol. Phys. 111, 607–630 (2013)](https://doi.org/10.1080/00268976.2013.779393).
The author-hosted abstracts and bibliography were accessible; the complete
publisher text of these three papers was not obtained during this review.
The option semantics above come from the official manual, not an invented
line-by-line interpretation of an inaccessible paper.

The successful five MRCI jobs contain 18 completed correlated roots. Every
selected energy agrees, within printed rounding (2e-8 Eh), with **relaxed**
Davidson plus its CABS correction. The fixed/rotated corrections are not read
by mistake. The parser and selector are unchanged. The failed far-point's
pre-CABS energy is not accepted as a finished MRCI-F12+Q result.

### CABS=2 outcome and the bounded CABS=3 test

| Far-point calculation | Failing call | NSTATI / NSTATR | Pre-CABS MRCI, Eh | Difference from original CABS=1 pre-CABS MRCI |
|---|---|---|---:|---:|
| CH3S C–S, CAS(5,4) | Fourth call; quartet root 1 | 3 / 3 | -437.405556015739 | -3.90e-10 Eh |
| CH3SS S–S, CAS(7,5) | First call; doublet root 1 | 6 / 6 | -835.141055331544 | +2.0e-12 Eh |

The first row rules out a blanket inability of CABS=2 to run any of these
calls; the second rules out a quartet-only restriction. Both have the full
declared reference manifold and the already repaired NSTATI/NSTATR equality.
Their preceding MRCI energies reproduce the original calculation. The error
therefore lies in the later correction stage, without evidence of a changed
CAS, failed MRCI convergence, or cluster interruption. The earlier CABS basis
construction completed with 273/348 vectors, respectively. Its metric
eigenvalues do not diagnose the unprinted matrix in the failing later routine.

The generic `strange..` message gives no condition, matrix, or offending
index. Neither the manual nor accessible literature identifies this guard.
An exact explanation requires the matching Molpro implementation or a
developer diagnosis. Calling it a specific linear-dependence, spin, or
reference-count bug would exceed the evidence. Repeating CABS=2 or increasing
MRCI iterations would not address the observed failure.

**Next controlled test:** `cabs_singles: 3`, documented as including couplings
between CABS singles and the MRCI terms. This changes the correction treatment;
it does not disable CABS or change the chosen CAS. The manual and the 2013
review cited above justify considering this method, not claiming that it
bypasses the failing routine. It may fail too. The
[two-point plan](../examples/anl/orbital_reviews/molpro_confirmation/cabs3_gate_plan.json)
uses the same original seeds, basis, geometry, root/reference lists, thresholds,
SCALF12=1 and CLUSTER=0. MAXDAV remains unchanged for these controls. It tests
all six CH3S and twelve S–S roots at 30 Å and contains no dependent points.
No general KinBot method or parser change is needed.

If both points complete, review every final selected energy, the native CABS
treatment, state manifold and relaxed-Davidson convention before preparing
the missing S–S interaction points and a matched CH3S interaction control.
CABS=3 includes coupling, so do not assume its energy is a pre-CABS number
plus the separately printed correction used for CABS=1. The old CABS=1
energy-addition review is not a validation of CABS=3 output. A completed
asymptote alone cannot establish accurate relative energies or spin splittings.
If either gate fails, retain the exact native reproducer for implementation
diagnosis; do not start another full scan or silently omit the correction.

## Waiting points and closure conditions

S–S MRCI `atz__center` requires `atz__far`; `atz__tilted` requires `atz__center`.
They had no submitted allocation/native output in `_r2`. Their waiting state
followed the failed predecessor. In the recovery, the CABS=2 S–S far point
also failed, so its center and tilted dependents remain blocked. The failed
CH3S CABS=2 control blocks its center too: **three dependency-blocked points**
in that supplement. Other waiting CASPT2 points follow their still-progressing
orbital chains. Neither waiting label itself represents a new Molpro failure.
The separate CABS=3 tests will not unfreeze or mutate those old chains.

The [supplemental plan](../examples/anl/orbital_reviews/molpro_confirmation/recovery_plan.json)
has 26 points in seven chains: four unresolved MRCI points, two matched CABS
controls, and 20 common-H0 CASPT2 points (two affected channels, two bases,
five points including the exact return). CH3SS C–S CASPT2 and its three
successful MRCI points are reused. There is no new large-CAS or IPEA series.
The original CASPT2 chains took roughly 3–32 minutes each; shared-reference
cost and scheduler delays can differ. No completion time is guaranteed.

Resolution requires actual completed outputs, not merely input preparation:

1. MAXDAV20 must converge all six requested tilted MRCI roots, retain the
   selected subspace/state manifold and complete the final CABS/ENERGD steps.
2. A CABS treatment must complete both all-root asymptote tests before
   preparing a consistent S–S far/center/tilted chain. The paired CH3S
   control must quantify changes in dissociation-referenced interaction and
   doublet–quartet gaps, not merely total-energy offsets. A material change
   requires adopting one consistent reference protocol across channels and
   rechecking its correction; it cannot be hidden by mixing old/new zeros.
3. Shared-reference CASPT2 must recover the sulfur-component degeneracy in
   both bases, retain the distinct CH3S fragment gap in S–S, and reproduce
   inward/return energies and active subspaces. It must use the specified
   same-spin reference weights in the **native output**. Compare the resulting
   angular interactions and spin gaps with the completed reference points.
4. Only then install the accepted basis-specific records/reference policy in
   production JSON. Numerical execution alone does not authorize a surface
   with the old 1.175 kcal/mol artifact. No independent per-root zero shifts
   are used to make that artifact disappear.

The current choice for the CASPT2 remedy is supported by completed control
data. CABS2 and MAXDAV20 remain specific native tests until the user executes
them on Blodgett; claiming that they have already succeeded would be false.
