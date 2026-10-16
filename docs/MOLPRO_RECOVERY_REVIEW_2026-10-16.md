# Native recovery review and MRCI root resolution

## Decision

The minimum active spaces remain **CH3S C–S (5,4), CH3SS C–S (7,5), and
CH3SS S–S (7,5)**. All 20 shared-reference CASPT2 recovery points and both
CABS=3 asymptote tests completed and passed independent native completion,
energy-selector and provenance checks. The CH3S tilted MAXDAV20 MRCI point
also completed numerically, but its separately requested roots contain strong
evidence of collapse. It is not an accepted distinct-state spectrum.

The seven-point MRCI supplement has stopped. **All four CH3S points have now
passed native review, including resolution of the tilted root collapse.**
The CH3SS S–S far point completed per the controller, but center failed at
iteration 1 with `dsygv failed` in `cif12solve`; tilted waits for center's
checkpoint. The new two-point [numerical trial](MOLPRO_MRCI_DSYGV_RECOVERY_2026-10-16.md)
preserves the completed far result. Full production awaits its independent
native review, the two interacting results and final model installation.

## Retained evidence and reproduction

The deposited archives remain in `examples/anl/orbital_reviews`:

- `minimum_active_space_molpro_confirmation_20261016_r2_recovery.evidence.tar.gz`:
  extracted to `local_runs/molpro_confirmation_recovery_review`; all **539**
  inventory files pass size/SHA-256 verification.
- `minimum_active_space_molpro_confirmation_20261016_r2_cabs3_gates.evidence.tar.gz`:
  extracted to `local_runs/molpro_confirmation_cabs3_review`; all **58**
  inventory files pass verification.

The recovery has **21 completed, two failed, three blocked** points. Twenty
completed points are CASPT2; the remaining one is the MAXDAV20 MRCI retry.
The two failed CABS=2 far points block S–S center/tilted and CH3S center.
Those unsubmitted dependents are not additional numerical failures. Both
separate CABS=3 gates completed. Historical failure statuses are retained.

Reproduce the full review from the repository root:

```bash
.venv/bin/python examples/anl/orbital_reviews/molpro_confirmation/summarize_resolution.py \
  local_runs/molpro_confirmation_recovery_review \
  local_runs/molpro_confirmation_cabs3_review
```

The committed [recovery audit](../examples/anl/orbital_reviews/molpro_confirmation/native_recovery_review.json),
[CABS audit](../examples/anl/orbital_reviews/molpro_confirmation/native_cabs3_review.json),
and [numerical findings](../examples/anl/orbital_reviews/molpro_confirmation/native_resolution_findings.json)
retain native paths/hashes, all energies, orbital tables, reference conventions,
return overlaps, error sentences and the duplicate-root evidence. No failed
pre-CABS energy is promoted to a completed result.

## CASPT2: the spatial splitting is resolved

IPEA remains **0.25 Eh**, real shift **0.3 Eh**. The change is the equal-weight
same-spin reference density: three references for CH3S C–S, six for CH3SS S–S.
The completed original CH3SS C–S checks already have one reference per spin and
remain applicable. These are separate from the CASSCF spin-block weights.

| Channel/basis | Maximum far-point sulfur-component span, kcal/mol | Largest center/return root difference, kcal/mol |
|---|---:|---:|
| CH3S C–S aDZ | 0.0000314 | 0.00000628 |
| CH3S C–S aTZ | 0.0000251 | 0.00000628 |
| CH3SS S–S aDZ | 0.0002887 | 0.0000126 |
| CH3SS S–S aTZ | 0.0002573 | 0.001908 |

For S–S, these spans are measured within each three-component sulfur group;
the two CH3S fragment components are **not** forced to be degenerate. Their
centroid gap remains about 1.300 kcal/mol in aDZ and 1.286 kcal/mol in aTZ.
The previous approximately 1.175 kcal/mol artificial CH3S atomic-component
split is absent. No new zero-IPEA calculation is required.

The largest return difference belongs to S–S aTZ quartet root 4. Its CASSCF
energy and active subspace reproduce; the difference appears in the final
CASPT2 energy, while the iterative energy reproduces much more closely. The
available output does not establish an internal implementation cause. The
observed repeatability is below 0.002 kcal/mol and is adequate for this
minimum-space screening; it must not be described as exact root-energy
invariance. No CAS enlargement or orbital-swap remedy is warranted by it.

All four center/return active-subspace comparisons have minimum singular value
greater than 0.9999999999999. The native orbital character remains consistent
with the selected fragment manifolds; no core intrusion is seen. Orbital labels
are Molpro C1, one-based:

| Channel | Closed / occupied | Natural active ordering |
|---|---|---|
| CH3S C–S | 10 / 14 | 11–13 sulfur 3p; 14 methyl radical |
| CH3SS C–S | 17 / 22 | 18–19 S2 pi; 20–21 pi-star; 22 methyl radical |
| CH3SS S–S | 17 / 22 | 18–19 CH3S sulfur pair; 20–22 atomic sulfur 3p |

The S–S canonical orbitals reorder the same active subspace: 18–20 atomic
sulfur, 21–22 CH3S pair. Rotations/reordering *within* the reviewed active
space are not active/inactive swapping. The original
[basis-specific orbital record](../examples/anl/orbital_reviews/molpro_confirmation/native_orbital_orderings.json)
and new audits supply the actual coefficients/occupations and checkpoint hashes.

Native eight-rank elapsed times for all roots at each recovery geometry were
24–34 s (CH3S aDZ), 67–89 s (CH3S aTZ), 204–212 s (S–S aDZ), and 533–693 s
(S–S aTZ). These are observed fixed-geometry job costs, not estimates of a
complete RotdPy Monte Carlo run.

## CABS: demonstrated completion, with a remaining interaction check

CABS=1 failed in the post-CI singles treatment; CABS=2 also failed after
converged MRCI, even with matching reference/internal-state counts. Their
messages identify `mrf12_singles.F90`, but do not expose the failed matrix or
the source-code guard. Claiming a specific linear-dependence defect or cluster
fault from these messages would be unsupported.

The documented coupled treatment, **CABS_SINGLES=3**, now completes all six
CH3S and twelve S–S roots. No native warnings/errors were detected. Every final
selected energy equals the printed relaxed-Davidson value of the coupled
calculation. Unlike the earlier separate-singles calculation, there is no
additional correction to add. The runtime selector remains `CLUSTER=0` and
the appropriate `ENERGD` array element; fixed/rotated Davidson alternatives
are not substituted. [Molpro CABS definitions](https://www.molpro.net/manual/doku.php?id=the_mrci_program#explicitly_correlated_mrci-f12).

At the CH3S asymptote, CABS=3 versus the completed CABS=1 calculation changes
the doublet–quartet gaps by less than 0.00002 kcal/mol. Spatial triplet spreads
are 0.0000314 kcal/mol. For S–S the spreads are below 0.000252 kcal/mol, while
the fragment centroid gap is approximately 1.235 kcal/mol.

MRCI still has nonzero asymptotic doublet–quartet residuals: approximately
−0.04009 kcal/mol for CH3S and −0.08079/−0.08071 for the S–S groups. Retain
these as method offsets in the documented, consistently dissociation-referenced
correction; they are not physical exchange interactions at infinite separation.
Do not combine CABS=1 far energies with CABS=3 interacting energies, or erase
the genuine fragment excitation gap by zeroing every electronic root separately.

## Root collapse: why MAXDAV20 alone was insufficient

Native file:
[CH3S tilted MAXDAV20 output](../local_runs/molpro_confirmation_recovery_review/ch3s_cs_atz_mrci_maxdav20/atz__tilted/atz__tilted.out).
Separate calls used `STATE,1,r` with the complete `REFSTATE` manifold. Increasing
MAXDAV improves the external iterative subspace, but does not require separate
calls to end at mutually orthogonal states.

| Requested independent calls | Final variational energies, Eh | Interpretation |
|---|---|---|
| Doublet 2 / 3 | −437.416361952422 / −437.416361952497 | Difference 7.5e−11 Eh; nearly identical dipoles |
| Quartet 1 / 3 | −437.400882693972 / −437.400882714922 | Difference 2.1e−8 Eh; nearly identical dipoles |

Doublet call 3 starts from a different internal state and drops by about
0.0148 Eh during the solve into the solution already reached by call 2.
This combined evidence is much stronger than near-equality of energies alone,
which can be legitimate degeneracy. A direct CI-vector overlap is not available
in this export, so the report does not claim such an overlap was measured.
The numerical completion marker cannot establish a distinct-state spectrum.

**Repair:** solve every declared root of a given spin simultaneously, using
`STATE,3,1,2,3` (or six for S–S), full matching `REFSTATE` and `NSTATI`,
and explicit `NOSWAP`. Molpro documents simultaneous states and warns that
separate excited-root solves can flip; including lower roots in the external
solve provides the protection absent from reference selection alone.
[Molpro STATE/REFSTATE](https://www.molpro.net/manual/doku.php?id=the_mrci_program#defining_state_numbers).
The multistate contraction framework is described in
[Knowles and Werner, Theor. Chim. Acta 84, 95–103 (1992)](https://doi.org/10.1007/BF01117405);
the executable syntax and indexing here follow the manual.

`NOSWAP` labels states by variational MRCI energy order. `ENERGD(i)` then selects
the relaxed-Davidson value corresponding to that calculated state, not the
ith separately initiated CASSCF root. The reconnaissance output explicitly
records this convention. It is not a claim of diabatic identity across crossings.
The new general JSON option is `correlation_options.state_policy="simultaneous"`;
it requires MRCI and `reference_policy="state_average"`. The default independent
policy and its `ENERGD(1)` selection remain unchanged. Numerical gradients
also select the matching array slot. Use a separate stage model for CASPT2,
which does not accept this MRCI-only policy.

## Minimal remaining submission and acceptance

[multistate_resolution_plan.json](../examples/anl/orbital_reviews/molpro_confirmation/multistate_resolution_plan.json)
contains seven points in three chains:

1. CH3S tilted, simultaneous roots and CABS=1: isolates root protection against
   the completed MAXDAV20 calculation at unchanged CABS and thresholds.
2. CH3S far/center/tilted, simultaneous roots and CABS=3: establishes a matched
   interaction curve and compares CABS=1/3 at the same tilted point.
3. CH3SS S–S far/center/tilted, simultaneous roots and CABS=3: finishes the
   previously blocked interacting references with protected states.

The far points are repeated because the simultaneous eigensolve is a new
native configuration; the completed CABS=3 gates provide its asymptotic
control. The CH3SS C–S one-root-per-spin checks and all accepted CASPT2 results
are reused. No old result is cancelled, deleted, overwritten, or silently
relabelled. Resources remain eight ranks, 1000 MW/rank and 104000 MB/node,
with all known exclusions and at most three simultaneous jobs.

Acceptance requires all requested simultaneous roots to complete, correct
relaxed-Davidson indexing, a distinct orthogonal eigenspectrum rather than
duplicate independently followed solutions, correct asymptotic grouping,
and matched far-referenced spin-gap/interaction comparisons. Review native
root coefficients and spectra, not just the job counts. The local tests
validate generation, parsing contracts and preservation; they cannot prove
that the installed Molpro executable supports this combination without error.
Commands and current deliverables are in the [handoff](../MOLPRO_CONFIRMATION_HANDOFF.md).

## Completed CH3S simultaneous-root review (October 16, afternoon)

The partial archive `multistate_ch3s_completed.evidence.tar.gz` has SHA-256
`6b31941179bd6c18763ef455898c3f4872dd0f71d76fbd072fd620467743897c`.
It contains the two completed CH3S chains, not the still-running S–S chain.
All **90 archived files** match the extracted bytes. This partial tar has no
producer inventory; the audit does not pretend it is a full `collect` export.
Each of the four jobs additionally passes manifest-bound input/output/WFU/XML
hashes, execution/native-input hashes, successful native termination, six
CASSCF and six correlated energy markers, and the declared simultaneous
three-root manifolds. There are no detected native warnings or errors.
All **24 selected energies** pass the independent relaxed-Davidson/CABS check.

The [machine-readable review](../examples/anl/orbital_reviews/molpro_confirmation/native_multistate_ch3s_review.json)
includes the native paths/hashes, energies, fixed-reference coefficient
matrices, printed state overlaps, orbital-space comparison and timing.
The [audit source](../examples/anl/orbital_reviews/molpro_confirmation/review_completed_chains.py)
is general to completed simultaneous-MRCI chains. Reproduce it with:

```bash
.venv/bin/python examples/anl/orbital_reviews/molpro_confirmation/review_completed_chains.py \
  local_runs/molpro_multistate_ch3s_review \
  examples/anl/orbital_reviews/multistate_ch3s_completed.evidence.tar.gz \
  --output examples/anl/orbital_reviews/molpro_confirmation/native_multistate_ch3s_review.json
```

**Root-collapse resolution.** At unchanged CABS=1 and convergence thresholds,
the tilted simultaneous calculation recovers three distinct states per spin.
For example, its third doublet variational energy is −437.397807322399 Eh,
instead of returning a duplicate of the lowest state at −437.416361952733 Eh.
The simultaneous CABS=3 calculation likewise has three distinct tilted states
per spin. Its maximum printed off-diagonal state overlap is 3.58e−6;
the fixed-reference projection matrices have smallest singular values above
0.9448. These are independent reference projections, not normalized full CI
overlaps. The explicit native overlaps and distinct reference mixtures provide
the complementary evidence. At the symmetry-axis center, the doublet pair
and quartet pair remain genuinely degenerate with independent reference
projections. Equal energies alone are not classified as collapse.

**CABS comparison.** The tilted CABS=1/3 controls have identical printed CASSCF
energies and the same active subspace (minimum singular value
0.99999999999997). CABS=3 changes individual relaxed-Davidson selected energies
by at most **0.002090 kcal/mol**, and the lowest doublet–quartet gap by
**0.001054 kcal/mol**. This supports coupled CABS=3 for this CH3S check without
enlarging CAS(5,4). The old independent CABS=3 far gate and new simultaneous
far point agree within 3e−8 Eh (0.000019 kcal/mol); their asymptotic spin
offset is unchanged. No failed or collapsed result is used as a reference.

| CABS=3 point | Lowest-doublet interaction | Lowest-quartet interaction | Lowest D−Q gap |
|---|---:|---:|---:|
| far | 0 | 0 | −0.040092 |
| center | −9.444162 | +2.149019 | −11.633273 |
| tilted | −7.770161 | +2.480771 | −10.291024 |

All entries are kcal/mol. Each interaction uses its own spin manifold's
lowest far energy; this diagnostic does not erase the far spin offset from
the raw gap, nor define final rate weights. The matched center-to-tilted
change is 1.674001 kcal/mol on the lowest doublet and 0.331752 on the lowest
quartet. These three points support the angular check, not an exhaustive PES.

**State labels and orbital ordering.** CH3S retains inactive orbitals 1–10 and
active natural orbitals **11.1–14.1**, five electrons in four orbitals. At the
tilted point, the lowest MRCI doublet has its largest fixed-reference
coefficient on reference 2; the lowest quartet on reference 3. `NOSWAP`
correctly labels these correlated states by variational energy order. This
is correlated-state reordering, not evidence of a changed active subspace.
Production installation must retain the method-specific root convention;
equal integer labels from different methods are not proof of equal state
character. In particular, the existing shared-reference CASPT2 center has
its lowest quartet in reference root 3. The lowest-per-spin diagnostic above
must not be implemented by blindly selecting root 1 from every method.

**Historical progress before the center failure.** The supplied S–S tail advances to iteration 6 of
the six-doublet simultaneous solve. Maximum displayed VAR(P) falls from
1.9e−7 to 9.7e−9, and energy changes from about 1.7e−5 to 9.5e−7 Eh.
Its 25,095,380 contracted configurations explain substantially greater cost
than CH3S. It is not yet converged to the requested 1e−10 residual threshold.
The quartet solve and dependent center/tilted points still have to finish.
The two waiting jobs are checkpoint dependencies, not additional failures.
That progress report was superseded when far completed and center failed;
see the current decision above. All original native inputs/results are retained.
