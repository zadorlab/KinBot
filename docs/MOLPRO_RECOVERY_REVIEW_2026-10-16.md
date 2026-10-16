# Native recovery review and MRCI root resolution

## Decision

The minimum active spaces remain **CH3S C–S (5,4), CH3SS C–S (7,5), and
CH3SS S–S (7,5)**. All 20 shared-reference CASPT2 recovery points and both
CABS=3 asymptote tests completed and passed independent native completion,
energy-selector and provenance checks. The CH3S tilted MAXDAV20 MRCI point
also completed numerically, but its separately requested roots contain strong
evidence of collapse. It is not an accepted distinct-state spectrum.

The next submission is seven MRCI-only points. It changes correlated root
protection and completes matched coupled-CABS interaction checks. It does not
repeat CASPT2, the active-space survey, IPEA tests, or optimizations. Full
KinBot production awaits the resulting MRCI state/reference acceptance and
installation of the reviewed models.

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
