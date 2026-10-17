# Native review: CASPT2-F12 real shift 0.1 Eh

Both matched CH3SS S–S CAS(7,5) calculations at 2.1 Å passed independent
review. Use IPEA 0.25 and real shift 0.1 Eh for the production CASPT2 route.
The minimum active spaces remain (5,4), (7,5), (7,5). This closes the requested
single-geometry numerical check, not the full kinetic validation.

## Native evidence

The deposited `minimum_active_space_molpro_confirmation_20261016_shift01.evidence.tar.gz`
contains 58 inventoried files; every size and SHA-256 matched. The previous
recovery export's 539 inventoried files were also reverified. Jobs 52908250
(aDZ-F12) and 52908251 (aTZ-F12) terminated normally with exit zero, all 12
correlated energies each, all required completion markers, and no native
warnings or convergence errors. Both use six doublet and six quartet roots,
IPEA 0.25, and six equal-weight same-spin references in each CASPT2 call.

The [machine-readable comparison](../examples/anl/orbital_reviews/molpro_confirmation/native_shift01_review.json)
records individual energies, reference energies, iteration counts, squared
norms, exact native-output/input hashes and orbital-overlap matrices. Native
outputs are retained under
`local_runs/molpro_confirmation_shift01_review/ch3ss_ss_{adz,atz}_caspt2_shift01/{adz,atz}__inner/`.
The comparison outputs are in the corresponding `local_runs/molpro_confirmation_recovery_review`
shared-reference chains. Both tar archives remain in `examples/anl/orbital_reviews`.

## What changed and what did not

The effective decks differ in real shift, title, and removal of `DIAB,2140.2`
by the recovery adapter. They are therefore not byte-identical except for
SHIFT. Initial checkpoint hash, geometry, basis, active space, SA weights,
roots, IPEA, reference policy and numerical tolerances match. The removed
DIAB directive has no measurable effect on these CASSCF results:

- Every printed CASSCF and correlated-reference energy is unchanged.
- The active natural-orbital overlap matrices are identity to numerical
  precision. Minimum singular values are 0.9999999999999987 (aDZ) and
  0.9999999999999714 (aTZ), using the checked same-geometry MO expansion.
- Within each spin manifold, correlated energy ordering is unchanged.
  This establishes no collapse or orbital substitution in this matched check;
  it is not a proof about every uncomputed sampling geometry.

## Relevant energy sensitivity

The production correction uses spin gaps, not absolute CASPT2 totals.

| Quantity, kcal/mol | aDZ-F12 | aTZ-F12 |
| --- | ---: | ---: |
| Lowest D–Q gap, shift 0.3 | -76.231721 | -76.929649 |
| Lowest D–Q gap, shift 0.1 | -75.976268 | -76.665211 |
| Change, 0.1 minus 0.3 | +0.255453 | +0.264439 |

The requested power-3.7 estimate `(3**3.7*gap_T - 2**3.7*gap_D)/(3**3.7-2**3.7)` changes by
**+0.267019 kcal/mol**. It is a two-basis estimate, not a calibrated
uncertainty bound for CASPT2-F12 gaps. Changes in each doublet's gap to the
common quartet-root-1 reference are retained individually in the JSON.
Absolute correlated energies decrease by 3.37–3.86 kcal/mol at aDZ and
3.57–4.09 at aTZ; those common shifts largely cancel in the gaps, but must
not be ignored when constructing dissociation zeros. Production must calculate
both reference asymptotes and interacting points at shift 0.1. No old
shift-0.3 interaction energy is relabelled or mixed into the new correction.

The largest final shifted squared norm increases from 1.151693 to 1.183351
(aDZ), and from 1.166607 to 1.200101 (aTZ). These are native `SQ.NORM` values,
including the zeroth-order norm. They are finite and converge, with no native
intruder/negative-norm warning. They do not alone establish perturbative
accuracy. The 0.267 kcal/mol gap sensitivity is small compared with this
short-distance spin splitting but is not zero and does not establish a
factor-of-two uncertainty in the final rate.

## Reproduce

From the repository, after extracting both collected archives into the
retained local review directories:

```bash
.venv/bin/python examples/anl/orbital_reviews/molpro_confirmation/review_shift_check.py \
  local_runs/molpro_confirmation_recovery_review \
  local_runs/molpro_confirmation_shift01_review \
  --output examples/anl/orbital_reviews/molpro_confirmation/native_shift01_review.json
```

No Molpro calculation is repeated. The review verifies all inventoried native
files before reading energies; unexpected native-input changes or inconsistent
root/reference/orbital evidence stop the comparison.
