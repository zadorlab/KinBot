# Soft motions and the next R8 characterization pass

This is an evidence record and restart scope, not an automatic well-merging
decision. The narrow sulfur-axis repair and restored harmonic midpoint optical
heuristic do not establish a lowest-energy path or an anharmonic partition
function. No new Gaussian calculations were run locally for this audit.

## Native evidence retained

The selected L2 database rows were compared with normally completed Gaussian
logs, including their full printed displacement vectors. Their geometries and
frequencies match the selected rows. The reproducible read-only script is
`examples/anl/ch3ss_anl1_f12/tools/audit_soft_motions.py`; the regression fixture
`tests/reference/soft_motion_native.json` retains the vectors, job names, row
IDs and native log SHA256 values. Local full report:
`local_runs/soft_motion_review_20261103/native_modes.json`.

| Saved MESS label | Worker | Evidence |
|---|---|---|
| w3 | 791671450790040000002 | L2 minimum, comparison case |
| w7 | 791671560700020000002-sedd8f7729b69b19b | L2 minimum |
| w14 | 791671560700020000002-s893815ff03c7331d | L2 minimum |
| w15 | 791671560700020000002-sdc19d1516ced03fc | L2 minimum |
| rxn22 | w7 R_Addition_MultipleBond_2_1_3 | L2 saddle connecting w7/w15 in saved network |
| rxn29 | w14 12_shift_S_R_2_1_3 | L2 saddle connecting w14/w7 in saved network |

The rounded saved MESS energies are w7 26.03, w14 26.68, w15 25.69,
rxn22 26.41 and rxn29 26.90 kcal/mol relative to initial CH3SS. These give
direction-dependent barriers of 0.38/0.72 and 0.22/0.87 kcal/mol. These are
MESS E0 differences; compare a HIR electronic barrier with electronic TS-minus-
minimum energies at the same theory before attributing a difference to path
choice. Keep the unrounded native quantities for the actual comparison.

The two saddle imaginary modes have approximately 98.7% and 96.9% squared
projection onto the combined local torsional subspace. This describes the
motion near those particular saddles; it does not exclude a lower bending
route. Minima also have soft predominantly nontorsional modes: w7 257.0,
w14 174.0 and 260.9, and w15 211.7 and 269.4 cm^-1. The w3 modes at 188.4,
256.0 and 304.8 cm^-1 have roughly 91%, 25% and 76% torsional character.
Residual motion is not automatically classified as an in-plane bend. Printed
vectors have limited precision; these percentages are character diagnostics.
Neither soft frequencies nor failure to discover another w3 minimum proves
that an analogous path exists or is absent.

For a Cartesian displacement d, first form u=M^(1/2)d and project out the
mass-weighted translations and rotations. Normalize the remaining u. Construct
physical infinitesimal top rotations about each eligible bridge axis, remove
the external motions, and orthonormalize their combined span by SVD to obtain T.
The reported torsional fraction is ||T u||^2; the residual is 1-||T u||^2.
Individual-axis overlaps are not additive because axes can couple. This is
not an overlap integral between vibrational wavefunctions, a finite path, or
a replacement for measured potentials and their kinetic metric.

## Restart scope and subsequent decision

Use the explicit `--reassess-thermochemistry` controller-refresh option for
both PESs. It keeps the cutoffs, accepted QC, attempt histories and live jobs,
and reopens characterization of every retained worker and accepted reaction.
Thus it includes w3/w7/w14/w15, their connecting and outgoing reactions, other
eligible wells/TSs/products, and the already staged reverse searches. Existing
homolysis checks resume under their existing settings. A naturally rejected
search is not retried merely because the rotor rule changed. Missing eligible
HIR/conformer characterization can run; valid compatible points are adopted.
The journal names the exact worker inputs and accepted reactions reassessed.

After native completion, compare measured geometries and paths, not labels:

1. Determine which scans remain in the declared population and which reach
   another minimum. Retain the native observations even when a scan is unusable
   as that population's HIR model.
2. Compare any alternate HIR connection with the saved TS/IRC route on a matched
   energy basis. A higher torsional path can coexist with a lower bend/wag path;
   it must not replace the lower path by assumption.
3. Check HIR/harmonic mode replacement, the optical midpoint decision, and
   population counting together. Do not count the same accessible domain both
   as separate wells and as an additional full-domain rotor.
4. Merge populations only if a consistent state-density representation is
   supported. Preserve/remap every distinct entrance and exit channel; an
   internal interconversion already represented in that state density cannot
   also be counted as independent flux.

The restart is a native test of the implemented repairs. It does not guarantee
that independent one-dimensional HIR scans will resolve coupled anharmonicity,
nor does it implement such a merger. The existing reverse-search integration
and unresolved homolysis export gates remain in force: controller completion
alone does not certify a final MESS model. Further model work, if needed, is
determined from those native results rather than presumed here.
