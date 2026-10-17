# Evidence for the three minimum active-orbital selections

The choices are CH3S C–S CAS(5,4), CH3SS C–S CAS(7,5), and CH3SS S–S
CAS(7,5). Selection combines fragment chemistry, occupations in larger parent
spaces, orbital coefficients, expansion comparisons at CASPT2, and path/return
checks. No occupation cutoff alone determines these choices. The S–S
simultaneous-MRCI interacting-point retry and its native state audit remain
pending; this document does not declare production complete.

## What the occupation numbers mean

The numbers below come from Molpro's explicit `NATURAL ORBITALS (STATE
AVERAGED)` tables in record `2140.2`. They are eigenvalues of the
state-averaged CASSCF one-particle density, between zero and two, not CASPT2
state populations. Equal averaging over an atomic sulfur 3p4 manifold gives
4/3 electrons per p orbital. Averaging the CH3S three-electron sulfur pair
over its two electronic components gives approximately 1.5 per orbital.
Fractional values from that averaging do not, by themselves, quantify the
correlation error or prove that individual states are degenerate.

Closed orbitals have occupation two by construction in the selected CAS.
Their resulting occupation is not independent evidence that closing them was
appropriate. The larger parent calculations and expansion controls supply
that evidence. Similarly, the printed AO coefficients identify dominant
atomic/angular character and relative phases; their squares are not population
percentages in a nonorthogonal AO basis.

## CH3S C–S: 11.1–14.1, CAS(5,4)

In the parent CAS(13,8) far calculation:

| Parent natural orbitals | Occupations | Character and decision |
|---|---|---|
| 7, 8, 10 | 2.00000 each | Methyl C–H bonding combinations; close |
| 9 | 2.00000 | Sulfur 3s; close |
| 11–13 | 1.33333 each | Complete atomic S 3p manifold; retain all three |
| 14 | 1.00000 | Methyl radical orbital; retain |

The parent printed coefficients identify 11, 12 and 13 as sulfur pz, px
and py, respectively. The methyl radical is carbon-dominated. The final
symmetrized aTZ confirmation mixes/rotates the sulfur p functions but keeps
the same three-dimensional sulfur subspace and methyl orbital. Thus the
chemical electron count is S 3p4 plus methyl p1 = five electrons in four
orbitals. A single bond/bond-antibond pair would omit sulfur components needed
for angular sampling and the three-root-per-spin asymptote.

Confirmed far occupations are `(1.33333, 1.33333, 1.33333, 1.00000)` in both
aDZ and aTZ. At the aTZ center they become
`(1.41637, 1.33331, 1.33331, 0.91701)` and at the inner point
`(1.49748, 1.49748, 1.47637, 0.52867)`. The redistribution is consistent with
bond formation in the retained subspace; it does not require preserving the
identity of each individual natural orbital within that subspace.

Adding S 3s, CAS(7,5), changes the tested lowest-doublet angular residual by
0.00039 kcal/mol and the angular D–Q gap residual by 0.01640 kcal/mol. The
larger (13,8) comparison gives 0.10803 and 0.03367 kcal/mol, respectively.
These are relative-potential diagnostics, not bounds on absolute accuracy.

## CH3SS C–S: reordered 18.1–22.1, CAS(7,5)

The larger parent CAS(19,12) shows why choosing the five highest occupied
orbitals, or simply applying an occupation threshold, is inadequate:

| Original parent indices | Occupations | Character and decision |
|---|---|---|
| 12–14 | 2.00000 each | Methyl C–H bonding combinations; close |
| 15, 16 | 1.99871, 1.99615 | Mainly sulfur 3s combinations; close |
| 17, 18 | 1.94578, 1.94578 | S2 pi bonding pair; retain |
| 19 | 1.94564 | S2 sigma bond; close in the minimum model |
| 20, 21 | 1.05334, 1.05334 | S2 pi antibonding pair; retain |
| 22 | 1.00000 | Methyl radical; retain |
| 23 | 0.06127 | S2 sigma antibonding partner; exclude from minimum CAS |

Parent 17 and 19 have almost identical occupations but different chemical
roles. The atomic coefficients distinguish the transverse pi manifold from
the sigma pair along the S–S bond. For example, parent 17 is the same-phase
sulfur pz combination; the other pi component is 18. The antibonding pi
partners 20/21 have opposite coefficients on the two sulfur atoms.

The desired original-parent set is therefore **17, 18, 20, 21, 22**. Swapping
17 and 19 places it contiguously at **18–22**, allowing `CLOSED,17; OCC,22`.
This initial choice is deliberate preparation, not an unintended swap during
a PES scan. The completed ordered checkpoints already contain it: applying
the original interchange again would corrupt the intended selection.

After minimum-space optimization, the aTZ far occupations of final indices
18–22 are `(1.95828, 1.95828, 1.04172, 1.04172, 1.00000)`. The corresponding
aDZ vector is `(1.95772, 1.95772, 1.04228, 1.04228, 1.00000)`. Both describe
six electrons in the complete S2 pi/pi-star manifold plus one methyl electron.
The aTZ native coefficients show equal-phase S contributions for 18/19,
opposite-phase contributions for 20/21, and carbon p character for 22.

The excluded sigma occupation 1.94564 and sigma-star occupation 0.06127
show real correlation; they are not exactly 2/0 in the parent. Closing them is
supported by the explicit (9,7) expansion comparison, which changes the tested
doublet angular residual by 0.00239 kcal/mol and the D–Q angular-gap residual
by 0.00095 kcal/mol. The S2 sigma bond remains intact in this C–S dissociation.
At the interacting aTZ center the active occupations become
`(1.95843, 1.95798, 1.13442, 1.04201, 0.90716)` and at the inner point
`(1.97064, 1.93072, 1.50723, 1.06898, 0.52243)`.

## CH3SS S–S: 18.1–22.1, CAS(7,5)

This channel uses a different seed and chemical subspace from CH3SS C–S,
despite the identical final index range. Its parent CAS(19,12) gives:

| Parent natural orbitals | Occupations | Character and decision |
|---|---|---|
| 12 | 2.00000 | Departing atomic sulfur 3s; close |
| 13–15 | 1.99981, 1.99945, 1.99943 | CH3S C–H bonding manifold; close |
| 16 | 1.99794 | Bound sulfur 3s; close |
| 17 | 1.97266 | C–S sigma bond in CH3S; close in minimum model |
| 18, 19 | 1.50018, 1.50017 | CH3S sulfur nonbonding/radical pair; retain |
| 20–22 | 1.33333 each | Complete departing atomic S 3p manifold; retain |
| 23 | 0.03036 | C–S sigma antibonding partner; exclude from minimum CAS |

The printed coefficients place 18/19 mainly on bound sulfur atom 2 and
20–22 on departing sulfur atom 3. Their seven electrons are three on the
CH3S pair and four on atomic sulfur. Both CH3S electronic components and all
three sulfur components are retained: six spatial roots per coupled spin,
not six degenerate states of atomic sulfur alone.

The reason for retaining the bound-sulfur pair is the two-component CH3S
ground electronic manifold, not a rule to include every lone pair. At C3v
the ground term is 2E; its spin-free components split when the symmetry is
lowered. The two sulfur-dominated orbitals support the complementary
three-electron arrangements `(2,1)` and `(1,2)`, and their mixtures. Calling
one a lone pair and the other the radical orbital describes one component,
not a permanent chemical distinction between them. Equal component averaging
explains the approximately `(1.5,1.5)` occupations. The 2E/Jahn–Teller framework
is supported by the primary study
[Jahn–Teller Effect versus Spin–Orbit Coupling in X2E CH3S](https://doi.org/10.1021/jp040443f);
assigning the native indices is based on our retained coefficients.

Freezing one of these two orbitals doubly occupied would remove the other
hole configuration from a common CAS reference. A single-component model
could describe a restricted geometry/state, but would not retain the complete
fragment manifold used for these angular and state-averaged calculations.
The deeper sulfur 3s lone pair is a different orbital and remains inactive.
This is the chemical reason for CAS(7,5) rather than a CAS(5,4) made from
atomic S 3p plus only one CH3S radical orbital.

The confirmed far occupations are `(1.50000, 1.50000, 1.33333, 1.33333,
1.33333)` in both bases. At the aTZ center they are `(1.52827, 1.50110,
1.33332, 1.33221, 1.30511)`; at the inner point `(1.71816, 1.63491,
1.49536, 1.36269, 0.78888)`. Canonical record `2141.2` reverses the order
of the fragment groups within the same CAS: atomic S occupies 18–20 and
the CH3S pair 21–22. This is internal reordering, not active-space loss.

Adding the C–S sigma/sigma-star pair, CAS(9,7), changes the tested doublet
angular residual by 0.04494 kcal/mol and the angular D–Q gap residual by
0.05703 kcal/mol. At 2.1 Å the interaction and D–Q differences reach 0.114
and 0.316 kcal/mol, respectively. Thus the evidence does not establish
uniform sub-0.1 kcal/mol agreement over every diagnostic quantity. It
supports the minimum space for the specified corrected-VRC construction;
absolute equality to the larger CAS is not its selection criterion.

For the user's factor-of-two kinetics target, the present evidence does not
require promotion to CAS(9,7). Retain (7,5) pending the native MRCI audit and
production rate checks. As an energy-scale comparison only, a uniform change
in activation free energy corresponding to a factor of two is `RT ln(2)`:
0.413 kcal/mol at 300 K, 1.377 at 1000 K, and 2.755 at 2000 K. The observed
angular residuals are smaller than these scales, but they are not measured
activation-free-energy errors. VRC rates depend on angular phase space,
variational bottlenecks and electronic-state treatment, and the full workflow
has additional thermochemical/master-equation uncertainties. No factor-of-two
rate bound follows from the sparse electronic comparisons, and choosing
(9,7) would not establish such a bound either. If a production sensitivity
check is needed, compare the corrected potentials/rates in the geometries
that control the flux before deciding whether a larger sampling CAS is useful.

## Continuity, states, and reproducibility

The minimum-space local forward/return all-root differences are 0.00000986,
0.00001097 and 0.00008911 kcal/mol for CH3S C–S, CH3SS C–S and CH3SS S–S.
Corresponding minimum active-span singular values are 0.9999999999986,
0.9999999999957 and 0.9999999998799. Native aDZ/aTZ return comparisons also
preserve the active spans; the four shared-reference recovery comparisons
have minimum singular values above 0.9999999999999. Native comparison uses
full MO expansions at identical geometry/basis and assumes native MO
normalization; it is not a cross-geometry CI-state overlap.

These checks support no detected active/inactive exchange or core intrusion
at the sampled points. They do not establish fixed state identity from a
root number. The CH3S independent MRCI calls previously collapsed; the
simultaneous calculations recover distinct tilted states and independent
reference projections, with maximum printed off-diagonal overlap 3.58e-6.
Correlated state reordering remains physical: the tilted lowest doublet and
quartet have dominant reference components 2 and 3. S–S requires the same
native audit after the retry, including its already-completed far result.

The quantitative expansion comparisons and their source outputs are in the
[minimum-space decision record](SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md).
State, selector and native continuity evidence is in the
[recovery review](MOLPRO_RECOVERY_REVIEW_2026-10-16.md). Final orbital labels
must be paired with the basis, geometry, record and WFU hash in
[native_orbital_orderings.json](../molpro_confirmation/native_orbital_orderings.json).

The [occupation evidence](../molpro_confirmation/native_occupation_evidence.json)
retains the three full parent active tables and all 30 completed CASPT2
confirmation/recovery occupation tables, with coefficients and file hashes.
The 30 output hashes were rechecked and their natural tables reparsed for
this explanation; all matched their retained audits and summed to the
requested active-electron counts within printed rounding.

The original parent outputs are retained in
`examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/ch3s_ch3ss_orbital_review_20261015.tar.gz`:

- CH3S: `02_recon/avt_f12__cs_r_30/avt_f12__cs_r_30.out`.
- CH3SS C–S: `05_recon/avt_f12__cs_r_30/avt_f12__cs_r_30.out`.
- CH3SS S–S: `08_recon/avt_f12__ss_r_30/avt_f12__ss_r_30.out`.

The new JSON records the archive and member hashes. Representative direct
native confirmation sources are the
[CH3S aTZ far table](../../../../../../local_runs/molpro_confirmation_recovery_review/ch3s_cs_atz_caspt2_shared_reference/atz__far/atz__far.out),
[CH3SS C–S aTZ far table](../../../../../../local_runs/molpro_confirmation_r2_review/ch3ss_cs_atz_caspt2/atz__far/atz__far.out),
and [CH3SS S–S aTZ far table](../../../../../../local_runs/molpro_confirmation_recovery_review/ch3ss_ss_atz_caspt2_shared_reference/atz__far/atz__far.out).
