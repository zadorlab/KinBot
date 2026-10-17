# CH3SS C–S and S–S: distinct CAS(7,5) subspaces

The main input is [../ch3ss.json](../ch3ss.json). It contains every editable
channel, active-electron/orbital count, basis, checkpoint hash/record/index,
state/root, weight and solver choice. The WFU files in [orbitals](orbitals/)
are native data, not another configuration. [orbital_inventory.json](orbital_inventory.json)
records their original source and hashes; production does not read that inventory.

The recommendation is the smallest chemically complete subspace supported by
relative interaction energies and spin gaps, angular/radial and return checks,
and expansion comparisons. Absolute full-valence total-energy error is not
the selection criterion. The evidence supports corrected VRC sampling, not a
factor-of-two rate guarantee. MRCI was closed by the user and is not a
production prerequisite. CASPT2 uses IPEA 0.25 and real shift 0.1 Eh; shared
same-spin references and the lowest-quartet correction selector remain explicit.

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
require promotion to CAS(9,7). Retain (7,5) for the CASPT2-only production run; the user closed
the optional MRCI audit. Production rate checks remain outstanding. As an energy-scale comparison only, a uniform change
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

## Retained evidence and interpretation

- [Occupation and AO-character records](occupation_evidence.json): parent-space
  occupations and confirmed points specific to this example. These are
  state-averaged natural occupations, not CASPT2 populations. Degenerate active
  orbitals may rotate internally without changing the active subspace.
- [Native CASPT2 inputs, outputs, natural-orbital XML and execution records](caspt2_native.evidence.tar.gz),
  copied byte-for-byte from the original exports; every member and source
  archive hash is listed in [native_inventory.json](native_inventory.json).
  Earlier confirmation points use shift 0.3; the CH3SS archive also includes
  the two matched shift-0.1 points. Production recomputes all energies and
  asymptotes at shift 0.1 rather than mixing these historical energies.
- [Joint literature/selection rationale](joint_history/docs/SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md)
  and [native shared-reference/return review](joint_history/docs/MOLPRO_RECOVERY_REVIEW_2026-10-16.md).
  The historical combined campaign is stored once and is not a runtime dependency.
- [Matched shift-0.1 review](joint_history/docs/CASPT2_SHIFT01_NATIVE_REVIEW_2026-10-16.md).
  Both aDZ/aTZ S–S checks completed all 24 energies; the exponent-3.7
  extrapolated spin-gap change was 0.267019 kcal/mol. This is a sensitive-point
  check, not an uncertainty estimate for all production geometries.

The tested paths showed no detected loss of the chosen active span or core
intrusion. They do not prove universal state tracking at uncomputed geometries.
Correlated quartet ordering does change: the input evaluates its complete
reviewed manifold and selects the lowest reference energy. Explicit doublet
roots still define the target surfaces and their geometry derivatives.
