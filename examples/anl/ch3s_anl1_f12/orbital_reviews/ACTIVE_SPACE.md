# CH3S C–S: CAS(5,4)

The main input is [../ch3s.json](../ch3s.json). It contains every editable
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
- [Joint literature/selection rationale](../../ch3ss_anl1_f12/orbital_reviews/joint_history/docs/SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md)
  and [native shared-reference/return review](../../ch3ss_anl1_f12/orbital_reviews/joint_history/docs/MOLPRO_RECOVERY_REVIEW_2026-10-16.md).
  The historical combined campaign is stored once and is not a runtime dependency.
- [Matched shift-0.1 review](../../ch3ss_anl1_f12/orbital_reviews/joint_history/docs/CASPT2_SHIFT01_NATIVE_REVIEW_2026-10-16.md).
  Both aDZ/aTZ S–S checks completed all 24 energies; the exponent-3.7
  extrapolated spin-gap change was 0.267019 kcal/mol. This is a sensitive-point
  check, not an uncertainty estimate for all production geometries.

The tested paths showed no detected loss of the chosen active span or core
intrusion. They do not prove universal state tracking at uncomputed geometries.
Correlated quartet ordering does change: the input evaluates its complete
reviewed manifold and selects the lowest reference energy. Explicit doublet
roots still define the target surfaces and their geometry derivatives.
