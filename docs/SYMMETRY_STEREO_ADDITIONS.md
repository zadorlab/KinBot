# Symmetry and stereochemistry: additions in this review

Scope: additions/corrections published in `776fecae1`, with the subsequent
restoration of the established harmonic midpoint approximation noted below.
These build on KinBot's existing symmetry and stereochemistry machinery; the
underlying concepts and original implementations are not claimed as new here.

- **Actual vibrational-mode comparison.** Read native Gaussian Cartesian
  displacement vectors, bind them to their own calculation's geometry, and
  compare mass-weighted modes using the same atom permutation and proper
  rotation as the geometry. Replaces the older geometric `matrix_corr` proxy
  on the L1/L2 saddle-comparison path. The earlier comparison implementation
  is recorded under Judit Zador, `2f4f62810` (2019).
- **Geometry-bound mode transfer.** Use calculated squared overlaps and a
  one-to-one assignment for VPT2 transfer, with ambiguous assignments and
  explicit isotope provenance checked. Reconstruct the L3 harmonic spectral
  model only from a complete, consistent native spectrum. Frequency proximity
  across calculations is not substituted for a calculated mode overlap.
- **Whole-structure external symmetry.** Test graph-permitted atom mappings
  against proper rotations of the entire geometry. Disconnected fragments
  cannot rotate independently to inflate the external symmetry number of one
  complex or transition state.
- **Consistent optical counting.** Share counting across MESS, MESSPF and
  thermochemistry; distinguish explicit mirror populations from mirror
  configurations already covered by conformer or hindered-rotor domains.
  Retain the established 4 kcal/mol harmonic-midpoint fallback for ambiguous
  puckering/mirror geometry, explicitly labeled as an approximation. Its
  removal in `776fecae1` was not a demonstrated bug fix and has been reversed.
  The shared counting extends the rotational/optical-weight work recorded
  under LukaDockx, `f0a1811b1` (2026); the midpoint approximation itself predates
  this review.
- **Evidence-based mirror-saddle flux accounting.** Calculate signed parity
  of the unstable native mode under each admissible mirror mapping. Share a
  single even self-mirror harmonic surface between its two represented mirror
  connections, conserving its total statistical weight. Missing evidence,
  mixed parity and PST/ROTD surfaces do not receive this rule. It does not
  calculate trajectories or dynamical branching.
  The rxn20 evidence supports this duplicate-surface correction, not a failure
  of the harmonic midpoint cutoff.
- **Consistent cached and network representations.** Refresh reconstructible
  legacy symmetry factors and stale geometry labels without changing raw QC;
  apply each correction once. Mirror connections use one consistent complete
  surface/energy reference. Mirror/conformer *work-selection* groups use exact
  first-passage union probabilities, without asserting physical lumpability.

Evidence: [math, native examples and limitations](SPATIAL_MODE_AND_MESS_REVIEW.md),
`tests/reference/mirror_even_transition.log`, and
`tests/reference/mess_symmetry_review.json`. The complete captured PES handoffs
were replayed twice with QC submission forbidden and original source hashes
unchanged. Native examples include the rxn20 parity (+0.9999996855), whole-body
sigma checks, and the twelve-point 360-degree rotor domain.

Git author records identify commit attribution, not whether a human or an
assistant authored individual lines. These additions do not replace the
existing connectivity, stereochemical identity, or stationary-point checks,
and do not establish coupled-rotor dynamics or electronic-wavefunction overlaps.

## Parallel transition states

Endpoint and stereochemical classes can contain several distinct saddles. PES
assembly, direct MESS output, PESViewer and pre-ANL selection retain these
parallel paths. MESS combines their complete statistical models in a `Union`,
including each member's energy and tunneling parameters.

Repeated discoveries are merged using an element/isotope-preserving atom map,
a proper rotation fitted to the geometry, and the mass-weighted imaginary-mode
overlap under that same transformation. The comparison uses 0.01 Å RMSD,
0.03 Å maximum atomic displacement, overlap at least 0.99, electronic-energy
agreement within 0.1 kcal/mol and frequency agreement within 5 cm⁻¹. Both
spectra must establish a first-order saddle. These are duplicate-comparison
tolerances. Accepted QC and stationary-point validation retain their own rules.

The first verified representative is retained, avoiding chains of approximate
equivalences. Mirror completion reflects both coordinates and displacements.
Full conformer-ensemble equivalence requires a one-to-one member comparison.
Missing mode or isotope-mass evidence keeps the observations with an explicit
unresolved-comparison warning. `saddle_selection.json` records the selected
names, native sources, mappings, measurements and decisions. Saved accepted
results supply this evidence during reassembly.

The native regression fixture `tests/reference/parallel_saddles.json` contains
two repeated observations of a loose saddle and a distinct compact saddle.
The latter survives selection in every discovery order, including reversed
IRC endpoint ordering.
