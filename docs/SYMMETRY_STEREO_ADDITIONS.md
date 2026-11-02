# Symmetry and stereochemistry: additions in this review

Scope: additions/corrections published in `776fecae1`, building on KinBot's
existing symmetry and stereochemistry machinery. This is not a claim that the
underlying concepts or original implementations were new to this review.

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
  Remove the midpoint harmonic-energy/4 kcal/mol heuristic as a counting
  decision, while retaining its diagnostic. This corrects/extends the
  rotational/optical-weight work recorded under LukaDockx, `f0a1811b1` (2026).
- **Evidence-based mirror-saddle flux accounting.** Calculate signed parity
  of the unstable native mode under each admissible mirror mapping. Share a
  single even self-mirror harmonic surface between its two represented mirror
  connections, conserving its total statistical weight. Missing evidence,
  mixed parity and PST/ROTD surfaces do not receive this rule. It does not
  calculate trajectories or dynamical branching.
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
