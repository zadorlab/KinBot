# Geometry, vibrational modes and MESS counting

This review changes shared KinBot routines. The CH3S/CH3SS export provides
regression evidence, not molecule-specific rules in installed code.

## What is actually calculated

For a chemically permitted atom mapping, center the two coordinate matrices
X and Y. With `X.T @ Y = U S V.T`, the proper Kabsch rotation for row vectors is
`R = U diag(1, 1, det(U @ V.T)) V.T`. Translation aligns the centroids.
The determinant constraint forbids an unnoticed reflection. Geometry RMSD
is `sqrt(sum_i |x_i R + t - y_i|^2 / N)`; checks that require every atom to
match also use the maximum atom error. Graph automorphisms only propose atom
permutations. They do not by themselves establish geometric equivalence.

Cartesian vibrational displacements transform with the geometry's R, with
the same atom mapping and **no translation**. Their normalized overlap is

```
C = sum_i m_i u_i . v_i / sqrt(sum_i m_i |u_i|^2 sum_i m_i |v_i|^2).
```

Independent eigenvectors have arbitrary signs, so mode identity uses `abs(C)`.
Mirror parity compares a mode with its own reflected/permuted/rotated copy
and retains the sign. Positive parity preserves the reaction coordinate;
negative parity reverses it. The vector arrays are never fitted independently
to maximize their agreement. `geometry.matrix_corr` is an older geometric
similarity function and is no longer used for the L1/L2 mode comparison.

Gaussian frequency outputs contain the vectors below each `Frequencies --`
line, in `Atom AN X Y Z` columns. The native reader checks ordered atoms, the
complete spectrum, normal completion, and the corresponding orientation.
Its geometry binding permits 1e-5 angstrom maximum atom error, accounting for
printed coordinate precision. A force-only output is not a mode source.
Printed vectors are finite-precision evidence; they are not an exact Hessian.

The L1/L2 saddle comparison reads each selected job's own geometry and mode.
It uses the existing 0.9 overlap/looser bond-geometry criterion, or the existing
tighter geometry-only criterion. Missing native vectors are recorded as
unavailable, not assigned a fictitious overlap or automatically retried.
Stationary-point order and stereochemical checks remain separate. Explicit
isotopes require matching mass provenance; the new mode path does not silently
substitute elemental masses for them.

VPT2 transfer rotates native modes into the accepted frame and uses squared
mass-weighted overlaps. Hungarian assignment maximizes the total one-to-one
overlap. It refuses low overlaps or ambiguous alternatives carrying materially
different shifts. A vector-free public helper now requires an explicit
reviewed mapping. Joining the harmonic and anharmonic tables **within one
native output** uses frequency consistency; unresolved degenerate assignments
with unequal shifts are rejected. Frequency proximity across different
calculations is not a substitute for mode vectors.

The L3 harmonic spectral reconstruction uses normalized mass-weighted native
modes and signed frequency squares, after checking their mass orthogonality.
Only printed roundoff is removed. The same established internal-rotation
projection is applied to that spectral model; accepted L2 geometries, raw
Hessians, rotor potentials and electronic energies are retained separately.

## Geometry, optical counting and dividing surfaces are separate

Rigid rotational symmetry uses proper rotations of the **whole assembly**.
Disconnected fragments cannot each rotate independently when assigning one
bound complex or saddle's external symmetry number. The surviving mappings
must form a group. The whole-geometry maximum atom tolerance is 0.05 angstrom.
Existing stable-molecule and internal-rotor conventions remain in place.

Optical counting distinguishes explicit mirror populations, distinct rigid
mirror geometries and mirror configurations already represented by a rotor
domain or conformer sum. These are domain/counting questions, not overlaps of
electronic wavefunctions. MESS, MESSPF and exported thermochemistry now use
the same parent counting routine instead of a stale `nopt` on one path.

The previous midpoint rule estimated `0.5 d.T H d` between a structure and an
aligned mirror, then used a 4 kcal/mol cutoff to change optical weight.
A local curvature at one structure cannot establish the intervening inversion
path or whether a torsional domain already covers it. That cutoff no longer
changes counting. The estimate remains a labeled diagnostic. Cases lacking
domain evidence remain explicitly unresolved; no new quantum overlap or
inversion barrier is claimed.

For a native self-mirror saddle, the unstable-mode parity adds information
that geometry alone lacks. The parity classifier uses a 0.01 angstrom
maximum atom error and signed overlaps separated at +/-0.99. All admissible
mirror mappings must agree. Its provenance contains geometry/mode hashes,
the mappings and overlaps, and the native output hash or selected Hessian
row/units. Missing or mixed-parity evidence does not authorize flux sharing.

When one even self-mirror harmonic surface is represented on two distinct
mirror connections, each copy receives half that surface's statistical
weight. Their summed flux is the original single-surface flux. This is not
an extra site degeneracy or a fitted trajectory branching fraction. Odd
parity, missing evidence, and PST/ROTD surfaces do not receive this rule.
The shared weight is attached to each individual RRHO member, once.

## Native R8 checks

* `rxn_20`: native L2 imaginary frequency -850.9602 cm-1. The maximum atom
  error for the mirror mapping is 0.002543411 angstrom; mass-weighted signed
  parity is +0.9999996855. This supports one even harmonic surface shared
  between the two represented mirror connections. It does not calculate
  post-transition-state trajectories or recrossing.
* `w_10` and `rxn_29`: whole-geometry external sigma is 1, not the former
  graph-only 2. With their retained optical division the factor is 0.5.
* `w_14`/`w_15`: actual proper C2 rotations support external sigma 2.
  Distinct mirror weights of two yield factor 1; the proposed C1 correction
  would be wrong for the exported geometries.
* `w_6`: twelve normally completed native rotor points cover 360 degrees.
  That represented rotor domain already covers its mirror; no extra factor
  two belongs in the partition function. This uses the declared independent
  rotor model, not an assumed coupled torsional dynamics calculation.
* `pr_6`/`pr_18`: the selected geometries are cyclic mirrors. `pr_11` is
  open-chain. Stale open-chain comments are refreshed from cached population
  geometries; they are not used to collapse these different products.

The exact native rxn20 excerpt is
`tests/reference/mirror_even_transition.log`, SHA-256
`f0423445a3f9e826947c7f805b9f1a4ee40de11bf3943a95b049a27726be20c2`.
Other selected native geometries and rotor rows are in
`tests/reference/mess_symmetry_review.json` with their original identifiers.

## Local replay and energy conventions

The local route is the owning `kinbot.pes ... no-kinbot` postprocessor on an
isolated export. It reuses the saved searches, databases and worker fragments;
it does not redo reaction discovery. Legacy single-RRHO symmetry corrections
are made in memory only when the old graph factor and optical division are
reconstructible. HIR/custom/multiple-member factors are not inferred from a
bare geometry. New fragments carry a version marker, preventing double
correction. Repeated complete replay must produce identical output.

Configured electronic levels and their declared scalar-reference correction
are kept together. Ordinary worker E0 is E + ZPE + that explicit correction;
accepted ANL E0 already includes its components and is not corrected again.
PES-wide L2/L3 energy assembly uses the same convention. Raw database energy
and ZPE remain unchanged. Cached level populations are refreshed against the
computed multiplicity and known charge; a different requested multiplicity
does not relabel existing QC. Reactive surfaces retain their own spin weight.
Only normally completed finite database results can supply PES energies.

The previously documented S2 state restoration is still a separate derived
model using native triplet QC. It is not accomplished by changing a singlet
label or degeneracy. The original archive remains intact.

The complete capture regression forbids QC/worker submission, runs both PES
handoffs twice, verifies original file hashes, checks the corrected native
factors and product labels, and exercises repeated preparation. Unit tests
also include the old independently-rotated-mode false positive, ambiguous
mode assignments, improper reflections, absent evidence and repeated energy
correction. Native MESS execution is recorded separately from these tests.

## Code history and limits

This review added the shared alignment/native-mode reader, calculated
mirror-mode parity, shared-surface counting and cached-model corrections.
`geometry.matrix_corr` predates this work (recorded author Judit Zador,
commit `2f4f62810`, 2019). The midpoint counting heuristic predates this pass
(`f0a1811b`, recorded author LukaDockx, September 2026). Earlier ANL thermal
and reuse commits record Hrishikesh Ram as author. Git author fields do not
establish whether a human or assistant wrote a particular line. The current
changes take responsibility for the paths they use, rather than relying on
authorship as validation.

RRHO/PST scouting, printed vector precision, independent hindered rotors and
statistical ergodicity remain model approximations. An ordinary Gaussian
termination or MESS exit alone cannot establish physical accuracy. Negative
pressure-dependent MESS coefficients are not clipped or used for pruning.
The collisionless/HPL scope and sensitivity of any shortlist must be stated.

## Reconciliation with the two independent review artifacts

The original artifacts were rechecked against the delivered inputs, rather
than treating every proposed change as an established bug.

| Flag | Resolution and scope |
| --- | --- |
| A: asymmetric mirror routes; unequal mirror PST thresholds | Shared mirror assembly applies one complete selected model consistently to both sides, including its threshold and endpoint references. Distinct classified pathways remain distinct. |
| B: split retention scores | Shared grouping and exact first-passage union probabilities retain or drop mirror/conformer work groups together. The kinetic wells remain separate. |
| C: merge rotamers and remove torsional saddles | Group them for work selection. Similar connectivity alone does not prove rapid equilibration or justify replacing separate kinetic populations by one hindered rotor. No such physical lumping was applied. |
| D: rotational/optical factors; rxn20 double counting | Whole-geometry rotations, documented rotor domains, and native reactive-mode parity resolve the audited cases. The proposed C1 assignment for w14/w15 was contradicted by their C2 geometries. |
| E: loose harmonic saddles | They are included in the 0.5/2 sensitivity scenarios. This is a scouting sensitivity check, not validation of a production VRC potential or proof that the RRHO error is bounded by two. |
| F: electronic levels | General writer, cached-population and energy-reference handling are native KinBot fixes. S, SH, CH3S and CH require explicit sourced models in the input; no molecule-specific level table is hidden in the code. |
| G: shallow w9/submerged rxn37 | The selected submerged-barrier approximation is explicit. A saddle below an optimized minimum at a different level does not by itself prove that minimum does not exist. Keep the full scout representation; collapsing it into one production channel needs its own model decision. |
| H: retained returning side branch | Exact branching includes returns. A visited branch need not produce an exit; retention is conservative and is not a recommendation to run every retained object at high level. Candidate scopes remain unapproved. |
| I: negative pressure rates | Unresolved numerical validation, excluded from all scope selection. Higher-precision native checks are separate from the KinBot fixes. |
| J: search-family names | Preserve these as provenance identifiers, not mechanistic assignments. Scope groups use endpoints and population identities, not the search-family label. |
| K: apparently identical product labels | Cached labels now follow the selected geometries: pr6/pr18 are cyclic mirrors, while pr11 is open-chain. They are not merged by their formerly stale comments. |
| L: CSH versus HCS product | The CH3S pr3 native fragment geometry has H bonded to S, as its H-S-C comment says. It is not replaced by the more stable connectivity merely because that would be expected. |

The follow-up input audit found that the first symmetry-review push retained
the previous CH omission. Both production and scout examples now explicitly
include CH X2Pi components at 0 and 28.14675 cm-1, with weight two each and a
single -A0/2 scalar-centroid correction. The value is from Colin and Bernath,
J. Mol. Spectrosc. 263 (2010), 120-122, Table 2,
https://doi.org/10.1016/j.jms.2010.06.013. This is the same effective electronic
component/RRHO approximation used for the scout's other diatomic terms; it
does not resolve spin-rotation mixing, Lambda doubling or hyperfine structure.
An explicit regression now checks coverage of all four flagged degenerate
species, cached CH rendering, unchanged raw geometry/frequency/energy fields,
and preservation of the reactive-surface spin weight.

## Completed release validation and updated scouting scope

The isolated published checkout plus only this review's files passed **538
tests, 97 subtests**, with one platform-specific skip. The exported pair was
replayed twice through the owning postprocessor with QC and worker launches
forbidden. Original inputs, databases, searches, screening decisions and
worker fragments remained byte-identical. This is local evidence, not a new
remote controller handoff.

The initial four full native MESS runs completed: both PESs with the symmetry
corrections, and both with the already declared electronic levels and
scalar-reference corrections restored together.
After the CH coverage follow-up, both full models and all six reduced
candidates were rerun successfully. Candidate membership is unchanged. The
follow-up passed 27 isolated tests and 17 subtests, including the new coverage
and cached-fragment regressions. The current full input SHA-256 values are:

* CH3S: `26e4703ea1bb0fd83baf08d33095e5ae4bf2a6b6718f066a2d4ec13cc3ffdcd9`
* CH3SS: `a336c86955e485c4125aa6a7e20090881b8c2e9086d6f80e3c5515ecf593ec08`

Both full models and all six reduced candidates also completed a native
1000 K / 760 torr direct-method check with energy-step/temperature 0.2.
Retained microcanonical state counts and off-diagonal HPL coefficients agree
exactly at the printed precision. This does not validate the pressure tables.

Work selection groups explicit mirrors and related population connectivity
families conservatively, while retaining separate wells in the kinetic
equations. The probability of using **any** channel in a group is computed as
one first-passage event, including returns; it is not a sum of individual
ever-hit probabilities. The independent 80-digit decimal solve checks the
baseline grouped probabilities.

The tested scenarios include internal energies 140, 148, 148.142152, 149 and
155 kcal/mol, native HPL rates over 300–2000 K, and factors 0.5/2 on individual
PST or loose-saddle groups. These are sensitivity tests, not rigorous error
bounds. All removals are then applied together and product branching compared
with the full model.

| PES / candidate threshold | Wells | Product sets | Connections | Maximum product total-variation shift |
| --- | ---: | ---: | ---: | ---: |
| CH3S full | 3 | 6 | 10 | — |
| CH3S 0.2% | 2 | 3 | 5 | 0.088 percentage points |
| CH3S 1% or 5% | 2 | 1 | 3 | 0.944 percentage points |
| CH3SS full | 17 | 18 | 91 | — |
| CH3SS 0.2% | 11 | 8 | 28 | 0.535 percentage points |
| CH3SS 1% | 7 | 3 | 10 | 3.380 percentage points |
| CH3SS 5% | 4 | 1 | 4 | 10.126 percentage points |

At 148.142152 kcal/mol the conditional isolated, ergodic CH3SS model gives
95.609% CH3 + S2, 2.087% CH2S + SH and 1.059% H + CH2SS. These are not
experimental photolysis yields. No R9 production scope is automatically
approved by these candidates.

The current full/reduced inputs, native outputs, numerical checks and state
restoration provenance are collected under
`local_runs/mess_review_20261102/deliverables/`. Earlier deliveries remain
historical. Native execution success is deliberately separate from the
pressure-dependent rate audit: some low-temperature off-diagonal coefficients
are grossly nonphysical, not merely tiny negative roundoff. Those tables are
excluded from every retention decision above.
