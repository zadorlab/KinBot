# ORCA active space validation for CH3S and CH3SS

This study selects a practical active space separately for CH3S C–S cleavage,
CH3SS C–S cleavage and CH3SS S–S cleavage. The target is the interfragment
potential used by VRC-TST, including its angular dependence and continuity.
A small, smooth calculation also needs the correct chemical manifold and
adequate dissociation-referenced **corrected** interaction energies. A larger
CAS is a diagnostic comparator, not an exact reference or an automatic reason
to expand the sampling space. The maintained scope is
[ORBITAL_RECON_DELIVERABLES.md](ORBITAL_RECON_DELIVERABLES.md).

**Current decision status (October 16, 05:18 UTC):** the leading minimum
candidates are CH3S C–S (5,4), CH3SS C–S (7,5), and CH3SS S–S (7,5), with
(9,7) retained as the S–S comparison. The minimum S–S radial/return checks are
still running. The earlier CH3S (13,8) recommendation is withdrawn: the
assistant-chosen 0.1 kcal/mol flag against that larger occupied-valence model
does not establish a need to expand (5,4). No native Molpro production model
has yet been approved. Both aDZ-F12 and aTZ-F12 transfer checks remain.

The dated entries below preserve the investigation history, including
superseded screening interpretations. They are not a list of additional
mandatory calculations. Current numerical comparisons and provenance are in
[minimum_vrc_evidence.json](../local_runs/qchem_casscf_20261015/minimum_vrc_evidence.json)
and [selection_decision.json](../local_runs/qchem_casscf_20261015/selection_decision.json).

## Decision-focused continuation, October 15 evening

At 02:20 UTC on October 16, the CH3S (5,4)/(13,8) comparisons at 30 Å,
8 Å and 3.1 Å are complete. After using one common dissociation zero per
model, the maximum difference across the six spin-free states is 0.00000668,
0.00011554 and 0.538522 kcal/mol, respectively. At 3.1 Å the lowest doublet
in (5,4) is 0.401647 kcal/mol more attractive than in (13,8). These are
differences between two approximations, not errors against an exact reference.
The finite-distance difference exceeds the declared 0.1 kcal/mol raw screening
target. Whether it can be removed adequately by a radial correction depends
on the matched off-axis calculation, which is still running at this snapshot.

Both 3.1 Å calculations converged all six requested CASPT2 roots. The minimum
same-geometry seed/final active-space principal overlaps are 0.999606 for
(5,4) and 0.999220 for (13,8); these support preservation of the intended
orbital span at this step but do not establish state tracking along a complete
path. The four additional orbitals of (13,8) have exported occupations
1.99997, 1.99997, 1.99991 and 1.99973. This supports the chemical motivation
for testing (5,4), without proving that dynamic-correlation differences are
negligible. The (5,4) quartet calculations print the previously observed
IJTU zero-eigenvalue-class diagnostics; they are retained as review flags,
not silently discarded or treated as validated physical information. The
lowest-doublet comparison above does not involve those quartet diagnostics.

The two CH3SS channels remain queued at that snapshot. No final or production
CAS selection is claimed for any bond. The immutable numerical assessment is
`decision_panel/assessments/20261016T022025776199Z.assessment.json`.

At 02:24 UTC the CH3S tilted pair also completed, and the CH3SS C–S chains
started. The lowest-doublet interaction energies at 3.1 Å are -9.702311 and
-9.300664 kcal/mol for (5,4) and (13,8), respectively; tilting by 25 degrees
raises them by 1.619074 and 1.511044 kcal/mol. The angular discrepancy after
the axial radial correction is thus 0.108030 kcal/mol for the lowest doublet,
and at most 0.174342 kcal/mol across the six sorted states. The small space is
a plausible sampling candidate but misses the declared 0.1 kcal/mol target
at this orientation. This is not proof that (13,8) is correct or that (5,4)
cannot yield adequate corrected kinetics. Its four orbitals are contained
in the larger calculation's active span to a minimum principal overlap of
0.99999914 at 3.1 Å. The two calculations therefore compare essentially the
same chemical radical/p manifold, with extra occupied orbitals in (13,8).

The queued `ch3s_followup_settings.json` combines two targeted questions:
an otherwise matched (7,5) test that adds only the S 3s pair, and inward
radial/return checks for (5,4) and (13,8). The radial distances are 4.0, 2.5,
1.9 and 3.1 Å, with the last point initialized from the bonded side and
using exactly the earlier 3.1 Å geometry. The 1.9 Å rigid-fragment point
tests chemical robustness; it is not a claim to have computed the relaxed
bonded minimum. The waiting radial-only coordinator was replaced before
launch by this combined follow-up; no running native job was interrupted.
It waits on the same execution lock as the decision panel and cannot add
concurrent resource usage. Every additional question and geometry remains
in its own specification JSON.

`watch_assessments.py assessment_watch_settings.json` observes completed
result records and automatically writes new immutable comparisons as the
decision and follow-up panels progress. It does not read unfinished energies
as results or make a chemical selection. `assessment_watch_status.json`
lists the newest report for each declared comparison. The observer exits
after both native coordinators finish.

The broad preoptimized-fragment panel was held before native energy jobs
started. Its prepared inputs and seeds remain intact. The active replacement
is `decision_panel_specification.json` / `decision_panel_settings.json`, run by
`run_parallel_panel.py`. It contains six candidate chains, with four points
each: 30 Å, 8 Å, the initial bonding-region reference cut, and a 25-degree
off-axis orientation at that same separation. Two four-rank workers request
4500 MB per rank, at most 36000 MB in total. Every geometry, orbital seed,
chemical question, state average and predecessor is explicit in JSON.

| Bond channel | First comparison | Chemical question |
| --- | --- | --- |
| CH3S C–S | (5,4) against (13,8) | Does adding occupied C–H valence orbitals materially alter the radical + complete S 3p interaction potential? |
| CH3SS C–S | (7,5) against (9,7) | Is the S2 sigma/sigma-star pair needed beyond its complete pi/pi-star manifold and the methyl radical? |
| CH3SS S–S | (7,5) against (9,7) | Does the C–S sigma/sigma-star pair matter in addition to the CH3S radical/lone-pair and atomic S 3p manifolds? |

These comparisons nominate candidates; the larger member is not automatically
a validated reference. An observed difference determines which orbital group
needs further testing. The initial panel does not establish a minimum-energy
path or sufficient radial/angular coverage by itself. Neighboring radial
points, reverse propagation and independent starts follow for the finalists,
as do reference/basis/shift and Molpro-transfer checks. Different choices for
the two CH3SS bonds remain supported by the per-bond surface-model interface.

The user clarified that tiny geometric symmetry deviations are acceptable for
this screening: qualitative chemical completeness is the starting criterion,
and real scans must test robustness. The prepared replacement happens to have
exact geometric C3v symmetry on the CH3S axial cut and Cs on the two CH3SS
reference cuts; explicit operation residuals are below 1e-8 Å. Changes from
the preoptimized-fragment construction are 0.000213, 0.000219 and 0.002228 Å,
respectively. All original geometries remain recorded. The off-axis tests are
deliberately unrestricted. ORCA still uses `NoSym` electronically; the nuclear
symmetry audit does not claim irrep-constrained orbitals or state tracking.
Approximate symmetry is not itself evidence of an orbital swap.

There is also no universal requirement that the *relaxed* CH3S path remain
C3v: spin-free Jahn–Teller surfaces and the spin–orbit coupled ground state
are different models. The primary [CH3S study](https://doi.org/10.1021/jp040443f)
finds that strong spin–orbit coupling quenches the Jahn–Teller distortion.
Our Cs distorted-fragment calculations are spin-free diagnostics, not a
claim about the final spin–orbit equilibrium geometry. The eventual surface
weights, asymptote and correction must state their spin–orbit treatment.

The literature basis is explicit. The user's
[Knepp et al. passage](https://doi.org/10.1039/B705934E) motivates retaining
complete radical, degeneracy and relevant pi manifolds, then testing whether
expansion changes the potential. Its alkyl + O2 (7,5) result is an analogy
for CH3 + S2, not a sulfur validation. The
[alkyl recombination study](https://doi.org/10.1039/B515914H) validates cheap
orientational sampling with higher-level one-dimensional corrections against
direct higher-level calculations. The recent
[NN-VRCTST methods, section 2.1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12573747/)
separates rigid-fragment sampling from geometry-relaxation and electronic
corrections; expensive optimization of every sampled orientation is not the
procedure. These precedents motivate the workflow, not a universal numerical
acceptance tolerance. The 0.1 kcal/mol screening target is our declared working
criterion and must ultimately be related to the requested kinetic sensitivity.

`assess_decision_panel.py` now calls KinBot's general `assess-vrc` comparison
code. Its immutable reports contain raw interaction errors, missing tests,
and the angular error remaining after a same-distance radial electronic
correction. One common energy zero per model preserves fragment-state gaps.
The adapter explicitly labels sorted spectra as provisional; no state
identity or final CAS choice is silently inferred. This closes a reporting
gap, not all remaining automation: reference-geometry relaxation, physical
state tracking and the final production handoff still require their evidence.

The completed fragment controls provide two quantitative diagnostics. Three
CASPT2 S2 points at the DFT distance and ±0.02 Å give a local quadratic minimum
near 1.89892 Å, compared with the PBE0 value 1.88389 Å. This is a local fit,
not a completed correlated optimization. The CH3S two-doublet vertical gap at
the distorted PBE0 geometry is 1.30848 kcal/mol, compared with 0.230142 kcal/mol
at the original frozen parent geometry. Geometry and the spin–orbit model
therefore matter for interpretation of the asymptotic manifold; the gap must
not be erased by separately zeroing each root. Existing native files and
`fragment_correlation_settings.json` reproduce these checks.

Reproduce a fresh assessment without rerunning energies:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=. \
  .venv/bin/python local_runs/qchem_casscf_20261015/assess_decision_panel.py \
  local_runs/qchem_casscf_20261015/decision_panel_specification.json
```

The assessment records input/output hashes and leaves the native results
unchanged. `decision_panel/assessments/latest.json` points to the newest report.

## Calculation record

The working directory is `local_runs/qchem_casscf_20261015/` under the repository.
The directory name is historical: the active calculations use ORCA 6.1.1.
Q-Chem was discontinued at the user's request. `validation_settings.json` is
the immutable specification of the original fixed-parent 33-job screen;
`optimized_fragment_settings.json` specifies its replacement using independently
preoptimized fragments. Each specifies every
geometry, electron count, orbital count, state manifold, orbital seed, process
count and numerical control. Example chemistry belongs in this JSON and its
seed metadata; the Python helpers do not choose spaces by molecular name.

Each job directory retains these artifacts:

| Artifact | Purpose |
| --- | --- |
| `panel.inp` | Exact native ORCA input |
| `panel.out` and `panel.launcher.stderr` | Native and launcher evidence |
| `execution.json` | Command, return code, runtime, solver settings and input, basis and seed hashes |
| `panel.gbw` and `panel.json` | Final orbitals and native overlap matrix |
| `seed/seed_copy.gbw`, `.json` and `seed.validation.json` | Exact-geometry initial orbitals and validation |
| `result.json` | Numerical convergence and orbital checks; chemical acceptance remains false |
| `independent_review.json` | Additional read-only reanalysis, when present |
| `attempts/` | Retained failed I/O attempts, if any |

`caspt2_validation/status.json` records controller progress.
`caspt2_validation/analysis.json` contains derived comparisons, not native data.
Timestamps use UTC; the experiment began on October 15 in Pacific time.

## Electronic structure definition

These are **conventional fully internally contracted CASPT2** calculations,
using the aug-cc-pVTZ-F12 *orbital basis*. They do not include an F12 correlation
operator. ORCA's implementation and its shift conventions must be checked
against Molpro before transferring the final production choice.

The local screening profile uses a real level shift of 0.3 hartree and IPEA shift zero,
equal state weights within each spin block and equal total weights for the two
blocks. On October 16 at 05:55 UTC the user declined additional long local
IPEA runs and chose native Molpro confirmation with IPEA 0.25 / real shift
0.3. The added controls were stopped. Shift insensitivity is not established;
consistent controls isolate CAS differences only at the tested settings. A
common orbital basis alone does not establish equivalence to Molpro CASPT2-F12.
See the [ORCA CASPT2 manual](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/caspt2.html).

All calculations are spin-free. CH3S has 25 electrons and uses three doublet
and three quartet roots. CH3SS has 41 electrons: its C–S channel uses one
doublet and one quartet root; its S–S channel uses six roots of each spin.
ORCA reorders multiplicity blocks internally. The parser reads the printed
multiplicity and root labels rather than assuming input order.

### State averages and their weights

The atomic sulfur ground term contributes **three spatial components of
S(3P)**, not six. Spin coupling and the electronic degeneracy of the other
fragment determine the molecular root count:

| Channel | Intended fragment manifold | Doublet/quartet roots | Normalized weight per root |
| --- | --- | --- | --- |
| CH3S C–S | CH3 ground state times the three S(3P) components | 3 / 3 | 1/6 |
| CH3SS C–S | CH3 ground state times the spatially nondegenerate ground triplet of S2 | 1 / 1 | 1/2 |
| CH3SS S–S | Two CH3S electronic components times the three S(3P) components | 6 / 6 | 1/12 |

Every current input gives equal total weight to the two spin blocks and equal
weights to roots within a block. These are **orbital-optimization weights**,
not thermal populations, spin statistical factors, weights in the composite
spin-splitting energy, or permission to count every averaged root as an
independent reactive surface. The doublet:quartet multiplicity ratio 2:4 does
not require that ratio in CASSCF orbital optimization. Keep the declared
averaging weights fixed along each diagnostic path.

Equal weights are the usual starting choice for an included manifold;
Molpro recommends equal state weights, while ORCA defaults to equal block
weights and equal roots within each block. These defaults agree for the equal
root counts used here; they need not agree for unequal root counts.
[Molpro state averaging](https://www.molpro.net/manual/doku.php?id=quickstart#state-averaged_mcscf),
[ORCA state averaging](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/CASSCF.html#basic-usage).

The two CH3S components are motivated by its parent X(2E) manifold and its
spin-free Jahn–Teller splitting. They are not two additional sulfur atomic
terms. This assignment is consistent with published CH3S electronic-structure
work, but the numerical splitting depends on geometry and method.
[CH3S electronic states](https://pmc.ncbi.nlm.nih.gov/articles/PMC7528163/).
At the present fixed Cs fragment geometry, the S–S (7,5) 30-A output has two
sets of three roots in each spin block, separated by about 1.316 kcal/mol;
the spread inside either set is below 0.0003 kcal/mol. Thus the six states
are **not all degenerate at that geometry**. See native
[S–S asymptote](../local_runs/qchem_casscf_20261015/decision_panel/ch3ss_ss_cas7_5_far/panel.out)
and the independent
[CH3S fragment calculation](../local_runs/qchem_casscf_20261015/fragment_correlation_checks/ch3s_cas5_4_optimized/panel.out).
Including both components in the diagnostic average preserves a balanced
description of this low-lying manifold; the production surface/population
model must separately state which components contribute to the requested rate.
Neither equal weights nor root counts alone prove state character at every
finite separation. CASPT2 energies remain state resolved, not averaged energies.

The frozen cores are six orbitals for CH3S and eleven for CH3SS, confirmed in
the native CASPT2 output. Subsequent checks must retain the same core definition.
“Inactive” in a CASSCF space is distinct from “frozen” in CASPT2: inactive valence
orbitals still contribute dynamic correlation.

### Candidate subspaces

Orbital counts alone do not define a CAS. These are the intended chemical
subspaces, verified against the seed populations and followed during optimization.
Indices in the seed metadata are zero based; Molpro orbital labels are not.

| Channel | CAS | Intended orbital content |
| --- | --- | --- |
| CH3S C–S | (5,4) | Methyl radical orbital and all three atomic sulfur 3p orbitals |
| CH3S C–S | (7,5) | Add sulfur 3s |
| CH3S C–S | (13,8) | Add the three occupied C–H valence bonds |
| CH3SS C–S | (7,5) | Methyl radical and S2 pi/pi-star manifold; S2 sigma is inactive and sigma-star external |
| CH3SS C–S | (9,7) | Methyl radical and the complete six-orbital S2 3p manifold |
| CH3SS C–S | (13,9) | Add both sulfur 3s-derived orbitals |
| CH3SS C–S | (19,12) | Add the three occupied C–H valence bonds |
| CH3SS S–S | (7,5) | CH3S sulfur radical/lone-pair manifold and free sulfur 3p; C–S sigma inactive |
| CH3SS S–S | (9,7) | Add C–S sigma and sigma-star |
| CH3SS S–S | (13,9) | Add both sulfur 3s orbitals, retaining three inactive C–H bonds |
| CH3SS S–S | (19,12) | Add the three occupied C–H valence bonds |

The (13,8) and (19,12) spaces are larger occupied-valence comparisons, **not
literal full-valence spaces**: they do not include all C–H antibonding partners.
Their suitability as references is itself checked through orbital character,
incremental space changes and the independent Molpro reference calculations.

## Geometry and orbital construction

The initial panel contains three geometries per space: a 30 Å far point, a
center point, and a tilted orientation at the same bond separation as the
center. The center separations are 3.1 Å for CH3S C–S, 3.3 Å for CH3SS C–S,
and 3.7 Å for CH3SS S–S. Exact coordinates and the source archive member/hash
are in each JSON job. `panel_geometries.json` additionally retains the neighboring
radial points. Internal fragment geometries are fixed for these comparisons.
Their relaxation correction is a separate part of the eventual potential.

### Geometry provenance correction and required optimization

The coordinates trace to `examples/anl/ch3s_anl1_f12/orbital_reviews/orbital_reconnaissance.json`
and `examples/anl/ch3ss_anl1_f12/orbital_reviews/orbital_reconnaissance.json`, through the archived
Molpro campaign. `kinbot.molpro_reconnaissance._anchor_geometry` translates a
declared fragment and optionally rotates it; it does not optimize nuclei.
No independently converged electronic-structure geometry optimization or
frequency record has been established for these input structures. Their
original geometry-generation method is not recorded in those JSON files.
Earlier references to "optimized" CASSCF results mean optimized orbitals and
CI coefficients, not optimized nuclear coordinates.

This is a material limitation of the initial screen. In the CH3S input the
three H–C–H angles are 109.47149 degrees, and carbon lies 0.36980 Å outside the
hydrogen plane. The CH3 fragment retained from CH3SS is also pyramidal, with
H–C–H angles 108.91895–109.81313 degrees. Isolated ground-state methyl is planar
with 120-degree angles ([NIST experimental geometry](https://cccbdb.nist.gov/expgeom2x.asp?casno=2229074&charge=0)).
The S2 fragment cut from CH3SS retains its parent S–S distance, 2.094169 Å.
The CH3S fragment cut from CH3SS retains C–S = 1.817793 Å and inequivalent
S–C–H angles. None of these inherited values establishes a relaxed fragment.

Consequently, the 33-job panel is a fixed-parent-geometry electronic diagnostic,
not the final equilibrium-fragment VRC-TST validation. Its energies remain
useful for testing the backend bridge, convergence and matched CAS differences.
The distorted-fragment far-point splitting must not be described as the
equilibrium asymptotic splitting. A constant asymptotic subtraction removes
the frozen fragment energy offset, but does not repair the change in the
interaction potential caused by the wrong fragment shape.

The required geometry work is:

| Channel | Isolated fragments whose structures must be established |
| --- | --- |
| CH3S C–S | CH3 doublet; sulfur atom triplet (no nuclear optimization) |
| CH3SS C–S | CH3 doublet; S2 triplet |
| CH3SS S–S | CH3S doublet; sulfur atom triplet (no nuclear optimization) |

First obtain and verify fragment minima, recording method, spin, optimized
coordinates, nuclear gradients and vibrational checks. An inexpensive geometry
preoptimization is only a starting point. Check the S2 bond length and CH3S
electronic/structural distortion at the correlated level; explicitly distinguish
a state-averaged reference geometry from a minimum on an individual electronic
surface. Do not impose exact CH3S degeneracy simply to remove a splitting.
Build a new common radial/angular panel from the reviewed fragment structures,
keeping its geometry specification distinct from the archived panel. All CAS
candidates must be compared at these same coordinates. Finally, use constrained
optimizations at selected approach separations to measure internal-fragment
relaxation and test whether an orientation-independent correction is adequate.
Fully unconstrained optimization of every angular sample would eliminate the
independent radial/angular coordinates that VRC-TST needs to sample.

The user confirmed that the correction potential is intended to include
geometry relaxation. Write this part explicitly as
`DeltaV_relax(R) = E_ref(R, q_relaxed(R)) - E_ref(R, q_frozen)` along the
specified reference approach, with identical method, state and energy zero
on both sides. Here `q` denotes fragment internal coordinates; the intermolecular
separation remains constrained. At infinite separation the correction must
approach zero when `q_frozen` uses the same method's equilibrium fragments;
otherwise record the asymptotic offset explicitly and use
`DeltaV_relax(R) - DeltaV_relax(infinity)` in an interaction potential whose
zero is relaxed separated fragments. Keep the absolute offset separately for
thermochemical bookkeeping. Add electronic-method and basis corrections
consistently, without counting relaxation twice. The radial
approximation is checked at representative off-axis orientations rather than
assumed universal. Geometry data and electronic labels must be transferred to
the correction-potential construction as well as to the rigid sampling input.

Three initial local jobs are prepared in
`local_runs/qchem_casscf_20261015/geometry_preoptimization_settings.json`:
CH3 doublet, S2 triplet and CH3S doublet. They use unrestricted PBE0, the exact
archived augmented cc-pVTZ-F12 orbital basis, no symmetry constraints,
TightOpt nuclear tolerances and analytic frequency calculations. These are
inexpensive **starting structures**, not CASPT2-validated equilibrium geometries.
The CH3S calculation may lower its symmetry; retain the resulting structure
and electronic evidence instead of forcing two roots to coincide. The settings
record the parent source coordinates, atom subsets and source hashes.
`run_geometry_preoptimization.py` preserves native inputs, outputs, Hessians,
optimized XYZ coordinates, execution metadata and frequencies; its results
explicitly leave correlated-geometry and production acceptance false. It uses
four MPI ranks with 2000 MB working memory per rank. Initially it waited for
the energy lock; after the last original far point completed, it was moved to
a separate geometry lock to run alongside the four-rank controls within the
local resource budget. Two retries are allowed only for explicit I/O failures.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  .venv/bin/python -u local_runs/qchem_casscf_20261015/run_geometry_preoptimization.py \
  local_runs/qchem_casscf_20261015/geometry_preoptimization_settings.json
```

The command serializes geometry jobs; do not start a second copy while the
launched process is active. It preserves an existing output
rather than silently repeating or overwriting an optimization.

Geometry optimization, electronic convergence and CAS adequacy are separate
acceptance requirements. No production CAS is approved until the geometry
requirements above have been satisfied. ORCA supports CASSCF nuclear gradients
and numerical gradients for methods without analytic derivatives; the selected
root and orbital state average must be explicit
([optimization manual](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html),
[CASSCF manual](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/CASSCF.html)).

The basis is reconstructed from archived Molpro XML, with the base cc-pVTZ-F12
functions checked against Basis Set Exchange. The archive contains the augmented
functions actually used by Molpro. `basis/provenance.json`, `molpro_source.xml`,
`cc-pvtz-f12.json` and `aug-cc-pvtz-f12.orca` preserve this definition. The resulting
spherical dimensions are H 26, C 65 and S 74, giving 217 and 291 functions for
CH3S and CH3SS respectively. This avoids substituting a similarly named basis.

`molpro_to_orca.py` converts the archived Cartesian expansion into ORCA's real
spherical functions, including normalization, ordering and angular phases.
The converter currently supports the s through f shells in these basis files;
unsupported basis representations require validation before use. Import into
GBW uses an RHF container for a single spatial coefficient matrix, preserving
the orbital data. Labeling that container CASSCF incorrectly changes ORCA's
interpretation of its dimensions.

`molpro_seeds/*.validation.json` records orbital permutations and source hashes.
For example, (7,5) CH3S swaps source indices 8 and 9; (7,5) CH3SS C–S swaps
16 and 18; (13,9) CH3SS S–S swaps 11 and 14. These are example-specific orbital
partitions recorded as data, not a general molecular rule.

`transport_seed.py` rigidly translates and rotates orbitals with the explicitly
listed fragments, obtains the target overlap matrix from native ORCA, then
orthonormalizes the core, inactive, active and external blocks in that order.
It rejects a non-rigid fragment fit, a rank-deficient orbital block or a failed
native overlap check. The final seed is `seed_copy.gbw`, not the intermediate
`transport_copy.gbw`. The 33 prepared seeds have native orthonormality errors
below 1e-7. This construction avoids projecting a distant sulfur atom onto its
old atomic position after a 27 Å displacement.

The initial Molpro-to-ORCA CH3S tilted check reproduced all six CASSCF energies
within 7e-9 hartree. Its active-space singular values were within 3e-9 of one;
the evidence is in `03_orca_ch3s_cas13_8/molpro_orca_subspace_comparison.json`
and the associated native outputs. This validates that particular orbital and
CASSCF bridge, not the CASPT2 implementation across programs.

### Excluded reference and unsuccessful trials

An archived CH3SS C–S center CAS(19,12) seed placed two sulfur core-like p
orbitals in the active space while retaining two C–H valence bonds inactive.
It is not used as an accuracy reference. The current panel starts from the
reviewed far seed and transports it to each geometry. A separate repaired seed
is retained for provenance, but is not substituted silently into prior results.

The old `02_orca_ch3ss_ss_cas13_9` trial stopped on energy with a gradient above
the target and had the wrong intended orbital partition. The
`04_orca_ch3ss_ss_refined_caspt2` trial was deliberately stopped during unstable
orbital optimization. The older `caspt2_panel/ch3s_cs_cas5_4_far` calculation
was stopped after poor long-distance orbital projection. Its replacement is in
`caspt2_validation/`. These attempts remain on disk and are excluded from the
accepted comparison set. The older small-space tilted pilot also exceeded the
gradient target despite normal termination.

## Numerical controls and execution

The current JSON uses CASSCF gradient 1e-5, CASSCF energy 1e-12 hartree, CI energy
1e-9 hartree and CI residual 1e-8, with 2000 maximum iterations. The tight CASSCF
energy value prevents premature energy-only stopping while the independent
gradient criterion remains practical. It is not Molpro's former squared
coefficient threshold of 1e-20. The optimizer starts with SuperCI_PT and switches
to NR below a gradient of 1e-4. Numerical tolerances are not method-error bounds.

Most jobs use 14 MPI ranks and 2000 MB MaxCore per rank. The (19,12) jobs use four
ranks and 6500 MB per rank to provide sufficient memory for each rank's tensors.
Only one energy job runs at a time on the 48 GiB machine. MaxCore is an algorithm
allocation parameter, not a strict bound on process RSS. The local ORCA binary
requires the MPI installation at `/opt/homebrew/opt/open-mpi-legacy/bin`; the
other installed MPI version caused launcher failures in preliminary tests.

Launch or resume from the repository root, after checking that another local
controller is not already active:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  .venv/bin/python -u local_runs/qchem_casscf_20261015/run_panel.py \
  local_runs/qchem_casscf_20261015/validation_settings.json
```

The controller takes an exclusive lock, runs one immutable input at a time,
retains failed attempts and permits two retries for explicit I/O errors.
Numerical failures are retained for review rather than classified as I/O.
Existing native outputs are never overwritten on resume. Resume acceptance
requires a completed execution record and the same hashes and orbital checks
as a new result. An interrupted record cannot become accepted merely because
some normal-looking output remains. Changes to a scientific input require a
new job directory. Do not rerun the older `prepare.py` over existing trials.

Monitor and regenerate the derived report without changing native files:

```bash
OMP_NUM_THREADS=1 .venv/bin/python \
  local_runs/qchem_casscf_20261015/summarize_panel.py \
  local_runs/qchem_casscf_20261015/validation_settings.json \
  --json local_runs/qchem_casscf_20261015/caspt2_validation/analysis.json
```

## Comparison and decision rules

Every root must show native CASPT2 convergence; a normal ORCA termination alone
is insufficient. CASSCF must satisfy the requested gradient and root coverage.
Orbital checks require an orthonormal export and maximum projection onto the
seed core span below 0.02. Same-geometry singular values between initial and
final active spaces identify subspace changes for chemical review. Rotations
within the active space do not constitute a change of active space.

Some small-space quartet calculations print `ERROR! All eigenvalues are <= 0`
while constructing the IJTU class, yet finish every CASPT2 root. The report
retains that diagnostic and its spin/root context. Its effect requires a
sensitivity check; numerical convergence alone does not dismiss it. Reference
weights, denominators and individual correlation classes must be inspected
alongside the total energy.

The first S–S (7,5) far-point attempt subsequently failed in the second quartet
CASPT2 root with `Diagonalization failed with error code: 71` in
`Tool-Linalg/qcrdiag2.cpp`. CASSCF had converged at gradient 6.68e-6. The ORCA
driver returned zero despite the printed error termination, so neither a zero
return code nor completion of the first CASPT2 root is sufficient for acceptance.
The result is excluded. `control_settings.json` repeats that point with four
MPI ranks at unchanged numerical settings, then tests contraction norm cutoffs
1e-5 and 1e-8 around the default 1e-6. The same cutoff tests for CH3S (5,4) assess
the nonfatal diagnostic there. These controls wait on the main panel's execution
lock and preserve their inputs in separate directories. They do not change the
physical Hamiltonian or active space.

For each candidate, the first screen subtracts one common far-point zero from
all its spin/root energies. This retains relative spin and electronic-state
splittings; there is no independent offset fitted for every root. At identical
geometries the report compares sorted same-spin spectra to the largest planned
space. Sorting is a screening device, not a substitute for tracking state
character or recognizing physical crossings.

The initial working targets are 0.1 kcal/mol for errors in relevant interaction
energies and angular differences, and 0.2 kcal/mol/Å for radial-slope differences.
They are explicit screening targets, not a universal accuracy guarantee or an
experimentally established rate uncertainty. Small offsets in absolute total
energy do not decide the CAS; changes in the interfragment potential do.
`selection_policy.json` makes these targets, proposed radial extensions and the
additional 0.01 kcal/mol numerical-sensitivity and propagation targets explicit.

A final recommendation requires:

0. Reviewed optimized fragment structures, matched candidate geometries, and
   a quantified constrained-relaxation correction where it matters.
1. Matched radial and angular comparisons against a chemically valid expanded
   space, including the approach to the bond, the capture region and asymptote.
2. Forward and reverse orbital propagation, and independent starts where
   needed, to identify hysteresis, core intrusion or an orbital exchange.
3. State assignment at the asymptote and along the scan; same-spin roots that
   correlate to different fragment states need not be degenerate.
4. Sensitivity to numerical thresholds, contraction cutoffs, level shift and
   basis size at representative points; unresolved native warnings remain
   attached to the affected candidate.
5. Matched Molpro CASPT2/CASPT2-F12 and available higher-level checks to test
   transferability, retaining core, state weights, shifts and orbital character.
6. Measured cost and a VRC-TST sampling check with the selected potential.

The smallest candidate meeting the accuracy and continuity requirements is
preferred. If two candidates cannot be distinguished within the reference
uncertainty, that uncertainty is reported rather than inventing a unique optimum.
Different CAS choices for the three channels are permitted. KinBot's general
per-bond channel and bottleneck parameter overrides already provide the routing;
installation of the selected spaces and the native CBH/ANL–ROTD–MESS production
run follow this validation.

## Reproduction snapshots

`capture_reproducibility.py` archives all input specifications, exact seeds,
orbital-conversion provenance, current helper code, documentation and finished
native evidence. It excludes running output files and transient integral files.
Every member has a SHA-256 hash in `reproducibility_manifest.json`; the archive
also has its own `.sha256` sidecar. The licensed ORCA executable is identified
by hash and is not redistributed. The first snapshot is
`reproducibility/orca_validation_snapshot_20261016T0115Z.tar.gz`; subsequent
snapshots retain later results under new filenames.

The 02:27:37 UTC snapshot contains 780 hash-verified files. Its isolated
extraction reproduced the CH3S interaction/angular comparison using the
archived general KinBot modules, without importing the working checkout.
The archive is `reproducibility/orca_validation_snapshot_20261016T022737Z.tar.gz`;
the `.verified.json` sidecar records that independent check. This is a
point-in-time snapshot; later observer reports and native jobs are retained
in the working directory and belong to the next snapshot.

After extracting a snapshot and providing a licensed ORCA installation, create
a fresh run without overwriting the captured calculations:

```bash
python make_replay.py validation_settings.json replay_settings.json \
  --directory replay --executable /absolute/path/to/orca \
  --mpi-bin /absolute/path/to/compatible/mpi/bin
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python run_panel.py replay_settings.json
```

Use Python with NumPy installed. The snapshot records the Python, NumPy and
platform versions used for the audit. Run the replay only when the machine has
resources available; a separate replay has its own execution lock. Native input
paths are regenerated for the extraction location, while the underlying basis
and seed hashes remain unchanged. The seed validation records retain their
original source paths as provenance.

To reproduce the current parallel decision panel, use `decision_panel_settings.json`
instead of `validation_settings.json`, then `run_parallel_panel.py` instead of
`run_panel.py`; the saved resource budget is checked before workers start.
The comparison adapter also requires ASE and the recorded RDKit release.
After extraction, `PYTHONPATH=source python assess_decision_panel.py
decision_panel_specification.json` recomputes the archived comparisons and
binds their evidence paths to the extracted files. Original reports retain
their original provenance paths.

## Blodgett record

The user handles SSH, transfers and remote execution. The restart on October 15
at 17:57 PDT preserved 17 validated pilots and prepared 51 unfinished pilots.
Its backup is `restart_metadata.20261016T005745823630Z`. The user subsequently
confirmed that Molpro jobs were running and that `scontrol` showed the exclusions.
No further restart is needed solely because historical retry counters are large.
The separate pilot-only manifests are not another required set of 68 calculations.
The remote numerical controls and restart procedure are retained in
[the remaining-checks record](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/SULFUR_REMAINING_CHECKS_2026-10-15.md).


## Geometry work resumed, 2026-10-16 01:35 UTC

The original parent-geometry controller was paused after the completed CH3SS
S–S CAS(19,12) far point. Ten of its eleven far-point calculations passed the
numerical/root-coverage checks; S–S CAS(7,5) did not. The first center calculation
had just started during the pause; its partial evidence is preserved and it is
marked `interrupted_for_reprioritization`. This is not a chemical failure.
`caspt2_validation/independent_review_summary.json` additionally records native
class warnings and active-span checks for every far point. A numerical pass
does not clear those warnings or approve a chemical model.

The fragment optimizations now have a separate four-rank queue, allowing them
to run alongside the four-rank numerical controls. Combined requested working
memory is 34,000 MB on the 48-GiB workstation. The subsequent interaction panel
waits for all fragment optimizations to finish and pass their frequency checks.
Its coordinator then acquires the shared energy lock and runs at most two
independent far→center→tilted chains, each with four MPI ranks and 4500 MB per
rank (eight ranks and 36,000 MB total requested working memory). Controls and
fragment checks cannot overlap this coordinator. The bounds, dependency ordering
and result aggregation were exercised with subprocess fixtures before launch.

CH3 completed its PBE0 optimization and harmonic frequency calculation with
normal native termination. The C–H distances are 1.079218–1.079235 Å; H–C–H
angles are 119.9805–120.0098 degrees. Carbon is only 0.000067 Å from the hydrogen
plane. All six vibrational frequencies are positive (lowest 531.77 cm⁻¹), and
<S²> is 0.754686. The exact XYZ, Hessian, frequency list, native output and input
hashes are retained under `fragment_geometry_preoptimization/ch3`. This verifies
the provisional geometry at the stated DFT level, not CASPT2 geometry accuracy.

`geometry_panel_specification.json` declares the fragment atom maps, reacting
bonds, source geometries and far→center→tilted orbital propagation. The generic
`build_fragment_panel.py` aligns optimized fragments around the reacting atoms,
preserves their internal distances and the requested separation, and builds
matched CAS comparisons. A changed internal geometry is initialized by
atom-centered coefficient transport and block orthonormalization; this is an
orbital guess and is never reported as a cross-geometry wavefunction overlap.

`correlated_fragment_specification.json` separately requests S2 CAS(8,6)
calculations at the DFT bond length and ±0.02 Å, plus a two-root CH3S CAS(5,4)
calculation at its optimized geometry. These test correlated bond-length
sensitivity and whether the parent-geometry asymptotic splitting persists after
relaxation. They do not replace a constrained interaction-region relaxation
curve. All cases, spin manifolds and displacements are explicit JSON data.

The remaining nuclear work is the paired frozen/relaxed correction on a
constrained approach path, with representative off-axis checks. Its correlated
energies must use the same state definition and asymptotic zero. No such
correction curve has yet been calculated; the fragment minima are the necessary
starting point.

S2 also completed normally, with an S–S distance of 1.8838925 Å, one positive
stretching frequency of 752.16 cm⁻¹, and <S²>=2.012492. The original 2.0941690 Å
S2 fragment was substantially elongated. CH3S subsequently completed both its
nuclear optimization and frequency calculation normally. Its C–S distance is
1.781110 Å; its three C–H distances are 1.089755, 1.089540 and 1.095494 Å.
All nine vibrations are positive (573.72–3130.56 cm⁻¹), and <S²>=0.755721.
All three isolated-fragment structures are now verified PBE0 minima. Their
`geometry_review.json` files retain pair distances and native-output hashes;
`relaxation_diagnostic.json` records the isolated strain-energy reduction,
which must not be confused with the interfragment relaxation correction.

The 33 optimized-fragment inputs are now prepared. Exact geometry equality was
checked across every CAS candidate at each channel/point, and all eleven initial
far-point seeds passed native import/orthonormality validation. The four
correlated-fragment checks are also prepared. Both controllers await the shared
energy lock held by the numerical-control queue.

The four-rank S–S CAS(7,5) control completed all twelve CASPT2 roots with full
native convergence. Six quartet contraction-class diagnostics remain. The
fourteen-rank failure is therefore not reproduced in this run, but the rank
comparison alone does not prove the cause. Contraction-cutoff sensitivity runs
remain in progress. Changing NThresh from 1e-6 to 1e-5 left all twelve reported
CASPT2 energies identical at their printed precision; the native output confirms
that 1e-5 was actually used. All six quartet warnings remained. This establishes
cutoff insensitivity for this pair only, and does not explain the warnings.
`summarize_controls.py` records root-by-root absolute changes and changes relative
to a common energy zero, native settings and evidence hashes in
`caspt2_controls/numerical_sensitivity.json`. Native diagnostics are retained even
when numerical convergence and root coverage pass.


## Final local decision and native handoff — 16 October UTC

All required local minimum-space calculations are finished. Recommend CH3S
C–S (5,4), CH3SS C–S (7,5), and CH3SS S–S (7,5), with different orbitals for
the two CH3SS channels. The final S–S return terminated normally at
2026-10-16T06:11:23Z, all twelve roots converged, maximum native-label return
error 0.0000891141 kcal/mol and minimum forward/return active-span singular
value 0.9999999998799036. No more local energy runs are needed for this decision.

The [final curated index](../local_runs/qchem_casscf_20261015/curated_evidence/20261016T0625Z/README.md)
contains 47 completed points and zero missing points, with retained warning
text and hashes. The [decision](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md)
replaces earlier tentative recommendations. The retained
[radial plot](../local_runs/qchem_casscf_20261015/curated_evidence/20261016T0625Z/radial_checks.pdf)
shows discrete checks, not an interpolated production PES.

Final archive path:
`local_runs/qchem_casscf_20261015/reproducibility/orca_minimum_space_final_20261016.tar.gz`.
The adjacent `.sha256` verifies the container; `.verified.json` records isolated
reparsing, archived-helper analysis replay, plots and fresh-input generation.
No new quantum-chemistry calculation is implied by analysis replay. The archive
captures documentation at its creation time; subsequent verification status is
recorded in that external report and the maintained checklist.

The Blodgett handoff preserves the old 59 validated results (subject to its
current re-audit), partial native outputs and seed records before retiring the
seven superseded expanded-space MRCI allocations. The exact counts can change
before the user runs it. The new 39-job confirmation uses native CASPT2-F12
in both requested bases, sparse aTZ MRCI-F12+Q and existing default CABS singles.
The relaxed-Davidson `CLUSTER=0` / `ENERGD(1)` selector remains unchanged.
The old CABS orthogonalization failure is a separate native issue, not evidence
that the converged CI equations or printed pre-CABS Davidson value constitute
a completed calculation. Diagnose any repeated failure individually.
