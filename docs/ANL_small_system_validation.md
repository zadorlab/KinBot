# Small-system CBH-ANL0-F12 validation target

## Selected system

Use acetaldehyde (`CC=O`) as the next complete validation system.  It is the
smallest target found that simultaneously exercises a nontrivial CBH rung,
ordinary saddle-point kinetics, barrierless VRC-TST, partition functions,
thermochemistry, and a pressure-dependent MESS network.

The requested network is deliberately restricted to:

1. `CH3CHO <=> CH3 + HCO`, the single barrierless C-C scission;
2. `CH3CHO <=> CH4 + CO`, the conventional tight molecular channel; and
3. `CH3CHO <=> CH2CO + H2`, the adjacent-H tight elimination channel.

No reaction search is launched from the product fragments.  Roaming is a
known physical contribution to the molecular products, so agreement of the
conventional tight channel with the total experimental molecular yield is not
an acceptance condition.  The barrierless dissociation/recombination and the
complete pressure-dependent network are compared with the published
VRC/RRKM work.

The connectivity-derived CBH1 reaction is

```text
CH3CHO + CH4 -> CH2O + C2H6
```

or, in the signed convention used by `kinbot.anl.cbh`,

```json
{"C": -1, "C=O": 1, "CC": 1, "CC=O": -1}
```

CBH0 is used only at terminal ATcT leaves.  Acetaldehyde's 0 K and 298.15 K
formation enthalpies are solved through the same CBH1 tree; no elemental
reference reaction is introduced for the finite-temperature result.

## Intended user input

The production interface should add only the following scientific choices to
an ordinary KinBot calculation input:

```json
{
  "title": "acetaldehyde_anl0_f12",
  "smiles": "CC=O",
  "charge": 0,
  "mult": 1,
  "barrier_threshold": 100.0,
  "conformer_search": 1,
  "high_level": 1,
  "rotor_scan": 1,
  "theory_preset": "uma-b2plyp-anl",
  "composite_method": "ANL0-F12",
  "vrc_tst_correction_model": "spin_splitting",
  "vrc_tst_scan_points": [2.5, 3.125, 3.75, 4.375, 5.0, 7.0, 10.0],
  "rotdpy_run": 1,
  "rotdpy_validation_profile": "production",
  "rotdpy_mess_mode": "production",
  "run_messpf": 1,
  "me": 1,
  "run_me": 1,
  "collider": "He",
  "epsilon": 300.0,
  "epsilon_unit": "K",
  "sigma": 4.5,
  "queuing": "slurm",
  "queue_name": "day-long-cpu",
  "ppn": 8,
  "queue_job_limit": 8
}
```

For a simultaneous automatic VRC and composite calculation, the positive
`queue_job_limit` is a total cap. KinBot freezes the accepted network, assigns
disjoint slot reservations to the two branches, and runs them concurrently.
The composite branch never reads the VRC branch's mutable reaction objects;
it consumes the frozen JSON plan. MESSPF and MESS remain behind a persisted
join record and cannot start from a partially completed branch.
For an automatic run, KinBot also scopes admission to the calculation tree,
serializes scheduler submissions across child controllers, and infers an
omitted scheduler username from the login environment. The portable example
therefore does not require a site-specific `username` field.

The spin-splitting preset keeps the original direct VRC correction as the
global compatibility default; selecting it is one extra scientific choice in
this validation input.  For this test its three reference branches all use
TZ/QZ two-point CBS values:

- RHF-UCCSD(T)-F12b/aug-cc-pV{T,Q}Z-F12 for the triplet interaction;
- CASPT2-F12b(2e,2o)/aug-cc-pV{T,Q}Z-F12 singlet and triplet energies, with
  real shift 0.20 and IPEA shift 0.25; and
- MRCI-F12b+Q(2e,2o)/aug-cc-pV{T,Q}Z-F12 singlet and triplet energies.

The explicit `b` labels remove a Molpro output ambiguity.  An unqualified
F12 command computes both variants: `ENERGY(1)`/`ENERGD(1)` are F12A and
`ENERGY(2)`/`ENERGD(2)` are F12B.  KinBot therefore selects `ENERGY(2)` for
CASPT2-F12b and `ENERGD(2)` for MRCI-F12b+Q.  The latter is generated with
`OPTION,CLUSTER=0`, so `+Q` is the relaxed-reference Davidson correction;
the F12 CABS-singles correction remains included.  The generated input binds
that value to `KBVRCSEL`, and the reader verifies that the public energy key
matches it.  The unqualified `caspt2-f12` and `mrci-f12+q` spellings default
to F12b.  Explicit F12a method names remain available when requested.

Every lower/upper pair uses the same ANL extrapolation:

```text
E_CBS = E_n + alpha_n (E_n - E_(n-1))
alpha_n = (n-1)^3.7 / (n^3.7 - (n-1)^3.7)
```

Thus `alpha_4 = 0.5265464668174269` for the requested TZ/QZ pairs.  The code
also accepts adjacent DZ/TZ, QZ/5Z, 5Z/6Z, or higher pairs and records the
power, coefficient, raw energies, and extrapolated energy in the correction
record.

At each separation the reference singlet energy is reconstructed as

```text
E_S(R) = E_T,CC,CBS(R)
       + w_CAS(R) [E_S,CAS,CBS(R) - E_T,CAS,CBS(R)]
       + w_MRCI(R)[E_S,MRCI,CBS(R) - E_T,MRCI,CBS(R)].
```

The default is now `arithmetic_mean`, the equal-average model used by Ram
et al.; method sensitivity must still be checked. The previously implicit
distance preference was specific to a validation example and has been removed.
`distance_switch` remains available only with explicit
`distance_knots_angstrom` and per-level `weights`. Archived, fully specified
correction records remain readable. A restart from an input that omitted its
weights now has a different correction definition and must regenerate that
correction rather than reuse an incompatible cache.

The correction is the reconstructed relaxed-path interaction minus the rigid
CASPT2/cc-pVDZ sampling interaction, after each curve is referenced to its
own 30-angstrom asymptote.  The production default relaxed path is optimized
with Molpro RS2-CASPT2(2e,2o)/aug-cc-pVTZ.  Conventional CASPT2 uses Molpro's
analytic gradient.  CASPT2-F12, MRCI+Q, and MRCI-F12+Q scan choices use
numerical gradients of their explicitly selected native energy.  In
particular, conventional MRCI+Q differentiates `ENERGD`, while MRCI-F12+Q
differentiates `ENERGD(1)` or `ENERGD(2)` for an explicitly requested F12A
or the default F12B surface.  A completed point is rejected if its constrained
distance changes
or an internal fragment bond breaks.

Gaussian relaxed scans remain available when a Gaussian method is explicitly
selected, preserving the original KinBot input path.  The separate default
UB3LYP/6-31+G(d) fragment calculations only construct HOMO-oriented ROTD_py
pivot points; they do not supply the relaxed scan geometries or energies.

Reaction search, supported families, and individual atom motifs are not
specified in this input. KinBot uses its ordinary defaults to enumerate all
supported unimolecular families and motifs, then accepts channels through its
normal optimization, IRC, connectivity, and barrier checks. Program paths,
module names, memory per node, and executable locations remain site
configuration rather than molecular input.

The user does not select a CBH rung. KinBot starts from the highest
structurally valid rung through CBH-3 for each connectivity, plans all lower
reference reactions needed to terminate at CBH-0/ATcT, and then selects the
highest rung having one common completed energy label within that reaction.
Species in different ladder reactions may use different ANL fallback tiers.
The checked input is
`examples/anl/acetaldehyde_cbh_anl0_f12/acetaldehyde.json`.

The automatic CBH and VRC planners consume the accepted reaction network;
neither requires a user-selected family, atom motif, or CBH rung.

## Automatic feasibility and fallback policy

`composite_method` is the requested electronic ceiling. For `ANL0-F12`, each
generated molecular state is assigned the highest tier allowed by its
connectivity metrics:

1. full ANL0-F12 for at most 16 atoms, 4 heavy atoms, and 4 heavy bonds;
2. L3 geometry and harmonic frequencies for at most 30 atoms, 8 heavy atoms,
   and 8 heavy bonds;
3. an L3 single point on the accepted L2 geometry for at most 50 atoms,
   12 heavy atoms, and 12 heavy bonds; or
4. the accepted L2 energy, geometry, and modes otherwise.

Every assigned state retains the ordered fallback chain
`ANL0-F12 -> ANL0 -> L3 -> L3-SP@L2 -> L2`, beginning at its allowed tier.
All structurally valid CBH reactions from CBH-3 down through CBH-0 and all
their lower-rung reference states are planned before licensed calculations
start. Assembly considers the highest structural rung first and requires one
exact completed method label for every species in that reaction.

At one electronic tier, the strict recipe is considered first. If an optional
additive family is unavailable, automatic assembly constructs explicitly
labeled reduced expressions from the common completed families before it
lowers the electronic tier. The mandatory reference CBS energy and ZPE are
never removed. The reducible ANL0 families are the complete higher-order
difference, DBOC, complete scalar-relativistic difference, and complete
core-valence CBS difference. Both sides of a difference are retained or
removed together. Profiles with fewer omissions sort first; equal-size
omissions drop the higher-order increment before DBOC, scalar relativistic,
and core-valence corrections. Every omission and the terminal task statuses
are written into the composite JSON.

Every completed species also exposes compatible subsets of its components.
Consequently, a species that completed more work can share the best common
reduced expression with one whose optional task failed without rerunning or
discarding the completed calculations. If no common expression is available,
selection lowers the electronic tier within that rung and then tries the next
lower structural rung. A lower-rung reference reaction independently selects
its own best common tier. Only the terminal CBH-0 boundary reads pinned ATcT
formation enthalpies.

If L3 geometry fails, only its dependent L3 harmonic calculation is removed;
independent single points use the accepted L2 geometry. Rejected or failed
VPT2 falls back to the same-level harmonic ZPE. Exhausting a declared retry
policy marks that component family unavailable and prevents its unstarted
partner jobs, after which assembly continues down the same ordered ladder.

## Implementation provenance

`kinbot/anl/cbh.py` is a direct RDKit/networkx graph implementation of the
generalized CBH definitions in Ramabhadran and Raghavachari, JCTC 2011,
DOI `10.1021/ct200279q`, with the radical-state extension tested against the
same atom, bond, charge, and radical invariants. KinBot does not import or
vendor AutoCBH, TCIT, or another public CBH reaction-generator library.
Repository history attributes this implementation and its tests to the
KinBot composite branch commits; public packages are not dependencies of the
module.

## Resource boundary

- Ordinary KinBot jobs written with the bundled Slurm template explicitly
  request the selected partition's discovered maximum wall time. Custom queue
  templates retain full control of this directive.
- Every individual Slurm job has a maximum wall time of 24 hours.
- The global live-job limit is eight.
- L1 uses the configured FairChem model.
- L2 geometry, frequencies, VPT2, and hindered rotors use the selected
  B2PLYP-D3(BJ)/cc-pVTZ profile.
- L3 and post-geometry ANL components start together once their preferred
  geometry is available. If L3 optimization becomes unavailable, only the L3
  harmonic calculation is skipped; every independent single point stages on
  the accepted L2 geometry.
- A component that exhausts its declared restart/resource policy is marked
  unavailable.  Unstarted tasks in the same correction family are skipped,
  and recipe selection first forms the best common reduced expression at the
  same tier, then follows the electronic and structural rung fallbacks above.
- Different CBH ladder nodes may use different best available composite
  levels.  Every species within one evaluated CBH reaction must use the same
  selected level.

## Acceptance checks

The final immutable report must prove:

1. the three requested channels, and no other channels, reached the network;
2. each tight transition state has one imaginary mode, a successful IRC, and
   the requested products;
3. one complete ROTD_py calculation contains the correction potential, pivot
   definitions, converged Monte Carlo surfaces, and native hashes;
4. every accepted ANL component is parsed from a hash-matching native output;
5. the CBH1 equation is balanced and recursively terminates only at pinned
   CBH0/ATcT records;
6. MESSPF covers all network and CBH-reference species, including accepted
   L2 hindered rotors;
7. Hf(0 K), Hf(298.15 K), and NASA fits reconstruct from the same CBH tree;
8. MESS completes and reports both pressure-dependent dissociation and the
   reverse association rates;
9. comparison tolerances are recorded rather than used as input data; and
10. restart, reparse, queue-reconciliation, retry exhaustion, component-family
    skipping, and ladder fallback are covered by deterministic tests.

The live chemistry run should not intentionally waste licensed calculations
to exercise every failure path.  Deterministic dispatcher tests cover those
paths, while the live run proves that their records and state transitions are
compatible with real Gaussian, Molpro, CFOUR, MRCC, ROTD_py, MESSPF, and MESS
outputs.

## Literature anchors

- S. J. Klippenstein, L. B. Harding, and B. Ruscic, *J. Phys. Chem. A* 2017,
  121, 6580-6602, DOI `10.1021/acs.jpca.7b05945`: ANL0, ANL0-F12, ANL1,
  geometries, component energies, and formation enthalpies for core
  combustion species.
- K. S. Gupte et al., *Proc. Combust. Inst.* 2007, 31, 167-174, DOI
  `10.1016/j.proci.2006.08.048`: direct CASPT2 VRC/RRKM acetaldehyde C-C
  dissociation and a 0.85 dynamical correction.
- L. B. Harding, Y. Georgievskii, and S. J. Klippenstein, *J. Phys. Chem. A*
  2010, 114, 765-777, DOI `10.1021/jp906919w`: pressure-dependent branching
  on the acetaldehyde decomposition surface.
- H. Ram et al., *J. Phys. Chem. A* 2025, DOI `10.1021/acs.jpca.4c07388`:
  RHF-UCCSD(T)-F12b plus CASPT2-F12 and MRCI-F12+Q spin splitting,
  multireference correction potentials, and VRC-TST statistical factors.
  KinBot makes the underlying F12b array selection explicit as documented
  above.
- The Molpro manual sections on explicitly correlated methods, multireference
  perturbation theory, and MRCI define the RHF-reference UCCSD(T)-F12,
  RS2-F12 `SHIFT`/`IPEA`, and Davidson-corrected `ENERGD` conventions used by
  the generated inputs.
- ATcT acetaldehyde record `75-07-0*0`: an independent thermochemical
  comparison, not a CBH1 input reference.
