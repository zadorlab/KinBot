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
an ordinary KinBot PES input:

```json
{
  "title": "acetaldehyde_anl0_f12",
  "smiles": "CC=O",
  "charge": 0,
  "mult": 1,
  "reaction_search": 1,
  "families": ["r12_insertion_R", "h2_elim", "hom_sci"],
  "selected_reactions": {
    "initial": {
      "r12_insertion_R": [[6, 1, 0]],
      "h2_elim": [[3, 0, 1, 6]],
      "hom_sci": [[0, 1]]
    }
  },
  "stop_product_reaction_search": 1,
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
matches it.  Explicit `f12a` method names remain available, while the older
unqualified `caspt2-f12` and `mrci-f12+q` spellings are compatible aliases for
the documented F12b convention.

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

The default distance switch uses CASPT2/MRCI weights `0.2/0.8`, `0.4/0.6`,
`0.5/0.5`, `0.6/0.4`, and `0.8/0.2` at 2.5, 3.125, 3.75, 4.375, and 5.0
angstrom, with endpoint weights held outside that range.  This gives the more
complete CI gap greater weight at short range and the size-consistent CASPT2
gap greater weight at long range.  Linear interpolation avoids a discontinuity.
Setting the weighting model to `arithmetic_mean` reproduces the equal average
used by Ram et al.

The correction is the reconstructed relaxed-path interaction minus the rigid
CASPT2/cc-pVDZ sampling interaction, after each curve is referenced to its
own 30-angstrom asymptote.  The MEP geometries remain those produced by the
configured KinBot VRC relaxed scan; the stored provenance distinguishes this
from the CASPT2-F12-relaxed MEP used in the published perfluoroalkyl work.

`selected_reactions` and `stop_product_reaction_search` are ordinary KinBot
parameters.  The selector uses the atom order generated from
the input SMILES and avoids a hard-coded connectivity ID.  Program paths,
module names, memory per node, and executable locations remain site
configuration rather than molecular input.

The user does not select a CBH rung. KinBot starts from the highest
structurally valid rung through CBH-3 for each connectivity, plans all lower
reference reactions needed to terminate at CBH-0/ATcT, and then selects the
highest rung having one common completed energy label within that reaction.
Species in different ladder reactions may use different ANL fallback tiers.
The checked input is
`examples/anl/acetaldehyde_cbh_anl0_f12/acetaldehyde.json`.

`selected_reactions` limits this deliberately small validation PES to three
published channels. It has no role in CBH balancing or rung selection.

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
exact completed method label for every species in that reaction. It lowers
the electronic tier within that rung before trying the next lower structural
rung. A lower-rung reference reaction independently selects its own best
common tier. Only the terminal CBH-0 boundary reads pinned ATcT formation
enthalpies.

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
  and recipe selection falls to the best complete level shared within that
  one CBH reaction.
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
