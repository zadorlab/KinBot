# Minimum active spaces for the sulfur VRC examples

## Decision and current boundary

The completed local evidence supports **CAS(5,4) for CH3S C–S and CAS(7,5)
for both CH3SS C–S and S–S**. All minimum-space inward and exact-return
checks are complete. The native shared-reference CASPT2 aDZ/aTZ confirmation
now supports these choices, including asymptotes and active-subspace returns.
The matched CASPT2 shift-0.1 check also passes. The user closed the optional
MRCI branch after Molpro 2023 reproduced the 2024 DSYGV failure. Its unresolved
S–S center/tilted calculations are not used in production. The separate
tilted CH3S MRCI root collapse was resolved with the same (5,4) space. See the
[native recovery review](MOLPRO_RECOVERY_REVIEW_2026-10-16.md). These are minimum-space recommendations for corrected
VRC sampling, not claims of full-valence total-energy accuracy or completed
production kinetics.

| Bond channel | Minimum candidate | Orbital content | SA roots (doublet / quartet) |
| --- | --- | --- | --- |
| CH3S C–S | (5,4) | Methyl radical orbital plus all three atomic S 3p orbitals | 3 / 3 |
| CH3SS C–S | (7,5) | Methyl radical orbital plus the complete S2 pi/pi-star set | 1 / 1 |
| CH3SS S–S | (7,5) | CH3S sulfur radical/lone-pair manifold plus atomic S 3p; C–S sigma remains inactive | 6 / 6 |

The two CH3SS spaces have equal dimensions but **different orbitals**. Install
separate models for each bond. Counting electrons/orbitals or aligning orbital
numbers alone cannot establish the selected chemical subspace.

## Production settings selected by the user — October 16, 2026

The final CASPT2-only inputs use:

| Stage/method | Requested setting |
| --- | --- |
| CASPT2-F12 sampling and correction | `IPEA=0.25`, real `SHIFT=0.1` Eh |
| Spin-gap basis estimate | aDZ/aTZ-F12, exponent 3.7 |
| High-spin anchor | ROHF-UCCSD(T)-F12b/aTZ-F12 ground quartet |

Historical MRCI denominator-shift trials remain in the
[diagnosis](MOLPRO_MRCI_DSYGV_RECOVERY_2026-10-16.md). No MRCI energies or
keywords enter the final production route. Same-spin CASPT2 reference
averaging remains explicit. All three quartet correction references evaluate
their complete declared manifold and select its lowest correlated energy;
this handles the observed CH3S root reordering without changing its H0.

This records the user's production choice; it does not alter running jobs,
the confirmation catalogue, or completed native inputs/results.
The original Molpro CASPT2 confirmation used real shift 0.3 Eh and IPEA 0.25.
The subsequent matched shift-0.1 check gives a 0.267019 kcal/mol change in
the exponent-3.7 spin-gap estimate at the sensitive S–S geometry; see the
[native review](CASPT2_SHIFT01_NATIVE_REVIEW_2026-10-16.md). This does not
establish shift insensitivity over the production PES. No additional
long CASPT2 parameter survey is requested. Do not mix 0.3 Eh asymptotic
energies with 0.1 Eh interacting energies when building a production potential.

## Identified Molpro orbital models and remaining state check

The [orbital-evidence explanation](SULFUR_ACTIVE_ORBITAL_EVIDENCE.md) provides
parent and final occupations, atomic character, the reason for the C–S seed
permutation, expansion comparisons and the limits of the continuity checks.

The native aDZ/aTZ records identify the following ordered C1 active spaces.
Labels are one-based `orbital.irrep`, in natural-orbital record `2140.2`:

| Channel | CAS | CLOSED / OCC | Active labels | Asymptotic chemical content |
| --- | --- | --- | --- | --- |
| CH3S C–S | (5,4) | 10 / 14 | 11.1, 12.1, 13.1, 14.1 | S 3p triplet, then methyl radical |
| CH3SS C–S | (7,5) | 17 / 22 | 18.1, 19.1, 20.1, 21.1, 22.1 | S2 pi pair, pi-star pair, then methyl radical |
| CH3SS S–S | (7,5) | 17 / 22 | 18.1, 19.1, 20.1, 21.1, 22.1 | CH3S sulfur pair, then atomic S 3p triplet |

These labels refer to the separate reviewed, basis-specific checkpoints in
[native_orbital_orderings.json](../molpro_confirmation/native_orbital_orderings.json),
which retains file paths, hashes, and native-output references. They are not
portable indices for a fresh RHF calculation. The original CH3SS C–S seed's
17.1/19.1 interchange has already been applied in these ordered checkpoints;
do not repeat it. S–S canonical record `2141.2` puts atomic S at 18–20 and the
CH3S pair at 21–22, preserving the same five-orbital active span.

The retained local inward/return tests and native aDZ/aTZ return tests show
no detected active-subspace loss or core intrusion at the checked points.
This supports these minimum chemical models; it does not require identical
individual orbital labels within degenerate active subspaces. Correlated
state reordering does occur. In particular, the CH3S tilted lowest MRCI
doublet and quartet have dominant reference components 2 and 3, respectively.
An identical integer root across methods therefore does not specify identical
state character. The completed CH3S simultaneous-root checks resolve the
earlier independent-call collapse; see the
[native root review](MOLPRO_RECOVERY_REVIEW_2026-10-16.md).

The remaining S–S acceptance requires review of far/center/tilted native
outputs, all requested roots and available overlap/reference-vector evidence,
energy selectors, and matched interaction/spin gaps. A successful process
exit alone is insufficient. If that audit passes, the remaining handoff work
is installing the reviewed orbital and state/surface models in the production
inputs, rather than another broad active-space search. Sparse validation does
not certify every later Monte Carlo geometry; production checks must continue
to detect loss of the intended subspace or state solution.

## Literature criterion

The selection target is the interaction potential relative to the appropriate
dissociation zero, its angular dependence, relevant spin gaps and repeatable
sampling. In corrected VRC-TST, a cheaper orientation potential is combined
with a higher-level reference-cut correction. This is the construction tested
in [Klippenstein, Georgievskii and Harding, PCCP 2006, 8, 1133–1147](https://doi.org/10.1039/B515914H).

[Knepp et al., PCCP 2007, 9, 4315–4331](https://doi.org/10.1039/B705934E)
describes the chemical reasons for minimal alkyl (2e,2o) references, enlargement
for fragment orbital degeneracy/resonance, and a complete O2 pi manifold for
alkyl + O2. This supports retaining complete relevant fragment manifolds,
not automatically using full valence. The sulfur assignments above are our
chemical inference tested by the retained calculations; they are not active
spaces prescribed for these molecules by that paper.

CASPT2 correlates inactive valence electrons outside the CAS reference, subject
to the frozen-core convention. Thus CASPT2(2,2) is not bare CASSCF(2,2), and
neither its absolute difference from full-valence CAS nor a universal error
bound follows from its success in corrected methyl-recombination calculations.
A radial correction can remove a radial bias; it cannot generally repair an
incorrect angular dependence, missing asymptotic states or orbital switching.
The previous preference for CH3S (13,8), based on an assistant-chosen angular
screening threshold, is withdrawn. Its results remain useful comparison data.

## State averaging and the VRC-TST objective: literature recheck, 16 October 2026

The relevant acceptance target remains the corrected interaction surfaces that
determine the variational flux. A general argument about XMS near conical
intersections is not, by itself, a requirement to replace the current sampling
method or repeat the active-space search. Conversely, VRC-TST does not make
electronic-state coupling irrelevant by definition. Changes of treatment should
be motivated by their effect on the kinetic bottleneck and contributing surfaces.

The NH2O precedent is
[Klippenstein, Mulvihill and Glarborg (2023), section 2.2](https://doi.org/10.1021/acs.jpca.3c05181)
([manuscript](https://backend.orbit.dtu.dk/ws/portalfiles/portal/337661919/nh2_o_revised_final.pdf)).
Its reported direct-sampling choices are:

| Channel/region | Active space | State average | Method/basis |
| --- | --- | --- | --- |
| NH2 + O, inner | (5,4) | two states | IPEA-shifted RS2/cc-pVQZ |
| NH2 + O, outer | (5,4) | three states | IPEA-shifted RS2/cc-pVQZ |
| NH + OH, inner | (3,3) | one target state | IPEA-shifted RS2-F12/cc-pVTZ-F12 |
| NH + OH, outer | (5,4) | two states | IPEA-shifted RS2/cc-pVQZ |

The article does not identify these as MS/XMS calculations or provide the
Molpro reference-density keywords. It supports the strategy, not an assertion
that its inputs exactly match our `REFSTATE`/`WEIGHT` implementation.

The separately inspected [publisher MESS input](https://doi.org/10.1021/acs.jpca.3c05181.s002)
uses unions of distinct electronic-surface contributions: three outer and two
inner NH2 + O tables, and two outer and one inner NH + OH tables. State averaging
of the electronic reference is therefore not an instruction to replace all
physical surfaces by their mean energy. State-average weights, electronic
degeneracies, and reactive surface selection are separate quantities. The local
SI audit remains [NH2O_SI_verification.json](../../../../../../docs/NH2O_SI_verification.json).

Two additional primary VRC-TST examples constrain the interpretation:

- [Glarborg et al., Proc. Combust. Inst. 41 (2025), 105942, section 2.1](https://doi.org/10.1016/j.proci.2025.105942)
  ([full text](https://backend.orbit.dtu.dk/ws/portalfiles/portal/416302663/1-s2.0-S1540748925001567-main.pdf))
  uses two-state-averaged CASPT2(3,2)/cc-pVQZ for the long-range NH3 + OH
  capture flux. This is another explicit use of a small fragment-manifold
  active space with state averaging for a VRC-TST observable.
- [Alarcon et al., Proc. Combust. Inst. 40 (2024), 105372, section 2.3](https://doi.org/10.1016/j.proci.2024.105372)
  ([methods text](https://www.sciencedirect.com/science/article/am/pii/S1540748924001809))
  explicitly reports MS-MR-CASPT2 for propargyl + O VRC-TST, with separate
  ground/excited contributions, active-space comparisons, and geometry/energy
  correction potentials. Thus a blanket claim that VRC-TST never uses or needs
  MS-MR is also incorrect. Its reported small differences between reference
  potentials occur outside the rate-controlling separation range. The section's
  MS-MR wording alone does not identify a particular XMS input option.

The [2006 alkyl-combination study](https://doi.org/10.1039/B515914H) directly
tested corrected inexpensive CASPT2 sampling against higher-level sampling for
the resulting kinetics. That is the appropriate motivation for our separation
between inexpensive angular sampling and the higher-level radial correction.
It is not a guarantee that any angular or state error can be corrected radially.

For these sulfur examples, the actual current CASPT2 calculations use separate
target roots with a common same-spin averaged reference density where multiple
roots are retained. They contain no `MIX`/`XMS` step. MRCI uses simultaneous roots
for the multiroot channels. Reference averaging and final MS effective-Hamiltonian
mixing are distinct in the
[Molpro documentation](https://www.molpro.net/manual/doku.php?id=multireference_rayleigh_schroedinger_perturbation_theory).
The current fixed reference manifolds also are not an exact copy of NH2O's
region-dependent state selection.

**Decision:** retain (5,4), (7,5), (7,5) and the evidence-based VRC-TST validation
scope. An XMS comparison is not a new automatic handoff prerequisite. Any further
comparison must address an observed ambiguity in the contributing interaction
surfaces or a material effect on the corrected flux. This literature review
does not complete the pending S-S MRCI retry or establish a numerical rate-error
bound. Production surface accounting and the native end-to-end kinetic run
remain required deliverables. No calculation settings were changed by this review.

## What was calculated

ORCA 6.1.1 **conventional FIC-CASPT2**, with the exact archived
aug-cc-pVTZ-F12 orbital basis, real shift 0.3 Eh and IPEA zero. This is not
CASPT2-F12. Six orbitals are frozen for CH3S and eleven for CH3SS. State
averages and all requested CASPT2 roots are retained explicitly.

The reference geometries use PBE0-preoptimized fragments and exact C3v (CH3S
axial cut) or Cs (CH3SS reference cuts) nuclear constructions. Tilted probes
deliberately break that reference symmetry. Electronic calculations use C1.
The user accepted approximate rigid geometries for this selection; they are
not relaxed MEPs. KinBot's production relaxed/frozen scans provide the later
geometry-relaxation contribution to the correction potential.

The retained [curated native index](../../../../../../local_runs/qchem_casscf_20261015/curated_evidence/20261016T0625Z/README.md)
reparses completed native outputs and checks input, basis, seed and result
provenance. It lists every native warning. Its final snapshot has **47 complete
points and no pending points**. Native inputs, outputs and orbital records are linked individually.
The [numerical comparison](../../../../../../local_runs/qchem_casscf_20261015/minimum_vrc_evidence.json)
separates raw interaction differences, angular residuals and spin gaps.

## Evidence for the minimum choices

Values below are differences between computational models, **not independent
accuracy bounds**. The angular residual is the candidate-minus-comparator
difference at the tilted geometry after subtracting that difference at the
matched reference-cut geometry. One common lowest-doublet asymptotic zero is
used for each model; genuine excited-fragment splittings remain present.

| Candidate / diagnostic expansion | Lowest-doublet angular residual (kcal/mol) | Lowest-doublet minus lowest-quartet angular gap residual (kcal/mol) |
| --- | ---: | ---: |
| CH3S (5,4) / (7,5), adding sulfur 3s | +0.00039 | +0.01640 |
| CH3S (5,4) / occupied-valence (13,8) | +0.10803 | +0.03367 |
| CH3SS C–S (7,5) / (9,7) | +0.00239 | +0.00095 |
| CH3SS S–S (7,5) / (9,7) | +0.04494 | +0.05703 |

CH3S (5,4) has the complete methyl + S(3P) manifold. Adding S 3s has a very
small effect on the tested doublet anisotropy. The occupied-valence comparator
changes the reference-cut interaction more than the angular spin gap; that is
not evidence that all those occupied orbitals must enter an expensive sampling
reference. The inward path and return are numerically repeatable. Native
[center](../../../../../../local_runs/qchem_casscf_20261015/decision_panel/ch3s_cs_cas5_4_center/panel.out),
[tilt](../../../../../../local_runs/qchem_casscf_20261015/decision_panel/ch3s_cs_cas5_4_tilted/panel.out),
[return](../../../../../../local_runs/qchem_casscf_20261015/ch3s_radial_checks/ch3s_cs_cas5_4/r3p1_reverse/panel.out).

CH3SS C–S (7,5) retains the S2 pi bonding/antibonding partners rather than just
the singly occupied orbitals. The (9,7) expansion does not materially alter the
tested anisotropy. The largest all-root return difference is approximately
0.000011 kcal/mol. Native
[center](../../../../../../local_runs/qchem_casscf_20261015/decision_panel/ch3ss_cs_cas7_5_center/panel.out),
[tilt](../../../../../../local_runs/qchem_casscf_20261015/decision_panel/ch3ss_cs_cas7_5_tilted/panel.out),
[return](../../../../../../local_runs/qchem_casscf_20261015/ch3ss_radial_checks/ch3ss_cs_cas7_5/r3p3_reverse/panel.out).

CH3SS S–S (7,5) retains the two relevant CH3S electronic components and the
atomic S 3p manifold. At 2.8 Å its lowest-doublet interaction differs from
(9,7) by about +0.065 kcal/mol and the D–Q gap by +0.008 kcal/mol. At 2.1 Å
these differences are +0.114 and +0.316 kcal/mol, respectively. Both spaces
describe strong doublet bonding on approach. The minimum-space inner-point
active-overlap singular values are 0.982 or larger, with maximum seed-core
projection below 5e-7. Those overlap diagnostics support continuity of the
active span, not a proof of every electronic state's identity. Native
[center](../../../../../../local_runs/qchem_casscf_20261015/decision_panel/ch3ss_ss_cas7_5_center/panel.out),
[tilt](../../../../../../local_runs/qchem_casscf_20261015/decision_panel/ch3ss_ss_cas7_5_tilted/panel.out),
[2.1 Å](../../../../../../local_runs/qchem_casscf_20261015/minimum_ss_radial/ch3ss_ss_cas7_5/r2p1/panel.out).

All minimum-space return calculations passed, comparing native root labels as
well as sorted spectra. Exact forward/return active-span overlaps were computed
from the native orbital coefficients and AO overlap matrix:

| Channel | Maximum all-root return difference (kcal/mol) | Minimum active-span singular value |
| --- | ---: | ---: |
| CH3S C–S (5,4) | 0.00000986 | 0.9999999999986 |
| CH3SS C–S (7,5) | 0.00001097 | 0.9999999999957 |
| CH3SS S–S (7,5) | 0.00008911 | 0.9999999998799 |

The final [S–S return output](../../../../../../local_runs/qchem_casscf_20261015/minimum_ss_radial/ch3ss_ss_cas7_5/r3p7_reverse/panel.out)
terminates normally and all twelve CASPT2 roots converge. The
[return assessment](../../../../../../local_runs/qchem_casscf_20261015/minimum_ss_radial/path_assessments/20261016T061242202509Z.json)
retains the numerical comparisons and file hashes.

These are sparse diagnostic cuts. Sorting energies within a spin block does
not track a physical state: for example the lowest correlated S–S doublet
changes root index along the inward path. A root reorder is not by itself an
active-orbital swap. Later production sampling must retain its reviewed state
definition, and material crossings or near-degeneracies require state review.

## State-average weights and sulfur degeneracy

Equal weights within each spin block and equal total doublet/quartet weights
are used throughout: 1/6 per CH3S C–S root, 1/2 per CH3SS C–S root, and 1/12
per CH3SS S–S root. This matches
[ORCA's default block-weight convention](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/CASSCF.html#basic-usage).
Weights are kept fixed along the tested path. These orbital-optimization
weights are **not** spin statistical factors, Boltzmann populations, weights
in the averaged spin-gap correction or an instruction to count every root
as an independent reactive surface.

Atomic sulfur contributes three spatial S(3P) components. CH3S C–S therefore
has three roots per coupled spin. CH3SS S–S has two CH3S components times those
three sulfur components: six roots per spin. At the current spin-free Cs
fragment geometry the six are two nearly degenerate triplets separated by
about 1.316 kcal/mol. They are not six degenerate sulfur states.
The independent fragment result gives about 1.308 kcal/mol.
See the native [asymptote](../../../../../../local_runs/qchem_casscf_20261015/decision_panel/ch3ss_ss_cas7_5_far/panel.out)
and [fragment](../../../../../../local_runs/qchem_casscf_20261015/fragment_correlation_checks/ch3s_cas5_4_optimized/panel.out)
and the [CH3S electronic-state literature](https://pmc.ncbi.nlm.nih.gov/articles/PMC7528163/).

## Limits that remain in the native evidence

- Small-space quartet calculations print `Class IJTU ... All eigenvalues
  are <= 0`. All requested CASPT2 roots converge, but the warning's origin has
  not been established. A contraction-cutoff test left the printed energies
  unchanged; that is not proof that the warning is harmless. Retain the
  independent Molpro check of the quartet spin anchor.
- Some excited S–S roots at the inner geometry have negative-denominator
  warnings in **both** (7,5) and (9,7). These are retained, not accepted as
  physical rate surfaces merely because the job converged. The minimum
  doublet and quartet at that point have reference weights about 0.860 and
  0.848. The (9,7) tilt also carries an excited-quartet warning.
- Absolute energies, the maximum over every averaged root, and TZ timings
  cannot by themselves decide the minimum sampling space. Native aDZ-F12
  cost and actual corrected-potential/flux sensitivity remain production
  checks. Earlier timings were sometimes affected by controller BLAS use.
- IPEA/real-shift insensitivity has **not** been established. All original
  active-space comparisons used IPEA 0 / real shift 0.3 consistently. On user
  instruction, added long IPEA controls were stopped. The native Molpro
  confirmation subsequently used IPEA 0.25 / real shift 0.3, with the
  shared-reference treatment resolving the artificial sulfur-component
  splitting. The user's final production choice is IPEA 0.25 / real shift
  0.1, as recorded above. These choices do not imply zero parameter dependence.

## Molpro and ordinary KinBot handoff

The [minimal Molpro package](../molpro_confirmation/README.md)
prepares 30 CASPT2-F12 points in aDZ/aTZ-F12 and nine sparse aTZ-F12 MRCI-F12+Q
points, using the same geometries. It tests native orbital ordering, basis
projection/reoptimization, asymptotes, angular and inward/return behavior.
It does not duplicate the production relaxation scans. Chemistry, orbital
choices and cluster exclusions are catalogue data, not shared-code branches.

The old Blodgett reference campaign is a standalone example study. The general
`kinbot.orbital_review`/`molpro_reconnaissance` machinery can generate and report
equivalent declared checks, including the native evidence used for selection.
It assists a chemical decision; it does not independently infer a unique CAS
from electron counts or natural occupations.

The [production package](../../../../ch3s_anl1_f12/orbital_reviews/RUNNING.md) now installs
separate basis-specific records for each bond and declares the physical
state/surface model. The intended launch is ordinary KinBot on the two supplied inputs:
CBH/ANL1-F12, constrained relaxed/frozen correction scans, aDZ-F12 orientation
sampling, aDZ/aTZ-F12 CASPT2-F12 power-3.7 spin gaps with the
RHF/ROHF-UCCSD(T)-F12b reference, RotdPy and MESS. This route is implemented
and locally tested; these examples have not yet completed it natively.
The user handles SSH and launch, with assistant support for fixes through
completion. See the maintained [twelve deliverables](../../../../../../docs/ORBITAL_RECON_DELIVERABLES.md)
and [deferred cleanup plan](../../../../../../docs/REPOSITORY_CLEANUP_PLAN.md). No cleanup is authorized
until both end-to-end examples finish.
