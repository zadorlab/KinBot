# Minimum active spaces for the sulfur VRC examples

## Decision and current boundary

The completed local evidence supports **CAS(5,4) for CH3S C–S and CAS(7,5)
for both CH3SS C–S and S–S**. All minimum-space inward and exact-return
checks are complete. The native shared-reference CASPT2 aDZ/aTZ confirmation
now supports these choices, including asymptotes and active-subspace returns.
MRCI simultaneous-root and matched coupled-CABS reference checks remain before
installation of the final production models. See the
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

The retained [curated native index](../local_runs/qchem_casscf_20261015/curated_evidence/20261016T0625Z/README.md)
reparses completed native outputs and checks input, basis, seed and result
provenance. It lists every native warning. Its final snapshot has **47 complete
points and no pending points**. Native inputs, outputs and orbital records are linked individually.
The [numerical comparison](../local_runs/qchem_casscf_20261015/minimum_vrc_evidence.json)
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
[center](../local_runs/qchem_casscf_20261015/decision_panel/ch3s_cs_cas5_4_center/panel.out),
[tilt](../local_runs/qchem_casscf_20261015/decision_panel/ch3s_cs_cas5_4_tilted/panel.out),
[return](../local_runs/qchem_casscf_20261015/ch3s_radial_checks/ch3s_cs_cas5_4/r3p1_reverse/panel.out).

CH3SS C–S (7,5) retains the S2 pi bonding/antibonding partners rather than just
the singly occupied orbitals. The (9,7) expansion does not materially alter the
tested anisotropy. The largest all-root return difference is approximately
0.000011 kcal/mol. Native
[center](../local_runs/qchem_casscf_20261015/decision_panel/ch3ss_cs_cas7_5_center/panel.out),
[tilt](../local_runs/qchem_casscf_20261015/decision_panel/ch3ss_cs_cas7_5_tilted/panel.out),
[return](../local_runs/qchem_casscf_20261015/ch3ss_radial_checks/ch3ss_cs_cas7_5/r3p3_reverse/panel.out).

CH3SS S–S (7,5) retains the two relevant CH3S electronic components and the
atomic S 3p manifold. At 2.8 Å its lowest-doublet interaction differs from
(9,7) by about +0.065 kcal/mol and the D–Q gap by +0.008 kcal/mol. At 2.1 Å
these differences are +0.114 and +0.316 kcal/mol, respectively. Both spaces
describe strong doublet bonding on approach. The minimum-space inner-point
active-overlap singular values are 0.982 or larger, with maximum seed-core
projection below 5e-7. Those overlap diagnostics support continuity of the
active span, not a proof of every electronic state's identity. Native
[center](../local_runs/qchem_casscf_20261015/decision_panel/ch3ss_ss_cas7_5_center/panel.out),
[tilt](../local_runs/qchem_casscf_20261015/decision_panel/ch3ss_ss_cas7_5_tilted/panel.out),
[2.1 Å](../local_runs/qchem_casscf_20261015/minimum_ss_radial/ch3ss_ss_cas7_5/r2p1/panel.out).

All minimum-space return calculations passed, comparing native root labels as
well as sorted spectra. Exact forward/return active-span overlaps were computed
from the native orbital coefficients and AO overlap matrix:

| Channel | Maximum all-root return difference (kcal/mol) | Minimum active-span singular value |
| --- | ---: | ---: |
| CH3S C–S (5,4) | 0.00000986 | 0.9999999999986 |
| CH3SS C–S (7,5) | 0.00001097 | 0.9999999999957 |
| CH3SS S–S (7,5) | 0.00008911 | 0.9999999998799 |

The final [S–S return output](../local_runs/qchem_casscf_20261015/minimum_ss_radial/ch3ss_ss_cas7_5/r3p7_reverse/panel.out)
terminates normally and all twelve CASPT2 roots converge. The
[return assessment](../local_runs/qchem_casscf_20261015/minimum_ss_radial/path_assessments/20261016T061242202509Z.json)
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
See the native [asymptote](../local_runs/qchem_casscf_20261015/decision_panel/ch3ss_ss_cas7_5_far/panel.out)
and [fragment](../local_runs/qchem_casscf_20261015/fragment_correlation_checks/ch3s_cas5_4_optimized/panel.out)
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
  instruction, added long IPEA controls were stopped. Proposed Molpro and
  production settings are IPEA 0.25 / real shift 0.3; standard settings do
  not imply zero parameter dependence. The chemical minimum is expected to
  be robust, but must be reviewed using the native confirmation results.

## Molpro and ordinary KinBot handoff

The [minimal Molpro package](../examples/anl/orbital_reviews/molpro_confirmation/README.md)
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

After native review, install separate basis-specific records for each bond
and sampling/correction stage, and finish the physical state/surface model.
Then the intended launch is ordinary KinBot on the two supplied inputs:
CBH/ANL1-F12, constrained relaxed/frozen correction scans, aDZ-F12 orientation
sampling, aTZ-F12 CASPT2-F12/MRCI-F12+Q averaged spin gaps with the
RHF/ROHF-UCCSD(T)-F12b reference, RotdPy and MESS. This route is implemented
and locally tested; these examples have not yet completed it natively.
The user handles SSH and launch, with assistant support for fixes through
completion. See the maintained [twelve deliverables](ORBITAL_RECON_DELIVERABLES.md)
and [deferred cleanup plan](REPOSITORY_CLEANUP_PLAN.md). No cleanup is authorized
until both end-to-end examples finish.
