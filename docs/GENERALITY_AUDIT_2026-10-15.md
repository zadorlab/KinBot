# Generality and production handoff audit — 15 October 2026

## Latest local verification (16 October UTC / 15 October PDT)

The final full local regression run passed **1,622 tests and 391 subtests**
in 73.45 seconds. The retained log is
`local_runs/qchem_casscf_20261015/reproducibility/final_handoff_regression_verified_20261016.log`.
This includes non-sulfur H2/N2 preparation, explicit target Cartesian geometry
validation, separate sampling/correction basis records, independent fixed
orbital-propagation chains, projected starting guesses and spin-block weights.
Native electronic
convergence is not simulated by those tests. The two warnings were the existing
deprecated `imp` import and denied optional Sella cache access. These checks do
not establish an executed production kinetics run.

Subsequent confirmation-export tests verify archive member hashes, original
seed retention, re-auditing of incomplete native results, refusal to overwrite,
and refusal to collect while scheduler state is live or unavailable. After
the final example-control edits, nine example/confirmation tests passed
(`local_runs/qchem_casscf_20261015/reproducibility/handoff_examples_20261016.log`).
The examples now explicitly request IPEA 0.25 and real shift 0.3 for both
sampling/scan CASPT2 and the CASPT2 spin-gap component. Historical recipes keep
their original shifts independently of changes to production input. Local
ORCA comparisons at IPEA zero do not establish parameter insensitivity.
The catalogue-control test also changes the SCF tolerance/iteration limit,
correlation reference policy and sparse point selection, then verifies the
rendered native input. Shared convergence values require decimal powers of ten.

Resolved recon plans now retain the requested ROTD Monte Carlo controls and
validation profile. Reports distinguish CAS-selection evidence from later
production relaxation, corrected-potential and flux checks. Backend
confirmation uses matched fixed geometries; it does not create a duplicate
relaxation/correction workflow. Geometry, fragment/state choices, candidate
orbital groups, thresholds and site settings remain input data. The two sulfur
examples are evidence and configuration under the example/local-run directories,
not runtime dispatch rules. No molecule-name or cluster-node special cases were
found in the recon/production runtime search; this limited static check is not
a line-by-line certification of every change relative to master.

The accuracy-preset audit also found an implicit chemical policy:
`anl1-f12` selected `full_valence` when the input omitted a CAS choice. That
default has been removed. The preset now requires an explicit orbital policy
or orbital/surface model; an intentional full-valence request remains valid.
The existing sulfur examples already supply a candidate policy. Their sampling
basis is now explicitly DZ-F12, as requested; relaxed-path and spin-splitting
levels remain TZ-F12. Historical comparison recipes explicitly retain their
TZ-F12 level. The focused example/review/ROTD/multistate checks passed **84
tests**, including distinct basis-specific WFU fixtures. Native DZ records
still need confirmation. Explicit source-basis provenance and
`project_orbitals: true` now permit a hash-bound checkpoint to supply a
projected starting guess in a different basis. This does not approve the new
basis model: its reoptimized orbitals and asymptote must be reviewed before
production installation. Unacknowledged cross-basis initialization is rejected.

The fixed-chain transfer package renders 39 jobs across nine chains without
duplicating the initial points as separate pilots. Tests verify one-time source
swaps, exact return geometries, requested bases, numerical thresholds and
immutable configuration checks. Its registry-aware monitor reuses the general
monitor. Those tests use synthetic checkpoint bytes and a mocked scheduler;
they are preparation/restart tests, not native Molpro confirmation. Chemistry,
geometries and cluster settings are confined to the example catalogue.

This follow-up examines the `composite` additions relative to local `master`
(`e8bec02b398235c0444e8ce922f2b711079f6768`), also checking the locally available
`origin/master` (`5c6525a7db09534d2e70cd689760eb50a9b5c618`). These are recorded
comparison bases, not a claim that the remote was freshly fetched. The branch
started at `b7eb5aa`. The diff includes orbital reconnaissance, CBH/ANL,
state-resolved VRC/ROTD/MESS, scheduler/runtime handling, and the earlier
stereochemical/thermochemical features. Static searches covered the added
runtime code for molecule names, identities, fixed site paths and special-case
acceptance. Execution/handoff paths were inspected directly and the complete
local regression suite was run. This is not a proof that every supported
molecule converges or that every native integration is correct.

## Findings corrected

| Finding | Change |
| --- | --- |
| ANL dispatcher permitted a historical UHF reproduction only for a hardcoded methyl benchmark tag/doublet | Explicit source-reproduction provenance now defines that diagnostic exception for any compatible open-shell state. Normal automatic ANL remains RHF/ROHF-UCC. |
| Empty rank-exact graphs required a particular workflow name | Acceptance now checks the actual molecular electron count and declared exact zero. |
| The hand-authored end-to-end verifier always compared rates to methyl recombination | It now audits hash-bound native rate tables for arbitrary networks. The published methyl fits remain an explicitly invoked literature diagnostic. Passing the general audit does not claim physical accuracy. |
| That verifier selected the only bimolecular model in the whole network and assumed a zero-energy flux boundary | It now selects the connected barrier endpoint and accepts the declared, hash-bound zero-flux threshold, including a negative native energy. |
| Default spin-gap weights had a fixed distance-dependent preference from an example | No implicit knots/weights remain. The default equal mean still requires validation; `distance_switch` requires explicit data. Old fully specified records remain readable. Inputs omitting weights now resolve differently and require a new correction calculation. |
| Automatic finalization expected one scalar ROTD file per reaction | It now follows selected scalar or multistate channels, validates every physical-state leaf and inner/outer stage, and checks the actual MESS flux-table hashes. Unselected optional exits follow the MESS selection policy. |
| Reviewed orbital selection discarded the optimizer/tolerances validated in recon | The state model can carry explicit `casscf_options` through selection, correction, gradient and sampling renderers. Defaults remain unchanged when omitted. |
| Correlated recon could reuse canonical orbitals from an older WFU | Each optimization now writes current canonical and natural records separately and selects the appropriate record explicitly in CASPT2/MRCI checks. Explicit correlated nonconvergence is rejected even after a normal exit. |
| Recon review did not expose the eventual ROTD request | Plans now retain resolved methods, bases, scan/surface grids, temperatures and available stage overrides. Reports supply model-specific potential data and secants without automatically approving accuracy/continuity. |

## What remains data or an explicit model choice

There are no CH3S/CH3SS names, campaign paths or proposed CAS counts in the
shared reconnaissance, VRC/ROTD or automatic CBH/ANL production logic. The new
example candidate spaces, chemical rationales and seed anchors live in
`examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/validation_candidates.json`. Its preparation
script derives methods and radial points from the linked KinBot requests.

Published energies, fits and geometries in `anl.literature`/`anl.kinetics` are
explicit diagnostic datasets, not calculated substitutes or automatic
production acceptance conditions. Element valence definitions, backend syntax,
physical constants, method recipes and configurable numerical defaults remain
part of the implementation. A preset is not an accuracy certificate. Explicit
`full_valence` and `radical_pair` policies remain user-requested model choices;
fragment multiplicities alone cannot justify a chemically truncated CAS.

The general orbital review supports user-defined AVAS targets, reviewed
natural-orbital partitions, independent channel/inner/outer spaces, arbitrary
supported spin/root averages, basis-specific records and optional new radial
points. Chemical assessment must consider reacting/radical orbitals, resonance,
fragment terms, correlation and consistent dissociation—not a universal
occupation cutoff or CAS size. Both accurate reference-relative potentials and
smooth state-consistent sampling are required before minimizing cost.

## Automatic route after selection

The production route is the ordinary KinBot input, not the historical
hand-authored `anl.end_to_end` verification manifest. `orbital_review
install-selection` installs the user's selected record and counts globally or
into an explicit bond/stage while preserving CBH/ANL, ROTD and MESS controls.
After the surface model's physical states/offsets/weights are reviewed, normal
`kinbot INPUT.json` handles reaction discovery, automatic CBH references and
ANL jobs, VRC correction scans, ROTD physical-state sampling, MESSPF/MESS and
`anl_automatic/final_report.json`. Different inner/outer spaces are supported
without inventing a vdW complex for `hom_sci`.

The sulfur examples now have completed local ORCA CASPT2 evidence, including
the minimum S–S radial/return check. Final orbital/surface approval, native
Molpro method/basis checks and the production run remain. Local regression
success cannot replace those calculations. The current transfer package is
[the minimal Molpro confirmation](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/molpro_confirmation/README.md).
The older [correlated validation notes](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/SULFUR_CORRELATED_VALIDATION.md)
retain the historical broad campaign and its recovery evidence.

## Native follow-up

The first correlated batch exposed a shared MR energy-selector bug that the
mocked tests had missed. The [failure and recovery audit](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/SULFUR_CORRELATED_VALIDATION.md#native-selector-failure-and-recovery)
corrects the method contract, documents the native evidence and provides an
in-place resume preserving 18 completed calculations. Production readers now
reject undefined-variable assignments. The supplied complete native output
is a regression fixture. The scheduler retry guard also handles requeues
using Slurm's restart count rather than rejecting every repeated job ID.

The subsequent native diagnostics distinguish CASSCF macroiteration exhaustion,
correlated reference-CI iteration exhaustion and missing/interrupted execution.
Generic state-model controls now preserve separately validated CASSCF and MR
iteration budgets through production selection. `prepare-incomplete` re-audits
finite independent validation pilots, preserves successful points, and prepares
only unfinished points with their original verified seeds. No molecule names
or job indices enter that repair. Accounting export records unavailable `sacct`
without blocking native evidence collection; an unavailable accounting service
does not establish why an allocation disappeared. The updated full regression
suite passed **1,451 tests and 391 subtests**, with the same two existing warnings.

## Verification

The complete local suite passed **1,434 tests and 391 subtests**. The subsequent
source-partition count guard has its own focused regression run. The only
warnings were the existing deprecated `imp` import and denied Sella cache
access outside the workspace.
Preparation was also exercised against the real exported source configurations
and geometries: 39 jobs, 260 correlated roots, one MULTI block per job, no
repeated AVAS or incompatible DIAB, and explicit current canonical records for
CASPT2. That preparation check mocked WFU bytes and scheduler snapshots; it did
not run Molpro or validate native correlated accuracy. The launch/export shell
blocks were syntax checked. Existing uploaded archives remain untouched.


## Final handoff audit

The final targeted set passed **135 tests** (12.20 seconds), covering confirmation,
examples, general VRC correction/state rendering, orbital review and solver controls.
The new transition tests verify seed hashes before cancellation, scheduler failure
before any signal, scoped shutdown, preserved native files/metadata, and a stop for
possible exact-geometry/space reuse. The final handoff retains the working relaxed-Davidson selector and default
CABS-singles CI (1), per user instruction. The old (9,7) CABS orthogonalization
failure is documented and remains a native risk to diagnose if repeated.
No speculative alternate-method control is added to the handoff.
Confirmation memory is explicitly propagated into the execution environment and
restored afterwards; hostile inherited SBATCH memory values are covered by tests.

The final full-suite rerun initially found six stale ROTD fixture expectations
for the former (13,8)/(9,7) example spaces. Expected native occupations were
updated to the selected (5,4)/(7,5) spaces (CH3S closed=10/occ=14;
CH3SS closed=17/occ=22). No production behavior was changed to satisfy those
expectations. The complete rerun then passed 1,622 tests and 391 subtests.
The Davidson reader test now includes all three Davidson/Pople variants and
a pre-CABS fatal-error case; the existing runtime selector remains unchanged.

## October 22 update: remove reference examples from the installed package

This supersedes the earlier allowance for diagnostic datasets inside
`anl.literature` and `anl.kinetics` above. All four historical ANL benchmark
geometries/energies, methyl recombination fits and the paper PES catalogue now
live under `validation/data/`. Source-reproduction and paper-comparison helpers
live under `validation/`, outside the installed package. The remaining generic
ANL comparator requires an explicit reference file; native rate auditing has no
implicit reference reaction. A fixed four-core benchmark cap and a fixed
molecule demonstration in `bond_combinations` were removed. Packaging includes
only `kinbot*`, excluding `kinbot.tests*`; validation and example directories are
not runtime packages. Regression fixtures retain explicit example evidence.

The source scan also flags these pre-existing chemistry rules for transparency:
`exceptions.py` contains legacy O2/CH2/OH identity exceptions; `stationary_pt.py`
contains elemental/fragment multiplicity defaults; `mess.py`, `product_complex.py`
and `tpl/mess_fragment_OH.tpl` implement the OH electronic-state model;
`pes.py` has the bare-proton zero-energy convention. These are existing physical
model rules, not benchmark results inferred from this campaign, and were not
removed as though they were example validation. Reaction-family patterns,
physical constants, method definitions and configurable solver limits likewise
remain general model implementation. This audit does not establish convergence
or physical accuracy for every possible molecule.

The new geometry/frequency refinement follows declared dependencies and native
convergence controls. VPT2 recovery derives point groups from native diagnostics;
no species, DOI, run path or observed example energies enter either decision.
