# Generality and production handoff audit — 15 October 2026

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
`examples/anl/orbital_reviews/validation_candidates.json`. Its preparation
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

The sulfur examples still need native correlated results, a user-approved CAS
and state/surface model, method/basis checks and the native production run.
Local regression success cannot replace those calculations. The next finite
39-job batch, 20-job cap, monitoring and export commands are in
[SULFUR_CORRELATED_VALIDATION.md](SULFUR_CORRELATED_VALIDATION.md).

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
