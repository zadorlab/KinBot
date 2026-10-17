# Orbital reconnaissance and KinBot production deliverables

Approved by the user on October 15, 2026, with the literature and handoff
amendments below. This is the maintained scope and completion checklist.
Update statuses and evidence links as work completes; do not replace outstanding
requirements with a new plan. The example-specific choices belong in example
data, not in general KinBot implementation.

**16 October scope clarification:** the user selected CASPT2-only production
spin gaps (IPEA 0.25, real shift 0.1 Eh), with aDZ/aTZ extrapolation and the
existing high-spin CC anchor. MRCI is a historical benchmark. The latest user instruction selects ordinary
single-bottleneck VRC for both examples. Existing optional inner/outer backend
hooks remain available but are not used here. The particular physical surfaces and statistical
assumptions must be documented separately; successful electronic convergence
does not establish them. The unsubmitted example models and new automatic VRC
spin-orbit offsets/population redistribution were removed. Retain and audit the
existing sulfur fine-structure treatment in thermochemistry separately. A
confirmation geometry called `inner` is not evidence of an inner kinetic
bottleneck. See the [current production scope](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/CASPT2_PRODUCTION_TRANSITION_2026-10-16.md).

## Current status: CASPT2 production handoff prepared

This table is the current status; the historical MRCI experiments below remain
provenance, not production prerequisites. The user has ended MRCI retries.

| Item | Status | Remaining work before or after handoff |
| --- | --- | --- |
| 1. Literature/minimum active spaces | Complete for handoff: (5,4), (7,5), (7,5) | Physical reactive components and distinct SA weights are explicit in the production models. |
| 2. Local ORCA evidence | Complete | Retain 47 curated native points and reproducible analysis. |
| 3. Molpro confirmation | CASPT2 complete, including matched shift 0.1 | MRCI branch closed by user; no further benchmark needed. |
| 4. General JSON controls | Complete for handoff; regression-tested | Single JSON per PES; per-bond indices, basis records, roots, weights, same-spin references and lowest-quartet selector are explicit general controls. |
| 5. General orbital reconnaissance | Complete for handoff; documented and tested | Chemical selection still requires review; it is not a unique-CAS oracle. |
| 6. Production inputs | Complete and locally checked | Six portable orbital records, ordinary VRC surfaces, full inputs and launch/preflight helper supplied. Native production remains pending. |
| 7. End-to-end production follow-through | r5 reached CBH-ANL and VRC, then failed; saved-results resume prepared | Four blockers repaired in code: WFU restored seed geometry after requested geometry, SH VPT2 framework-label rejection, H2S native framework mismatch, HSS reference bond-order label mismatch. Native confirmation and ANL/RotdPy/MESS completion remain pending. |
| 8. Integration/general defaults | Profile routing and minimum recovery audited | Current coverage and exact full-suite results are in the [DFT restart validation](../examples/anl/tools/dft_restart_validation.json). Native production remains pending. |
| 9. Evidence/documentation | Complete for handoff | Parameter-by-parameter provenance, native evidence, limitations and commands linked below. |
| 10. Cleanup inventory | Example/configuration cleanup complete; native purge deferred | Removed obsolete external model JSONs, caspt2_production and acetaldehyde example. Each PES has one root input with its own checkpoints/evidence/docs below orbital_reviews. Original mixed evidence is retained once; do not purge before production finishes. |
| 11. Consolidated production handoff | Downstream repair and bounded monitor prepared | Same DFT preset, 1% RotdPy target, tight L1/L2/L3 criteria and CASPT2 settings. Both PES summaries appear first; full records remain available in details/JSON views. |
| 12. Previous-run retirement | User selected in-place r5 resume with external failure archive | Helper requires both PES inactive, archives obsolete attempts/logs outside the run, preserves accepted QC and rejected reaction decisions, then resumes both controllers. No unrelated work or shared scratch is removed. |

The [production handoff](../examples/anl/tools/README.md) and its linked full
parameter reference are authoritative for the new runs. “Complete for handoff”
means the requested preparation is done, not that the two PES calculations,
their scientific acceptance, or deferred cleanup have already happened.

The [17 October downstream repair record](PRODUCTION_DOWNSTREAM_REPAIR_2026-10-17.md)
contains native errors, the documented input changes, validation limits and
the decision to retain the broad barrier window for this repair. The 30 A
seed-geometry evaluations are not valid 2.5 A production scan points.

The [species spectroscopy review](SULFUR_ELECTRONIC_LEVELS_2026-10-17.md)
records the complete supplied inventory. Experimental CH3S and SH level pairs
are in both main JSONs for MESS/MESSPF. S, SH and CH3S now explicitly request
one additive experimental-centroid SO correction in ANL (CH3S uses the
approved effective approximation). Saved reference models receive the current
input correction during assembly. All three CAS root manifolds are unchanged. CH3S JT coupling is
a plausible contributor to the earlier negative mode, not a proven diagnosis;
the independent negative curvature required reoptimization in any case.

The user also completed the [27-entry Molpro basis-library check](MOLPRO_BASIS_LIBRARY_CHECK_2026-10-16.md)
on Blodgett 2024.1: all H/C/S lookups loaded, including augmented 5Z-F12.
This resolves the basis-availability concern without changing the ANL pair.

The user's launch transcript reports successful preparation and native backend
preflight (`5 passed in 9.58s`), nine archived campaign directories, and new
CH3S/CH3SS controller PIDs 1008942/1008943 at 18:31:10 PDT on 16 October.
The immediate monitor reported one process and zero Slurm allocations per
PES. This confirms launch only; it is not evidence of completed QC work or
later controller health. The run root is
`/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016`.

Subsequent native logs show both Gaussian initial L2 jobs rejected
`scf(xqc,(xqc,conver=10))` before any SCF iteration. This regression is
resolved in the generic merger/writer and covered by actual rendered-input
tests. PES failure propagation and monitoring are also repaired. See
[the full repair record](GAUSSIAN_PRODUCTION_REPAIR_2026-10-16.md). User requested
a clean restart and ran `_r2`. That attempt passed Gaussian syntax but exposed
the two minimum rejections documented in the [DFT restart plan](MINIMUM_RECOVERY_2026-10-16.md).
That handoff replaced `_r2` with `_r3_dft`. The latest clean restart replaces
`_r3_dft` with `_r4_dft` following the [sulfur QC-SCF failure](GAUSSIAN_SCF_RECOVERY_2026-10-16.md).
The user explicitly chose to stop and archive both runs and restart fresh.
No scientific completion claim follows from local tests.

Later `_r4_dft` evidence revealed that Gaussian can return exit 0 with an
explicit SCF nonconvergence warning. The reader now rejects this case.
The [atomic symmetry correction and continuation assessment](GAUSSIAN_SCF_RECOVERY_2026-10-16.md)
documents the actual converged single-atom check and the later accepted
CH3S channel. No new clean restart is required by the 22:17 snapshot;
an optional r4-to-r5 replacement is available if requested.

The new shift export verifies 58/58 inventoried files and 24/24 correlated
energies. The matched CASSCF/reference energies and active orbital subspaces
are unchanged. The lowest D–Q gap changes by 0.255453 kcal/mol at aDZ,
0.264439 at aTZ and 0.267019 in the power-3.7 estimate. See the
[native shift review](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/CASPT2_SHIFT01_NATIVE_REVIEW_2026-10-16.md).
This is acceptance of the requested matched-point check, not proof of
shift-independent rates or all-geometry root continuity.

## Objective and handoff

Determine the **minimum chemically adequate active space** for the requested
corrected VRC-TST calculation on CH3S C–S, CH3SS C–S and CH3SS S–S. Establish
defensible local evidence, confirm the Molpro implementation in both requested
bases, and provide complete general KinBot inputs and reproducible evidence.
After handoff, **the user launches the actual end-to-end KinBot calculations
on Blodgett, and the assistant helps diagnose and fix issues through completion**.
The user handles all SSH and transfers. Preparation and unit/integration tests
must never be described as completed native production calculations.

## 1. Literature-grounded minimum active-space decisions — supported locally and by native CASPT2

- [x] Recommend a space separately for all three bond channels, including its
  actual orbitals, chemical role, state average, and physical surfaces.
- [x] Review and cite the VRC-TST literature; distinguish supported procedures
  from our inference for the sulfur examples.
- [x] Judge dissociation-referenced interaction energies, relevant spin gaps,
  angular dependence, correct dissociation chemistry, and stable sampling.
  Assess the **corrected potential**, including what the intended correction
  actually recovers. Absolute total-energy agreement with a full-valence or
  larger CAS is **not** the selection objective.
- [x] Start with (5,4) for CH3S C–S and (7,5) for each CH3SS channel; test (9,7)
  or other expansions only where evidence warrants them. These are candidates,
  not predetermined answers. Prefer the smallest adequate space; do not
  promote a larger CAS just because its TZ timing is similar.
- [x] Report native sampling cost at the sampling basis after Molpro confirmation.
  Recovery timings: CH3S aDZ 24–34 s; CH3SS S–S aDZ 204–212 s for all
  requested roots on eight ranks. These are point costs, not Monte Carlo totals.
- [x] Separate numerical convergence, model adequacy, reference uncertainty,
  and statistical Monte Carlo uncertainty. An assistant-chosen screening flag
  is not a user-approved universal accuracy requirement.

**Literature basis and interpretation.** Klippenstein, Georgievskii and Harding
validated inexpensive orientation sampling with a one-dimensional higher-level
correction for alkyl-radical combination kinetics. This supports validating the
resulting interaction potential and kinetics rather than demanding identical
absolute electronic energies from sampling and reference calculations.
[PCCP 2006, 8, 1133–1147](https://doi.org/10.1039/B515914H).

The passage supplied by the user from Knepp et al. reports (2e,2o) for alkyl
radical recombinations, expansion for orbital degeneracy and resonance, and
(7e,5o) for alkyl + O2, where the smaller (3e,3o) model produced unphysical
maxima. It also separates inexpensive orientation sampling from basis and
relaxation corrections. These are chemical selection principles, not a rule
assigning the same electron/orbital count to every sulfur channel.
[PCCP 2007, 9, 4315–4331](https://doi.org/10.1039/B705934E).

For an explicitly chosen state/surface convention, write the sampling
interaction as `V_sample(R,Omega) = E_sample(R,Omega) - E_sample(infinity)`.
A radial correction can remove a predominantly radial bias. It cannot in
general remove an orientation-dependent residual, repair an incorrect
asymptote, or establish which electronic state is being followed. Therefore
check a small number of matched off-axis and return points as well as the
reference cut. Preserve genuine fragment excitation/splitting energies; do
not force every root to zero independently. The detailed averaged-spin
correction must use its declared state and energy conventions.

The adequacy of (2e,2o) in a corrected methyl-recombination model establishes
neither equality to full-valence absolute energies nor a numerical absolute
error for that comparison. We will not invent such an error or assume that
any sufficiently smooth small CAS is accurate. Cost and accuracy are assessed
together for the requested observable.

Current evidence: [minimum comparison](../local_runs/qchem_casscf_20261015/minimum_vrc_evidence.json),
[selection policy](../local_runs/qchem_casscf_20261015/selection_policy.json),
[decision record](../local_runs/qchem_casscf_20261015/selection_decision.json).
The previous CH3S (13,8) recommendation was withdrawn; its numerical comparison
is retained as diagnostic evidence.

The [16 October state-averaging literature recheck](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md#state-averaging-and-the-vrc-tst-objective-literature-recheck-16-october-2026)
compares the NH2O protocol with additional VRC-TST applications, including an
explicit MS-MR example. Keep reference averaging, multistate mixing and physical
surface flux accounting separate. XMS is not an automatic additional handoff
requirement; the criterion remains the accuracy of the corrected VRC-TST flux.

## 2. Complete the necessary local ORCA checks — complete

- [x] Finish minimum-space radial and forward/reverse checks, with targeted
  comparison calculations only for unresolved questions.
- [x] Audit native convergence, chemical active-orbital character, core
  intrusion, state ordering, warnings and restart dependence.
- [x] Document the final choices and residual uncertainty. Conventional ORCA
  CASPT2 in an F12 orbital basis is not CASPT2-F12.

The final local choices are **(5,4), (7,5), (7,5)**, with different orbital
subspaces for the two CH3SS channels. All local inward/return checks are complete;
47 curated points have been reparsed from native outputs. Maximum native-label
return differences are 0.00000986, 0.00001097 and 0.00008911 kcal/mol, respectively.
See [the decision and caveats](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md).
No extra full-valence expansion is needed merely to match absolute energies.

## 3. Minimal one-time Molpro confirmation — CASPT2 accepted; MRCI branch closed

- [x] Supply reproducible aDZ-F12/aTZ-F12 CASPT2-F12 checks of the selected spaces.
- [x] Establish basis-specific orbital initialization, states/asymptotes,
  active-space continuity and representative radial/angular/return behavior
  for CASPT2, reusing the original CH3SS C–S results.
- [x] Supply sparse aTZ-F12 MRCI-F12+Q checks and explicit acceptance review.
- [x] Complete the requested matched shift-0.1 CASPT2 check in both bases.
- MRCI root/CABS investigations are retained as historical benchmarks. The user
  removed MRCI from the production route after the second-build DSYGV failure.

The original run completed 35/39 points. The subsequent recovery completed
20 CASPT2 points plus the MAXDAV20 MRCI point, with two CABS=2 failures and
three blocked dependents. Both separate CABS=3 gates then completed all 18
requested roots. All 597 files in the two new export inventories verify.

Shared-reference CASPT2 retains IPEA 0.25/shift 0.3, restores equivalent sulfur
components to within 0.000289 kcal/mol, preserves physical fragment gaps,
and reproduces return active spaces. The largest root-energy return difference
is 0.001908 kcal/mol. Minimum spaces remain (5,4), (7,5), (7,5).

MAXDAV20 resolved iterative stagnation but separate MRCI calls reached duplicate
solutions. Four of the seven simultaneous-root supplement points are now
reviewed: CH3S tilted CABS=1 and far/center/tilted CABS=3. All 24 selected
energies pass, the collapsed tilted roots are resolved, and the CABS control
changes the lowest spin gap by only 0.001054 kcal/mol. S–S far subsequently
completed per the controller; its full native export awaits review. Center
failed in the first simultaneous MRCI iteration with `dsygv failed`, leaving
tilted dependent on center. The 0.1 Eh numerical preconditioning trial also
failed at the same first-iteration solve (job 52908243); tilted remains unrun.
The center-only comparison with the installed 2023 executable then reproduced
the same DSYGV failure. The user closed this MRCI branch; its unresolved
points are excluded from production. All five completed multistate points are preserved;
no CASPT2 or active-space survey is repeated. The
[diagnosis and controlled comparison](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/MOLPRO_MRCI_DSYGV_RECOVERY_2026-10-16.md)
distinguish documented controls from the unproven crash mechanism. See the
[current native review](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/MOLPRO_RECOVERY_REVIEW_2026-10-16.md),
[historical diagnosis](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/MOLPRO_CONFIRMATION_RESOLUTION_2026-10-16.md), and
[submission handoff](../MOLPRO_CONFIRMATION_HANDOFF.md). The three S–S results
are not relabelled as successes by the CASPT2-only production decision.

## 4. General user controls in KinBot JSON — implemented and checked for handoff

- [x] Verify per-reaction and per-stage CAS counts, exact orbital selection or
  reviewed ordered records, basis-specific initialization, state averages,
  selected roots/surfaces, weights and energy conventions.
- [x] Expose necessary numerical convergence, iteration, memory, retry and
  scheduler controls with validation and documented defaults.
- [x] Verify that the chosen controls reach native inputs and survive restarts.
  No molecule, bond, node name or local path may be hardcoded in shared logic.

Explicit MRCI `correlation_options.state_policy="simultaneous"` computes the
full same-spin manifold, records energy-order labels and selects the matching
relaxed-Davidson array/gradient slot. The default independent policy is unchanged.
- [x] Preserve explicit root/reference conventions and native-confirmed shared-
  reference CASPT2 settings. Do not infer equivalent character from integer
  root labels, or impose multiroot work on a one-root case. Historical MRCI
  simultaneous-root controls remain available but are unused in production.

## 5. General orbital reconnaissance procedure — revised and locally tested

- [x] Make planning, chemical candidate review, checks, interpretation and
  production handoff usable with minimal explicit user input for other systems.
- [x] Explain which operations are automated and which physical/chemical
  choices require review; fail clearly when required evidence is absent.
- [x] Distinguish the standalone Blodgett/example studies from the general
  tool. Demonstrate that equivalent declared checks can provide the evidence
  used here, without claiming a unique active-space oracle.

Procedure: [ORBITAL_REVIEW_WORKFLOW.md](ORBITAL_REVIEW_WORKFLOW.md).

## 6. Complete production inputs for both examples — prepared and locally checked

- [x] Configure all three hom_sci channels with the selected per-bond CAS counts.
- [x] Record the user's October 16 production settings: CASPT2-F12
  `IPEA=0.25`, real `SHIFT=0.1` Eh for sampling/correction. The requested
  matched check now passes. MRCI is no longer used in production.
- [x] Install their native-confirmed, basis-specific orbital and surface models.
- [x] Provide the requested CASPT2-F12/aDZ-F12 sampling and aTZ-F12 relaxed scan,
  CASPT2-F12 aDZ/aTZ power-3.7 spin-gap estimate referenced to
  RHF-UCCSD(T)-F12b/aTZ-F12, including geometry-relaxation correction.
  Retain the TZ-only gap and extrapolation difference for review.
- [ ] Complete CBH-ANL at the requested ANL1-F12 levels, RotdPy to the requested
  Monte Carlo uncertainty, and MESS with checked formation enthalpies,
  thermochemistry and fluxes.
- [x] State the actual level/root used for relaxed geometry scans. A composite
  energy correction must not be described as a composite-gradient optimization.

## 7. Honest readiness and production follow-through — continuing

- [x] Report separately: configured, statically checked, tested, native
  confirmed, and completed end to end.
- [ ] After handoff the user launches both PES calculations on Blodgett; assist
  with bug fixes and scientific/output checks throughout those runs.
- [ ] Do not promise convergence, a completion date, or zero cluster failures.
  Retry recoverable I/O/infrastructure errors with provenance while preserving
  completed valid results; do not disguise scientific failures as successes.

## 8. Generality and integration audit relative to master — complete for handoff

- [x] Audit added features for molecule-specific assumptions, configuration
  propagation, failure/restart behavior, and connections between workflow stages.
- [x] Fix identified defects and run appropriate regression/integration checks.
- [x] Record the audit's actual coverage and remaining native verification;
  do not claim every changed line or external program has been certified.

Earlier recovery regression: 1,652 tests and 391 subtests pass. The current
root-protection change passes 400 targeted tests, including seven new native-input
transformations, energy/gradient indexing and simultaneous-call monitoring. This includes
all 39 original and 26 supplemental native-input transformation paths, not
new native calculations. See
[the audit](GENERALITY_AUDIT_2026-10-15.md) for coverage and native limitations.
The final CASPT2-only handoff pass adds 571 passing focused tests and 31 subtests,
plus five separate native RotdPy/MESS/MESSPF checks. These counts are different
test selections and must not be added to the historical counts as unique tests.
The [parameter reference](../examples/anl/ch3s_anl1_f12/orbital_reviews/PARAMETERS.md)
records explicit input provenance against master. New lowest-manifold energy
selection is regression-tested; its constituent native root calculations,
rather than a newly submitted wrapper calculation, supply the electronic evidence.

## 9. Current documentation and retained native evidence — verified

- [x] Curate conclusions with direct references to retained native inputs,
  outputs, orbital records, geometries, settings, versions and hashes.
- [x] Provide reproducible generation, execution, analysis and transfer steps.
- [x] Mark superseded conclusions and retain their provenance. Archive failed
  attempts needed to explain solver/state decisions without mixing them into
  the accepted production evidence.
- [x] Verify the final archive and replay of analysis from retained results.

All 1,469 archived members passed checksum verification. Isolated replay with
archived analysis source reproduced all 47 curated native observations exactly,
with no missing points, and recomputed the return/orbital-overlap analyses.
The root [handoff record](../MOLPRO_CONFIRMATION_HANDOFF.md) identifies the
committed transfer bundle and the separately retained full native archive.

## 10. Deferred cleanup inventory of the full repository — initial inventory complete

- [x] Inventory material to keep, archive outside version control, ignore,
  or eventually remove, with reasons and dependencies.
- [x] Include generated runs/scratch, duplicate exports, historical documents,
  caches, temporary files and misplaced artifacts; untracked source/tests are
  not automatically extraneous.
- [ ] **Do not perform this cleanup until both examples complete end to end
  in KinBot on Blodgett.** Retain native evidence and reproducibility records.

Initial plan: [REPOSITORY_CLEANUP_PLAN.md](REPOSITORY_CLEANUP_PLAN.md), with
[full top-level inventory](REPOSITORY_CLEANUP_INVENTORY_2026-10-15.json).
Refresh after native production; no cleanup has been performed.

## 11. Consolidated handoff — CASPT2 production launch prepared

- [x] Give all three CAS decisions and their evidence/limitations.
- [x] Give the minimal Molpro validation package, commands and native acceptance checks.
- [x] After native review, give the final production launch/configuration instructions.
- [x] Answer explicitly whether the standalone study and general tool supply
  equivalent evidence, and what remains before a normal end-to-end run.
- [x] Include the cleanup inventory and link this maintained checklist.

## 12. Preserve and retire previous Blodgett work — subsequent native runs retained

- [x] Account for the current old campaign: user snapshot reports 59 validated,
  two failed and seven running expanded-space S–S MRCI jobs.
- [x] Prepare scoped inventory, metadata backup, seed-hash checks, scheduler
  checks, cancellation and post-stop audit without deleting native files.
- [x] Identify complete matching geometry/basis/CAS points for reuse review;
  reuse hash-bound reviewed seed records. Changed geometries or equal CAS
  dimensions alone do not justify reusing old numerical energies.
- [x] Diagnose both reported failures: reference-CI nonconvergence and CABS
  singles orthogonalization after converged MRCI. Do not accept partial energies.
- [x] User ran the transition and one-time confirmation; native exports and
  subsequent recovery/gate outputs are retained with hashes.
- [x] Supply a scoped transition that inventories all named previous campaigns,
  verifies the new seeds and environment first, then stops and archives the old
  directories while preserving native files. Its scope/preservation/duplicate
  protections are tested. The user executes the live inventory and retirement;
  local exports do not certify the current scheduler state.

Historical MRCI defaults remain CABS=1 and relaxed Davidson. Both CABS=3 gates
pass native review. The MRCI supplement explicitly chooses simultaneous roots,
energy-ordered `ENERGD(i)` and matched coupled-CABS interactions; no old
calculation is silently assigned those new conventions. MAXDAV20 completed
numerically but did not protect separate roots from collapse. That
historical submission cancelled/deleted nothing. The new production helper
retires only the explicit superseded campaign list after its preflight succeeds.

## Maintenance record

- Native CAS/MRCI evidence: general recon/production renderers now share documented orbital/CI printing. Per-output `.orbital_audit.json` files extract native occupations, AO coefficients, state/reference information, input declarations and hashes; 32 retained native CASPT2 outputs are regression fixtures. These reports record evidence and inconsistencies, not a universal no-swapping verdict. See [the diagnostic contract](MOLPRO_NATIVE_ORBITAL_DIAGNOSTICS.md).

- Final restart audit: self-contained inputs retain all existing settings and six
  checkpoint hashes. Added explicit bounded ANL L3 minimum recovery before
  dependent energies; tests cover recovery and bounded failure. Current full
  regression: 1,856 tests plus 424 subtests, zero failures/skips. Native recovery
  and the end-to-end acceptance remain pending. Obsolete untracked example
  folders are archived, not deleted, by the restart helper.

- 2026-10-16, production L2 follow-up: the corrected Gaussian input reaches a
  small-force CH3S stationary point, but its fresh Hessian has a −587.334 cm⁻¹
  mode. Independent energy/force displacements confirm negative curvature.
  Added general bounded minimum recovery in the normal worker, detailed
  rejection reporting and a two-PES monitor regression. The
  [evidence and retry plan](MINIMUM_RECOVERY_2026-10-16.md) retain the L2
  method. CH3SS subsequently rejected its CH3S L1 UMA product (−53.062 cm⁻¹).
  At the user's request both examples now select the general
  `b3lyp-b2plyp-anl` preset, with B3LYP-D3(BJ)/cc-pVDZ at L1 and unchanged
  B2PLYP-D3(BJ)/cc-pVTZ at L2. Profile propagation repairs cover preliminary
  reaction paths, scans and IRC options. The current scoped restart archives
  both `_r2` trees and launches `_r3_dft`, monitoring both. Native recovery and
  end-to-end production remain unverified. No active-space or VRC model changes.

- 2026-10-15: user approved the eleven deliverables, added the literature-based
  minimum-space criterion, and clarified that handoff is followed by user-run
  Blodgett production with continued assistant bug-fix support.
- Checklist created while local S–S validation is still running. Unchecked
  items are not claimed complete by the existence of this file.
- 2026-10-16 05:55 UTC: user declined additional long local IPEA checks. The
  added controllers were stopped, preserving the original S–S return check.
  Existing CAS comparisons consistently use real shift 0.3/IPEA 0; this is
  not proof of shift insensitivity. Proposed Molpro/production controls are
  IPEA 0.25/real shift 0.3, with native confirmation before approval.
  The production shift choice in this historical entry is superseded by
  the user's later October 16 decision below; completed evidence is unchanged.

- 2026-10-16 06:11 UTC: final local S–S return terminated normally with all
  twelve CASPT2 roots converged. Added deliverable 12 from the user's request
  to preserve, inspect and retire superseded Blodgett work before confirmation.
- 2026-10-16: final native archive passed isolated replay; user requested a
  committed repository-root transfer bundle and pull/submit handoff. Cleanup
  remains deferred; native Molpro and end-to-end production remain pending.
- 2026-10-16: verified all 833 files in the `_r2` export and independently
  reviewed 35 completed points. Six return subspaces are preserved. User-supplied
  earlier controls resolve the CH3S splitting by equal same-spin reference
  averaging at unchanged IPEA/shift; the previous claim that the large CH3S
  splitting persisted at zero IPEA is withdrawn. Prepared a targeted 26-point
  supplement: four unresolved MRCI points, two CABS controls, and twenty
  corrected-reference CASPT2 points. Existing results are not discarded.

- 2026-10-16: verified all 597 files in the recovery and CABS3 exports. Accepted
  the twenty shared-reference CASPT2 checks and both CABS3 asymptote tests;
  identified duplicate roots in the completed MAXDAV20 MRCI output. Prepared
  seven MRCI-only simultaneous-root checks. Minimum CAS choices are unchanged;
  complete production configuration remains pending native reference closure.

- 2026-10-16, after preparation of the S–S denominator-shift retry: user
  selected final production IPEA 0.25, CASPT2 real shift 0.1 Eh, and MRCI
  iteration denominator shifts 0.1 Eh. Record these for final installation;
  do not change running jobs or relabel historical 0.3 Eh evidence. Modest
  CASPT2 relative-energy sensitivity is a working expectation, not a proven
  result for the full PES. No added long parameter survey was requested.
  The [decision record](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md) now lists all
  three checkpoint-bound active orbital sets, the already-applied C–S seed
  permutation, and the distinction between preserved active spans, physical
  state reordering, and the remaining S–S collapse/solver acceptance checks.
- 2026-10-16: on the user's request, added the detailed
  [orbital-selection evidence](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/docs/SULFUR_ACTIVE_ORBITAL_EVIDENCE.md). Re-extracted
  three larger-parent occupation/character tables from the retained archive,
  rechecked 30 completed native output hashes and their natural occupations,
  and recorded all tables and source hashes in `native_occupation_evidence.json`.
  No new electronic-structure jobs or production-approval changes were made.

- 2026-10-16 17:34 PDT: shift-0.1 export independently reviewed (58 files, 24 energies). User closed MRCI after Molpro23 reproduced DSYGV. CASPT2-only production remains the active scope; final integration and handoff are still in progress.

- 2026-10-16, final scope clarification: ordinary single-bottleneck VRC only in these examples. Added a general, explicit lowest-correlated-manifold energy selector for the CH3S quartet correction reference; shared CASPT2 reference weights remain fixed. It does not claim universal state tracking or permit a mismatched gradient. The aDZ/aTZ exponent remains 3.7. Final handoff integration is still being checked.

- 2026-10-16, handoff preparation complete: installed the same general ground-
  quartet minimum policy in all three models, retained explicit target roots,
  and completed inputs, parameter provenance, native backend checks and scoped
  transition tests. The user launches both PES runs; production acceptance and
  repository cleanup remain future work. No claim of universal root continuity
  or successful end-to-end kinetics is made.
