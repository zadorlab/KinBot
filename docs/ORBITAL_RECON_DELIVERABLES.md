# Orbital reconnaissance and KinBot production deliverables

Approved by the user on October 15, 2026, with the literature and handoff
amendments below. This is the maintained scope and completion checklist.
Update statuses and evidence links as work completes; do not replace outstanding
requirements with a new plan. The example-specific choices belong in example
data, not in general KinBot implementation.

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

The [16 October state-averaging literature recheck](SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md#state-averaging-and-the-vrc-tst-objective-literature-recheck-16-october-2026)
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
See [the decision and caveats](SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md).
No extra full-valence expansion is needed merely to match absolute energies.

## 3. Minimal one-time Molpro confirmation — CASPT2 and CH3S MRCI accepted; S–S MRCI pending

- [x] Supply reproducible aDZ-F12/aTZ-F12 CASPT2-F12 checks of the selected spaces.
- [x] Establish basis-specific orbital initialization, states/asymptotes,
  active-space continuity and representative radial/angular/return behavior
  for CASPT2, reusing the original CH3SS C–S results.
- [x] Supply sparse aTZ-F12 MRCI-F12+Q checks and explicit acceptance review.
- [ ] Finish correlated MRCI root protection and matched coupled-CABS interactions;
  then install final reference-stage settings.

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
One center-only comparison with the installed 2023 executable is prepared,
with identical electronic input and original seed. It is a build diagnostic,
not an established cure. All five completed multistate points are preserved;
no CASPT2 or active-space survey is repeated. The
[diagnosis and controlled comparison](MOLPRO_MRCI_DSYGV_RECOVERY_2026-10-16.md)
distinguish documented controls from the unproven crash mechanism. See the
[current native review](MOLPRO_RECOVERY_REVIEW_2026-10-16.md),
[historical diagnosis](MOLPRO_CONFIRMATION_RESOLUTION_2026-10-16.md), and
[submission handoff](../MOLPRO_CONFIRMATION_HANDOFF.md). The three S–S results
are not accepted before execution and review.

## 4. General user controls in KinBot JSON — audited; simultaneous MRCI confirmed for CH3S

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
- [ ] After the S–S check, apply the confirmed root-protection guidance to the
  general multiroot workflow and production stage models. Preserve explicit
  state conventions; do not assume identical CASSCF/CASPT2/MRCI integer root
  labels mean identical character, or impose multiroot work on a one-root case.

## 5. General orbital reconnaissance procedure — revised and locally tested

- [x] Make planning, chemical candidate review, checks, interpretation and
  production handoff usable with minimal explicit user input for other systems.
- [x] Explain which operations are automated and which physical/chemical
  choices require review; fail clearly when required evidence is absent.
- [x] Distinguish the standalone Blodgett/example studies from the general
  tool. Demonstrate that equivalent declared checks can provide the evidence
  used here, without claiming a unique active-space oracle.

Procedure: [ORBITAL_REVIEW_WORKFLOW.md](ORBITAL_REVIEW_WORKFLOW.md).

## 6. Complete production inputs for both examples — pending native choices

- [x] Configure all three hom_sci channels with the selected per-bond CAS counts.
- [x] Record the user's October 16 production settings: CASPT2-F12
  `IPEA=0.25`, real `SHIFT=0.1` Eh for sampling/correction; MRCI Davidson
  iteration denominator shifts 0.1 Eh for pairs, singles and internals.
  Validate the latter in the current S–S retry. These are recorded choices,
  not a claim that final production inputs are already installed or validated.
- [ ] Install their native-confirmed, basis-specific orbital and surface models.
- [ ] Provide the requested CASPT2-F12/aDZ-F12 sampling and aTZ-F12
  CASPT2-F12/MRCI-F12+Q averaged-spin-splitting correction referenced to
  RHF-UCCSD(T)-F12b, including the configured geometry-relaxation correction.
- [ ] Complete CBH-ANL at the requested ANL1-F12 levels, RotdPy to the requested
  Monte Carlo uncertainty, and MESS with checked formation enthalpies,
  thermochemistry and fluxes.
- [ ] State the actual level/root used for relaxed geometry scans. A composite
  energy correction must not be described as a composite-gradient optimization.

## 7. Honest readiness and production follow-through — continuing

- [x] Report separately: configured, statically checked, tested, native
  confirmed, and completed end to end.
- [ ] After handoff the user launches both PES calculations on Blodgett; assist
  with bug fixes and scientific/output checks throughout those runs.
- [ ] Do not promise convergence, a completion date, or zero cluster failures.
  Retry recoverable I/O/infrastructure errors with provenance while preserving
  completed valid results; do not disguise scientific failures as successes.

## 8. Generality and integration audit relative to master — local audit complete

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

## 11. Consolidated handoff — Molpro launch prepared; production handoff follows native review

- [x] Give all three CAS decisions and their evidence/limitations.
- [x] Give the minimal Molpro validation package, commands and native acceptance checks.
- [ ] After native review, give the final production launch/configuration instructions.
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
- [ ] Reconcile the original broad campaign retirement inventory at production
  handoff; the later confirmation exports do not by themselves count every
  earlier allocation or preserved job.

General defaults remain CABS=1 and relaxed Davidson. Both CABS=3 gates now
pass native review. The MRCI supplement explicitly chooses simultaneous roots,
energy-ordered `ENERGD(i)` and matched coupled-CABS interactions; no old
calculation is silently assigned those new conventions. MAXDAV20 completed
numerically but did not protect separate roots from collapse. This new
submission cancels/deletes nothing and preserves independently running work.

## Maintenance record

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
  The [decision record](SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md) now lists all
  three checkpoint-bound active orbital sets, the already-applied C–S seed
  permutation, and the distinction between preserved active spans, physical
  state reordering, and the remaining S–S collapse/solver acceptance checks.
- 2026-10-16: on the user's request, added the detailed
  [orbital-selection evidence](SULFUR_ACTIVE_ORBITAL_EVIDENCE.md). Re-extracted
  three larger-parent occupation/character tables from the retained archive,
  rechecked 30 completed native output hashes and their natural occupations,
  and recorded all tables and source hashes in `native_occupation_evidence.json`.
  No new electronic-structure jobs or production-approval changes were made.
