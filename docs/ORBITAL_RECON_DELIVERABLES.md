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

## 1. Literature-grounded minimum active-space decisions — local decision complete; native transfer pending

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
- [ ] Report native sampling cost at the sampling basis after Molpro confirmation.
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

## 3. Minimal one-time Molpro confirmation — 35/39 validated; targeted resolution running

- [x] Supply reproducible commands/configuration for **aDZ-F12 and aTZ-F12**
  CASPT2-F12 confirmation of the selected spaces on the same fixed geometries.
- [ ] Establish basis-specific orbital initialization and ordering, correct
  states/asymptotes, active-orbital continuity, and representative radial,
  angular and return behavior. Reuse valid evidence where possible.
- [x] Include only the sparse requested **aTZ-F12 MRCI-F12+Q** checks needed
  for the reference/spin-splitting comparison, plus explicit acceptance review.
- [ ] Record every method, shift, frozen-core, state and basis convention;
  resolve material differences between local screening and production.

This is a transfer check, not another broad campaign, geometry optimization,
or duplicate of the production relaxed scans. Do not reuse a TZ orbital file
as an unverified DZ orbital model. A count such as (7,5) alone is insufficient
to identify the intended orbital subspace. Catalogue:
[Molpro seeds and geometries](../examples/anl/orbital_reviews/molpro_confirmation/geometry_and_seed_catalog.json).
The [prepared transfer recipe](../examples/anl/orbital_reviews/molpro_confirmation/README.md)
contains 30 CASPT2-F12 points across both bases and nine sparse aTZ-F12
MRCI-F12+Q points. Preparation, input propagation and resume checks pass;
the `_r2` run has now completed 35 points. The first submission exposed
cross-basis record collisions, inherited RHF state inconsistencies, and an
MRCI-F12 excited-reference restriction. The [manual/input audit](MOLPRO_CONFIRMATION_INPUT_AUDIT_2026-10-16.md)
records repairs, native evidence, and the preserved full restart. Local tests
alone did not establish handoff readiness. The
[current resolution record](MOLPRO_CONFIRMATION_RESOLUTION_2026-10-16.md)
verifies all 35 outputs, gives basis-specific natural/canonical orbital
orderings, confirms the existing relaxed-Davidson selector, and distinguishes
the two failures from their two unsubmitted dependents. Existing control data
demonstrate that shared-reference CASPT2 removes the large CH3S splitting at
IPEA 0.25/shift 0.3. The supplemental plan tests that policy on the current
CH3S/S–S geometries in both bases and completes targeted MRCI recovery.
These behaviours require demonstrated closure, not an indefinitely recorded warning.
At 09:45:56 PDT on October 16 the user reports 10/26 supplemental points
validated, five running, ten waiting and one failed CABS=2 CH3S control.
The subsequent native report confirms **both** CABS=2 far points failed:
CH3S quartet root 1 and S–S doublet root 1, after converged MRCI with matching
NSTATI/NSTATR. CABS=2 is rejected for production in these conditions. Exact
internal cause of `strange..` remains unavailable; no cluster fault is shown.
The MAXDAV20 retry had progressed to root 2. Existing completed points and
other running chains remain in place. A separate two-point CABS=3 test is
prepared using the same spaces/references/thresholds, all 18 roots and no
dependent scans. Native success and relative-energy checks remain required;
neither the final CABS settings nor production readiness is claimed.
At 10:12:32 PDT, the new monitor reports 19/26 validated, two failed, four
waiting and one running (S–S aTZ CASPT2 return). The plan/count accounting
indicates MAXDAV20 completed; native energies await export review. Three
waiting points are blocked CABS dependents; the fourth is the last CASPT2
tilted successor. Allow that chain to finish.

## 4. General user controls in KinBot JSON — locally audited and tested; native pending

- [x] Verify per-reaction and per-stage CAS counts, exact orbital selection or
  reviewed ordered records, basis-specific initialization, state averages,
  selected roots/surfaces, weights and energy conventions.
- [x] Expose necessary numerical convergence, iteration, memory, retry and
  scheduler controls with validation and documented defaults.
- [x] Verify that the chosen controls reach native inputs and survive restarts.
  No molecule, bond, node name or local path may be hardcoded in shared logic.

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

Latest recovery regression: 1,652 tests and 391 subtests pass. This includes
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

## 12. Preserve and retire previous Blodgett work — adapter tested; user execution pending

- [x] Account for the current old campaign: user snapshot reports 59 validated,
  two failed and seven running expanded-space S–S MRCI jobs.
- [x] Prepare scoped inventory, metadata backup, seed-hash checks, scheduler
  checks, cancellation and post-stop audit without deleting native files.
- [x] Identify complete matching geometry/basis/CAS points for reuse review;
  reuse hash-bound reviewed seed records. Changed geometries or equal CAS
  dimensions alone do not justify reusing old numerical energies.
- [x] Diagnose both reported failures: reference-CI nonconvergence and CABS
  singles orthogonalization after converged MRCI. Do not accept partial energies.
- [ ] User runs the transition and one-time confirmation on Blodgett. Check
  the actual final preservation counts and queued controls from that run.

The default MRCI convention remains `CABS_SINGLES=1`, `CLUSTER=0` (relaxed
Davidson), selecting the final `ENERGD(1)` after CABS processing. The repeated
minimum-space CABS failure persisted in the CABS=2 recovery/control trial;
the two-point CABS=3 test now isolates the remaining documented coupled
treatment. The stalled external CI has a separate MAXDAV=20 trial. Neither
changes the Davidson selector. All original run/results are preserved.
Actual native completion and energy review are required before approval.

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
