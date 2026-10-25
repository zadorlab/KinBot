# Chemical orbital review and production handoff

New spin-splitting selections now require a hash-bound native CC/MR
reference review in addition to the orbital review. See the
[review gate and receipt commands](VRC_DIRECT_CBS_HANDOFF_R7.md#future-spin-splitting-review).
Direct MR selections do not require a CC reference. A scientific controller
migration can stage an additive orbital proposal and install it together with
a correction-model change, preserving unchanged scan graphs.


The general entry point is `python -m kinbot.orbital_review`. It uses the
existing Molpro/Slurm execution machinery. There are no molecule-name branches
or prescribed sulfur active spaces in the implementation.

Production review requests also have installed, molecule-independent tools:

```bash
python -m kinbot.requested_orbital_recon prepare REQUEST.json CASE_INPUT.json RECIPE.json SCREEN
python -m kinbot.requested_orbital_recon candidates SCREEN CANDIDATES.json COMPARISON
python -m kinbot.vrc_relaxation_recon prepare PRODUCTION_ROOT PANEL_PLAN.json RELAXATIONS
python -m kinbot.requested_orbital_recon review SCREEN CASE_INPUT.json SELECTION.json
```

Recipes supply the molecular identity, geometry, states, explicit candidate
spaces, distances and native seed evidence. Preparation does not infer or approve
a production model. The `review --apply` and `review-batch --apply` commands
publish an explicitly reviewed selection through live controls. Rejection uses
the shared reaction-decision controls described in `LIVE_DRIVER_CONTROLS.md`.
The legacy `examples/anl/tools/requested_orbital_recon.py` and
`vrc_relaxation_recon.py` paths are compatibility entry points to these modules.

New active-space preparation and planning require an explicit `probe_definition`
with counted electrons/orbitals, chemical rationale and electronic-state model.
An intentionally broad valence diagnostic instead requires
`"probe_space": {"mode": "full_valence"}`. Merely setting
`probe_active_space=true` no longer authorizes a valence CI. Legacy manifests
remain readable and already prepared jobs can be monitored/resumed; reproducing
or expanding an implicit legacy recipe requires a new explicit choice. Do not
edit a running campaign's hash-bound input to add this field.

Requested-review recipes and batch plans may carry `status: withdrawn` or
`status: proposed`; preparation rejects either before creating a campaign.
An absent status retains the existing explicit-recipe interface. `status: ready`
only permits preparation, and is not scientific approval or production acceptance.

Independent prepared campaigns can run together with
`python -m kinbot.review_campaign REVIEW/inventory.json --queue-all`:

```json
{
  "schema": 1,
  "campaigns": [
    {"name": "orbital_comparison", "kind": "orbitals", "directory": "comparison"},
    {"name": "geometry_panel", "kind": "relaxations", "directory": "relaxations"}
  ]
}
```

Directories are relative to the inventory file. This queues each finite
inventory while respecting native orbital-predecessor dependencies; it does not
accept results or promise immediate allocations. Without `--queue-all`, the
workers use the configured reserved review capacity. Explicit `--allow-node`
options affect the review workers only. Repeating after owners exit adopts
their completed and live calculations. A concurrent duplicate cannot take the
existing workers' locks. Run-specific example launchers supply plans and legacy
owner identities; the screen preparation, relaxation, scientific validation,
execution, supervision and numerical recovery implementations reside in KinBot.

Shared admission retains uncertain `sbatch` outcomes. A purged known allocation
can release its reservation only with a healthy scheduler query, an explicit
invalid-ID response, a matching submission receipt, and verified current or
archived task/execution evidence. Reconciliation is journaled. The legacy
write-before-admission state is repaired only when the old archived allocation
still blocks the exact submission and no newer claim, receipt, worker artifact
or live allocation exists.

Relaxation panels use the shared bounded scratch-allocation and Molpro `hsdel`
solver retries. The solver sequence changes P-space, then internal expansion
workspace and numerical convergence, with at most two replacements. It requires
the native singularity signature and retains each attempt. It does not change
the CAS, state weights, selected roots, Hamiltonian, separation constraint or
optimizer-step budget. Unknown or exhausted failures remain visible for review.

Recon and production share [native orbital/state printing and extracted
per-output evidence](MOLPRO_NATIVE_ORBITAL_DIAGNOSTICS.md). Each CAS/MRCI output
gets an adjacent `.orbital_audit.json` with occupations, labeled AO coefficients,
state information and source hashes. Evidence completeness is distinguished
from orbital/state continuity.

The target is the **requested eventual ROTD calculation**: chemically adequate
states, accurate interfragment potentials and spin gaps, continuity throughout
the contributing geometries, and affordable sampling. These requirements are
joint. A smooth small-space curve can be inaccurate; a larger space or a lower
total energy is not automatically a better reference. Compare chemically
motivated alternatives and justified correlated references before accepting
cost savings. Orbital occupations and CASCI projections guide hypotheses;
they do not certify an error bound or replace reaction-specific chemistry.

The optimization target is the **minimum chemically adequate active space for
that corrected VRC calculation**. Start with the reacting orbitals and complete
required fragment manifolds. Add an orbital group only when it repairs a
demonstrated deficiency in dissociation chemistry, relevant spin splittings,
the corrected angular potential, or robust sampling. Agreement in absolute
total energy with an occupied-valence CAS is not the selection objective.
A comparison tolerance is a declared diagnostic, not proof that its comparator
is exact. Exceeding an assistant-chosen cutoff at one geometry does not, by
itself, establish that an expansion is necessary.

Report the requested reactive surface and correction's spin anchor separately
from the maximum across every state-averaged root. State averaging can be needed
to retain a balanced fragment manifold without making every averaged root an
independent flux surface. Conversely, do not silently omit excited surfaces
that the physical kinetic model includes. The state/flux model must state this
choice explicitly; a CAS size cannot determine it.

The sulfur examples and their native evidence are documented separately in
[the local ORCA validation record](SULFUR_ORCA_VALIDATION_2026-10-15.md).
Molecule names, candidate orbital groups, bonds, geometry constraints, state
averages, tolerances and site resources belong in input data, not runtime
branches. The [generality audit](GENERALITY_AUDIT_2026-10-15.md) records the
scope and limits of the earlier implementation review.

`vrc_tst_accuracy_preset: "anl1-f12"` specifies energy levels and correction
settings. It requires an explicit active-space policy or an orbital/surface
model; it no longer implicitly chooses full valence. A user can still request
`full_valence` deliberately. Installing a selected model or per-bond/stage
selection supplies this requirement without introducing a separate launch path.

## Decisions and completion criteria

Keep these decisions separate in the review record:

| Decision | Evidence needed | What it authorizes |
| --- | --- | --- |
| Chemical candidate | Fragment terms, radical/bond orbitals, resonance partners, required degeneracies and orbital character | A finite comparison calculation |
| Numerically usable model | Native convergence of every requested root, correct active/core partition, independent starts and return checks | Comparing its correlated energies |
| Locally recommended CAS | Matched interaction and angular tests against justified expansions or independent references, cost and limitations | Selecting a candidate for the target VRC calculation |
| Production backend confirmation | Actual method/basis/core/shift conventions, native orbital ordering, state identity and scan behavior | Installing the reviewed backend records |
| Completed kinetics | Converged corrected scans and Monte Carlo flux, audited thermochemistry, executed MESS | Reporting rates and their tested uncertainty |

A local recommendation need not wait for a complete kinetics calculation.
Conversely, a local CASPT2 result does not approve CASPT2-F12 automatically.
An F12 orbital basis used with conventional CASPT2 does not make the
correlation treatment explicitly correlated. Different implementations must
be compared first at matched geometry, active orbital content, state averages,
frozen core, shifts and correlation definition; then test the production
method and basis changes explicitly.

## What accuracy means for this workflow

Reconnaissance is an aid to model development, not an automatic discovery of
the chemically correct CAS. AVAS generates candidate orbital subspaces from
specified atomic valence targets. The user reviews their chemical content,
state averages and dissociation limits; occupations alone cannot make that
decision. Native convergence and agreement between initial guesses establish
numerical reproducibility, not physical accuracy. A larger CASSCF space is a
comparison model, not an exact reference for the correlated potential.

The relevant quantity is the interaction potential for each physical state,
relative to a consistent separated-fragment reference. For corrected sampling,
compare the **final** potential, schematically
`V_final(R,Omega) = V_sampling(R,Omega) + Delta(R)`, with independent correlated
checks. A radial correction can remove radial bias, but cannot generally repair
an incorrect angular dependence, a missing electronic manifold or discontinuous
root/orbital switching. Physical crossings must not be smoothed away. The aim
is accurate state-consistent surfaces over geometries contributing to the flux,
then the cheapest numerical representation that preserves those surfaces.

CAS(2,2) for methyl + methyl has a chemical justification: the two radical
orbitals span the forming bond and its antibond. In CASPT2 or MRCI this does not
mean only two electrons receive correlation; external excitations recover
dynamic correlation beyond the CAS reference, subject to the chosen core policy.
Bare CASSCF(2,2), correlated CAS(2,2), and a corrected VRC sampling model therefore
have different accuracy. The [alkyl-recombination study](https://doi.org/10.1039/B515914H)
validated corrected CASPT2/cc-pVDZ sampling against direct higher-level
CAS+1+2/aug-cc-pVTZ sampling. A [recent implementation](https://pmc.ncbi.nlm.nih.gov/articles/PMC12573747/)
explicitly uses (2e,2o) for methyl + methyl, together with a high-level correction.
Neither result assigns a universal kcal/mol error to CAS(2,2), or establishes
that it suffices for resonance-stabilized radicals or open-shell atomic terms.

The user-supplied passage from
[Knepp et al., PCCP 2007, 9, 4315–4331](https://doi.org/10.1039/B705934E)
provides a particularly direct statement of the chemical selection strategy.
It reports CAS(2,2) for alkyl recombination, expansion to CAS(4,3) and state
averaging for CH3 + OH, and inclusion of the relevant complete pi manifold
for resonance-stabilized radicals. For alkyl + O2, its CAS(7,5) contains the
alkyl radical orbital plus all four O2 pi/pi* orbitals; CAS(3,3) produced
unphysical maxima. These are examples supporting chemically complete small
spaces, not a molecule-to-CAS lookup table or a universal full-valence rule.

The same passage describes inexpensive double-zeta orientational sampling
with an orientation-independent correction from higher-level calculations.
Its conventional TZ/QZ CBS correction removes a modest spurious saddle near
2.2 angstrom in ethyl + O2. Therefore validate the **topology and magnitude of
the corrected interaction potential**, including angular residuals and
fragment relaxation, before judging a basis or space adequate. Conventional
Dunning-basis results do not establish that F12 CBS is required, or that
TZ-F12 automatically suffices. A sparse basis check can decide whether the
correction is converged; separate fragment deformation from a frozen-fragment
interaction comparison. Corrections transferred between chemically similar
systems still require validation for the target event.

Inner and outer models can require different spaces and averages. For example,
the [NH2O study, section 2.2](https://doi.org/10.1021/acs.jpca.3c05181)
uses (3e,3o) for inner NH + OH addition and (5e,4o) outside, where the additional
orbital represents the near-degenerate OH states. These are literature examples
of chemical reasoning, not defaults copied into other reactions. Stationary
well/reaction energies remain the responsibility of the CBH/ANL energy framework;
the VRC calculation must also deliver accurate relative interaction energies.

## Procedure

The protocol must distinguish four kinds of evidence: chemical orbital/state
content, numerical robustness, accuracy of the corrected interfragment
potential, and the resulting kinetics. A completed electronic calculation
establishes only part of the second item. A useful selection report must say
which comparison remains missing and which omitted orbital group that
comparison tests. Expanding every space indiscriminately is not the procedure.

`python -m kinbot.orbital_review assess-vrc OBSERVATIONS.json REPORT.json`
now evaluates explicit, backend-independent potential comparisons. The
observations are hash-bound to native input/output and parsed results. The
schema contains `points` (coordinates and separation), `models` (per-point
`energies_hartree` keyed by matched state labels, coordinates, and `evidence`
path/SHA-256 pairs), and `comparisons` (candidate, comparator, chemical question,
required points, asymptote, angular-reference mapping, tolerance and state
assignment). Model identities and orbital content are input data.

Each model uses one lowest included asymptotic energy as the common zero for
all roots and spins. The report retains raw interaction-potential differences
and computes the residual after applying the electronic correction measured
on the reference cut at the **same separation**. This residual is the error
in angular dependence that a correction depending only on separation cannot
remove. Missing points stay incomplete; a different geometry, missing state
mapping or changed evidence causes an error. The result reports a measured
angular test, not automatic approval of a reference, CAS or production run.
This comparison does not yet calculate the separate geometry-relaxation
correction. Qualitative orbital/state analysis and the actual scan checks
remain essential to the user's choice.

The local ORCA adapter `assess_decision_panel.py` supplies observations to this
same interface. Until cross-geometry electronic-state identities are checked,
it labels its sorted same-spin energies as provisional spectra. Root-index
agreement, or a symmetry-preserving nuclear path, must not be advertised as
proof that orbitals or physical states did not switch. General declared
Cartesian symmetry operations can be audited with `kinbot.orbital_geometry`;
the module contains no molecular point-group choices.

Native Molpro review reports also include `state_order_diagnostics`: the
correlated energy order within each method/control/spin group, the selected
reference-root label, and its energy above the lowest computed root. A root
that was lowest at CASSCF can cease to be lowest after correlation. This
diagnostic does not silently change the production root or assign flux weights.

The same `assess-vrc` input can include `path_tests`. Each entry supplies `id`,
`model`, `state_assignment`, an ordered `path`, `asymptote`, `return_pairs` and
`return_tolerance_kcal_mol`. Optional `spin_pairs` explicitly names `label`,
`target` and `reference` state keys. Return pairs must have exactly identical
nuclear coordinates. Reports retain every state's return error and the
declared spin gaps, and compute secants only between adjacent available path
points. Missing points remain incomplete. These tests measure hysteresis;
they do not infer physical root identities or certify unsampled smoothness.

1. `plan CONFIG PLAN.json` derives fragment/spin context from the existing
   reconnaissance configuration and linked KinBot input. It exposes all
   radial/orientation anchors and selects a small finite comparison panel.
   The initial model is the explicitly supplied chemical candidate, not a
   production choice. A valence reference requires deliberate opt-in; it is not
   the mandatory starting point for a VRC review.
   `rotd_target` records resolved sampling/correction methods and bases, the
   actual scan and dividing-surface grids, temperature controls and available
   surface-model stage overrides. Missing/unreviewed stage choices stay explicit.
   Scan distances are reactive atom-pair distances; outer COM surface distances
   are not interchangeable with them. The linked surface-model hash is checked
   before preparing the plan.
2. Edit/duplicate models to state the chemical purpose, AVAS targets, electron
   and orbital counts, state averages/weights, channel, and region (`diagnostic`,
   `single`, `inner`, `outer`). Retain radical orbitals, the forming/breaking
   bond manifold, relevant resonance partners and complete required degenerate
   fragment manifolds. Orbital indices or near-2 occupations alone do not
   identify a chemically adequate subspace. `available_anchors` lists valid
   names; each model's `anchors` selects its actual validation geometries.
3. `prepare PLAN.json ROOT` creates a fresh campaign with three comparisons:
   CASCI with fresh AVAS, optimized CASSCF with fresh AVAS, and optimized CASSCF
   with checkpoint initialization. Fresh AVAS explicitly starts from converged
   RHF/ROHF record 2100.2, rather than Molpro's standalone atomic-density guess.
   Checkpoint calculations optimize orbitals; they are not frozen CASCI scans.
4. `--reuse-from OLD_CAMPAIGN` seeds matching checkpoint comparisons from
   completed old WFU records. Basis, atoms, coordinates, electron/active-space
   counts, AVAS targets, input hash and WFU hash must agree. Unavailable seeds
   fall back to radial propagation or a fresh start. Old energies are never
   relabeled as optimized results. Keep the old campaign at its original path.
   In particular, the variant name `casscf_propagated` does not guarantee
   neighboring-geometry continuation: a matching external CASCI checkpoint
   overrides its radial predecessor. Reports expose the actual initialization
   as `fresh`, `external_checkpoint` (with path/hash), or `predecessor`.
5. Run with `orbital_campaign run`. Independent systems share one Slurm job
   limit; pilots run concurrently, then all ready full-panel probes can launch
   together after the pilots pass. Full panels do not repeat a redundant
   one-job startup check. Genuine checkpoint predecessors must still finish
   and validate before dependent probes launch. Live allocations, newly
   submitted IDs not yet visible in Slurm, and missing allocations awaiting
   validated output all consume capacity. The finite panels do not trigger automatic spatial
   expansion. Legacy adaptive refinement remains an explicit opt-in; both CLI
   defaults are now zero passes. For a stopped mixed-method campaign,
   `run ROOT --probe-method casscf --max-jobs 20 --max-refinements 0` resumes
   only the optimized references, rechecking/reusing their completed pilots.
   It refuses live or uncertain unselected allocations, re-audits unselected
   controls without submitting them, and reports
   `selected_complete_review_required` rather than whole-campaign completion.
   `progress` labels each system with `selected`; a later run without the
   method filter explicitly resumes the entire original campaign.
   `stop-controller ROOT` stops only the verified Linux controller process
   and waits for its campaign lock to release, preserving running Slurm jobs.
   A new controller adopts their recorded job IDs. This differs from
   `run_cleanup --stop`, which also cancels the scoped calculations.
   If a submitted probe fails native execution or scientific validation,
   the shared controller stops new submissions, continues monitoring all
   already submitted campaign jobs (including other phases), then reports
   all observed failures together. `progress` shows `draining_after_failure`
   while this is happening, with failed job names, IDs and full directories.
   The top-level `error` gives the first failure while draining;
   `progress ROOT --failures-only` provides a compact list without the full
   system inventory.
   Successful remaining jobs are retained; failed energies remain rejected.
   Unsubmitted jobs and descendants of failed checkpoints are not launched.
   The final campaign status is `failed`, requiring review before a repair.
   Uncertain submissions and scheduler errors still stop immediately, since
   allocation ownership cannot safely be inferred in those cases.
   A separate controller-side recovery handles a retained `scontrol` record
   identifying `FAILED`, `ExitCode=0:53` (Slurm batch I/O initialization), with
   the matching job ID/work directory and a specific batch host. After the
   60-second publication grace, it retries at most twice, only when there is
   no native input/output, execution record, nonempty stdio, or current
   retry-guard journal. It journals the old ID and scheduler evidence before
   resubmission. The identified host is excluded across this campaign's
   manifests and later phases, preserving site/user exclusions. No node name
   is built into the runtime. Progress exposes `startup_retries` and
   `slurm_excluded_nodes`. Live/requeued allocations retain their IDs and
   consume the same global cap. Unknown disappearance, cancellation, native
   convergence failure and exhausted retries still require review. This
   mechanism needs no `sacct`; the batch-script retry guard cannot handle
   failures that occur before the shell starts.
   Retained scheduler records are also saved for **started** workers, before
   expensive output/checkpoint audits, so termination reasons can survive
   `scontrol` expiry even when accounting is disabled. Missing records are not
   evidence of a particular termination cause.

   `run --allocation-retries 2` additionally opts into at most two replacement
   allocations per probe for finalized native I/O failures or interrupted known
   workers whose execution record remains `running`. The latter is labelled
   `interrupted_allocation`, not assumed to be I/O. Interrupted workers wait for
   the publication grace; finalized I/O failures can be retried after queue exit.
   Recovery checks queue/scheduler ownership again, archives partial
   products, and restores original seeds. History and limits survive controller
   restarts. Native nonconvergence, known cancellation/time/memory limits and
   uncertain submissions are not retried by this policy. The existing batch
   guard separately allows two in-allocation I/O retries. With both limits at
   two there can be at most nine native starts across three allocations for
   one probe; reducing either limit is supported. Old `attempts/` archives
   are retained but no longer copied into new retry seed snapshots.

   `run --keep-going` finishes other independent manifests after a failure;
   the failed manifest and its dependent points remain blocked. Progress is
   `running_with_failures`, and unresolved failures still produce a failed
   final status/nonzero exit. Full panels remain gated on successful pilots.
   `retry-failed --include-interrupted --job-id ID` can explicitly recover a
   particular stopped allocation while preserving live peers. Repeat the ID
   flag to select several; replacement IDs are not affected by repeating the
   same command.

   To add comparisons without a second competing job cap, prepare an addition
   **inside** the parent tree, stop only its controller, and use
   `orbital_campaign append ROOT ROOT/ADDITION`. The addition must be unstarted
   with distinct titles. Adoption is idempotent, carries its review metadata,
   and prevents launching the child as a separate campaign. Resume the parent
   to share its cap with preserved live jobs.

   Default Slurm memory now includes Molpro's 300-MW/rank program-overhead
   estimate in addition to the requested stack, then adds 25% headroom.
   Eight ranks with `memory_mw=1000` request `104000M`; explicit `memory_mb`
   remains supported. This estimate does not bound all GA/auxiliary allocations.
   Existing batch files are unchanged; Slurm's `SBATCH_MEM_PER_NODE` can change
   the resource request for their next submission without changing chemistry.
   See [Molpro parallel memory](https://www.molpro.net/manual/doku.php?id=running_molpro_on_parallel_computers)
   and [Slurm environment precedence](https://slurm.schedmd.com/sbatch.html#SECTION_INPUT-ENVIRONMENT-VARIABLES).
   For native convergence failures, `orbital_review repair-convergence OLD NEW`
   prepares only rejected anchors, preserving old evidence. Its default small-CAS
   repair uses quadratic RHF and exact CI without L-BFGS. If that also stalls,
   `--solver so_sci --probe-method casscf` explicitly selects Molpro's combined
   first/second-order optimizer with L-BFGS for RHF and CASSCF, excluding optional
   CASCI controls. Run the new campaign with `--pilot-only` first. Basis, active
   space and state averages are unchanged; success
   still requires native convergence. Repeating the same failed repair controls
   is refused. Direct recon configurations can select `so_sci` or `second_order`
   in `scf_options.solver` and `casscf_options.solver` (`diis` is also available
   for SCF). See the [Molpro SCF](https://www.molpro.net/manual/doku.php?id=the_scf_program)
   and [MULTI](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi) manuals.
   New CASSCF decks explicitly set `ACCURACY,ENERGY` inside MULTI: otherwise
   SO-SCI in Molpro 2024.1 can derive a target 100 times tighter than the global
   `GTHRESH,ENERGY` value. `casscf_options.energy_convergence` sets that local
   target; it defaults to `scf_convergence` (1e-10 when omitted). The optional
   `casscf_options.gradient_convergence` separately sets the orbital-gradient
   criterion; omitting it retains Molpro's gradient default. These options do
   not change the global SCF/correlation thresholds. Explicit targets apply
   equally to fresh and checkpoint-initialized CASSCF. A repair from an older
   implicit-target SO-SCI campaign is allowed once; an identical explicit-target
   repair is still refused. Old input files and validation results are preserved.
   Failure reports include the native MULTI threshold line when available.
   An explicitly reviewed change to the orbital-search tolerance can be made with
   `repair-convergence ... --energy-convergence 1e-8`. This changes only local
   CASSCF `ACCURACY,ENERGY`; SCF/correlation targets, gradient tolerance, basis,
   active space and state average remain unchanged. Molpro documents 1e-8 as
   its first-order energy default, but the target is a stopping criterion,
   not an error bound. Use a pilot and compare converged state energies/character
   before propagating the choice. `--verbosity 4` enables detailed MULTI output
   (`casscf_options.verbosity` in a direct configuration). Increasing print level
   alone does not bypass the guard against repeating failed numerical controls.
6. `report ROOT OUTPUT.json` supplies geometry, chemical rationale, requested
   state averages, occupations, complete native orbital/CI text, root energies,
   elapsed times and hash-bound artifact paths. It compares geometrically
   equivalent calculations across initializations and reports radial energies
   relative to the farthest completed radial point. That point is labeled as a
   finite reference, not automatically certified as an asymptote. Native XML
   exports and the full outputs are included by `orbital_campaign collect`;
   binary WFU files remain on the cluster.
   Reports include both pilots and full-panel attempts. Repeat
   `--additional-campaign OTHER_ROOT` to compare the original panel with repair
   campaigns, without rerunning a successful pilot just to put it in a panel.
   Counts are attempts, not unique geometries. Failed and unsubmitted attempts
   remain listed, and campaign/phase, solver settings and checkpoint provenance
   are retained. Equivalent-geometry comparisons require matching active-space
   definitions and state averages as well as atoms, basis, method and root labels.
   Radial profiles remain separate by campaign and phase; combining evidence
   does not silently splice unlike runs into a potential curve.
   `potential_profiles` additionally joins explicitly named comparison models
   across single-geometry systems, keeping basis, stage, orientation, space,
   state average and correlation controls separate. It supplies finite-reference
   energies, root-labelled secant slopes and timings. Duplicate distances or a
   single isolated point receive no manufactured interaction reference. These
   are inspection data: slopes alone neither certify continuity nor accuracy.
   Collect each stopped campaign first to re-audit historical completion labels,
   then generate the combined report. Collection and reporting hold the campaign
   and manifest locks, reject uncertain submissions, and check recorded job IDs
   against Slurm before reading native evidence. Scheduler-query failure aborts
   the export instead of treating it as an empty queue.
   A failed initialization comparison need not prevent inspection of a converged
   alternative at that geometry. It does leave initialization robustness
   unresolved; it is not a reason to approve a failed result or a production CAS.
7. After orbital/initialization review, duplicate a plan into a fresh small
   correlated validation campaign. Each model may request `correlation_checks`,
   for example `[{"method":"caspt2-f12(5,4)","shift":0.3,"ipea":0.25},
   {"method":"mrci-f12+q(5,4)"}]`. Counts must match that model. Checks run
   only with optimized CASSCF, for every requested spin/root, and use the
   production selected-state and Davidson-energy selectors. They are disabled by default
   on the larger reference spaces. Inspect interaction energies, spin gaps,
   intruder/shift sensitivity and rate sensitivity before a final selection.
   For an explicit partition of reviewed natural orbitals, use
   `prepare-partitions SPEC.json ROOT`. Each `models` entry supplies `name`,
   `model` (comparison group), `source_directory`, `region`, `chemical_note`,
   target `active_electrons`/`active_orbitals`, `close_lowest`, `discard_highest`,
   and `correlation_checks`; optional `state_average` defines a separate average.
   Each `state_average` entry accepts a positive `block_weight` (default 1)
   in addition to `multiplicity`, `roots`, and `weights`. The `weights` sum
   to one **within** that multiplicity; the effective normalized weight of
   root i in block b is `block_weight[b] * weights[b][i] / sum(block_weight)`.
   Omitted block weights retain equal total weight per spin block. A declared
   `reference_policy: "state_average"` in a single-spin correlated call uses
   the within-block weights only. These numerical orbital/reference averages
   do not set VRC surface populations or spin statistical factors.
   To preserve these reviewed roots and weights in production, select with
   `state_average_policy: "reviewed"` (the selection command's default).
   The alternative `channel_derived` policy deliberately rebuilds equal-weight
   averages from the accepted fragment channel and does not retain custom
   within-block or block weights. It is not appropriate when transferring a
   specifically validated weighted average.
   Optional `orbital_swaps: [[i,j], ...]` exchanges explicitly reviewed orbitals
   **before** closing/discarding the ends of the source-active block. Indices
   are one-based within that source block, not absolute molecule-wide indices.
   Exchanges are applied in order using MULTI's documented
   `ROTATE,orb1.sym,orb2.sym,0` following `START`. Optional
   `source_wfu_sha256` binds that identity decision to the previously reviewed
   checkpoint and rejects replacement sources. Neither occupations nor a CAS
   count establishes orbital character; the user supplies these choices.
   See [MULTI initial-orbital exchanges](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#rotating_pairs_of_initial_orbitals).
   A direct `molpro_reconnaissance` input can instead supply fixed `anchors`
   with `name`, `geometry_angstrom` (N-by-3, in the original atom order), and
   `orbital_chain`. With `orbital_initialization: "propagated"`, points with
   the same chain name consume the previous point's audited target-space
   checkpoint, separately for each basis. This supports exact angular and
   forward/return comparisons without reconstructing geometry from a radial
   coordinate. Explicit coordinates cannot also request bond displacement or
   rotation; these fixed chains cannot use adaptive `scan_policy`.
   Every chain starts with its own hash-bound `initial_records` entry, keyed
   `BASIS_NAME__ANCHOR_NAME`. A source partition's swaps apply only at that
   initial point. Subsequent points reuse the target CAS and do not repeat
   the source swaps. Campaign pilots cover each chain and basis.

   An initial record normally contains `path` and `sha256`. Optional `basis`
   declares its source orbital basis. If that differs from the destination,
   `project_orbitals: true` is required and emits `START,...,OVL`; an explicit
   projection always requires a source-basis declaration. The checkpoint is
   then a **projected starting guess**, not an approved destination-basis
   model. Reoptimize, inspect the native destination-basis orbitals and state
   manifold, and bind production selection to that new checkpoint. The
   source declaration is provenance supplied by the user, not inferred from
   the opaque WFU. Legacy records without a basis declaration remain accepted
   for existing workflows; this does not certify cross-basis portability.
   Foreign-basis seeds are staged separately on FILE,3 and projected into a
   fresh FILE,2; a dump record cannot be overwritten with another basis dimension.
   Generated C1 RHF references explicitly reset spin and determinant occupations
   rather than inheriting an AVAS/MULTI restart state.
   See [Molpro orbital projection](https://www.molpro.net/manual/doku.php?id=general_program_structure#orbital_directives).
   Optional `casscf_options` overrides the partition preparation's solver defaults;
   optional `correlation_options` sets `max_iterations`,
   `max_reference_iterations`, `max_reference_vectors`, and/or `max_external_vectors` for correlated calls. These limits must be
   positive integers. `max_reference_iterations` is Molpro's internal/reference-CI `MAXITI`,
   distinct from CASSCF macroiterations and from the correlated `MAXIT` limit.
   `max_reference_vectors` is Molpro's internal Davidson expansion capacity
   (`OPTION,MAXVI`). `max_external_vectors` controls the separate external
   Davidson subspace (`OPTION,MAXDAV`); increasing MAXVI does not increase MAXDAV.
   Both are independent of the iteration ceiling. These controls change solver budgets without changing convergence
   tolerances, the Hamiltonian, or the selected states. Omitted correlation
   controls retain native defaults. See the
   [Molpro MRCI iteration controls](https://www.molpro.net/manual/doku.php?id=the_mrci_program#maximum_number_of_iterations),
   which also apply to [RS2](https://www.molpro.net/manual/doku.php?id=multireference_rayleigh_schroedinger_perturbation_theory).
   MRCI alone additionally accepts `correlation_options.denominator_shifts`,
   e.g. `{"pairs":0.1,"singles":0.1,"internal":0.1}` in hartree, rendering
   SHIFTP/SHIFTS/SHIFTI. Only supplied keys are overridden; finite nonnegative
   numbers are required. These alter Davidson update denominators, not the
   CASPT2 level shift or the intended converged MRCI eigenproblem. They are
   optional numerical experiments, never automatic infrastructure retries.
   No default changes. See [Molpro denominator shifts](https://www.molpro.net/manual/doku.php?id=the_mrci_program#denominator_shifts).
   Separately, `correlation_options.reference_policy` accepts `selected_root`
   (the default) or `state_average`. The latter adds `REFSTATE` for all roots
   of the requested spin in the reviewed CASSCF average; CASPT2 also receives
   that spin's weights, so its Fock operator uses the reference-averaged
   density. A check can override this policy with its own `reference_policy`.
   For MRCI-F12, selecting an excited root requires `state_average`: the
   renderer supplies the full same-spin reference manifold and matching
   `NSTATI`. A lone excited reference is rejected before submission because
   the installed F12 solver aborts when its reference/internal counts differ.
   This changes the correlation treatment, not just the solver effort. It
   does **not** enable XMS-CASPT2 or guarantee rotational/fragment invariance.
   Test its asymptotes and interaction potentials independently. The same
   explicit policy is retained through selection into production.
   `correlation_options.state_policy` accepts `independent` (unchanged default)
   or `simultaneous`. The latter requires MRCI and `reference_policy="state_average"`;
   all declared same-spin roots are calculated together with matching STATE,
   REFSTATE and NSTATI. REFSTATE alone does not prevent separately optimized
   roots from reaching the same solution. Explicit NOSWAP gives variational
   MRCI energy order; selected root i reads relaxed-Davidson `ENERGD(i)` and
   uses the matching numerical-gradient slot. This is correlated eigenstate
   rank, not guaranteed continuity of a CASSCF root's character through a crossing.
   Recon records that convention and computes one common call per spin.
   Validate the native spectrum and root character before production.
   Use separate stage models: this MRCI-only policy is rejected for CASPT2.
   For MRCI-F12 only, `correlation_options.cabs_singles` accepts integer 1
   (singles CI, unchanged default), 2 (perturbational), or 3 (coupled).
   This is a method choice requiring consistent reference curves and matched
   relative-energy checks, not an automatic I/O retry or numerical override.
   It is rejected for CASPT2 and non-F12 MRCI. The relaxed-Davidson selector
   still reads the final `ENERGD(1)` for an independent call, or the requested
   `ENERGD(i)` for simultaneous roots, after the entire command completes.
   An optional `{"method":"uccsd(t)-f12b","multiplicity":4}` check supplies
   a separate, explicitly requested ground-state spin anchor. The value 4 is
   only an example. Such checks must follow the MR checks; they run a new
   restricted HF reference in a separate orbital record, preserving the
   CASSCF orbitals. They validate native F12b plus unscaled triples energies
   and reject failed SCF/CC calculations. They cannot supply arbitrary
   excited roots or replace a multireference low-spin surface. Their
   usefulness depends on the reference character at the geometry examined.
   Source CASSCF convergence, hashes, geometry and basis are checked. Partitions
   reoptimize all closed and active orbitals and retain explicit source provenance.
   They cannot change basis. Optional `target_distance_angstrom` generates a new
   point on the source's declared dissociation coordinate, retaining its fragment
   orientation and reoptimizing from the reviewed seed. Source/target geometry
   and the checkpoint hash are retained separately. This permits small stencils
   in the actual requested scan region without repeating a large-space recon.
   Alternatively, `target_geometry_angstrom` supplies an explicit N-by-3 array
   in the source atom order and coordinate frame. It is mutually exclusive with
   `target_distance_angstrom`; no source-anchor translation or rotation is
   reapplied. This supports a fixed-geometry comparison with another backend,
   including changed internal fragment coordinates, while keeping the original
   seed provenance. The new CASSCF must converge and its orbital character must
   be checked; transferring a seed never transfers acceptance of the old result.
   Each entry is a single-geometry system, so run this
   finite campaign with `--pilot-only` to avoid repeating its panel duplicate.
   Correlated checks write/select current pseudo-canonical orbitals explicitly
   for CASPT2 and natural orbitals for MRCI; stale WFU orbital sets are not used.
8. The user chooses the final space, orbital ordering and states. Register the
   actual reviewed WFU with `molpro_reconnaissance select --basis ...`; supply
   the reviewed counts, average and root choices. The command does not reorder
   orbitals or certify chemical accuracy. New selections carry
   `state_average_policy="reviewed"`, so production preserves their averages
   and rejects conflicting count overrides. Old files without this field retain
   the legacy channel-derived policy; recreate them to preserve a reviewed
   average explicitly.
   `select --casscf-options-json '{...}'` records the solver and tolerances
   validated during recon in the selected state model. The general model field
   `casscf_options` supports `solver` (`so_sci`/`second_order`),
   `energy_convergence`, `gradient_convergence`, `max_iterations`,
   `exact_diagonalization`, `quasi_newton`, and `verbosity`. Production correction,
   sampling and geometry calculations share this renderer, so the selection need
   not silently revert to numerical controls that failed during recon. Omitting
   the field retains legacy defaults; choosing tolerances still needs numerical
   validation, independently of active-space accuracy.
   `select --correlation-options-json '{"max_reference_iterations":200}'`
   similarly retains explicitly validated MR iteration limits in the general
   `state_model.correlation_options` field. This is an example of the syntax,
   not a universal convergence prescription. Native solver failure remains
   an error even when energies are printed; `IGNORE_ERROR` is not enabled.

For failed jobs, `execution.json` is published before checkpoint staging, so
missing/changed seeds and input-validation failures are retained. If a job
never produced native output, audit reports the missing output and available
batch/retry stderr rather than diagnosing missing occupations as CAS failures.
Only the recorded allocation's latest restart journal is used. A scheduler
interruption needs retained `scontrol` evidence or Slurm accounting to establish its cause; a stale
`status="running"` execution record is not proof of native nonconvergence.
Campaign `collect` includes `scheduler_accounting.json`, queried with `sacct`
for only the recorded allocation IDs and their prior attempts. This preserves
states, exit codes/signals and restart counts alongside native diagnostics.
An unavailable accounting service is recorded without blocking the export.

`orbital_campaign monitor ROOT --details` reads one live `squeue` snapshot and
bounded tails of current native/launcher/guard files without auditing, hashing,
locking or modifying the run. It shows actual Slurm exclusions and requested
memory, current stage/iteration, native CI parameters, output age and detected
diagnostics. `--json` exports the same snapshot. A log match is not a terminal
failure verdict; quiet output alone is not a stall. Requeued guard journals are
selected by restart count; old archived attempts are not scanned as new errors.

Exclusions are passed as actual `sbatch --exclude` arguments. `SBATCH_EXCLUDE`
is not a Slurm input variable; KinBot translates that older instruction for
compatibility. After stopping only the controller, use
`orbital_campaign refresh-exclusions ROOT` to apply the campaign's recorded
exclusions to its pending allocations, preserving existing exclusions and
verifying Slurm's response. Running allocations remain in place.
Repeat `--node HOST` to persist additional user-specified exclusions across
all campaign panels. After repairing the infrastructure cause,
`orbital_campaign retry-startup-failures ROOT` explicitly recovers stopped
pre-shell I/O failures, including exhausted ones. It requires the same
identity, signal, missing-worker-artifact and publication-grace evidence as
automatic recovery. Old histories remain intact; recovered jobs receive two
new automatic startup retries. Native solver failures, unknown disappearances,
cancellations and running jobs are not selected by this operation.

For a diagnosed reference-CI failure, an explicit in-place numerical retry is
`orbital_campaign retry-failed ROOT --job-id ID --max-reference-vectors N --max-reference-iterations M`.
The controller must be stopped and the selected allocation must have ended.
It archives the old attempt, preserves the prepared input and original seed,
and records the numerical override in the manifest and hash-bound executed
native input. Completed and unrelated running jobs are preserved. The override
changes internal Davidson capacity and/or iteration ceiling; convergence requirements, active
space, states, shifts and reference policy are unchanged. Acceptance still
requires native convergence. Carry a successfully validated value into the
reviewed production model's `correlation_options`.

Campaign execution settings can be supplied in `ROOT/campaign_settings.json`
(loaded automatically) or with `--settings FILE` on `run`, `refresh-exclusions`,
`retry-startup-failures`, and `retry-failed`. Example:

```json
{
  "schema": 1,
  "max_jobs": 20,
  "max_refinements": 0,
  "pilot_only": true,
  "probe_method": "casscf",
  "keep_going": true,
  "allocation_retries": 2,
  "startup_retries": 2,
  "io_retries": 2,
  "memory_mb": 104000,
  "excluded_nodes": [],
  "solver_options": {
    "scf": {"energy_convergence": 1e-8, "max_iterations": 2000},
    "casscf": {"energy_convergence": 1e-8, "gradient_convergence": 1e-5, "max_iterations": 999},
    "correlation": {
      "energy_convergence": 1e-7, "coefficient_convergence": 1e-5,
      "max_iterations": 2000, "max_reference_iterations": 2000,
      "max_reference_vectors": 60
    }
  }
}
```

These are editable execution settings, separate from the immutable scientific
input snapshots. Molpro MULTI limits CASSCF `max_iterations` to 999; larger
values are rejected before submission. This limit is separate from SCF and
correlated/reference-CI controls. Unknown keys, noninteger limits, unsafe node names, and retry
counts outside 0–10 are rejected. Explicit CLI values override JSON; boolean
flags also have `--no-pilot-only` and `--no-keep-going` forms. Without settings,
existing defaults remain unchanged. `memory_mb` is the per-node Slurm request;
it must cover all Molpro rank stacks, and the user must keep it within the site's
node capacity. Optional `scratch_root` is an absolute directory path.

`startup_retries` bounds confirmed pre-shell I/O retries; a manual
`retry-startup-failures` grants that many additional automatic retries only to
the explicitly recovered stopped attempts. Existing manual allowances and all
attempt histories persist. `io_retries` controls the shared worker I/O guard;
`allocation_retries` controls stopped-worker/allocation recovery. These limits
do not authorize accepting unconverged scientific results. A controller restart
does not erase consumed retry counts.

`solver_options` applies explicit numerical tolerances and iteration limits to
new submissions and explicit `retry-failed --job-id` targets. Completed results
retain their original settings; running or queued jobs are not rewritten.
Energy thresholds are in hartree. The coefficient tolerance is a norm: `1e-5`
means the literal Molpro `THRVAR=1e-10`, not `1e-20`. MR calculations also receive
local `THRDEN=1e-7`, independent of SCF settings. Fresh generated MR inputs and
production VRC MR inputs use these two defaults unless overridden in the state
model's `correlation_options`. The legacy coupling of global COEFF to the SCF
energy tolerance was a bug. No convergence checks are bypassed. A numerical
tolerance is not an estimate of electronic-structure model error; final choices
still need relative-energy and PES validation. See
[Molpro's THRDEN/THRVAR definitions](https://www.molpro.net/manual/doku.php?id=the_mrci_program).

The older `reference_solver_options` applies to **new submissions with MR correlation
checks**, and to explicitly selected `retry-failed --job-id` targets. It does not
alter running, queued or completed scientific inputs. Explicit per-job retry
limits take precedence over campaign defaults. Prepared inputs remain immutable;
the actual transformed input and its limits are hash-bound in `execution.json`.
Refinement retains these execution settings. The run records effective settings
in `campaign.json`; review archives include the local settings file.

The scientific JSON already provides these controls when preparing a new panel:

| JSON location | Controls |
| --- | --- |
| Recon input `scf_options` | `solver`, `max_iterations` |
| Recon input `scf_convergence` | SCF energy threshold, independent of CI coefficients |
| Recon input `casscf_options` | `solver`, `max_iterations`, `exact_diagonalization`, `quasi_newton`, `energy_convergence`, `gradient_convergence`, `verbosity` |
| Recon input `probe_definition.state_model.correlation_options` | `max_iterations`, `max_reference_iterations`, `max_reference_vectors`, `max_external_vectors`, `energy_convergence`, `coefficient_convergence`, `reference_policy`, MRCI `state_policy`/`denominator_shifts`, MRCI-F12-only `cabs_singles` |
| Each recon `correlation_checks` entry | `method`, `shift`, `ipea`, optional `reference_policy` |
| Recon input resources | `ppn`, `memory_mw` per Molpro rank, `memory_mb` per node, `partition`, `walltime`, `scratch_min_mb` |
| Partition recipe model | `casscf_options`, `correlation_options` |
| Final KinBot `vrc_tst_state_model` or a reviewed stage model | `casscf_options`, `correlation_options` |

For example, the same correlated limits in a new recon or final production state
model are `"correlation_options": {"max_reference_iterations": 2000,
"max_reference_vectors": 60}`. Chemical definitions, states and accuracy thresholds
belong in these reviewed scientific inputs. They cannot be changed through the
campaign execution file. Successful convergence with these limits is necessary,
but does not itself validate the active-space accuracy or the sampling potential.

For finite independent single-geometry validation panels,
`orbital_review prepare-incomplete OLD_ROOT NEW_ROOT` prepares only incomplete
pilots, preserving completed results in `OLD_ROOT`. Optional
`--casscf-max-iterations` and `--reference-max-iterations` raise iteration ceilings
without lowering existing limits or altering tolerances, chemistry, or seeds.
It verifies that source allocations are stopped and uses the original input
configurations and verified seeds, never interrupted output WFUs. It rejects
multigeometry/dependent reconnaissance chains, which must use their normal
resume machinery. Run the new panel with `--pilot-only`, then combine evidence
using `orbital_review report NEW_ROOT REPORT --additional-campaign OLD_ROOT`.

AVAS/CASCI occupation inspection follows the
[Molpro recommendation](https://www.molpro.net/manual/doku.php?id=automated_construction_of_atomic_valence_active_spaces).
The complete panel is a custom validation aid. The
[NH2O paper, section 2.2](https://backend.orbit.dtu.dk/ws/portalfiles/portal/337661919/nh2_o_revised_final.pdf)
motivates choosing chemical subspaces and state averages by channel and
inner/outer region, rather than applying one global diagnostic average to
every production surface.

## Inner/outer routing and multiple active spaces

The surface model already has `channels[].bond` and
`bottlenecks[].parameters`. Each bottleneck can override
`vrc_tst_active_space` and inline `vrc_tst_state_model`. Thus different channels
and inner/outer regions can use different CAS counts, orbital records and
state averages inside one main input. Each inline model can contain several
basis-specific WFU records. Legacy external state-model files remain readable.

Install each user-selected model into a **new** KinBot input using:

```text
python -m kinbot.orbital_review install-selection INPUT.json SELECTED.json NEXT.json \
    --surface-file SURFACES.json --bond 1 2 --stage inner
```

The command writes only `NEXT.json`, embedding the selection and preserving other
bottlenecks, checks that the selected average can supply that stage's surfaces,
and prevents inherited global root overrides from changing the selection.
Use `NEXT.json` for the next stage selection; omit `--surface-file` to update
its inline surface model. It leaves the
surface model's `reviewed` flag unchanged: selecting orbitals does not approve
electronic statistical weights, energy offsets or equilibration assumptions.
Without a surface model the command installs a global selection instead.

KinBot creates a separate correction/sampling calculation for each selected
stage and spin/root. Inner dividing surfaces use reactive pivots; outer
surfaces use fragment centers of mass. ROTD_py integrates each scalar surface;
KinBot sums the physical electronic contributions with the reviewed energy
offsets and weights. MESS combines inner/outer fluxes using the selected
statistical or minimum model. Correction variants are not extra electronic
states. Correction scans retain their atom-pair coordinate even when ROTD's
dividing surfaces use COM separation.

This route applies to accepted `hom_sci` channels without `do_vdW`. Two
variational bottlenecks do not require an explicit stationary vdW well.
The IRC-based vdW identifier is therefore not extended for this workflow;
its homolytic placeholder is the bonded parent, not a candidate complex.
No radial-minimum classifier or automatic well insertion was added.
An explicitly retained intermediate-well network is a distinct kinetic model
from MESS's reduced two-bottleneck treatment.

The existing high-level examples retain CBH/ANL, ROTD, MESSPF and MESS requests.
Installing a selection preserves those controls. A licensed native end-to-end
run is still required to validate the entire chain on Blodgett; local tests
cannot establish native convergence, active-space accuracy or rate accuracy.

## Geometry and correction protocol

Backend confirmation and production correction scans have different purposes.
After selecting a space locally, confirm it in the production backend on the
same fixed geometries, with the same chemical orbital subspace, state average,
frozen core and numerical settings. A new geometry optimization would obscure
this comparison. Include an asymptote, representative interaction-region
points, an angular perturbation and a return-path check; add reference-space
comparisons only where needed to confirm the selection. Native CASPT2-F12
confirmation is not a second campaign to build the production correction.

KinBot already owns the production constrained relaxed scans and the matching
frozen-fragment geometries (`VTS._do_molpro_multireference_scan`). The geometry
objective is the configured `vrc_tst_scan_method` / `vrc_tst_scan_basis` and
selected electronic root. With `spin_splitting`, KinBot evaluates the CC
reference and each multireference target/reference pair on those relaxed
geometries, assembles the corrected target surface, and subtracts the sampling
potential evaluated with frozen fragments. RotdPy consumes that correction
potential during sampling. This includes geometry relaxation at the configured
scan level; it does not optimize the composite spin-splitting energy itself.
See `kinbot/vrc_correction.py:correction_components` and
`assemble_spin_splitting_reference` for the explicit energy composition.

Use a matched frozen-fragment geometry panel to isolate the electronic model.
When preparing a new reference geometry, an alternative to separate fragment
optimizations is to optimize the
separated pair at a sufficiently large fixed separation, constraining the
intended nuclear symmetry and relative orientation. Verify that residual
interaction is negligible for the required accuracy; a nominal distance alone
does not establish this. Freeze those optimized internal geometries and vary
separation for the rigid reference cut. Use the same geometries for every
candidate/reference comparison. Add angular perturbations and radial neighbors
in the region contributing to the requested VRC flux.

This separated-pair optimization does **not** supply the finite-separation
relaxation correction. For a method H and sampling method L, retain separately

```
Delta_electronic(R) = V_H(R, frozen) - V_L(R, frozen)
Delta_relaxation(R) = V_H(R, relaxed) - V_H(R, frozen)
V_corrected(R, Omega) = V_L(R, Omega, frozen)
                      + Delta_electronic(R) + Delta_relaxation(R)
```

All potentials use the declared consistent asymptotic references. The production
KinBot scan supplies the finite-R relaxation contribution; computing extra
relaxed points during active-space reconnaissance is optional and should answer
a specific unresolved selection question.
Check whether each radial correction transfers to relevant off-axis
orientations. If angular residuals are too large, improve the sampling model
or use a correction that resolves that dependence. Do not hide the discrepancy
by assigning an independent zero to each electronic root.

State explicitly which energy defines each geometry optimization: a physical
root or a deliberately chosen manifold average. A state-averaged orbital
optimization does not by itself imply a state-averaged geometry objective.
Preserving nuclear C3v/Cs/etc. with tied coordinates is distinct from imposing
electronic irreducible representations. The latter and electronic state
identity must be checked separately. A spin-free symmetry-constrained geometry
also need not equal a spin-orbit vibronic equilibrium geometry.

For Molpro, use a Z-matrix with tied variables and an explicit `ACTIVE` list
when preserving a constrained nuclear path. With the RF/Z-matrix optimizer,
only the listed variables change. An arbitrary energy objective, including
Davidson-corrected MRCI, requires `OPTG,VARIABLE=...` and numerical gradients.
PQS optimization of a Z-matrix can replace tied parameters, so it must not be
used as a substitute for these restrictions. See the
[Molpro OPTG manual](https://www.molpro.net/manual/doku.php?id=geometry_optimization_optg).
Use the documented iteration limits and verify final geometry, objective,
convergence and orbital character before materializing the frozen scans.

## Reproducible evidence and production handoff

The review record should contain the target bond/stage, chemical orbital
content, state labels and weights, original/native input, program version,
basis coefficients, core policy, shifts, solver controls, seeds and orbital
ordering, geometry provenance, output, parsed energies and comparison reports.
Bind artifacts by SHA-256 and retain unsuccessful attempts as rejection or
recovery evidence in an archive. Production tables should contain accepted
observations only, with links to any limitations affecting their interpretation.
Do not delete the sole copy of source checkpoints when cleaning a review tree.

After native confirmation, install each bond/stage selection using the command
above. Keep the requested CBH/ANL thermochemistry, correction, sampling and
MESS settings in the parent input. The spin-splitting correction combines
method estimates of the **same physically assigned target/reference gap** with
a validated single-reference high-spin anchor; it does not average arbitrary
spatial roots. A ground-state coupled-cluster calculation cannot stand in for
an excited high-spin root. Orbital identity, electronic statistical weights,
spin-orbit offsets and equilibration assumptions are separate review items.

### Shared CASPT2 references and a ground-state high-spin anchor

Two independent choices must be recorded. `correlation_options.reference_policy`
set to `"state_average"` emits `REFSTATE` for the entire declared same-spin
manifold and its normalized `WEIGHT` values in every CASPT2 call. This gives
the roots a common reference density in the zeroth-order Hamiltonian. CASSCF
state averaging alone does not make this CASPT2 choice. The native sulfur
controls show why it matters: the apparent atomic-component splitting with
separate IPEA references disappears with equal same-spin reference weights,
while genuine fragment excitation gaps remain. This is shared-reference,
independent-root CASPT2, not an implicit switch to MS/XMS-CASPT2.

For a ground-state high-spin CC anchor, a CASPT2 root number need not remain
the lowest correlated energy at every geometry. The general state-model
option `"energy_selection": {"4": "lowest"}` evaluates every declared quartet
root with those same shared-reference settings and returns their minimum.
Multiplicity 4 is an example, not a hardcoded spin. Each value is retained
in the native `KBGAPS` vector; `KBVRCSEL=MIN(KBGAPS)` is the returned energy.
The array is reset for each evaluation, so values restored from a checkpoint
cannot enter the minimum. A failure in any requested root invalidates the
evaluation; the code does not select only the successfully computed subset.
This policy requires CASPT2, shared references, and `selected_roots` equal to
1 for that multiplicity. Remove the `energy_selection` entry to request a
particular fixed root instead. Geometry derivatives require a fixed root:
the calculator rejects lowest-manifold gradients rather than differentiating
the last root while reporting a different energy.

Use this option when reconnaissance establishes the relevant manifold and
the intended reference is its ground state. It handles energy reordering
without a molecule-specific root switch. It does not discover missing states,
prove the CC solution's character, repair a changing active subspace, or
guarantee smoothness through every crossing. Confirm those separately with
the reconnaissance checks. All three example production models enable it for
the quartet correction reference. CH3S demonstrates why: the native center
point has its lowest quartet CASPT2 energy at root 3; other reviewed points
have root 1. A one-root manifold incurs no extra correlation call.

These controls are available for any reviewed active space. They are not
silently imposed on arbitrary excited-state or diabatic calculations.
Native syntax: [Molpro CASPT2 reference controls](https://www.molpro.net/manual/doku.php?id=multireference_rayleigh_schroedinger_perturbation_theory)
and [indexed variables and MIN](https://www.molpro.net/manual/doku.php?id=variables).

The sampling basis may differ from the correction basis. Specify it explicitly
and validate the corrected potential; selecting an accuracy preset alone does
not establish that it matches the user's requested sampling cost. Likewise,
state the requested `rotdpy_flux_parameters.flux_rel_err` and sample ceilings:
a sample ceiling is not evidence that the uncertainty target was achieved.
Production mode must retain completion/thermochemistry/flux checks through
MESS. Passing parser and orchestration tests is not a native end-to-end result.

## Single-input production selection (current interface)

After reviewing orbitals, supply their **one-based C1 indices in the exact
source record**, not indices from another basis or from an unbound SCF run.
`molpro_reconnaissance select` accepts `--active-indices 3 6 8`, for example,
with `--active-electrons 4 --active-orbitals 3`, the checkpoint, record, basis,
reviewed state averages and solver options. The numbers here illustrate an
interface; they are not a chemical recommendation for an arbitrary molecule.
`--orbital-type natural` identifies the intended set in a record that may hold
multiple orbital types. Each basis has its own source record and index list.

`orbital_review install-selection INPUT SELECTION OUTPUT` now embeds the
whole selection in `OUTPUT` as `vrc_tst_state_model`; the selection JSON is
not a runtime dependency. For distinct reaction spaces, supply `--bond i j
--stage single` and, for a scoped channel, `--reaction EXACT_REACTION_NAME`:
the model is embedded under
`vrc_tst_surface_model.channels[].bottlenecks[].parameters`. An existing
external surface definition can be imported with `--surface-file`, but the
resulting output is still one JSON. No second surface JSON is written.
The user can also edit these fields directly in the final input.

Within each inline `basis_models[BASIS]`, `orbital_selection` contains:

```json
{
  "mode": "molpro_record",
  "path": "orbitals/reviewed.wfu",
  "sha256": "REPLACE_WITH_THE_SHA256_OF_THE_REVIEWED_BINARY",
  "record": "2140.2",
  "orbital_type": "natural",
  "active_indices": [3, 6, 8],
  "align_active_orbitals": false
}
```

This illustrative fragment must be paired with its actual electron/orbital
counts and state model. Positive distinct indices and list length are checked;
file bytes are checked against the SHA-256. Basis lookup, molecular atom order
and charge are checked when resolving the selection. Relative binary paths
belong to the containing input and become absolute in generated worker inputs.
The binary is required electronic-structure data, not another configuration
file. It cannot be replaced by naked MO indices or a different molecule's WFU.

For N electrons and CAS(n,m), KinBot computes `CLOSED=(N-n)/2` and
`OCC=CLOSED+m`. It packs the selected indices into that active block using
[Molpro's documented zero-angle ROTATE exchanges](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#rotating_pairs_of_initial_orbitals).
The input list gives the initial order. The remaining orbitals follow the
resulting permutation; choosing orbital character and the appropriate closed
space still requires review. Repacking with `align_active_orbitals=true` is
rejected because DIAB would refer to the source's different active block.
Use the packed, validated checkpoint before enabling DIAB. The `select`
command defaults alignment to false; legacy direct models retain their prior
default when the field is absent. Production records explicitly set it.

After a completed optimization, the next checkpoint model records the packed
indices, natural type, new hash, and the same alignment choice. Source swaps
are therefore not repeated. Method labels are materialized with the selected
CAS counts. Changing counts alone does not choose orbitals: the index list,
checkpoint and state manifold must agree. `state_average` roots are required;
omitted weights default equally within a spin block and block weights default
to one. `selected_roots` defaults to root one, and all numerical options remain
user overrides. The sulfur examples explicitly retain their validated solver
options and shared correlation references; those choices are data, not
molecule-name branches in the implementation.

The two production examples now show the full interface in their main JSON
files, including all six aDZ/aTZ index lists. Their old duplicate state/surface
configuration files were removed. Historical native evidence is retained.

## Distinct reactions sharing a separation coordinate

A channel may include `"reaction": "EXACT_REACTION_NAME"` alongside `bond`.
Use this when different accepted reactions use the same atom pair as their
separation coordinate. For example, a homolysis and a rearrangement's product
complex can have the same two anchor atoms but different fragment compositions.
The reaction name is an exact routing key, not a pattern or a chemistry rule.

Two channels on one bond require distinct explicit reaction names. Mixing a
bond-only channel with scoped channels on that bond is refused. Unambiguous
legacy bond-only models remain valid; a collision in accepted reactions or
current review requests requires explicit scoping. Additive publication cannot
silently rename or rescope an installed scientific model.

Selection, publication, worker acknowledgement, scan dispatch and reference
probes use the same reaction identity. Both endpoint fragments must also match
accepted products in connectivity, stereochemistry, charge and spin before
fragment QC. A fragment-identity hold requires the owning scan driver to verify
the endpoint; installing an orbital model cannot acknowledge it away.

## Finite character panels and a code-refresh handoff

SCF-only recipes can set `orbital_character_only: true`,
`orbital_print_virtuals: 8`, and `orbital_print_threshold: 0.05`. They require
at least two distinct distances per dissociation coordinate, retain native
checkpoints/exports, and do not trigger potential-coverage refinement from SCF
energies. Coordinate-generated CAS radial screens require at least four distances;
explicit candidate anchors can define a smaller diagnostic comparison. Product-exit
`geometry_source` may set `include_reference_geometry: true` to include the exact
accepted product-complex geometry alongside the explicit translated anchors.
Fragment internals are preserved during translation.

`requested_orbital_recon batch ROOT PLAN OUTPUT --shared-pool --preflight`
validates the complete plan in temporary storage without staging/submitting it.
Add `--resume` to verify an identical saved campaign. The `character OUTPUT`
command reports hash-checked native orbital coefficients and occupations;
SCF character is not a multistate natural-orbital selection and cannot be
published through the reviewed-selection publisher. Assign indices, construct
explicit compact probes, then review their optimized orbitals and states.

When new installed code must be loaded as well as a new panel started:

```bash
python -m kinbot.orbital_review_handoff ROOT PLAN OUTPUT \
  --case-input ROOT/CASE/input.json --expected-max-jobs 23
python -m kinbot.orbital_review_handoff ROOT PLAN OUTPUT \
  --case-input ROOT/CASE/input.json --expected-max-jobs 23 \
  --apply --wait-idle 600 --wait-owner 120
```

The first command previews the panel. The second preflights it, uses the existing
controller refresh under its lock, waits for replacement owners and pending input
revisions, reconciles accepted models with live controls, stages/resumes the panel,
and runs it. Repeat `--case-input` for more cases. No scientific model or resource
policy is changed by this orchestration. The optional expected limit checks the
configured pool size, not scheduler availability. Review jobs use shared admission.
The existing refresh preserves running/pending QC and independent RotdPy owners;
repeating it with the same code is a no-op. The panel/handoff locks prevent duplicate
owners. Owner waits expire with a resumable error; scientific or ownership conflicts
are retained. New active-space publication is a separate explicit review action.

## Compact CASSCF probes from explicit SCF indices

A completed SCF character screen may seed a compact CASSCF comparison without
an intervening all-valence CI. This is an explicit review choice: the code does
not infer orbital character, active counts or state averages from occupations.
The candidate definition includes `"orbital_source": "scf"`, the chosen
`active_electrons`, `active_orbitals`, `state_model`, and `chemical_note`.
Do not combine this with an AVAS selector or a natural-orbital partition.

Each candidate anchor specifies the target `source` anchor, an optional different
`seed` anchor from the same source screen, and basis-specific maps:

```json
{
  "source": "near",
  "seed": "far",
  "active_indices": {"dz": [3, 5, 6], "tz": [3, 5, 7]},
  "seed_sha256": {"dz": "REVIEWED_DZ_WFU_SHA256", "tz": "REVIEWED_TZ_WFU_SHA256"}
}
```

The above indices and hashes are illustrative, not a molecular selection.
The source screen must be complete, and the declared WFU hashes must match its
validated records. Preparation also binds the source input/output hashes, checks
that every selected index is printed, requires the selected occupations to give
the declared active electron count, and includes every singly occupied orbital.
Each basis needs its own matching native seed. Passing these guards establishes
consistent initialization, not the correctness of the molecular state or CAS.

The native runner copies the original canonical record `2100.2` onto FILE,3,
starts the new RHF guess and MULTI explicitly from `2100.3`, and writes new SCF
and natural orbitals on a fresh FILE,2. Source-index exchanges occur once in
MULTI. This prevents the preliminary RHF calculation from replacing the reviewed
canonical indices before packing the active space. Native output, checkpoint,
XML exports and execution provenance remain available for review.

`candidates SOURCE RECIPE OUTPUT --preflight` checks without preparing output;
`--resume` verifies an identical existing comparison. A shared `batch` can mix
ordinary screen recipes and completed-source candidate recipes. Add
`"source": "RELATIVE/PATH/TO/COMPLETED/SCREEN"` to each candidate entry alongside
its `case`, `recipe` and `name`. Candidate sources must be inside the same
production policy tree. One campaign owner applies the plan cap and the existing
shared pool policy across all child manifests. No production model is published.

Preparation preflights the whole plan before creating the campaign. A later
failure removes only the newly created, unsubmitted staging; existing source
screens and production work remain intact. Atomic directory ownership prevents
a competing preparer from deleting an existing campaign. Repeated run commands
use the campaign/manifest locks and job journals. For a new panel, a production
controller restart is unnecessary: the new review process and workers load the
current installed code.

### Same-basis SCF seed restart and diagnosed rechecks

SCF-seeded compact probes now use normal `START,2100.3,TYPE=CANONICAL`
restarts in the matching basis. The source checkpoint remains on FILE,3 and the
new results use fresh FILE,2. The executor records this as
`native_xml_scf_seed_same_basis_v8`. The prepared input and its hash remain
unchanged; the transformed native input and its separate hash identify what ran.
Source indices, their packing exchanges, geometry, CAS, spin states, weights and
convergence requirements are retained.

The previous `native_xml_scf_seed_file_v7` executor forced `OVL`, even for the
same basis. Molpro's [orbital restart documentation](https://www.molpro.net/manual/doku.php?id=general_program_structure#selecting_orbitals_and_density_matrices_orbital_density)
describes `OVL` as overlap projection, normally used when basis dimensions
change. It can lose occupied rank across a large fragment displacement. A new
review owner automatically archives and retries a native VIRTORB orbital-count
failure once only when the old v7 transform, same-basis seed, job ID, source and
input hashes, and failure before MULTI are verified. Live jobs and uncertain
submissions are protected. A failure of the replacement setting is not retried
with the same setting. All original files and Slurm IDs remain in attempt history.
Neither `NOCHECK` nor `IGNORE_ERROR` is added.

A single-state calculation may print the bare heading `NATURAL ORBITALS`.
The completion reader accepts that table only for a declared single root,
retaining occupation bounds, electron-count, energy, convergence, checkpoint and
export checks. A state-specific table cannot replace a multi-state average.
An existing native success can therefore be re-audited without recomputation.
Computational completion is still separate from acceptance of orbital character.

If human inspection rejects the character of a *completed* legacy SCF-seeded
diagnostic, an explicit recheck can use the corrected restart:

```bash
python -m kinbot.requested_orbital_recon run PANEL --report-failures \
  --wait-owner 86400 --recheck-seed-job JOB_ID
```

Repeat `--recheck-seed-job` for additional inspected jobs. This option is limited
to intact, unpublished v7 completions with independent SCF seeds; it is not a
general permission to replace accepted results. It rejects orbital dependency
chains, live/ambiguous allocations and changed evidence. Every requested target
is checked before archival. The old ID is recorded as `reviewed_seed_recheck`;
repeating the handoff, including after partial application, does not repeat it.
Ordinary submission then enforces the same panel and shared-pool limits and
exclusions. The monitor displays inherited shared-pool exclusions separately.

`--wait-owner` waits for the current review owner to finish, without sending
signals or stopping allocations. It expires with an error after the specified
number of seconds. The replacement review process and subsequently launched
workers load the updated code; production controllers need no restart for this
review-only change. These orchestration and parser paths are locally tested
with simulated scheduling and retained native excerpts. Native recovery and the
intended orbital subspace must be verified from the new outputs.
