# Chemical orbital review and production handoff

The general entry point is `python -m kinbot.orbital_review`. It uses the
existing Molpro/Slurm execution machinery. There are no molecule-name branches
or prescribed sulfur active spaces in the implementation.

The target is the **requested eventual ROTD calculation**: chemically adequate
states, accurate interfragment potentials and spin gaps, continuity throughout
the contributing geometries, and affordable sampling. These requirements are
joint. A smooth small-space curve can be inaccurate; a larger space or a lower
total energy is not automatically a better reference. Compare chemically
motivated alternatives and justified correlated references before accepting
cost savings. Orbital occupations and CASCI projections guide hypotheses;
they do not certify an error bound or replace reaction-specific chemistry.

The [October 15 native review](SULFUR_ORBITAL_REVIEW_2026-10-15.md) identifies
solver failures hidden by the original completion checks, gives candidate
spaces from the actual orbital data, and now combines the successful repairs
with the original panel. All eight geometries in each sulfur channel have
converged fresh and checkpoint-initialized references. Use that combined
evidence for candidate selection; no further reference-panel restart is needed.
The [sulfur correlated validation](SULFUR_CORRELATED_VALIDATION.md) now gives
the finite next batch, 20-job launch/monitor/export commands and the remaining
checks required to complete the examples.

See also the [generality and automatic handoff audit](GENERALITY_AUDIT_2026-10-15.md).

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

Inner and outer models can require different spaces and averages. For example,
the [NH2O study, section 2.2](https://doi.org/10.1021/acs.jpca.3c05181)
uses (3e,3o) for inner NH + OH addition and (5e,4o) outside, where the additional
orbital represents the near-degenerate OH states. These are literature examples
of chemical reasoning, not defaults copied into other reactions. Stationary
well/reaction energies remain the responsibility of the CBH/ANL energy framework;
the VRC calculation must also deliver accurate relative interaction energies.

## Procedure

1. `plan CONFIG PLAN.json` derives fragment/spin context from the existing
   reconnaissance configuration and linked KinBot input. It exposes all
   radial/orientation anchors and selects a small finite comparison panel.
   The initial model is a diagnostic valence reference, not a production choice.
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
   Successful remaining jobs are retained; failed energies remain rejected.
   Unsubmitted jobs and descendants of failed checkpoints are not launched.
   The final campaign status is `failed`, requiring review before a repair.
   Uncertain submissions and scheduler errors still stop immediately, since
   allocation ownership cannot safely be inferred in those cases.
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
   Optional `casscf_options` overrides the partition preparation's solver defaults;
   optional `correlation_options` sets `max_iterations` and/or
   `max_reference_iterations` for the correlated calls. Both limits must be
   positive integers. The latter is Molpro's internal/reference-CI `MAXITI`,
   distinct from CASSCF macroiterations and from the correlated `MAXIT` limit.
   These controls change iteration budgets without changing convergence
   tolerances, the Hamiltonian, or the selected states. Omitted correlation
   controls retain native defaults. See the
   [Molpro MRCI iteration controls](https://www.molpro.net/manual/doku.php?id=the_mrci_program#maximum_number_of_iterations),
   which also apply to [RS2](https://www.molpro.net/manual/doku.php?id=multireference_rayleigh_schroedinger_perturbation_theory).
   Source CASSCF convergence, hashes, geometry and basis are checked. Partitions
   reoptimize all closed and active orbitals and retain explicit source provenance.
   They cannot change basis. Optional `target_distance_angstrom` generates a new
   point on the source's declared dissociation coordinate, retaining its fragment
   orientation and reoptimizing from the reviewed seed. Source/target geometry
   and the checkpoint hash are retained separately. This permits small stencils
   in the actual requested scan region without repeating a large-space recon.
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
interruption still needs Slurm accounting to establish its cause; a stale
`status="running"` execution record is not proof of native nonconvergence.
Campaign `collect` includes `scheduler_accounting.json`, queried with `sacct`
for only the recorded allocation IDs and their prior attempts. This preserves
states, exit codes/signals and restart counts alongside native diagnostics.
An unavailable accounting service is recorded without blocking the export.
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
`vrc_tst_active_space` and `vrc_tst_state_model_file`. Thus different channels
and inner/outer regions can use different CAS counts, orbital records and
state averages. A state-model file may contain several basis-specific WFU
records for the same space; use a different state-model file for a different
space or state average.

Install each user-selected model into a **new** KinBot input using:

```text
python -m kinbot.orbital_review install-selection INPUT.json SELECTED.json NEXT.json \
    --surface-file SURFACES.json --bond 1 2 --stage inner
```

The command writes `NEXT.json` and `NEXT.surfaces.json`, preserves other
bottlenecks, checks that the selected average can supply that stage's surfaces,
and prevents inherited global root overrides from changing the selection.
Use those new files as the inputs for the next stage selection. It leaves the
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

## Blodgett: stop, update, reuse and restart

Run this block in the existing Blodgett SSH session. `set -e` is confined to a
child shell, so a failed command does not exit the SSH login shell. It stops
only the specified campaign trees. It does not remove old results or WFU files.
Pull first to install the corrected cleanup command. Cancellation waits up to
five minutes by default: Slurm can retain a job in `COMPLETING` while its
processes and node epilog finish. The command prints progress and saves
`cleanup_status.json`; a timeout reports the remaining jobs/processes, returns
nonzero and leaves the directory in place. A repeated `--stop` is safe once
the queue has drained. Use a longer `--timeout SECONDS` if needed; do not
restart while scoped jobs remain.
The fixed new campaign path makes a repeated launch resume the same work;
use a new path if you edit the comparison plan.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
OLD=/home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525
NEW=/home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_20261015
PLANS=/home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_20261015_inputs

env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite
"$KB_PY" -m kinbot.run_cleanup "$OLD" --stop --timeout 300
if [ -f "$NEW/campaign.json" ]; then
    "$KB_PY" -m kinbot.run_cleanup "$NEW" --stop --timeout 300
fi
"$KB_PY" -m pip install -e . --no-build-isolation --no-deps
"$KB_PY" -m kinbot.orbital_review --help
"$KB_PY" -m kinbot.orbital_campaign archive-excluded "$OLD"
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
export KINBOT_MOLPRO_SCRATCH_MIN_MB=65536
mkdir -p /home/hram/.cache/kinbot/molpro

mkdir -p "$PLANS"
if [ ! -f "$NEW/campaign.json" ]; then
    STAMP=$(date +%Y%m%d_%H%M%S)
    "$KB_PY" -m kinbot.orbital_review plan \
        /home/hram/KinBot/examples/anl/ch3s_anl1_f12/orbital_reconnaissance.json \
        "$PLANS/ch3s_$STAMP.json"
    "$KB_PY" -m kinbot.orbital_review plan \
        /home/hram/KinBot/examples/anl/ch3ss_anl1_f12/orbital_reconnaissance.json \
        "$PLANS/ch3ss_$STAMP.json"
    "$KB_PY" -m kinbot.orbital_review prepare "$PLANS/ch3s_$STAMP.json" "$NEW" \
        --additional-plan "$PLANS/ch3ss_$STAMP.json" --reuse-from "$OLD"
else
    "$KB_PY" -m kinbot.orbital_campaign retry-failed "$NEW" --include-interrupted
fi
"$KB_PY" -m kinbot.orbital_campaign preflight "$NEW"
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$NEW" \
    --max-jobs 20 --max-refinements 0 \
    >> "$NEW/controller.log" 2>&1 < /dev/null &
printf '%s\n' "$!" > "$NEW/controller.pid"
"$KB_PY" -m kinbot.orbital_campaign progress "$NEW"
tail -n 30 "$NEW/controller.log"
BASH
```

This first run resolves initialization and optimized reference-orbital
behavior. It does not assume CAS(5,4), CAS(9,6), a particular inner/outer
boundary or a production average. We can use its export with the existing
CASCI archive to choose the smaller candidate comparisons hands-on. Completed
new jobs are reused on resume; an interrupted allocation is explicitly reset.
A native failure still needs its reported diagnostic addressed before retry.
With the current examples, preparation produces 72 finite-panel jobs and
nine native pilots, under one twenty-job cap. Every requested basis is retained;
the examples currently request only TZ-F12. No correlated checks are enabled
in this initial reference-orbital campaign.

Read status at any time:

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
/home/hram/KinBot/.venv/bin/python -m kinbot.orbital_campaign progress \
    /home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_20261015
tail -n 50 /home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_20261015/controller.log
squeue -u hram -o '%.18i %.16P %.48j %.10T %.10M %.6D %R'
BASH
```

After the controller finishes, export:

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
NEW=/home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_20261015
STAMP=$(date +%Y%m%d_%H%M%S)
"$KB_PY" -m kinbot.orbital_review report "$NEW" "$NEW/review_report_$STAMP.json"
"$KB_PY" -m kinbot.orbital_campaign collect "$NEW" \
    "/home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_$STAMP.tar.gz"
BASH
```

Copy the printed archive to the local checkout for review. Keep both campaign
directories on Blodgett because the selection and restart provenance refer
to their WFU paths. Production launch commands follow the user's actual
space/surface selections; the current example inputs deliberately reference
state/surface model files that have not yet been scientifically selected.
