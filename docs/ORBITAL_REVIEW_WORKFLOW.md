# Chemical orbital review and production handoff

The general entry point is `python -m kinbot.orbital_review`. It uses the
existing Molpro/Slurm execution machinery. There are no molecule-name branches
or prescribed sulfur active spaces in the implementation.

The [October 15 native review](SULFUR_ORBITAL_REVIEW_2026-10-15.md) identifies
solver failures hidden by the original completion checks, gives candidate
spaces from the actual orbital data, and provides a selective repair command.
Use that follow-up for the completed sulfur campaign rather than rerunning
the full reference panel.

## Procedure

1. `plan CONFIG PLAN.json` derives fragment/spin context from the existing
   reconnaissance configuration and linked KinBot input. It exposes all
   radial/orientation anchors and selects a small finite comparison panel.
   The initial model is a diagnostic valence reference, not a production choice.
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
6. `report ROOT OUTPUT.json` supplies geometry, chemical rationale, requested
   state averages, occupations, complete native orbital/CI text, root energies,
   elapsed times and hash-bound artifact paths. It compares geometrically
   equivalent calculations across initializations and reports radial energies
   relative to the farthest completed radial point. That point is labeled as a
   finite reference, not automatically certified as an asymptote. Native XML
   exports and the full outputs are included by `orbital_campaign collect`;
   binary WFU files remain on the cluster.
7. After orbital/initialization review, duplicate a plan into a fresh small
   correlated validation campaign. Each model may request `correlation_checks`,
   for example `[{"method":"caspt2-f12b(5,4)","shift":0.3,"ipea":0.25},
   {"method":"mrci-f12b+q(5,4)"}]`. Counts must match that model. Checks run
   only with optimized CASSCF, for every requested spin/root, and use the
   production F12B and Davidson-energy selectors. They are disabled by default
   on the larger reference spaces. Inspect interaction energies, spin gaps,
   intruder/shift sensitivity and rate sensitivity before a final selection.
8. The user chooses the final space, orbital ordering and states. Register the
   actual reviewed WFU with `molpro_reconnaissance select --basis ...`; supply
   the reviewed counts, average and root choices. The command does not reorder
   orbitals or certify chemical accuracy. New selections carry
   `state_average_policy="reviewed"`, so production preserves their averages
   and rejects conflicting count overrides. Old files without this field retain
   the legacy channel-derived policy; recreate them to preserve a reviewed
   average explicitly.

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
