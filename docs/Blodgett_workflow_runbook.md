# Blodgett orbital reconnaissance checkpoint — 14 October 2026

Run all commands on Blodgett yourself. No SSH or cluster job has been started
by this development session. Use the existing `composite` checkout and its
`.venv`; no SCP or replacement source directory is required. The orbital
campaign can run while the remaining MESS audit continues. Do not start the
full KinBot chemistry campaign until the reconnaissance has been reviewed.

## Pull and install

The checkout confirmed by the user's Blodgett output is `/home/hram/KinBot`.
Run each block in a child Bash shell if using `set -e`; an optional command
must not terminate the SSH login shell. Commit `4793a0f` contains the tested
reconnaissance implementation; checking only `4bd1d25` is insufficient.

```bash
cd /home/hram/KinBot
KB_REPO=/home/hram/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite

RECON_COMMIT=4793a0f
git merge-base --is-ancestor "$RECON_COMMIT" HEAD || {
    echo 'Required reconnaissance checkpoint is missing'
    exit 1
}

.venv/bin/python -m pip install -e . --no-build-isolation --no-deps
KB_PY="$KB_REPO/.venv/bin/python"
"$KB_PY" -c 'import kinbot, kinbot.orbital_campaign; print(kinbot.__file__)'
```

Preparation discovers Molpro from the current environment, configured
installation paths, or available site modules. It writes a versioned module
load and executable PATH into each `site_setup.sh` when available. Every batch
job explicitly selects `KINBOT_BACKEND=molpro`, initializes modules, and sources
that setup. Modules that supply dependencies are loaded by the site's Molpro
module. Single-node Intel MPI defaults to shared-memory transport unless you
have explicitly selected another fabric. Each probe has isolated scratch.
No Gaussian, CFOUR, MESS, ROTD_py, or UMA calculation is needed for this phase.

If automatic discovery cannot expose Molpro, the preflight below fails before
submission. Load the site's actual Molpro module, then prepare a **fresh**
campaign so that it records that environment. If the installation needs an
explicit executable override, set `KINBOT_MOLPRO_EXECUTABLE` to its absolute
launcher path before preparation. Do not invent a module version name or add
an unrelated Conda library directory to `LD_LIBRARY_PATH`.

## Preserve the acetaldehyde failure evidence and clean its active path

The cluster restart explains the dead driver/jobs and stale `running` labels.
The archived execution record identifies Molpro exit 255 at the first force
evaluation; its native stdout/stderr/output are needed for the cause. HIR
harmonic fallbacks and the optional vdW-well failure were nonfatal decisions.

```bash
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
RUN=/home/hram/kinbot-runs/acetaldehyde_anl0_f12_day
POINT="$RUN/vrctst/molpro_scan/441041030570000000001_hom_sci_1_2_vts_pt00_5914412a7cd7e609"
if [ -d "$POINT/attempts/relaxed_scan/1" ]; then
    cat "$POINT/attempts/relaxed_scan/1/execution.json"
    find "$POINT/attempts/relaxed_scan/1" -maxdepth 2 -type f -print
fi

cd /home/hram/kinbot-runs
"$KB_PY" -m kinbot.run_cleanup "$RUN"
"$KB_PY" -m kinbot.run_cleanup "$RUN" --stop --archive
```

The cleanup targets current owned Linux processes in that run tree and Slurm
jobs whose current WorkDir lies inside it. It checks PID identity, protects
the invoking shell, cancels newly observed scoped jobs, verifies that work
has stopped, and moves the entire tree to a dated sibling. The printed
`archive` contains all attempts and `cleanup_inventory.json`. It refuses to
archive during a scheduler-query failure or while work remains alive. It
does not delete shared scratch or terminate unrelated jobs. Keep this archive
for diagnosis of the earlier transient failure.

Blodgett reports `Slurm accounting storage is disabled`. No `sacct` command
is required. Cleanup and monitoring query the live scheduler with `squeue`.
Do not place optional accounting queries unguarded under `set -e`.

## Prepare and launch both systems with ten nodes total

This block initializes its own paths. Preparation is local input generation,
not submission. Use a fresh path.

```bash
cd /home/hram/KinBot
KB_REPO=/home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
CAMPAIGN="/home/hram/kinbot-runs/ch3s_ch3ss_orbital_$(date +%Y%m%d_%H%M%S)"
"$KB_PY" -m kinbot.orbital_campaign prepare "$CAMPAIGN" \
    "$KB_REPO/examples/anl/ch3s_anl1_f12/orbital_reconnaissance.json" \
    "$KB_REPO/examples/anl/ch3ss_anl1_f12/orbital_reconnaissance.json"

"$KB_PY" -m kinbot.orbital_campaign preflight "$CAMPAIGN" || {
    echo 'Site setup failed; do not launch this campaign'
    exit 1
}
printf '%s\n' "$CAMPAIGN" > "$HOME/kinbot-runs/last_sulfur_orbital_campaign.txt"

nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$CAMPAIGN" \
    --max-jobs 10 --max-refinements 3 \
    > "$CAMPAIGN/controller.log" 2>&1 < /dev/null &
printf '%s\n' "$!" > "$CAMPAIGN/controller.pid"
"$KB_PY" -m kinbot.orbital_campaign progress "$CAMPAIGN"
tail -n 40 "$CAMPAIGN/controller.log"
```

Each fresh manifest first submits one probe and waits for its validated
completion before widening submission. The controller then completes all
native pilots: 8 CH3S jobs and 16 CH3SS jobs,
covering all eight bases for each dissociation. Any invalid native basis,
AVAS/CASCI block, root output, or missing orbital artifact stops progression
before bulk scans. Successful pilots are followed by CH3S then CH3SS scans;
the controller waits for allocations to leave before switching systems.
**Do not launch a second controller:** ten is a shared campaign limit, not
an allowance for ten nodes per concurrent system.

Each job requests one node, eight MPI ranks, one thread per rank, 1000 MW
Molpro stack per rank, 80000 MB Slurm memory including headroom, and 24 hours
on `day-long-cpu`. These are explicit example resource settings. If that
partition has changed, edit the example resource fields **before preparing a
fresh campaign**. A rejected pilot submission stops before bulk work.

There are initially 248 full CH3S probes and 600 full CH3SS probes, plus the
24 pilots: **872 native jobs**, never all submitted simultaneously. Dense
radial grids, rigid-fragment orientation probes and a 30 Å initial endpoint
replace the old three-distance approach. All requested root energies and
natural occupations drive interval/range refinement, up to three passes.
An unresolved plateau is reported explicitly; it is not silently accepted.
The eight native bases are augmented conventional D/T/Q/5Z and augmented
F12 D/T/Q/5Z, with Molpro's default definitions and no rev2 override. Actual
availability is checked by the installed Molpro in the pilots.

## Monitor, resume, and collect

In a later shell, start in the same checkout and restore these variables:

```bash
KB_REPO=$(git rev-parse --show-toplevel)
KB_PY=/home/hram/KinBot/.venv/bin/python
CAMPAIGN=$(cat "$HOME/kinbot-runs/last_sulfur_orbital_campaign.txt")
"$KB_PY" -m kinbot.orbital_campaign progress "$CAMPAIGN"
tail -n 60 "$CAMPAIGN/controller.log"
squeue -u "$USER" -o '%.18i %.12P %.35j %.8T %.10M %.6D %R'
```

A disconnected interactive shell does not stop the `nohup` controller.
After a controller interruption, the **same** `run` command resumes completed
and already submitted work without duplicating it; a manifest lock rejects
another live controller. An unsuccessful native job needs diagnosis before
retry; do not erase its journal or blindly resubmit it. Preserve stdout,
stderr, Molpro output and Slurm accounting. The controller stops rather than
repeating a known failed job indefinitely.

Once the controller has finished (including an explicit unresolved-coverage
status), collect review evidence:

```bash
"$KB_PY" -m kinbot.orbital_campaign collect "$CAMPAIGN" "$CAMPAIGN.review.tar.gz"
```

The archive includes immutable input/context snapshots, manifests, input
decks, native outputs, stdout/stderr and Molden files. Binary WFU records stay
on Blodgett and must be retained: the eventual reviewed active-space models
reference them by path and SHA-256. Share the review archive/results before
launching the full KinBot campaigns. No automatic minimal-CAS or reactive
state selection is inferred merely from a population cutoff.

## Procedure after reconnaissance

Review occupation and orbital character across every basis, radial interval,
orientation and channel; resolve coverage failures and state crossings;
compare larger diagnostic spaces where needed. Select consistent active
spaces, averaged/physical roots and basis-specific WFU records. Establish
spin–orbit/vibronic energy references and count physical states once. Then
validate native correction, ROTD sampling, Ne staging and MESS per channel
before full production. Homolytic vdW-complex discovery is deferred under the
current scope; existing inherited optional-well behavior is unchanged.

The ongoing MESS work and actual NH2O SI verification are recorded in
[the audit](ANL_workflow_audit.md) and
[the SI-derived contract report](NH2O_SI_verification.json).

## First native pilot failure and correction

The 19:26 campaign ended at Molpro input parsing: its generated ROHF line
had an extra closing brace (`UNMATCHED }`). This was a generator error,
not evidence that the requested basis was unavailable. The correction removes
the extra brace and structurally validates every generated deck before it
can be submitted. The native output also confirmed the default WFU repository
was `/home/hram/wfu`; execution now passes `-W` with each job's own directory
so checkpoint staging and artifact auditing use the same location.

The controller now writes a native `execution.json`, includes the final
Molpro output in failure diagnostics, and checks a single first probe before
widening a fresh manifest. Preserve the failed campaign using
`kinbot.run_cleanup --stop --archive`, and prepare a fresh campaign from the
corrected code. Its changed input hashes must not be written over the old
submission history. The initial 872-job count is unchanged.
