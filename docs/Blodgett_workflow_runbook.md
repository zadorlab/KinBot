# Blodgett orbital reconnaissance checkpoint — 14 October 2026

**Current procedure:** see [implementation and acceptance status](WORKFLOW_READINESS_2026-10-14.md).
Preserve the existing campaign. Its reduced D/T/Q basis selection is retained.
For the latest CH3SS scratch I/O failure, use the final recovery section below.
Conventional-5Z RHF completed
using about 138 GB of scratch, then Molden rejected h functions. Use the
published native XML/WFU recovery below after all campaign allocations finish.

Run all commands on Blodgett yourself. No SSH or cluster job has been started
by this development session. Use the existing `composite` checkout and its
`.venv`; no SCP or replacement source directory is required. The orbital
campaign can run independently of MESS development. Do not start the
full KinBot chemistry campaign until the reconnaissance has been reviewed.

## Pull and install

The checkout confirmed by the user's Blodgett output is `/home/hram/KinBot`.
Run each block in a child Bash shell if using `set -e`; an optional command
must not terminate the SSH login shell. Commit `c14913f` contains the tested
native export recovery; checking only `4bd1d25` is insufficient.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_REPO=/home/hram/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite

RECON_COMMIT=c14913f
git merge-base --is-ancestor "$RECON_COMMIT" HEAD || {
    echo 'Required reconnaissance checkpoint is missing'
    exit 1
}

.venv/bin/python -m pip install -e . --no-build-isolation --no-deps
KB_PY="$KB_REPO/.venv/bin/python"
"$KB_PY" -c 'import kinbot, kinbot.orbital_campaign; print(kinbot.__file__)'
BASH
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

For the current campaign use the final restart section instead. This fresh
campaign block applies the recommended six-basis selection explicitly.
This block initializes its own paths. Preparation is local input generation,
not submission. Use a fresh path.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_REPO=/home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
CAMPAIGN="/home/hram/kinbot-runs/ch3s_ch3ss_orbital_$(date +%Y%m%d_%H%M%S)"
export KINBOT_MOLPRO_SCRATCH="/scratch/$USER/kinbot-molpro"
"$KB_PY" -m kinbot.orbital_campaign prepare "$CAMPAIGN" \
    "$KB_REPO/examples/anl/ch3s_anl1_f12/orbital_reconnaissance.json" \
    "$KB_REPO/examples/anl/ch3ss_anl1_f12/orbital_reconnaissance.json"

"$KB_PY" -m kinbot.orbital_campaign exclude-bases "$CAMPAIGN" \
    --basis aug-cc-pV5Z --basis av5z-f12

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
BASH
```

Each fresh manifest first submits one probe and waits for its validated
completion before widening submission. The controller then completes all
native pilots: 6 CH3S jobs and 12 CH3SS jobs,
covering the six selected bases for each dissociation. Any invalid native basis,
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

There are initially 186 full CH3S probes and 450 full CH3SS probes, plus the
18 pilots: **654 native jobs**, never all submitted simultaneously. Dense
radial grids, rigid-fragment orientation probes and a 30 Å initial endpoint
replace the old three-distance approach. All requested root energies and
natural occupations drive interval/range refinement, up to three passes.
An unresolved plateau is reported explicitly; it is not silently accepted.
The selected bases are augmented conventional D/T/QZ and augmented
F12 D/T/QZ, with Molpro's default definitions and no rev2 override. The
original example configurations retain eight bases; omitting `exclude-bases`
in a fresh campaign explicitly requests the 872-job comparison. Actual
availability is checked by the installed Molpro in the pilots.

## Monitor, resume, and collect

In a later shell, start in the same checkout and restore these variables:

```bash
cd /home/hram/KinBot
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
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
CAMPAIGN=$(cat /home/hram/kinbot-runs/last_sulfur_orbital_campaign.txt)
"$KB_PY" -m kinbot.orbital_campaign collect "$CAMPAIGN" "$CAMPAIGN.review.tar.gz"
```

The archive includes immutable input/context snapshots, manifests, input
decks, native outputs, stdout/stderr, XML/sidecars and legacy Molden files. Binary WFU records stay
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

## Successful native pilot with a parser rejection

The replacement CH3S augmented-DZ pilot completed normally, returned zero,
saved both Molden files and its job-local WFU, and printed all three doublet
and three quartet energies. Its variable assignments use the unit label
`HARTREE`, while the initial reconnaissance reader required `AU`. The reader
now accepts both equivalent labels. The same unit handling is corrected in
the ASE gradient/sampling parser and the ANL stored-F12-energy validation;
repeated assignment/SHOW values must agree with the native total energy.

This is a reader-only correction: **preserve and resume this campaign in
place** after pulling it. Do not archive or regenerate the successful pilot,
and do not resubmit its job manually. The controller re-audits its existing
output, accepts the completed probe, and continues the remaining pilots.

## Recover the existing scratch-failed campaign

This older recovery sequence retains all bases. Use the final section for the
latest requested cost-reduced restart.
The compute-node inventory supplied by the user shows `/scratch` has about
5.4 TB free; `/tmp` has about 37 GB. The login node has no `/scratch` mount.
Set `KINBOT_MOLPRO_SCRATCH="/scratch/$USER/kinbot-molpro"` before starting the
controller. The worker creates/checks it on its allocated node, and submission
explicitly exports that environment. No login-node scratch mkdir is needed.

After pulling the final checkpoint and reinstalling in `/home/hram/KinBot`,
with all campaign allocations gone, run:

```bash
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
CAMPAIGN=$(cat /home/hram/kinbot-runs/last_sulfur_orbital_campaign.txt)
export KINBOT_MOLPRO_SCRATCH="/scratch/$USER/kinbot-molpro"
"$KB_PY" -m kinbot.orbital_campaign retry-failed "$CAMPAIGN"
"$KB_PY" -m kinbot.orbital_campaign preflight "$CAMPAIGN"
cd "$CAMPAIGN"
if [ -f controller.log ]; then
    mv controller.log "controller.previous.$(date +%Y%m%d_%H%M%S).log"
fi
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$CAMPAIGN" \
    --max-jobs 10 --max-refinements 3 \
    > controller.log 2>&1 < /dev/null &
printf '%s\n' "$!" > controller.pid
```

`retry-failed` re-audits completed outputs, preserves successful probes,
archives only known failed native attempts, and clears only their old job IDs.
It refuses live allocations or uncertain submissions. Do not resubmit `.sbatch`
files manually or delete journals. A native failure after retry needs its own
inspection; there is no unlimited automatic retry loop.

## Molpro high-angular-momentum orbital exports

The conventional 5Z RHF completed using about 138 GB of scratch, then its
Molden export failed because Molpro's interface does not support h functions.
This is not a 5Z-energy/VRC limitation. Execution now replaces only visualization
exports with Molpro XML in a separately hashed native deck, preserving the
original prepared input and all completed results. Native natural occupations,
XML and WFU evidence are audited. An export-only failed attempt can reuse its
converged RHF WFU if the native output, exact input hash and archived checkpoint
establish that it completed; other failed SCF jobs cannot use this recovery.
Wait until all pilot allocations leave before `retry-failed`. Keep using the
existing campaign and the compute-node scratch override. No basis downgrade
or removal of high-angular-momentum functions is performed.

The default VRC production sampling level is now single-basis augmented TZ-F12.
Correction-potential reference and spin-gap energies retain TZ/QZ-F12 CBS.
Explicit legacy sampling-CBS requests still work but are no longer the preset.
The existing sulfur reconnaissance retains all eight requested basis families,
including conventional 5Z and 5Z-F12. Changing the production sampling preset
does not change a campaign's snapshotted reconnaissance inputs.

Recovery release validation: 1236 tests and 391 subtests pass. This includes
native XML execution/auditing, old Molden acceptance, archived SCF recovery,
and scalar/CBS sampling compatibility. Licensed Molpro XML execution still
requires the next native retry; local tests verify the execution contract.

## Recommended cost-reduced restart of the existing campaign

The last reported live allocation was `52907305`. Stop it with the following
command if it is still running, then wait for it to leave `squeue`. A job
already gone need not be cancelled. Do not delete or archive the campaign.

```bash
cd /home/hram/KinBot
scancel 52907305
```

After pulling the release containing `exclude-bases`, this block preserves the
six validated CH3S pilots and the original 5Z output/WFU evidence. It selects
six bases in all four manifests and preserves that choice through adaptive
refinement. The remaining initial counts are 186/450 full probes and 6/12
pilots, or 654 jobs total before refinements. The command refuses live
campaign allocations, uncertain submissions, unknown bases or an empty plan.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot

env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite

git merge-base --is-ancestor c14913f HEAD
.venv/bin/python -m pip install -e . --no-build-isolation --no-deps

KB_PY=/home/hram/KinBot/.venv/bin/python
CAMPAIGN=$(cat /home/hram/kinbot-runs/last_sulfur_orbital_campaign.txt)
export KINBOT_MOLPRO_SCRATCH="/scratch/$USER/kinbot-molpro"

"$KB_PY" -m kinbot.orbital_campaign exclude-bases "$CAMPAIGN" \
    --basis aug-cc-pV5Z --basis av5z-f12
"$KB_PY" -m kinbot.orbital_campaign retry-failed "$CAMPAIGN"
"$KB_PY" -m kinbot.orbital_campaign preflight "$CAMPAIGN"

cd "$CAMPAIGN"
if [ -f controller.log ]; then
    mv controller.log "controller.previous.$(date +%Y%m%d_%H%M%S).log"
fi
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$CAMPAIGN" \
    --max-jobs 10 --max-refinements 3 \
    > controller.log 2>&1 < /dev/null &
printf '%s\n' "$!" > controller.pid
printf 'Campaign: %s\nController PID: %s\n' "$CAMPAIGN" "$!"
BASH
```

Each allocated job loads the preserved Molpro site setup automatically.
There is no new guessed module name or basis data to install. An optional
positive integer `KINBOT_MOLPRO_SCRATCH_MIN_MB` raises the compute-node free
capacity floor if desired; it cannot weaken a task's requested minimum.

Final resumed checkpoint validation: **1251 tests and 391 subtests pass**;
all **159 non-template Python modules** parse and `git diff --check` passes.
The actual NH2O SI verification also passes again. Tests cover direct/PES
scalar and multistate assembly with UQ energy/mode changes, retained optional
exit selection, reference-root contracts, native export rendering, runtime
failure records, cost-reduced restart and preservation through refinement.
Native Molpro/Slurm/MESS execution remains a cluster acceptance step.

## Recover the CH3SS XML title overflow in the same campaign

The reported `strings::string_insert` failure happened at the first XML export,
after RHF and before AVAS/CASCI. The generated title exceeded the exporter's
internal string capacity. Native recon and production correction decks now use
short deterministic ASCII titles. This is a general execution safeguard, not
a molecule-specific or basis-specific exception. Prepared recon inputs and
their hashes remain unchanged. Completed legacy Molden results remain valid.

After pulling this title-fix release, use the following block once all campaign
allocations have left. The scheduler-checked retry preserves six completed
CH3S pilots, the existing D/T/Q basis selection, and all failure evidence. It
reuses the failed CH3SS pilot's RHF only if exact input, a finite native RHF
energy, the RHF program summary and its archived WFU support recovery.
Otherwise the failed probe starts again from its normal initialization.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite
.venv/bin/python -m pip install -e . --no-build-isolation --no-deps
KB_PY=/home/hram/KinBot/.venv/bin/python
"$KB_PY" -c 'from kinbot.molpro import compact_molpro_titles; assert len(compact_molpro_titles("***," + "x" * 100).split(",", 1)[1]) <= 32'
CAMPAIGN=$(cat /home/hram/kinbot-runs/last_sulfur_orbital_campaign.txt)
export KINBOT_MOLPRO_SCRATCH="/scratch/$USER/kinbot-molpro"
"$KB_PY" -m kinbot.orbital_campaign retry-failed "$CAMPAIGN"
"$KB_PY" -m kinbot.orbital_campaign preflight "$CAMPAIGN"
cd "$CAMPAIGN"
if [ -f controller.log ]; then
    mv controller.log "controller.previous.$(date +%Y%m%d_%H%M%S).log"
fi
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$CAMPAIGN" \
    --max-jobs 10 --max-refinements 3 \
    > controller.log 2>&1 < /dev/null &
printf '%s\n' "$!" > controller.pid
printf 'Campaign: %s\nController PID: %s\n' "$CAMPAIGN" "$!"
BASH
```

Preserved `site_setup.sh` files still load Molpro for each allocation.
No new campaign, module guess, basis exclusion, or manual `.sbatch` submission
is needed. The native retry is the acceptance check for this Molpro export
path; local regression fixtures do not replace licensed execution.

Title-fix release validation: **1260 tests and 391 subtests pass**, including
reported-title execution, immutable-input preservation, SCF recovery from
both known export failures, rejection of unproven SCF recovery, long production
channel names, Unicode titles and padded user titles. All 159 non-template
Python modules parse; all nine runbook Bash blocks pass `bash -n`.

## Recover the CH3SS QZ launcher scratch I/O failure

Five CH3SS pilots completed after the title fix. The conventional-QZ C–S
pilot then failed before creating native output because its Molpro launcher
received `Input/output error` from `mktemp` under `/scratch`. This is distinct
from the earlier export failure. The log establishes a failed filesystem
operation; it does not identify the underlying disk, mount or node fault.
Free space from `df` does not establish that writes or metadata operations
work. Code cannot repair that filesystem.

The runtime now tests nested directories and file I/O and records rejected
roots. An explicit override is strict unless `KINBOT_MOLPRO_SCRATCH_FALLBACKS`
permits alternatives, in order, separated by the platform's path separator
(`:` on Blodgett). The same capacity requirement applies to every candidate.
Reconnaissance records its host and launcher stdout/stderr, including stderr
in controller errors even when Molpro produced no `.out` file. The shared
runtime check also covers production Molpro paths. There is no automatic
native retry or guarantee against a filesystem failing after its probe.

Let the remaining job `52907312` finish. Wait for all campaign allocations to
leave before running the following block. It preserves successful pilots and
the reduced basis selection. New jobs try `/scratch2` and may explicitly use
shared home scratch if that local filesystem is unavailable or fails checks.
The 65536 MB floor is an operational free-space threshold, not a forecast or
reservation of a particular calculation's peak storage.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite
.venv/bin/python -m pip install -e . --no-build-isolation --no-deps
KB_PY=/home/hram/KinBot/.venv/bin/python
"$KB_PY" -c 'from kinbot.anl.runtime import _check_scratch_io'
CAMPAIGN=$(cat /home/hram/kinbot-runs/last_sulfur_orbital_campaign.txt)
export KINBOT_MOLPRO_SCRATCH="/scratch2/$USER/kinbot-molpro"
export KINBOT_MOLPRO_SCRATCH_FALLBACKS="$HOME/.cache/kinbot/molpro"
export KINBOT_MOLPRO_SCRATCH_MIN_MB=65536
"$KB_PY" -m kinbot.orbital_campaign retry-failed "$CAMPAIGN"
"$KB_PY" -m kinbot.orbital_campaign preflight "$CAMPAIGN"
cd "$CAMPAIGN"
if [ -f controller.log ]; then
    mv controller.log "controller.previous.$(date +%Y%m%d_%H%M%S).log"
fi
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$CAMPAIGN" \
    --max-jobs 10 --max-refinements 3 \
    > controller.log 2>&1 < /dev/null &
printf '%s\n' "$!" > controller.pid
printf 'Campaign: %s\nController PID: %s\n' "$CAMPAIGN" "$!"
BASH
```

Each compute job evaluates scratch health on its allocated node. The login
preflight checks modules and scheduler access; it does not validate compute
filesystems. Existing Molpro site setup remains in effect. Leave archived
attempts and completed WFU records intact.

The follow-up execution metadata records 5,564,398 MB free for failed job
`52907311`, corroborating that a simple capacity shortfall does not explain
its `mktemp` EIO. A second failed job, `52907317` (QZ-F12 C–S), has an empty
runtime record and requires its own setup error to be examined. Slurm no
longer retains either job in `scontrol show job`; identify their nodes in
controller/job-completion logs before investigating mount/kernel/storage
errors. Do not infer the failed node from another job's live allocation.

Scratch-fix release validation: **1267 tests and 391 subtests pass**; all 159
non-template Python modules parse and all ten runbook Bash blocks pass syntax
checks. Regression coverage includes failed directory allocation, fsync and
read-back, explicit override preservation, capacity floors and launcher
stderr/host retention when native output is absent.
