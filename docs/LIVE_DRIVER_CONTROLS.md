# Live driver controls

These opt-in Slurm controls change scheduling and request another calculation
attempt without changing chemistry or editing source code. A controller must
already include this implementation. An older running Python process needs one
controller refresh; subsequent control commands do not require refreshes.
Running allocations keep their immutable inputs and code snapshot.

## Enable and inspect

Use the common parent directory for cooperating PES trees:

```bash
KB_RUN=/absolute/path/to/run
python -m kinbot.run_control enable "$KB_RUN"
python -m kinbot.run_control list "$KB_RUN"
python -m kinbot.run_control show "$KB_RUN"
```

`list` reports dispatcher targets and registered ordinary QC/search targets.
`show` reports the configuration and each operator request's receipt. The regular
`kinbot.run_monitor` shows a compact control summary. Reading these reports does
not apply requests, submit calculations, or cancel allocations.

Defaults can also live in the ordinary single KinBot JSON input:

```json
{
  "queuing": "slurm",
  "queue_live_policy": {"max_jobs": 8, "manage_pending": true},
  "queue_live_control": {
    "hold_failures": false,
    "resource_rules": [
      {"name": "long-qc", "selectors": ["*/tasks/*"],
       "resources": {"walltime": "7-00:00:00"}}
    ]
  }
}
```

These are initial defaults. Later live edits take precedence, including after a
normal restart. Empty `queue_live_control` retains the existing restart and
failure behavior. Nested control roots are rejected rather than splitting
ownership. PBS/local operation is unchanged; these new allocation controls are
for Slurm.

## Shared capacity and nodes

```bash
python -m kinbot.run_control capacity "$KB_RUN" --max-jobs 8
python -m kinbot.run_control capacity "$KB_RUN" --max-jobs 8 --apply
python -m kinbot.run_control capacity "$KB_RUN" --max-jobs 0 --apply
python -m kinbot.run_control capacity "$KB_RUN" --max-jobs 12 --apply
python -m kinbot.run_control capacity "$KB_RUN" --exclude nodeA,nodeB --apply
python -m kinbot.run_control capacity "$KB_RUN" --allow nodeA --apply
```

Omitting `--apply` previews capacity changes. Lowering capacity retains running
allocations and holds excess pending jobs. Zero pauses admissions and pending
backfill, without stopping active QC. Increasing capacity releases KinBot's
recorded holds and admits more ready work on the next controller poll. Slurm
availability and priority still determine start times; this is not a promise to
start immediately. Holds and releases are verified. A race with job startup
retains that running calculation. Pre-existing user/administrator holds are not
released. Do not independently manage holds on allocations already held by
this control interface.

Every cooperating branch shares the same serialized admission ledger, including
uncertain submissions. Pending and dependency-pending allocations reserve
capacity conservatively; they are not discounted as hypothetical free nodes.
The pool still enforces single-node exclusive allocations and existing scratch
policy. Node exclusions apply to new submissions; replacing selected pending
jobs is explicit, as described below. These commands do not remove unrelated
exclusions, change chemistry, or change scratch selection.

## CPU, memory and walltime rules

```bash
python -m kinbot.run_control resources "$KB_RUN" --name l3 \
  --select '*/tasks/l3_geometry*' --cores 8 --memory-mb 96000 --walltime 10-00:00:00
python -m kinbot.run_control resources "$KB_RUN" --name one-job \
  --select 'case/well/conf/well_0003' --cores 4 --memory-mb 24000
python -m kinbot.run_control resources "$KB_RUN" --name all-jobs \
  --select '*' --walltime 7-00:00:00
python -m kinbot.run_control resources "$KB_RUN" --name one-job --remove
```

Rules are published immediately; their revision/history is recorded. Selectors
are shell-style globs relative to the control root (`*` may span directories),
and can be repeated. Later matching rules override individual resource fields.
Inspect matching targets with `list` before publishing a broad rule. Removing a
rule restores base settings for unsubmitted dispatcher work; it never edits a
completed/running allocation. Set explicit desired resources when retrying an
ordinary worker whose previous retry already used an override.

Memory is the total allocation in MB, not Molpro memory per rank. Native memory
budgets and thread/rank counts are recalculated together with the Slurm request.
Configured per-node limits and backend memory checks still apply. Unsupported
fixed custom native decks are rejected rather than receiving a scheduler-only
edit. Resource rules support Gaussian/ORCA ordinary workers, Molpro/Gaussian
optimizations, and Molpro/MRCC/CFOUR/Gaussian external dispatcher tasks. Generated
Molpro ROTD sampling and MESS batches also read rules before submission. ROTD
freezes its native resource settings into that allocation without altering
sampling inputs, cached energies, active spaces, or the pending sample pickle.

Already-submitted jobs keep their resources. For a pending job, set the rule
and request a pending-only replacement. Running jobs need an explicit
`--cancel-running` request to use a new allocation. Ordinary QC reads the rule
again when shared admission becomes available. Dispatcher edits of unsubmitted
inputs are archived, and resource-only changes retain scientific VRC identity.
If a graph still contains old workers without the new manifest protocol,
resource edits are deferred until those workers finish.

An existing shared `minimum_walltime` remains a floor for new submissions. Use
`python -m kinbot.run_policy set "$KB_RUN" --clear-minimum-walltime` deliberately
if that floor should no longer apply. An explicit pending-time update is not
subject to this submission floor.
Confirmed automatic timeout escalations are also preserved when resource rules
are reapplied or removed; a shorter stale rule must not shrink a retry budget.

## Pending walltime without replacing the calculation

```bash
python -m kinbot.run_control pending-time "$KB_RUN" --job JOB_ID --walltime 10-00:00:00
python -m kinbot.run_control pending-time "$KB_RUN" --job JOB_ID --walltime 10-00:00:00 --apply
python -m kinbot.run_control pending-time "$KB_RUN" --all-pending --delta=+02:00:00 --apply
python -m kinbot.run_control pending-time "$KB_RUN" --job JOB_ID --delta=-01:00:00 --apply
```

Repeat `--job` to select several allocations. All limits must be positive whole
minutes and fit their current partition. Each job is held and rechecked before
`scontrol update`; the result and release are verified. A job that starts first
is retained. Scheduler permission denial is reported, not treated as success.
The receipt path is printed; `resume-scheduler RECEIPT.json` reconciles an
interrupted update using its recorded absolute limit, without applying a delta
twice. Running-job extensions remain subject to Slurm permissions and are not
attempted by this command.

## Targeted retry requests

```bash
python -m kinbot.run_control retry "$KB_RUN" \
  --target 'case/well/anl_automatic/states/state_ID/composite/tasks/l3_geometry' \
  --reason 'Reviewed native failure; try again with the revised resources'
# Add --apply to publish the request after inspecting its generation.
```

Repeat `--target` for multiple explicit targets. Each request is independent and
has its own receipt. Kinds and target names are:

* `--kind task` (default): `GRAPH/tasks/TASK_ID`, including ANL and VRC tasks.
* `--kind qc`: worker basename without `.py`, including L1/L2 and rotor workers.
* `--kind reaction`: rejected reaction's exact name inside the active well.
* `--kind conformer`: rejected conformer basename, including ring conformers.

Use `--pending-only --apply` to replace only a still-pending attempt, or
`--cancel-running --apply` to explicitly stop and retry a running allocation.
Without either option, a live allocation is retained. Attempts are fenced by
saved input hashes, allocation IDs and/or dispatcher attempt numbers. A request
for an obsolete generation cannot restart its replacement. An uncertain
scheduler state never authorizes another submission. Successful native results
that arrive during cancellation are retained for normal reconciliation.

Manual retries permit another attempt after an otherwise ineligible failure or
exhausted automatic retry budget; they do not reset that budget. Old artifacts
are archived before restaging. Ordinary QC restores verified pre-allocation
seeds, fences old database rows, and uses durable submission receipts. These
are clean retries, not a promise of native partial-iteration restart. Existing
native restart and automatic timeout/I/O behavior remains available separately.
All geometry, curvature, connectivity, electronic-state and result-provenance
checks run again. A retry never marks a rejected calculation as successful.

Task retries refuse to invalidate already staged/running/completed dependent
work. A rejected reaction can also be retried after its search finishes while
its revision-aware PES owner remains in exploration/postprocessing. The owner
first quiesces downstream dispatch across that PES, then reopens exploration.
Running native allocations continue. Closed conformer searches and finalized
PES/ROTD/MESS campaigns still require dependency-aware recovery.
Reaction retry retains shared product wells and other reaction namespaces.
Missing allocation ownership or missing ordinary-worker retry seeds require
reconciliation, not a guessed resubmission. Custom ordinary backends outside
Gaussian/ORCA can use walltime rules but do not receive guessed CPU/memory edits.

## Keep a retry window open

```bash
python -m kinbot.run_control enable "$KB_RUN" --hold-failures
python -m kinbot.run_control release "$KB_RUN" --kind task \
  --target 'GRAPH/tasks/TASK_ID' --reason 'Reviewed; continue normal failure handling' --apply
python -m kinbot.run_control enable "$KB_RUN" --no-hold-failures
```

With `hold_failures`, drivers remain available after failures, and unstarted
dependent fallback work waits. Independent ready work continues. Ordinary QC
and search rejections also wait for retry or release; this includes expected
unproductive reaction guesses, so enabling it can intentionally delay normal
network completion. `release` restores normal handling for that exact failed
generation; it does not waive required L3 success. The default is off.

A dead/exited controller cannot consume a live request. Refresh it once using
the existing recovery tools. For the already-running two-PES salvage layout,
`run_refresh RUN --apply` replaces controllers and preserves allocations; do not
add pending replacement or VRC migration flags simply to enable these controls.
Configuration/history, request receipts, scheduler receipts and archived
attempts remain in the run directory for review. Nothing in this interface
assumes a username, chemistry example, hostname, partition, or scratch path.

## Exploration and reviewed network extensions

`pes_phase_policy` defaults to `auto`: automatic CBH/ANL and VRC/ROTD PES runs
finish exploration of all eligible wells before downstream dispatch begins.
`simultaneous_kinbot` limits exploring wells, not parked controllers. Thus a
value of one permits serial exploration followed by concurrent downstream
work. Ordinary runs without these branches retain their prior behavior.
Explicit `legacy` or `exploration_first` policies remain available. A live
revisit is supported only by an exploration-first owner.

```bash
python -m kinbot.run_control revisit "$KB_RUN" --worker 'case/WELL' \
  --vdw-detection 0 --rejected --reason 'Reassess complex cutoff and rejected channels'
# Repeat with --apply to request it.
```

This grants one additional attempt to each currently rejected channel. It
replays accepted results through the ordinary validation path, retaining shared
product calculations. It does not waive connectivity, curvature, stereochemical
or scientific acceptance checks. Zero is a well-depth cutoff, not an instruction
to accept unbound structures. A scientifically changed previously accepted
geometry requires explicit dependency invalidation; the interface refuses to
silently reuse its old single points.

The ANL plan may grow with new species, TSs, CBH references and accepted vdW
complexes. Unchanged graph identities keep complete/submitted QC. Derived
profiles, formations and handoff reports are archived and rebuilt. A reference
promoted to a network species with a new L2 geometry gets a separate graph;
single points from the old geometry cannot become its results. Disconnected
complexes use same-method binding energies plus fragment formation enthalpies,
not an invented covalent CBH decomposition.

## Resubmit the same JSON without a second driver

New parent wells cannot inherit the initial well's orbital indices. Missing
surface models create `vrc_model_review.json` and per-channel
`vrctst/orbital_review/REACTION/request.json` records with parent geometry,
products, state and bond indices. The unreviewed channel starts no fragment/scan
QC. Other reviewed channels and ANL work proceed. These records request orbital
reconnaissance and human review; they are not automatically approved active
spaces or proof that native reconnaissance has run. Final network assembly and
PES sampling wait for all required corrections, including these holds.

After reviewing the native orbitals, add `vrc_tst_surface_models` to the original
PES `input.json`. It maps parent routing names (recommended) or SMILES to the
existing schema-1 reviewed surface model format. Each new parent's bottleneck
must carry its own `parameters.vrc_tst_state_model` (including basis models,
active orbital indices, state-average weights and seed provenance), or an
explicit reviewed state-model file. Inline models keep the scientific settings
in the single input. Required native WFU seed files remain explicit hashed
artifacts. Choose scan/sampling/reference methods, bases, extrapolation and spin
settings through that bottleneck's existing parameter overrides.

```bash
python -m kinbot.run_control resubmit /absolute/path/to/case/input.json
python -m kinbot.run_control resubmit /absolute/path/to/case/input.json --apply
python -m kinbot.run_control resubmit /absolute/path/to/case/input.json --status
```

An input stored outside the run directory remains supported. Run resubmission
from its PES directory, or add `--run-directory /absolute/path/to/case`.
Accepted snapshots retain the original directory for resolving orbital paths.

An unchanged resubmission is a no-op. Re-running `python -m kinbot.pes INPUT`
also detects its live owner and uses this interface instead of starting another
driver. Explicit `run_control resubmit` is preferable for preview and receipts.
The existing owner consumes a durable batch, acknowledges each worker, and
re-enters downstream dispatch. It does not replace controller PIDs. The
receipt distinguishes `queued`, `published`, `applied` and `blocked`; publication
alone does not establish successful worker application.

Only `vrc_tst_surface_models` additions and explicit `vdW_detection` updates are
accepted by this scientific file interface. Existing reviewed channels must
remain unchanged. Live CPU/memory/walltime/capacity edits use the resource
commands above. Unsupported scientific edits are refused with an explanation;
they never silently invalidate expensive work. Resubmit after exploration has
completed and while its PES owner remains alive; a finished controller cannot
consume a request. Editing a file without resubmission has no effect on the
running owner or new wells: they read the last accepted snapshot.

Input publication is serialized with the PES closing/reopening barrier and
shared control requests. Partial batches are invisible until publication
completes. Receipts are idempotent and preserve exact input/prepared-manifest
hashes. The review hold, exploration epoch and task results are visible in the
ordinary monitor. This interface does not start native QC on the login node.
