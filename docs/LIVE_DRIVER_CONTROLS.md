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
work. Search retries are accepted only while that search still owns the
network; once results have gone downstream, use dependency-aware recovery.
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
