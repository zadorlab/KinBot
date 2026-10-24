# Live driver controls

These opt-in Slurm controls change scheduling and request another calculation
attempt, including operator-specified numerical fallbacks, without editing source code. A controller must
already include this implementation. An older running Python process needs one
controller refresh; subsequent control commands do not require refreshes.
Running allocations keep their immutable inputs and code snapshot.

After an outage, `kinbot.run_refresh ROOT --apply` archives the previous PES
coordination records after verifying controller shutdown and acquiring the
worker locks. This includes an empty or truncated `pes_exploration.json`:
its exact bytes remain in the controller archive, while replacement owners
create their new generation. Discovery inventories, saved search decisions,
accepted chemistry, dispatch ledgers and live allocations stay in place.
Independent ROTD campaigns retain their own ownership: resume an absent owner
with `kinbot.rotdpy_campaign retry INPUT.py`. The command verifies the immutable
input/correction contract and retains an existing owner; valid completed
results are reused and existing sampling databases remain available to ROTD.
PES and sampling ownership tolerate a hostname changing between its short name
and fully qualified name only when the recorded boot ID, PID and process start
time all match. Receipts are not rewritten to hide a mismatch. Different boot
IDs, reused PIDs, unrelated hosts and changed scientific contracts still stop
adoption or signalling.

For a complete retained-run handoff from a checkout, use
`bash examples/anl/tools/refresh_retained_run.sh ROOT MAX_JOBS`.
It adds currently down/unknown/draining/failed nodes to the existing exclusions,
updates the shared admission cap, refreshes PES/KB owners, resumes eligible
independent ROTD owners, then reconciles saved dispatcher results and applicable
automatic retries without submitting QC itself. It checks PES owner health
before and after reconciliation and writes a timestamped reconciliation report.
Review holds and scientific validation remain in force. The cap is a maximum,
not a guarantee of a fixed number of idle nodes.
The handoff also verifies each generated RotdPy setup and routes future sample
workers through the new immutable code snapshot. Retained owners, input files,
sampling databases and already executing workers are unchanged. The cached
Slurm template sources this setup on each allocation. The route is recorded in
`rotdPy/.kinbot_worker_code.json`, and prior setup files are archived. A later
runtime-only update can use `python -m kinbot.rotdpy_worker_code ROOT PIN` without
replacing sampling or PES owners; `PIN` must be a verified snapshot in ROOT.

Native CC iteration exhaustion now has a bounded automatic ladder in ANL
dispatch: Molpro CCSD(T)/CC-F12 optimization and single-point jobs use MAXIT
100 then 200; CFOUR VCC uses CC_MAXCYC 600 then 1000. An existing larger limit
is never reduced. Native CC-specific failure evidence and the failed attempt's
identity are required. Failed attempts are archived, blocked descendants are
reconciled, and geometry, theory, state definitions and convergence thresholds
stay fixed. MRCC's `ccprog=mrcc` higher-order engine has a compiled iteration
ceiling, so an ineffective `ccmaxit` edit is not attempted there.

The Molpro sampling bridge shares the dispatcher's two-step `hsdel` singularity
recovery: enlarge the numerical P-space, then reduce the iterative workspace
and tighten numerical convergence. It keeps the same sample geometry, selected
roots, state-average weights, active space, level shift and IPEA shift. Native
attempts and restart records are archived under `.kinbot_sample_attempts` in
the worker directory; exhausted settings stay exhausted across worker restarts.
Normal termination and the energy parser remain mandatory. Unknown failures
remain failures, and failed samples are never assigned replacement energies.
Sampling ownership is checked independently of current scan readiness: replay
of saved scan points can temporarily publish `running` during replacement.
A live, hash-bound campaign is retained through that interval. An absent owner
waits for a complete, resolved scan and its correction binding before launch;
the handoff records `resume_requested` so normal VRC owner staging consumes
that request once readiness returns. A real review hold remains held and is
reported in `sampling_resume.json`; it does not abort other campaigns or the
reconciliation and health checks. Completed sampling files stay in place even
when their use is held. Input/correction changes and ownership mismatches still
stop the handoff, rather than authorize reuse of different science.
If an earlier handoff is still active, pass its PID and exact stdout-log path as
the third and fourth arguments. The script waits for the verified process and
any existing reconciliation writer; it does not terminate an I/O-blocked writer.
Use `nohup` with an explicit log to keep this handoff running after logout.
An optional `KB_DECISIONS_FILE=/absolute/path/decisions.json` adds owner-applied
reaction decisions after loading the replacement code. The JSON is a list of
objects with exactly `worker` (relative to ROOT), `reaction` (exact accepted
name), `decision` (`reject` or `include`), `reason`, and boolean `cancel_running`.
The whole batch is validated before policy changes or controller replacement.
The handoff waits for these receipts before resuming absent sampling owners and
running the full reconciliation. Re-running a completed decision is idempotent.
The complete handoff gives local child processes, pool admission, graph staging
and submission bookkeeping up to 30 minutes to settle before replacement.
`run_refresh --wait-idle SECONDS` exposes this bounded wait directly. Each retry
inventories all owners again and logs the blocking process or task. It never
signals a child, resets a submitting task, or bypasses a scientific/ownership
validation failure. A race during quiescence resumes the original controllers
and releases the admission lock before waiting, so the owner can finish its
write. An unresolved submission ends the wait with its exact ledger path;
refresh never guesses its allocation ID or submits another copy.

## Persistent reaction selection

Once these controls are loaded, a user decision requires no controller restart:

```bash
python -m kinbot.run_control decide-reaction "$KB_RUN" \
  --worker case/worker --reaction exact_accepted_reaction_name \
  --decision reject --cancel-running --reason 'Reviewed relaxed path changes fragments' --apply
```

Omit `--apply` for a read-only preview. The network owner waits for its branch
barrier, archives the previous selection and affected dispatch ledgers, cancels
only verified dedicated VRC tasks/sampling allocations, then excludes the
channel from VRC and final PES selection. Shared fragment and ANL calculations
remain available. Native scan evidence and raw accepted L2 chemistry stay in
place. A partial correction batch retires only that channel's tasks.
An interrupted cancellation remains an applying request until its allocations
finish. Search checkpoints, including rejected searches, are not reset.

`--decision include` removes this operator exclusion and restores the saved
task ledger for reconciliation. It does not accept failed or unknown chemistry,
waive scan validation, or invent a converged result. Previously cancelled work
may still need an ordinary explicit task retry. Operator-excluded sampling
owners cannot restart through the campaign launcher or retained-run refresh.

## VPT2 geometry recovery and correction-energy reuse

The shared Gaussian VPT2 framework recovery exhausts its applicable symmetry
settings using native evidence and saved attempt history. For a declared minimum
with a supported generated Gaussian route, it then makes one isolated tighter
L2 geometry branch at the same method, basis, dispersion and integration grid.
The force criterion is at most `1e-5` eV/angstrom (at least ten times tighter than
an explicitly declared source criterion), and SCF convergence is at least 10.
A separate ordinary L2 harmonic job must produce the complete positive spectrum
before regular VPT2 starts at that geometry. The fresh geometry generation gets
its own full symmetry ladder, including applicable native-group/Loose, C1,
NoSymm and mass-centered alternatives. Exhaustion then marks VPT2 unavailable
even when the general task-failure hold is enabled. Unknown failures and failed
refinement/Hessian validation remain held for diagnosis.

The old L2 geometry and all L3 geometries/single points keep their inputs and
ownership. The L2-only fallback uses the refined electronic energy and harmonic
ZPE together; existing L3-SP@L2 components retain their original geometry and
harmonic reference. ANL/CBH profile assembly sees an unavailable anharmonic
component explicitly and uses its existing common harmonic fallback machinery.
No anharmonic term is silently replaced by zero in an anharmonic recipe.

VRC assembly can reuse a completed relaxed-scan central energy for any matching
correction component. Matching includes geometry/atom order, method, basis,
charge, multiplicity, electronic model, shift/IPEA, and the hash-bound propagated
orbital lineage. Native final input/output, final geometry and execution records
are retained as the component provenance. A different frozen geometry requires
its own sampling-level calculation. Already dispatched single points are
adopted, and existing completed outputs take precedence over reuse.

For a direct correction the arithmetic remains
`[E_high,relaxed(R) - E_high,relaxed(infinity)] -
[E_sample,frozen(R) - E_sample,frozen(infinity)]`.
Both curves use the same separation grid; each has its own asymptotic reference.
The sampling level may itself contain multiple CBS components. Relaxed-scan
readiness and immutable correction/input provenance gate new ROTD sampling and
use of sampling results; an
additional manual approval of the assembled correction is not required.
Actual converged sampling results are still required before rate integration.

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

The monitor labels the driver's `vrc_tst` branch as `vrc_scans/corrections`:
completion there means preparation and sampling handoff have finished. It lists
RotdPy campaign statuses separately, linked to each worker by the saved correction
source. Missing campaign evidence is not reported as completed sampling. An
operator-rejected scan's old errors and unresolved points remain in `--details`
and `--json`, but are excluded from active scan failure/review summaries. These
display changes require only an updated monitor, not replacement controllers.

Defaults can also live in the ordinary single KinBot JSON input:

```json
{
  "queuing": "slurm",
  "queue_live_policy": {"max_jobs": 8, "manage_pending": true},
  "queue_live_control": {
    "hold_failures": false,
    "hold_task_failures": true,
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

New jobs default to `partition_max`, resolved from the selected Slurm
partition's maximum time. Explicit finite limits in the input, task resources,
or the rules below override that default. The fresh r7 campaign has no minimum
walltime floor. Removing a rule restores the prepared task's base allocation;
it does not alter existing scheduler allocations or silently reprepare QC.

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

## A manual fallback for one attempt or an explicit batch

Write a JSON proposal and pass `--fallback FILE.json` to `retry`. Its contents
are embedded in the request, so editing the proposal afterward cannot change
an already queued retry. This does not reset automatic retry budgets. Repeat
`--target` for several jobs of the same kind; all proposals are validated before
publication, and each owner then rechecks its individual attempt. A batch is
not an all-or-nothing cancellation transaction: a target whose generation
changed is retained and reported separately.

For an ASE optimization task (ANL geometry, MR relaxed point or asymptote):

```json
{
  "schema": 1,
  "calculator_kwargs": {"scf_convergence": 1e-8},
  "optimizer": {"steps": 300, "sella_kwargs": {"internal": false}},
  "policy": {"optimization_step_retries": 2},
  "resources": {"cores": 8, "memory_mb": 96000, "walltime": "10-00:00:00"},
  "environment": {"CUSTOM_RUNTIME_OPTION": "value"}
}
```

For an ordinary Gaussian L1/L2 worker (`--kind qc`), a rejected conformer,
or the newly generated workers of a rejected reaction:

```json
{
  "schema": 1,
  "worker_parameters": {
    "kwargs": {"scf": "conver=8,CDIIS,NoVarAcc,NoIncFock,MaxCycle=512", "guess": null},
    "steps": 400
  },
  "environment": {"GAUSS_MEMDEF": null, "GAUSS_PDEF": null},
  "resources": {"cores": 8, "walltime": "7-00:00:00"}
}
```

Use only variables present in the actual generated worker. Supported literal
assignments are `kwargs`, `kwargs1`, etc., `prod_kwargs`, `sella_kwargs`, `steps`
and `fmax`. Dictionary values merge; `null` removes a keyword. Missing or
ambiguous assignments are refused before cancellation. A reaction-wide
fallback applies to its subsequently generated namespace; shared product wells
are unchanged and need their own explicitly targeted fallback if necessary.

External native tasks accept `input_edits: [{"old": "exact text", "new":
"replacement text"}]`; each old string must match exactly once. They also
accept `command` as an argument list, preserving the task's native filenames
and placeholders. ASE tasks instead accept a command string through
`calculator_kwargs.command`. For example, an installed replacement executable
can be selected without changing the KinBot driver. Exact native edits and
custom executables are operator-reviewed numerical repairs: KinBot cannot
prove arbitrary native text has identical chemistry. Do not use them to
change basis, method, states, geometry, or the meaning of a parsed energy.
Those changes require a scientific migration and dependent-result invalidation.

Both kinds accept `setup`, a list of single shell command lines, and
`environment`, mapping names to string values or `null` to unset. These commands
run **inside the replacement allocation only**, after site setup. Manual task
runtime setup is reported as `runtime_environment_deferred` by preflight;
its executable availability is established in the allocation. Operator setup
is executable shell code, not a sandbox. Environment values are quoted as data.
Scheduler/account identity variables and Python import-path overrides are
rejected. Backend `KINBOT_*_EXECUTABLE`, `_ROOT` and `_SETUP` overrides are
supported. Scratch locality and native-result validation still apply.

```bash
python -m kinbot.run_control retry "$KB_RUN" \
  --target 'case/well/anl_automatic/states/state_ID/composite/tasks/l3_geometry' \
  --fallback /absolute/path/to/fallback.json --reason 'Reviewed numerical fallback'
# Inspect the preview, then repeat with --apply.
```

Task options are validated before cancellation. An unsubmitted staged task can
also receive a fallback through this command; no allocation is cancelled.
Running and pending tasks retain the explicit stop rules described above.
The archive contains the old attempt; manifests record the original logical
task, the patched task hashes and the embedded fallback. Unrelated allocations
and completed files remain untouched. Later VRC model refresh recognizes the
approved fallback as the existing task rather than submitting a duplicate.
Further manual retries use the effective previous worker, including prior
fallbacks, rather than silently reverting to its initial source.

The per-attempt `optimization_step_retries` override governs additional chunks
inside that replacement worker. Campaign timeout count/factor, minimum
walltime and orbital-review reserve are live `kinbot.run_policy set` options;
they are distinct from native SCF/CC iteration limits. A new executable must
remain output-compatible with the selected parser. Successful completed tasks,
uncertain submissions, and changed sources with existing dependent work are
not overwritten through a numerical retry.

## Keep a retry window open

```bash
python -m kinbot.run_control enable "$KB_RUN" --hold-failures
python -m kinbot.run_control enable "$KB_RUN" --no-hold-failures --hold-task-failures
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

`hold_task_failures` keeps only dispatcher ANL/VRC failures open for manual
retry/release. Ordinary scientifically rejected reaction searches can finish,
so they do not indefinitely block outer-PES selection. Fresh r7 enables this
task-only option. Independent ready tasks continue while a failed task waits.

A dead/exited controller cannot consume a live request. Refresh it once using
the existing recovery tools. For the already-running two-PES salvage layout,
`run_refresh RUN --apply` replaces controllers and preserves allocations; do not
add pending replacement or VRC migration flags simply to enable these controls.
Configuration/history, request receipts, scheduler receipts and archived
attempts remain in the run directory for review. Nothing in this interface
assumes a username, chemistry example, hostname, partition, or scratch path.

Each well now saves reaction decisions and constrained-path progress in
`.kinbot_search_state.json`. A successful intermediate QC calculation cannot
reopen a rejected search after a controller refresh. Accepted reactions rebuild
their in-memory objects from retained QC, with replacement submissions blocked;
the monitor labels this as `restoring accepted QC`. The outer PES can briefly
report `exploring` while those objects and its generation-bound L2 selection
are reconstructed. That does not authorize new reaction attempts.

Path submissions record their next step at scheduler acknowledgement, and scan
energy histories survive replacement without counting a point twice. Explicit
reaction retries archive the previous decision and reset only their target.
Input/theory or reaction-definition mismatches require an explicit revision.
An authorized network-policy revision can revalidate accepted results and
perform the additional work that policy requires.

When upgrading an older run, `run_refresh` captures the existing search monitor
while controllers are frozen, archives it, and installs the checkpoint after
the controllers stop. Admission receipts distinguish the most recent path
submission from a rendered input still waiting for capacity. An old partial
scan without its energy history, or a path with ambiguous step provenance,
requires review rather than guessing. Jobs already submitted by an older
controller remain owned and are retained. Loading this fix does not undo
earlier extra attempts or waive scientific review holds.

PES startup and rediscovery retain existing worker inputs verbatim, including
the initial well's starting coordinates. Re-embedding the same SMILES during
startup must not rewrite a completed search's input. Checkpoints retain an input
snapshot, and refresh checks compatibility before replacing controllers. For
older checkpoints affected by an overwritten starting structure, refresh can
restore the structure from the archived worker parameter log only when doing so
reproduces the full saved input SHA-256 and the theory fingerprint still matches.
The overwritten input is archived; search decisions, step history, native files
and job IDs are retained. Other mismatches still require an explicit revision.

## Exploration and reviewed network extensions

`pes_phase_policy` defaults to `auto`: automatic CBH/ANL PES runs finish
exploration of all eligible wells and outer L2 pathway selection before new
ANL planning/dispatch begins. Reviewed VRC channels may start once their own
well finishes exploration, while other wells continue exploring. Final network
assembly still waits for the entire exploration and both branches.
`simultaneous_kinbot` limits exploring wells, not parked controllers. Thus a
value of one permits serial exploration with independent VRC work. If that
one exploring well is waiting on its last QC job, the next well still waits;
raising a shared node cap does not increase well concurrency. Ordinary runs
without these branches retain their prior behavior.
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
they never silently invalidate expensive work. Additive VRC reviews can be
resubmitted while other wells are still exploring: ready branches receive
generation-fenced requests, and searching/newly discovered wells inherit the
accepted models at their ready boundary. A worker becoming ready between
submission and consumption is included; conflicting search revisions block
publication. `vdW_detection` file revisions still require the global
postprocessing phase; use an explicit network revisit for exploration changes.
A finished controller cannot consume a request. Editing a file without resubmission has no effect on the
running owner or new wells: they read the last accepted snapshot.

Input publication is serialized with the PES closing/reopening barrier and
shared control requests. Partial batches are invisible until publication
completes. Receipts are idempotent and preserve exact input/prepared-manifest
hashes. The review hold, exploration epoch and task results are visible in the
ordinary monitor. This interface does not start native QC on the login node.
An `applied` model-only receipt acknowledges the currently targeted branches;
still-exploring wells inherit that accepted input when they become ready.

These controls cover configuration-driven retries, not arbitrary changes to
already imported Python code. A future controller-code defect can still require
a controller replacement; unsafe native checkpoints and contradictory scientific
provenance still require review. No live interface can guarantee every backend
will converge or implement an unsupported native method.

## Reconcile saved results and known failures without replacing controllers

```bash
python -m kinbot.run_control reconcile "$KB_RUN"
python -m kinbot.run_control reconcile "$KB_RUN" --apply
```

The first command previews current dispatcher graphs. `--select` accepts a
repeatable graph-path glob relative to the run root, for example
`'ch3s/*/anl_automatic/states/*/composite'`. The applied command takes each
unowned graph's existing drive lock, reconciles native results, applies known
bounded recoveries, and stages eligible tasks. It **does not submit or cancel
QC jobs**, launch drivers, or release the outer-PES L2 admission gate. An
actively owned graph is reported as `owned_by_controller` and left to that
controller. A failed graph is reported as `blocked`; other graphs are still
processed. Retry archives and journals are excluded from the inventory.

All Gaussian VPT2 dispatcher tasks now use the same native-diagnostic recovery
rules, regardless of task name, molecule, ANL variant or driver:

- A named old/new framework-group mismatch can receive `Symmetry=None`.
- A subsequent unassigned old group under `Symmetry=None` can receive one
  `Symmetry=(PG=<native New group>)` retry. The October 22 CHS case reports
  `New: CS - Old:` and therefore selects `Symmetry=(PG=CS)`.
- For verified linear geometries, new VPT2 input rendering retains point-group
  detection with `NoSymm`. An older unsubmitted input is restaged; a failed
  old input is retried only for the native `Molecule is non-linear but NTrRo=5`
  diagnostic. This interface repair is recorded once per task.
- Simple and compound anharmonic frequency keywords are recognized consistently
  by input validation, rendering, recovery and the native result reader.

These are finite recovery paths, not a generic segmentation-fault retry.
Native failure evidence, input/geometry hashes and allocation termination are
checked. Existing completed/allocated descendants prevent automatic source
replacement. Interrupted archival resumes through the live-control journal;
unchanged allocations retain their manifests and files. Legacy allocations
that lack the manifest read lock defer a workflow edit until they finish.
Unknown errors or unsupported custom symmetry combinations remain available
for a reviewed `retry --fallback FILE.json` request. The geometry, electronic
method, basis and result acceptance criteria are not relaxed by these repairs.

`NoSymm`, `Symmetry=None` and `PG` have different meanings in the
[Gaussian symmetry manual, distributor-hosted translation](https://conflex.co.jp/gaussian_support/symmetry.php).
`PG` caps the computational group; it does not prove convergence or impose a
new accepted chemical structure. Native completion and VPT2 result checks
remain required after recovery.

For nonlinear framework-group failures, the established sequence includes
`NoSymm` after `PG=C1`, even when the log lacks the additional mass/Eckart
diagnostics. It retains group detection without initial reorientation. If the
complete native log still shows a different initial full point group and
mass-based framework at the Eckart handoff, one additional
`Symmetry=(PG=<native New group>,Loose,CenterOfMass)` trial is applicable.
This selects the documented mass center instead of the default charge center.
Both require the specific framework mismatch, preserve the submitted geometry
and theory, and retain native validation. These are bounded recovery attempts,
not a guarantee that Gaussian will complete. No particular molecule or group
is used to choose them.

Explicit task retries also support a narrow lost-before-output case when
Slurm has purged the allocation and accounting is unavailable. The owner
requires the original admission receipt, unchanged staged files, expiry of
the result-publication grace, no worker or scheduler output, an explicit
`scontrol` invalid-ID response, and a successful queue query with that job
absent. An empty best-effort query or a scheduler outage is insufficient.
The request retains this evidence and records the original exit as unknown;
it does not synthesize a failed native result or infer a node failure. This
path requires an explicit generation-fenced request and does not create an
automatic repeated-resubmission loop. Partial native work needs its own
diagnosis. Existing descendant and archival safeguards still apply.

## Preserve and restore the discovered-well inventory

PES startup preserves `chemids` instead of resetting it to the initial well.
Controller handoff also archives `chemids` and `pes_exploration.json`. When
upgrading a run affected by the older truncation bug, `run_refresh` restores
verified missing workers from the preceding handoff's `pes_summary.txt` before
relaunch. `controller_refresh.json: discovery_restored` records the additions.
Existing worker inputs must match the well name and current run format.
Discovery restoration does not manufacture L2 completion or submit QC.

An existing exploring owner can also accept an explicit saved inventory:

```bash
python -m kinbot.pes_inventory "$KB_RUN/ch3s" \
  --summary /absolute/path/to/saved/pes_summary.txt
# Inspect added/retained names, then repeat with --apply.
```

The merge is append-only, atomic and phase-locked. Normal exploration limits
admit the added workers, and the full L2 selection gate remains enforced.
Restoring additional wells after a live owner entered downstream planning
requires the existing exploration-revisit procedure or a controller handoff.
Stray numeric directories are never treated as accepted discoveries.

Installing these code fixes into an already pinned old controller needs one
controller replacement. Once installed, known recovery policies, numerical
keyword/command/environment overrides, resources, explicit retries, and VRC
review/resubmission use the live interfaces above. A genuinely new software
defect can still require a software upgrade; arbitrary scientific changes
must invalidate affected dependents rather than silently reuse them.

## Paired geometry and frequency refinement

A failed geometry, harmonic calculation, or dependent single point can be
retried with tighter geometry and electronic convergence through its existing
owner. This applies to TZ and QZ stages, minima and transition states, without
molecule, reaction, task-name or basis-name rules. The present native adapter
supports generated Molpro ASE geometry/harmonic graphs; unsupported custom
frequency decks are refused before any job is stopped.

```bash
python -m kinbot.run_control refine "$KB_RUN" \
  --target RELATIVE_GRAPH/tasks/harmonic_tz \
  --fmax 0.0001 --scf 1e-11 --cc 1e-11 \
  --reason 'Refine geometry and numerical Hessian after curvature validation'
# Inspect the affected tasks, then repeat with --apply.
```

`--target` can instead identify the QZ geometry, its frequency task, or a
single point consuming that geometry. The command resolves the geometry source
and the entire dependency closure (`geometry_from`, `geometry_fallback` and
`depends_on`). All existing results in that closure, including completed
single points, are archived and recomputed. Tasks outside it retain their files
and allocation IDs. A QZ refinement does not reset its TZ ancestors. A TZ
refinement resets QZ geometry and everything downstream of it.

The default waits for affected allocations to finish. Add `--cancel-running`
when publishing the request to stop only the reviewed affected allocations.
Job ownership, task/geometry hashes and request generations are checked before
cancellation. The controller journals archival and replacement, resumes an
interrupted transaction, and blocks a stale request instead of affecting a
newer attempt. Pending refinement fences new submissions using the old source.
Verified converged geometry can seed its own refinement; descendant geometry
optimizations always consume the replacement parent, not their previous seed.

The force target is in eV/angstrom. `--scf` and `--cc` tighten the native
calculations used for both geometry derivatives and the numerical Hessian;
stricter existing thresholds are preserved. The Hessian finite-difference step
and imaginary-mode acceptance rule are unchanged. Molpro's default step is
0.01 bohr; making it smaller is not automatically more accurate when energies
are noisy. The retry is a new calculation, not a claim that extra negative
curvature is numerical error. See the [Molpro frequency manual](https://www.molpro.net/manual/doku.php?id=harmonic_vibrational_frequencies_frequencies).

New requests use a distinct `refinement` kind so older controllers cannot
misinterpret them as single-task retries. After the one code upgrade, subsequent
requests require no controller replacement. A control-enabled ANL graph owner
now remains available after hard failures instead of exiting once other work
finishes. `live_retry_status.json` reports the failure; no downstream composite
handoff occurs while it remains unresolved. Final handoff is serialized with
request publication; already handed-off graphs require reopening the owning
workflow, because changing a consumed result cannot safely update external
consumers by itself. Requests still respect PES admission and shared capacity.

## Native-derived VPT2 recovery

The existing automatic framework-error ladder now also handles this situation:
Gaussian was asked for the group printed as `New`, but its initial `Full point
group` agrees with `Old`, and `Framework definition based on atomic masses`
agrees with `New`. It retries once with `Symmetry=(PG=<native New>,Loose)`.
The group is read from the native log, not selected from molecule-specific or
point-group-specific cases. Attempt history prevents cycling. Existing geometry,
method, basis and VPT2 validation are retained. This is the documented symmetry
detection tolerance, not permission to accept an unconverged frequency job.

Unknown failures and contradictory native group assignments remain visible.
A segmentation trace by itself is not a retry trigger. This procedure cannot
infer a reliable numerical remedy for an arbitrary new Gaussian error.

An explicit live retry can also recognize a purged timeout allocation without
`execution.json` when the saved scheduler ownership/timestamps and the exact
job's TIME LIMIT diagnostic prove termination and a successful queue query
confirms inactivity. Missing files alone are never sufficient. This permits
operator-budgeted retries of deliberately bounded probes without turning on
unlimited automatic retries.

## Numerical-gradient stagnation recovery

Molpro/Sella tasks that exhaust their step budget can receive one automatic
accuracy retry when the final twenty optimizer steps all have trust radii at
or below 0.001, forces above the requested target, and no material improvement
in the best force. Ordinary slow progress, native execution errors, analytic
gradients and deliberately bounded fragment probes do not trigger this rule.
Already exhausted accuracy settings are not repeated.

The graph owner verifies the attempt and trajectory, archives the failed task,
and resumes its last evaluated positions with a fresh optimizer and recomputed
forces. It tightens SCF/CC convergence to at least 1e-11/1e-10 and, for MR
models, CASSCF energy/gradient to 1e-11/1e-7 and correlation energy/coefficient
to 1e-11/1e-8. Tighter existing values stay tighter. It selects Molpro's
four-point numerical gradient without changing displacement size, selected
energy, active space, roots, state weights, reference policy, shift or IPEA.
The stationary-point order, original constraint targets, force criterion and
step budget remain unchanged. Completed or allocated descendants prevent
automatic source replacement; unrelated tasks retain their inputs and jobs.
`state.json` records the proof in `molpro_gradient_retries`. A second failure
remains available for review; the history survives controller refreshes.

The same controls are available for an explicitly reviewed task retry:
`"calculator_kwargs": {"numerical_gradient": {"stencil": "fourpoint"}}`.
An optional positive `step_bohr` sets both distance and symmetric-displacement
steps; it is not automatically reduced by this recovery. See the
[Molpro FORCE documentation](https://www.molpro.net/manual/doku.php?id=energy_gradients)
for four-point accuracy and cost. Local tests verify recovery and input
generation, not convergence of a future native calculation.
