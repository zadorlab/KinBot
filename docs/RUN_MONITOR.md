# Read-only PES monitor and log tracker

Use either one PES directory or a production handoff parent containing
`prepared.json`:

```bash
python -m kinbot.run_monitor /path/to/run
python -m kinbot.run_monitor /path/to/run --all
python -m kinbot.run_monitor /path/to/run --all --details
python -m kinbot.run_monitor /path/to/run --all --raw-logs --lines 12
python -m kinbot.run_monitor /path/to/run --follow
python -m kinbot.run_monitor /path/to/run --json
```

`examples/anl/tools/run_examples.py monitor RUN` accepts the same
options. `watch -n 10 'python -m kinbot.run_monitor RUN'` refreshes a snapshot.
`--follow` shows a snapshot then follows newly appended controller and native
QC text, discovering new logs and handling rotation/truncation. Ctrl-C stops
only the viewer. `--lines N` changes recent progress lines; `--interval N`
sets the follower poll interval (1–60 seconds; default 3).

The dashboard separates **KINBOT LOGS** (including each `kinbot.log` and
`pes.log` tail) from **QC OUTPUT** (native solver and Sella progress/errors).
Input-parameter dictionaries are excluded from the progress view. Coordinate
tables and Molpro dataset record numbers are not treated as solver iterations.
Shared exclusions appear once;
Slurm exclusions are grouped by their actual scheduler value. Paths within a
PES are shortened to ANL state/task or VRC channel/task labels; adopted jobs
with an old run alias are resolved for display. Process rows show module names/PIDs. Long
diagnostics wrap to the terminal width without truncation.

Both PES overviews appear before either expanded section. Routine successful
result records are omitted from the default view. Dispatch task counts are
aggregated across graphs; current failures, pending retries and dependencies
blocked by failure are expanded. Completed reaction rejections are counted
as **not-retained searches**, with individual rows under `--all` or `--details`.
This is the driver's decision about that search, not proof that the physical
pathway cannot exist. Unsuccessful TS guesses, incorrect IRC connections and
discarded conformers can be ordinary exploration outcomes.

Nonzero QC exits and unsuccessful `.result.json`/`.scf_recovery.json` records
are counted separately from task/controller errors. Where a later reaction
monitor identifies the owning search, counts distinguish not-retained,
accepted and still-active searches. Shared fragments, unmatched jobs and QC
evidence newer than the search decision remain unclassified. These are counts
of attempt records, not counts of broken channels or recommendations to retry.
Only a same-worker, exact reaction-name prefix with an underscore boundary is
matched; records from different wells are never combined by job name alone.

Default QC tails omit older attempts under terminal searches. `--all`,
`--details`, `--raw-logs` and JSON retain their diagnostics within the log read
limits. Explicit input-syntax, missing-module and I/O diagnostics remain
visible even under a rejected search, as do controller/task errors. Those
diagnostics identify technical problems; they do not certify that a native
retry will succeed. Scientific rejection alone is not a repair request.
A waiting task is not a failed QC calculation. Completed reference-worker
log tails are available under `--details`, keeping the main `kinbot.log` and
`pes.log` visible by default.

The normal display shows four controller/four QC tails per PES. It reads one
tail per dispatch task before older force steps, and shows one representative
healthy log per task; current error logs retain priority. `--all` shows every reaction row and up
to six controller/six QC tails. `--details` expands every task/record and full process commands,
work directories, duplicate status records and recovery metadata. These views
cover current records, not archived attempts. `--raw-logs`
shows the unfiltered last `--lines` lines of each log instead of selected
progress. Use it to inspect output not recognized by the progress filter.
`--json` retains original paths, commands, status records and bounded raw tails.
`--follow` continues to stream newly appended raw text, including new logs.

A monitor-only code update can be pulled while calculations run. An existing
`watch` invocation loads the updated viewer on its next refresh. A persistent
`--follow` viewer needs to be stopped with Ctrl-C and started again; that stops
only the viewer. No restart of the PES controllers or QC jobs is needed.

The snapshot shows:

- Live Slurm job IDs, states, hosts, elapsed/wall time, CPU/memory requests,
  work directories, pending reasons and actual scheduler exclusions.
- Local controller/worker commands and PIDs, explicit completion/failure
  records, PES worker counts, and nonzero QC exit markers.
- Actual recent `pes.log`/`kinbot.log` events and native Gaussian/Molpro
  progress/errors with file paths and ages; driver license banners are omitted.
- Signed frequency rejection reasons, latest Sella progress, minimum-recovery
  attempts/selected displacement and final `.result.json` status. Archived
  rejected Hessians are skipped so they do not masquerade as new failures.
  Failed QC exit markers prioritize that job's native/Sella logs ahead of
  newer routine output; old worker frequency-rejection messages are recognized.
  ANL `minimum_recovery.json` reports show curvature checks, retry controls,
  selected descent energies and explicit failure reasons before dependent
  geometry tasks are released.
- After the explicit saved-results migration, hash-matched rejected QC exit
  markers are shown as cached outcomes, with their original native attempts
  in the external archive. A marker written by a new failed attempt is shown
  as a new attempt outcome. The exact matching `exited with worker status` message
  is also collapsed in the default controller tail, but remains in `--details`,
  raw tails and JSON. Other errors are not hidden based on age or job name alone.
  The same message is collapsed when a later owning search is terminal.
  This does not turn a rejected reaction into a success.
- Dispatcher `attempts/` folders beside `workflow.json` are historical and are
  excluded from the dashboard, JSON and log follower. Their paths are listed in
  the JSON snapshot; all archived files stay untouched on disk. A retry's current
  status and any new failure under `tasks/` remain visible.
- Reaction-stage counts, pending optimizations and failures, ANL/VRC dispatch
  task counts and dependency blockers, post-reaction/MESS execution gates.
- RotdPy execution records and surface-output file counts. File counts do not
  measure accepted Monte Carlo convergence or rate accuracy.

No `sacct` database is needed. Slurm outages and unreadable records are shown
explicitly. Absence of jobs/processes is not success; old `pes_summary.txt`
“Finished” entries only report that workers exited. New runs additionally
record `kinbot_status.json` and `pes_status.json`; a completed controller is
still distinct from scientific acceptance of its outputs.

The view reads bounded tails (64 KiB/file, normally up to 40 logs, up to 200
with `--all` or `--details`). It prioritizes main controller logs, failed tasks
and live allocations. Log coverage is reported in JSON and `--details`; a
bounded read is not a read warning. Use `--all --details` for all
reaction/task rows, including accepted/completed entries. Save this view to a
file instead of trying to capture a long `watch` screen:

```bash
python -m kinbot.run_monitor /path/to/run --all --details > monitor-report.txt
```

It skips archives,
scratch, checkpoint payloads and per-sample RotdPy trees, and does not follow
external symlinks. Log warnings may describe earlier attempts; quiet output
alone does not prove a stall. The monitor never advances a dispatcher,
modifies an ASE database, submits jobs, or retries calculations.

For detailed native CAS/MRCI orbital evidence, workers separately write
`*.orbital_audit.json` alongside their outputs. The heavier recursive
`python -m kinbot.molpro_diagnostics RUN` command can rebuild/list those reports
for both PES trees without changing calculations. It is deliberately separate
from this bounded live dashboard; see the [native evidence contract](MOLPRO_NATIVE_ORBITAL_DIAGNOSTICS.md).
# Gaussian SCF recovery

`<job>.scf_recovery.json` records a bounded Link 508 recovery: original and retry
SCF options, preserved native artifacts and outcome (`retrying`, `complete`, or
`failed`). These records appear with the other workflow artifacts. Logs under
`restart_gaussian_scf/` are historical and do not count as current failures.
See the [SCF recovery record](GAUSSIAN_SCF_RECOVERY_2026-10-16.md).
