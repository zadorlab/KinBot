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
tables are not treated as solver iterations. Shared exclusions appear once;
Slurm exclusions are grouped by their actual scheduler value. Paths within a
PES are relative to that PES and process rows show module names/PIDs. Long
diagnostics wrap to the terminal width without truncation.

Both PES overviews appear before either expanded section. Routine successful
result records and completed reference workers are counted together; dependent
waiting tasks are summarized by their blockers. A waiting task is not a failed
QC calculation. The normal display shows at most eight detailed records and
four controller/four QC tails per PES. `--all` shows every reaction row and up
to six controller/six QC tails. `--details` expands every task/record and full process commands,
work directories, duplicate status records and recovery metadata. `--raw-logs`
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
  as a current failure. This does not turn a rejected reaction into a success.
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
with `--all`). It prioritizes controller logs and errors. Use `--all --details` for all
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
