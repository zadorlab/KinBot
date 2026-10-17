# Read-only PES monitor and log tracker

Use either one PES directory or a production handoff parent containing
`prepared.json`:

```bash
python -m kinbot.run_monitor /path/to/run
python -m kinbot.run_monitor /path/to/run --all
python -m kinbot.run_monitor /path/to/run --follow
python -m kinbot.run_monitor /path/to/run --json
```

`examples/anl/tools/run_examples.py monitor RUN` accepts the same
options. `watch -n 10 'python -m kinbot.run_monitor RUN'` refreshes a snapshot.
`--follow` shows a snapshot then follows newly appended controller and native
QC text, discovering new logs and handling rotation/truncation. Ctrl-C stops
only the viewer. `--lines N` changes recent progress lines; `--interval N`
sets the follower poll interval (1–60 seconds; default 3).

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
with `--all`). It prioritizes controller logs and errors. Use `--all` for all
reaction/task rows, including accepted/completed entries. It skips archives,
scratch, checkpoint payloads and per-sample RotdPy trees, and does not follow
external symlinks. Log warnings may describe earlier attempts; quiet output
alone does not prove a stall. The monitor never advances a dispatcher,
modifies an ASE database, submits jobs, or retries calculations.
