# Recovery after the interrupted r6 handoff

The user reported that the earlier restart command hung while Blodgett was
having problems. Its actual stopping point is unknown. The approved action is
to **archive r6 and start fresh in r7**, not reuse r6 results. The r7 inputs
retain the exact scientific settings committed in `29432951`.

## Command

Run on the same Blodgett login host used for the previous launch:

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
printf 'Updating KinBot (120-second network timeout)...\n'
timeout 120 git pull --ff-only origin composite
bash examples/anl/tools/recover_on_blodgett.sh
BASH
```

The command prints a detached handoff PID and a unique log such as
`/home/hram/kinbot-runs/r7-restart.XXXXXX.log`. Follow the printed `tail -F`
command; Ctrl-C stops only that log viewer. The detached helper and the eventual
PES drivers do not depend on the SSH connection staying open. A PID alone is
not confirmation that the PES has launched.

The fresh parent is
`/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261017_r7_dft_scf8`.
It has `restart_status.json` with the current phase/error and
`restart_report.json` with preparation, archives and driver information.

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python -m kinbot.run_monitor /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261017_r7_dft_scf8'
```

## Interruption handling

- A nonblocking handoff lock prevents two recovery helpers from launching
  duplicates. Preparation now hashes only declared inputs/checkpoints; runtime
  reference JSON/WFU files no longer contaminate its immutable manifest.
- Backend checks precede retirement. An unavailable scheduler is an error,
  never evidence of an empty queue. The helper records the failing phase and
  exits; rerun the same command after the cluster recovers.
- If r5/r6 already disappeared during the prior archive step, the helper records
  that absence and leaves existing `.stopped.*` archives intact. It does not
  require an old directory to exist merely to launch a new one.
- If a package preparation/launch helper for r6 is still alive outside the run
  directory, its exact script/action/root arguments identify it. Its PID and
  process start time enter the existing scoped stop/descendant tracking. No
  user-wide process kill, PID-file-only kill, or shared-scratch deletion occurs.
- Existing r5/r6 PES trees are stopped, checked inactive, and moved as complete
  dated sibling archives. Retirement interrupted between the two PES can be
  repeated. A partially prepared old tree with valid case markers can also be
  archived. Unrecognized directories are rejected rather than guessed about.
- If only one r7 driver started, a repeat detects its actual process and starts
  only the unstarted case. If both run, both are retained. A valid completed PES
  is also retained. A previously started r7 case with no live driver is reported
  for diagnosis, not blindly relaunched over scientific results or orphaned QC.
- The helper checks recorded driver hosts before relying on local `/proc`.
  It cannot terminate controllers on a different login host, revive a down
  filesystem, or force Slurm to finish cancellation. Those conditions must
  clear before the handoff can safely finish.

Resources, exclusions, thresholds, electronic states, seeds, and all chemical
settings are unchanged. The monitor still covers both PES. None of these local
tests establishes native production completion.

## Verification

**85 focused tests passed**, including the existing cleanup, monitor,
saved-result repair and production-handoff regressions. Shell syntax and diff
whitespace checks also passed. The chemical input files were not changed.

The recovery tests cover a scheduler outage, interrupted retirement, interruption
between driver launches, repeated launch, preservation of native evidence,
completed/stale new runs, handoff locking, exact launcher selection and
partially prepared old trees. Tests run with simulated cluster/process state;
no Blodgett jobs are controlled locally. Commands and results are recorded in
`examples/anl/tools/validation/handoff_recovery.log`.
