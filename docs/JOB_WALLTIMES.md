# Job walltimes and read-only activity checks

ANL MRCC CCSDT(Q)/CCSDTQ(P) tasks now default to seven days, including the
previously 12-hour DZ and 24-hour TZ quadruples jobs. TZ L3 optimization and
numerical harmonic tasks also default to seven days; ANL1 QZ geometries,
5Z/6Z reference energies, and basis-matched harmonic tasks already did so.
`vrc_tst_walltime` now defaults to `7-00:00:00` for Molpro VRC calculations.
Explicit task/input overrides remain supported. Automatic ANL jobs respect an
explicit partition's maximum time through the existing resource resolution.
These defaults affect newly built graphs, not persisted graphs or allocations.

For an existing shared pool, `queue_live_policy.minimum_walltime` (or its live
CLI counterpart) can extend *future requests* from already-staged scripts:

```bash
python -m kinbot.run_policy set /path/to/run --minimum-walltime 7-00:00:00
```

The pool's current capacity and exclusions are preserved. Existing pinned
controllers must first be refreshed onto code supporting this option:

```bash
python -u -m kinbot.run_refresh /path/to/run --apply
```

This adopts running and pending jobs and archives/replaces only controllers.
Do not use `--replace-pending` for a walltime-only change. No native job is
cancelled or requeued. Scratch selection, exclusivity, cores, memory,
convergence, active spaces, and accepted results are unchanged.

## Extend existing allocations

```bash
python -m kinbot.run_walltime /path/to/run
python -m kinbot.run_walltime /path/to/run --minimum 7-00:00:00 --apply
```

The first command is read-only. The second selects the current account's jobs
whose resolved work directories lie under the run, checks partition maxima,
rechecks current limits immediately before the request, and issues a single
`scontrol update JobId=... TimeLimit=7-00:00:00`. Longer limits are retained.
Array aggregates and unknown work directories fail closed. A scheduler timeout
is followed by verification, never an automatic repeat of the update.

The command prints verified increases and unresolved IDs separately. A nonzero
exit indicates unverified jobs or partition-limit skips; no restart is needed
to recover from a rejected update. If a job leaves the queue during the check,
its extension cannot be verified. Coordinate with other operators: Slurm has
no compare-and-set time-limit update. Do not concurrently lower/change the same
limits from another session.

Slurm ordinarily allows only a privileged user to increase a running or
suspended job's limit. On `Access/permission denied`, send the printed command
to the cluster administrator; do not cancel/requeue the computation. Pending
jobs may be updated even when running-job requests fail, so the command checks
each job after the request. See [Slurm scontrol](https://slurm.schedmd.com/scontrol.html).

## Read-only progress

```bash
KB_REPORT=$(mktemp -d /tmp/kinbot-timecheck.XXXXXX)
python -m kinbot.run_walltime /path/to/run --progress --output "$KB_REPORT/report.json"
```

This samples Slurm CPU/I/O counters and bounded native-log content twice, 20
seconds apart. It prints remaining walltime and optimizer, curvature-probe, or
native iteration evidence. The only written file is the explicitly requested
new report outside the run. Native logs, inputs, task state and scratch are
untouched. Permission/malformed-read errors are reported, not treated as success.

The scratch mirror can refresh a `.live.log` mtime without adding calculation
output. This check compares content hashes instead. CPU/I/O or output activity
supports “still working,” not a finish-time guarantee. An unchanged 20-second
sample does not establish a stall. Curvature probes are extra work after force
convergence; iterative coupled cluster can have a substantial noniterative
correction afterward. A residual trend or optimizer force is not a whole-job
ETA, and neither estimates the full remaining VRC/ANL network.
