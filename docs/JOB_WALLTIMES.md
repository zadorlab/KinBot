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

## Automatic walltime retries

Slurm timeout recovery is opt-in and independent of the existing I/O,
scratch-node, L3 geometry, and VPT2 symmetry retry budgets. It covers ordinary
KinBot L1/L2 workers (including Gaussian geometry, frequency and single-point
workers) and all external/ASE ANL and VRC dispatcher tasks: Gaussian, Molpro,
MRCC and CFOUR. It does not implement an allocation retry for MESS or ROTD
sampling, or for PBS/local execution.

Enable one replacement allocation per logical calculation, with a seven-day
minimum and a factor of two applied to the previous allocation's actual limit:

```bash
python -m kinbot.run_policy set /path/to/run \
  --minimum-walltime 7-00:00:00 --timeout-retries 1 --timeout-walltime-factor 2
python -u -m kinbot.run_refresh /path/to/run --apply
```

Refresh once to install this capability in older pinned controllers. Running
and pending allocations are retained; their current limits do not change.
Later policy edits are live. Setting `--timeout-retries 0` disables creating
new timeout retries. Prepared transactions still finish consistently; pausing
admission with `--max-jobs 0` also pauses their submissions.

The same controls belong in the single KinBot JSON input, for example:

```json
"queue_live_policy": {
  "max_jobs": 8,
  "excluded_nodes": [],
  "minimum_walltime": "7-00:00:00",
  "timeout_retries": 1,
  "timeout_walltime_factor": 2.0
}
```

The node count here is illustrative, not a program default. Existing policy
capacity and exclusions are unchanged by a timeout-only live edit. Defaults
remain zero timeout retries and factor two; allowed retries are 0–10. A retry
uses the same partition, capped at its current maximum. A 24-hour failure
therefore receives seven days with these settings; a seven-day failure receives
14 days if the partition permits it. No longer allocation available means a
visible hard failure, not another identical attempt or a lower scientific tier.

Recovery requires the exact job to be absent from the live queue and confirmed
`TIMEOUT` in its matching `scontrol` record. When that record has expired, an
exact job-ID Slurm `DUE TO TIME LIMIT` diagnostic plus the saved allocation,
work directory, start time and limit can establish exhaustion without `sacct`.
Unknown scheduler status, cancellation, arbitrary exit codes, and generic
nonconvergence are not timeout evidence. Complete verified results are accepted
before retry; missing result publication receives the existing grace period.

**These automatic replacements are clean retries, not continuation of the
interrupted optimizer, Hessian probe sequence or coupled-cluster iteration.**
Only the failed task is repeated. Previously completed graph tasks and VRC scan
points remain available. Dispatcher retries restage the same task using its
accepted upstream geometry, preserve input/provenance checks, and block dependent
work until completion. Ordinary QC workers restore hash-verified pre-allocation
seeds from the existing I/O guard, preserving required `Guess=Read` inputs while
rejecting checkpoints changed during the failed allocation. The worker's
scientific settings and random seed are retained. Legacy/custom workers without
those seed records require review rather than guessed recovery.

Failed attempts are archived under the graph's `attempts/timeout_*` or the
ordinary worker's `.kinbot_timeout_retries/`. Journals survive controller
replacement, and profiled L1/L2 route ledgers adopt replacement IDs. Shared
admission, exclusivity, exclusions and local scratch policies still apply.
Uncertain `sbatch` submissions are never blindly repeated. Timeout exhaustion
does not consume or reactivate the independent L3 geometry retry budget.

## Existing continuation and other recovery

- Slurm I/O guard: by default two retries **inside the current allocation**;
  no extra walltime. It restores saved pre-allocation seeds.
- Scratch startup: by default two allocation retries for verified local
  scratch faults, with failed-node exclusions.
- ANL L3: configurable `l3_max_retries` (normally one) gives a clean retry of
  execution/scheduler failures. It does not resume Sella curvature probes.
- VPT2: bounded recovery for the recognized framework-symmetry diagnostic.
- Direct MRCC: `kinbot.anl.dispatch resume-mrcc GRAPH TASK` is an explicit
  amplitude restart using `rest=1`, not an automatic timeout policy. The
  checkpoint is now hash-bound and staged into node-local execution; `fort.17`
  is also staged if present. The [MRCC manual](https://www.mrcc.hu/MRCC/manual/pdf/manual.pdf)
  requires `fort.16`, and additionally `fort.17` for multiple roots or geometry
  optimization. Missing required auxiliary files need review before using this
  manual route. A nonempty/hash-matching file alone does not prove its native
  consistency after an interrupted write.

New MRCC workers preserve amplitude files once after a handled native failure
for manual review. They do not repeatedly copy amplitudes to shared storage
while running or add them to every successful job's outputs. An abrupt Slurm
kill can prevent publication; no automatic partial-checkpoint salvage from old
node-local scratch is promised. Gaussian `Opt=Restart`/`Freq=Restart`, Molpro
WFU continuation and CFOUR checkpoint continuation are **not** enabled by this
change. Their timeout replacements use the clean-retry route described above.

Validation uses scheduler mocks, interruption tests and native-program shims;
it does not substitute for a real Gaussian/Molpro/MRCC/CFOUR timeout on the site.

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
