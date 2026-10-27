# Durable restart implementation status

## Accepted ASE geometry progress

Dispatcher ASE optimizations publish immutable evaluated points below
`tasks/<id>/.kinbot_checkpoints/optimization/<contract>/`. The pointer is
published only after the point has been flushed. The contract binds the complete
task record, original starting coordinates, atom order, and the original step
budget. A fresh execution of that identical task automatically restores the
last published accepted coordinates. In Sella, the accepted trajectory index
comes from the optimizer log; the last trajectory frame can be a Hessian probe.

The working trajectory, log, step outputs and referenced wavefunction files are
archived before counters restart. Slurm's retry wrapper excludes the durable
checkpoint directory from its rollback snapshot. Restaging an identical task
through a recovery transaction copies the last hash-bound point from the old
attempt; the old attempt remains intact. Changed task definitions do not inherit
the checkpoint through this path.

This is geometry continuation. It recomputes the force and optimizer Hessian;
it does not claim reuse of interrupted native electronic iterations or frequency
displacements. The original total step budget survives process restarts. Existing
geometry, stationary-point and electronic-state validation still runs.

## Persistent allocation retries

Shared `queue_live_policy` options, also available through `kinbot.run_policy set`:

```json
{
  "max_jobs": 3,
  "allocation_retries": 2,
  "allocation_retry_delay": 300,
  "allocation_retry_max_delay": 3600
}
```

`allocation_retries` defaults to zero. Enabled dispatcher recovery requires a
matching Slurm job ID and working directory with `NODE_FAIL`, `BOOT_FAIL` or
`PREEMPTED`, plus an independent healthy query proving the allocation is not
live. It verifies the staged inputs and any execution record before preparing a
retry. Completed results are reconciled first. An active or complete descendant
prevents resetting its parent.

The `retry_wait` entry saves the evidence and absolute retry deadline. Delays
double between attempts up to the configured maximum. A surviving driver waits;
a fresh driver resumes the same wait. Slurm communication failures do not retire
allocations or reset the retry count. If Slurm requeues the original ID, KinBot
reattaches to it. Once due, the ordinary transactional recovery archives the old
attempt and restages the same task; admission outages leave that one replacement
staged. Exhausting the allocation budget is a hard infrastructure failure, not
scientific unavailability of an optional component.

Unknown disappearance, `CANCELLED`, generic `FAILED`, out-of-memory and walltime
are not reclassified as node failures. Walltime has its existing separate policy.
Unexplained/admin cancellations still require attributable evidence and an
explicit recovery decision. Native restart tasks and deliberately bounded VRC
diagnostics retain their existing separate restrictions.

This implementation applies to dispatcher graphs. Ordinary PES QC and orbital
recon controllers have separate existing recovery paths; this document does not
claim that they already share this new persistent allocation procedure.

## Read-only seed evidence

`python -m kinbot.run_seed_audit SOURCE --output-parent OUTSIDE_SOURCE` exports
shallow PES/search metadata, full `kinbot.db` files and SQLite sidecars, current
graph/task/execution records, and orbital request inventories. `--task RELATIVE`
adds the selected task's native text, optimizer trajectories and archived
attempts. `--evidence PRESERVATION_DIRECTORY` adds outage/preservation receipts
and their copied native text evidence. Binary checkpoints are inventoried in the
receipts, not recopied into this audit.

Every captured file is rechecked before completion. Missing, changed, oversized
or symlinked files produce an explicit incomplete status. Default limits are
2 GiB total and 256 MiB per file. The source is never opened for writing and no
driver, QC executable, submission or cancellation is run. This is a sequential
audit of a quiescent source, not a live filesystem snapshot or a full run backup.

## Periodic native checkpoint candidates

The shared local-scratch workspace can periodically publish named native files
to `.kinbot_checkpoints/native/<label>/<attempt>/` in the result directory.
This applies to callers of `qc_scratch.workspace`; it does not change native
input theory, geometry, state definitions, or restart keywords.

```json
{
  "qc_checkpoint_interval": 600,
  "qc_checkpoint_max_mb": 16384,
  "qc_checkpoint_generations": 2
}
```

The interval defaults to zero (disabled); `_r8`'s staged overrides enable a
ten-minute interval, a 16 GiB ceiling per candidate set, and two retained
generations per native invocation. The copy worker polls every 30 seconds.
The ceiling is not a total run storage quota: each force evaluation or other
native invocation has its own attempt directory. Older attempts remain intact.

The publisher copies recognized `.wfu`, `.chk`, `.rwf`, `.gbw` and named MRCC/
CFOUR restart files from the invocation's private scratch directory (at most two
directory levels below it). It records exact staged-input and task hashes,
Slurm identity, file hashes, and bounded native output tails. It rejects symlinks,
nonregular files, oversized sets, and files that change or grow during copying.
The entire recognized set must remain stable across the copy. File bytes and
directory entries are flushed before publishing the generation pointer. Readers
must verify that pointer and every file hash before considering native reuse.

Publication/storage failures retain the prior complete generation and retry
with backoff up to an hour. Interrupted pointer publication can leave one extra
uncommitted generation; subsequent attempts prune old orphans while preserving
the pointed-to generation. A failed final copy retains local scratch. The live
output mirror also retries after transient filesystem failures. Normal
exit joins the copy thread before cleanup; partial copies never replace a good
pointer. The Slurm rollback wrapper excludes `.kinbot_checkpoints`.

These files are labelled `restart_candidate_requires_native_validation`.
Stable bytes are not proof of a consistent native restart phase, especially
when several native files must agree. A hard kill before the first successful
publication, or a shared filesystem unavailable for the entire invocation,
can still leave no durable candidate. Backend-specific restart validation and
automatic selection/restoration are separate work and are not claimed here.

## Validation and work still required

Local tests exercise actual ASE BFGS and Sella with synthetic force calculators,
abrupt interruption, preservation of accepted points and step budgets, archive
transaction interruption, fresh drivers, scheduler/requeue/admission races,
completed/live sibling preservation and repeated invocation. They also exercise
SQLite WAL and failed-search-history capture without changing the source files.
The periodic publisher tests include a real subprocess hard kill before normal
workspace finalization followed by removal of its local scratch, changing files,
interrupted pointer writes, copy ceilings, and shared-storage outage/recovery.
These are local simulated-site checks, not a verified native QC restart.

Still required for the complete migration: a general independent seed importer
with compatibility and path-rebinding proofs; retention of terminal PES failures;
reassessment of accepted paths under a changed L1 policy; backend-specific native
restart validation and restoration; consistent retry integration into
the ordinary QC and orbital owners; and remote native verification. The new run
must not be labelled ready based only on the features documented above.
