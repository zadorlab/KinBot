# Node-local QC storage, exclusive nodes, and outage recovery

This replaces the HOME scratch recommendations in earlier handoffs. On October
19 the supplied login session could not read its home directory or shell
profiles (`Input/output error`), while Slurm responded. That establishes a
filesystem-access failure. It does not establish that these calculations caused
the storage outage; that requires the administrators' storage/kernel evidence.

## Source of the regression

Commit `05fd09dd950a59893f337db75537735e0d86dace` added HOME and result-directory
fallbacks when local scratch lacked capacity. Commit
`c340cc8961ad8f3d2caf25c4425537530fee9a9b` explicitly exported the HOME cache
in the production handoff. Later restart/resume scripts inherited it. This
violated the requested node-local policy. Capacity alone cannot establish
locality, and a larger shared filesystem is not an acceptable fallback.

Molpro also has three separate locations to control: scratch (`-d`), its
integral repository (`-I`), and wavefunction repository (`-W`). Its documented
default for the last is `$HOME/wfu`. Setting only `-d` leaves that default in
effect. Command-line choices outrank inherited Molpro options. See the
[Molpro execution documentation](https://www.molpro.net/manual/doku.php?id=running_molpro).

## Reconciled requirements and implementation

| Requirement | Behavior |
|---|---|
| General, automatic account selection | Resolve the executing UID through the operating system; try `/scratch/<account>` first on the compute node. No username is embedded in the implementation. |
| No inherited HOME override | HOME and paths resolving into HOME are rejected. Old Molpro overrides cannot outrank the preferred account scratch. |
| Verified alternatives | Check the containing filesystem/mount, free capacity, private-directory creation, write/fsync/rename/read-back. Reject network, unknown, and memory-backed filesystem types. If no candidate works, stop before native QC. |
| All native work in the active example workflow | Gaussian, Molpro force jobs, ANL Molpro/Gaussian/CFOUR/MRCC tasks, orbital reconnaissance, legacy Molpro and ROTD sampling use a private local working directory. ORCA and Q-Chem adapters also use it. |
| Wavefunction and integral placement | Pass explicit local `-d`, `-I`, and `-W` for Molpro; set local `GAUSS_SCRDIR`/`TMPDIR` for Gaussian. Stage the explicitly referenced input checkpoints. Absolute native checkpoint/output paths cannot bypass this policy. |
| Results and diagnostics | Publish native outputs, required results and restart checkpoints back to the calculation directory. Do not publish integral/intermediate scratch. A bounded live tail (up to three 16-KiB tails every 30 seconds) supports the monitor. |
| Copy-back failure | Retain local files, include host and path in the error, and never mark the calculation successful. The runtime JSON records the local path before launch. |
| Exclusive QC nodes | Ordinary Slurm L1/L2/reference/pivot jobs default to `queue_exclusive=1`, including custom headers. ANL, reconnaissance and legacy Molpro jobs request exclusivity; ROTD defaults to `rotdpy_exclusive=1`. Conflicting sharing directives are rejected. |
| Resource separation | Reserving a node does not increase native ranks/threads, memory requests, or the two examples' separate ten-allocation limits. Drivers stay on the login host. |
| Preserve scientific settings | No changes to methods, selected active spaces, orbital indices, spin/root handling, CBS exponent, IPEA/level shifts, convergence targets, VPT2 symmetry recovery, or L3 geometry dependencies. |
| Salvage conservatively | The read-only recap distinguishes active ownership, verified outputs, recorded successes, failed attempts, missing outcomes and invalid geometry. It neither submits nor cancels jobs. |

The single KinBot JSON accepts `qc_scratch_fallbacks` (default `[]`, extra
absolute local paths) and `qc_scratch_min_mb` (default 4096; task-specific larger
requirements still win). These settings propagate to nested workers. Candidate
selection checks `/scratch/<account>`, then explicit candidates and alternatives,
site/scheduler candidates and `/tmp`. Every candidate must pass the same checks.
Previous explicit I/O failures may exclude that root for a recorded retry.
These workspaces require one node per native allocation.

The kernel mount check cannot prove how a site's block device is physically
backed. Administrators should confirm the intended local mount. Here the user
identified `/scratch/<account>` as actual node-local disk. Private scratch
survives failed publication only while the node/site retains it; node loss or a
site epilog can still remove it.

Slurm `--exclusive` requests a node without other jobs, subject to the
partition's sharing policy. It allocates all CPUs but only the requested memory;
`--mem=0` is a separate request for all memory. This distinction matters when
interpreting `squeue`: 96 allocated CPUs does not mean 96 Molpro MPI ranks.
See [Slurm's sbatch documentation](https://slurm.schedmd.com/sbatch.html).

Shared run storage still holds controller state, scheduler logs, published
results and restart checkpoints. It is not used as computational scratch.
An unavailable home directory still prevents loading the installed Python/code
or updating controller state. This change cannot make an ongoing filesystem
outage safe to ignore.

## Read-only recovery report after storage is available

Do not run cleanup/restart against an inaccessible filesystem. Once normal
reads work, pull the fix and inventory the existing run:

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
git pull --ff-only origin composite
KB_ROOT=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage
KB_REPORT=$(mktemp /tmp/kinbot-r6-recap.XXXXXX)
.venv/bin/python -m kinbot.run_recap "$KB_ROOT" --json "$KB_REPORT" --limit 6
BASH
```

The compact report includes current failure sentences. The JSON preserves all
task and branch details. A missing scheduler job is never proof of completion.
`recorded_success` means a completion record exists; it is not yet a salvage
verdict. Once storage is stable, add `--verify` to validate native results,
artifact hashes and geometry provenance. This reads the declared artifacts,
including potentially large checkpoints, so it is deliberately not the first
step during a storage incident. Ordinary L1/L2 success markers/model files are
reported separately and do not claim composite-level hash verification.

Pulling does not replace pinned running controllers or already submitted jobs.
Review this report before choosing which jobs to adopt, retry or archive. The
recap changes no run files, does no cross-run result reuse, and preserves
expensive active allocations. Unknown errors must be diagnosed before restart;
no command in this section silently treats them as repaired.

### Exporting the native evidence after the first recap

The October 19 recap found no local controllers, one live allocation, 48 task
execution records reporting success, and unresolved submission/geometry/native
failures. The first compact display grouped failed geometry calculations with
their affected descendants; its 40 CH3S `invalid_geometry` entries therefore
did not mean 40 independent native failures. The recap now shows the primary
failures first, summarizes the geometry cascade separately, and identifies the
remaining live job. A recorded `running` branch with no local controller is
labelled explicitly. This does not establish the absence of processes on other
login hosts, so run the recap on the original controller host.

Use `--evidence` to include the recent native failure JSON records (including
saved rank-local diagnostics), input files, submission scripts, task/execution
records, geometry provenance, and bounded native output beginnings/tails in
the JSON. It excludes large wavefunctions and integral files. The terminal
output remains compact. Tasks with ambiguous submission outcomes remain
unresolved even if an execution record exists; the recap never resets or
resubmits them.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
git pull --ff-only origin composite
KB_ROOT=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage
KB_REPORT=$(mktemp /tmp/kinbot-r6-evidence.XXXXXX)
.venv/bin/python -m kinbot.run_recap "$KB_ROOT" \
  --evidence --json "$KB_REPORT" --limit 8
KB_EXPORT="$HOME/$(basename "$KB_REPORT").json.gz"
gzip -c "$KB_REPORT" > "$KB_EXPORT"
printf '\nEvidence export: %s\n' "$KB_EXPORT"
BASH
```

Transfer that small compressed JSON for diagnosis before preparing a mutation
or restart plan. This export does not hash large checkpoint files (`--verify`
is separate), alter run files, cancel the surviving allocation, or launch jobs.
The report is a point-in-time observation; jobs can finish while it is collected.

A compute-node storage-only probe, when the cluster is healthy, can use the
same resolver without starting Gaussian or Molpro:

```bash
srun --partition=short-cpu --nodes=1 --ntasks=1 --exclusive --time=00:02:00 \
  --exclude=node1,node9,node10,node13,node14,node16,node21,node23,node24,node30,node33 \
  /home/hram/KinBot/.venv/bin/python -m kinbot.qc_scratch --backend molpro --probe
```

It prints the actual node, selected root, filesystem and I/O check result and
removes its private test directory. It is not a native electronic-structure
validation or a guarantee against later disk exhaustion.

## Validation scope

Regression coverage includes runtime selection, inherited overrides, Linux
mount parsing and symlinks, EIO failures, checkpoint staging, copy-back failure
retention, native Molpro force execution through a simulated executable,
reconnaissance, VRC profiles, Slurm headers, recovery ownership and read-only
recaps. Existing L3 salvage, controller refresh, automatic ANL and handoff tests
remain part of the verification. Licensed cluster executables and the restored
cluster filesystem are not available to the local test suite.

Final targeted suite: 630 passed, 10 subtests passed, one Linux-only pidfd test
skipped on macOS. Python 3.11 syntax, shell syntax, and a comparison confirming
unchanged scientific parameters in both example JSONs also passed.


## Reviewed r6 evidence and in-place recovery

The local export `kinbot-r6-evidence.L3YrXs.json.gz` was collected at
2026-10-19 19:15:16 UTC. It contains no read errors and reports no surviving
local controllers. It is a snapshot, not a current scheduler query.

* CH3SS Gaussian VPT2 allocation **52909915**, on node27, was still writing
  second-derivative output. Its task is `state_6a45048aedc75231/gaussian_vpt2`
  for `[S]S=C`. Preserve this allocation and its attempt-2 `Symmetry=(PG=C1)`
  input if it remains live. Its existing process does not acquire new runtime
  policy merely because the repository is updated.
* CH3S `state_074fcc5ef5188035/l3_geometry_tz` completed the electronic and
  numerical-gradient output, then failed during wavefunction/file handling.
  The native launcher reported a missing source WFU in shared task-directory
  scratch. This is not evidence that the CC iterations failed to converge.
* CH3SS `relaxed_scan` failed before MULTI with a Fortran I/O error:
  `File already opened in another unit` in `iow.F90`. The recorded native
  working/scratch paths were long shared task-directory paths. The corrected
  private, short local workspace removes that configuration; the export does
  **not** prove which internal Molpro condition produced this error. Native
  retry is still required to establish that it is resolved.
* The other failed CH3S L3/VRC attempts show MPI aborts. Several central
  electronic evaluations converged before the abort. Their excerpts do not
  establish the underlying cause of every abort; do not relabel them as SCF
  convergence failures or claim the storage outage caused every failure.
* Seven geometry allocations disappeared without final execution records.
  Treat them as interrupted, not complete. Retain accepted TZ geometries as
  QZ seeds; archive incomplete optimization attempts and restart them from
  their accepted upstream geometry. Partial Sella iterations are not promoted
  to accepted geometries.
* H2 `state_3551007628dd3b3d/anl1_5z` is recorded as `submitting`, but has an
  `executed` record and artifact hashes. Verify those files before deciding
  to rerun it. The snapshot has 47 other recorded successes across both cases;
  these counts alone are not salvage approval.

The high invalid-geometry count includes descendants of failed geometry tasks,
not just independently failed optimizations. A single point made on an obsolete
or fallback geometry cannot be reused for the requested L3-based composite.

### Recovery command and guarantees

`python -m kinbot.run_recover RUN` is read-only: it queries the queue, verifies
staged-input and result hashes, parses completed external results, and computes
geometry dependency invalidations. `--apply` archives and prepares the restart;
`--apply --launch` additionally starts the two login-host controllers. There are
no species, job IDs, usernames, or cluster paths in the recovery implementation.
It operates on prepared salvage runs and discovers ANL graphs plus VRC scan,
Molpro pivot, and correction graphs.

Before mutation it requires healthy scheduler queries, no live local
controllers, consistent input provenance, and one unambiguous owner for every
live allocation. It refuses unrecognized live work or live calculations on
incompatible geometries. It never calls `scancel` or signals a process.

A successful result with a lost submission response can be accepted after full
artifact/native validation. Missing results are never synthesized from a queue
exit or a native energy line. Submitted jobs still in the queue keep their task
directories and job IDs; completed results remain in place. Failed, interrupted,
wrong-geometry, and old unsubmitted attempts move to a timestamped sibling
`.outage_archive.*` directory. Controller logs are archived there too. Deferred
single-point retries retain their attempt/archive history and are staged only
when their required geometry has completed. VRC also handles an archived scan
whose replacement has not yet been staged.

Recovery is journaled in `outage_recovery.json`. Repeating the same command can
finish an interrupted archive/input migration or partial two-controller launch.
It will not repeatedly relaunch replacements that have already failed. An
unknown new allocation, changed ownership, or a live controller during migration
stops recovery for review. The full journal records every retained and reset
task; the terminal summary is compact.

Input migration enables exclusive QC allocations and supplies missing local
scratch configuration defaults. It preserves the scientific input settings,
node exclusions, and existing per-example concurrency limits. Replacement
controllers use a fresh immutable code snapshot; already running native jobs
continue with their existing snapshot.

On Blodgett, after storage reads work normally:

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
git pull --ff-only origin composite
KB_ROOT=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage
KB_PY=/home/hram/KinBot/.venv/bin/python

# Storage-only compute-node check; no electronic-structure calculation.
# Every actual job independently repeats locality, capacity, and I/O checks.
srun --partition=short-cpu --nodes=1 --ntasks=1 --exclusive --time=00:02:00 \
  --exclude=node1,node9,node10,node13,node14,node16,node21,node23,node24,node30,node33 \
  "$KB_PY" -m kinbot.qc_scratch --backend molpro --probe --min-mb 63000

"$KB_PY" -u -m kinbot.run_recover "$KB_ROOT" --apply --launch
BASH
```

Keep the old r5 compatibility alias, run directories, checkpoint repositories,
and saved code snapshots until the recovery verifies their dependencies. No
blanket HOME cleanup is part of this command. The archive move is on the same
filesystem; it does not copy native scratch or delete unreviewed files.

Use the same monitor root afterward:

```bash
watch -n 10 /home/hram/KinBot/.venv/bin/python -m kinbot.run_monitor \
  /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage --lines 3
```

The storage and archive tests cannot run licensed Gaussian/Molpro or verify the
restored cluster filesystem locally. Retrying the affected native calculations
is the remaining check for the unresolved MPI and Fortran I/O failures.

Recovery validation: 258 tests passed; one Linux-only process-handle test was
skipped on macOS. Coverage includes lost submission responses, live-job
ownership, invalid geometry descendants, artifact/native rejection, deferred
attempt history, VPT2 input preservation, VRC restaging, interrupted archives
and input migration, and partial two-controller launches. Python 3.11 syntax
and every documented Bash block also passed. The test environment emitted a
non-fatal Sella/JAX cache-permission warning; it did not affect these checks.

## Read-only salvage and scheduling audit while drivers run

`kinbot.run_audit` lists every materialized ANL/VRC task, its molecular identity,
recovery provenance, current recorded status, current Slurm allocation, required
geometry, dependencies, conditional fallback triggers and resources. The complete
list of jobs salvaged at recovery comes from `outage_recovery.json`; these are
results verified at recovery, not a new hash/chemistry validation. It separately
identifies subsequent recorded completions and retained pre-update allocations.
It does not read large wavefunctions or native output files.

The policy report distinguishes:

* User/worker input values, branch slot reservations, and policy source files in
  the recorded pinned controller snapshot compared with the audit's code.
* On-disk `job.slurm` directives and their committed hashes.
* Current `scontrol` evidence for full-node exclusivity and exclusions, including
  Slurm's compressed hostlists. A missing record is unknown, not a pass.
* Native `*.runtime.json` evidence: the compute host, selected scratch directory,
  filesystem, successful I/O probe, available capacity and local working path.
  No scratch record means unobserved. Legacy incomplete records are unknown;
  an alternate verified local root is explicitly distinguished from the preferred
  `/scratch/<account>`. No login-node filesystem probe substitutes for compute-node
  evidence. Retained pre-update jobs keep their original runtime.

There is no fixed global submission order. ANL rotates among graph keys and
selects staged tasks; its in-memory cursor is not persisted. VRC is concurrent,
and its later scan points/corrections, ROTD sampling and MESS work may not yet
have task graphs. The report lists dependency levels and eligibility, not an
invented future queue or completion time. Slurm controls actual start order.

The auditor calls only read operations and scheduler queries. It never calls
dispatcher `advance` (even `advance(submit=False)` can mutate files), preflight,
submission, cancellation or native QC. It tolerates a moving live snapshot by
reporting read/provenance errors instead of interpreting missing data as success.
Optional report files must be outside the run and cannot overwrite existing
report files. No controller refresh is required to use a newly pulled auditor.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
git pull --ff-only origin composite
KB_AUDIT=$(mktemp -d /tmp/kinbot-r6-audit.XXXXXX)
/home/hram/KinBot/.venv/bin/python -B -m kinbot.run_audit \
  /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage \
  --output "$KB_AUDIT" --limit 2
printf '\nRead every task with: less %s/audit.txt\n' "$KB_AUDIT"
BASH
```

The terminal summary is bounded per live/ready group; `audit.txt` and
`audit.json` contain all task identities, dependencies and policy evidence.
Other live allocations outside materialized dispatch graphs are explicitly
reported, rather than silently omitted.
