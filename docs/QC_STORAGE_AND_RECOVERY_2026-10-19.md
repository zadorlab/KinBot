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
