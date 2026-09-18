# CH4 composite dispatcher: first HPC smoke test

This is a single-molecule **dispatch and interface** test. It submits Gaussian
ASE/Sella geometry, Molpro ASE/Sella geometry, then independent Molpro, CFOUR,
and Gaussian jobs. MRCC is deferred for this first offsite test. A green
dispatcher status means that the processes and declared files passed checks;
the Molpro geometry calculator also parses
energy and forces. Full scientific parsing and ANL energy assembly are later
gates.

## 1. Clone the tested branch on a Slurm login node

After `composite` has been pushed to a remote you can read:

```bash
git clone --branch composite --single-branch https://github.com/zadorlab/KinBot.git
cd KinBot
git rev-parse HEAD
git status --short --branch
```

Compare the printed commit with the tested commit supplied alongside this
runbook. If the branch is in a fork, replace the URL with the fork URL.

For an existing HPC checkout, update the branch before generating the CH4
specification:

```bash
cd ~/KinBot  # replace if the checkout is elsewhere
git switch composite
git pull --ff-only origin composite
git rev-parse --short HEAD
```

On Blodgett, do this before loading Python or QC modules. If Git's HTTPS
helper reports a `libhogweed.so.6` / `__gmpn_cnd_sub_n` symbol error after
modules have been loaded, retry the pull with the module library overrides
removed for that command only:

```bash
env -u LD_LIBRARY_PATH -u LD_PRELOAD git pull --ff-only origin composite
git rev-parse --short HEAD
```

This keeps the current shell's modules loaded. The [Linux dynamic loader](https://man7.org/linux/man-pages/man8/ld.so.8.html)
uses `LD_LIBRARY_PATH` when resolving shared libraries, and
[`env -u`](https://www.gnu.org/s/coreutils/manual/html_node/env-invocation.html)
removes a variable from the environment of the command it starts.

If you already generated `ch4_dispatch.json` from an older commit, regenerate
it after the pull. Prepare a new run directory; prepared task graphs are
immutable and older ones may still contain the deferred MRCC task.

If the cluster cannot access GitHub, transfer a Git bundle from the development
machine and clone it on the cluster:

```bash
# Development machine, in the KinBot checkout:
git bundle create composite.bundle composite
scp composite.bundle YOUR_USER@YOUR_HPC:~/

# Cluster login node:
git clone --branch composite ~/composite.bundle KinBot
cd KinBot
git rev-parse HEAD
```

The bundle must be made after the final test commit.

## 2. Install the Python environment on shared storage

On Blodgett, `module load python/3.7.3` exposes Conda 22.9.0. This is only
the bootstrap: the new `.venv` below uses Python 3.11. System curl reports
`/etc/ssl/certs/ca-certificates.crt` as its CA file. Use that **literal path**;
`/path/printed/after/CAfile` was a placeholder in earlier instructions and
must not be entered. If an earlier attempt saved that placeholder in
`~/.condarc`, the guarded replacement below repairs it while saving a backup.

```bash
module load python/3.7.3
test -r /etc/ssl/certs/ca-certificates.crt
if [ -f "$HOME/.condarc" ]; then
  cp -p "$HOME/.condarc" "$HOME/.condarc.before-kinbot"
  sed -i 's|/path/printed/after/CAfile|/etc/ssl/certs/ca-certificates.crt|g' "$HOME/.condarc"
fi
conda config --set ssl_verify /etc/ssl/certs/ca-certificates.crt
conda config --show ssl_verify
```

Continue only when `ssl_verify` prints that exact file path. Then create the
environment:

```bash
conda create -y -c conda-forge -p "$PWD/.venv" python=3.11 pip numpy scipy ase=3.29.0 networkx rmsd pytest
```

Continue only when Conda finishes successfully and `.venv/bin/python` exists.
Unload the bootstrap Python module so its `/opt/anaconda3` libraries do not
affect later programs, then install and check KinBot with the new interpreter:

```bash
module unload python/3.7.3
PIP_CERT=/etc/ssl/certs/ca-certificates.crt .venv/bin/python -m pip install 'sella==2.6.0'
PIP_CERT=/etc/ssl/certs/ca-certificates.crt .venv/bin/python -m pip install -e . --no-deps
.venv/bin/python -m pytest -q tests/test_molpro_ase.py tests/test_anl_dispatch.py tests/test_theory_profiles.py
.venv/bin/python -c 'import sys, ase, importlib.metadata; print(sys.executable); print("ASE", ase.__version__); print("Sella", importlib.metadata.version("sella"))'
```

### If environment creation still fails

The earlier `CustomValidationError` came from entering the placeholder path.
If the repair block above does not clear it, inspect `~/.condarc` for any
remaining `/path/printed/after/CAfile` text, then replace it with the literal
path above. Conda 22.9.0 does not support `ssl_verify: truststore`, but it does
support a CA bundle path. `conda config --set` writes to `~/.condarc` by
default. See the [Conda 22.9 configuration command](https://docs.conda.io/projects/conda/en/22.9.x/commands/config.html)
and [certificate settings](https://docs.conda.io/projects/conda/en/22.9.x/user-guide/configuration/use-condarc.html).

If Conda still reports an SSL error with
`/etc/ssl/certs/ca-certificates.crt`, the cluster may require a different CA
chain; ask site support for it. Do not disable certificate verification. The
missing `.venv/bin/python` messages are a consequence of the failed `conda
create`; run the pip and test commands only after environment creation succeeds.
The `PIP_CERT` assignments use the same CA bundle for pip's HTTPS requests,
as [pip documents](https://pip.pypa.io/en/stable/topics/https-certificates/).

Run preparation with this same `.venv/bin/python` on the cluster. The batch
scripts store its absolute path. The licensed QC executables need to be
available on compute nodes through modules or the site setup script below.
Sella 2.6.0 is installed from [PyPI](https://pypi.org/project/sella/)
because the [conda-forge `sella` package](https://anaconda.org/conda-forge/sella)
currently provides 2.1.0, below KinBot's 2.6.0 requirement. If the HPC login
node cannot reach PyPI, transfer a wheel for the same release and install it
with pip there.

## 3. Set resources and prepare a fresh run

For a fresh CH4 run on Slurm, generate the fixture with automatic resources:

```bash
.venv/bin/python examples/anl/ch4_dispatch.py ch4_dispatch.json --auto-resources
```

Set only `limits.max_nodes` to the number of exclusive nodes you may use
concurrently (the fixture starts at 3). `prepare` reads `sinfo`, selects the
shortest available partition that fits each walltime, and sizes the job from
the smallest eligible node in that partition. An explicit
`resources.partition` takes priority. Automatic jobs request `--mem=0` so
Slurm allocates all memory on their exclusive node. Molpro ranks receive at
least 1024 MW of stack plus the documented 200 MW per-rank overhead, with 15%
node headroom; preparation fails if even one rank cannot meet the minimum.
The default performance cap is 16 ranks, chosen from 1, 2, 4, 8, 12, or 16;
`resources.max_cores` and `resources.min_stack_mw` can adjust a particular
method after a timing benchmark. Gaussian/CFOUR use one Slurm task with the
selected shared-memory core count. The CH4 fixture's original fixed-resource
mode remains available by omitting `--auto-resources` for repeatable local
tests. An existing prepared run keeps its original resource specification.
On the reported Blodgett layout, `day-long-cpu` should fit the CH4 walltimes;
review the choice and resolved resources in `ch4_run/workflow.json` and each
staged `job.slurm`. Automatic sizing requires `sinfo` on the target site.
The two geometry jobs are sequential; only after
the Molpro ASE/Sella geometry is accepted can up to `max_nodes` independent jobs run
at once. Each task's Slurm script requests `--exclusive`.
Molpro `FORCE,NUMERICAL` is run once per Sella geometry step; harmonic
frequencies also use numerical derivatives. Check their walltimes against
the site's expected CH4 performance. Each Sella step writes its own Molpro
`.inp`, `.out`, `.log`, and gradient `.xyz`. The calculator passes Molpro's
`-g` option to request the detailed `.log` for every force step.
Sella writes `optimization.traj`, `optimization.log`, and `final.xyz`.
Molpro requests `--ntasks=<ranks> --cpus-per-task=1` from Slurm and passes the
same count with `-n`; `OMP_NUM_THREADS=1` and `MKL_NUM_THREADS=1` avoid
multiplying threads per rank. For Molpro 2024's default single-node disk
mode, the generated commands use per-process `-m`; check that the site's
`.molprorc` does not add `-M` or `-G`.

Load the QC modules you intend to use **before** preparation. Gaussian is
already visible as `g16` on the reported Blodgett login node, even though it
has no listed module. These commands are specific to that site; the dispatcher
does not contain these paths or module versions:

```bash
module load molpro/molpro24 cfour/2.1
command -v g16 molpro xcfour
.venv/bin/python -m kinbot.anl.dispatch prepare ch4_dispatch.json ch4_run
cat ch4_run/site_setup.sh
cat ch4_run/tasks/l2_geometry/job.slurm
```

Use a new run-directory name if you change chemistry, task resources, or
the JSON after preparation. The dispatcher hashes its prepared specification.

## 4. Inspect discovered setup and preflight

`prepare` writes `ch4_run/site_setup.sh` from the executables and loaded
modules visible on the login node. It records the exact Molpro and CFOUR
module names found in `LOADEDMODULES` and pins their executable directories
on `PATH`. From `g16`, it finds and sources the adjacent
`bsd/g16.profile`, sets `g16root`, and gives Gaussian a scratch directory.
It uses a writable `GAUSS_SCRDIR` first, then `SLURM_TMPDIR`, `SCRATCH`,
`TMPDIR`, or the task directory. The reported login-node `/scratch` is not
writable, so preflight may use a fallback even when compute nodes can use
scratch. Keep the persistent run directory under your chosen working area
(for example under `$HOME`); Gaussian scratch is temporary.
The generated setup sources `g16.profile` in a checked conditional because
profile initialization probes may return nonzero under the batch script's
`set -e`. A nonzero final profile status is reported explicitly.
From `xcfour`, it looks for `../basis/GENBAS` and sets `CFOUR_GENBAS`; the
dispatcher copies that file beside `ZMAT` when the CFOUR task runs. An
existing `CFOUR_GENBAS` takes priority. The generated script is sourced by
each batch job after `set -euo pipefail` and can be edited for unusual sites.
The reported `/opt/cfour/2.1/basis/GENBAS` and
`/opt/gaussian/g16/bsd/g16.profile` match these discovery rules. Check the
generated script for the paths actually present on your system.

If `site_setup.sh` lacks a program or needs special license setup, edit that
file and rerun preflight. If you change chemistry, resources, or partition,
edit the JSON and prepare a fresh run directory. MRCC is omitted from this
test because no MRCC module was listed.
[Molpro's MRCC interface](https://www.molpro.net/manual/doku.php?id=the_mrcc_program_of_m._kallay_mrcc)
also requires a separate MRCC installation with its executables on `PATH`.

Do not run the three QC programs on the login node. This check only verifies
module setup and executable/file discovery there; compute-node execution is
the actual test:

```bash
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run
```

Preflight checks each task's own module environment and must report `g16`,
`molpro`, and `xcfour`, a valid
`CFOUR_GENBAS`, and an exclusive Slurm directive for each staged job. It also
imports ASE, Sella, and KinBot with the Python path pinned for batch jobs and
uses Slurm's `sbatch --test-only` to validate currently staged job scripts
without submitting them ([Slurm `sbatch` manual](https://slurm.schedmd.com/sbatch.html)).
Partition discovery uses the `%P`, `%c`, `%m`, `%l`, and `%a` fields of
[`sinfo`](https://slurm.schedmd.com/sinfo.html). It cannot determine account
or QoS access; `sbatch --test-only` is the final scheduler check.
If your cluster exposes these programs only on compute nodes, use a short
site-approved test allocation to check the same setup and executable paths
before starting this workflow.

## 5. Submit and monitor from the login node

The first call submits only the Gaussian ASE/Sella geometry job:

```bash
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run --once
.venv/bin/python -m kinbot.anl.dispatch status ch4_run
squeue -u "$USER"
```

Wait for `ch4_run/tasks/l2_geometry/execution.json` and for that job to leave
`squeue`, then accept its geometry
and stage Molpro without submitting it yet. Check the newly staged Slurm script
and submit Molpro:

```bash
.venv/bin/python -c 'from kinbot.anl.dispatch import advance; advance("ch4_run", submit=False)'
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run --once
```

Wait for `ch4_run/tasks/l3_geometry/execution.json` and for that job to leave
`squeue`. Before releasing the
other programs, inspect the first force evaluation and Sella's final geometry:

```bash
cat ch4_run/tasks/l3_geometry/execution.json
cat ch4_run/tasks/l3_geometry/l3_geometry_step_0001.inp
tail -n 40 ch4_run/tasks/l3_geometry/l3_geometry_step_0001.out
tail -n 40 ch4_run/tasks/l3_geometry/l3_geometry_step_0001.log
ls -lh ch4_run/tasks/l3_geometry/l3_geometry_step_0001.xyz \
  ch4_run/tasks/l3_geometry/final.xyz \
  ch4_run/tasks/l3_geometry/optimization.traj
.venv/bin/python - <<'PY'
from pathlib import Path
from ase.io import read
from kinbot.ase_modules.calculators.molpro import parse_numerical_gradient, parse_output
task = Path('ch4_run/tasks/l3_geometry')
atoms = read(task / 'geometry.xyz')
print('first CCSD(T) energy, eV:', parse_output(task / 'l3_geometry_step_0001.out'))
print('first max |force|, eV/Angstrom:',
      abs(parse_numerical_gradient(task / 'l3_geometry_step_0001.out',
                                   task / 'l3_geometry_step_0001.xyz', atoms)).max())
PY
```

Confirm the first `.out` prints CCSD(T)/cc-pVTZ and normal termination.
Inspect the `.log` as well. In this standalone Molpro 2024
`FORCE,NUMERICAL` calculation the complete `Numerical gradient for
KB_GEOM_ENERGY` table is in `.out`; it appeared in the second Blodgett attempt.
The calculator reads that table in Hartree/Bohr, changes its sign to return
ASE forces in eV/Å, and maps the rows through Molpro's `PUT,XYZ` geometry.
The second live attempt reported `mppx mode, nproc=11` despite the requested
eight processes. Setting `MPPX=0` for the successful retry did not change the
total: its `.out` still reports 11 compute processes plus one helper. The
Blodgett Molpro launch script discards its own `-np` launcher argument under
Slurm; Intel Hydra then uses the Slurm task layout. New Molpro scripts request
one Slurm task per intended MPI rank. The default Molpro mppx path is restored
for numerical gradients. A Molpro output reporting more total MPI processes
than declared ranks now fails the task. Verify this fix with the small DZ
single point before releasing larger jobs.
Check the force sign and units against the native table and confirm that the
Sella geometry has the same atom order. If the installed Molpro
output differs from the documented format, stop and return that input/output
pair for a parser fix. The completed L3 result may be accepted now. Stage the
independent jobs, preflight, and submit only the small Molpro DZ rank probe:

```bash
.venv/bin/python -c 'from kinbot.anl.dispatch import advance; advance("ch4_run", submit=False)'
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run
grep -E '^#SBATCH --(ntasks|cpus-per-task|mem|exclusive)' ch4_run/tasks/molpro_dz_sp/job.slurm
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run --once --only molpro_dz_sp
```

When that job leaves `squeue`, verify the native `.out` reports no more than
the intended total MPI count and normal termination. Reconcile the completed
task, then release the other ready jobs:

```bash
grep -E 'Distribution of processes|Memory per process|Molpro calculation terminated' ch4_run/tasks/molpro_dz_sp/molpro_dz_sp.out
.venv/bin/python -c 'from kinbot.anl.dispatch import advance; advance("ch4_run", submit=False)'
.venv/bin/python -m kinbot.anl.dispatch status ch4_run
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run --once
```

To keep checking and releasing ready jobs automatically, run the driver in a
permitted login or workflow session:

```bash
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run --interval 30
```

You can stop the driver and restart it. If persistent login-node processes
are disallowed, rerun `drive ch4_run --once` periodically. The driver records
Slurm job IDs, and a completed geometry is checked before downstream inputs
are created. `status` prints `waiting`, `staged`, `submitted`, `complete`, or
`failed` for every task. The parent `drive` command exits successfully only
when every task is `complete`.

Inspect a task later with, for example:

```bash
cat ch4_run/tasks/l3_geometry/execution.json
tail -n 40 ch4_run/tasks/l3_geometry/slurm.stderr
tail -n 40 ch4_run/tasks/l3_geometry/l3_geometry_step_0001.out
cat ch4_run/tasks/l3_geometry/l3_geometry_step_0001.inp
cat ch4_run/tasks/l3_geometry/l3_geometry_step_0001.xyz
ls -lh ch4_run/tasks/l3_geometry/final.xyz \
  ch4_run/tasks/l3_geometry/optimization.traj
```

If a task fails, read its `execution.json`, `slurm.stderr`, launcher stderr,
and native output. After correcting site setup, a transient cluster issue,
or the calculator implementation while keeping the task chemistry and
resources fixed, archive the failed attempt and restage that one task:

```bash
.venv/bin/python -c 'from kinbot.anl.dispatch import advance; advance("ch4_run", submit=False)'
.venv/bin/python -m kinbot.anl.dispatch retry ch4_run FAILED_TASK_ID
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run --once
```

The first Blodgett Molpro attempt exited before creating `.out` or `.log`:
its step `.stderr` reported `PSM3 can't open nic unit` and an OFI failure
during `PMPI_Init`. For this one-node job, setting
`I_MPI_FABRICS=shm` in the `molpro)` branch of the editable
`site_setup.sh` asks Intel MPI to communicate within the node without
initializing the network fabric.
The second attempt with this setting completed the CH4 CCSD(T) numerical
gradient, confirming that the MPI workaround reached the calculation. It
then failed at `PUT,XYZGRAD` with Molpro 2024's message "gradient is not
available ... for saving". The input now uses `PUT,XYZ`; the calculator
reads the numerical-gradient table in `.out` and verifies atom order using
the XYZ geometry. The third attempt completed the Sella optimization. See
[Intel MPI fabric control](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-reference-linux/2021-14/communication-fabrics-control.html).
The dispatcher treats Slurm's `Invalid job id specified` response for a
finished job as inactive, so the `advance` call can mark its recorded
execution failed and permit `retry`.

The original files remain under `ch4_run/attempts/FAILED_TASK_ID/1/` (and
`2/` for the second attempt). Keep the single-node `I_MPI_FABRICS=shm` line
in the editable `site_setup.sh` for this Blodgett run.
If a submission is marked `submitting`, reconcile its Slurm job ID manually
before any new submission; the driver deliberately stops.

For the live `ch4_run_auto3`, first pull the updated `composite` branch and
continue from its completed L3 geometry. Do not retry the completed task or
regenerate its prepared workflow:

```bash
cd ~/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD git pull --ff-only origin composite
git rev-parse --short HEAD
grep -n I_MPI_FABRICS ch4_run_auto3/site_setup.sh
.venv/bin/python -c 'from kinbot.anl.dispatch import advance; advance("ch4_run_auto3", submit=False)'
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run_auto3
grep -E '^#SBATCH --(ntasks|cpus-per-task|mem|exclusive)' ch4_run_auto3/tasks/ccsdt_dz/job.slurm
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run_auto3 --once --only ccsdt_dz
```

This submits only the small DZ CCSD(T) job. Its historical fixture requests
four ranks and 16 GB; the new Slurm layout should be `--ntasks=4` and
`--cpus-per-task=1`. `--mem=0` applies to fresh auto-sized runs, not this
already prepared fixed-resource workflow. Wait for the probe to leave
`squeue`, then check:

```bash
cat ch4_run_auto3/tasks/ccsdt_dz/execution.json
grep -E 'Distribution of processes|Memory per process|Molpro calculation terminated' ch4_run_auto3/tasks/ccsdt_dz/ccsdt_dz.out
.venv/bin/python -c 'from kinbot.anl.dispatch import advance; advance("ch4_run_auto3", submit=False)'
.venv/bin/python -m kinbot.anl.dispatch status ch4_run_auto3
```

The total process count should be at most four, the task should become
`complete`, and the `.out` should show normal termination. If so, submit the
other independent jobs with `drive ch4_run_auto3 --once`; the driver admits
only `limits.max_nodes` exclusive jobs concurrently. If the process count is
still wrong, retain the DZ input/output and Slurm stderr for diagnosis; the
automatic count check will fail that task and hold the rest.

The Blodgett DZ probe passed: Molpro 2024.1 reported four total processes
(three compute, one helper), 225 MW per compute process, a
`-40.387076267138` Hartree RHF CCSD(T)/cc-pVDZ energy, and normal termination.
The input coordinates match the accepted L3 CH4 geometry. The historical task
name `ccsdt_dz` and variable `KB_CCSDT` are misleading: this calculation is
**CCSD(T), not CCSDT**, and must not enter a CCSDT correction. Fresh CH4
fixtures name this task `molpro_dz_sp`; the prepared workflow retains its
original task ID. With the process-count gate passed, continue the live run:

```bash
.venv/bin/python -c 'from kinbot.anl.dispatch import advance; advance("ch4_run_auto3", submit=False)'
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run_auto3
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run_auto3 --once
.venv/bin/python -m kinbot.anl.dispatch status ch4_run_auto3
```

The first driver pass can submit up to three ready tasks. Repeat `drive
ch4_run_auto3 --once` after jobs finish until all tasks complete, or use the
bounded polling driver described above. Review the actual F12b, harmonic,
DBOC, and VPT2 outputs before using any numerical result in an ANL expression.

If the Slurm queue is empty but a task still reads `submitted`, the older
`status` command is showing the saved state. Reconcile without submitting:

```bash
.venv/bin/python -c 'from kinbot.anl.dispatch import advance; advance("ch4_run_auto3", submit=False)'
.venv/bin/python -m kinbot.anl.dispatch status ch4_run_auto3
```

From the subsequent dispatcher update, `status` performs that reconciliation
itself under the run lock; `status --cached` retains the saved-state view for
offline inspection. A job that left Slurm without `execution.json` becomes
`failed`, so inspect its Slurm stderr and native output before retrying.

### CFOUR runtime failure observed in the first CH4 run

The first `cfour_dboc` attempt exited 127 before reading `ZMAT` because
`libgfortran.so.4` was missing from the executable's loader path. The
compute-node probe found a usable copy in the site's Anaconda installation.
Current dispatch code discovers and checks a matching Fortran runtime for
the CFOUR child automatically; normal users do not need to find or export a
library path. It leaves the Python/Slurm environment and other programs'
library paths alone. If a site has no usable copy, preflight reports the
missing SONAME and the `CFOUR_LIBGFORTRAN` override. A passing loader check
does not prove the DBOC calculation will complete.

On the already prepared `ch4_run_auto3`, update the branch and retry only the
failed leaf task. The immutable CH4 geometry and seven successful siblings
remain accepted:

```bash
cd ~/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD git pull --ff-only origin composite
git rev-parse --short HEAD
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run_auto3
.venv/bin/python -m kinbot.anl.dispatch retry ch4_run_auto3 cfour_dboc
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run_auto3
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run_auto3 --once --only cfour_dboc
.venv/bin/python -m kinbot.anl.dispatch status ch4_run_auto3
```

Run `status` again after the job leaves Slurm. If CFOUR is still failed,
review `tasks/cfour_dboc/execution.json`, `cfour.err`, `cfour.out`, and
`slurm.stderr` before retrying. `retry` archives each previous attempt under
`attempts/cfour_dboc/<attempt>`. Do not hand the DBOC to the composite expression
until the final native output and numerical value have been reviewed.

The next CFOUR attempt reached `xjoda` but exposed a second input problem:
CFOUR 2.1 read only the first 80 columns of the 123-character keyword line
and rejected the truncated `MU` keyword. The current branch fixes this at
staging time, including for the older prepared workflow. To retry that failed
leaf again, pull the new commit, then run:

```bash
cd ~/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD git pull --ff-only origin composite
.venv/bin/python -m kinbot.anl.dispatch retry ch4_run_auto3 cfour_dboc
sed -n '1,30p' ch4_run_auto3/tasks/cfour_dboc/ZMAT
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run_auto3
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run_auto3 --once --only cfour_dboc
```

The new `ZMAT` should contain separate `CALC`, `BASIS`, `DBOC`, `COORD`,
`UNITS`, `CHARGE`, `MULTIPLICITY`, `MEM_UNIT`, and `MEMORY_SIZE` lines. After
the job leaves Slurm, `status` reconciles it. Inspect the native DBOC line
and numerical value even if the dispatch gate is `complete`.

The corrected CH4 run now has all eight tasks `complete`. The CFOUR native
output prints both HF and MP1 DBOC summaries. The requested HF value is
`0.0025887093` Hartree (`568.156016 cm-1`); MP1 is
`0.0026718675` Hartree and is a separate level. To inspect the existing
completed output with the method-aware parser after pulling the current
branch, run:

```bash
.venv/bin/python -m kinbot.anl.results cfour-dboc \
  ch4_run_auto3/tasks/cfour_dboc/cfour.out --level HF
```

This command only reads `cfour.out`. It does not change the completed run or
submit another job. The complete text archive was subsequently reviewed.
After pulling the current branch, the other native results can be inspected
without resubmitting or rewriting the completed run:

```bash
cd ~/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD git pull --ff-only origin composite
.venv/bin/python -m kinbot.anl.results molpro-energy \
  ch4_run_auto3/tasks/f12_tz/f12_tz.out --method 'CCSD(T)-F12b' --basis cc-pVTZ-F12
.venv/bin/python -m kinbot.anl.results molpro-energy \
  ch4_run_auto3/tasks/f12_qz/f12_qz.out --method 'CCSD(T)-F12b' --basis cc-pVQZ-F12
.venv/bin/python -m kinbot.anl.results molpro-energy \
  ch4_run_auto3/tasks/ccsdt_dz/ccsdt_dz.out --method 'CCSD(T)' --basis cc-pVDZ
.venv/bin/python -m kinbot.anl.results molpro-harmonic \
  ch4_run_auto3/tasks/harmonic/harmonic.out --basis cc-pVTZ
.venv/bin/python -m kinbot.anl.results gaussian-vpt2 \
  ch4_run_auto3/tasks/gaussian_vpt2/vpt2.log --method B3LYP --basis cc-pVTZ
```

The Gaussian parser reports `review_required: true` because the native VPT2
output contains rotor/framework and unreliable cubic-force warnings. It
reports the correction `-137.32185 cm-1`; treat that value as provisional
until a convergence comparison is complete. The Molpro F12b values are
`-40.454906199189` and `-40.456608306474` Hartree; the harmonic ZPE is
`0.04479801` Hartree. The historical `ccsdt_dz` directory contains
conventional CCSD(T), `-40.387076267138` Hartree. Fresh prepared CH4
fixtures attach these parser kinds to their external tasks automatically;
the old run remains unchanged.

## 6. Save results for the next implementation pass

Once all tasks finish, retain `workflow.json`, `state.json`, each task's
`task.json`, `execution.json`, `job.slurm`, input, stdout/stderr, native
output, and Slurm logs. In particular, retain every Molpro geometry step's
`.inp`, `.out`, `.log`, and gradient `.xyz`, plus Sella's final `.xyz` and
trajectory. Record `git rev-parse HEAD`, `module list`, and version
banners from Gaussian, Molpro, and CFOUR. Keep the full run directory on the
HPC. CFOUR's binary scratch files made the first full archive 99 MB even
though its text outputs were small. For transfer and parser review, package
only text inputs, outputs, and state records:

```bash
find ch4_run_auto3 -type f \( -name '*.json' -o -name '*.out' \
  -o -name '*.log' -o -name '*.inp' -o -name '*.com' \
  -o -name '*.xyz' -o -name '*.xml' -o -name '*.err' \
  -o -name '*.stderr' -o -name '*.stdout' -o -name '*.slurm' \
  -o -name 'ZMAT' -o -name 'site_setup.sh' \) -print0 | \
  tar --null -czf ch4_review.tgz -T -
sha256sum ch4_review.tgz
```

This includes matching files under `attempts` while excluding large binary
scratch such as CFOUR's `IIII`, `MOINTS`, and `GAMLAM`. Review the archive under
your site's sharing rules before transferring it. The first CH4 native-output
audit is documented in [composite QC validation](composite_qc_validation.md).
The VPT2 convergence check and missing ANL higher-order components come next;
no ANL composite energy or MESS input has been accepted from this smoke test.
