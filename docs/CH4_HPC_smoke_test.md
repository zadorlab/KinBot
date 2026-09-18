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

Generate the CH4-only task graph, then edit its JSON before preparation:

```bash
.venv/bin/python examples/anl/ch4_dispatch.py ch4_dispatch.json
```

In `ch4_dispatch.json`, set `limits.max_nodes` to the number of exclusive
nodes you may use concurrently (the example is 3). Check every task's
`resources.cores`, `resources.memory_mb`, and `resources.walltime` against the
site's limits. When `resources.partition` is absent, `prepare` reads `sinfo`
and selects the available partition with the shortest time limit that fits
that task's cores, memory, and walltime. An explicit `"partition": "NAME"`
inside a task's `resources` takes priority. On the reported Blodgett layout,
the default `short-cpu` partition permits only 30 minutes, so the CH4 tasks
need a longer partition; discovery should select `day-long-cpu`. Review the
actual choice in `ch4_run/workflow.json` and each staged `job.slurm`. If
`sinfo` is unavailable during preparation, set the partition explicitly.
The two geometry jobs are sequential; only after
the Molpro ASE/Sella geometry is accepted can up to `max_nodes` independent jobs run
at once. Each task's Slurm script requests `--exclusive`.
Molpro `FORCE,NUMERICAL` is run once per Sella geometry step; harmonic
frequencies also use numerical derivatives. Check their walltimes against
the site's expected CH4 performance. Each Sella step writes its own Molpro
`.inp`, `.out`, `.log`, and gradient `.xyz`. The calculator passes Molpro's
`-g` option to request the detailed `.log` for every force step.
Sella writes `optimization.traj`, `optimization.log`, and `final.xyz`.
Molpro runs `-n` MPI processes with `OMP_NUM_THREADS=1` and
`MKL_NUM_THREADS=1` inside its exclusive node to avoid multiplying threads
per process. For Molpro 2024's default single-node disk mode, the generated
commands use per-process `-m`; check that the site's `.molprorc` does not add
`-M` or `-G`.

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
from kinbot.ase_modules.calculators.molpro import parse_output, parse_xyzgrad
task = Path('ch4_run/tasks/l3_geometry')
atoms = read(task / 'geometry.xyz')
print('first CCSD(T) energy, eV:', parse_output(task / 'l3_geometry_step_0001.out'))
print('first max |force|, eV/Angstrom:',
      abs(parse_xyzgrad(task / 'l3_geometry_step_0001.xyz', atoms)).max())
PY
```

Confirm the first `.out` prints CCSD(T)/cc-pVTZ and normal termination.
Inspect the `.log` as well; in the supplied Molpro 2024
`OPTG` cases the detailed numerical gradient table is in `.log`, not `.out`.
Check that the XYZGRAD force sign and units agree with any native gradient table,
and that the Sella geometry has the same atom order. If the installed Molpro
output differs from the documented format, stop and return that input/output
pair for a parser fix. Then accept the L3 result, preflight the newly staged
jobs, and release them:

```bash
.venv/bin/python -c 'from kinbot.anl.dispatch import advance; advance("ch4_run", submit=False)'
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run
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
and native output. After correcting only the site module setup or a transient
cluster issue, archive the failed attempt and restage that one task:

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
`site_setup.sh` is the next site-specific retry to test; it asks Intel MPI
to communicate within the node without initializing the network fabric.
This setting is not yet a validated Molpro result. See
[Intel MPI fabric control](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-reference-linux/2021-14/communication-fabrics-control.html).
The dispatcher treats Slurm's `Invalid job id specified` response for a
finished job as inactive, so the `advance` call can mark its recorded
execution failed and permit `retry`.

The original files remain under `ch4_run/attempts/FAILED_TASK_ID/1/`. For
an input-method or resource change, edit the source JSON and prepare a fresh
run instead. If a submission is marked `submitting`, reconcile its Slurm job
ID manually before any new submission; the driver deliberately stops.

## 6. Save results for the next implementation pass

Once all tasks finish, retain `workflow.json`, `state.json`, each task's
`task.json`, `execution.json`, `job.slurm`, input, stdout/stderr, native
output, and Slurm logs. In particular, retain every Molpro geometry step's
`.inp`, `.out`, `.log`, and gradient `.xyz`, plus Sella's final `.xyz` and
trajectory. Record `git rev-parse HEAD`, `module list`, and version
banners from Gaussian, Molpro, and CFOUR. One way to package the run is:

```bash
tar -czf ch4_run_results.tgz -C ch4_run workflow.json state.json tasks
sha256sum ch4_run_results.tgz
```

Also retain `ch4_run/attempts` if there were retries. Review the archive under
your site's sharing rules before transferring it. The next pass will compare
printed methods, energies, frequencies, and DBOC/VPT2 values, implement
versioned parsers, and then assemble the ANL expression and MESS handoff.
