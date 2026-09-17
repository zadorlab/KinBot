# CH4 composite dispatcher: first HPC smoke test

This is a single-molecule **dispatch and interface** test. It submits Gaussian
ASE/Sella geometry, Molpro native geometry, then independent Molpro, CFOUR,
MRCC, and Gaussian jobs. A green dispatcher status means that the processes
and declared files passed checks; scientific output parsing and ANL energy
assembly are later gates.

## 1. Clone the tested branch on a Slurm login node

After `composite` has been pushed to a remote you can read:

```bash
git clone --branch composite --single-branch https://github.com/zadorlab/KinBot.git
cd KinBot
git rev-parse HEAD
git status --short --branch
```

Compare the printed commit with the tested commit supplied alongside this
runbook. If the branch is in a fork, replace the URL with the fork URL. If the
cluster cannot access GitHub, transfer a Git bundle from the development
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

Load or activate the site's Miniforge installation first. The commands below
assume `mamba` is available. Replace `mamba` with `conda` if needed.

```bash
mamba create -y -c conda-forge -p "$PWD/.venv" \
  python=3.11 pip numpy scipy ase=3.29.0 networkx rmsd pytest
.venv/bin/python -m pip install 'sella==2.6.0'
.venv/bin/python -m pip install -e . --no-deps
.venv/bin/python -m pytest -q tests/test_anl_dispatch.py tests/test_theory_profiles.py
.venv/bin/python -c 'import sys, ase, importlib.metadata; print(sys.executable); print("ASE", ase.__version__); print("Sella", importlib.metadata.version("sella"))'
```

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
target partition. If required, add `"partition": "YOUR_PARTITION"` inside
each task's `resources`. The two geometry jobs are sequential; only after
the Molpro geometry is accepted can up to `max_nodes` independent jobs run
at once. Each task's Slurm script requests `--exclusive`.
Molpro CCSD(T) geometry and harmonic frequencies use numerical derivatives,
so check those two walltimes against the site's expected CH4 performance.

```bash
.venv/bin/python -m kinbot.anl.dispatch prepare ch4_dispatch.json ch4_run
cat ch4_run/tasks/l2_geometry/job.slurm
```

Use a new run-directory name if you change chemistry, task resources, or
the JSON after preparation. The dispatcher hashes its prepared specification.

## 4. Configure the site's licensed programs

Edit `ch4_run/site_setup.sh`. It is sourced by each batch job **after**
`set -euo pipefail`. For example, replace the names and paths below with the
actual site installations:

```bash
#!/usr/bin/env bash
# If the module function is unavailable in non-login batch shells, source
# your site's module initialization file here.
module load gaussian/YOUR_VERSION
module load molpro/YOUR_VERSION
module load cfour/YOUR_VERSION
module load mrcc/YOUR_VERSION
export CFOUR_GENBAS=/absolute/path/to/cfour/basis/GENBAS
```

Do not run the four QC programs on the login node. This check only verifies
module setup and executable/file discovery there; compute-node execution is
the actual test:

```bash
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run
```

Preflight must report `g16`, `molpro`, `xcfour`, and `dmrcc`, a valid
`CFOUR_GENBAS`, and an exclusive Slurm directive for each staged job. It also
imports ASE, Sella, and KinBot with the Python path pinned for batch jobs and
uses Slurm's `sbatch --test-only` to validate currently staged job scripts
without submitting them ([Slurm `sbatch` manual](https://slurm.schedmd.com/sbatch.html)).
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

To keep releasing ready jobs automatically, run the driver in a permitted
login or workflow session:

```bash
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run --interval 30
```

You can stop the driver and restart it. If persistent login-node processes
are disallowed, rerun `drive ch4_run --once` periodically. The driver records
Slurm job IDs, and a completed geometry is checked before downstream inputs
are created. `status` prints `waiting`, `staged`, `submitted`, `complete`, or
`failed` for every task. The parent `drive` command exits successfully only
when every task is `complete`.

Inspect a task with, for example:

```bash
cat ch4_run/tasks/l3_geometry/execution.json
tail -n 40 ch4_run/tasks/l3_geometry/slurm.stderr
tail -n 40 ch4_run/tasks/l3_geometry/l3_geometry.out
ls -lh ch4_run/tasks/l3_geometry/l3_geometry.log \
  ch4_run/tasks/l3_geometry/l3_geometry.xyz
```

If a task fails, read its `execution.json`, `slurm.stderr`, launcher stderr,
and native output. After correcting only the site module setup or a transient
cluster issue, archive the failed attempt and restage that one task:

```bash
.venv/bin/python -m kinbot.anl.dispatch retry ch4_run FAILED_TASK_ID
.venv/bin/python -m kinbot.anl.dispatch preflight ch4_run
.venv/bin/python -m kinbot.anl.dispatch drive ch4_run --once
```

The original files remain under `ch4_run/attempts/FAILED_TASK_ID/1/`. For
an input-method or resource change, edit the source JSON and prepare a fresh
run instead. If a submission is marked `submitting`, reconcile its Slurm job
ID manually before any new submission; the driver deliberately stops.

## 6. Save results for the next implementation pass

Once all tasks finish, retain `workflow.json`, `state.json`, each task's
`task.json`, `execution.json`, `job.slurm`, input, stdout/stderr, native
output, and Slurm logs. In particular, retain Molpro geometry `.out`, `.log`,
and final `.xyz`. Record `git rev-parse HEAD`, `module list`, and version
banners from all four codes. One way to package the run is:

```bash
tar -czf ch4_run_results.tgz -C ch4_run workflow.json state.json tasks
sha256sum ch4_run_results.tgz
```

Also retain `ch4_run/attempts` if there were retries. Review the archive under
your site's sharing rules before transferring it. The next pass will compare
printed methods, energies, frequencies, and DBOC/VPT2 values, implement
versioned parsers, and then assemble the ANL expression and MESS handoff.
