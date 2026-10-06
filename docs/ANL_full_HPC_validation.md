# Profiled KinBot and ANL external-site validation

## Readiness decision

The completed v5 run validates profiled KinBot reaction handling, its VRC
correction dispatch, ROTD_py execution, Gaussian L2, Gaussian VPT2, Molpro
F12, and CFOUR DBOC interfaces. It does **not** validate the conventional
Molpro CCSD(T) base. Those inputs contained `UHF_UCCSD=1`; Molpro 2024.1
printed an `RHF-UCCSD(T)` total label for ethane while its native `(T)`
contribution was exactly zero and the total equaled UCCSD. Consequently the
v5 L3 geometry, harmonic calculation, conventional single points, and every
post-geometry result evaluated on that L3 geometry are chemistry-invalid.
Keep them as failure and restart evidence, but do not assemble an ANL result
from them.

The corrected path generates `rhf` followed by bare `uccsd(t)`, requires a
native triples contribution, and verifies CCSD + `(T)` against the selected
total. First run the source-geometry ethane CCSD(T)/cc-pVDZ probe and compare
it with the published `-79.582320541811` hartree value. Then make a fresh
current-base graph from the already accepted v5 L2 geometry. Only its new L3
geometry may seed a new post-geometry graph. No KinBot reaction search or
ROTD_py sampling needs to be repeated for this correction.

The corrected ethane continuation performs these operations:

1. KinBot uses FairChem UMA at L1 for the initial structure, conformer search,
   and a restricted C-C homolytic reaction search.
2. KinBot refines accepted stationary points with
   B2PLYP-D3(BJ)/cc-pVTZ through ASE/Sella at L2 and evaluates hindered rotors
   on that same L2 surface.
3. The accepted homolytic scission bypasses stationary-saddle frequency
   validation. KinBot prepares the VRC asymptote with Gaussian, dispatches the
   sampling/high-level Molpro corrections on exclusive nodes, and writes a
   runnable rotdPy input and executes its reduced Slurm/Molpro sampling.
4. A machine-readable gate requires the accepted channel, both methyl product
   entries, four normal parent hindered-rotor points, a consistent VRC
   correction record, and a completed, hash-verified rotdPy result.
5. The accepted L2 parent geometry enters the exclusive-node dispatcher.
6. Molpro performs the CCSD(T)/cc-pVTZ ASE/Sella L3 geometry calculation.
7. After that geometry succeeds, Molpro harmonic, F12/TZ, F12/QZ, and
   CCSD(T)/DZ jobs, CFOUR DBOC, and a frequency-only Gaussian VPT2 job are
   allowed to run concurrently up to the requested node limit.
8. Every native output is hash checked and reparsed. The F12 pair is CBS
   extrapolated as a verified component.

The methyl product minimum still needs a frequency calculation for its MESS
partition function. `hom_sci` itself has no stationary transition state and
therefore has no transition-state Hessian or one-imaginary-frequency test.
The generated rotdPy smoke input is intentionally small: one dividing-surface
distance, small temperature/energy/angular grids, at most eight samples, and
one concurrent Molpro sampling job. It validates execution, parsing, and
restart records; it does not
establish production VRC convergence.

The final audit must say:

```text
interface_complete_recipe_incomplete
```

That status means the requested programs, scheduler path, geometry handoff,
parsers, restart logic, and one CBS operation worked. It deliberately does
not label the partial result ANL1-F12.

## Why this continuation is profiled ANL0-F12

The original literature defines canonical ANL0, ANL0-F12, and ANL1. The
branch does not invent an ANL1-F12 equation from the method name. Its
implemented assembler follows ANL0-F12 and records the selected
`SCALE_TRIP=1` and B2PLYP-D3(BJ) VPT2 choices in a profiled recipe label. The
following remain outside this first assembled result:

- a pinned, citable equation and label for the requested F12/ANL1 variant;
- the distinct ANL1 a'5Z/a'6Z reference, QZ geometry, and TZ/QZ ZPE graph;
- a citable definition for any requested ANL1-F12 extension;
- automatic full-network L3 graph construction for every reference species,
  well, product, and stationary TS;
- uncertainty propagation and the production MESS handoff;
- an external stationary-TS/IRC/restart validation.

Until those gates pass, a run may request the `ANL1-F12` ladder head to select
the higher B2PLYP L2 profile, but it must fall back to the highest complete,
provenance-matched recipe. It must not rename a partial sum ANL1-F12.

## Bugs fixed for this validation

- Profiled Gaussian L2 calculations now honor
  `frequency_mode: native_hessian` for wells and transition states:
  ASE/Sella performs the geometry optimization, then one Gaussian frequency
  calculation runs at that fixed optimized geometry without a Gaussian
  optimization keyword. The previous path incorrectly launched ASE
  finite-difference force displacements.
- Completion polling for ASE/Sella jobs now reads the optimizer completion
  log, while native Gaussian frequency-recovery jobs continue to use their
  Gaussian log. This prevents a completed Slurm job from leaving the KinBot
  driver polling indefinitely.
- Gaussian/Sella hindered-rotor points retain their per-angle checkpoint when
  later force evaluations use `Guess=Read`. Constraint failures are written
  to the Sella log instead of being silently converted to a generic failed
  scan.
- Repeated identical products, such as the two methyl fragments from ethane
  homolysis, share one optimization state while remaining two entries in the
  product channel for correct stoichiometry.

- Profiled inputs no longer stop in `QuantumChemistry` with an unimplemented
  router. A persistent router sends reaction/conformer work to L1 and L2
  refinements and hindered rotors to L2.
- The router records each job's backend and scheduler ID in
  `.kinbot_theory_jobs.json`, so result readers and restarts return to the
  correct backend.
- A preset backend change no longer inherits an empty or unrelated executable;
  Gaussian defaults to `g16` unless the profile supplies a command.
- The Blodgett example explicitly selects the additive
  `slurm_partition.tpl`, so its `queue_name` is emitted as `--partition`.
  The shared legacy `slurm.tpl` retains master KinBot's `-q` behavior.
- PBS and Slurm jobs now invoke the exact Python interpreter that submitted
  them. A `.venv` installation therefore remains active on compute nodes even
  when the login shell was not activated.
- The validation graph keeps Gaussian VPT2 on the accepted L2 geometry and
  does not add a Gaussian optimization.
- The Molpro L3 geometry is an ASE/Sella optimization. Molpro supplies energy
  and numerical forces and retains `.out`, `.log`, `.xml`, and `.xyz`
  evidence for every step.
- All dispatcher jobs request exclusive nodes. Molpro rank counts are capped
  by the method-specific scaling limit and reduced when the node cannot meet
  the configured minimum memory per rank.
- Every Molpro invocation gets an isolated temporary repository. KinBot checks
  `KINBOT_MOLPRO_SCRATCH`, `SLURM_TMPDIR`, `SCRATCH`, `TMPDIR`, the user's
  KinBot cache, the task filesystem, and `/tmp` in that order, skipping a
  candidate that lacks the task's minimum free space. The default minimum is
  half the requested node memory or 4096 MB, whichever is larger. Set
  `KINBOT_MOLPRO_SCRATCH` when a site has a preferred high-capacity filesystem;
  the selected root and capacity are recorded in `execution.json` and the
  per-invocation directory is removed after Molpro exits.
- A completed scheduler job is reconciled from its immutable execution
  record, including the Slurm `Invalid job id` case after a job leaves the
  queue.

## First installation on the HPC

Clone the branch directly for a fresh checkout:

```bash
cd ~
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  git clone --branch composite --single-branch \
  https://github.com/zadorlab/KinBot.git KinBot
cd ~/KinBot
```

For an existing checkout, pull before loading QC modules. This avoids the
previously observed Git and `libhogweed` collision from a module-modified
`LD_LIBRARY_PATH`.

```bash
cd ~/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD git fetch origin composite
env -u LD_LIBRARY_PATH -u LD_PRELOAD git switch composite
env -u LD_LIBRARY_PATH -u LD_PRELOAD git pull --ff-only origin composite
git rev-parse --short HEAD
```

Initialize the user's Miniforge installation. Adjust the first path if
Miniforge is installed elsewhere.

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
command -v conda
conda --version
```

If `.venv` does not exist, create it. Use the cluster CA bundle that `curl -vI`
reports; on Blodgett it was `/etc/ssl/certs/ca-certificates.crt`.

```bash
cd ~/KinBot
export KINBOT_CA_FILE=/etc/ssl/certs/ca-certificates.crt
export SSL_CERT_FILE="$KINBOT_CA_FILE"
export REQUESTS_CA_BUNDLE="$KINBOT_CA_FILE"
export PIP_CERT="$KINBOT_CA_FILE"

conda config --set ssl_verify "$KINBOT_CA_FILE"
conda create -y -c conda-forge -p "$PWD/.venv" \
  python=3.11 pip numpy scipy ase=3.29.0 networkx rmsd pytest openbabel
```

Install or update KinBot and the optional FairChem dependency in the same
environment:

```bash
PIP_CERT="$KINBOT_CA_FILE" .venv/bin/python -m pip install --upgrade pip
PIP_CERT="$KINBOT_CA_FILE" .venv/bin/python -m pip install -e '.[fc]'
```

ATcT's public `atct` client is now a KinBot dependency and is installed by the
preceding command. ROTD_py remains private and is pinned as a KinBot
submodule. Verify that the HPC account's SSH key is registered with a GitHub
account that can read `zadorlab/ROTD_py`:

```bash
ssh -T git@github.com
```

The expected message names the GitHub user and says authentication succeeded.
GitHub deliberately returns status 1 because it provides no shell. If an older
HTTPS attempt was interrupted, remove only its incomplete submodule checkout:

```bash
git submodule deinit -f -- external/ROTD_py 2>/dev/null || true
rm -rf external/ROTD_py .git/modules/external/ROTD_py
```

Initialize and install the pinned revision into the KinBot environment:

```bash
env -u LD_LIBRARY_PATH -u LD_PRELOAD git submodule sync --recursive
env -u LD_LIBRARY_PATH -u LD_PRELOAD git submodule update --init --recursive
PIP_CERT="$KINBOT_CA_FILE" .venv/bin/python -m pip install -e external/ROTD_py
git -C external/ROTD_py rev-parse HEAD
.venv/bin/python - <<'PY'
import atct
from rotd_py.flux.fluxbase import FluxBase
from rotd_py.new_multi import Multi
from rotd_py.sample.multi_sample import MultiSample
from kinbot.rotdpy import ensure_available
print('ATcT and rotdPy imports OK')
print(ensure_available())
PY
```

The import name is `rotd_py`. KinBot fails before licensed VRC jobs when the
pinned package or one of its runtime dependencies is unavailable. The executor
records the installed package location and Git revision.

Prime and verify a small pinned ATcT API cache on the networked login node:

```bash
mkdir -p "$HOME/.cache/kinbot/atct"
.venv/bin/python -m kinbot.anl.atct \
  1.220 "$HOME/.cache/kinbot/atct" C '[H][H]' --refresh
```

KinBot temporarily ignores a nonstandard `socks://` value in `ALL_PROXY`
during ATcT requests while preserving the site's `HTTP_PROXY` and
`HTTPS_PROXY` settings, then restores the original environment.

The command must report version `1.220` and a SHA-256 digest. If the live API
has advanced, KinBot fails the version check so the new table can be reviewed
and deliberately pinned rather than silently changing CBH reference data.

FAIR Chemistry's UMA checkpoint is gated on Hugging Face. Request access to
`facebook/UMA`, create a token with read access to public gated repositories,
then log in without placing the token in a shell history:

```bash
.venv/bin/python -m pip install --upgrade huggingface_hub
unset HF_HUB_OFFLINE
.venv/bin/hf auth login
.venv/bin/hf auth whoami
```

A successful `whoami` does not itself grant access to UMA. A
`403 GatedRepoError` stating that the user is not in the authorized list means
the network and token authentication succeeded, but that individual Hugging
Face account has not been granted access. While logged into the same account
shown by `hf auth whoami`, visit <https://huggingface.co/facebook/UMA>, submit
or accept the repository access terms, and wait if the page reports a pending
manual review. For a fine-grained token, enable **Read access to contents of
all public gated repos you can access**, then run `hf auth login --force` with
the updated token. Gated access can only be requested through the browser.

Do not use Hugging Face's standalone installer on Blodgett: it discovers
`/opt/anaconda3/bin/python` (Python 3.7) instead of KinBot's Python 3.11
environment. If the CLI reports an unknown `socks://` proxy scheme, first
retain the site's HTTP/HTTPS proxies while omitting the incompatible generic
proxy for the login and model-download commands:

```bash
env | grep -i '_proxy'
env -u ALL_PROXY -u all_proxy .venv/bin/hf auth login
```

If the site exposes only a SOCKS proxy, install HTTPX's SOCKS support and use
the `socks5://` scheme required by HTTPX. For the Blodgett value observed in
the traceback:

```bash
.venv/bin/python -m pip install 'httpx[socks]'
export ALL_PROXY=socks5://proxy.ca.sandia.gov:80
export all_proxy="$ALL_PROXY"
.venv/bin/hf auth login
```

Download the model once on a networked login node. Omitting Blodgett's
incompatible generic SOCKS proxy still leaves its `HTTP_PROXY` and
`HTTPS_PROXY` settings in place. Save the returned path; FairChem otherwise
uses a separate `~/.cache/fairchem` Hub cache and will not find a checkpoint
downloaded into the CLI's default `~/.cache/huggingface/hub` cache while
offline:

```bash
unset HF_HUB_OFFLINE
env -u ALL_PROXY -u all_proxy \
  .venv/bin/hf download facebook/UMA checkpoints/uma-s-1p2.pt
```

Set the path printed by that command and verify the molecular task directly
from the checkpoint. The production batch jobs use the same path with
`HF_HUB_OFFLINE=1`, which the supplied run script sets by default:

```bash
export KINBOT_FAIRCHEM_MODEL=/path/printed/by/hf/download
test -r "$KINBOT_FAIRCHEM_MODEL"

HF_HUB_OFFLINE=1 .venv/bin/python - "$KINBOT_FAIRCHEM_MODEL" <<'PY'
import sys

from fairchem.core import FAIRChemCalculator
from ase import Atoms
from kinbot.fairchem_utils import load_predictor

predictor = load_predictor(sys.argv[1], 'cpu')
atoms = Atoms('H2', positions=[[0., 0., 0.], [0., 0., 0.74]])
atoms.info.update({'charge': 0, 'spin': 1})
atoms.calc = FAIRChemCalculator(predictor, task_name='omol')
print('FairChem UMA H2 energy (eV):', atoms.get_potential_energy())
PY
```

Using the local checkpoint path avoids Hub cache lookup and proxy handling on
compute nodes. KinBot still removes a nonstandard `socks://` generic proxy
automatically during an offline named-model cache lookup and restores the
environment immediately afterward.

Official FairChem installation and model-access instructions are at
<https://github.com/FAIR-Chem/fairchem/blob/main/docs/core/install.md>; its
current quick start identifies `uma-s-1p2` and `omol` at
<https://github.com/FAIR-Chem/fairchem/blob/main/README.md>.

## Update an existing checkout

```bash
cd ~/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  git -c submodule.recurse=false pull --ff-only origin composite
export KINBOT_CA_FILE=/etc/ssl/certs/ca-certificates.crt
export PIP_CERT="$KINBOT_CA_FILE"
PIP_CERT="$KINBOT_CA_FILE" .venv/bin/python -m pip install -e '.[fc]'
ssh -T git@github.com
env -u LD_LIBRARY_PATH -u LD_PRELOAD git submodule sync --recursive
env -u LD_LIBRARY_PATH -u LD_PRELOAD git submodule update --init --recursive
PIP_CERT="$KINBOT_CA_FILE" .venv/bin/python -m pip install -e external/ROTD_py
```

Run the local suite before spending licensed-code allocation:

```bash
cd ~/KinBot
mkdir -p .mpl-cache
MPLCONFIGDIR="$PWD/.mpl-cache" \
  .venv/bin/python -m pytest -q --ignore=tests/test_kinbot.py
```

`tests/test_kinbot.py` is a pre-existing collection/import harness rather than
the runnable unit suite. Every collected test must pass before the HPC run;
the exact count changes as this branch adds regressions.

## Load external programs and run

Load programs only after Git operations. Gaussian has no module on Blodgett,
so source its installed profile.

```bash
cd ~/KinBot
module load molpro/molpro24 cfour/2.1
source /opt/gaussian/g16/bsd/g16.profile
command -v g16 molpro xcfour sbatch squeue sinfo
```

The example defaults FairChem to CPU because the reported partitions did not
advertise GPUs. The driver watches and advances dependent jobs. Run it from a
persistent login session, or use `nohup` when `tmux` is unavailable:

```bash
cd ~/KinBot
export KINBOT_PYTHON="$PWD/.venv/bin/python"
bash examples/anl/ethane_profiled_hpc/run.sh \
  day-long-cpu 3 "$KINBOT_FAIRCHEM_MODEL"
```

The equivalent detached launch is:

```bash
nohup bash examples/anl/ethane_profiled_hpc/run.sh \
  day-long-cpu 3 "$KINBOT_FAIRCHEM_MODEL" \
  > ethane_profiled_hpc_driver.log 2>&1 < /dev/null &
echo $! > ethane_profiled_hpc_driver.pid
disown
```

The first argument is the Slurm partition. The second is the maximum number of
simultaneous exclusive nodes for both the VRC Molpro correction stage and the
ANL dispatcher after the L3 geometry succeeds. The
optional third argument is the FairChem registered model name or local
checkpoint path; the local path is required for the documented offline HPC
run. The test directory defaults to `~/KinBot/ethane_profiled_hpc_run_v5`;
override it with `KINBOT_PROFILED_TEST_DIR`.

Monitor either layer with:

```bash
cd ~/KinBot
squeue -u "$USER"
.venv/bin/python -m kinbot.anl.dispatch status \
  ethane_profiled_hpc_run_v5/anl_interface
.venv/bin/python -m kinbot.anl.dispatch status \
  ethane_profiled_hpc_run_v5/vrctst/molpro/dispatch
tail -f ethane_profiled_hpc_run_v5/kinbot.log
```

Before a downstream dispatcher is created, its status is
`{"workflow": "not_prepared"}`. Inspect the generated input and completed
ROTD_py execution record with:

```bash
ls -lh ethane_profiled_hpc_run_v5/rotdPy/*.py
cat ethane_profiled_hpc_run_v5/rotdPy/*.execution.json
cat ethane_profiled_hpc_run_v5/rotdPy/*.rotdpy.json
```

The command is restartable. Rerun the same `run.sh` command after an
interruption. It reuses `kinbot.db`, the persistent L1/L2 route manifest, and
the dispatcher's immutable state.

The v5 failure that ended with `Molecular charge must be an integer` happened
while validating the first VRC dispatcher specification, before the dispatch
directory or a Molpro job was created. After pulling the fix, resume that same
v5 directory with `run.sh`; retain the completed conformer, rotor, product,
and Gaussian VRC records. No dispatcher retry or run-directory cleanup is
needed for this specific failure.

The following v5 attempt created both VRC tasks but exposed a single-node
Intel MPI transport failure (`PSM3`, `PMPI_Init`, and OFI endpoint errors).
After updating KinBot, retry both failed tasks explicitly and then restart the
top-level driver. The child runtime now selects shared-memory MPI for a
single-node Molpro allocation unless the site or user already selected a
fabric. The same default is present in new dispatcher setup and rotdPy sample
scripts.

Both VRC retries then completed. The first rotdPy handoff exited before sample
submission because `vrc_tst_noscan` passed its single asymptotic point to a
cubic correction spline. An asymptote-only input now omits the radial
correction while retaining its sampling-level asymptotic energy. A real scan
must provide at least three scan points plus the asymptote. The deferred
`qu.tpl` also now escapes the braces in its shell `I_MPI_FABRICS` default
through rotdPy's formatting pass. Pull the correction and rerun the same v5
driver; do not retry or remove the already complete VRC dispatcher.

The first ANL fanout then completed L3 geometry, harmonic frequencies,
F12/TZ, and conventional CCSD(T)/DZ. F12/QZ failed while Molpro was writing
distributed integral files below `/tmp`: the compute node exposed only 40 GB
there, and rank three reported `I/O error. Perhaps full disk?`. Updated KinBot
capacity-checks its Molpro repository before launch. On Blodgett, where the
batch environment defines neither `SCRATCH` nor `SLURM_TMPDIR`, it selects
`$HOME/.cache/kinbot/molpro` instead of the undersized node `/tmp`. Retry only
`f12_qz`; the completed sibling calculations remain reusable.

For a failed dispatcher task, inspect its `execution.json`, `slurm.stderr`,
and native output before retrying:

```bash
.venv/bin/python -m kinbot.anl.dispatch retry \
  ethane_profiled_hpc_run_v5/anl_interface TASK_ID
.venv/bin/python -m kinbot.anl.dispatch drive \
  ethane_profiled_hpc_run_v5/anl_interface --interval 20
```

Use `retry` only when the native calculation failed. If the native program
completed and only KinBot's result parser rejected the output, update KinBot
and reparse the preserved output without submitting another licensed job:

```bash
.venv/bin/python -m kinbot.anl.dispatch reparse \
  ethane_profiled_hpc_run_v5/anl_interface TASK_ID
```

This operation requires the task's declared native success marker, unchanged
staged input and geometry, and an `execution.json` traceback showing that the
failure occurred in `parse_result`. It preserves the original failure as
`execution.failed.json`, hashes that record with the accepted artifacts, and
marks the task complete only after all native output checks and the corrected
parser pass.

Review the final machine-readable report:

```bash
cat ethane_profiled_hpc_run_v5/anl_interface_audit.json
cat ethane_profiled_hpc_run_v5/kinbot_gate.json
```

The report is acceptable for this validation only when every task is
`complete`, the native parsers are listed, a finite F12 CBS value is present,
and the top-level status remains `interface_complete_recipe_incomplete`.
The KinBot gate must separately report `kinbot_reaction_complete`.

Do not resume an older run directory that already recorded failed HIR or
`hom_sci` state. Preserve it for diagnosis and select a fresh directory after
pulling this fix, for example:

```bash
export KINBOT_PROFILED_TEST_DIR="$PWD/ethane_profiled_hpc_run_v5"
```

The completed v5 B2PLYP-D3(BJ)/cc-pVTZ VPT2 job exposed a parser-only defect.
Gaussian terminated normally after the complete anharmonic analysis, but its
short route near the beginning omitted the long `EmpiricalDispersion`
keyword. The keyword and `GD3BJ` value were retained in Gaussian's archive
record near the end of the 1.7 MB log. The parser now searches the complete
normally terminated output for the dispersion model while continuing to
check the requested method, basis, frequency mode, named ZPE components, and
mode table. After pulling this correction, recover that exact calculation
with:

```bash
cd ~/KinBot
.venv/bin/python -m kinbot.anl.dispatch reparse \
  ethane_profiled_hpc_run_v5/anl_interface gaussian_vpt2
.venv/bin/python -m kinbot.anl.dispatch status \
  ethane_profiled_hpc_run_v5/anl_interface

cd ethane_profiled_hpc_run_v5
../.venv/bin/python -m kinbot.anl.validation audit anl_interface \
  | tee anl_interface_audit.json
```

No Gaussian resubmission is part of this recovery. The parsed result will
retain Gaussian's VPT2 warnings and set `review_required`; this interface
validation records that scientific review flag without treating a normally
terminated job as an execution failure. `reparse` is idempotent: when this
task is already complete, it verifies and returns the existing record.

The v5 Molpro jobs predate the unrestricted-command policy. Its audit must
therefore report `interface_complete_legacy_recipe_incompatible`, list the
legacy Molpro tasks, and expose the old F12 CBS only as a diagnostic. That
status is expected and prevents those electronic/ZPE values from entering a
current ANL recipe. The accepted L2 geometry and Gaussian VPT2 result remain
reusable.

## Continue the completed v5 ethane calculation

Update KinBot before loading QC modules, reinstall the editable packages, and
run the local regression suite:

```bash
cd ~/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  git -c submodule.recurse=false pull --ff-only origin composite
git submodule sync --recursive
git submodule update --init --recursive

export KINBOT_CA_FILE=/etc/ssl/certs/ca-certificates.crt
export PIP_CERT="$KINBOT_CA_FILE"
.venv/bin/python -m pip install -e . --no-deps --no-build-isolation
.venv/bin/python -m pip install -e external/ROTD_py \
  --no-deps --no-build-isolation
MPLCONFIGDIR="$PWD/.mpl-cache" \
  .venv/bin/python -m pytest -q --ignore=tests/test_kinbot.py
```

`--no-build-isolation` is required on an offline or proxy-restricted login
node. Without it, pip may try to download the `pyproject.toml` build
requirements even when `--no-deps` is present.

Make all native programs visible during preparation. The Blodgett path below
is site configuration rather than a repository default. `KINBOT_MRCC_ROOT`
may be replaced by a module, normal `PATH`, `MRCC_ROOT`, `MRCC_HOME`,
`EBROOTMRCC`, or `KINBOT_MRCC_COMMAND`.

```bash
module load molpro/molpro24 cfour/2.1
source /opt/gaussian/g16/bsd/g16.profile
export KINBOT_MRCC_ROOT=/opt/mrcc

command -v molpro xcfour g16 sbatch squeue sinfo
test -x "$KINBOT_MRCC_ROOT/dmrcc"
test -x "$KINBOT_MRCC_ROOT/scf"
test -x "$KINBOT_MRCC_ROOT/mrcc"
```

Before continuing the profiled geometry, run the inexpensive conventional
CCSD(T)/cc-pVDZ calculation at the published ethane TZ geometry. This catches
method-selection and parsing errors without mixing in a geometry difference:

```bash
cd ~/KinBot
.venv/bin/python -m kinbot.anl.literature prepare-higher-order \
  ethane-tz-2017 ethane_tz_ccsdt_dz_449cfbd \
  --task ccsdt_dz --max-nodes 1 --partition day-long-cpu
.venv/bin/python -m kinbot.anl.dispatch preflight \
  ethane_tz_ccsdt_dz_449cfbd
.venv/bin/python -m kinbot.anl.dispatch drive \
  ethane_tz_ccsdt_dz_449cfbd --interval 20
.venv/bin/python -m kinbot.anl.literature compare-run \
  ethane-tz-2017 ethane_tz_ccsdt_dz_449cfbd \
  | tee ethane_tz_ccsdt_dz_449cfbd_literature.json
```

The comparison must pass near `-79.582320541811` hartree, and the generated
input must contain `rhf` followed by bare `uccsd(t)`.

The Blodgett validation at revision `449cfbd` returned
`-79.582320352486` hartree, an error of `1.89325e-7` hartree, and passed the
pinned tolerance.

Next prepare a current base graph from the accepted and hash-verified v5 L2
geometry. It reruns the L3 Sella optimization through Molpro's normal
RHF-referenced unrestricted CCSD(T) path, then fans out the current F12 TZ/QZ,
harmonic, and CFOUR DBOC tasks. It does not repeat the L2 optimization or
Gaussian VPT2 calculation.

```bash
cd ~/KinBot
base="$PWD/ethane_profiled_hpc_run_v5"

.venv/bin/python -m kinbot.anl.validation \
  prepare-current-base-from-run \
  "$base/anl_interface" "$base/anl_current_base_ethane_449cfbd" \
  --geometry-task l2_geometry \
  --max-nodes 3 --partition day-long-cpu

.venv/bin/python -m kinbot.anl.dispatch preflight \
  "$base/anl_current_base_ethane_449cfbd"
```

Run the current base first. Only after its new L3 geometry is complete can the
post-geometry graph be prepared from that exact geometry. With `--anl0-only`,
the graph contains the ANL0 higher-order term and every common core-valence and
scalar-relativistic single point. One dispatcher therefore fans them out
together while enforcing a shared limit of three exclusive nodes. No `tmux`
installation is needed:

```bash
nohup bash -lc '
  set -euo pipefail
  cd "$HOME/KinBot"
  base="$PWD/ethane_profiled_hpc_run_v5"
  .venv/bin/python -m kinbot.anl.dispatch drive \
    "$base/anl_current_base_ethane_449cfbd" --interval 20
  .venv/bin/python -m kinbot.anl.validation audit-current-base \
    "$base/anl_current_base_ethane_449cfbd" \
    > "$base/anl_current_base_ethane_449cfbd_audit.json"
  .venv/bin/python -m kinbot.anl.validation \
    prepare-post-geometry-from-run \
    "$base/anl_current_base_ethane_449cfbd" \
    "$base/anl_post_geometry_ethane_449cfbd" \
    --max-nodes 3 --anl0-only
  .venv/bin/python -m kinbot.anl.dispatch preflight \
    "$base/anl_post_geometry_ethane_449cfbd"
  .venv/bin/python -m kinbot.anl.dispatch drive \
    "$base/anl_post_geometry_ethane_449cfbd" --interval 20
  .venv/bin/python -m kinbot.anl.validation audit-anl0-post-geometry \
    "$base/anl_post_geometry_ethane_449cfbd" \
    > "$base/anl_post_geometry_ethane_449cfbd_anl0_audit.json"
' > ethane_anl_449cfbd.log 2>&1 < /dev/null &
echo $! > ethane_anl_449cfbd.pid
disown
```

Monitor without mutating either graph:

```bash
base="$PWD/ethane_profiled_hpc_run_v5"
squeue -u "$USER" -o "%.18i %.30j %.2t %.10M %.10l %R"
.venv/bin/python -m kinbot.anl.dispatch status \
  "$base/anl_current_base_ethane_449cfbd"
.venv/bin/python -m kinbot.anl.dispatch status \
  "$base/anl_post_geometry_ethane_449cfbd"
tail -f ethane_anl_449cfbd.log
```

The completed `449cfbd` continuation produced a common L3 geometry hash of
`3b499b5022d557e6af706cca1a68949b57c9b6ade480832c60a57e464b26919b`.
The current-base audit completed the F12 TZ/QZ reference, harmonic ZPE, and
CFOUR DBOC. The post-geometry audit completed the unrestricted
CCSD(T)/cc-pVDZ and direct-MRCC CCSDT(Q)/cc-pVDZ pair, both core-valence
pairs, and both scalar-relativistic tasks. Against the pinned ethane source,
the harmonic ZPE, DBOC, and CCSDT(Q)-CCSD(T) increment differed by
`3.0e-8`, `-6.9e-9`, and `-9.9318e-8` hartree, respectively.

### Superseded v5 recovery record

The old `fb3eab2` current-base graph and its dependent
`anl_post_geometry_ethane` graph used `UHF_UCCSD=1`. Their ethane output
reported zero for `(T)`, so their conventional coupled-cluster energies and
all recipes depending on those energies are invalid. Keep the directories as
diagnostic provenance, but do not reparse, resume, or use them in a recipe.

An early CFOUR `xncc` CCSDT(Q) route and the first direct-MRCC restart attempts
also remain only as provenance. Fresh graphs use direct MRCC with RHF or
semicanonical ROHF orbitals for high-order terms. The `--anl0-only` graph omits
the ANL1-only CCSDT(Q)/cc-pVTZ and CCSDTQ(P)/cc-pVDZ jobs. Chemistry validation
still requires comparison of the affordable ethane ANL0-F12 components and
the final 0 K heat of formation with the 2017 source workbook.

The literature command compares only tasks present in a run. Apply it to the
base/interface run as well to cover harmonic, F12, and DBOC components. Exact
source comparison requires the source geometry and reference convention; a
result from a deliberately changed production profile must be reported as
such instead of relaxing the tolerance.

Prepare the short source-matched MRCC checks from the published QZ geometries
after the ethane correction fan-out is under control:

```bash
cd ~/KinBot
.venv/bin/python -m kinbot.anl.literature prepare-higher-order \
  methane-qz-2017 methane_qz_higher_2017 --max-nodes 2
.venv/bin/python -m kinbot.anl.dispatch preflight methane_qz_higher_2017
.venv/bin/python -m kinbot.anl.dispatch drive \
  methane_qz_higher_2017 --interval 20
.venv/bin/python -m kinbot.anl.validation audit-higher-order \
  methane_qz_higher_2017 | tee methane_qz_higher_2017_audit.json
.venv/bin/python -m kinbot.anl.literature compare-run \
  methane-qz-2017 methane_qz_higher_2017 \
  | tee methane_qz_higher_2017_literature.json
```

For methyl, first run only the source-comparable conventional component:

```bash
.venv/bin/python -m kinbot.anl.literature prepare-higher-order \
  methyl-qz-2017 methyl_qz_ccsdt_dz_2017 \
  --task ccsdt_dz --max-nodes 1
.venv/bin/python -m kinbot.anl.dispatch preflight methyl_qz_ccsdt_dz_2017
.venv/bin/python -m kinbot.anl.dispatch drive \
  methyl_qz_ccsdt_dz_2017 --interval 20
.venv/bin/python -m kinbot.anl.literature compare-run \
  methyl-qz-2017 methyl_qz_ccsdt_dz_2017 \
  | tee methyl_qz_ccsdt_dz_2017_literature.json
```

The published methyl CCSDT(Q) and CCSDTQ(P) targets used a UHF determinant,
as stated by the paper. The current profile deliberately uses semicanonical
ROHF, so those absolute values are retained as source provenance but are not
tight numerical acceptance targets for the changed reference. Stage the
modern `ccsdtq_dz` and `ccsdtqp_dz` pair separately and use
`audit-higher-order` to validate input, execution, parsing, and the common
geometry correction. Preserve one real interrupted MRCC attempt and
successful `resume-mrcc` continuation as restart evidence; do not manufacture
interruptions for every species.

```bash
.venv/bin/python -m kinbot.anl.literature prepare-higher-order \
  methyl-qz-2017 methyl_qz_rohf_higher \
  --task ccsdtq_dz --task ccsdtqp_dz --max-nodes 2
.venv/bin/python -m kinbot.anl.dispatch preflight methyl_qz_rohf_higher
.venv/bin/python -m kinbot.anl.dispatch drive \
  methyl_qz_rohf_higher --interval 20
.venv/bin/python -m kinbot.anl.validation audit-higher-order \
  methyl_qz_rohf_higher | tee methyl_qz_rohf_higher_audit.json
```

The completed ethane VPT2 parser recorded native Gaussian warnings. Assembly
therefore stops until those warnings and the mode table are reviewed. First
print the warnings and exact native-output hash:

```bash
base="$PWD/ethane_profiled_hpc_run_v5"
.venv/bin/python - <<'PY'
import json
from pathlib import Path

path = Path('ethane_profiled_hpc_run_v5/anl_interface/tasks/gaussian_vpt2/execution.json')
record = json.loads(path.read_text())
print('native_output_sha256 =', record['artifacts']['vpt2.log'])
print('warnings =')
for warning in record['details']['parsed_result']['warnings']:
    print(' -', warning)
PY
```

If the exact output is accepted after inspection, create
`$base/ethane_vpt2_review.json` with this schema, substituting the printed
hash, reviewer, and concrete scientific rationale:

```json
{
  "schema": 1,
  "task_id": "gaussian_vpt2",
  "native_output_sha256": "PASTE_THE_PRINTED_SHA256",
  "decision": "accept",
  "reviewer": "REVIEWER_NAME",
  "rationale": "DESCRIBE_WHICH_WARNINGS_AND_MODES_WERE_REVIEWED"
}
```

Then assemble the profiled ANL0-F12 value. The zero spin-orbit entry below is
an explicit state policy for this nondegenerate closed-shell validation and
is recorded in the result; change it if the selected production convention
uses another state-specific value.

```bash
base="$PWD/ethane_profiled_hpc_run_v5"
.venv/bin/python -m kinbot.anl.validation assemble-anl0-f12 \
  "$base/anl_interface" \
  "$base/anl_post_geometry_ethane_449cfbd" \
  "$base/anl_post_geometry_ethane_449cfbd" \
  "$base/ethane_profiled_anl0_f12_449cfbd.json" \
  --state-id ethane-singlet \
  --spin-orbit-hartree 0.0 \
  --spin-orbit-backend known_zero \
  --spin-orbit-source "nondegenerate closed-shell ethane validation policy" \
  --base-run "$base/anl_current_base_ethane_449cfbd" \
  --vpt2-review "$base/ethane_vpt2_review.json"
cat "$base/ethane_profiled_anl0_f12_449cfbd.json"
```

### Generate the CBH-0 reference graphs from SMILES

Reference fragments do not need a pre-existing KinBot database row. The
SMILES entry point creates a deterministic initial Cartesian structure and
then runs the same Gaussian/Sella L2 optimization, Molpro/Sella L3
optimization, and native component graph. The embedded coordinates are never
used as a final energy geometry. Diatomics use Cartesian Sella coordinates.

Prepare methane and hydrogen with one exclusive node available to each graph:

```bash
cd ~/KinBot
base="$PWD/ethane_profiled_hpc_run_v5"
refs="$base/cbh0_reference_graphs"
mkdir -p "$refs"

.venv/bin/python -m kinbot.anl.validation prepare-from-smiles \
  C "$refs/methane_interface" \
  --charge 0 --multiplicity 1 --max-nodes 1 --partition day-long-cpu
.venv/bin/python -m kinbot.anl.validation prepare-from-smiles \
  '[H][H]' "$refs/hydrogen_interface" \
  --charge 0 --multiplicity 1 --max-nodes 1 --partition day-long-cpu

.venv/bin/python -m kinbot.anl.dispatch preflight "$refs/methane_interface"
.venv/bin/python -m kinbot.anl.dispatch preflight "$refs/hydrogen_interface"
```

Run the independent graphs concurrently, then prepare each post-geometry
ANL0 fan-out from its exact accepted L3 geometry:

```bash
nohup bash -lc '
  set -uo pipefail
  cd "$HOME/KinBot"
  base="$PWD/ethane_profiled_hpc_run_v5"
  refs="$base/cbh0_reference_graphs"

  run_reference() (
    set -euo pipefail
    name=$1
    .venv/bin/python -m kinbot.anl.dispatch drive \
      "$refs/${name}_interface" --interval 20
    .venv/bin/python -m kinbot.anl.validation \
      prepare-post-geometry-from-run \
      "$refs/${name}_interface" "$refs/${name}_post" \
      --max-nodes 1 --partition day-long-cpu --anl0-only
    .venv/bin/python -m kinbot.anl.dispatch preflight "$refs/${name}_post"
    .venv/bin/python -m kinbot.anl.dispatch drive \
      "$refs/${name}_post" --interval 20
    .venv/bin/python -m kinbot.anl.validation audit \
      "$refs/${name}_interface" > "$refs/${name}_interface_audit.json"
    .venv/bin/python -m kinbot.anl.validation audit-anl0-post-geometry \
      "$refs/${name}_post" > "$refs/${name}_post_audit.json"
  )

  run_reference methane &
  methane_pid=$!
  run_reference hydrogen &
  hydrogen_pid=$!
  methane_rc=0
  hydrogen_rc=0
  wait "$methane_pid" || methane_rc=$?
  wait "$hydrogen_pid" || hydrogen_rc=$?
  if (( methane_rc || hydrogen_rc )); then
    printf "reference failures: methane=%s hydrogen=%s\\n" \
      "$methane_rc" "$hydrogen_rc" >&2
    exit 1
  fi
' > cbh0_reference_graphs.log 2>&1 < /dev/null &
echo $! > cbh0_reference_graphs.pid
disown
```

Inspect each `gaussian_vpt2` parsed result after the driver exits. Supply a
hash-bound review only when that record has `review_required: true`; an
unflagged result must not receive a review file. Assemble each accepted
reference with its interface run as the base and its post run as both the
higher-order and correction provider.

Once the methane and hydrogen composite JSON files exist, generate and solve
the ethane CBH-0 reaction directly from the three records:

```bash
base="$PWD/ethane_profiled_hpc_run_v5"
refs="$base/cbh0_reference_graphs"
mkdir -p "$HOME/.cache/kinbot/atct"

.venv/bin/python -m kinbot.anl.cbh generate CC --rung 0
.venv/bin/python -m kinbot.anl.cbh solve-records \
  CC "$base/ethane_cbh0_anl0_f12.json" --rung 0 \
  --energy "CC=$base/ethane_profiled_anl0_f12_449cfbd.json" \
  --energy "C=$refs/methane_profiled_anl0_f12.json" \
  --energy "[H][H]=$refs/hydrogen_profiled_anl0_f12.json" \
  --atct-version 1.220 \
  --atct-cache "$HOME/.cache/kinbot/atct"
cat "$base/ethane_cbh0_anl0_f12.json"
```

The solver recomputes the electronic, zero-point, and zero-K sums from every
input component, requires the exact same recipe label and state across the
CBH reaction, resolves only pinned 0 K gas-phase ATcT records, and stores all
source hashes in the result.

The methyl minimum can then be staged from the already accepted KinBot row,
using `--multiplicity 2`, and sent through the same base, higher-order, and
common-correction sequence. Before doing that, list the exact database job
name instead of assuming it:

```bash
base="$PWD/ethane_profiled_hpc_run_v5"
.venv/bin/python - <<'PY'
from ase.db import connect
for row in connect('ethane_profiled_hpc_run_v5/kinbot.db').select():
    if '150390060000000000002' in getattr(row, 'name', ''):
        print(row.name, row.data.get('status'))
PY
```

## Gates after this continuation

1. Validate the native CFOUR and MRCC method echoes and total-energy lines in
   both audit reports before accepting the assembled ethane value.
2. Complete the methyl doublet graph and compare ethane and methyl values to
   a pinned literature/ATcT release. Do not compare unlike 0 K and 298 K
   quantities.
3. Add the monatomic/diatomic reference-species graphs needed by the methyl
   CBH-0 reaction before claiming an automated CBH result. The current CBH
   solver and local table tests are complete, but the general QC graph still
   needs zero-mode handling for atomic H and the appropriate small-species
   component policy.
4. Use the existing completed ROTD_py surface as the interface result. A
   production CH3 + CH3 rate comparison needs converged sampling grids and the
   published CASPT2/cc-pVDZ surface plus its higher-level one-dimensional
   correction, followed by a MESS run with the accepted formation enthalpies.
5. Run a separate stationary-TS/IRC case before declaring the ANL/MESS path
   validated for ordinary transition-state kinetics.
