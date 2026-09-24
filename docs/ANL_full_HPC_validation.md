# Profiled KinBot and ANL external-site validation

## Readiness decision

The branch is ready for a **profiled workflow and QC-interface validation** on
the HPC. It is not yet ready to claim a complete ANL1-F12 energy, CBH heat of
formation, or production MESS result.

The supplied ethane run performs these real operations:

1. KinBot uses FairChem UMA at L1 for the initial structure, conformer search,
   and a restricted C-C homolytic reaction search.
2. KinBot refines accepted stationary points with
   B2PLYP-D3(BJ)/cc-pVTZ through ASE/Sella at L2 and evaluates hindered rotors
   on that same L2 surface.
3. The accepted L2 parent geometry enters the exclusive-node dispatcher.
4. Molpro performs the CCSD(T)/cc-pVTZ ASE/Sella L3 geometry calculation.
5. After that geometry succeeds, Molpro harmonic, F12/TZ, F12/QZ, and
   CCSD(T)/DZ jobs, CFOUR DBOC, and a frequency-only Gaussian VPT2 job are
   allowed to run concurrently up to the requested node limit.
6. Every native output is hash checked and reparsed. The F12 pair is CBS
   extrapolated as a verified component.

The final audit must say:

```text
interface_complete_recipe_incomplete
```

That status means the requested programs, scheduler path, geometry handoff,
parsers, restart logic, and one CBS operation worked. It deliberately does
not label the partial result ANL1-F12.

## Why a complete ANL1-F12 test is not available yet

The original literature defines canonical ANL0, ANL0-F12, and ANL1. The
current branch does not invent an ANL1-F12 equation from the method name. A
complete higher rung also needs calculations absent from the current external
graph:

- a pinned, citable equation and label for the requested F12/ANL1 variant;
- all-electron and frozen-core CCSD(T) TZ/QZ CBS calculations;
- scalar-relativistic DKH and nonrelativistic reference calculations;
- closed-shell CFOUR CCSDT(Q), with native output validation;
- CCSDTQ(P)/cc-pVDZ, currently assigned to MRCC and intentionally disabled;
- state-specific spin-orbit data or a validated zero policy;
- automatic L3 graph construction for every accepted well, product, and TS;
- assembly of each complete E0 value, CBH reaction selection, ATcT reference
  resolution, uncertainty propagation, and the production MESS handoff;
- an external stationary-TS/IRC/restart validation.

Until those gates pass, a run may request the `ANL1-F12` ladder head to select
the higher B2PLYP L2 profile, but it must fall back to the highest complete,
provenance-matched recipe. It must not rename a partial sum ANL1-F12.

## Bugs fixed for this validation

- Profiled inputs no longer stop in `QuantumChemistry` with an unimplemented
  router. A persistent router sends reaction/conformer work to L1 and L2
  refinements and hindered rotors to L2.
- The router records each job's backend and scheduler ID in
  `.kinbot_theory_jobs.json`, so result readers and restarts return to the
  correct backend.
- A preset backend change no longer inherits an empty or unrelated executable;
  Gaussian defaults to `g16` unless the profile supplies a command.
- Slurm's `queue_name` is emitted as `--partition`, rather than `-q` (QOS).
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
env -u LD_LIBRARY_PATH -u LD_PRELOAD git pull --ff-only origin composite
export KINBOT_CA_FILE=/etc/ssl/certs/ca-certificates.crt
export PIP_CERT="$KINBOT_CA_FILE"
PIP_CERT="$KINBOT_CA_FILE" .venv/bin/python -m pip install -e '.[fc]'
```

Run the local suite before spending licensed-code allocation:

```bash
cd ~/KinBot
mkdir -p .mpl-cache
MPLCONFIGDIR="$PWD/.mpl-cache" \
  .venv/bin/python -m pytest -q --ignore=tests/test_kinbot.py
```

`tests/test_kinbot.py` is a pre-existing collection/import harness rather than
the runnable unit suite. The expected branch result at the time this runbook
was written is 236 passed, 2 skipped, and 190 subtests passed.

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
advertise GPUs. Run it from a persistent login session such as `tmux` because
the driver watches and advances dependent jobs:

```bash
cd ~/KinBot
export KINBOT_PYTHON="$PWD/.venv/bin/python"
bash examples/anl/ethane_profiled_hpc/run.sh \
  day-long-cpu 3 "$KINBOT_FAIRCHEM_MODEL"
```

The first argument is the Slurm partition. The second is the maximum number of
simultaneous exclusive dispatcher nodes after the L3 geometry succeeds. The
optional third argument is the FairChem registered model name or local
checkpoint path; the local path is required for the documented offline HPC
run. The test directory defaults to `~/KinBot/ethane_profiled_hpc_run`;
override it with `KINBOT_PROFILED_TEST_DIR`.

Monitor either layer with:

```bash
cd ~/KinBot
squeue -u "$USER"
.venv/bin/python -m kinbot.anl.dispatch status \
  ethane_profiled_hpc_run/anl_interface
tail -f ethane_profiled_hpc_run/kinbot.log
```

The command is restartable. Rerun the same `run.sh` command after an
interruption. It reuses `kinbot.db`, the persistent L1/L2 route manifest, and
the dispatcher's immutable state. For a failed dispatcher task, inspect its
`execution.json`, `slurm.stderr`, and native output before retrying:

```bash
.venv/bin/python -m kinbot.anl.dispatch retry \
  ethane_profiled_hpc_run/anl_interface TASK_ID
.venv/bin/python -m kinbot.anl.dispatch drive \
  ethane_profiled_hpc_run/anl_interface --interval 20
```

Review the final machine-readable report:

```bash
cat ethane_profiled_hpc_run/anl_interface_audit.json
```

The report is acceptable for this validation only when every task is
`complete`, the native parsers are listed, a finite F12 CBS value is present,
and the top-level status remains `interface_complete_recipe_incomplete`.

## Remaining decisions before the production end-to-end test

1. Pin the exact equation and publication label for the desired F12/ANL1
   ladder head. If it is a Ram-style ladder rather than a canonical recipe,
   keep the `ANL-LADDER` provenance label.
2. Decide whether the first production reaction test remains ethane C-C
   homolysis or uses a stationary TS. The stationary-TS test is required
   before kinetics can be considered validated.
3. Enable MRCC for the CCSDTQ(P) increment, or explicitly define the highest
   allowed fallback when MRCC is unavailable.
4. Select the ATcT release/update policy and freeze a reviewed local snapshot
   with hashes for reproducible CBH results.
5. Provide the site MESS/MESSPF commands and compare the internal NASA7 fit
   against PAC99 before production thermochemistry is released.
