# Full ethane / methyl CBH-ANL, ROTD_py, MESSPF, and MESS validation

This runbook is the current external-site acceptance procedure. Historical
smoke-test transcripts and superseded recovery commands are retained in Git
history rather than repeated here.

## What the completed calculation must prove

The gate reuses the already completed electronic-structure and production
ROTD_py calculations. It verifies:

- FairChem UMA L1 search and B2PLYP-D3(BJ)/cc-pVTZ L2 refinement;
- one accepted ethane L2 hindered rotor;
- ANL0-F12 component assembly for ethane, methane, H2, CH3, and H;
- published source-geometry higher-order values for CH4 and CH3;
- ethane and methyl CBH-0 Hf(0 K);
- the 24-surface production ROTD_py result from
  CASPT2(2e,2o)/cc-pVDZ sampling and MRCI+Q(2e,2o)/cc-pVTZ corrections;
- MESSPF partition functions for every live and CBH-reference species;
- CBH-reaction Hf(298.15 K), NASA-7 anchors, and continuity;
- a MESS `Core Rotd` source hash matching the selected ROTD_py result;
- the CH3 + CH3 high-pressure rate table against the configured literature
  envelope.

The final pass must not submit ANL single points, VRC correction points, or
ROTD_py samples again. It may repeat only missing L2 or hindered-rotor work
needed to reconstruct the current KinBot/MESS network.

## Scientific contracts

For every signed CBH reaction, `nu_target = -1`:

```text
Delta_r E(0) = sum_i nu_i E_i(0)
Hf_target(0) = sum_(j != target) nu_j Hf_j(0) - Delta_r E(0)
```

MESSPF supplies:

```text
H_i(T)-H_i(0) = R T^2 d ln(Q_i)/dT + R T
```

and the finite-temperature solve uses the same CBH equation:

```text
Delta_r H(T) = Delta_r E(0) + sum_i nu_i [H_i(T)-H_i(0)]
Hf_target(T) = sum_(j != target) nu_j Hf_j(T) - Delta_r H(T)
```

Only CBH-0 may read direct pinned ATcT values. Higher rungs recursively use
computed lower-rung formation enthalpies at both 0 and 298.15 K. Elemental
thermal increments are not used.

## Prerequisites

The checkout must contain:

- the completed `ethane_profiled_hpc_run_v5` calculation tree;
- the initialized private `external/ROTD_py` submodule;
- the cached FairChem model used by the completed search;
- Gaussian, Molpro, CFOUR, and MRCC available through site modules or user
  configuration;
- working `mess` and `messpf` executables.

The final gate itself does not download models or rerun completed licensed
work.

## Update the checkout

```bash
cd ~/KinBot

env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  git -c submodule.recurse=false pull --ff-only origin composite

.venv/bin/python -m pip install -e . \
  --no-build-isolation --no-deps

.venv/bin/python -m pytest -q \
  tests/test_anl_cbh_atct.py \
  tests/test_anl_cbh_table_handoff.py \
  tests/test_anl_mess_thermochemistry.py \
  tests/test_anl_end_to_end.py
```

Do not clean or reset the completed calculation directory.

## Verify MESS

```bash
cd ~/KinBot

export KINBOT_MESS_COMMAND="$HOME/.local/kinbot-mess/bin/mess"
export KINBOT_MESSPF_COMMAND="$HOME/.local/kinbot-mess/bin/messpf"

test -x "$KINBOT_MESS_COMMAND" || {
  echo "MESS executable missing"
  exit 1
}

test -x "$KINBOT_MESSPF_COMMAND" || {
  echo "MESSPF executable missing"
  exit 1
}
```

These variables accept another site installation without changing KinBot.

## Prepare the immutable inputs

```bash
cd ~/KinBot

base="$PWD/ethane_profiled_hpc_run_v5"
manifest="$base/full_verification_manifest.json"

export KINBOT_MESS_COMMAND="$HOME/.local/kinbot-mess/bin/mess"
export KINBOT_MESSPF_COMMAND="$HOME/.local/kinbot-mess/bin/messpf"

.venv/bin/python \
  examples/anl/ethane_profiled_hpc/full_verification_manifest.py \
  "$base" \
  --output "$manifest" \
  --partition day-long-cpu \
  --ppn 8 \
  --max-jobs 8

set -o pipefail

.venv/bin/python -m kinbot.anl.end_to_end prepare "$manifest" \
  | tee "$base/full_verification_preparation.console.json"
```

Preparation must finish with `"status": "prepared"`. It performs no Slurm
submission. It reparses and hash-checks native QC data, assembles the
composites, solves CBH-0, verifies literature checks, writes the strict
handoff, and writes `cbh_thermochemistry_models.json` for CH4, H2, and H.

Inspect the compact preparation result:

```bash
.venv/bin/python - <<'PY'
import json
from pathlib import Path

base = Path('ethane_profiled_hpc_run_v5')
record = json.loads((base / 'full_verification_preparation.json').read_text())
print('status:', record['status'])
print('ROTD surfaces:', record['rotdpy']['surface_count'])
print('handoff species:', record['handoff_species'])
print('CBH reference models:', record['cbh_thermochemistry_models'])
for name, entry in record['cbh'].items():
    formation = entry['formation']
    print(name, formation['method'], formation['formation_0k_kcal_mol'])
PY
```

## Run MESSPF and MESS

```bash
cd ~/KinBot

base="$PWD/ethane_profiled_hpc_run_v5"
manifest="$base/full_verification_manifest.json"

export KINBOT_MESS_COMMAND="$HOME/.local/kinbot-mess/bin/mess"
export KINBOT_MESSPF_COMMAND="$HOME/.local/kinbot-mess/bin/messpf"

nohup .venv/bin/python -m kinbot.anl.end_to_end run "$manifest" \
  > "$base/full_verification.console.log" 2>&1 < /dev/null &

echo $! > ethane_full_verification.pid
echo "driver PID: $(cat ethane_full_verification.pid)"
```

The existing production ROTD_py files are inputs to this command. No
`rotd-*` jobs should be submitted.

## Short monitoring command

```bash
cd ~/KinBot

base="$PWD/ethane_profiled_hpc_run_v5"
pid=$(cat ethane_full_verification.pid)

kill -0 "$pid" 2>/dev/null \
  && echo "driver running: $pid" \
  || echo "driver exited"

squeue -h -u "$USER" -o '%.2t %.8M %.30j %R' \
  | grep -E 'hir|well_high|mess|kb-' || true

tail -n 25 "$base/full_verification.console.log"
```

If `rotd-*` jobs appear, stop the run and inspect the generated KinBot input;
this continuation must use `rotdpy_run=0` and the completed production result.

## Audit after completion

```bash
cd ~/KinBot

base="$PWD/ethane_profiled_hpc_run_v5"
manifest="$base/full_verification_manifest.json"

.venv/bin/python -m kinbot.anl.end_to_end audit "$manifest" \
  | tee "$base/full_verification_audit.console.json"
```

The final file is:

```text
ethane_profiled_hpc_run_v5/full_verification_report.json
```

The top-level status must be `passed`. The report must include:

- exact native and derived record hashes;
- the selected ANL recipe and signed component contributions;
- CBH reactions and recursive reference provenance;
- MESSPF input/output hashes for network and CBH-reference species;
- independently reconstructed Hf(298.15 K) CBH traces;
- NASA-7 H/S anchors and range continuity;
- the selected ROTD_py number-of-states hash;
- MESS high-pressure rates and literature comparison.

## Failure diagnostics

Preparation failures are deterministic and should be fixed before submitting
anything. For a MESSPF failure, inspect:

```bash
run="$PWD/ethane_profiled_hpc_run_v5/vrc_production_run"
find "$run/me/partition_functions" -maxdepth 1 \
  -type f \( -name '*.log' -o -name '*.stderr' \) -size +0 -print
```

For a KinBot continuation failure:

```bash
tail -n 80 ethane_profiled_hpc_run_v5/vrc_production_run/kinbot.log
tail -n 80 ethane_profiled_hpc_run_v5/full_verification.console.log
```

Do not delete native outputs. The dispatcher and audits rely on their hashes,
and parser fixes should reparse them without repeating calculations.
