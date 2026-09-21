#!/usr/bin/env bash
set -euo pipefail

partition=${1:-day-long-cpu}
max_nodes=${2:-3}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd "$script_dir/../../.." && pwd)
run_dir=${KINBOT_PROFILED_TEST_DIR:-$repo_dir/ethane_profiled_hpc_run}
python_bin=${KINBOT_PYTHON:-$repo_dir/.venv/bin/python}
export HF_HUB_OFFLINE=${HF_HUB_OFFLINE:-1}

if [ ! -x "$python_bin" ]; then
    echo "KinBot Python is not executable: $python_bin" >&2
    exit 2
fi
for command in sbatch squeue sinfo g16 molpro xcfour; do
    command -v "$command" >/dev/null || {
        echo "Missing required command: $command" >&2
        exit 2
    }
done

mkdir -p "$run_dir"
"$python_bin" - "$script_dir/ethane.json" "$run_dir/ethane.json" "$partition" <<'PY'
import json
import sys
from pathlib import Path

source, target = Path(sys.argv[1]), Path(sys.argv[2])
partition = sys.argv[3]
data = json.loads(source.read_text())
data['queue_name'] = partition
target.write_text(json.dumps(data, indent=2) + '\n')
PY

"$python_bin" -c \
    'import ase, sella; from fairchem.core import FAIRChemCalculator; from kinbot.fairchem_utils import load_predictor; FAIRChemCalculator(load_predictor("uma-s-1p2", "cpu"), task_name="omol"); print("Python, ASE, Sella, and cached FairChem UMA model are ready")'

cd "$run_dir"
"$python_bin" -m kinbot.kb ethane.json

parent_job=301020900180000000001_well_high
if [ ! -f anl_interface/workflow.json ]; then
    "$python_bin" -m kinbot.anl.validation prepare-from-db \
        kinbot.db "$parent_job" anl_interface --max-nodes "$max_nodes" \
        --partition "$partition"
fi
"$python_bin" -m kinbot.anl.dispatch preflight anl_interface
"$python_bin" -m kinbot.anl.dispatch drive anl_interface --interval 20
"$python_bin" -m kinbot.anl.validation audit anl_interface | tee anl_interface_audit.json

echo "Validation completed. Review $run_dir/anl_interface_audit.json"
