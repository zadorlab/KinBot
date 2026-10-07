#!/usr/bin/env bash
set -euo pipefail

partition=${1:-day-long-cpu}
max_nodes=${2:-8}
fairchem_model=${3:-${KINBOT_FAIRCHEM_MODEL:-uma-s-1p2}}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd "$script_dir/../../.." && pwd)
evidence_dir=${KINBOT_PROFILED_TEST_DIR:-$repo_dir/ethane_profiled_hpc_run_v5}
# PR 108 calculation directories carry an immutable format marker.  The
# completed profiled smoke run predates that format and must remain read-only;
# a production VRC continuation therefore gets its own current-format run.
run_dir=${KINBOT_VRC_PRODUCTION_DIR:-$evidence_dir/vrc_production_run}
python_bin=${KINBOT_PYTHON:-$repo_dir/.venv/bin/python}
reaction=301020900180000000001_hom_sci_1_2

test -x "$python_bin" || {
    echo "KinBot Python is not executable: $python_bin" >&2
    exit 2
}
test -d "$evidence_dir" || {
    echo "Profiled validation directory is missing: $evidence_dir" >&2
    exit 2
}
for command in sbatch squeue sinfo g16 molpro; do
    command -v "$command" >/dev/null || {
        echo "Missing required command: $command" >&2
        exit 2
    }
done

# Fail before reaction work if a compute-node dependency or the private
# sampler is unavailable.  Loading the predictor here also proves that the
# supplied path names a usable checkpoint while the login-node environment is
# still visible.
"$python_bin" - "$fairchem_model" <<'PY'
import sys

from fairchem.core import FAIRChemCalculator
from kinbot.fairchem_utils import load_predictor
from kinbot.rotdpy import ensure_available

model = sys.argv[1]
FAIRChemCalculator(load_predictor(model, 'cpu'), task_name='omol')
rotdpy = ensure_available()
print(f'FairChem model {model!r} and rotdPy are ready')
print(f"rotdPy revision: {rotdpy.get('git_revision')}")
PY

mkdir -p "$run_dir"
"$python_bin" - "$script_dir/ethane_vrc_production.json" \
    "$run_dir/ethane_vrc_production.json" "$partition" \
    "$fairchem_model" "$max_nodes" <<'PY'
import json
import sys
from pathlib import Path

source, target = Path(sys.argv[1]), Path(sys.argv[2])
data = json.loads(source.read_text())
data['queue_name'] = sys.argv[3]
data['queue_template'] = str(
    source.resolve().parents[3] / 'kinbot' / 'tpl' / 'slurm_partition.tpl')
data['fc_model_path'] = sys.argv[4]
data['vrc_tst_max_nodes'] = int(sys.argv[5])
data['rotdpy_max_jobs'] = int(sys.argv[5])
target.write_text(json.dumps(data, indent=2) + '\n')
PY

# A reduced result in this child directory has a restart database containing
# a different set of surfaces and no correction potential. Preserve it, but
# never mix those samples into this calculation.
if test -d "$run_dir/rotdPy"; then
    profile=interface
    if test -f "$run_dir/rotdPy/$reaction.rotdpy.json"; then
        profile=$("$python_bin" - "$run_dir/rotdPy/$reaction.rotdpy.json" <<'PY'
import json
import sys
from pathlib import Path
record = json.loads(Path(sys.argv[1]).read_text())
print((record.get('validation') or {}).get('name', 'interface'))
PY
        )
    elif test -f "$run_dir/rotdPy/$reaction.py" \
        && grep -q "'name': 'production'" \
            "$run_dir/rotdPy/$reaction.py"; then
        # Preserve an interrupted production database so ROTD_py can restart.
        profile=production
    fi
    if test "$profile" != production; then
        stamp=$(date +%Y%m%dT%H%M%S)
        mv "$run_dir/rotdPy" "$run_dir/rotdPy_interface_$stamp"
    fi
fi

cd "$run_dir"
"$python_bin" -m kinbot.kb ethane_vrc_production.json
"$python_bin" -m kinbot.rotdpy check \
    "rotdPy/$reaction.py" --profile production \
    | tee rotdpy_production_gate.json

echo "Production ROTD_py result passed: $run_dir/rotdpy_production_gate.json"
