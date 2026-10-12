#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd "$script_dir/../../.." && pwd)
python_bin=${KINBOT_PYTHON:-$repo_dir/.venv/bin/python}
base=${KINBOT_PROFILED_TEST_DIR:-$repo_dir/ethane_profiled_hpc_run_v5}
partition=${1:-day-long-cpu}
max_jobs=${2:-8}
ppn=${3:-8}
manifest=${KINBOT_FULL_VERIFICATION_MANIFEST:-$base/full_verification_manifest.json}

test -x "$python_bin" || {
    echo "KinBot Python is not executable: $python_bin" >&2
    exit 2
}

"$python_bin" "$script_dir/full_verification_manifest.py" \
    "$base" --output "$manifest" --partition "$partition" --ppn "$ppn" \
    --max-jobs "$max_jobs"

# Preparation reparses every native QC record, assembles each profiled
# composite, solves both CBH reactions, validates the production ROTD result,
# and writes the strict two-species handoff before a MESS job is submitted.
"$python_bin" -m kinbot.anl.end_to_end run "$manifest" \
    | tee "$base/full_verification.log"

echo "Full verification passed: $base/full_verification_report.json"
