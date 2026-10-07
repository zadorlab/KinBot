#!/usr/bin/env bash
set -u

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd "$script_dir/../../.." && pwd)
evidence_dir=${KINBOT_PROFILED_TEST_DIR:-$repo_dir/ethane_profiled_hpc_run_v5}
run_dir=${KINBOT_VRC_PRODUCTION_DIR:-$evidence_dir/vrc_production_run}
python_bin=${KINBOT_PYTHON:-$repo_dir/.venv/bin/python}
reaction=301020900180000000001_hom_sci_1_2

echo "===== time ====="
date
echo "===== queue ====="
squeue -u "$USER" -o "%.18i %.32j %.2t %.10M %.10l %R" 2>&1

echo "===== production driver ====="
pid_file="$repo_dir/ethane_vrc_production.pid"
if test -f "$pid_file"; then
    pid=$(cat "$pid_file")
    if kill -0 "$pid" 2>/dev/null; then
        echo "running: PID $pid"
    else
        echo "exited: PID $pid"
    fi
else
    echo "PID file not found: $pid_file"
fi

pointer="$run_dir/vrctst/molpro/current_dispatch.json"
if test -f "$pointer"; then
    dispatch=$($python_bin - "$pointer" <<'PY'
import json
import sys
from pathlib import Path
print(json.loads(Path(sys.argv[1]).read_text())['run_dir'])
PY
)
    echo "===== correction dispatcher: $dispatch ====="
    "$python_bin" -m kinbot.anl.dispatch status "$dispatch" 2>&1
else
    echo "===== correction dispatcher ====="
    echo "not prepared"
fi

echo "===== correction record ====="
if test -f "$run_dir/vrctst/corr_$reaction.json"; then
    "$python_bin" - "$run_dir/vrctst/corr_$reaction.json" <<'PY'
import json
import sys
from pathlib import Path
record = json.loads(Path(sys.argv[1]).read_text())
print(json.dumps({
    'point_count': len(record.get('dist', [])),
    'levels': record.get('levels'),
}, indent=2))
PY
else
    echo "not complete"
fi

echo "===== ROTD_py ====="
manifest="$run_dir/rotdPy/$reaction.rotdpy.json"
input="$run_dir/rotdPy/$reaction.py"
if test -f "$manifest"; then
    "$python_bin" -m kinbot.rotdpy check "$input" --profile production 2>&1
elif test -f "$input"; then
    find "$run_dir/rotdPy/kb_$reaction" -maxdepth 2 -type f \
        \( -name 'surface_*.dat' -o -name 'Ne_*.out' \) -print 2>/dev/null \
        | sort | tail -n 40
    tail -n 40 "$run_dir/rotdPy/$reaction.rotdpy.stderr" 2>/dev/null
else
    echo "input not generated"
fi

echo "===== recent worker errors ====="
find "$run_dir/perm" -maxdepth 1 -type f -name '*.err' -size +0c \
    -print 2>/dev/null | sort | tail -n 4 | while read -r error_file; do
        echo "--- $error_file ---"
        tail -n 25 "$error_file"
    done

echo "===== latest KinBot messages ====="
tail -n 35 "$run_dir/kinbot.log" 2>/dev/null
