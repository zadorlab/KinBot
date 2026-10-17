#!/usr/bin/env bash
# User-run diagnostic; all original results and other jobs are preserved.
set -euo pipefail
KB_PACKAGE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
KB_REPO=$(cd -- "$KB_PACKAGE/../../../.." && pwd)
cd "$KB_REPO"
export PYTHONPATH="$KB_REPO${PYTHONPATH:+:$PYTHONPATH}"
KB_PY=${KB_PY:-"$KB_REPO/.venv/bin/python"}
KB_SOURCE=${1:-/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_solver_retry}
KB_RUN=${2:-/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_molpro23_probe}
if [ -e "$KB_RUN" ]; then
    printf 'Diagnostic destination already exists; inspect it before resubmitting: %s\n' "$KB_RUN" >&2
    exit 1
fi

# Initialize the same module machinery as the source, then change only the
# Molpro module in this script's process. Do not purge unrelated modules.
export KINBOT_BACKEND=molpro
set +u
source "$KB_SOURCE/ch3ss_ss_atz_mrci_preconditioned/site_setup.sh"
KB_OLD_EXECUTABLE=$(command -v molpro)
if ! type module >/dev/null 2>&1; then
    printf 'The source setup did not initialize the module command.\n' >&2
    exit 1
fi
IFS=: read -r -a KB_LOADED_MODULES <<< "${LOADEDMODULES:-}"
for KB_LOADED_MODULE in "${KB_LOADED_MODULES[@]}"; do
    case "$KB_LOADED_MODULE" in molpro/*) module unload "$KB_LOADED_MODULE" ;; esac
done
# The source setup prepends an executable directory separately from the module.
# Remove that directory as well, or it can hide the newly loaded version.
PATH=$("$KB_PY" - "$PATH" "$KB_OLD_EXECUTABLE" <<'PY'
import os, sys
old = os.path.dirname(sys.argv[2])
print(os.pathsep.join(p for p in sys.argv[1].split(os.pathsep) if p != old))
PY
)
export PATH
module load molpro/molpro23
set -u
export KINBOT_MOLPRO_EXECUTABLE
KINBOT_MOLPRO_EXECUTABLE=$(command -v molpro)
KB_MODULE_PATH=$(module show molpro/molpro23 2>&1)
export KB_MODULE_PATH
"$KB_PY" - "$KB_SOURCE" "$KB_RUN" "$KB_PACKAGE" <<'PY'
import json, os, re, subprocess, sys
from pathlib import Path
from kinbot.execution import atomic_json, file_sha256

source, root, package = map(Path, sys.argv[1:])
sys.path.insert(0, str(package))
from prepare_recovery import prepare
from run_confirmation import validated_record

executable = os.environ['KINBOT_MOLPRO_EXECUTABLE']
version = subprocess.run([executable, '--version'], text=True, capture_output=True,
                         check=True, timeout=30)
version_text = version.stdout + version.stderr
if not re.search(r'\b2023[.]\d+', version_text):
    raise RuntimeError('Expected a 2023 Molpro build; got: ' + version_text)
record = prepare(source, root, package/'mrci_build_probe_plan.json', None)
old = validated_record(source)
old_system = next(s for s in old['systems'] if s['title'] == record['systems'][0]['title'])
old_manifest = json.loads(Path(old_system['manifest']).read_text())
prior = next(j for j in old_manifest['jobs'] if j['anchor'] == 'center')
system = record['systems'][0]
manifest = json.loads(Path(system['manifest']).read_text())
job, = manifest['jobs']
if job['input_sha256'] != prior['input_sha256']:
    raise RuntimeError('Build probe input differs from the failed trial; no job submitted.')
if job['initial_record']['sha256'] != prior['initial_record']['sha256']:
    raise RuntimeError('Build probe starting checkpoint differs; no job submitted.')
if prior.get('retry_solver_options', old_manifest.get('solver_options', {})) != manifest['solver_options']:
    raise RuntimeError('Source has different effective solver overrides; no job submitted.')
if any(prior.get(key) is not None for key in
       ('retry_max_reference_vectors', 'retry_max_reference_iterations', 'recovered_scf')):
    raise RuntimeError('Source has additional retry transformations; review before comparing builds.')
native = Path(prior['directory'])/(prior['name']+'.native.inp')
execution = json.loads((Path(prior['directory'])/'execution.json').read_text())
native_digest = file_sha256(native)
if native_digest != execution.get('native_input_sha256'):
    raise RuntimeError('Failed native input no longer matches its execution record.')
# Test the generated batch setup itself, not merely the login-shell binary.
setup = Path(system['manifest']).parent/'site_setup.sh'
probe = subprocess.run(['bash', '-c',
    'set -eo pipefail; export KINBOT_BACKEND=molpro; source "$1"; molpro --version',
    'probe', str(setup)], text=True, capture_output=True, check=True, timeout=30)
if not re.search(r'\b2023[.]\d+', probe.stdout + probe.stderr):
    raise RuntimeError('Generated site setup does not select Molpro 2023; no job submitted.')
atomic_json(root/'build_comparison.json', {
    'schema': 1, 'production_approved': False, 'source_directory': str(source),
    'source_job_id': prior.get('job_id'), 'module': 'molpro/molpro23',
    'executable': executable, 'version': version_text.strip(),
    'site_setup_sha256': file_sha256(setup), 'input_sha256': job['input_sha256'],
    'initial_checkpoint_sha256': job['initial_record']['sha256'],
    'source_native_input_sha256': native_digest,
    'input_and_starting_checkpoint_identical': True,
    'scope': 'one center point; both six-root spin blocks; no tilted submission',
    'module_show': os.environ['KB_MODULE_PATH']})
print(version_text.strip())
print('Identical electronic input and seed verified:', root)
PY
nohup "$KB_PY" -u "$KB_PACKAGE/run_confirmation.py" run "$KB_RUN" \
    >> "$KB_RUN/controller.log" 2>&1 < /dev/null &
KB_CONTROLLER_PID=$!
sleep 3
if ! kill -0 "$KB_CONTROLLER_PID" 2>/dev/null; then
    tail -n 60 "$KB_RUN/controller.log"
    exit 1
fi
printf 'Build comparison only: %s\nController PID: %s\n' "$KB_RUN" "$KB_CONTROLLER_PID"
"$KB_PY" "$KB_PACKAGE/run_confirmation.py" monitor "$KB_RUN" --all --roots
