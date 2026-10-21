"""Build comparison must preserve the failed physical problem and evidence."""
import tempfile
import json
from pathlib import Path
import re
import os
import subprocess
import sys
from types import SimpleNamespace

from kinbot import molpro_reconnaissance as recon
from kinbot.execution import file_sha256
from test_molpro_confirmation_package import PACKAGE
from test_mrci_preconditioner import stopped_multistate


def test_build_probe_uses_identical_deck_and_pristine_seed(tmp_path, monkeypatch):
    adapter, multistate, retry = stopped_multistate(tmp_path, monkeypatch)
    record = adapter.prepare(multistate, retry, PACKAGE/'mrci_preconditioner_plan.json', None)
    manifest_path = Path(record['systems'][0]['manifest'])
    manifest = json.loads(manifest_path.read_text())
    for item in manifest['jobs']:
        item['status'] = 'prepared'
    prior = manifest['jobs'][0]
    directory = Path(prior['directory'])
    prior['status'] = 'incomplete'
    native = directory/(prior['name']+'.native.inp')
    native.write_text('Synthetic failed native-input provenance fixture')
    (directory/'execution.json').write_text(json.dumps({
        'status': 'failed', 'returncode': 255, 'native_input_sha256': file_sha256(native)}))
    (directory/(prior['name']+'.out')).write_text('? dsygv failed\n? The problem occurs in cif12solve\n')
    (directory/(prior['name']+'.wfu.failed')).write_bytes(b'Do not restart from this failed output')
    manifest_path.write_text(json.dumps(manifest))
    record['status'] = 'failed'
    (retry/'confirmation.json').write_text(json.dumps(record))
    before = {p: file_sha256(p) for p in retry.rglob('*') if p.is_file()}
    target = tmp_path/'build_probe'
    result = adapter.prepare(retry, target, PACKAGE/'mrci_build_probe_plan.json', None)
    assert result['execution']['max_jobs'] == 1
    assert len(result['systems']) == 1
    job, = json.loads(Path(result['systems'][0]['manifest']).read_text())['jobs']
    assert job['anchor'] == 'center'
    assert job['input_sha256'] == prior['input_sha256']
    assert job['initial_record']['sha256'] == prior['initial_record']['sha256']
    assert Path(job['initial_record']['path']).read_bytes() == Path(prior['initial_record']['path']).read_bytes()
    assert all(file_sha256(p) == digest for p, digest in before.items())
    # Execute the launcher's Python preflight with synthetic version responses;
    # it must bind the already-tested plan to the failed native input as well.
    script = (PACKAGE/'probe_mrci_build_on_blodgett.sh').read_text()
    code = script.split('"$KB_PY" - "$KB_SOURCE" "$KB_RUN" "$KB_PACKAGE" <<\'PY\'\n', 1)[1].split('\nPY\n', 1)[0]
    with monkeypatch.context() as context:
        context.setitem(sys.modules, 'prepare_recovery', adapter)
        context.setattr(sys, 'argv', ['probe', str(retry), str(target), str(PACKAGE)])
        context.setattr(sys, 'path', list(sys.path))
        context.setenv('KINBOT_MOLPRO_EXECUTABLE', '/synthetic/molpro23')
        context.setenv('KB_MODULE_PATH', 'synthetic module definition')
        context.setattr(subprocess, 'run', lambda *a, **k: SimpleNamespace(
            stdout='Molpro 2023.2\n', stderr='', returncode=0))
        exec(compile(code, 'build-probe-preflight', 'exec'), {})
    proof = json.loads((target/'build_comparison.json').read_text())
    assert proof['source_native_input_sha256'] == file_sha256(native)
    assert proof['input_and_starting_checkpoint_identical'] is True
    assert proof['production_approved'] is False
    assert all(file_sha256(p) == digest for p, digest in before.items())
    monkeypatch.setattr(recon, 'qc_runtime_environment',
                        lambda *a, **k: ({}, {'scratch_directory': tempfile.mkdtemp(dir=tmp_path)}))
    monkeypatch.setattr(recon, 'cleanup_qc_runtime', lambda *a: None)
    native_decks = []
    def launch(command, **kwargs):
        native_decks.append((kwargs['cwd']/command[-1]).read_text())
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', launch)
    for item in (prior, job):
        item['retry_solver_options'] = manifest['solver_options']
        recon._execute(item)
    assert native_decks[0] == native_decks[1]
    assert native_decks[0].count('state,6,1,2,3,4,5,6;') == 4  # STATE + REFSTATE, both spins
    assert len(re.findall(r'^reconc\d+=', native_decks[0], re.M)) == 12


def test_module_switch_removes_pinned_old_path(tmp_path):
    old, new = tmp_path/'old', tmp_path/'new'
    old.mkdir()
    new.mkdir()
    for directory, version in ((old, '2024.1'), (new, '2023.2')):
        p = directory/'molpro'
        p.write_text(f'#!/bin/bash\nprintf "Molpro {version}\\n"\n')
        p.chmod(0o700)
    source = tmp_path/'source'
    setup = source/'ch3ss_ss_atz_mrci_preconditioned/site_setup.sh'
    setup.parent.mkdir(parents=True)
    # Emulate both module-inserted and separately pinned old PATH entries.
    setup.write_text('''export PATH="$OLD_BIN:$OLD_BIN:$PATH"
export LOADEDMODULES=molpro/molpro24
module() {
    case "$1" in
        unload) export LOADEDMODULES= ;;
        load) export PATH="$NEW_BIN:$PATH"; export LOADEDMODULES=molpro/molpro23 ;;
        show) printf 'module definition\\n' ;;
        *) return 1 ;;
    esac
}
''')
    script = (PACKAGE/'probe_mrci_build_on_blodgett.sh').read_text()
    start = script.index('export KINBOT_BACKEND=molpro')
    end = script.index('"$KB_PY" - "$KB_SOURCE" "$KB_RUN" "$KB_PACKAGE"')
    code = 'set -euo pipefail\n' + script[start:end]
    code += '\nprintf "%s\\n%s\\n" "$KINBOT_MOLPRO_EXECUTABLE" "$PATH"\n'
    environment = {**os.environ, 'KB_SOURCE': str(source), 'KB_PY': sys.executable,
                   'OLD_BIN': str(old), 'NEW_BIN': str(new)}
    result = subprocess.run(['bash', '-c', code], env=environment,
                            capture_output=True, text=True, check=True)
    executable, search_path = result.stdout.strip().splitlines()
    assert executable == str(new/'molpro')
    assert str(old) not in search_path.split(os.pathsep)
