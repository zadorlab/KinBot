"""Actual outer + inner shell sequencing with simulated finite QC campaigns."""
import json
import shutil
import subprocess
import sys
from pathlib import Path

from tests.test_retained_run_handoff import handoff, RUNNER


DIAGNOSTICS = r'''
from types import ModuleType
screens = ModuleType('requested_orbital_recon')
relaxation = ModuleType('vrc_relaxation_recon')
def prepare_panel(production, plan, output):
    event('stage_relaxations')
    if os.environ.get('KB_TEST_BAD_STAGE'):
        raise ValueError('diagnostic identity mismatch')
    output.mkdir(parents=True, exist_ok=True)
    return dict(graphs=['same-graph'])
def prepare_screens(production, plan, output):
    assert not output.exists()
    event('stage_screens')
    output.mkdir()
    atomic_json(output/'review_launch.json', dict(systems=['same-system']))
    # Canonical JSON must compare equal to the differently formatted source.
    atomic_json(output/'plan.json', json.loads(plan.read_text()))
    return dict(tasks=24)
relaxation.prepare = prepare_panel
screens.batch = prepare_screens
screens._saved = lambda output: json.loads((output/'review_launch.json').read_text())
sys.modules['requested_orbital_recon'] = screens
sys.modules['vrc_relaxation_recon'] = relaxation
'''

SCRIPT_BRANCH = r'''
elif args[0].endswith('.py'):
    assert (root/'controller_refresh.json').exists()
    assert os.environ['PYTHONPATH'] == str(root/'.salvage_code.test')
    code = json.loads((root/'first/rotdPy/.kinbot_worker_code.json').read_text())
    assert code['code'] == str(root/'.salvage_code.test')
    name = Path(args[0]).stem + ':' + args[1]
    if args[1] == 'run':
        marker = Path(args[2])/'diagnostic.complete'
        if marker.exists():
            event('retained:' + name)
        else:
            marker.touch()
            event(name)
    else:
        print('diagnostic native evidence')
'''


def progress(handoff):
    root, decisions, retained, run = handoff
    checkout = root.parent/'checkout'
    source = Path(__file__).resolve().parents[1]
    script = checkout/'examples/anl/tools/refresh_progress_20261024.sh'
    shutil.copyfile(source/'examples/anl/tools/refresh_progress_20261024.sh', script)
    plan = checkout/'examples/anl/orbital_review_batches/r7_progress_20261024'
    plan.mkdir(parents=True)
    (plan/'screens.json').write_text('{"z":1,"a":2}')
    shutil.copyfile(decisions, plan/'decisions.json')
    python = checkout/'.venv/bin/python'
    runner = RUNNER.replace('args = sys.argv[1:]', DIAGNOSTICS+'\nargs = sys.argv[1:]')
    runner = runner.replace('else:\n    source = sys.stdin', SCRIPT_BRANCH+'else:\n    source = sys.stdin')
    python.write_text('#!'+sys.executable+'\n'+runner)
    git = root.parent/'bin/git'
    git.write_text('#!/bin/sh\nprintf "%s\\n" "${KB_TEST_REVISION:-reviewed-commit}"\n')
    git.chmod(0o755)
    def invoke(**environment):
        return subprocess.run(['bash', str(script), str(root)], env=dict(run.env, **environment),
                              capture_output=True, text=True, timeout=45)
    return root, retained, invoke


def test_progress_handoff_stages_before_one_replacement_and_retains_repeat(handoff):
    root, retained, invoke = progress(handoff)
    result = invoke()
    assert result.returncode == 0, result.stdout+result.stderr
    events = lambda: [json.loads(s) for s in (root/'events.jsonl').read_text().splitlines()]
    assert events() == ['stage_relaxations', 'stage_screens', 'refresh', 'decision_applied',
        'resume:absent', 'completed:complete', 'reconcile', 'monitor',
        'requested_orbital_recon:run', 'vrc_relaxation_recon:run']
    result = invoke()
    assert result.returncode == 0, result.stdout+result.stderr
    assert events().count('refresh') == events().count('stage_screens') == 1
    assert events()[-2:] == ['retained:requested_orbital_recon:run','retained:vrc_relaxation_recon:run']
    assert all(p.read_bytes() == old for p,old in retained.items())
    before = events()
    result = invoke(KB_TEST_REVISION='another-commit')
    assert result.returncode != 0 and 'another checkout revision' in result.stderr
    assert events() == before


def test_diagnostic_mismatch_prevents_controller_or_job_changes(handoff):
    root, retained, invoke = progress(handoff)
    result = invoke(KB_TEST_BAD_STAGE='1')
    assert result.returncode != 0 and 'diagnostic identity mismatch' in result.stderr
    assert not (root/'controller_refresh.json').exists()
    assert all(p.read_bytes() == old for p,old in retained.items())
