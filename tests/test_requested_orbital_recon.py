"""Exercise the exact new-well screen recipe without submitting native work."""
import importlib.util
import json
from pathlib import Path

import pytest

from kinbot.execution import atomic_json

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('requested_orbital_recon',
    ROOT/'examples/anl/tools/requested_orbital_recon.py')
helper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(helper)
RECIPE = ROOT/'examples/anl/ch3s_anl1_f12/orbital_reviews/ch2sh_cs_screen.json'


@pytest.fixture
def request_tree(tmp_path):
    pool = tmp_path/'production'
    case = pool/'case'
    directory = case/'471030810040000000002/vrctst/orbital_review/reaction'
    directory.mkdir(parents=True)
    request = directory/'request.json'
    recipe = json.loads(RECIPE.read_text())
    data = dict(recipe['request_match'], schema=1, status='needs_review',
                smiles='[H][S]=[C]([H])[H]', positions=[
        [1.03440581452718, -.3719293385761376, .03988588054708078],
        [.3381038973372055, .7267766509223542, 1.1679876585701665],
        [.7607169349022845, 1.8652195195547316, .6015953344491213],
        [.8725603978962576, -1.4191794072599508, .22000490140932294],
        [1.6514884746631548, -.030599788522718163, -.7700430150262825]])
    atomic_json(request, data)
    atomic_json(case/'input.json', {'queue_name': 'test-partition'})
    atomic_json(pool/'.kinbot_resource_policy.json', dict(schema=1, max_jobs=23,
        excluded_nodes=['badnode'], manage_pending=True, orbital_review_reserve=3))
    return request, case/'input.json', tmp_path/'screen'


def test_exact_request_prepares_eight_independent_probes_and_no_production_writes(request_tree, monkeypatch):
    request, case_input, output = request_tree
    original = {p: p.read_bytes() for p in case_input.parent.parent.rglob('*') if p.is_file()}
    monkeypatch.setattr(helper.subprocess, 'run', lambda *a, **k: pytest.fail('No launch during preparation'))
    result = helper.prepare(request, case_input, RECIPE, output)
    assert result['tasks'] == 8 and result['max_jobs'] == 3
    assert result['active_space']['active_electrons'] == 13
    assert result['active_space']['active_orbitals'] == 11
    manifest = json.loads((output/'native/reconnaissance.json').read_text())
    assert manifest['slurm_excluded_nodes'] == ['badnode']
    for job in manifest['jobs']:
        assert job['orbital_predecessor'] is None
        assert [(s['multiplicity'], s['root']) for s in job['root_states']] == [
            (4, 1), (4, 2), (2, 1), (2, 2), (2, 3), (2, 4)]
        d = Path(job['directory'])
        deck = (d/(job['name']+'.inp')).read_text()
        assert '{rhf,so-sci,maxit=256;wf,25,1,3;occ,14;closed,11;' in deck
        assert '{multi,cas-ci;' in deck and 'avas,nela=13' in deck
        assert 'center,5,1s' in deck and 'start,2100.2' in deck
        assert 'wf,25,1,1\nstate,4' in deck
        assert 'rs2' not in deck.lower()
        script = (d/(job['name']+'.sbatch')).read_text()
        assert '#SBATCH --exclusive' in script and '#SBATCH -p test-partition' in script
        assert '# KINBOT_WALLTIME=partition_max' in script
    assert all(p.read_bytes() == value for p, value in original.items())
    with pytest.raises(FileExistsError):
        helper.prepare(request, case_input, RECIPE, output)


@pytest.mark.parametrize('field,value', [('parent', 'wrong'), ('multiplicity', 4),
                                        ('status', 'model_accepted')])
def test_wrong_or_accepted_request_cannot_be_prepared(request_tree, field, value):
    request, case_input, output = request_tree
    data = json.loads(request.read_text())
    data[field] = value
    atomic_json(request, data)
    with pytest.raises(ValueError):
        helper.prepare(request, case_input, RECIPE, output)
    assert not output.exists()


def test_refuse_to_write_into_production_tree(request_tree):
    request, case_input, _ = request_tree
    with pytest.raises(ValueError, match='outside the production'):
        helper.prepare(request, case_input, RECIPE, case_input.parent/'screen')


def test_launch_uses_current_exclusions_reserved_cap_and_original_manifest(request_tree, monkeypatch):
    request, case_input, output = request_tree
    helper.prepare(request, case_input, RECIPE, output)
    policy_file = case_input.parent.parent/'.kinbot_resource_policy.json'
    policy = json.loads(policy_file.read_text())
    policy.update(orbital_review_reserve=2, excluded_nodes=['badnode', 'newbadnode'])
    atomic_json(policy_file, policy)
    monkeypatch.setattr(helper.subprocess, 'run', lambda *a, **k: None)
    from kinbot import run_walltime
    monkeypatch.setattr(run_walltime, 'partition_walltime', lambda _: '31-00:00:00')
    calls = []
    def launch(paths, **kwargs):
        calls.append((paths, kwargs))
        state = json.loads(paths[0].read_text())
        assert state['slurm_excluded_nodes'] == ['badnode', 'newbadnode']
        return [dict(status='complete_review_required', jobs=state['jobs'])]
    monkeypatch.setattr(helper.recon, 'run_many', launch)
    monkeypatch.chdir(output)  # Restore cwd automatically after this test.
    result = helper.run(output)
    assert len(calls) == 1 and calls[0][1] == dict(max_jobs=2, allocation_retries=2)
    assert not result['production_model_accepted']


def test_screen_monitor_reads_registry_without_pilot_or_manifest_writes(request_tree, monkeypatch):
    request, case_input, output = request_tree
    helper.prepare(request, case_input, RECIPE, output)
    from kinbot import orbital_monitor
    monkeypatch.setattr(orbital_monitor, '_queue', lambda: ({}, None))
    before = {p: p.read_bytes() for p in output.rglob('*') if p.is_file()}
    data = orbital_monitor.snapshot(output, registry=output/'review_launch.json', include_all=True)
    assert data['counts']['manifest']['WAITING'] == 8
    assert '8/8' not in orbital_monitor.render(data, details=True)
    assert all(p.read_bytes() == value for p, value in before.items())
