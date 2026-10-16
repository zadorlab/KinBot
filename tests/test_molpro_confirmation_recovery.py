"""Recovery provenance and rendering tests; fixtures are NOT native evidence."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from kinbot.execution import file_sha256
from test_molpro_confirmation_package import module, build_fixture, PACKAGE


def fixture(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(PACKAGE))
    adapter = module('prepare_recovery')
    source = build_fixture(tmp_path).parent
    record = json.loads((source/'confirmation.json').read_text())
    reviewed = []
    unresolved = {('ch3s_cs_atz_mrci', 'tilted'), ('ch3ss_ss_atz_mrci', 'far'),
                  ('ch3ss_ss_atz_mrci', 'center'), ('ch3ss_ss_atz_mrci', 'tilted')}
    for s in record['systems']:
        p = Path(s['manifest'])
        m = json.loads(p.read_text())
        rows = []
        for j in m['jobs']:
            d = Path(j['directory'])
            j['status'] = 'incomplete' if (s['title'], j['anchor']) in unresolved else 'complete'
            row = {'name': j['name'], 'status': j['status'], 'input_sha256': j['input_sha256'],
                   'output_sha256': None, 'wfu_sha256': None}
            if j['status'] == 'complete':
                out = d/(j['name']+'.out')
                out.write_text('Synthetic completed output fixture')
                wfu = d/(j['name']+'.wfu')
                wfu.write_bytes(b'Synthetic selected-CAS checkpoint fixture')
                j.update(wfu=str(wfu), wfu_sha256=file_sha256(wfu))
                (d/'execution.json').write_text(json.dumps({'status':'complete','returncode':0}))
                row.update(output_sha256=file_sha256(out), wfu_sha256=file_sha256(wfu))
            rows.append(row)
        p.write_text(json.dumps(m))
        reviewed.append({'title': s['title'], 'jobs': rows})
    review = tmp_path/'review.json'
    review.write_text(json.dumps({'systems': reviewed}))
    monkeypatch.setattr(adapter.run_cleanup, 'inspect', lambda root: {'jobs': [], 'processes': []})
    return adapter, source, tmp_path/'recovery', PACKAGE/'recovery_plan.json', review


def test_recovery_preserves_source_and_renders_requested_changes(tmp_path, monkeypatch):
    adapter, source, target, plan, review = fixture(tmp_path, monkeypatch)
    before = {p: file_sha256(p) for p in source.rglob('*') if p.is_file()}
    record = adapter.prepare(source, target, plan, review)
    assert sum(s['jobs'] for s in record['systems']) == 26
    assert len(record['systems']) == 7 and record['execution']['memory_mb'] == 104000
    for system in record['systems']:
        m = json.loads(Path(system['manifest']).read_text())
        cfg = json.loads(Path(system['configuration']).read_text())
        assert cfg['probe_definition']['state_model']['correlation_options']['reference_policy'] == 'state_average'
        for i, j in enumerate(m['jobs']):
            text = (Path(j['directory'])/(j['name']+'.inp')).read_text()
            if 'mrci' in system['title']:
                assert 'option,cluster=0' in text and 'KBVRCSEL=energd(1)' in text
            else:
                assert 'shift=0.3,ipea=0.25' in text
                assert 'refstate,' in text and 'weight,' in text
                assert 'KBVRCSEL=energy(1)' in text and 'cabs_singles=' not in text
            assert 'maxvi=60' in text and 'thrvar=1.d-10' in text
            if 'maxdav20' in system['title']:
                assert 'maxdav=20' in text and 'cabs_singles=1' in text
                assert cfg['probe_definition']['orbital_partition']['close_lowest'] == 0
                assert cfg['probe_definition']['orbital_partition']['orbital_swaps'] == []
            elif 'mrci' in system['title']:
                assert 'cabs_singles=2' in text
            assert (j.get('orbital_predecessor') is not None) == (i > 0)
        seed = next(iter(cfg['initial_records'].values()))
        assert target in Path(seed['path']).parents
        assert file_sha256(seed['path']) == seed['sha256']
    assert all(file_sha256(p) == sha for p, sha in before.items())
    assert adapter.prepare(source, target, plan, review) == record
    first = json.loads(Path(record['systems'][0]['manifest']).read_text())['jobs'][0]
    inp = Path(first['directory'])/(first['name']+'.inp')
    inp.write_text(inp.read_text()+'\nchanged\n')
    with pytest.raises(ValueError, match='input changed'):
        adapter.prepare(source, target, plan, review)


@pytest.mark.parametrize('plan_name,point_count', [('recovery_plan.json', 26), ('cabs3_gate_plan.json', 2)])
def test_all_recovery_native_transformations_keep_method_choices(tmp_path, monkeypatch, plan_name, point_count):
    from types import SimpleNamespace
    import re
    from kinbot import molpro_reconnaissance as recon
    adapter, source, target, plan, review = fixture(tmp_path, monkeypatch)
    plan = PACKAGE/plan_name
    record = adapter.prepare(source, target, plan, review)
    monkeypatch.setattr(recon, 'qc_runtime_environment', lambda *a, **k: ({}, {'scratch_directory': str(tmp_path)}))
    monkeypatch.setattr(recon, 'cleanup_qc_runtime', lambda *a: None)
    checked = []
    def inspect_launch(command, **kwargs):
        directory = Path(kwargs['cwd'])
        text = (directory/command[-1]).read_text()
        assert 'maxit,999;' in text
        for block in re.findall(r'\{(?:mrci-f12|rs2-f12)[^{}]*\}', text):
            assert 'thrvar=1.d-10' in block and 'maxvi=60' in block
            if 'mrci-f12' in block:
                assert 'option,cluster=0' in block
                nref = int(re.search(r'refstate,(\d+)', block)[1])
                assert f'option,nstati={nref};' in block
                if 'maxdav20' in directory.parent.name:
                    assert 'maxdav=20' in block and 'cabs_singles=1' in block
                elif 'cabs3_gate' in directory.parent.name:
                    assert 'cabs_singles=3' in block and 'maxdav=' not in block
                else:
                    assert 'cabs_singles=2' in block
            else:
                assert 'shift=0.3,ipea=0.25' in block
                assert 'refstate,' in block and 'weight,' in block
        checked.append(directory)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', inspect_launch)
    for s in record['systems']:
        m = json.loads(Path(s['manifest']).read_text())
        for j in m['jobs']:
            j['retry_solver_options'] = m['solver_options']
            if j.get('orbital_predecessor'):
                (Path(j['directory'])/(j['name']+'.wfu')).write_bytes(b'Synthetic propagation fixture')
            recon._execute(j)
    assert len(checked) == point_count


def test_cabs3_gates_only_change_cabs_treatment_and_launch_no_dependents(tmp_path, monkeypatch):
    import re
    adapter, source, target, _, review = fixture(tmp_path, monkeypatch)
    before = {p: file_sha256(p) for p in source.rglob('*') if p.is_file()}
    original = json.loads((source/'confirmation.json').read_text())
    source_cfgs = {s['title']: json.loads(Path(s['configuration']).read_text()) for s in original['systems']}
    record = adapter.prepare(source, target, PACKAGE/'cabs3_gate_plan.json', review)
    assert record['execution']['max_jobs'] == 2
    assert record['execution']['memory_mb'] == 104000
    assert record['execution']['excluded_nodes'] == original['execution']['excluded_nodes']
    assert record['production_approved'] is False
    assert len(record['systems']) == 2
    for s in record['systems']:
        cfg = json.loads(Path(s['configuration']).read_text())
        prior = source_cfgs[s['title'].removesuffix('_cabs3_gate')]
        expected_space = deepcopy(prior['probe_definition'])
        expected_space['state_model']['correlation_options']['cabs_singles'] = 3
        assert cfg['probe_definition'] == expected_space
        assert cfg['correlation_checks'] == prior['correlation_checks']
        assert cfg['bases'] == prior['bases']
        far = next(a for a in prior['anchors'] if a['name'] == 'far')
        assert cfg['anchors'][0]['geometry_angstrom'] == far['geometry_angstrom']
        seed = next(iter(cfg['initial_records'].values()))
        assert seed['sha256'] == next(iter(prior['initial_records'].values()))['sha256']
        assert file_sha256(seed['path']) == seed['sha256']
        jobs = json.loads(Path(s['manifest']).read_text())['jobs']
        assert len(jobs) == 1 and jobs[0]['anchor'] == 'far'
        assert jobs[0].get('orbital_predecessor') is None
        text = (Path(jobs[0]['directory'])/(jobs[0]['name']+'.inp')).read_text()
        blocks = re.findall(r'\{mrci-f12[^{}]*\}', text)
        roots = 3 if s['title'].startswith('ch3s_cs_') else 6
        assert len(blocks) == 2*roots
        for spin in (1, 3):
            for root in range(1, roots+1):
                assert any(f'wf,{25 if roots == 3 else 41},1,{spin};state,1,{root};' in b for b in blocks)
        for block in blocks:
            assert 'cabs_singles=3' in block
            assert f'refstate,{roots},'+','.join(map(str, range(1, roots+1)))+';' in block
            assert f'option,nstati={roots};' in block
            assert 'option,cluster=0' in block and 'maxdav=' not in block
        assert text.count('KBVRCSEL=energd(1)') == 2*roots
    assert all(file_sha256(p) == sha for p, sha in before.items())


@pytest.mark.parametrize('defect', ['live', 'changed_output', 'changed_checkpoint', 'unknown_scheduler'])
def test_recovery_refuses_stale_or_live_source_before_creating_run(tmp_path, monkeypatch, defect):
    adapter, source, target, plan, review = fixture(tmp_path, monkeypatch)
    if defect == 'live':
        monkeypatch.setattr(adapter.run_cleanup, 'inspect', lambda root: {'jobs': ['123'], 'processes': []})
    elif defect == 'unknown_scheduler':
        def unavailable(root):
            raise RuntimeError('scheduler unavailable')
        monkeypatch.setattr(adapter.run_cleanup, 'inspect', unavailable)
    else:
        suffix = '*.out' if defect == 'changed_output' else '*.wfu'
        next(source.rglob(suffix)).write_bytes(b'changed')
    with pytest.raises((RuntimeError, ValueError)):
        adapter.prepare(source, target, plan, review)
    assert not target.exists()


def test_same_geometry_subspace_review_allows_internal_rotation_but_detects_leakage(monkeypatch):
    review = module('review_export')
    a = {'atoms': ['same'], 'basis': 'same', 'coefficients': np.eye(4)}
    b = deepcopy(a)
    theta = .6
    b['coefficients'][1:3, 1:3] = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    monkeypatch.setattr(review, 'orbital_data', lambda p: a if p == 'a' else b)
    result = review.same_geometry_overlap('a', 'b', 1, 2)
    assert min(result['active_singular_values']) == pytest.approx(1.)
    b['coefficients'] = np.eye(4)[:, [0, 3, 2, 1]]
    assert min(review.same_geometry_overlap('a', 'b', 1, 2)['active_singular_values']) == pytest.approx(0.)
    b['atoms'] = ['changed']
    with pytest.raises(ValueError, match='identical coordinates'):
        review.same_geometry_overlap('a', 'b', 1, 2)
