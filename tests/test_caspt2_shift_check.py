"""The production shift check must preserve the reviewed numerical experiment."""
from copy import deepcopy
import json
from pathlib import Path

from kinbot.execution import file_sha256
from test_molpro_confirmation_package import PACKAGE
from test_molpro_confirmation_recovery import fixture


def test_shift_check_changes_shift_only_and_preserves_predecessor(tmp_path, monkeypatch):
    adapter, source, target, _, review = fixture(tmp_path, monkeypatch)
    plan = json.loads((PACKAGE/'caspt2_shift01_plan.json').read_text())
    assert plan['max_jobs'] == 2
    # Synthetic source fixture has the original names. The shared-reference
    # controls in the plan must still be rendered, not assumed from a name.
    for item in plan['chains']:
        assert item['source_chain'].endswith('_shared_reference')
        item['source_chain'] = item['source_chain'].removesuffix('_shared_reference')
    p = tmp_path/'shift_plan.json'
    p.write_text(json.dumps(plan))
    before = {f: file_sha256(f) for f in source.rglob('*') if f.is_file()}
    record = adapter.prepare(source, target, p, review)
    assert len(record['systems']) == 2
    assert record['execution']['max_jobs'] == 2
    assert record['execution']['memory_mb'] == 104000
    original = json.loads((source/'confirmation.json').read_text())
    for system, item in zip(record['systems'], plan['chains']):
        cfg = json.loads(Path(system['configuration']).read_text())
        prior = next(s for s in original['systems'] if s['title'] == item['source_chain'])
        old_cfg = json.loads(Path(prior['configuration']).read_text())
        assert cfg['bases'] == old_cfg['bases']
        assert cfg['anchors'][0]['geometry_angstrom'] == next(
            a['geometry_angstrom'] for a in old_cfg['anchors'] if a['name'] == 'inner')
        old_manifest = json.loads(Path(prior['manifest']).read_text())
        center = next(j for j in old_manifest['jobs'] if j['anchor'] == 'center')
        seed = next(iter(cfg['initial_records'].values()))
        assert seed['sha256'] == center['wfu_sha256']
        assert cfg['probe_definition']['orbital_partition']['orbital_swaps'] == []
        jobs = json.loads(Path(system['manifest']).read_text())['jobs']
        assert len(jobs) == 1 and jobs[0]['anchor'] == 'inner'
        deck = (Path(jobs[0]['directory'])/(jobs[0]['name']+'.inp')).read_text()
        assert deck.count('{rs2-f12,shift=0.1,ipea=0.25') == 12
        assert deck.count('refstate,6,1,2,3,4,5,6') == 12
        assert deck.count('KBVRCSEL=energy(1)') == 12
        assert 'mrci-f12' not in deck and 'shift=0.3' not in deck
    assert all(file_sha256(f) == sha for f, sha in before.items())
