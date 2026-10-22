"""Restart and live discovery restoration must retain existing QC ownership."""
import json
import os

import pytest

from kinbot import pes_inventory as inventory, pes_phase
from kinbot.execution import atomic_json
from kinbot.run_format import ensure_current_run


@pytest.fixture
def root(tmp_path):
    ensure_current_run(tmp_path, create=True)
    for name in ('123', '456', '789'):
        directory = tmp_path/name
        ensure_current_run(directory, create=True)
        atomic_json(directory/(name+'.json'), dict(title=name, pes=1))
        (directory/'kinbot.db').write_bytes(b'accepted and allocated work')
    return tmp_path


def summary(root, names):
    p = root/'saved_summary.txt'
    p.write_text(f'Total\t\t{len(names)}\nRunning\t{len(names)}\nFinished\t0\n\n'
                 'Running:\n' + ''.join('\t'+n+'\n' for n in names) +
                 '\nFinished:\n\nWaiting:\n')
    return p


def test_restart_retains_all_discovered_wells_and_order(root):
    (root/'chemids').write_text('123\n789\n456\n')
    before = {p: p.read_bytes() for p in root.glob('[0-9]*/*') if p.is_file()}
    for _ in range(2):
        assert inventory.initialize(root, '123') == ['123', '789', '456']
    assert before == {p: p.read_bytes() for p in before}
    assert (root/'chemids').read_text() == '123\n789\n456\n'


def test_fresh_inventory_contains_only_initial_well(root):
    # Stray prepared directories are not discoveries.
    assert inventory.initialize(root, '123') == ['123']


def test_live_restore_is_append_only_idempotent_and_does_not_start_workers(root):
    inventory.initialize(root, '123')
    pes_phase.initialize(root)
    pes_phase.register(root, '123', os.getpid())
    phase = (root/pes_phase.NAME).read_bytes()
    saved = summary(root, ['456', '123', '789'])
    report = inventory.restore(root, saved)
    assert report['added'] == ['456', '789']
    assert (root/'chemids').read_text() == '123\n'
    inventory.restore(root, saved, apply=True)
    assert (root/'chemids').read_text() == '123\n456\n789\n'
    assert inventory.restore(root, saved, apply=True)['added'] == []
    assert (root/pes_phase.NAME).read_bytes() == phase
    # Owner's stale pre-tick list cannot close the L2 gate after the merge.
    record = json.loads(phase)
    record['workers']['123']['status'] = 'ready'
    atomic_json(root/pes_phase.NAME, record)
    assert pes_phase.tick(root, ['123'], max_searches=1, known_jobs=['123'])[1] == 'exploring'


@pytest.mark.parametrize('defect', ['missing', 'wrong_input', 'bad_summary', 'symlink'])
def test_restore_validates_every_candidate_before_writing(root, defect):
    inventory.initialize(root, '123')
    p = summary(root, ['123', '456', '789'])
    if defect == 'missing':
        (root/'789/789.json').unlink()
    elif defect == 'wrong_input':
        atomic_json(root/'789/789.json', dict(title='another', pes=1))
    elif defect == 'bad_summary':
        p.write_text('Total 3\nRunning:\n123\n456\n')
    else:
        (root/'789/789.json').unlink()
        (root/'789/789.json').symlink_to(root/'456/456.json')
    with pytest.raises((OSError, ValueError)):
        inventory.restore(root, p, apply=True)
    assert (root/'chemids').read_text() == '123\n'


def test_live_restore_cannot_expand_a_finalized_selection(root, monkeypatch):
    inventory.initialize(root, '123')
    record = pes_phase.initialize(root)
    record['phase'] = 'postprocessing'
    atomic_json(root/pes_phase.NAME, record)
    monkeypatch.setattr(pes_phase, 'owner_alive', lambda _: True)
    with pytest.raises(RuntimeError, match='downstream planning'):
        inventory.restore(root, summary(root, ['123', '456']), apply=True)
    assert (root/'chemids').read_text() == '123\n'


def test_handoff_restores_preceding_summary_without_qc_changes(root, monkeypatch):
    from kinbot import run_refresh
    campaign = root/'campaign'; campaign.mkdir()
    case = campaign/'ch3s'
    ensure_current_run(case, create=True)
    for n in ('123', '456'):
        ensure_current_run(case/n, create=True)
        atomic_json(case/n/(n+'.json'), dict(title=n, pes=1))
    inventory.initialize(case, '123')
    old = root/'campaign.controllers.old'; (old/'ch3s').mkdir(parents=True)
    summary(old/'ch3s', ['123', '456']).rename(old/'ch3s/pes_summary.txt')
    new = root/'campaign.controllers.new'; new.mkdir()
    atomic_json(new/'previous_refresh.json', dict(phase='launched', archive=str(old)))
    monkeypatch.setattr(run_refresh.run_salvage, '_cases', lambda _: ({}, [case/'input.json']))
    report = run_refresh._restore_discovery(campaign, new)
    assert report['ch3s']['added'] == ['456']
    assert (case/'chemids').read_text() == '123\n456\n'
