import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import pes_phase as phase, live_control as control
from kinbot.execution import atomic_json


@pytest.fixture
def root(tmp_path, monkeypatch):
    monkeypatch.setattr(phase, '_ROOT', None)
    monkeypatch.setattr(phase, '_WORKER', None)
    monkeypatch.setattr(phase, '_EPOCH', None)
    return tmp_path


def state(root, **values):
    m = json.loads((root/phase.NAME).read_text())
    m.update(values)
    atomic_json(root/phase.NAME, m)
    return m


def worker_status(root, name, status):
    m = state(root)
    m['workers'][name]['status'] = status
    state(root, workers=m['workers'])


def test_serial_exploration_does_not_serialize_downstream(root):
    phase.initialize(root)
    phase.register(root, 'first', os.getpid())
    busy, status = phase.tick(root, ['first','second'], max_searches=1)
    assert busy == ['first'] and status == 'exploring'
    worker_status(root,'first','ready')
    busy, status = phase.tick(root, ['first','second'], max_searches=1)
    assert busy == [] and status == 'exploring'  # Free slot for second well.
    phase.register(root,'second', os.getpid())
    worker_status(root,'second','ready')
    assert phase.tick(root,['first','second'],max_searches=1) == ([], 'postprocessing')
    # Both may postprocess even though exploration concurrency is one.
    worker_status(root,'first','postprocessing')
    worker_status(root,'second','postprocessing')
    assert phase.tick(root,['first','second'],max_searches=1)[1] == 'postprocessing'
    worker_status(root,'first','assembled')
    assert phase.tick(root,['first','second'],max_searches=1)[1] != 'closing'
    worker_status(root,'second','assembled')
    assert phase.tick(root,['first','second'],max_searches=1)[1] == 'closing'


def test_discovered_well_extends_barrier(root):
    phase.initialize(root)
    phase.register(root,'a',os.getpid())
    worker_status(root,'a','ready')
    assert phase.tick(root,['a','new'],max_searches=1)[1] == 'exploring'


def test_live_reopen_quiesces_both_branches_before_retry(root, monkeypatch):
    control.configure(root)
    record = phase.initialize(root)
    for name in ('a','b'):
        (root/name).mkdir()
        phase.register(root,name,os.getpid())
        worker_status(root,name,'postprocessing')
    state(root,phase='postprocessing')
    target = root/'a'/'failed_reaction'
    control.search_scope(target,'reaction',-999,close=True)
    request = control.request(root,target,kind='reaction',reason='retry rejected channel')
    assert phase.tick(root,['a','b'],max_searches=1)[1] == 'exploring'
    assert state(root)['epoch'] == 1 and state(root).get('grants') == []
    monkeypatch.setattr(phase,'_ROOT',root)
    monkeypatch.setattr(phase,'_WORKER','a')
    monkeypatch.setattr(phase,'_EPOCH',0)
    monkeypatch.setenv('KINBOT_PES_PHASE_TOKEN',record['token'])
    with pytest.raises(phase.ExplorationReopened):
        phase.checkpoint()
    worker_status(root,'a','ready')
    phase.tick(root,['a','b'],max_searches=1)
    assert state(root).get('grants') == []
    worker_status(root,'b','ready')
    assert phase.tick(root,['a','b'],max_searches=1)[0] == ['a']
    control.acknowledge(root/'.kinbot_controls/requests'/(request['id']+'.json'),request,'applied')
    assert phase.tick(root,['a','b','new'],max_searches=1)[1] == 'exploring'


def test_closed_search_refuses_dead_owner_or_finalized_pes(root):
    control.configure(root)
    phase.initialize(root)
    (root/'a').mkdir()
    phase.register(root,'a',os.getpid())
    target = root/'a'/'reaction'
    control.search_scope(target,'reaction',-999,close=True)
    state(root,phase='closing')
    with pytest.raises(ValueError,match='handed off'):
        control.request(root,target,kind='reaction',reason='late')


def test_checkpoint_noop_outside_managed_pes(root):
    phase.checkpoint()
    phase.wait_for_exploration()
    phase.wait_for_selection()
    phase.assembled()


def test_early_vrc_releases_search_slot_without_opening_anl(root, monkeypatch):
    record = phase.initialize(root)
    phase.register(root, 'a', os.getpid())
    monkeypatch.setattr(phase, '_ROOT', root)
    monkeypatch.setattr(phase, '_WORKER', 'a')
    monkeypatch.setenv('KINBOT_PES_PHASE_TOKEN', record['token'])
    phase.wait_for_exploration(allow_vrc=True)
    phase.checkpoint()  # Exploring globally is legal for VRC in this epoch.
    assert state(root)['workers']['a']['status'] == 'vrc'
    assert phase.tick(root, ['a','b'], max_searches=1) == ([], 'exploring')
    assert state(root)['waiting_workers'] == ['b']
    phase.register(root, 'b', os.getpid())
    worker_status(root, 'b', 'ready')
    assert phase.tick(root, ['a','b'], max_searches=1) == ([], 'postprocessing')
    phase.wait_for_selection()
    assert state(root)['workers']['a']['status'] == 'postprocessing'


def test_revision_during_early_vrc_waits_for_branch_acknowledgement(root, monkeypatch):
    control.configure(root)
    record = phase.initialize(root)
    for name in ('a','b'):
        (root/name).mkdir()
        phase.register(root, name, os.getpid())
        worker_status(root, name, 'vrc')
    phase.register(root, 'still_exploring', os.getpid())
    target = root/'a'/'failed_reaction'
    control.search_scope(target, 'reaction', -999, close=True)
    control.request(root, target, kind='reaction', reason='review while VRC runs')
    monkeypatch.setattr(phase, '_ROOT', root)
    monkeypatch.setattr(phase, '_WORKER', 'a')
    monkeypatch.setattr(phase, '_EPOCH', 0)
    monkeypatch.setenv('KINBOT_PES_PHASE_TOKEN', record['token'])
    for _ in range(3):
        phase.tick(root, ['a','b','still_exploring'], max_searches=2)
        assert state(root)['epoch'] == 1  # One invalidation, not one per tick.
        assert state(root)['grants'] == []
    with pytest.raises(phase.ExplorationReopened):
        phase.wait_for_selection()
    worker_status(root, 'a', 'ready')
    phase.tick(root, ['a','b','still_exploring'], max_searches=2)
    assert state(root)['grants'] == []
    worker_status(root, 'b', 'ready')
    assert phase.tick(root, ['a','b','still_exploring'], max_searches=2)[0] == ['still_exploring','a']
    assert not state(root)['quiescing']


def test_concurrent_discovery_publishes_one_input_without_overwriting(root,monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    import threading
    from kinbot import pes
    (root/'chemids').write_text('old\n')
    monkeypatch.setattr(pes,'routing_name',lambda s:'new')
    writes=[]
    monkeypatch.setattr(pes,'write_input',lambda *a:writes.append(a))
    barrier=threading.Barrier(2)
    def publish():
        barrier.wait()
        return pes.publish_well('input.json',object(),10,20,root,2)
    with ThreadPoolExecutor(2) as pool:
        futures=[pool.submit(publish) for _ in range(2)]
        assert sorted(f.result(timeout=5) for f in futures)==[False,True]
    assert len(writes)==1 and (root/'chemids').read_text()=='old\nnew\n'


def test_new_well_published_between_inventory_and_tick_keeps_barrier_closed(root):
    phase.initialize(root)
    phase.register(root,'a',os.getpid())
    worker_status(root,'a','ready')
    (root/'chemids').write_text('a\nnew\n')
    assert phase.tick(root,['a'],max_searches=1,known_jobs=['a'])[1]=='exploring'


def test_real_serial_workers_all_explore_before_any_postprocessing(root):
    import subprocess
    import sys
    import time
    record=phase.initialize(root)
    child='''from pathlib import Path
import time
from kinbot import pes_phase as p
assert p.attach()
Path('explored').write_text(str(time.time_ns()))
p.wait_for_exploration()
Path('postprocessing').write_text(str(time.time_ns()))
p.assembled()
'''
    processes={}
    env={**os.environ,'PYTHONPATH':str(Path(__file__).parents[1]),
         'KINBOT_PES_PHASE_ROOT':str(root),'KINBOT_PES_PHASE_TOKEN':record['token']}
    try:
        deadline=time.monotonic()+15
        while time.monotonic()<deadline:
            busy,status=phase.tick(root,['a','b'],max_searches=1)
            if not busy and len(processes)<2:
                name=('a','b')[len(processes)]
                (root/name).mkdir()
                proc=subprocess.Popen([sys.executable,'-c',child],cwd=root/name,env=env)
                processes[name]=proc
                phase.register(root,name,proc.pid)
            if len(processes)==2 and all(p.poll() is not None for p in processes.values()):
                break
            time.sleep(.03)
        assert len(processes)==2
        assert all(p.poll()==0 for p in processes.values())
        assert max(int((root/n/'explored').read_text()) for n in processes) < min(
            int((root/n/'postprocessing').read_text()) for n in processes)
    finally:
        for p in processes.values():
            if p.poll() is None:
                p.terminate()
            p.wait(timeout=5)


@pytest.mark.parametrize('slots', [1, 8])
def test_real_workers_run_vrc_early_and_anl_during_vrc_hold(root, slots):
    """Use the real branch join and owner loop, including a single-slot input."""
    import subprocess
    import sys
    import time
    record = phase.initialize(root)
    child = '''from pathlib import Path
from types import SimpleNamespace
import time
from kinbot import pes_phase as p, kb
from kinbot.anl import automatic
assert p.attach()
root = Path.cwd().parent
if Path.cwd().name == 'a':
    def prepare(*args, **kwargs):
        Path('planned').touch()
        return {}
    def vrc(*args):
        Path('vrc_started').touch()
        while not (root/'release_vrc').exists():
            p.checkpoint()
            time.sleep(.02)
        return {'status':'complete'}
    def anl(*args):
        Path('anl_done').touch()
        return {'status':'complete'}
    automatic.prepare_automatic_cbh_anl = prepare
    automatic.validate_mess_environment = lambda *a: None
    kb._execute_vrc_branch = vrc
    kb._execute_automatic_branch = anl
    p.wait_for_exploration(allow_vrc=True)
    well = SimpleNamespace(reac_obj=[], reac_ts_done=[])
    kb._post_reaction_branches(well, dict(queue_job_limit=SLOTS,
        vrc_tst_scan={'a':['scission']}, vrc_tst_max_nodes=SLOTS,
        rotdpy_max_jobs=SLOTS), object(), 'input.json', True)
else:
    # The second well starts while the first is running/reviewing VRC.
    assert (root/'a/vrc_started').exists()
    assert not (root/'a/planned').exists()
    Path('explored').touch()
    p.wait_for_exploration()
p.wait_for_selection()
p.assembled()
'''.replace('SLOTS', str(slots))
    processes = {}
    env = {**os.environ, 'PYTHONPATH':str(Path(__file__).parents[1]),
           'OPENBLAS_NUM_THREADS':'1', 'OMP_NUM_THREADS':'1',
           'KINBOT_PES_PHASE_ROOT':str(root), 'KINBOT_PES_PHASE_TOKEN':record['token']}
    try:
        deadline = time.monotonic()+25
        while time.monotonic() < deadline:
            busy, status = phase.tick(root, ['a','b'], max_searches=1)
            if not busy and len(processes) < 2 and (not processes or (root/'a/vrc_started').exists()):
                name = ('a','b')[len(processes)]
                (root/name).mkdir()
                proc = subprocess.Popen([sys.executable,'-c',child], cwd=root/name, env=env)
                processes[name] = proc
                phase.register(root, name, proc.pid)
            if (root/'a/anl_done').exists() and not (root/'release_vrc').exists():
                assert (root/'b/explored').exists()
                assert status == 'postprocessing'
                (root/'release_vrc').touch()
            join = root/'a/.kinbot_post_reaction_branches.json'
            if join.exists():
                branches = json.loads(join.read_text())['branches']
                assert all(b['status'] != 'failed' for b in branches.values()), branches
            if len(processes)==2 and all(p.poll() is not None for p in processes.values()):
                break
            time.sleep(.03)
        assert len(processes)==2 and all(p.poll()==0 for p in processes.values())
        assert (root/'a/planned').exists() and (root/'a/anl_done').exists()
    finally:
        for proc in processes.values():
            if proc.poll() is None:
                proc.terminate()
            proc.wait(timeout=5)
