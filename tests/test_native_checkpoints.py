"""Durable candidates survive interrupted copies without becoming restart approvals."""
import json
import os
import threading
from pathlib import Path

import pytest

from kinbot import native_checkpoints as nc
from kinbot.execution import file_sha256


@pytest.fixture
def setup(tmp_path):
    result=tmp_path/'results';result.mkdir()
    scratch=tmp_path/'scratch';work=scratch/'work';work.mkdir(parents=True)
    (work/'MINP').write_text('calc=CCSDT(Q)\nrest=0\n')
    (result/'task.json').write_text('{"task":"original"}')
    runtime={'scratch_directory':str(scratch)}
    pub=nc.Publisher(result,work,runtime,'mrcc','test',{'MINP'},
                     config={'qc_checkpoint_interval':1,'qc_checkpoint_max_mb':1})
    return pub,work,result,runtime


def latest(pub):
    p=json.loads((pub.root/'latest.json').read_text())
    d=pub.root/p['generation'];m=json.loads((d/'manifest.json').read_text())
    assert file_sha256(d/'manifest.json')==p['sha256']
    for f in m['files']:assert file_sha256(d/f['path'])==f['sha256']
    return p,d,m


def test_generations_atomic_bounded_and_native_identity(setup):
    pub,work,result,runtime=setup
    (work/'fort.16').write_bytes(b'amplitudes1')
    (work/'huge.integrals').write_bytes(b'ignored')
    (work/'test.out').write_text('Iteration 14\n')
    assert pub.poll(force=True)
    p,d,m=latest(pub)
    assert m['contract']['task_sha256']==file_sha256(result/'task.json')
    assert m['contract']['inputs']=={'MINP':file_sha256(work/'MINP')}
    assert [x['source'] for x in m['files']]==[str(work/'fort.16')]
    assert m['status']=='restart_candidate_requires_native_validation'
    assert (d/'diagnostics/MINP').read_bytes()==(work/'MINP').read_bytes()
    assert pub.poll(force=True)==p
    for i in range(3):
        (work/'fort.16').write_bytes(f'amplitudes{i+2}'.encode());pub.poll(force=True)
    assert len(list(pub.root.glob('generation-*')))==2
    assert not d.exists()
    assert latest(pub)[2]['files'][0]['sha256']==file_sha256(work/'fort.16')
    other=nc.Publisher(result,work,runtime,'mrcc','test',{'MINP'},config=pub.config)
    other.poll(force=True)
    assert other.root!=pub.root and (pub.root/'latest.json').exists()


@pytest.mark.parametrize('change',['during_copy','between_files','added_file','stop','input'])
def test_invalid_snapshot_keeps_previous_generation(setup,monkeypatch,change):
    pub,work,_,_=setup
    p=work/'fort.16';p.write_bytes(b'first');pub.poll(force=True)
    before=(pub.root/'latest.json').read_bytes()
    p.write_bytes(b'next')
    native=nc._copy
    stop=threading.Event()
    def copying(source,target,expected,event):
        if change=='during_copy':source.write_bytes(b'changed-before-open')
        result=native(source,target,expected,event)
        if change=='between_files':source.write_bytes(b'changed-after-copy')
        if change=='added_file':(work/'fort.17').write_bytes(b'new')
        if change=='stop':stop.set()
        if change=='input':(work/'MINP').write_text('different theory')
        return result
    monkeypatch.setattr(nc,'_copy',copying)
    assert pub.poll(force=True,stop=stop) is None
    assert (pub.root/'latest.json').read_bytes()==before
    assert not list(pub.root.glob('*.preparing'))
    assert len(list(pub.root.glob('generation-*')))==1
    latest(pub)


def test_copy_cannot_chase_file_growth_or_follow_replacement(setup,monkeypatch):
    pub,work,result,_=setup
    source=work/'fort.16';source.write_bytes(b'a'*100)
    expected=nc.signature(source)
    original_open=nc.os.open
    def growing(path,flags,*args,**kwargs):
        if Path(path)==source:source.write_bytes(b'b'*200)
        return original_open(path,flags,*args,**kwargs)
    monkeypatch.setattr(nc.os,'open',growing)
    with pytest.raises(OSError,match='replaced'):nc._copy(source,result/'copy',expected,None)
    assert (result/'copy').stat().st_size==0


@pytest.mark.parametrize('kind',['symlink','fifo','oversized'])
def test_unsafe_or_oversized_candidates_not_published(setup,kind):
    pub,work,result,_=setup
    if kind=='symlink':(work/'fort.16').symlink_to(work/'MINP')
    elif kind=='fifo':os.mkfifo(work/'fort.16')
    else:(work/'fort.16').write_bytes(b'a'*(1024**2+1))
    assert pub.poll(force=True) is None
    assert not (pub.root/'latest.json').exists()
    assert json.loads((result/'test.checkpoint.json').read_text())['error']


def test_storage_outage_retries_on_deadline_without_resetting_good_copy(setup,monkeypatch):
    pub,work,_,runtime=setup
    (work/'fort.16').write_bytes(b'old');pub.poll(force=True)
    prior=(pub.root/'latest.json').read_bytes()
    (work/'fort.16').write_bytes(b'new')
    original=nc.publish;clock=[0.]
    monkeypatch.setattr(nc.time,'monotonic',lambda:clock[0])
    def offline(*a,**kw):raise OSError('shared filesystem unavailable')
    monkeypatch.setattr(nc,'publish',offline)
    pub.poll(force=True);assert pub.failures==1 and pub.next_at==60
    assert (pub.root/'latest.json').read_bytes()==prior
    monkeypatch.setattr(nc,'publish',original)
    clock[0]=59;assert pub.poll() is None
    clock[0]=60;assert pub.poll()
    assert pub.failures==0 and 'native_checkpoint_error' not in runtime
    latest(pub)


def test_pointer_publication_failure_keeps_previous(setup,monkeypatch):
    pub,work,_,_=setup
    (work/'fort.16').write_bytes(b'old');pub.poll(force=True)
    prior=(pub.root/'latest.json').read_bytes()
    (work/'fort.16').write_bytes(b'new')
    original=nc.publish
    def interrupted(path,data):
        if path.name=='latest.json':raise OSError('interrupted pointer update')
        original(path,data)
    monkeypatch.setattr(nc,'publish',interrupted)
    pub.poll(force=True)
    assert (pub.root/'latest.json').read_bytes()==prior
    latest(pub)
    monkeypatch.setattr(nc,'publish',original)
    pub.poll(force=True)
    assert len(list(pub.root.glob('generation-*')))==2
    latest(pub)


@pytest.mark.parametrize('backend,checkpoint',[('molpro','job.wfu'),('gaussian','job.chk'),('mrcc','fort.16')])
def test_workspace_preserves_checkpoint_on_native_failure(setup,monkeypatch,backend,checkpoint):
    from kinbot.qc_scratch import workspace
    from kinbot.anl.runtime import cleanup_qc_runtime
    _,_,result,runtime=setup
    # Workspace allocates its own work subdirectory.
    root=result.parent/'native';root.mkdir();runtime={'scratch_directory':str(root)}
    (result/'MINP').write_text('unchanged input')
    monkeypatch.setenv('KINBOT_QC_CHECKPOINT_INTERVAL','600')
    with pytest.raises(RuntimeError,match='native interrupted'):
        with workspace(result,runtime,backend=backend,inputs=['MINP'],label='job') as work:
            (work/checkpoint).write_bytes(b'restart candidate')
            (root/'Gaussian.rwf').write_bytes(b'Gaussian scratch')
            raise RuntimeError('native interrupted')
    cleanup_qc_runtime(runtime)
    assert not root.exists()
    marker=json.loads((result/'job.checkpoint.json').read_text())
    assert len(marker['checkpoint']['files'])==2
    for item in marker['checkpoint']['files']:
        assert file_sha256(Path(marker['checkpoint']['directory'])/item['path'])==item['sha256']
    assert (result/'MINP').read_text()=='unchanged input'


def test_final_failed_publication_keeps_scratch(setup,monkeypatch):
    from kinbot.qc_scratch import workspace
    from kinbot.anl.runtime import cleanup_qc_runtime
    _,_,result,_=setup
    root=result.parent/'native';root.mkdir();runtime={'scratch_directory':str(root)}
    (result/'MINP').write_text('unchanged input')
    monkeypatch.setenv('KINBOT_QC_CHECKPOINT_INTERVAL','600')
    monkeypatch.setattr(nc.Publisher,'capture',lambda *a,**kw: (_ for _ in ()).throw(OSError('offline')))
    with workspace(result,runtime,backend='mrcc',inputs=['MINP']) as work:
        (work/'fort.16').write_bytes(b'survives locally')
    cleanup_qc_runtime(runtime)
    assert runtime['retain_scratch'] and (work/'fort.16').exists()


def test_default_disabled_and_validation(setup):
    _,work,result,runtime=setup
    pub=nc.Publisher(result,work,runtime,'mrcc','job',[],config={})
    (work/'fort.16').write_bytes(b'a')
    assert pub.poll(force=True) is None and not pub.root.exists()
    for values in ({'qc_checkpoint_interval':True},{'qc_checkpoint_generations':1},
                   {'qc_checkpoint_max_mb':0},{'qc_checkpoint_interval':-1}):
        with pytest.raises(ValueError):nc.validate(values)


def test_periodic_workspace_checkpoint_survives_process_kill(tmp_path):
    """Exercise the real mirror thread and workspace, not normal finalization."""
    import subprocess
    import sys
    import time
    result=tmp_path/'results';result.mkdir();(result/'MINP').write_text('unchanged native input')
    scratch=tmp_path/'scratch';scratch.mkdir()
    driver=tmp_path/'driver.py'
    driver.write_text('''
import sys,time,threading
from pathlib import Path
from types import SimpleNamespace
from kinbot import qc_scratch
from kinbot.native_checkpoints import Publisher
# Speed the timer only. Publication and workspace code run unchanged.
real_event=threading.Event()
class FastEvent:
 def wait(self,seconds):return real_event.wait(.02)
 def set(self):real_event.set()
 def is_set(self):return real_event.is_set()
qc_scratch.threading=SimpleNamespace(Thread=threading.Thread,Event=FastEvent)
original=Publisher.__init__
def init(self,*a,**kw):
 original(self,*a,**kw);self.next_at=0
Publisher.__init__=init
result,scratch=map(Path,sys.argv[1:])
with qc_scratch.workspace(result,{'scratch_directory':str(scratch)},backend='mrcc',inputs=['MINP'],label='native') as work:
 (work/'fort.16').write_bytes(b'durable amplitudes')
 while True:time.sleep(.1)
''')
    env=dict(os.environ,KINBOT_QC_CHECKPOINT_INTERVAL='1',PYTHONPATH=str(Path.cwd()))
    with (tmp_path/'stderr').open('w') as err:
        child=subprocess.Popen([sys.executable,str(driver),str(result),str(scratch)],env=env,stderr=err)
        try:
            deadline=time.monotonic()+10
            pointers=[]
            while time.monotonic()<deadline and child.poll() is None:
                pointers=list(result.glob('.kinbot_checkpoints/native/native/*/latest.json'))
                if pointers:break
                time.sleep(.02)
            assert pointers,(tmp_path/'stderr').read_text()
            child.kill();child.wait(timeout=5)
        finally:
            if child.poll() is None:child.kill();child.wait(timeout=5)
    # Destroy ephemeral storage after the simulated node/process loss.
    import shutil
    shutil.rmtree(scratch)
    pointer=json.loads(pointers[0].read_text());generation=pointers[0].parent/pointer['generation']
    assert file_sha256(generation/'manifest.json')==pointer['sha256']
    m=json.loads((generation/'manifest.json').read_text())
    assert (generation/m['files'][0]['path']).read_bytes()==b'durable amplitudes'
    assert not (result/'fort.16').exists()  # normal exit copy-back never ran


def test_persistent_pointer_failure_cannot_accumulate_unbounded_generations(setup,monkeypatch):
    pub,work,_,_=setup
    (work/'fort.16').write_bytes(b'original');pub.poll(force=True)
    previous=(pub.root/'latest.json').read_bytes();original=nc.publish
    def offline_pointer(path,data):
        if path.name=='latest.json':raise OSError('pointer publication unavailable')
        return original(path,data)
    monkeypatch.setattr(nc,'publish',offline_pointer)
    for i in range(6):
        (work/'fort.16').write_bytes(str(i).encode());pub.poll(force=True)
        assert len(list(pub.root.glob('generation-*')))<=3  # two retained plus newly uncommitted
        assert (pub.root/'latest.json').read_bytes()==previous
        latest(pub)
