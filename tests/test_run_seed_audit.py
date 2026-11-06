import json
from pathlib import Path
import sqlite3
import tarfile

import pytest

from kinbot.execution import atomic_json, file_sha256
from kinbot import run_seed_audit as audit


def case(root):
    p=root/'example';p.mkdir(parents=True)
    atomic_json(p/'input.json',{'smiles':'test'})
    w=p/'worker';w.mkdir();atomic_json(w/'worker.json',{'smiles':'test'})
    atomic_json(w/'.kinbot_search_state.json',{'failed':{'status':-2},'accepted':{'status':-1}})
    return p,w


def test_full_database_with_wal_and_terminal_search_history_source_unchanged(tmp_path):
    root=tmp_path/'old';p,w=case(root)
    db=sqlite3.connect(w/'kinbot.db');db.execute('pragma journal_mode=wal')
    db.execute('create table evidence (status text)');db.execute('insert into evidence values (?)',('normal',));db.commit()
    native=w/'rotdPy/samples';native.mkdir(parents=True);(native/'huge.out').write_text('not a search ledger')
    before={f:file_sha256(f) for f in root.rglob('*') if f.is_file()}
    result=audit.capture(root,tmp_path/'audit')
    assert result['status']=='captured'
    output=Path(result['directory'])
    with sqlite3.connect(output/'run/example/worker/kinbot.db') as copied:
        assert copied.execute('select * from evidence').fetchall()==[('normal',)]
    assert json.loads((output/'run/example/worker/.kinbot_search_state.json').read_text())['failed']['status']==-2
    assert before=={f:file_sha256(f) for f in before}
    with tarfile.open(result['archive']) as archive:
        assert 'run/example/worker/kinbot.db-wal' in archive.getnames()
        assert not any('huge.out' in n for n in archive.getnames())
    db.close()


def test_selected_task_archives_and_checkpoint_receipts_without_binary_copy(tmp_path):
    root=tmp_path/'old';p,w=case(root)
    task=w/'anl_automatic/states/id/composite/tasks/geom';task.mkdir(parents=True)
    atomic_json(task/'task.json',{'id':'geom'});atomic_json(task.parent.parent/'workflow.json',{})
    old=task/'.kinbot_slurm_io_retries/attempt/interrupted';old.mkdir(parents=True)
    (old/'optimization.traj').write_bytes(b'original frames')
    (old/'large.wfu').write_bytes(b'not part of audit')
    evidence=tmp_path/'preserved';node=evidence/'checkpoints/nodes/node-a';node.mkdir(parents=True)
    (node/'native.out').write_text('native diagnostic')
    atomic_json(node/'preservation.json',{'files':[{'target':str(node/'native.out')}]})
    result=audit.capture(root,tmp_path/'audit',tasks=[str(task.relative_to(root))],evidence=evidence)
    manifest=json.loads((Path(result['directory'])/'manifest.json').read_text())
    names=[x['path'] for x in manifest['files']]
    assert any(n.endswith('optimization.traj') for n in names)
    assert any(n.endswith('native.out') for n in names)
    assert not any(n.endswith('.wfu') for n in names)


def test_changed_source_is_not_labelled_complete(tmp_path,monkeypatch):
    root=tmp_path/'old';p,w=case(root);original=audit.read_file
    def changing(source,*a,**k):
        result=original(source,*a,**k)
        if source.name=='worker.json':source.write_text('{"changed":true}')
        return result
    monkeypatch.setattr(audit,'read_file',changing)
    result=audit.capture(root,tmp_path/'audit')
    assert result['status']=='incomplete' and result['issues']==1


def test_escape_symlinks_and_capacity_are_explicit_omissions(tmp_path):
    root=tmp_path/'old';p,w=case(root)
    (w/'kinbot_monitor.out').symlink_to(tmp_path/'outside');(tmp_path/'outside').write_text('private')
    (w/'kinbot.db').write_bytes(b'0'*(2*1024**2))
    result=audit.capture(root,tmp_path/'audit',max_mb=1)
    assert result['status']=='incomplete' and result['issues']>=2
    with pytest.raises(ValueError,match='outside'):audit.capture(root,root/'audit')
    with pytest.raises(ValueError,match='relative'):audit.capture(root,tmp_path/'audit',tasks=['../other'])


def test_mess_handoff_roundtrip_replays_native_worker_artifacts(tmp_path):
    from types import SimpleNamespace
    from tests.test_mess_population_identity import (
        stage_collision_network,assemble_collision_network,assert_separate_complexes)
    from tests.test_mess_pes_artifacts import write_worker
    root=tmp_path/'source';scope=root/'case'
    cases=stage_collision_network(scope)
    atomic_json(scope/'input.json',cases[0].par)
    (scope/'chemids').write_text('\n'.join(c.worker.name for c in cases)+'\n')
    (scope/'molpro').mkdir()  # Empty backend directories are part of the handoff.
    for c in cases:write_worker(c)
    expected=assemble_collision_network(cases)
    for c in cases:
        (c.worker/'native.chk').write_bytes(b'unrelated checkpoint')
        huge=c.worker/'rotdPy/samples';huge.mkdir(parents=True)
        (huge/'ignored.out').write_text('unrelated sampling')
    before={p.relative_to(root):file_sha256(p) for p in root.rglob('*') if p.is_file()}
    result=audit.capture(root,tmp_path/'exports',mess_handoff=True)
    assert result['status']=='captured'
    dest=tmp_path/'replay';dest.mkdir()
    with tarfile.open(result['archive']) as archive:archive.extractall(dest,filter='data')
    manifest=json.loads((dest/'manifest.json').read_text())
    assert manifest['purpose']=='mess_handoff_replay' and not manifest['submits_jobs']
    assert (dest/'run/case/molpro').is_dir()
    assert all(file_sha256(dest/item['path'])==item['sha256'] for item in manifest['files'])
    assert not any(item['path'].endswith(('.chk','ignored.out')) for item in manifest['files'])
    copied=[SimpleNamespace(**dict(vars(c),root=dest/'run/case',
                  worker=dest/'run/case'/c.worker.name)) for c in cases]
    for _ in range(2):
        models=assemble_collision_network(copied)
        assert models==expected
        assert_separate_complexes(copied,models)
    assert before=={p.relative_to(root):file_sha256(p) for p in root.rglob('*') if p.is_file()}


def test_mess_handoff_records_late_artifact_and_symlink_as_incomplete(tmp_path,monkeypatch):
    root=tmp_path/'source';p,w=case(root)
    (w/'one.mess').write_text('artifact')
    (w/'two.mess').symlink_to(tmp_path/'outside')
    (tmp_path/'outside').write_text('off-root')
    original=audit.read_file
    def arriving(source,*a,**k):
        result=original(source,*a,**k)
        if source.name=='one.mess':(w/'late.mess').write_text('late worker output')
        return result
    monkeypatch.setattr(audit,'read_file',arriving)
    result=audit.capture(root,tmp_path/'exports',mess_handoff=True)
    manifest=json.loads((Path(result['directory'])/'manifest.json').read_text())
    assert result['status']=='incomplete'
    assert any('inventory changed' in issue['error'] for issue in manifest['issues'])
    assert any('symlink' in issue['error'] for issue in manifest['issues'])


def test_homolysis_native_captures_full_attempts_without_driving_or_mutating(tmp_path,monkeypatch):
    root=tmp_path/'source';p,w=case(root)
    search='worker_dissociation_saddle_1_2';parent='worker_hom_sci_1_2'
    atomic_json(w/'homolysis_checks.json',dict(checks=[dict(parent=parent,search=search,
                outcome='unresolved',scan_outcome='saddle_submitting',saddles=[])]))
    atomic_json(w/'.kinbot_search_state.json',dict(reactions={search:dict(status=-999,step=40,job_id='123')}))
    atomic_json(w/(search+'.dissociation_scan.json'),dict(outcome='saddle_submitting',
                points=[dict(step=0,job_id='122',energy_hartree=-1.)]))
    atomic_json(w/'.kinbot_theory_jobs.json',dict(jobs={search+'_IRC_F':dict(job_id='124'),
                                                       'unrelated':dict(job_id='999')}))
    natives={search+'.log':b'full native output\n'*5000,search+'.traj':b'trajectory',
             search+'.chk':b'checkpoint',search+'_high_freq.out':b'native frequencies',
             search+'_IRC_F.log':b'IRC', 'fragment_well.log':b'fragment state evidence',
             'other_tight_TS.log':b'alternative tight saddle',
             '.kinbot_timeout_retries/token/attempt/files/native.log':b'old attempt',
             'restart_gaussian_scf/archive/native.com':b'old SCF input',
             'perm/'+search+'.err':b'worker stderr'}
    for name,data in natives.items():
        path=w/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
    excluded=w/'vrctst/sampling/big.wfu';excluded.parent.mkdir(parents=True);excluded.write_bytes(b'not needed')
    other=p/'resolved';other.mkdir();atomic_json(other/'resolved.json',{})
    atomic_json(other/'homolysis_checks.json',dict(checks=[dict(outcome='parent_not_retained')]))
    (other/'native.log').write_text('resolved worker native output')
    before={f:file_sha256(f) for f in root.rglob('*') if f.is_file()}
    calls=[]
    def query(args):
        calls.append(args)
        assert args[0]=='sacct'
        return dict(argv=args,status='observed',returncode=0,stdout='retained accounting',stderr='')
    monkeypatch.setattr('kinbot.run_outage_audit.command',query)
    result=audit.capture(root,tmp_path/'exports',homolysis_native=True)
    assert result['status']=='captured' and result['inconclusive_checks']==1
    output=Path(result['directory'])
    for name,data in natives.items():assert (output/'run/example/worker'/name).read_bytes()==data
    assert not (output/'run/example/resolved/native.log').exists()
    assert not (output/'run/example/worker/vrctst/sampling/big.wfu').exists()
    assert before=={f:file_sha256(f) for f in before}
    assert len(calls)==1 and calls[0][calls[0].index('-j')+1]=='122,123,124'
    inventory=json.loads((output/'homolysis_inventory.json').read_text())
    assert inventory['holds'][0]['search_state']['status']==-999
    with tarfile.open(result['archive']) as archive:
        assert archive.extractfile('run/example/worker/'+search+'.log').read()==natives[search+'.log']
        assert {'summary.txt','scheduler.json','homolysis_inventory.json'}<=set(archive.getnames())


def test_all_ordinary_native_includes_unheld_workers_and_reverse_models(tmp_path, monkeypatch):
    root=tmp_path/'source';p,w=case(root)
    files = {
        w/'hir/point.log': b'full HIR native output',
        w/'conf/selected.log': b'normal selected conformer',
        w/'thermochemistry_model.json': b'{"source":"complete model"}',
        p/'reverse_abstractions/pair/tight_TS.log': b'reverse saddle',
        p/'me/mess_0000.inp': b'combined MESS',
    }
    for path, data in files.items():
        path.parent.mkdir(parents=True, exist_ok=True);path.write_bytes(data)
    excluded=w/'scratch/raw.rwf';excluded.parent.mkdir();excluded.write_bytes(b'raw scratch')
    before={f:file_sha256(f) for f in root.rglob('*') if f.is_file()}
    monkeypatch.setattr('kinbot.run_outage_audit.command',lambda args: pytest.fail('No homolysis job IDs to query'))
    result=audit.capture(root,tmp_path/'exports',ordinary_native=True)
    assert result['status']=='captured' and result['inconclusive_checks']==0
    with tarfile.open(result['archive']) as archive:
        for path, data in files.items():
            assert archive.extractfile('run/'+str(path.relative_to(root))).read()==data
        assert 'run/'+str(excluded.relative_to(root)) not in archive.getnames()
    assert before=={f:file_sha256(f) for f in before}
    assert json.loads((Path(result['directory'])/'manifest.json').read_text())['purpose']=='ordinary_native_audit'


def test_homolysis_live_changes_are_reported_not_reset(tmp_path,monkeypatch):
    root=tmp_path/'source';p,w=case(root)
    name='worker_dissociation_saddle_1_2'
    atomic_json(w/'homolysis_checks.json',dict(checks=[dict(search=name,outcome='unresolved')]))
    (w/(name+'.log')).write_text('running')
    original=audit.read_file
    def changing(source,*a,**k):
        data=original(source,*a,**k)
        if source.name==name+'.log':source.write_text('running: more progress')
        return data
    monkeypatch.setattr(audit,'read_file',changing)
    monkeypatch.setattr('kinbot.run_outage_audit.command',lambda args:pytest.fail('No job IDs to query'))
    result=audit.capture(root,tmp_path/'exports',homolysis_native=True)
    assert result['status']=='incomplete' and result['issues']==1
    assert result['inconclusive_checks'] is None
    assert json.loads((w/'homolysis_checks.json').read_text())['checks'][0]['outcome']=='unresolved'


def test_homolysis_scheduler_batches_all_recorded_ids_and_retains_unavailability(tmp_path,monkeypatch):
    from kinbot.homolysis_evidence import scheduler
    worker=tmp_path/'run/case/worker';worker.mkdir(parents=True)
    name='worker_dissociation_saddle_1_2'
    atomic_json(worker/(name+'.dissociation_scan.json'),dict(points=[dict(job_id=str(i)) for i in range(201)]))
    calls=[]
    def unavailable(args):
        calls.append(args)
        return dict(argv=args,status='unknown',error='scheduler unavailable')
    monkeypatch.setattr('kinbot.run_outage_audit.command',unavailable)
    result=scheduler(tmp_path,[dict(worker='case/worker',check=dict(search=name),search_state={})])
    assert len(result['job_ids'])==201 and len(calls)==3
    assert all(q['status']=='unknown' for q in result['accounting'])


def test_native_inventory_walks_nested_scopes_once_and_keeps_attempts(tmp_path, monkeypatch):
    from kinbot import homolysis_evidence as native
    root=tmp_path/'source';p,w=case(root)
    (w/'native.log').write_text('evidence')
    attempt=w/'restart_gaussian_scf/attempt';attempt.mkdir(parents=True)
    (attempt/'native.com').write_text('old settings')
    walked=[];real=native.os.walk
    def walk(path, **kwargs):
        walked.append(path)
        yield from real(path, **kwargs)
    monkeypatch.setattr(native.os,'walk',walk)
    found=native.sources([w,p,w],all_workers=True)
    assert walked==[p]
    assert {w/'native.log',attempt/'native.com'}<=found


def test_capture_reports_all_phases_and_fast_archive_is_verified_before_publish(tmp_path, monkeypatch):
    root=tmp_path/'source';p,w=case(root)
    records=[];levels=[];real_open=audit.tarfile.open;real_verify=audit._verify_archive
    def opened(*args,**kwargs):
        if len(args)>1 and args[1]=='w:gz':levels.append(kwargs.get('compresslevel'))
        return real_open(*args,**kwargs)
    def verified(partial,*args):
        assert partial.name.endswith('.partial')
        assert not partial.with_suffix('').exists()
        return real_verify(partial,*args)
    monkeypatch.setattr(audit.tarfile,'open',opened)
    monkeypatch.setattr(audit,'_verify_archive',verified)
    result=audit.capture(root,tmp_path/'exports',ordinary_native=True,progress_callback=records.append)
    assert levels==[1]
    assert {r['phase'] for r in records}=={
        'inventory','checking_size','copying','rechecking_source_hashes','rechecking_inventory',
        'summarizing_copies','scheduler_accounting','compressing_archive',
        'verifying_archive','hashing_archive','finished'}
    assert Path(result['archive']).is_file()
    assert not Path(result['archive']+'.partial').exists()
    assert json.loads((Path(result['directory'])/'result.json').read_text())==result


def test_failed_archive_verification_keeps_partial_and_never_publishes_final(tmp_path, monkeypatch):
    root=tmp_path/'source';case(root)
    def bad(archive, expected, progress):
        expected=dict(expected);expected[next(iter(expected))]='not the checksum'
        real(archive,expected,progress)
    real=audit._verify_archive
    monkeypatch.setattr(audit,'_verify_archive',bad)
    with pytest.raises(RuntimeError,match='Archive verification mismatch'):
        audit.capture(root,tmp_path/'exports',ordinary_native=True)
    assert len(list((tmp_path/'exports').glob('*/seed-evidence.tar.gz.partial')))==1
    assert not list((tmp_path/'exports').glob('*/seed-evidence.tar.gz'))
    assert not list((tmp_path/'exports').glob('*/result.json'))


def test_detached_cli_finishes_after_launcher_exits_and_ignores_hangup(tmp_path):
    import os
    import signal
    import subprocess
    import sys
    import time
    root=tmp_path/'source';p,w=case(root)
    # Enough real files to signal during collection, not an already dead PID.
    for i in range(1200):(w/f'point_{i}.log').write_text('native output\n'*40)
    before={f:file_sha256(f) for f in root.rglob('*') if f.is_file()}
    command=[sys.executable,'-B','-m','kinbot.run_seed_audit',str(root),
             '--ordinary-native','--mess-handoff','--output-parent',str(tmp_path/'exports'),'--detach']
    launched=subprocess.run(command,capture_output=True,text=True,timeout=15,check=True)
    receipt=json.loads(launched.stdout);job=Path(receipt['job_directory'])
    deadline=time.monotonic()+30;sent=False
    while time.monotonic()<deadline:
        state=json.loads((job/'status.json').read_text())
        if state['state']=='running' and not sent:
            assert state['session_id']==receipt['pid']==state['pid']
            os.kill(receipt['pid'],signal.SIGHUP);sent=True
        if state['state'] in {'complete','incomplete','failed'}:break
        time.sleep(.01)
    assert sent
    assert state['state']=='complete',state
    result=json.loads((job/'result.json').read_text())
    assert Path(result['archive']).is_file() and result['issues']==0
    assert file_sha256(result['archive'])==result['sha256']
    assert before=={f:file_sha256(f) for f in before}
    # The read-only status command works from a new login without old variables.
    status=subprocess.run([sys.executable,'-B','-m','kinbot.run_seed_audit',
        '--status',str(tmp_path/'exports')],capture_output=True,text=True,check=True,timeout=15)
    assert str(job) in status.stdout and '"state": "complete"' in status.stdout


def test_detached_child_failure_is_saved_and_does_not_claim_completion(tmp_path):
    import subprocess
    import sys
    import time
    root=tmp_path/'unrecognized';root.mkdir()
    launched=subprocess.run([sys.executable,'-B','-m','kinbot.run_seed_audit',str(root),
        '--output-parent',str(tmp_path/'exports'),'--detach'],
        capture_output=True,text=True,check=True,timeout=15)
    job=Path(json.loads(launched.stdout)['job_directory'])
    deadline=time.monotonic()+15
    while time.monotonic()<deadline:
        state=json.loads((job/'status.json').read_text())
        if state['state']=='failed':break
        time.sleep(.01)
    assert state['state']=='failed' and 'No recognized seed evidence' in state['error']
    assert not (job/'result.json').exists()
    assert not list(job.glob('*/seed-evidence.tar.gz'))


def test_rapid_repeat_launch_keeps_child_result_and_latest_status(tmp_path, monkeypatch, capsys):
    from types import SimpleNamespace
    root=tmp_path/'source';case(root)
    args=SimpleNamespace(root=root,output_parent=tmp_path/'exports',evidence=None,max_mb=1,
                         task=[],mess_handoff=False,homolysis_native=False,ordinary_native=False)
    def quick(command, **kwargs):
        assert kwargs['start_new_session'] and kwargs['close_fds']
        assert kwargs['stdin']==audit.subprocess.DEVNULL
        job=Path(command[command.index('--job-directory')+1])
        atomic_json(job/'status.json',dict(state='complete',test='child finished first'))
        return SimpleNamespace(pid=1234)
    monkeypatch.setattr(audit.subprocess,'Popen',quick)
    first=audit.detach(args);second=audit.detach(args)
    assert first['job_directory']!=second['job_directory']
    for result in (first,second):
        assert json.loads(Path(result['status_file']).read_text())['state']=='complete'
    audit.show_status(args.output_parent)
    assert second['job_directory'] in capsys.readouterr().out


def test_periodic_checkpoint_store_is_never_walked_but_current_and_failed_attempts_are_captured(tmp_path, monkeypatch):
    import os
    root=tmp_path/'source';p,w=case(root)
    # Reproduce the exported paths, including checkpoint stores nested in HIR.
    stores=[w/'.kinbot_checkpoints',w/'hir/.kinbot_checkpoints',
            p/'reverse_abstractions/pair/.kinbot_checkpoints']
    for store in stores:
        old=store/'native/job/attempt/generation-old/files/work';old.mkdir(parents=True)
        (old/'job.chk').write_bytes(b'restart generation')
    files={w/'job.chk':b'current checkpoint',w/'job.log':b'full native output',
           w/'job.checkpoint.json':b'{"generation":"retained"}',
           w/'hir/point.log':b'rotor',w/'conf/point.log':b'conformer',
           w/'restart_gaussian_scf/attempt/job.log':b'failed SCF output',
           p/'reverse_abstractions/pair/job.log':b'reverse native output'}
    for path,data in files.items():
        path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
    name='worker_dissociation_saddle_1_2'
    atomic_json(w/'homolysis_checks.json',dict(checks=[dict(search=name,outcome='no_validated_saddle')]))
    before={f:file_sha256(f) for f in root.rglob('*') if f.is_file()}
    original=os.scandir
    def bounded(path):
        assert not any(Path(path).is_relative_to(store) for store in stores),path
        return original(path)
    monkeypatch.setattr(os,'scandir',bounded)
    result=audit.capture(root,tmp_path/'exports',ordinary_native=True,mess_handoff=True)
    assert result['status']=='captured' and result['inconclusive_checks']==1
    output=Path(result['directory'])
    manifest=json.loads((output/'manifest.json').read_text())
    assert set(map(str,stores))<=set(manifest['excluded_trees'])
    with tarfile.open(result['archive']) as archive:
        assert not any('.kinbot_checkpoints' in n for n in archive.getnames())
        for path,data in files.items():
            assert archive.extractfile('run/'+str(path.relative_to(root))).read()==data
    assert before=={f:file_sha256(f) for f in before}


def test_native_size_budget_fails_before_any_file_copy(tmp_path, monkeypatch):
    root=tmp_path/'source';p,w=case(root)
    (w/'native.log').write_bytes(b'0'*(2*1024**2))
    monkeypatch.setattr(audit,'read_file',lambda *a,**k:pytest.fail('Budget must be checked before copying'))
    with pytest.raises(ValueError,match='no native files copied'):
        audit.capture(root,tmp_path/'exports',ordinary_native=True,max_mb=1)
    plan=json.loads(next((tmp_path/'exports').glob('*/selection.json')).read_text())
    assert plan['selected_bytes']>plan['budget_bytes']
    assert plan['largest_files'][0]['path'].endswith('native.log')
    assert not list((tmp_path/'exports').glob('*/run'))
    assert not list((tmp_path/'exports').glob('*/seed-evidence.tar.gz*'))


def test_native_inventory_limit_is_explicit_and_does_not_copy(tmp_path, monkeypatch):
    from kinbot import homolysis_evidence as native
    root=tmp_path/'source';case(root)
    monkeypatch.setattr(native,'MAX_NATIVE_FILES',1)
    monkeypatch.setattr(audit,'read_file',lambda *a,**k:pytest.fail('Inventory must finish before copying'))
    with pytest.raises(ValueError,match='inventory exceeds file limit'):
        audit.capture(root,tmp_path/'exports',ordinary_native=True)
    assert not list((tmp_path/'exports').glob('*/run'))
