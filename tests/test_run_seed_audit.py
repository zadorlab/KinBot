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
