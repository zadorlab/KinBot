import json
from pathlib import Path
import sqlite3
import tarfile

from ase import Atoms
from ase.db import connect
import pytest

from kinbot import qc_decision_audit as audit


def put(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data) if isinstance(data, dict) else data)


def fixture(tmp_path):
    root = tmp_path/'run'
    case = root/'case'
    worker = case/'123'
    put(case/'input.json', dict(high_level=1, imagfreq_threshold=50.))
    put(worker/'123.json', dict(high_level=1, imagfreq_threshold=50.))
    put(worker/'kinbot_monitor.out', '-1 5 rx_good 456\n-1 5 rx_bad 456\n-999 3 rx_rejected\n')
    put(worker/'kinbot.log', 'rx_rejected Both IRCs lead to the initial well, identical reaction found\n')
    mol = Atoms('HOH', positions=[[0,0,0], [1,0,0], [1,1,0]])
    db = connect(worker/'kinbot.db')
    for name, freq in [('123_well_high', [100,200,300]), ('456_well_high',[100,200,300]),
                       ('rx_good_high',[-100,200,300]), ('rx_bad_high',[100,200,300])]:
        db.write(mol, name=name, data=dict(status='normal',energy=-100.,zpe=.01,frequencies=freq))
    db.write(mol,name='rx_good_conf_0001',data=dict(status='error'))
    db.write(mol,name='rx_rejected',data=dict(status='error'))
    return root, worker, db


def test_capture_is_read_only_and_does_not_equate_rejection_with_retry(tmp_path):
    root, worker, db = fixture(tmp_path)
    # A known SCF failure and its old input are evidence, not automatic retry permission.
    old = worker/'restart_gaussian_scf/old/probe.com'
    put(old, '#p B3LYP/cc-pVDZ force\n')
    put(worker/'rx_rejected.scf_recovery.json', dict(status='failed',
        retry_scf='CDIIS,conver=8',preserved_files={'input': {'archive': str(old.relative_to(worker))}}))
    put(worker/'rx_rejected.log', 'SCF Done: -100\nSCFChk: convergence 1e-6 required 1e-8\nError termination l701.exe\n')
    put(worker/'rx_good_high.log', 'Frequencies -- -100 200 300\nNormal termination\n')
    put(worker/'rx_rejected.py', 'raise AssertionError("Never execute this worker")\n')
    put(worker/'rx_rejected.optimization_budget.json', dict(status='exhausted', chunks_used=3,chunk_steps=250))
    put(worker/'rx_rejected_minimum_recovery/attempt_1_x/recovery.json', dict(status='no_lower_energy_trial'))
    put(worker/'rx_rejected_minimum_recovery/attempt_1_x/plus/probe.log', 'SCF Done: -100\n')
    put(root/'.kinbot_policy_history/000001.json', dict(optimization_step_retries=2))
    put(root/'case/me/mess_0000.inp', 'TemperatureList[K] 300 1000\n')
    before = {p.relative_to(root): p.read_bytes() for p in root.rglob('*') if p.is_file()}
    result = audit.capture(root, tmp_path/'export')
    assert before == {p.relative_to(root): p.read_bytes() for p in root.rglob('*') if p.is_file()}
    report = json.loads((tmp_path/'export/report.json').read_text())
    w = report['cases'][0]['workers'][0]
    checks = {c['job']:c for c in w['accepted_checks']}
    assert checks['rx_good_high']['status'] == 'cache_data_gates_pass'
    assert not checks['rx_good_high']['native_acceptance_verified']
    assert checks['rx_bad_high']['status'] == 'inspect_evidence'
    assert 'transition state' in str(checks['rx_bad_high']['questions'])
    assert all('rx_rejected' != c['job'] for c in w['accepted_checks'])
    assert w['rejected'][0]['category'] == 'returned_to_initial_well'
    assert not w['rejected'][0]['proves_channel_absent']
    assert 'SCFChk:' in str(w['native'])
    assert not report['retry_success_inferred'] and not report['global_snapshot']
    assert w['metadata']['rx_rejected.optimization_budget.json']['chunks_used'] == 3
    assert w['metadata']['rx_rejected_minimum_recovery/attempt_1_x/recovery.json']['status']=='no_lower_energy_trial'
    assert report['root_metadata']['.kinbot_policy_history/000001.json']['optimization_step_retries']==2
    assert 'mess_0000.inp' in report['cases'][0]['mess_inputs']
    assert result['issues'] == 0
    with tarfile.open(result['archive']) as tar:
        assert 'databases/case/123/kinbot.db' in tar.getnames()
        assert 'files/case/123/restart_gaussian_scf/old/probe.com' in tar.getnames()
        assert 'files/case/123/rx_rejected_minimum_recovery/attempt_1_x/plus/probe.log' in tar.getnames()
    with pytest.raises(ValueError, match='empty'):
        audit.capture(root,tmp_path/'export')
    with pytest.raises(ValueError, match='outside'):
        audit.capture(root,root/'export')


def test_newer_failed_row_is_not_replaced_by_older_normal_row(tmp_path):
    root, worker, db = fixture(tmp_path)
    db.write(Atoms('HOH'), name='rx_good_high',data={'status':'error'})
    result = audit.capture(root,tmp_path/'export')
    report = json.loads((tmp_path/'export/report.json').read_text())
    w = report['cases'][0]['workers'][0]
    assert w['database']['rx_good_high']['total_rows'] == 2
    assert [r['status'] for r in w['database']['rx_good_high']['latest_rows']] == ['error','normal']
    check = next(c for c in w['accepted_checks'] if c['job']=='rx_good_high')
    assert check['status']=='inspect_evidence'
    assert 'latest expected cache row is not normal' in str(check['questions'])


def test_database_backup_includes_committed_wal_and_makes_no_source_db_edits(tmp_path):
    source = tmp_path/'source.db'
    db = connect(source)
    db.write(Atoms('H'), name='atom', data={'status':'normal','energy':-1.,'zpe':0.})
    con = sqlite3.connect(source)
    try:
        con.execute('PRAGMA journal_mode=WAL')
        con.execute('PRAGMA wal_autocheckpoint=0')
        con.execute("UPDATE systems SET mtime=mtime+1 WHERE id=1")
        con.commit()
        before = source.read_bytes()
        found = audit.database_snapshot(source,tmp_path/'export/snapshot.db')
        assert found['atom'][0]['mtime_ase'] == con.execute('SELECT mtime FROM systems').fetchone()[0]
        assert source.read_bytes() == before
    finally:
        con.close()


def test_capture_handles_unstarted_worker_and_rejects_external_sources(tmp_path):
    root, worker, db = fixture(tmp_path)
    put(root/'case/999/999.json', {})
    outside = tmp_path/'outside.log'
    outside.write_text('SECRET OUTSIDE RUN')
    (worker/'external.log').symlink_to(outside)
    result = audit.capture(root,tmp_path/'export')
    assert result['issues'] == 1
    assert 'SECRET OUTSIDE RUN' not in (tmp_path/'export/evidence.txt').read_text()
    assert not (tmp_path/'export/files/case/123/external.log').exists()


def test_controller_symlink_cannot_supply_another_runs_rejection(tmp_path):
    root, worker, db = fixture(tmp_path)
    outside = tmp_path/'outside.log'
    outside.write_text('rx_rejected PRIVATE OTHER RUN Neither IRC leads to the initial well\n')
    (worker/'kinbot.log').unlink()
    (worker/'kinbot.log').symlink_to(outside)
    result = audit.capture(root, tmp_path/'export')
    assert result['issues'] > 0
    assert 'PRIVATE OTHER RUN' not in (tmp_path/'export/evidence.txt').read_text()
