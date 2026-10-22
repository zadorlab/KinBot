import json

from kinbot.execution import atomic_json
from kinbot.native_failure_report import snapshot, excerpt


def test_archived_abort_report_is_read_only_and_separate_from_current_job(tmp_path):
    graph = tmp_path/'graph'
    task = graph/'tasks/relaxed_scan'
    atomic_json(graph/'workflow.json', {})
    atomic_json(graph/'state.json', {'tasks':{'relaxed_scan':{'status':'submitted','attempt':3,'job_id':'200'}}})
    atomic_json(task/'task.json', {'task':{'resources':{'cores':1}}})
    archive = graph/'attempts/live_retry/graphs/tasks/relaxed_scan'
    atomic_json(archive/'step.failure.json', {'returncode':255,'host':'node31',
        'files':{'scratch/rank.log_5':{'tail':'MPI_Abort from rank 5'}}})
    atomic_json(archive/'execution.json', {'status':'failed','slurm_job_id':'100'})
    (archive/'step.log').write_text('Error in displaced calculation\n'+'bookkeeping\n'*10000)
    before = {p:p.read_bytes() for p in graph.rglob('*') if p.is_file()}
    report = snapshot(task, 'step')
    assert report['current']['job_id'] == '200'
    assert report['resources']['cores'] == 1
    assert report['attempts'][0]['execution']['slurm_job_id'] == '100'
    assert 'Error in displaced' in report['attempts'][0]['files']['step.log']['error_context']
    assert report['attempts'][0]['files']['scratch/rank.log_5']['tail']=='MPI_Abort from rank 5'
    assert before == {p:p.read_bytes() for p in graph.rglob('*') if p.is_file()}


def test_excerpt_declares_truncation(tmp_path):
    path = tmp_path/'native.out'
    path.write_text('too old\n' + 'ERROR\n'*200)
    evidence = excerpt(path, limit=1000)
    assert evidence['scan_truncated'] and evidence['error_context_truncated']
    assert evidence['scan_start_byte'] == path.stat().st_size-1000
