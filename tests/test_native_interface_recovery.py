"""Known interface repairs preserve scientific inputs, live work and lineage."""
import json
from pathlib import Path

import pytest

from kinbot.anl import dispatch
from kinbot import run_recover
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution
from tests.test_molpro_rank_exact_interfaces import H2_ENERGY


def fail(graph, kind, *, job='42', status='failed'):
    d = graph/'tasks/bad'
    state = dispatch._load(graph)[2]
    entry = state['tasks']['bad']
    entry.update(status=status, job_id=job)
    if kind == 'cfour':
        zmat = d/'ZMAT'
        zmat.write_text(zmat.read_text().replace('AUG-PVTZ', 'aug-cc-pVTZ'))
        entry['input_sha256'] = file_sha256(zmat)
        (d/'cfour.out').write_text('@GTFLGS-W, Option aug-cc-pVTZ for keyword BASIS\n'
                                  'not known. Default used.\n@MKDECK-F, I/O error on file MOL.\n')
        error = 'Success marker absent from cfour.out.'
    else:
        (d/'bad_step_0001.out').write_text(H2_ENERGY +
            ' Numerical gradient for KB_GEOM_ENERGY\n'
            ' Atom dE/dx dE/dy dE/dz d2E/dx2 d2E/dy2 d2E/dz2\n'
            ' 1 0.0 0.0 0.01 0 0 0\n 2 0.0 0.0 -0.01 0 0 0\n'
            ' Molpro calculation terminated\n')
        (d/'bad_step_0001.xyz').write_bytes((d/'geometry.xyz').read_bytes())
        error = 'Molpro output lacks the requested undisplaced F12b energy.'
    atomic_json(d/'execution.json', dict(schema=1, task_id='bad', status='failed',
        geometry_sha256=entry['geometry_sha256'], error=error, slurm_job_id=job,
        failure_phase='validation' if kind == 'cfour' else 'execution'))
    atomic_json(graph/'state.json', state)


@pytest.fixture(params=['cfour', 'molpro'])
def graph(tmp_path, monkeypatch, request):
    kind = request.param
    spec = dispatch_spec()
    spec['molecule'] = dict(symbols=['H', 'H'], positions=[[0, 0, 0], [0, 0, .74]],
                            charge=0, multiplicity=1)
    bad = next(t for t in spec['tasks'] if t['id'] ==
               ('cfour_dboc' if kind == 'cfour' else 'l3_geometry'))
    bad.update(id='bad', geometry_from='initial')
    if kind == 'cfour':
        bad['input_template'] = bad['input_template'].replace('cc-pVTZ', 'aug-cc-pVTZ')
        bad['result_parser']['basis'] = 'aug-cc-pVTZ'
    else:
        bad['profile'].update(method='CCSD(T)-F12b', basis='aug-cc-pVTZ-F12',
                              calculator_kwargs={'effective_method': 'CCSD-F12b'})
    spec['tasks'] = [bad] + [molpro_task(n, 'rhf\n', geometry_from='initial')
                             for n in ('live', 'done', 'child')]
    spec['tasks'][-1]['depends_on'] = ['bad']
    root = dispatch.prepare(spec, tmp_path/'graph')
    fail(root, kind)
    _mock_execution(root, 'done')
    state = dispatch._load(root)[2]
    state['tasks']['live'].update(status='submitted', job_id='43')
    state['tasks']['child'] = dict(status='blocked', blocked_by=['bad'])
    atomic_json(root/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda job: str(job) == '43')
    return root, kind


def test_repair_once_preserves_work_and_bound_scientific_spec(graph):
    root, kind = graph
    original = (root/'workflow.json').read_bytes()
    unchanged = {p: file_sha256(p) for name in ('live', 'done')
                 for p in (root/'tasks'/name).iterdir() if p.is_file()}
    state = dispatch.advance(root)
    repaired = state['tasks']['bad']
    assert repaired['status'] == 'staged' and repaired['attempt'] == 2
    assert (Path(repaired['previous_attempt'])/'execution.json').is_file()
    assert state['tasks']['live']['job_id'] == '43'
    assert state['tasks']['live']['status'] == 'submitted'
    assert state['tasks']['done']['status'] == 'complete'
    assert (root/'workflow.json').read_bytes() == original
    assert unchanged == {p: file_sha256(p) for p in unchanged}
    if kind == 'cfour':
        assert 'BASIS=AUG-PVTZ' in (root/'tasks/bad/ZMAT').read_text()
    assert dispatch.advance(root) == state
    fail(root, kind, job='44')
    state = dispatch.advance(root)
    assert state['tasks']['bad']['status'] == 'failed'
    assert state['tasks']['bad']['attempt'] == 2


def test_allocated_intermediate_failure_is_untouched(graph, monkeypatch):
    root, kind = graph
    fail(root, kind, status='submitted')
    monkeypatch.setattr(dispatch, '_job_active', lambda job: True)
    assert dispatch.advance(root)['tasks']['bad']['status'] == 'submitted'
    assert not (root/'.live_control_transaction.json').exists()


@pytest.mark.parametrize('defect', ['wrong_provenance', 'native_failure', 'downstream_complete'])
def test_repair_does_not_mask_real_errors_or_invalidate_dependents(graph, defect):
    root, kind = graph
    p = root/'tasks/bad/execution.json'
    if defect == 'wrong_provenance':
        value = json.loads(p.read_text()); value['geometry_sha256'] = 'wrong'
        atomic_json(p, value)
        with pytest.raises(RuntimeError, match='execution record does not match'):
            dispatch.advance(root)
    else:
        if defect == 'native_failure':
            native = root/'tasks/bad'/('cfour.out' if kind == 'cfour' else 'bad_step_0001.out')
            native.write_text('SCF did not converge.\nERROR EXIT\n')
        else:
            state = dispatch._load(root)[2]
            state['tasks']['child']['status'] = 'complete'
            # Exercise recovery directly; the fabricated descendant has no
            # accepted output and must never be invalidated by the repair.
            from kinbot.anl.interface_recovery import recover
            assert recover(root, dispatch._load(root)[1], state) == state
            return
        assert dispatch.advance(root)['tasks']['bad']['status'] == 'failed'
    assert not (root/'.live_control_transaction.json').exists()


def test_interrupted_archival_resumes_without_double_retry(graph, monkeypatch):
    root, _ = graph
    original = run_recover.atomic_json
    def interrupted(path, data):
        if Path(path) == root/'state.json':
            raise OSError('interrupted transaction')
        return original(path, data)
    monkeypatch.setattr(run_recover, 'atomic_json', interrupted)
    with pytest.raises(OSError, match='interrupted transaction'):
        dispatch.advance(root)
    monkeypatch.setattr(run_recover, 'atomic_json', original)
    state = dispatch.advance(root)
    assert state['tasks']['bad']['attempt'] == 2
    assert state['tasks']['live']['job_id'] == '43'
    assert len(list((root/'attempts').glob('live_*'))) == 1


def test_old_unsubmitted_cfour_input_is_repaired_before_allocation(graph):
    root, kind = graph
    if kind != 'cfour':
        return
    state = dispatch._load(root)[2]
    state['tasks']['bad']['status'] = 'staged'
    state['tasks']['bad'].pop('job_id')
    (root/'tasks/bad/execution.json').unlink()
    atomic_json(root/'state.json', state)
    state = dispatch.advance(root)
    assert state['tasks']['bad']['status'] == 'staged'
    assert state['tasks']['bad']['attempt'] == 1
    assert 'BASIS=AUG-PVTZ' in (root/'tasks/bad/ZMAT').read_text()


@pytest.mark.parametrize('evidence', ['valid', 'active', 'queue_unknown', 'wrong_job', 'changed_input', 'wrong_owner'])
def test_explicit_retry_with_purged_scheduler_record_requires_verified_quiescence(graph, monkeypatch, evidence):
    from kinbot.anl import live_control as owner
    root, _ = graph
    directory = root/'tasks/bad'
    monkeypatch.setattr(owner, 'slurm_job_record', lambda job: {})
    def active(job):
        if evidence == 'queue_unknown':
            raise dispatch.SchedulerUnavailableError('Slurm unreachable')
        return evidence == 'active'
    monkeypatch.setattr(dispatch, '_job_active', active)
    if evidence == 'wrong_job':
        p = directory/'execution.json'; record = json.loads(p.read_text())
        record['slurm_job_id'] = '99'; atomic_json(p, record)
    if evidence == 'changed_input':
        from ase.io import read, write
        atoms = read(directory/'geometry.xyz')
        atoms.positions[0, 0] += .1
        write(directory/'geometry.xyz', atoms)
    if evidence == 'wrong_owner':
        monkeypatch.setattr(owner, 'slurm_job_record', lambda job:
            dict(JobId='42', JobState='FAILED', WorkDir=str(root/'elsewhere')))
    if evidence == 'valid':
        assert owner.stop_attempt(directory, '42', root/'unused.json', {})
    else:
        with pytest.raises((RuntimeError, ValueError)):
            owner.stop_attempt(directory, '42', root/'unused.json', {})
