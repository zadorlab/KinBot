"""Actual dispatch staging/restart for the bounded example verification panel."""
from copy import deepcopy
import fcntl
import importlib.util
import json
from pathlib import Path
import tarfile

from ase.units import Hartree, invcm
import pytest

from kinbot.anl import dispatch
from kinbot.anl.results import parse_gaussian_harmonic, parse_result
from kinbot.execution import atomic_json
from tests.test_anl_dispatch import _mock_execution

PANEL = Path(__file__).resolve().parents[1]/'examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114'


def load(name):
    spec = importlib.util.spec_from_file_location('validation_'+name, PANEL/(name+'.py'))
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def site(spec):
    for task in spec['tasks']:
        task['resources'] = dict(cores=8, memory_mb=128000, walltime='1-00:00:00', partition='test')


def spectrum(frequencies):
    zpe = .5*sum(x for x in frequencies if x > 0)*invcm/Hartree
    return (' #P B2PLYP/cc-pVTZ Freq=Numer EmpiricalDispersion=GD3BJ\n\n'
            ' Frequencies -- '+' '.join(map(str, frequencies))+'\n'
            f' Zero-point correction= {zpe:.9f}\nNormal termination of Gaussian\n')


def test_native_saddle_hessian_excludes_imaginary_from_zpe():
    request = dict(kind='gaussian_harmonic', method='B2PLYP', basis='cc-pVTZ', dispersion='GD3BJ',
                   natoms=3, stationary_point='transition_state')
    output = spectrum([-350., 800., 1000.])
    parsed = parse_result(output, request)
    assert parsed['frequencies_cm_inverse'] == [-350., 800., 1000.]
    assert parsed['zpe']['hartree'] == pytest.approx(900*invcm/Hartree, abs=1e-9)
    for invalid in ([50.,800.,1000.],[-350.,-20.,1000.],[-.2,800.,1000.],[-350.,800.]):
        with pytest.raises(ValueError):parse_result(spectrum(invalid), request)
    # The existing VPT2 minimum gate has not been weakened.
    with pytest.raises(ValueError):
        parse_gaussian_harmonic(output, method='B2PLYP', basis='cc-pVTZ', dispersion='GD3BJ', natoms=3)


def test_entire_frozen_panel_validates_and_has_bounded_scope(tmp_path):
    panel = load('campaign'); graphs = panel.graphs(tmp_path)
    assert len(graphs) == 69
    assert sum(len(g['tasks']) for g in graphs.values()) == 220
    for name, graph in graphs.items():
        spec = deepcopy(graph); site(spec); dispatch.validate_spec(spec)
        assert spec['intent']['scientific_approval'] is False
        for task in spec['tasks']:
            text = task['input_template']
            if task['backend'] == 'gaussian':
                assert 'B2PLYP/cc-pVTZ' in text and 'B3LYP' not in text
            if '%oldchk=' in text:
                assert str(tmp_path/name/'tasks/hessian/hessian.chk') in text
                assert 'MaxPoints=' in text and 'Geom=AllCheck' in text
                assert 'GradientOnly,EulerPC,MassWeighted' in text
            assert task['resources']['cores'] == 8
    evidence = json.loads((PANEL/'evidence.json').read_text())
    assert len(evidence['screening_checks']) == 24
    assert len(evidence['screening_points']) == 24
    assert all(c['resolution_staged'] in graphs for c in evidence['coverage_notes'])
    recheck = json.loads((PANEL/'remaining_l2.json').read_text())
    assert len(recheck['cases']) == 44
    assert all(c['mode'] == 'reuse' for c in recheck['cases'])


def test_all_graphs_stage_through_real_dispatcher(tmp_path,monkeypatch):
    panel=load('campaign')
    monkeypatch.setattr(dispatch,'assign_partitions',site)
    monkeypatch.setattr(dispatch,'render_site_setup',lambda _:'# simulated site\n')
    panel.stage(tmp_path/'native')
    assert len(list((tmp_path/'native').glob('*/state.json')))==69
    for path in (tmp_path/'native').glob('*/state.json'):
        _, spec, state=dispatch._load(path.parent)
        for task in spec['tasks']:
            if task['id'] in state['tasks']:
                dispatch._verify_stage_files(path.parent,task,state['tasks'][task['id']])
    panel.stage(tmp_path/'native')


@pytest.fixture
def small_panel(tmp_path, monkeypatch):
    panel = load('campaign'); root = tmp_path/'native'
    point = dict(symbols=['O','H','H'],positions=[[0.,0,0],[1.,0,0],[0,1.,0]],charge=0,multiplicity=1)
    monkeypatch.setattr(panel, 'graphs', lambda r: {'example':panel.native_irc('example',point,r/'example')})
    monkeypatch.setattr(dispatch, 'assign_partitions', site)
    monkeypatch.setattr(dispatch, 'render_site_setup', lambda _: '# simulated site\n')
    return panel, root


def test_interrupted_preparation_preserves_and_recovers_without_submitting(small_panel, monkeypatch):
    panel, root = small_panel
    prepare = dispatch.prepare
    def fail(spec, directory):
        directory.mkdir(); (directory/'workflow.json').write_text(json.dumps(spec))
        raise OSError('simulated interrupted write')
    monkeypatch.setattr(dispatch, 'prepare', fail)
    with pytest.raises(OSError):panel.stage(root)
    monkeypatch.setattr(dispatch, 'prepare', prepare)
    panel.stage(root)
    assert len(list((root/'interrupted_preparations').glob('*/workflow.json'))) == 1
    _, _, state = dispatch._load(root/'example')
    assert state['tasks']['hessian']['status'] == 'staged'
    # Also recover publication of state before publication of requested.json.
    (root/'example/requested.json').unlink()
    panel.stage(root)
    assert (root/'example/requested.json').is_file()


def test_full_dispatch_adopts_running_work_and_repeats_without_retry(small_panel, monkeypatch):
    panel, root = small_panel; panel.stage(root)
    graph = root/'example'
    _, _, state = dispatch._load(graph)
    state['tasks']['hessian'].update(status='submitted', job_id='123')
    atomic_json(graph/'state.json',state)
    monkeypatch.setattr(dispatch, '_job_active', lambda _:True)
    real_advance = dispatch.advance; calls = []
    def tick(path, submit):
        calls.append(str(path))
        # Exercise actual reconciliation, hashes and dependent task staging;
        # submission itself is represented by the live allocation above.
        return real_advance(path, submit=False)
    monkeypatch.setattr(dispatch, 'advance', tick)
    class Interrupted(Exception):pass
    monkeypatch.setattr(panel.time, 'sleep', lambda _:(_ for _ in ()).throw(Interrupted()))
    with pytest.raises(Interrupted):panel.run(root)
    assert dispatch._load(graph)[2]['tasks']['hessian']['job_id'] == '123'
    checkpoint = graph/'tasks/hessian/hessian.chk'
    _mock_execution(graph,'hessian')
    execution = graph/'tasks/hessian/execution.json'; e=json.loads(execution.read_text())
    e['details']['parsed_result'] = parse_gaussian_harmonic(spectrum([-350.,800.,1000.]),
        method='B2PLYP',basis='cc-pVTZ',dispersion='GD3BJ',natoms=3,stationary_point='transition_state')
    atomic_json(execution,e)
    # An output file arriving before the allocation exits must not be adopted.
    state = real_advance(graph)
    assert len(state['tasks']) == 1
    monkeypatch.setattr(dispatch, '_job_active', lambda _:False)
    state = real_advance(graph)
    assert len(state['tasks']) == 5
    before = checkpoint.read_bytes()
    for ident in ('forward_10','reverse_10','forward_5','reverse_5'):
        _mock_execution(graph,ident)
    panel.run(root); panel.run(root)
    assert checkpoint.read_bytes() == before
    assert all(e['status']=='complete' for e in dispatch._load(graph)[2]['tasks'].values())
    completion=json.loads((root/'completion.json').read_text())
    assert not completion['scientific_approval'] and not completion['mess_updated']


def test_second_controller_cannot_touch_owned_graph(small_panel, monkeypatch):
    panel,root=small_panel; panel.stage(root)
    before=(root/'example/state.json').read_bytes()
    monkeypatch.setattr(dispatch,'advance',lambda *a,**k:pytest.fail('duplicate controller'))
    with (root/'controller.lock').open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX | fcntl.LOCK_NB)
        panel.run(root)
        with pytest.raises(BlockingIOError):panel.stage(root)
    assert (root/'example/state.json').read_bytes()==before


def test_scheduler_outage_and_changed_plan_preserve_native_state(small_panel,monkeypatch):
    panel,root=small_panel;panel.stage(root)
    before=(root/'example/state.json').read_bytes()
    monkeypatch.setattr(dispatch,'advance',lambda *a,**k:(_ for _ in ()).throw(dispatch.SchedulerUnavailableError('unavailable')))
    class Interrupted(Exception):pass
    monkeypatch.setattr(panel.time,'sleep',lambda _:(_ for _ in ()).throw(Interrupted()))
    with pytest.raises(Interrupted):panel.run(root)
    progress=json.loads((root/'progress.json').read_text())
    assert progress['infrastructure']=={'example':'unavailable'}
    assert not progress['execution_finished']
    assert (root/'example/state.json').read_bytes()==before
    monkeypatch.setattr(panel,'fingerprint',lambda:{'different':'definition'})
    with pytest.raises(ValueError,match='definitions changed'):panel.stage(root)
    assert (root/'example/state.json').read_bytes()==before


def test_failed_hessian_blocks_paths_without_new_attempt(small_panel, monkeypatch):
    panel,root=small_panel; panel.stage(root); graph=root/'example'
    _,_,state=dispatch._load(graph); state['tasks']['hessian']['status']='failed'
    atomic_json(graph/'state.json',state)
    real=dispatch.advance
    monkeypatch.setattr(dispatch,'advance',lambda path,submit:real(path,submit=False))
    panel.run(root); before=(graph/'state.json').read_bytes();panel.run(root)
    assert (graph/'state.json').read_bytes()==before
    assert {e['status'] for e in dispatch._load(graph)[2]['tasks'].values()}=={'failed','blocked'}


def test_collector_does_not_walk_checkpoint_trees(tmp_path):
    collector=load('collect');root=tmp_path/'validation';root.mkdir()
    work=root/'remaining_l2/test';work.mkdir(parents=True)
    (root/'remaining_l2/campaign.json').write_text(json.dumps({'cases':[{'name':'test'}]}))
    (work/'source.log').write_text('Full native output\nNormal termination\n')
    ignored=work/'.kinbot_checkpoints/a/b';ignored.mkdir(parents=True)
    (ignored/'large.log').write_text('not selected')
    (work/'large.chk').write_bytes(b'not selected')
    out=collector.collect(root,tmp_path/'exports')
    with tarfile.open(out/'validation-evidence.tar.gz') as tar:
        assert set(tar.getnames())=={'remaining_l2/campaign.json','remaining_l2/test/source.log'}
    assert json.loads((out/'manifest.json').read_text())['issues']==[]
