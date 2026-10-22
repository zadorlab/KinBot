from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.io import write

from kinbot import vrc_branch_probe as probe, run_policy
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from kinbot.vrc_scan_recovery import digest


@pytest.fixture
def source(tmp_path):
    well = tmp_path/'ch3s/123'
    graph = well/'vrctst/molpro_scan/channel_vts_fragments_asymptote_original'
    directory = graph/'tasks/relaxed_scan'
    # Make the mocked checkpoint before staging the native contract.
    checkpoint = tmp_path/'seed.wfu'; checkpoint.write_bytes(b'initial orbitals')
    model = dict(state_average=[dict(multiplicity=m, roots=3, weights=[1/3]*3) for m in (4,2)],
        selected_roots={'2':1,'4':1}, canonical_orbitals=True,
        orbital_selection=dict(mode='molpro_record', path=str(checkpoint),
            sha256=file_sha256(checkpoint), record='2140.2', orbital_type='natural',
            active_indices=[12,13,14,15], align_active_orbitals=False))
    atoms = Atoms(['C','S','H','H','H'], positions=[[0,0,0],[30,0,0],[31.33,0,0],
                                                               [30,1.48,0],[30,-1.48,0]])
    molecule = dict(symbols=atoms.get_chemical_symbols(), positions=atoms.positions.tolist(),
                    charge=0, multiplicity=2)
    task = dict(id='relaxed_scan', kind='ase_optimize', geometry_from='initial', geometry_output='final.xyz',
        profile=dict(calculator='molpro', method='caspt2-f12(3,4)', basis='avdz-f12',
            command='molpro', optimizer='sella', frequency_mode='ase_forces',
            calculator_kwargs=dict(state_model=model, caspt2_shift=.1, caspt2_ipea=.25,
                gradient_mode='numerical', scf_convergence=1.e-8, correlation_convergence=1.e-7)),
        optimizer=dict(steps=100, fmax=.0005, order=0, sella_kwargs={'internal':False},
                       constraints=dict(fixed_bonds=[[0,1]], minimum_bonds=[])),
        resources=dict(cores=8, memory_mb=104000, min_stack_mw=1000, partition='cpu',
                       walltime='31-00:00:00', requested_walltime='partition_max'))
    spec = dict(schema=1, name='asymptote', molecule=molecule, limits={'max_nodes':1}, tasks=[task])
    dispatch.prepare(spec, graph)
    checkpoint = directory/'final.wfu'; checkpoint.write_bytes(b'final orbitals')
    model = deepcopy(model)
    model['orbital_selection'].update(path=str(checkpoint), sha256=file_sha256(checkpoint))
    write(directory/'final.xyz', atoms)
    execution = dict(schema=1, status='executed', details={'orbital_checkpoint_model':model},
        artifacts={p.name:file_sha256(p) for p in [directory/'task.json', directory/'final.xyz', checkpoint]})
    atomic_json(directory/'execution.json', execution)
    state = probe.read(graph/'state.json'); state['tasks']['relaxed_scan'].update(status='complete',job_id='101')
    atomic_json(graph/'state.json', state)
    context = dict(symbols=molecule['symbols'], charge=0, multiplicity=2,
        maps=[[0],[1,2,3,4]], coordinate=[0,1], equiv=[[0],[1]],
        fragment_geometry=atoms.positions.tolist(), separation_angstrom=30., state_model=model,
        parameters=dict(vrc_tst_scan_method='caspt2-f12(3,4)', vrc_tst_scan_basis='avdz-f12',
            vrc_tst_scan_shift=.1, vrc_tst_scan_ipea=.25, vrc_tst_scan_gradient='numerical',
            vrc_tst_scf_convergence=1.e-8, vrc_tst_correlation_convergence=1.e-7, vrc_tst_scan_fmax=.0005))
    review = well/'vrctst/scan_review/channel/scan_review.json'
    atomic_json(review, dict(reaction='channel', context=context, context_sha256=digest(context)))
    atomic_json(tmp_path/run_policy.NAME, dict(schema=1, max_jobs=23, excluded_nodes=[]))
    return tmp_path, graph, review


def test_probe_preserves_science_and_internal_geometry_without_touching_production(source):
    root, original, review = source
    before = {p:file_sha256(p) for p in (root/'ch3s').rglob('*') if p.is_file()}
    graph = probe.prepare(review, 6.)
    assert graph.parent == root/'diagnostics/vrc_branch_probes'
    spec = probe.read(graph/'workflow.json')
    task = spec['tasks'][0]
    original_task = probe.read(original/'tasks/relaxed_scan/task.json')['task']
    model = probe.read(review)['context']['state_model']
    expected_profile = deepcopy(original_task['profile'])
    expected_profile['calculator_kwargs']['state_model'] = model
    assert task['profile'] == expected_profile
    assert task['optimizer']['fmax'] == .0005
    assert task['optimizer']['steps'] == 50 and task['optimizer']['step_retries'] == 0
    assert task['resources']['cores'] == 8 and task['resources']['walltime'] == '01:00:00'
    assert 'requested_walltime' not in task['resources']
    assert spec['intent']['vrc_fragment_retry'] and not spec['intent']['production_model_accepted']
    x = np.array(spec['molecule']['positions'])
    assert np.linalg.norm(x[0]-x[1]) == pytest.approx(6.)
    assert np.linalg.norm(x[2]-x[1]) == pytest.approx(1.33)
    assert np.linalg.norm(x[3]-x[1]) == pytest.approx(1.48)
    assert run_policy.find(graph) == root/run_policy.NAME
    assert before == {p:file_sha256(p) for p in before}
    # Re-running after a terminal failure must not allocate a replacement.
    state = probe.read(graph/'state.json')
    state['tasks']['relaxed_scan'].update(status='failed',job_id='102')
    atomic_json(graph/'state.json', state)
    assert probe.prepare(review, 6.) == graph
    assert probe.status(graph)['task']['job_id'] == '102'
    assert len(list(graph.parent.glob('probe_*/workflow.json'))) == 1


@pytest.mark.parametrize('change,match', [('checkpoint','checkpoint'), ('final','artifact'),
    ('task','artifact'), ('context','context hash'), ('rejected','checkpoint'),
    ('active','not complete'), ('method','method/basis'), ('threshold','numerical contract'),
    ('map','partition'), ('state','electronic-state')])
def test_probe_refuses_changed_or_unverified_source_before_staging(source, change, match):
    root, graph, review = source
    directory = graph/'tasks/relaxed_scan'
    if change == 'checkpoint': (directory/'final.wfu').write_bytes(b'changed')
    elif change == 'final': (directory/'final.xyz').write_text('changed')
    elif change == 'task': (directory/'task.json').write_text('{}')
    elif change == 'rejected': atomic_json(directory/'scan_rejection.json', {})
    elif change == 'active':
        s = probe.read(graph/'state.json'); s['tasks']['relaxed_scan']['status']='submitted'
        atomic_json(graph/'state.json',s)
    else:
        s = probe.read(review)
        if change in ('context','method'): s['context']['parameters']['vrc_tst_scan_method']='caspt2-f12(5,4)'
        if change == 'threshold': s['context']['parameters']['vrc_tst_correlation_convergence']=1.e-9
        if change == 'map': s['context']['maps']=[[0],[1,2,3,3]]
        if change == 'state': s['context']['state_model']['selected_roots']['2']=2
        if change != 'context': s['context_sha256']=digest(s['context'])
        atomic_json(review,s)
    with pytest.raises((ValueError,RuntimeError), match=match): probe.prepare(review, 6.)
    assert not (root/'diagnostics').exists()


@pytest.mark.parametrize('distance,options', [(30.,{}), (float('nan'),{}), (-1.,{}),
    (6.,{'steps':101}), (6.,{'walltime':'03:00:00'})])
def test_probe_budget_is_bounded(source, distance, options):
    with pytest.raises(ValueError, match='interior distance'):
        probe.prepare(source[2], distance, **options)


def test_status_reports_transfer_without_accepting_production(source, capsys):
    root, original, review = source
    graph = probe.prepare(review, 6.)
    directory = graph/'tasks/relaxed_scan'
    atoms = Atoms(['C','S','H','H','H'], positions=[[0,0,0],[6,0,0],[7.33,0,0],
                                                               [1.11,0,0],[6,-1.48,0]])
    write(directory/'final.xyz', atoms)
    atomic_json(directory/'execution.json', dict(status='executed',
        artifacts={'final.xyz':file_sha256(directory/'final.xyz')}, details={'energy_hartree':-437.}))
    before = review.read_bytes()
    report = probe.status(graph)
    distances = {tuple(r['atoms']):r['final'] for r in report['pair_distances_angstrom']}
    assert distances[(1,4)] == pytest.approx(1.11)
    assert distances[(2,4)] == pytest.approx(4.89)
    assert report['production_model_accepted'] is False
    assert probe.main(['status',str(root)]) == 0
    assert '1.11' in capsys.readouterr().out
    assert review.read_bytes() == before
