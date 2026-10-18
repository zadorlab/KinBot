"""The example transition must preserve results and refuse duplicate launches."""
import importlib.util
import json
import os
from pathlib import Path

import pytest

PACKAGE = Path(__file__).parents[1]/'examples/anl/tools'


@pytest.fixture
def handoff(monkeypatch):
    monkeypatch.syspath_prepend(str(PACKAGE))
    spec = importlib.util.spec_from_file_location('caspt2_handoff', PACKAGE/'run_examples.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_prepared_package_is_portable_immutable_and_records_resolved_settings(handoff, tmp_path):
    root = handoff.prepare(tmp_path/'production')
    assert handoff.prepare(root) == root
    record = handoff.verify(root)
    for case in handoff.CASES:
        resolved = json.loads((root/case/'resolved_parameters.json').read_text())
        assert resolved['epsilon_unit'] == 'cm-1'
        assert resolved['queue_job_limit'] == 10
        assert resolved['vrc_tst_surface_model_file'] == ''
        assert resolved['vrc_tst_surface_model']['channels']
        assert resolved['rotdpy_flux_parameters']['flux_rel_err'] == 1
        assert resolved['sella_fmax'] == resolved['anl_geometry_fmax'] == 5e-4
        assert resolved['vrc_tst_scan_fmax'] == 5e-4
        assert resolved['anl_scf_convergence'] == resolved['vrc_tst_scf_convergence'] == 1e-8
        assert resolved['anl_cc_convergence'] == resolved['vrc_tst_correlation_convergence'] == 1e-7
        assert resolved['barrier_threshold_L2'] == 400
        assert resolved['barrier_threshold'] == 410
        assert max(resolved['TemperatureList']) == 5000
        # Match the actual RotdPy grid generator used in the rendered worker.
        def grid(values):
            start, interval, factor, count = values
            points = [start]
            for _ in range(1, count):
                start += interval
                points.append(start)
                interval *= factor
            return points
        assert grid(resolved['rotdpy_temperature_grid']) == resolved['TemperatureList']
        gas_constant = .00198720425864083  # kcal mol-1 K-1
        assert resolved['ModelEnergyLimit'] > (resolved['barrier_threshold'] +
            resolved['ExcessEnergyOverTemperature'] * gas_constant * 5000)
        assert grid(resolved['rotdpy_energy_grid'])[-1] * gas_constant > resolved['ModelEnergyLimit']
        assert f'{case}/resolved_parameters.json' in record['sha256']
    model = root/'ch3s/input.json'
    model.write_text('{}')
    with pytest.raises(ValueError, match='Prepared artifact changed'):
        handoff.verify(root)
    with pytest.raises(ValueError, match='Existing prepared settings differ'):
        handoff.prepare(root)


def test_retirement_inventories_first_and_preserves_unrelated_runs(handoff, monkeypatch, tmp_path):
    from kinbot import run_cleanup, orbital_campaign
    root = tmp_path/'production'
    root.mkdir()
    parent = tmp_path/'old'
    parent.mkdir()
    expected = [parent/name for name in handoff.OLD_NAMES[:2]]
    for old, marker in zip(expected, ('campaign.json', 'confirmation.json')):
        old.mkdir()
        (old/marker).write_text('{}')
        (old/'native.out').write_text('retained evidence')
    unrelated = parent/'unrelated_running_calculation'
    unrelated.mkdir()
    (unrelated/'native.out').write_text('unrelated evidence')
    events = []
    monkeypatch.setattr(run_cleanup, 'inspect', lambda path: events.append(('inspect', path)) or {})
    monkeypatch.setattr(orbital_campaign, 'stop_controller', lambda path: events.append(('controller', path)))

    def stop(path, *, archive, timeout):
        assert archive and timeout == 300
        events.append(('stop', path))
        destination = path.with_name(path.name+'.stopped.test')
        path.rename(destination)
        return {'archive': str(destination)}

    monkeypatch.setattr(run_cleanup, 'stop', stop)
    reports = handoff.retire(root, parent)
    assert events[:2] == [('inspect', path) for path in expected]
    assert [path for action, path in events if action == 'stop'] == expected
    assert (unrelated/'native.out').read_text() == 'unrelated evidence'
    for report in reports:
        assert (Path(report['archive'])/'native.out').read_text() == 'retained evidence'
    assert len(list(root.glob('retirement_inventory.*.json'))) == 1
    assert len(list(root.glob('retirement_results.*.json'))) == 1


def test_unrecognized_old_run_aborts_before_any_stop(handoff, monkeypatch, tmp_path):
    from kinbot import run_cleanup
    root = tmp_path/'production'
    root.mkdir()
    bad = tmp_path/handoff.OLD_NAMES[0]
    bad.mkdir()
    monkeypatch.setattr(run_cleanup, 'stop', lambda *a, **k: pytest.fail('must not stop'))
    with pytest.raises(ValueError, match='Unrecognized old campaign'):
        handoff.retire(root, tmp_path)


def test_previous_launch_refuses_duplicate_even_if_queue_empty(handoff, monkeypatch, tmp_path):
    from kinbot import run_cleanup
    case = tmp_path/'ch3s'
    case.mkdir()
    (case/'launch.json').write_text('{}')
    monkeypatch.setattr(handoff, 'verify', lambda root: {'inputs': [str(case/'input.json')]})
    monkeypatch.setattr(run_cleanup, 'inspect', lambda path: {'jobs': [], 'processes': []})
    with pytest.raises(RuntimeError, match='already started'):
        handoff.assert_not_started(tmp_path)


def test_rendered_preflight_uses_both_real_production_profiles(handoff, tmp_path, monkeypatch):
    monkeypatch.setenv('SBATCH_EXCLUDE', 'existing-exclusion')
    root = handoff.prepare(tmp_path/'production')
    result = handoff.render_gaussian_inputs(root, handoff.verify(root))
    assert os.environ['SBATCH_EXCLUDE'] == 'existing-exclusion'
    assert len(result['files']) == 72
    for item in result['files']:
        text = (root/item['path']).read_text()
        l2 = Path(item['path']).name.startswith('l2_')
        assert ('B2PLYP/cc-pVTZ' if l2 else 'B3LYP/cc-pVDZ') in text
        assert 'scf(xqc,conver=8)' in text
        assert 'scf(xqc,(' not in text
        assert 'EmpiricalDispersion(GD3BJ)' in text
        assert 'integral(UltraFine)' in text
        assert '\nopt(' not in text and '\nirc(' not in text
        if item['path'].endswith('_energy.com'):
            assert '\nfreq' not in text
    assert not (root/'ch3s/launch.json').exists()


def test_current_inputs_propagate_thresholds_into_reference_and_all_l3_backends(handoff, tmp_path):
    from ase.build import molecule
    from kinbot.parameters import Parameters
    from kinbot.theory import resolve_profiles
    from kinbot.anl.automatic import _reference_input
    from kinbot.anl.validation import automatic_composite_spec
    from kinbot.ase_modules.calculators.molpro import render_input
    root = handoff.prepare(tmp_path/'production')
    for case in handoff.CASES:
        par = Parameters(str(root/case/'input.json')).par
        directory = root/case/'reference'
        directory.mkdir()
        ref = _reference_input(root/case/'input.json',
            {'smiles': '[SH]', 'charge': 0, 'multiplicity': 2},
            directory, directory/'thermochemistry_model.json')
        ref_par = Parameters(str(ref)).par
        assert ref_par['anl_scf_convergence'] == 1e-8
        assert ref_par['sella_fmax'] == 5e-4
        for profile in resolve_profiles(ref_par).values():
            assert profile.calculator_kwargs['scf'] == 'xqc,conver=8'
        atoms = molecule('CH4')
        spec = automatic_composite_spec({
            'symbols': atoms.get_chemical_symbols(), 'positions': atoms.positions.tolist(),
            'charge': 0, 'multiplicity': 1,
            'source': {'geometry_sha256': '1'*64, 'record_sha256': '2'*64,
                       'energy_hartree': -40., 'harmonic_zpe_hartree': .04}},
            highest_tier=par['composite_method'], scf_convergence=par['anl_scf_convergence'],
            cc_convergence=par['anl_cc_convergence'], geometry_fmax=par['anl_geometry_fmax'])
        for task in spec['tasks']:
            if task.get('kind') == 'ase_optimize':
                assert task['optimizer']['fmax'] == 5e-4
                kw = task['profile']['calculator_kwargs']
                assert (kw['scf_convergence'], kw['correlation_convergence']) == (1e-8, 1e-7)
                deck = render_input(atoms, basis='cc-pVTZ', geometry_name='step.xyz',
                    scf_convergence=kw['scf_convergence'], correlation_convergence=kw['correlation_convergence'])
                assert 'gthresh,energy=1.d-8,' in deck
                assert 'thrden=1.d-7,thrvar=1.d-10' in deck
            template = task.get('input_template', '')
            if task.get('backend') == 'gaussian':
                assert 'SCF=(XQC,Conver=8)' in template
            elif task.get('backend') == 'molpro':
                assert 'gthresh,energy=1.d-8,' in template
                if ',thrden=' in template:
                    assert 'thrden=1.d-7,thrvar=1.d-10' in template
            elif task.get('backend') == 'mrcc':
                assert 'sctol=8\n' in template and 'cctol=7\n' in template
            elif task.get('backend') == 'cfour':
                assert 'SCF_CONV=8\n' in template
                if task['id'] == 'cfour_dboc':
                    assert 'CALC=SCF\n' in template
                else:
                    assert 'CC_CONV=7\n' in template


def test_profile_threshold_overrides_remain_available(handoff, tmp_path):
    from kinbot.parameters import Parameters
    from kinbot.theory import resolve_profiles
    root = handoff.prepare(tmp_path/'production')
    par = Parameters(str(root/'ch3s/input.json')).par
    par['anl_scf_convergence'] = 1e-9
    assert all(p.calculator_kwargs['scf'] == 'xqc,conver=9'
               for p in resolve_profiles(par).values())
    par['l1_profile'] = {'calculator_kwargs': {'scf': 'diis,conver=7'}}
    profiles = resolve_profiles(par)
    assert profiles['l1'].calculator_kwargs['scf'] == 'diis,conver=7'
    assert profiles['l2'].calculator_kwargs['scf'] == 'xqc,conver=9'


@pytest.mark.parametrize('stop_fails', [False, True])
def test_clean_restart_stops_both_before_archiving_parent(handoff, monkeypatch, tmp_path, stop_fails):
    from kinbot import run_cleanup
    old = handoff.prepare(tmp_path/'previous')
    root = handoff.prepare(tmp_path/'replacement')
    (old/'ch3s/evidence.log').write_text('failed native input evidence')
    events = []
    def inspect(path):
        events.append(('inspect', path.name))
        return {'directory': str(path), 'jobs': [], 'processes': []}
    def stop(path, **kwargs):
        events.append(('stop', path.name))
        assert kwargs == {'archive': False, 'timeout': 300.}
        if stop_fails and path.name == 'ch3ss':
            raise RuntimeError('scheduler unavailable')
        assert old.exists()
        return {'status': 'stopped'}
    monkeypatch.setattr(run_cleanup, 'inspect', inspect)
    monkeypatch.setattr(run_cleanup, 'stop', stop)
    if stop_fails:
        with pytest.raises(RuntimeError, match='scheduler unavailable'):
            handoff.retire_production(root, old)
        assert old.exists()
        assert not (root/'ch3s/launch.json').exists()
    else:
        result = handoff.retire_production(root, old)
        assert not old.exists()
        assert (Path(result['archive'])/'ch3s/evidence.log').read_text() == 'failed native input evidence'
        assert events == [('inspect', 'ch3s'), ('inspect', 'ch3ss'),
                          ('stop', 'ch3s'), ('stop', 'ch3ss'),
                          ('inspect', 'ch3s'), ('inspect', 'ch3ss')]


def test_cleanup_refuses_unrecognized_or_overlapping_production(handoff, monkeypatch, tmp_path):
    from kinbot import run_cleanup
    old = handoff.prepare(tmp_path/'previous')
    root = handoff.prepare(tmp_path/'replacement')
    monkeypatch.setattr(run_cleanup, 'stop', lambda *a, **k: pytest.fail('must not stop'))
    with pytest.raises(ValueError, match='separate directories'):
        handoff.retire_production(old, old)
    data = json.loads((old/'prepared.json').read_text())
    data['inputs'][0] = '/unrelated/input.json'
    (old/'prepared.json').write_text(json.dumps(data))
    with pytest.raises(ValueError, match='Unrecognized previous'):
        handoff.retire_production(root, old)


def test_inline_vrc_preflight_renders_all_reviewed_channels_and_roots(handoff, tmp_path):
    root = handoff.prepare(tmp_path/'production')
    result = handoff.render_vrc_inputs(root, handoff.verify(root))
    assert len(result['files']) == 20  # Four physical surfaces, five MR components each.
    for entry in result['files']:
        assert entry['method'] in ('caspt2-f12(5,4)', 'caspt2-f12(7,5)')
        text = (root/entry['path']).read_text()
        assert 'shift=0.1,ipea=0.25' in text and 'refstate,' in text
        assert 'rotate,' not in text and 'diab,' not in text
        assert entry['active_indices'] == (list(range(11,15)) if 'ch3s/' in entry['path'] else list(range(18,23)))
    assert not list((root/'assets').glob('*state_model.json'))
    assert not list(root.rglob('surfaces.json'))


def test_obsolete_layout_archive_preserves_untracked_data_and_ignores_other_directories(handoff, tmp_path, monkeypatch):
    repo = tmp_path/'repo'
    old = repo/'examples/anl/caspt2_production/orbitals'
    old.mkdir(parents=True)
    (old/'old.wfu').write_bytes(b'retained native data')
    unrelated = repo/'examples/anl/unrelated'
    unrelated.mkdir()
    (unrelated/'input.json').write_text('{}')
    root = tmp_path/'new_run'; root.mkdir()
    monkeypatch.setattr(handoff, 'REPO', repo)
    calls = []
    def untracked(command, **kwargs):
        calls.append(command)
        assert kwargs['cwd'] == repo
        assert command[:4] == ['git', 'ls-files', '-z', '--']
        return b''
    monkeypatch.setattr(handoff.subprocess, 'check_output', untracked)
    result = handoff.archive_legacy_example_layout(root)
    assert len(result) == len(calls) == 1
    assert (Path(result[0]['archive'])/'orbitals/old.wfu').read_bytes() == b'retained native data'
    assert not old.parent.exists()
    assert (unrelated/'input.json').is_file()
    assert json.loads((root/'repository_layout_archive.json').read_text()) == result
    assert handoff.archive_legacy_example_layout(root) == []


def test_layout_archive_refuses_all_moves_if_any_candidate_is_tracked(handoff, tmp_path, monkeypatch):
    repo = tmp_path/'repo'
    first = repo/'examples/anl/caspt2_production'; first.mkdir(parents=True)
    second = repo/'examples/anl/orbital_reviews'; second.mkdir()
    root = tmp_path/'new_run'; root.mkdir()
    monkeypatch.setattr(handoff, 'REPO', repo)
    monkeypatch.setattr(handoff.subprocess, 'check_output',
        lambda command, **kwargs: b'source.py\0' if command[-1].endswith('orbital_reviews') else b'')
    with pytest.raises(ValueError, match='tracked example files'):
        handoff.archive_legacy_example_layout(root)
    assert first.is_dir() and second.is_dir()
    assert not (root/'repository_layout_archive.json').exists()


def test_two_pes_salvage_resumes_interrupted_preparation_without_losing_inputs(handoff, tmp_path, monkeypatch):
    import errno
    from kinbot import run_salvage
    from kinbot.anl import salvage, dispatch
    from tests.test_anl_salvage import small_spec
    from kinbot.execution import atomic_json, file_sha256
    from kinbot.parameters import Parameters
    from kinbot.qc import QuantumChemistry
    old = handoff.prepare(tmp_path/'old')
    new = tmp_path/'new'
    for case in handoff.CASES:
        worker = old/case/'123'
        graph = worker/'anl_automatic/states/example/composite'
        graph.parent.mkdir(parents=True)
        dispatch.prepare(small_spec(), graph)
        # Existing L1/L2 routes were generated at the original tighter SCF.
        par = Parameters(str(old/case/'input.json')).par
        par['anl_scf_convergence'] = 1e-10
        with monkeypatch.context() as context:
            context.chdir(worker)
            QuantumChemistry(par)._record('accepted_high', 'l2', '201')
        atomic_json(worker/'anl_automatic/network_request.json', {'resources': {}})
        (worker/'kinbot.log').write_text('old failure evidence')
        ref = worker/'anl_automatic/reference_l2/failed_reference'
        ref.mkdir(parents=True)
        (ref/'driver.stderr').write_text('failed SCF')
    monkeypatch.setattr(run_salvage, '_inventory', lambda root: {'jobs': [], 'processes': []})
    monkeypatch.setattr(run_salvage, '_exchange', lambda *a: (_ for _ in ()).throw(OSError(errno.EOPNOTSUPP, 'test shared filesystem')))
    monkeypatch.setattr(run_salvage, '_pin_code', lambda root: root/'.salvage_code')
    repair = salvage.repair_graph
    calls = []
    def interrupted(graph, archive, **options):
        calls.append(graph)
        if len(calls) == 2:
            raise RuntimeError('interrupted preparation')
        return repair(graph, archive, **options)
    monkeypatch.setattr(salvage, 'repair_graph', interrupted)
    with pytest.raises(RuntimeError, match='interrupted preparation'):
        run_salvage.prepare(old, new)
    monkeypatch.setattr(salvage, 'repair_graph', repair)
    report = run_salvage.prepare(old, new)
    assert report['status'] == 'prepared'
    assert run_salvage.prepare(old, new) == report
    manifest = handoff.verify(old)
    assert len(manifest['inputs']) == 2
    for filename in manifest['inputs']:
        values = json.loads(Path(filename).read_text())
        assert values['vrc_tst_pivot_backend'] == 'molpro'
        assert values['anl_require_l3_geometry'] == 1
        assert values['vrc_tst_surface_model']['channels']
        assert values['queue_job_limit'] == 10
    assert len(list(Path(report['archive']).glob('*/123/kinbot.log'))) == 2
    assert len(list(Path(report['archive']).glob('*/123/anl_automatic/reference_l2/failed_reference/driver.stderr'))) == 2
    for case in handoff.CASES:
        payload = json.loads((old/case/'123/.kinbot_theory_jobs.json').read_text())
        assert payload['profiles']['l2']['calculator_kwargs']['scf'] == 'xqc,conver=8'
        assert payload['jobs']['accepted_high']['job_id'] is None
        assert (Path(report['archive'])/case/'123/.kinbot_theory_jobs.json').exists()


def test_salvage_resumes_after_cancel_and_rehome_with_live_qz(handoff, tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot import run_cleanup, run_salvage
    from kinbot.anl import salvage, dispatch
    from kinbot.execution import atomic_json, file_sha256
    from tests.test_anl_salvage import small_spec
    from tests.test_anl_dispatch import _mock_execution

    old = handoff.prepare(tmp_path/'old')
    new = tmp_path/'new'
    for case in handoff.CASES:
        base = old/case/'123/anl_automatic/states'
        for name in ('a_accepted', 'b_fallback') if case == 'ch3s' else ('a_accepted',):
            graph = base/name/'composite'
            graph.parent.mkdir(parents=True)
            dispatch.prepare(small_spec(), graph)
            _mock_execution(graph, 'accepted_geometry', geometry='final.xyz')
            _mock_execution(graph, 'l2_property')
            dispatch.advance(graph)
    invalid = old/'ch3s/123/anl_automatic/states/b_fallback/composite'
    state = dispatch._load(invalid)[2]
    state['tasks']['l3_geometry']['status'] = 'unavailable'
    atomic_json(invalid/'state.json', state)
    state = dispatch.advance(invalid)
    state['tasks']['energy'].update(status='submitted', job_id='102')
    atomic_json(invalid/'state.json', state)
    live = old/'ch3s/123/anl_automatic/states/a_accepted/composite'
    state = dispatch._load(live)[2]
    state['tasks']['l3_geometry'].update(status='submitted', job_id='101')
    atomic_json(live/'state.json', state)
    live_hashes = {str(p.relative_to(live)): file_sha256(p)
                   for p in (live/'tasks/l3_geometry').iterdir() if p.is_file()}
    jobs = [{'job_id': '101', 'directory': str(live/'tasks/l3_geometry')},
            {'job_id': '102', 'directory': str(invalid/'tasks/energy')}]
    cancelled = []
    def cancel(root, selected):
        cancelled.extend(j['job_id'] for j in selected)
        jobs[:] = [j for j in jobs if j not in selected]
        # scancel need not leave an execution.json result behind.
    monkeypatch.setattr(run_cleanup, '_cancel', cancel)
    monkeypatch.setattr(run_cleanup, '_jobs', lambda root: list(jobs))
    monkeypatch.setattr(run_cleanup, '_processes', lambda root: [])
    monkeypatch.setattr(dispatch, '_job_active', lambda job: any(j['job_id'] == job for j in jobs))
    def exchange(a, b):
        temporary = tmp_path/'exchange'
        a.rename(temporary)
        b.rename(a)
        temporary.rename(b)
    monkeypatch.setattr(run_salvage, '_exchange', exchange)
    monkeypatch.setattr(run_salvage, '_pin_code', lambda root: root/'.salvage_code')
    repair = salvage.repair_graph
    def old_failure(graph, archive, **kwargs):
        if graph.parent.name == 'b_fallback':
            state = dispatch.advance(graph, submit=False)
            assert state['tasks']['energy']['status'] == 'submitted'
            assert 'missing_since' in state['tasks']['energy']
            raise RuntimeError('old stale submitted guard')
        return repair(graph, archive, **kwargs)
    monkeypatch.setattr(salvage, 'repair_graph', old_failure)
    with pytest.raises(RuntimeError, match='old stale submitted guard'):
        run_salvage.prepare(old, new)
    assert cancelled == ['102']
    assert old.is_symlink() and old.resolve() == new
    assert not list(new.glob('*/launch.json'))
    first_migration = (live/'salvage.json').read_bytes()

    # Resume the exact partial layout left on Blodgett, during publication
    # grace, with the original old-path live job still in the queue.
    monkeypatch.setattr(salvage, 'repair_graph', repair)
    report = run_salvage.prepare(old, new)
    assert report['status'] == 'prepared'
    assert cancelled == ['102']  # no second cancellation / no valid job stopped
    assert (live/'salvage.json').read_bytes() == first_migration
    assert live_hashes == {str(p.relative_to(live)): file_sha256(p)
                           for p in (live/'tasks/l3_geometry').iterdir() if p.is_file()}
    records = {Path(r['graph']).parent.name: r for r in report['graphs']
               if '/ch3s/' in r['graph']}
    assert records['a_accepted']['adopted_jobs'] == {'l3_geometry': '101'}
    assert records['b_fallback']['retired_jobs'] == {'energy': '102'}
    assert (Path(records['b_fallback']['archive'])/'energy/task.json').exists()
    assert not (invalid/'tasks/energy').exists()
    assert dispatch._load(invalid)[2]['tasks']['l3_geometry']['status'] == 'staged'
    assert len(handoff.verify(new)['inputs']) == 2

    # Both detached drivers start exactly once after successful preparation.
    monkeypatch.setattr('kinbot.anl.automatic.validate_automatic_environment', lambda par: None)
    launches = []
    def start(command, **kwargs):
        launches.append((command, kwargs))
        return SimpleNamespace(pid=700+len(launches))
    monkeypatch.setattr(run_salvage.subprocess, 'Popen', start)
    launched = run_salvage.launch(new)
    assert launched['status'] == 'launched'
    assert run_salvage.launch(new) == launched
    assert len(launches) == 2
    assert {kwargs['cwd'].name for _, kwargs in launches} == {'ch3s', 'ch3ss'}
    assert all(kwargs['start_new_session'] for _, kwargs in launches)
    assert len(list(new.glob('*/launch.json'))) == 2


def test_resume_failed_theory_startup_preserves_qc_and_preflights_real_router(handoff, tmp_path, monkeypatch):
    from copy import deepcopy
    from types import SimpleNamespace
    from ase import Atoms
    from kinbot import run_cleanup, run_salvage
    from kinbot.anl import dispatch
    from kinbot.parameters import Parameters
    from kinbot.qc import QuantumChemistry
    from kinbot.execution import atomic_json, file_sha256
    from tests.test_anl_salvage import small_spec
    from tests.test_anl_dispatch import _mock_execution
    root = handoff.prepare(tmp_path/'run')
    archive = tmp_path/'archive'; archive.mkdir()
    record = json.loads((root/'prepared.json').read_text())
    old_ledgers, databases = {}, {}
    error = 'ValueError: Theory profiles changed since jobs were submitted; archive or remove .kinbot_theory_jobs.json before a new run.'
    for filename in record['inputs']:
        path = Path(filename)
        worker = path.parent/'123'; worker.mkdir()
        par = Parameters(filename).par
        before = deepcopy(par)
        before['anl_scf_convergence'] = 1e-10
        with monkeypatch.context() as context:
            context.chdir(worker)
            qc = QuantumChemistry(before)
            qc._record('accepted_well', 'l1', '201')
            qc.l1.db.write(Atoms('H'), name='accepted_well',
                           data={'status': 'normal', 'energy': -10., 'zpe': 0., 'frequencies': []})
        ledger = worker/'.kinbot_theory_jobs.json'
        old_ledgers[ledger] = ledger.read_bytes()
        databases[worker/'kinbot.db'] = file_sha256(worker/'kinbot.db')
        (worker/'kinbot.log').write_text(error)
        atomic_json(worker/'kinbot_status.json', {'status': 'failed', 'error': error})
        atomic_json(path.parent/'pes_status.json', {'status': 'failed', 'error': error})
        atomic_json(path.parent/'launch.json', {'pid': 1234})
    graph = root/'ch3s/123/anl_automatic/states/valid/composite'
    graph.parent.mkdir(parents=True)
    dispatch.prepare(small_spec(), graph)
    _mock_execution(graph, 'accepted_geometry', geometry='final.xyz')
    state = dispatch.advance(graph)
    state['tasks']['l3_geometry'].update(status='submitted', job_id='101')
    atomic_json(graph/'state.json', state)
    graph_hashes = {str(p.relative_to(graph)): file_sha256(p)
                    for p in graph.rglob('*') if p.is_file()}
    atomic_json(root/'salvage.json', {'status': 'launched', 'archive': str(archive),
                                    'inputs': record['inputs'], 'code': str(root/'.salvage_code')})
    pin = root/'.salvage_code'; pin.mkdir()
    (pin/'old_code.py').write_text('# original pinned source\n')
    atomic_json(pin/'manifest.json', {'old_code.py': file_sha256(pin/'old_code.py')})
    monkeypatch.setattr(run_cleanup, '_jobs', lambda root: [
        {'job_id': '101', 'directory': str(graph/'tasks/l3_geometry')}])
    monkeypatch.setattr(run_cleanup, '_processes', lambda root: [])
    with pytest.raises(ValueError, match='Theory profiles changed'):
        run_salvage._preflight_worker_routes(root)
    pin_code = run_salvage._pin_code
    def interrupted_pin(root, **kwargs):
        raise RuntimeError('interrupted before new code snapshot')
    monkeypatch.setattr(run_salvage, '_pin_code', interrupted_pin)
    with pytest.raises(RuntimeError, match='interrupted before new code snapshot'):
        run_salvage.resume_failed_controllers(root)
    monkeypatch.setattr(run_salvage, '_pin_code', pin_code)
    result = run_salvage.resume_failed_controllers(root)
    assert result['status'] == 'prepared'
    assert result['code'] != str(pin) and (pin/'old_code.py').exists()
    assert Path(result['code']).is_dir()
    assert run_salvage.resume_failed_controllers(root) == result
    assert not list(root.glob('*/launch.json'))
    saved = Path(result['controller_restart']['archive'])
    for ledger, content in old_ledgers.items():
        assert (saved/ledger.relative_to(root)).read_bytes() == content
        revised = json.loads(ledger.read_text())
        assert revised['profiles']['l1']['calculator_kwargs']['scf'] == 'xqc,conver=8'
        assert revised['jobs']['accepted_well']['job_id'] is None
        assert revised['jobs']['accepted_well']['source_profile_sha256'] == json.loads(content)['profile_sha256']
    assert databases == {p: file_sha256(p) for p in databases}
    assert graph_hashes == {str(p.relative_to(graph)): file_sha256(p)
                            for p in graph.rglob('*') if p.is_file()}
    run_salvage._preflight_worker_routes(root)
    monkeypatch.setattr('kinbot.anl.automatic.validate_automatic_environment', lambda par: None)
    launches = []
    def start(command, **kwargs):
        launches.append((command, kwargs))
        return SimpleNamespace(pid=800+len(launches))
    monkeypatch.setattr(run_salvage.subprocess, 'Popen', start)
    assert run_salvage.launch(root)['status'] == 'launched'
    assert run_salvage.launch(root)['status'] == 'launched'
    assert len(launches) == 2
    assert all(kwargs['env']['PYTHONPATH'].startswith(result['code']) for _, kwargs in launches)
    with pytest.raises(RuntimeError, match='already launched'):
        run_salvage.resume_failed_controllers(root)


def test_failed_controller_recovery_refuses_live_controller_or_other_failure(handoff, tmp_path, monkeypatch):
    from kinbot import run_cleanup, run_salvage
    from kinbot.execution import atomic_json
    root = handoff.prepare(tmp_path/'run')
    atomic_json(root/'salvage.json', {'status': 'launched'})
    monkeypatch.setattr(run_cleanup, '_jobs', lambda root: [])
    monkeypatch.setattr(run_cleanup, '_processes', lambda root: [{'pid': 123}])
    with pytest.raises(RuntimeError, match='still alive'):
        run_salvage.resume_failed_controllers(root)
    monkeypatch.setattr(run_cleanup, '_processes', lambda root: [])
    atomic_json(root/'ch3s/pes_status.json', {'status': 'failed', 'error': 'different failure'})
    with pytest.raises(RuntimeError, match='expected the reviewed'):
        run_salvage.resume_failed_controllers(root)
