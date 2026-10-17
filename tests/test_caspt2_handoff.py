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
                assert 'thresh,energy=1.d-7,coeff=1.d-10}' in deck
            template = task.get('input_template', '')
            if task.get('backend') == 'gaussian':
                assert 'SCF=(XQC,Conver=8)' in template
            elif task.get('backend') == 'molpro':
                assert 'gthresh,energy=1.d-8,' in template
                if ';thresh,energy=' in template:
                    assert 'thresh,energy=1.d-7,coeff=1.d-10}' in template
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
