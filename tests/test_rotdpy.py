"""General rotdPy input generation, execution, and restart contracts."""

import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pytest

from kinbot.parameters import Parameters
from kinbot.pes import _rotdpy_correction_block, create_rotdpy_inputs
from kinbot.rotdpy import (ensure_available, main, number_of_states_file,
                           read_result, run, sampling_progress)
from kinbot.molpro import _vrc_multireference_method


def test_generated_input_uses_configured_sampling_and_portable_scratch():
    with TemporaryDirectory() as temporary:
        previous = Path.cwd()
        os.chdir(temporary)
        try:
            Path('input.json').write_text(json.dumps({
                'barrier_threshold': 50.,
                'queuing': 'slurm',
                'queue_name': 'short-cpu',
                'vrc_tst_noscan': {'parent': ['parent_hom_sci_1_2']},
                'rotdpy_dist': [3.0],
                'rotdpy_processors': 2,
                'rotdpy_max_jobs': 3,
                'rotdpy_walltime': '01:30:00',
                'rotdpy_temperature_grid': [100., 50., 1., 2],
                'rotdpy_energy_grid': [0., 5., 1., 3],
                'rotdpy_angular_grid': [0., 1., 1., 4],
                'rotdpy_flux_parameters': {
                    'pot_smp_max': 8, 'pot_smp_min': 4,
                    'tot_smp_max': 8, 'tot_smp_min': 4,
                    'flux_rel_err': 5, 'smp_len': 1,
                },
            }))
            par = Parameters('input.json', show_warnings=False).par
            Path('vrctst').mkdir()
            reaction = 'parent_hom_sci_1_2'
            Path(f'vrctst/corr_{reaction}.json').write_text(json.dumps({
                'dist': [30.], 'e_samp': [0.], 'e_high': [0.],
                'scan_ref': [[0, 0]], 'ra': [[0], [0]],
                'e_inf_samp': -1., 'e_inf_high': -1.,
                'frags_atom': [['H'], ['H']],
                'frags_geom': [[[0., 0., 0.]], [[1., 0., 0.]]],
                'frags_mult': [2, 2], 'unique': [[[0]], [[0]]],
            }))
            created = create_rotdpy_inputs(
                par, [['parent', reaction, ['h', 'h'], 0.]], [],
                correction_root='.')
            assert created == [f'rotdPy/{reaction}.py']
            generated = Path(created[0]).read_text()
            compile(generated, created[0], 'exec')
            assert "'method': 'caspt2(2,2)'" in generated
            assert "'basis': 'vdz'" in generated
            assert "'processors': 2" in generated
            assert "'max_jobs': 3" in generated
            assert "os.environ.get('SCRATCH')" in generated
            assert 'corrections = None' in generated
            assert "'KinBot 1D scan'" not in generated
            assert 'temperature = generate_grid(*[100.0, 50.0, 1.0, 2])' in generated
            assert "'pot_smp_min': 4" in generated
            scheduler = Path('rotdPy/qu.tpl').read_text()
            assert '#SBATCH --partition=short-cpu' in scheduler
            assert '#SBATCH --ntasks={procs}' in scheduler
            assert '#SBATCH --time=01:30:00' in scheduler
            assert 'export I_MPI_FABRICS="${{I_MPI_FABRICS:-shm}}"' in scheduler
            sample_job = scheduler.format(
                surf_id=0, face_id=1, samp_id=2, procs=2, mem=4096)
            assert 'export I_MPI_FABRICS="${I_MPI_FABRICS:-shm}"' in sample_job
            assert '#SBATCH --ntasks=2' in sample_job
            assert str(Path(sys.executable)) in scheduler
        finally:
            os.chdir(previous)


def test_multipoint_vrc_scan_retains_cubic_correction_contract():
    correction = {
        'dist': [3., 4., 5., 30.],
        'e_samp': [4., 2., 1., 0.],
        'e_high': [5., 2.5, 1.2, 0.],
        'scan_ref': [[0, 1]],
    }
    rendered = _rotdpy_correction_block(correction, noscan=False)
    assert "'KinBot 1D scan'" in rendered
    assert "'r_sample' : [3.0, 4.0, 5.0, 30.0]" in rendered

    with pytest.raises(ValueError, match='three scan points plus'):
        _rotdpy_correction_block({
            **correction,
            'dist': [3., 4., 30.],
            'e_samp': [2., 1., 0.],
            'e_high': [2.5, 1.2, 0.],
        }, noscan=False)


def test_vrc_mrci_q_uses_davidson_corrected_energy():
    method, energy = _vrc_multireference_method('MRCI+Q(2,2)', 18)
    assert 'occ,10' in method
    assert 'closed,8' in method
    assert 'wf,18,1,0' in method
    assert '{mrci}' in method
    assert energy == 'energd(1)'


def test_production_result_requires_converged_numeric_multipoint_sampling(
        tmp_path):
    input_file = tmp_path / 'channel.py'
    input_file.write_text('# production input\n')
    result_root = tmp_path / 'kb_channel'
    outputs = [result_root / 'Ne_0.out'] + [
        result_root / 'output' / f'surface_{index}.dat'
        for index in range(3)]
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text('0.0 1.0\n')
    relatives = [str(path.relative_to(tmp_path)) for path in outputs]
    hashes = {relative: __import__('hashlib').sha256(
        path.read_bytes()).hexdigest()
              for relative, path in zip(relatives, outputs)}
    validation = {
        'name': 'production',
        'sampling_level': {'method': 'caspt2(2,2)', 'basis': 'vdz'},
        'trusted_correction_level': {
            'method': 'mrci+q(2,2)', 'basis': 'vtz'},
        'correction': {'kind': 'one_dimensional', 'point_count': 4,
                       'source_sha256': 'a' * 64},
        'dividing_surfaces': {
            'requested_distances_angstrom': [3., 4., 5.],
            'generated_count': 3},
    }
    statistics = {
        str(index): {'converged': True, 'accepted_samples': 100}
        for index in range(3)}
    manifest = {
        'schema': 3, 'status': 'complete', 'reaction': 'channel',
        'surface_count': 3,
        'input_sha256': __import__('hashlib').sha256(
            input_file.read_bytes()).hexdigest(),
        'validation': validation, 'surface_statistics': statistics,
        'result_files': relatives, 'result_sha256': hashes,
    }
    (tmp_path / 'channel.rotdpy.json').write_text(json.dumps(manifest))
    assert read_result(input_file, required_profile='production')[
        'surface_count'] == 3

    manifest['surface_statistics']['1']['converged'] = False
    (tmp_path / 'channel.rotdpy.json').write_text(json.dumps(manifest))
    with pytest.raises(RuntimeError, match='not converged'):
        read_result(input_file, required_profile='production')


def test_production_selection_rejects_interface_smoke_manifest(tmp_path):
    input_file = tmp_path / 'channel.py'
    input_file.write_text('# smoke input\n')
    surface = tmp_path / 'kb_channel' / 'output' / 'surface_0.dat'
    number = tmp_path / 'kb_channel' / 'Ne_0.out'
    surface.parent.mkdir(parents=True)
    surface.write_text('0.0 1.0\n')
    number.write_text('0.0 1.0\n')
    relatives = [str(number.relative_to(tmp_path)),
                 str(surface.relative_to(tmp_path))]
    hashes = {relative: __import__('hashlib').sha256(
        (tmp_path / relative).read_bytes()).hexdigest()
              for relative in relatives}
    (tmp_path / 'channel.rotdpy.json').write_text(json.dumps({
        'schema': 2, 'status': 'complete', 'reaction': 'channel',
        'surface_count': 1, 'result_files': relatives,
        'result_sha256': hashes,
    }))
    with pytest.raises(RuntimeError, match='not production validated'):
        number_of_states_file(input_file, required_profile='production')


def test_executor_records_success_and_reuses_matching_result():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        input_file = root / 'channel.py'
        (root / 'qu.tpl').write_text('# test scheduler template\n')
        surface = root / 'kb_channel' / 'output' / 'surface_0.dat'
        number = root / 'kb_channel' / 'Ne_0.out'
        surface.parent.mkdir(parents=True)
        surface.write_text('surface flux\n')
        number.write_text('0.0 1.0\n')
        surface_hash = __import__('hashlib').sha256(surface.read_bytes()).hexdigest()
        number_hash = __import__('hashlib').sha256(number.read_bytes()).hexdigest()
        input_file.write_text(
            "import json\n"
            "from pathlib import Path\n"
            "Path('channel.rotdpy.json').write_text(json.dumps({"
            "'schema': 2, 'status': 'complete', 'reaction': 'channel', "
            "'surface_count': 2, "
            "'result_files': ['kb_channel/Ne_0.out', "
            "'kb_channel/output/surface_0.dat'], "
            f"'result_sha256': {{'kb_channel/Ne_0.out': '{number_hash}', "
            f"'kb_channel/output/surface_0.dat': '{surface_hash}'}}"
            "}))\n"
            "print('sampling complete')\n")
        with patch('kinbot.rotdpy.ensure_available', return_value={
                'version': 'test', 'location': '/test',
                'git_revision': 'abc'}):
            record = run(input_file, python=sys.executable)
            reused = run(input_file, python=sys.executable)
        assert record['status'] == 'complete'
        assert record['surface_count'] == 2
        assert reused == record
        assert (root / 'channel.rotdpy.stdout').read_text() == (
            'sampling complete\n')


def test_number_of_states_selection_is_hash_verified():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        input_file = root / 'channel.py'
        input_file.write_text('# generated input\n')
        result_root = root / 'kb_channel'
        surface = result_root / 'output' / 'surface_0.dat'
        first = result_root / 'Ne_0.out'
        final = result_root / 'Ne_1.out'
        surface.parent.mkdir(parents=True)
        surface.write_text('surface flux\n')
        first.write_text('0.0 1.0\n')
        final.write_text('0.0 2.0\n')
        files = [first, final, surface]
        relatives = [str(path.relative_to(root)) for path in files]
        hashes = {
            relative: __import__('hashlib').sha256(path.read_bytes()).hexdigest()
            for relative, path in zip(relatives, files)
        }
        (root / 'channel.rotdpy.json').write_text(json.dumps({
            'schema': 2, 'status': 'complete', 'reaction': 'channel',
            'surface_count': 1, 'result_files': relatives,
            'result_sha256': hashes,
        }))

        selected, provenance = number_of_states_file(input_file)
        assert selected == final.resolve()
        assert provenance['energy_index'] == 1
        assert provenance['selection'] == 'highest_available_correction'
        explicit, provenance = number_of_states_file(input_file, 0)
        assert explicit == first.resolve()
        assert provenance['selection'] == 'explicit'
        with pytest.raises(RuntimeError, match='no Ne_2.out'):
            number_of_states_file(input_file, 2)
        final.write_text('changed\n')
        with pytest.raises(RuntimeError, match='hash changed'):
            number_of_states_file(input_file)


def test_executor_fails_early_when_rotdpy_is_not_installed():
    with TemporaryDirectory() as temporary:
        input_file = Path(temporary) / 'channel.py'
        input_file.write_text('raise AssertionError("must not execute")\n')
        with patch('kinbot.rotdpy.importlib.util.find_spec',
                   return_value=None):
            with pytest.raises(RuntimeError, match='same environment'):
                run(input_file, python=sys.executable)


def test_executor_reports_captured_rotdpy_stderr():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        input_file = root / 'channel.py'
        input_file.write_text(
            'import sys\n'
            'sys.stderr.write("scheduler template field failure\\n")\n'
            'raise SystemExit(3)\n')
        (root / 'qu.tpl').write_text('# scheduler fixture\n')
        with patch('kinbot.rotdpy.ensure_available', return_value={
                'version': 'test', 'location': '/test',
                'git_revision': 'abc'}):
            with pytest.raises(RuntimeError,
                               match='status 3.*scheduler template field'):
                run(input_file, python=sys.executable)
        record = json.loads((root / 'channel.execution.json').read_text())
        assert record['returncode'] == 3
        assert 'scheduler template field failure' in record['error']
        assert (root / record['stderr']).read_text() == \
            'scheduler template field failure\n'


def test_rotdpy_revision_must_match_pinned_submodule():
    with patch('kinbot.rotdpy.dependency_provenance', return_value={
            'version': '0.1.0', 'location': '/different',
            'git_revision': 'wrong'}):
        with pytest.raises(RuntimeError, match='pinned'):
            ensure_available()


def test_command_line_check_reports_completed_manifest(tmp_path, capsys):
    input_file = tmp_path / 'channel.py'
    input_file.write_text('# generated input\n')
    surface = tmp_path / 'kb_channel' / 'output' / 'surface_0.dat'
    number = tmp_path / 'kb_channel' / 'Ne_0.out'
    surface.parent.mkdir(parents=True)
    surface.write_text('surface flux\n')
    number.write_text('0.0 1.0\n')
    files = [number, surface]
    relative = [str(path.relative_to(tmp_path)) for path in files]
    hashes = {name: __import__('hashlib').sha256(path.read_bytes()).hexdigest()
              for name, path in zip(relative, files)}
    (tmp_path / 'channel.rotdpy.json').write_text(json.dumps({
        'schema': 2, 'status': 'complete', 'reaction': 'channel',
        'surface_count': 1, 'result_files': relative,
        'result_sha256': hashes,
    }))
    assert main(['check', str(input_file)]) == 0
    assert json.loads(capsys.readouterr().out)['status'] == 'complete'


def test_live_progress_counts_converged_surfaces_and_active_samples(tmp_path):
    input_file = tmp_path / 'channel.py'
    input_file.write_text('# generated input\n')
    root = tmp_path / 'kb_channel'
    for index in range(3):
        (root / f'Surface_{index}' / 'jobs').mkdir(parents=True)
    output = root / 'output'
    output.mkdir()
    (output / 'surface_0.dat').write_text('complete\n')
    (root / 'Surface_1/jobs/surf1_face0_samp843.pkl').write_bytes(b'x')

    with patch('kinbot.rotdpy.time.time',
               return_value=input_file.stat().st_mtime + 120.):
        progress = sampling_progress(input_file)
    assert progress['status'] == 'running'
    assert progress['surface_count'] == 3
    assert progress['converged_surfaces'] == 1
    assert progress['remaining_surfaces'] == 2
    assert progress['active_samples'] == [
        {'surface': 1, 'face': 0, 'sample': 843}]
    assert progress['elapsed_seconds'] == 120
    assert progress['linear_eta_seconds'] == 240
