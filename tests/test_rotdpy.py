"""General rotdPy input generation, execution, and restart contracts."""

import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pytest

from kinbot.parameters import Parameters
from kinbot.pes import create_rotdpy_inputs
from kinbot.rotdpy import ensure_available, run


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
            assert 'temperature = generate_grid(*[100.0, 50.0, 1.0, 2])' in generated
            assert "'pot_smp_min': 4" in generated
            scheduler = Path('rotdPy/qu.tpl').read_text()
            assert '#SBATCH --partition=short-cpu' in scheduler
            assert '#SBATCH --ntasks={procs}' in scheduler
            assert '#SBATCH --time=01:30:00' in scheduler
            assert str(Path(sys.executable)) in scheduler
        finally:
            os.chdir(previous)


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


def test_executor_fails_early_when_rotdpy_is_not_installed():
    with TemporaryDirectory() as temporary:
        input_file = Path(temporary) / 'channel.py'
        input_file.write_text('raise AssertionError("must not execute")\n')
        with patch('kinbot.rotdpy.importlib.util.find_spec',
                   return_value=None):
            with pytest.raises(RuntimeError, match='same environment'):
                run(input_file, python=sys.executable)


def test_rotdpy_revision_must_match_pinned_submodule():
    with patch('kinbot.rotdpy.dependency_provenance', return_value={
            'version': '0.1.0', 'location': '/different',
            'git_revision': 'wrong'}):
        with pytest.raises(RuntimeError, match='pinned'):
            ensure_available()
