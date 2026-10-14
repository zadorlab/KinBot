import json
from pathlib import Path

import numpy as np
import pytest

from kinbot.execution import atomic_json, file_sha256, scheduler_unavailable
from kinbot.reader_qchem import read_hess
from kinbot.utils import queue_command
from kinbot.zmatrix import get_minimum_pathway


def test_invalid_queue_name_reports_the_supplied_value():
    with pytest.raises(ValueError, match='Unexpected value for queueing: invalid'):
        queue_command('invalid')


def test_disconnected_rotors_use_the_imported_motif_search():
    bond = np.zeros((4, 4), dtype=int)
    for first, second in ((0, 1), (1, 2), (2, 3)):
        bond[first, second] = bond[second, first] = 1
    distance, path = get_minimum_pathway(
        bond, [0], [3], 4, ['C', 'C', 'C', 'C'])
    assert distance == 2
    assert path == [0, 1, 2, 3]
    reverse_distance, reverse_path = get_minimum_pathway(
        bond, [3], [0], 4, ['C', 'C', 'C', 'C'])
    assert reverse_distance == distance
    assert reverse_path == [3, 2, 1, 0]


def test_qchem_hessian_reader_has_its_runtime_dependencies(tmp_path):
    job = tmp_path / 'qchem'
    job.with_suffix('.fchk').write_text(
        'Cartesian Force Constants\n'
        '1.0 0.0 2.0 0.0 0.0 3.0\n')
    assert read_hess(job, 1) == pytest.approx(np.diag([1., 2., 3.]))


def test_hessian_reader_invokes_formchk_without_shell_parsing(
        tmp_path, monkeypatch):
    job = tmp_path / 'job with spaces'
    job.with_suffix('.chk').write_text('checkpoint')

    def formchk(command, **kwargs):
        assert command == [
            'formchk', str(job.with_suffix('.chk')),
            str(job.with_suffix('.fchk'))]
        assert kwargs['check'] is True
        job.with_suffix('.fchk').write_text(
            'Cartesian Force Constants\n'
            '1.0 0.0 2.0 0.0 0.0 3.0\n')

    monkeypatch.setattr('kinbot.reader_gauss.subprocess.run', formchk)
    assert read_hess(job, 1) == pytest.approx(np.diag([1., 2., 3.]))


def test_shared_execution_helpers(tmp_path):
    destination = tmp_path / 'record.json'
    atomic_json(destination, {'value': 3})
    assert json.loads(destination.read_text()) == {'value': 3}
    assert not Path(f'{destination}.tmp').exists()
    assert file_sha256(destination) == (
        'a3c2f33d714cb8590e407c83bd199c5351b3e26897acfaaf4264353ede52dfaf')
    assert scheduler_unavailable('slurm_load_jobs error: controller down')
    assert not scheduler_unavailable('job exited with status 1')
