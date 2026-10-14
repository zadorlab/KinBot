import numpy as np
import pytest

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
