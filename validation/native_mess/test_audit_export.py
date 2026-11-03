import importlib.util
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location('mess_export_audit', Path(__file__).with_name('audit_export.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

VALID = '''BEGIN 1 300 1 0.1
STATE 0 w_1 0 1 1 1
STATE 1 w_1 1 0 2 4
COLLISION_DIAGONAL 0 4
COLLISION_DIAGONAL 1 1
EDGE 0 1 4 collision
EDGE 1 0 1 collision
SINK 0 pr_1 3 barrier
END
'''


def run(tmp_path, text):
    path = tmp_path/'generator.txt'
    path.write_text(text)
    return module.audit(path)


def test_balanced_native_rates_with_physical_sink(tmp_path):
    result = run(tmp_path, VALID)
    assert result['states'] == 2 and result['sink_records'] == 1
    assert result['max_detailed_balance_relative_error'] == 0
    assert result['max_collision_diagonal_relative_error'] == 0


@pytest.mark.parametrize('old,new', [
    ('SINK 0 pr_1 3', 'SINK 0 pr_1 -3'),
    ('EDGE 0 1 4', 'EDGE 0 1 nan'),
    ('0 1 1 1', '0 1 inf 1'),
    ('END\n', ''),
    ('END\n', 'END\nBEGIN 1 300 1 .1\n'),
    ('STATE 1', 'STATE 0'),
    ('EDGE 0 1', 'EDGE 0 2'),
])
def test_invalid_native_evidence_is_not_accepted(tmp_path, old, new):
    with pytest.raises(ValueError):
        run(tmp_path, VALID.replace(old, new))


def test_reports_balance_error_without_repair(tmp_path):
    result = run(tmp_path, VALID.replace('EDGE 0 1 4', 'EDGE 0 1 8'))
    assert result['max_collision_diagonal_relative_error'] == .5
    assert result['max_detailed_balance_relative_error'] == pytest.approx(.5)


def test_large_equilibrium_weights_do_not_overflow(tmp_path):
    result = run(tmp_path, VALID.replace('0 1 1 1', '0 1 1e200 1').replace('1 0 2 4', '1 0 2e200 4'))
    assert result['max_detailed_balance_relative_error'] < 1e-12
