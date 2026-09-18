"""Native result labels and units are part of the scientific acceptance gate."""

from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import pytest

from kinbot.anl.dispatch import _run_external
from kinbot.anl.results import parse_cfour_dboc


# The values and surrounding labels are from the completed Blodgett CH4 job.
_CFOUR_CH4 = """Total SCF energy including DBOC  -40.210746196622267
Summary of diagonal Born-Oppenheimer correction at Hartree-Fock level
The total diagonal Born-Oppenheimer correction (DBOC) is: 0.0025887093 a.u.
The total diagonal Born-Oppenheimer correction (DBOC) is: 568.156016 cm-1
The total diagonal Born-Oppenheimer correction (DBOC) is: 6.797 kJ/mole
Summary of diagonal Born-Oppenheimer correction at MP1 level
The total diagonal Born-Oppenheimer correction (DBOC) is: 0.0026718675 a.u.
The total diagonal Born-Oppenheimer correction (DBOC) is: 586.407129 cm-1
The total diagonal Born-Oppenheimer correction (DBOC) is: 7.015 kJ/mole
The final electronic energy is -40.210746196622267 a.u.
This computation required 31.93 seconds (walltime).
"""


def test_cfour_selects_hf_dboc_without_confusing_mp1_or_final_energy():
    result = parse_cfour_dboc(_CFOUR_CH4)
    assert result['selected_level'] == 'HF'
    assert result['selected']['hartree'] == pytest.approx(0.0025887093)
    assert result['selected']['cm_inverse'] == pytest.approx(568.156016)
    assert result['reported_levels']['MP1']['hartree'] == pytest.approx(0.0026718675)
    assert result['final_energy_including_dboc_hartree'] == pytest.approx(
        -40.210746196622267)
    assert parse_cfour_dboc(_CFOUR_CH4, level='MP1')['selected']['hartree'] == \
        pytest.approx(0.0026718675)


def test_cfour_rejects_incomplete_and_inconsistent_dboc_output():
    with pytest.raises(ValueError, match='completion'):
        parse_cfour_dboc(_CFOUR_CH4.replace('This computation required',
                                           'Computation did not complete'))
    with pytest.raises(ValueError, match='a.u./cm-1'):
        parse_cfour_dboc(_CFOUR_CH4.replace('568.156016', '500.000000'))
    with pytest.raises(ValueError, match='repeated the HF'):
        parse_cfour_dboc(_CFOUR_CH4.replace(
            'Summary of diagonal Born-Oppenheimer correction at MP1 level',
            'Summary of diagonal Born-Oppenheimer correction at Hartree-Fock level'))


def test_external_cfour_task_stores_selected_native_result():
    with TemporaryDirectory() as temporary:
        directory = Path(temporary)
        (directory / 'ZMAT').write_text('*CFOUR(CALC=SCF\nDBOC=ON)\n')
        task = {
            'id': 'dboc', 'kind': 'external', 'backend': 'cfour',
            'input_name': 'ZMAT',
            'resources': {'cores': 1, 'memory_mb': 1024},
            'command': [sys.executable, '-c',
                        'import sys; sys.stdout.write(' + repr(_CFOUR_CH4) + ')'],
            'stdout': 'cfour.out', 'stderr': 'cfour.err',
            'required_outputs': ['cfour.out'],
            'success_marker': {'file': 'cfour.out',
                               'contains': 'The total diagonal Born-Oppenheimer correction'},
            'result_parser': {'kind': 'cfour_dboc', 'file': 'cfour.out',
                              'level': 'HF'},
        }
        details = _run_external(directory, {'task': task})
        assert details['returncode'] == 0
        assert details['parsed_result']['selected']['hartree'] == \
            pytest.approx(0.0025887093)
