"""CFOUR keyword formatting must survive its fixed-width ZMAT reader."""

import pytest

from kinbot.anl.cfour import normalize_cfour_zmat


def test_wraps_prepared_single_line_without_changing_keyword_values():
    old = ('CH4\nC 0 0 0\n\n'
           '*CFOUR(CALC=SCF,BASIS=cc-pVTZ,DBOC=ON,COORD=CARTESIAN,'
           'UNITS=ANGSTROM,CHARGE=0,MULTIPLICITY=1,MEM_UNIT=MB,'
           'MEMORY_SIZE=11200)\n')
    wrapped = normalize_cfour_zmat(old)
    assert '*CFOUR(CALC=SCF\nBASIS=cc-pVTZ\nDBOC=ON\n' in wrapped
    assert 'CHARGE=0\nMULTIPLICITY=1\nMEM_UNIT=MB\nMEMORY_SIZE=11200)\n' in wrapped
    assert max(map(len, wrapped.splitlines())) <= 72
    assert normalize_cfour_zmat(wrapped) == wrapped


def test_nested_method_parentheses_remain_in_one_keyword():
    wrapped = normalize_cfour_zmat('*CFOUR(CALC=CCSD(T),BASIS=cc-pVTZ)\n')
    assert wrapped == '*CFOUR(CALC=CCSD(T)\nBASIS=cc-pVTZ)\n'


def test_oversize_single_keyword_fails_before_submission():
    with pytest.raises(ValueError, match='72-column'):
        normalize_cfour_zmat('*CFOUR(BASIS=' + 'X' * 80 + ')\n')
