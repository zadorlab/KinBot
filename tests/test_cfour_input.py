"""CFOUR keyword formatting must survive its fixed-width ZMAT reader."""

import pytest

from kinbot.anl.cfour import normalize_cfour_zmat, basis_keyword


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


@pytest.mark.parametrize('canonical,native', [
    ('aug-cc-pVDZ', 'AUG-PVDZ'), ('aug-cc-pVTZ', 'AUG-PVTZ'),
    ('aug-cc-pVQZ', 'AUG-PVQZ'), ('aug-cc-pV5Z', 'AUG-PV5Z'),
    ('aug-cc-pV6Z', 'AUG-PV6Z'), ('aug-cc-pCVTZ', 'AUG-PCVTZ'),
    ('cc-pVTZ', 'cc-pVTZ'), ('SPECIAL', 'SPECIAL'), ('custom-basis', 'custom-basis')])
def test_documented_cfour_basis_spelling(canonical, native):
    assert basis_keyword(canonical) == native
    zmat = normalize_cfour_zmat(f'Basis test\nH 0 0 0\n\n*CFOUR(CALC=CCSD,BASIS={canonical},DBOC=ON)\n')
    assert f'BASIS={native}\n' in zmat
    assert normalize_cfour_zmat(zmat) == zmat


def test_cfour_basis_parser_accepts_alias_but_not_a_different_basis():
    from tests.test_anl_results import _CFOUR_CH4, _CFOUR_CCSDTQ
    from kinbot.anl.results import parse_cfour_dboc, parse_cfour_energy
    text = _CFOUR_CH4.replace('cc-pVTZ', 'AUG-PVTZ')
    assert parse_cfour_dboc(text, basis='aug-cc-pVTZ')['selected_level'] == 'HF'
    for basis in ('cc-pVTZ', 'aug-cc-pVDZ', 'aug-cc-pCVTZ'):
        with pytest.raises(ValueError, match='echo basis'):
            parse_cfour_dboc(text, basis=basis)
    kwargs = dict(method='CCSDT(Q)', basis='aug-cc-pVTZ', reference='RHF',
                  correlation='unrestricted', core='frozen', program='cfour', driver='VCC')
    assert parse_cfour_energy(_CFOUR_CCSDTQ.replace('cc-pVTZ', 'AUG-PVTZ'), **kwargs)['basis'] == 'aug-cc-pVTZ'


def test_normalized_native_input_remains_consistent_with_canonical_parser():
    from kinbot.anl.results import validate_result_parser
    request = dict(kind='cfour_dboc', file='cfour.out', level='CCSD',
                   reference='UHF', basis='aug-cc-pVTZ')
    zmat = normalize_cfour_zmat('*CFOUR(CALC=CCSD,BASIS=aug-cc-pVTZ,DBOC=ON,REFERENCE=UHF)\n')
    validate_result_parser(request, backend='cfour', template=zmat, outputs=['cfour.out'])
    for changed in (zmat.replace('AUG-PVTZ', 'AUG-PVDZ'),
                    zmat.replace('AUG-PVTZ', 'AUG-PVTZ\nBASIS=AUG-PVTZ')):
        with pytest.raises(ValueError):
            validate_result_parser(request, backend='cfour', template=changed, outputs=['cfour.out'])
