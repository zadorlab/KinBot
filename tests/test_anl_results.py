"""Native result labels and units are part of the scientific acceptance gate."""

from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import pytest

from kinbot.anl.dispatch import _run_external
from kinbot.anl.results import (
    parse_cfour_dboc, parse_gaussian_vpt2, parse_molpro_energy,
    parse_molpro_harmonic, validate_result_parser,
)
from kinbot.energy import attach_anharmonic_frequencies


# The values and surrounding labels are from the completed Blodgett CH4 job.
_CFOUR_CH4 = """BASIS=cc-pVTZ
Total SCF energy including DBOC  -40.210746196622267
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
    assert parse_cfour_dboc(_CFOUR_CH4, basis='cc-pVTZ')['selected_level'] == 'HF'
    with pytest.raises(ValueError, match='basis'):
        parse_cfour_dboc(_CFOUR_CH4, basis='cc-pVQZ')
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


# Short excerpts preserve the exact labels and units seen in the completed
# native CH4 outputs; the licensed program files are not bundled with tests.
_MOLPRO_F12 = """basis=cc-pVTZ-F12
ccsd(t)-f12,scale_trip=1
kb_f12b=energy(2)
 !CCSD(T)-F12a total energy -40.458434986468
 !CCSD(T)-F12b total energy -40.454906199189
 CCSD(T)-F12/cc-pVTZ-F12 energy= -40.454906199189
 Molpro calculation terminated
"""

_MOLPRO_HARMONIC = """basis=cc-pVTZ
ccsd(t)
frequencies,numerical
 PROGRAM * FREQUENCIES (Calculation of harmonic vibrational spectra for CCSD(T))
   Low Vibration      Wavenumber
        Nr             [1/cm]
        1                0.00
     Vibration        Wavenumber
        Nr             [1/cm]
        1             1343.29
        2             1343.90
        3             1344.76
        4             1570.56
        5             1570.88
        6             3033.49
        7             3151.55
        8             3152.19
        9             3153.44

 Zero point energy:  0.04479801 [H]     9832.03 [1/CM]      117.62 [KJ/MOL]
 Molpro calculation terminated
"""

_GAUSSIAN_VPT2 = """#p B3LYP/cc-pVTZ Opt=(Tight,CalcFC) Freq=Anharmonic NoSymm
 WARNING: Unreliable CUBIC force constant i= 2,j= 1,k= 4
 Fundamental Bands
 -----------------
 Mode(n) Status      E(harm)   E(anharm)      I(harm)   I(anharm)
    1(1) active      3153.44     3010.18        0.10       0.12
 H  2(1) active      3152.19     3009.34        0.11       0.13
    3(1) active      3151.55     3008.91        0.12       0.14
    4(1) active      3033.49     2915.22        0.13       0.15
    5(1) active      1570.88     1534.77        0.14       0.16
    6(1) active      1570.56     1534.46        0.15       0.17
    7(1) active      1344.76     1312.63        0.16       0.18
    8(1) active      1343.90     1311.85        0.17       0.19
    9(1) active      1343.29     1311.24        0.18       0.20
 Overtones
 Anharmonic Zero Point Energy
 ----------------------------
 Harmonic       : cm-1 =  9783.68667 ; Kcal/mol =  27.973
 Anharmonic Pot.: cm-1 =  -137.81811 ; Kcal/mol =  -0.394
 Watson+Coriolis: cm-1 =     0.49625 ; Kcal/mol =   0.001
 Total Anharm   : cm-1 =  9646.36482 ; Kcal/mol =  27.580
 Normal termination of Gaussian 16 at Thu Sep 17 19:56:44 2026.
"""


def test_molpro_f12b_selects_exact_total_energy():
    result = parse_molpro_energy(_MOLPRO_F12, method='CCSD(T)-F12b',
                                 basis='cc-pVTZ-F12')
    assert result['energy_hartree'] == pytest.approx(-40.454906199189)
    assert result['method'] == 'CCSD(T)-F12b'
    without_fixture_variable = _MOLPRO_F12.replace('kb_f12b=energy(2)\n', '')
    assert parse_molpro_energy(without_fixture_variable, method='CCSD(T)-F12b',
                               basis='cc-pVTZ-F12')['energy_hartree'] == \
        pytest.approx(-40.454906199189)
    with pytest.raises(ValueError, match='summary'):
        parse_molpro_energy(_MOLPRO_F12.replace(
            'CCSD(T)-F12/cc-pVTZ-F12 energy= -40.454906199189',
            'CCSD(T)-F12/cc-pVTZ-F12 energy= -40.458434986468'),
            method='CCSD(T)-F12b', basis='cc-pVTZ-F12')
    with pytest.raises(ValueError, match='basis'):
        parse_molpro_energy(_MOLPRO_F12, method='CCSD(T)-F12b',
                            basis='cc-pVQZ-F12')


def test_molpro_harmonic_ignores_zero_modes_and_crosschecks_zpe():
    result = parse_molpro_harmonic(_MOLPRO_HARMONIC, basis='cc-pVTZ')
    assert len(result['wavenumbers_cm_inverse']) == 9
    assert result['low_modes_cm_inverse'] == [0.0]
    assert result['review_required'] is False
    assert result['wavenumbers_cm_inverse'][0] == 1343.29
    assert result['zpe']['hartree'] == pytest.approx(0.04479801)
    with pytest.raises(ValueError, match='modes disagree'):
        parse_molpro_harmonic(_MOLPRO_HARMONIC.replace('3153.44', '3253.44'),
                              basis='cc-pVTZ')
    with pytest.raises(ValueError, match='nonreal'):
        parse_molpro_harmonic(_MOLPRO_HARMONIC.replace('1343.29', '1343.29i'),
                              basis='cc-pVTZ')
    with pytest.raises(ValueError, match='imaginary low mode'):
        parse_molpro_harmonic(_MOLPRO_HARMONIC.replace(
            '1                0.00', '1              -20.00'),
            basis='cc-pVTZ')


def test_gaussian_vpt2_named_zpe_and_warnings():
    result = parse_gaussian_vpt2(_GAUSSIAN_VPT2, method='B3LYP',
                                 basis='cc-pVTZ')
    assert result['anharmonic_correction_cm_inverse'] == pytest.approx(-137.32185)
    assert result['optimized_in_job'] is True
    assert result['review_required'] is True
    assert len(result['warnings']) == 1
    assert len(result['anharmonic_fundamentals_cm_inverse']) == 9
    assert result['harmonic_fundamentals_cm_inverse'][0] == pytest.approx(3153.44)
    assert result['anharmonic_fundamentals_cm_inverse'][1] == pytest.approx(3009.34)
    with pytest.raises(ValueError, match='components disagree'):
        parse_gaussian_vpt2(_GAUSSIAN_VPT2.replace('9646.36482', '9600.36482'),
                             method='B3LYP', basis='cc-pVTZ')
    with pytest.raises(ValueError, match='final normal termination'):
        parse_gaussian_vpt2(_GAUSSIAN_VPT2.split(' Normal termination')[0],
                             method='B3LYP', basis='cc-pVTZ')
    with pytest.raises(ValueError, match='route'):
        parse_gaussian_vpt2(_GAUSSIAN_VPT2.replace('B3LYP/cc-pVTZ',
                                                  'B3LYP/cc-pVTZ-F12'),
                             method='B3LYP', basis='cc-pVTZ')
    frequency_only = _GAUSSIAN_VPT2.replace(
        'B3LYP/cc-pVTZ Opt=(Tight,CalcFC)',
        'B2PLYP/cc-pVTZ EmpiricalDispersion=GD3BJ')
    result = parse_gaussian_vpt2(frequency_only, method='B2PLYP',
                                 basis='cc-pVTZ', dispersion='GD3BJ')
    assert result['dispersion'] == 'GD3BJ'
    assert result['optimized_in_job'] is False
    with pytest.raises(ValueError, match='dispersion'):
        parse_gaussian_vpt2(frequency_only, method='B2PLYP',
                             basis='cc-pVTZ')


def test_gaussian_vpt2_fundamentals_are_an_explicit_mess_handoff():
    clean = _GAUSSIAN_VPT2.replace(
        ' WARNING: Unreliable CUBIC force constant i= 2,j= 1,k= 4\n', '')
    result = parse_gaussian_vpt2(clean, method='B3LYP', basis='cc-pVTZ')
    species = type('Species', (), {
        'reduced_freqs': list(reversed(
            result['harmonic_fundamentals_cm_inverse'])),
        'anl_thermochemistry_frequencies': None,
        'anl_thermochemistry_frequency_source': None,
    })()
    fundamentals = attach_anharmonic_frequencies(species, result)
    assert fundamentals[0] == pytest.approx(1311.24)
    assert fundamentals[-1] == pytest.approx(3010.18)
    source = species.anl_thermochemistry_frequency_source
    assert source['model'] == 'harmonic_plus_vpt2_mode_shifts'
    assert source['method'] == 'B3LYP'
    assert source['mode_map'] == list(range(9, 0, -1))
    # A rotor-projected target can contain fewer modes; the unmatched VPT2
    # mode is omitted rather than serialized alongside the hindered rotor.
    species.reduced_freqs = species.reduced_freqs[1:]
    assert len(attach_anharmonic_frequencies(species, result)) == 8
    species.reduced_freqs = [5000.] * 10
    with pytest.raises(ValueError, match='KinBot needs 10'):
        attach_anharmonic_frequencies(species, result)
    species.reduced_freqs = [5000.] * 9
    with pytest.raises(ValueError, match='mode match exceeds'):
        attach_anharmonic_frequencies(species, result)
    with pytest.raises(ValueError, match='quality review'):
        attach_anharmonic_frequencies(species, parse_gaussian_vpt2(
            _GAUSSIAN_VPT2, method='B3LYP', basis='cc-pVTZ'))


def test_gaussian_vpt2_rejects_incomplete_or_flags_nonpositive_fundamentals():
    with pytest.raises(ValueError, match='not consecutive'):
        parse_gaussian_vpt2(_GAUSSIAN_VPT2.replace('    9(1) active',
                                                   '   10(1) active'),
                            method='B3LYP', basis='cc-pVTZ')
    result = parse_gaussian_vpt2(
        _GAUSSIAN_VPT2.replace('1311.24', '-1.00'),
        method='B3LYP', basis='cc-pVTZ')
    assert result['review_required'] is True
    assert 'Nonpositive anharmonic fundamental mode(s): 9' in result['warnings']


def test_gaussian_linear_mode_labels_keep_degenerate_components():
    linear = (_GAUSSIAN_VPT2
              .replace('    5(1) active', '    5(1,-1) active')
              .replace('    6(1) active', '    5(1,+1) active')
              .replace('    7(1) active', '    6(1,-1) active')
              .replace('    8(1) active', '    6(1,+1) active')
              .replace('    9(1) active', '    7(1) active'))
    result = parse_gaussian_vpt2(linear, method='B3LYP', basis='cc-pVTZ')
    assert result['fundamental_mode_labels'][4:8] == [
        '5(1,-1)', '5(1,+1)', '6(1,-1)', '6(1,+1)']


def test_parser_declaration_must_match_method_and_input():
    request = {'kind': 'molpro_energy', 'file': 'sp.out',
               'method': 'CCSD(T)-F12b', 'basis': 'cc-pVTZ-F12'}
    template = 'basis=cc-pVTZ-F12\nccsd(t)-f12,scale_trip=1\n'
    validate_result_parser(request, backend='molpro', template=template,
                           outputs=['sp.out'])
    with pytest.raises(ValueError, match='invalid result_parser'):
        validate_result_parser(request, backend='molpro',
                               template=template.replace('scale_trip=1', 'scale_trip=0'),
                               outputs=['sp.out'])
