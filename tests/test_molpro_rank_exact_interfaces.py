"""Regression across geometry, harmonic and energy readers for small species."""
from ase import Atoms
from ase.io import write
from ase.units import Bohr, Hartree, invcm, kJ, mol
import numpy as np
import pytest

from kinbot.ase_modules.calculators.molpro import (
    render_input, parse_output, parse_numerical_gradient)
from kinbot.anl.results import parse_molpro_harmonic, parse_molpro_energy


# Energy lines copied from Blodgett job 52917398. Its ! marks the final
# CCSD-F12b result; the CCSD(T)-F12b run prints the base without that marker.
H2_ENERGY = ''' !RHF-UCCSD-F12b energy                -1.174052556428
 SETTING KB_GEOM_ENERGY =        -1.17405256  AU
 KB_GEOM_ENERGY   =             -1.1740525564275448 AU
'''


@pytest.mark.parametrize('bang', ['!', ''])
def test_h2_final_f12b_label_and_precise_undisplaced_value(tmp_path, bang):
    p = tmp_path/'h2.out'
    p.write_text(H2_ENERGY.replace('!RHF', bang+'RHF') +
                 ' Computing numerical gradient using startcommand=RHF\n' +
                 H2_ENERGY.replace('-1.17405256', '-1.17415256') +
                 ' Molpro calculation terminated\n')
    assert parse_output(p, method='CCSD(T)-F12b', effective_method='CCSD-F12b',
                        require_precise_energy=True) == pytest.approx(
                            -1.1740525564275448 * Hartree, abs=1.e-12, rel=0)


@pytest.mark.parametrize('bad', ['F12a', 'displaced_only', 'no_gradient'])
def test_wrong_or_displaced_f12b_is_not_accepted(tmp_path, bad):
    text = H2_ENERGY
    if bad == 'F12a':
        text = text.replace('F12b', 'F12a')
    elif bad == 'displaced_only':
        text = text.splitlines(keepends=True)[1:] + text.splitlines(keepends=True)[:1]
        text = ''.join(text)
    else:
        # Incomplete energy calculation must not be accepted either.
        text += ' ERROR EXIT\n'
    p = tmp_path/'h2.out'
    p.write_text(text + ' Molpro calculation terminated\n')
    with pytest.raises(ValueError):
        parse_output(p, method='CCSD(T)-F12b', effective_method='CCSD-F12b',
                     require_precise_energy=True)


@pytest.mark.parametrize('f12', [False, True])
@pytest.mark.parametrize('symbols,charge,mult,electrons', [
    ('H', 0, 2, 1), ('H2', 1, 2, 1), ('H2', 0, 1, 2), ('He', 0, 1, 2)])
def test_small_electron_deck_and_energy_force_contract(tmp_path, f12, symbols, charge, mult, electrons):
    atoms = Atoms(symbols)
    if len(atoms) == 2:
        atoms[1].position[2] = .74
    method = 'CCSD(T)-F12b' if f12 else 'CCSD(T)'
    effective = 'HF' if electrons == 1 else 'CCSD-F12b' if f12 else 'CCSD'
    basis = 'aug-cc-pVTZ-F12' if f12 else 'cc-pVTZ'
    deck = render_input(atoms, basis=basis, geometry_name='point.xyz',
                        charge=charge, mult=mult, method=method, effective_method=effective)
    assert 'kb_geom_energy=energy\nshow[1,f28.16],kb_geom_energy\n' in deck
    assert 'forces,numerical,variable=kb_geom_energy,startcmd=rhf' in deck
    assert 'uccsd(t)' not in deck and 'scale_trip' not in deck
    assert ('{uccsd' in deck) == (electrons == 2)
    label = 'RHF STATE 1.1' if electrons == 1 else 'RHF-U'+effective
    native = f' !{label} energy -1.174052556428\n'
    p = tmp_path/'point.out'
    p.write_text(native + '\n'.join(H2_ENERGY.splitlines()[1:]) +
                 '\n Numerical gradient for KB_GEOM_ENERGY\n'
                 ' Atom dE/dx dE/dy dE/dz d2E/dx2 d2E/dy2 d2E/dz2\n' +
                 ''.join(f' {i+1} 0.0 0.0 0.0001 0.0 0.0 0.0\n' for i in range(len(atoms))) +
                 ' Molpro calculation terminated\n')
    write(tmp_path/'point.xyz', atoms, format='xyz')
    assert np.isfinite(parse_output(p, method=method, effective_method=effective,
                                   require_precise_energy=True))
    forces = parse_numerical_gradient(p, tmp_path/'point.xyz', atoms)
    np.testing.assert_allclose(forces[:, 2], -.0001*Hartree/Bohr)
    with pytest.raises(ValueError, match='electron count'):
        render_input(Atoms('CH4'), basis=basis, geometry_name='point.xyz',
                     method=method, effective_method=effective)


@pytest.mark.parametrize('f12', [False, True])
@pytest.mark.parametrize('electrons', [1, 2])
def test_rank_exact_frequency_and_single_point_readers(f12, electrons):
    method = 'CCSD(T)-F12b' if f12 else 'CCSD(T)'
    effective = 'HF' if electrons == 1 else 'CCSD-F12b' if f12 else 'CCSD'
    driver = 'RHF' if electrons == 1 else 'U'+effective.upper()
    command = '' if electrons == 1 else 'u'+effective.lower()+'\n'
    basis = 'aug-cc-pVTZ-F12' if f12 else 'cc-pVTZ'
    label = 'RHF STATE 1.1' if electrons == 1 else 'RHF-U'+effective
    prefix = (f'basis={basis}\nrhf\n{command}'
              f' PROGRAMS * TOTAL {driver} RHF-SCF INT\n'
              f' !{label} energy -1.174052556428\n')
    kwargs = dict(method=method, basis=basis, reference='ROHF' if electrons == 1 else 'RHF',
                  rank_exact_electrons=electrons, effective_method=effective)
    assert np.isfinite(parse_molpro_energy(prefix+' Molpro calculation terminated\n',
                                           **kwargs)['energy_hartree'])
    zpe = 2000*invcm/Hartree
    native = prefix + (f'frequencies,numerical\n'
        f' PROGRAM * FREQUENCIES (Calculation of harmonic vibrational spectra for {driver})\n'
        ' Vibration Wavenumber\n 1 4000.00\n\n'
        f' Zero point energy: {zpe:.10f} [H] 2000.00 [1/CM] {zpe*Hartree*mol/kJ:.5f} [KJ/MOL]\n'
        ' Molpro calculation terminated\n')
    assert parse_molpro_harmonic(native, **kwargs)['wavenumbers_cm_inverse'] == [4000.]
    with pytest.raises(ValueError, match='excitation rank'):
        parse_molpro_harmonic(native, **{**kwargs, 'rank_exact_electrons': 8})
    if f12 and electrons == 2:
        with pytest.raises(ValueError, match='F12b energy'):
            parse_molpro_harmonic(native.replace('!RHF-UCCSD-F12b', '!RHF-UCCSD-F12a'), **kwargs)
