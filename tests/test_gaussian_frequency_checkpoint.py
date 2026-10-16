"""Exercise frequency displacements without launching Gaussian."""
import numpy as np
import pytest
from ase.build import molecule
from ase.calculators.calculator import Calculator, all_changes

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot.frequencies import calc_vibrations


@pytest.mark.parametrize('guess, remaining', [
    ('Read', None), ('(Read,Mix,Always)', 'Mix,Always'), ('Mix', 'Mix')])
def test_displacements_do_not_read_the_removed_checkpoint(tmp_path, monkeypatch,
                                                         guess, remaining):
    monkeypatch.chdir(tmp_path)
    mol = molecule('H2')
    equilibrium = mol.positions.copy()
    mol.calc = Gaussian(label='optimization', method='AM1', chk='optimization',
                        guess=guess, force='',
                        command='/profile/g16 < PREFIX.com > PREFIX.log')
    calls = []

    def harmonic_forces(self, atoms=None, properties=None, system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        calls.append(dict(self.parameters))
        assert 'chk' not in self.parameters
        assert self.parameters.get('guess') == remaining
        assert self.parameters['method'] == 'AM1'
        assert self.command.startswith('/profile/g16 ')
        displacement = atoms.positions - equilibrium
        self.results = {'energy': .5 * np.sum(displacement ** 2),
                        'forces': -displacement}

    monkeypatch.setattr(Gaussian, 'calculate', harmonic_forces)
    freqs, zpe, hessian = calc_vibrations(mol, 'hydrogen')
    assert len(calls) > 1
    assert len(freqs) == 1 and np.isfinite(freqs).all()
    assert np.isfinite(zpe)
    assert hessian.shape == (6, 6)
    assert np.allclose(mol.positions, equilibrium)
    assert tmp_path.samefile('.')
