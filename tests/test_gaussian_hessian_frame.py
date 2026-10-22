import numpy as np
import pytest
from ase.build import molecule
from ase.units import Bohr

from kinbot.gaussian_hessian import read_formatted_hessian


def checkpoint(path, atoms, hessian):
    fields = {'Atomic numbers': atoms.numbers,
              'Current cartesian coordinates': (atoms.positions/Bohr).ravel(),
              'Cartesian Force Constants': hessian[np.tril_indices(len(hessian))]}
    path.write_text('\n'.join(f'{k} R N= {len(v)}\n' + ' '.join(map(str, v))
                              for k, v in fields.items())+'\n')


@pytest.mark.parametrize('formula', ['CH3', 'H2O', 'CO2', 'H2', 'H'])
def test_native_hessian_rotates_with_selected_geometry(tmp_path, formula):
    atoms = molecule(formula)
    # A rotation-covariant pair-spring Hessian also covers the arbitrary
    # transverse frame of a linear molecule and the atom's absent rotation.
    h = np.zeros((3*len(atoms),)*2)
    for i in range(len(atoms)):
        for j in range(i):
            v = atoms.positions[i]-atoms.positions[j]
            block = np.outer(v, v)/(v@v)
            for a, b, sign in [(i,i,1),(j,j,1),(i,j,-1),(j,i,-1)]:
                h[3*a:3*a+3,3*b:3*b+3] += sign*block
    rotation, _ = np.linalg.qr(np.random.default_rng(7).normal(size=(3,3)))
    target = atoms.positions@rotation + [4.,-3.,2.]
    transform = np.kron(np.eye(len(atoms)), rotation)
    p = tmp_path/'native.fchk'; checkpoint(p, atoms, h)
    actual = read_formatted_hessian(p, atoms.get_chemical_symbols(), target)
    np.testing.assert_allclose(actual, transform.T@h@transform, atol=1e-10)
    if len(atoms)>1:
        target[0,0]+=.01
        with pytest.raises(ValueError, match='geometry differs'):
            read_formatted_hessian(p, atoms.get_chemical_symbols(), target)


def test_incomplete_or_wrong_atom_checkpoint_refused(tmp_path):
    a=molecule('H2O'); p=tmp_path/'bad.fchk'; checkpoint(p,a,np.eye(9))
    with pytest.raises(ValueError, match='atom order'):
        read_formatted_hessian(p, ['S','H','H'], a.positions)
    p.write_text(p.read_text().replace('Current cartesian coordinates','Missing coordinates'))
    with pytest.raises(ValueError, match='coordinates'):
        read_formatted_hessian(p, a.get_chemical_symbols(), a.positions)
