"""Native cube layouts and coordinate frames must give the same VRC pivots."""

import numpy as np
import pytest
from ase.units import Bohr

from kinbot.orbital_cube import read_orbital_cube, lobe_directions


def write_cube(path, *, gaussian=False, phase=1, rotation=None, wrap=6,
               atom_number=16, atoms=None):
    rotation = np.eye(3) if rotation is None else rotation
    spacing = 0.2 * rotation
    origin = np.array([-1., -1., -1.]) @ rotation
    indices = np.indices((11, 11, 11)).reshape(3, -1).T
    local = -1. + indices * 0.2
    values = phase * local[:, 0] * np.exp(-8 * (local**2).sum(axis=1))
    atoms = [(atom_number, np.zeros(3))] if atoms is None else atoms
    lines = ['Orbital cube', 'Native Gaussian or Molpro layout',
             f'{-len(atoms) if gaussian else len(atoms)} '
             + ' '.join(map(str, origin / Bohr)) + (' 1' if gaussian else '')]
    lines += ['11 ' + ' '.join(map(str, v / Bohr)) for v in spacing]
    lines += [f'{number} 0 ' + ' '.join(map(str, xyz / Bohr))
              for number, xyz in atoms]
    if gaussian:
        lines.append('1 9')
    lines += [' '.join(map(str, values[i:i+wrap]))
              for i in range(0, len(values), wrap)]
    path.write_text('\n'.join(lines) + '\n')
    return path


@pytest.mark.parametrize('gaussian,wrap', [(True, 6), (False, 6), (False, 11)])
@pytest.mark.parametrize('phase', [1, -1])
def test_native_layout_and_orbital_phase_do_not_change_pivots(tmp_path, gaussian, wrap, phase):
    cube = read_orbital_cube(write_cube(
        tmp_path/'homo.cube', gaussian=gaussian, wrap=wrap, phase=phase))
    assert np.allclose(lobe_directions(cube, np.zeros(3)), [[-1, 0, 0], [1, 0, 0]])


def test_rotated_grid_and_other_elements(tmp_path):
    rotation = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
    cube = read_orbital_cube(write_cube(
        tmp_path/'homo.cube', rotation=rotation, atom_number=17))
    assert cube['atoms'].get_chemical_symbols() == ['Cl']
    assert np.allclose(lobe_directions(cube, np.zeros(3)), [[0, 0, -1], [0, 0, 1]])


def test_skew_grid_and_clipped_search_box():
    spacing = np.array([[.15, .05, 0.], [0., .15, 0.], [0., 0., .15]])
    shape = (5, 5, 5)
    center = np.array([2, 2, 2]) @ spacing
    xyz = np.indices(shape).reshape(3, -1).T @ spacing
    orbital = (xyz[:, 2] - center[2]).reshape(shape)
    cube = {'data': orbital, 'origin': np.zeros(3), 'spacing': spacing}
    directions = lobe_directions(cube, center)
    assert np.allclose(directions, [[0, 0, -1], [0, 0, 1]])
    with pytest.raises(ValueError, match='outside'):
        lobe_directions(cube, [100., 0., 0.])


def test_incomplete_and_nonfinite_cube_rejected(tmp_path):
    path = write_cube(tmp_path/'homo.cube')
    full = path.read_text()
    path.write_text(full.rsplit('\n', 2)[0] + '\n')
    with pytest.raises(ValueError, match='incomplete'):
        read_orbital_cube(path)
    path.write_text(full.replace('0.0', 'nan', 1))
    with pytest.raises(ValueError, match='invalid'):
        read_orbital_cube(path)


def test_fragment_atom_reordering_uses_cube_atom_center(tmp_path, monkeypatch):
    from kinbot.fragments import Fragment
    from kinbot.species_routing import routing_name
    # CO written O,C in the cube, but C,O in the fragment. Carbon lies on the
    # origin and is the requested reactive atom after reordering.
    fragment = Fragment(0, 1, ['C', 'O'], [[0., 0., 0.], [1.13, 0., 0.]],
                        [0], {'pp_orient': 'homo'}, str(tmp_path), [[0]],
                        charge=0, mult=1)
    directory = tmp_path/'vrctst'
    directory.mkdir()
    write_cube(directory/f'{routing_name(fragment)}_vts.cube',
               atoms=[(8, np.array([1.13, 0., 0.])), (6, np.zeros(3))])
    assert np.allclose(fragment.pp_from_homo(0), [[-1, 0, 0], [1, 0, 0]])
