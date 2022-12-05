from typing import Dict, List, Tuple, Union, Optional
from numbers import Real
from collections import namedtuple
import re
from string import digits
import numpy as np
from ase import Atoms
from ase.units import Angstrom, Bohr, nm


# split on newlines or semicolons
_re_linesplit = re.compile(r'\n|;')
# split definitions on whitespace or on "=" (possibly also with whitespace)
_re_defs = re.compile(r'\s*=\s*|\s+')


_ZMatrixRow = namedtuple(
    '_ZMatrixRow', 'ind1 dist ind2 a_bend ind3 a_dihedral',
)


ThreeFloats = Union[Tuple[float, float, float], np.ndarray]


class _ZMatrixToAtoms:
    known_units = dict(
        distance={'angstrom': Angstrom, 'bohr': Bohr, 'au': Bohr, 'nm': nm},
        angle={'radians': 1., 'degrees': np.pi / 180},
    )

    def __init__(self, dconv: Union[str, Real], aconv: Union[str, Real],
                 defs: Optional[Union[Dict[str, float],
                                str, List[str]]] = None) -> None:
        self.dconv = self.get_units('distance', dconv)  # type: float
        self.aconv = self.get_units('angle', aconv)  # type: float
        self.set_defs(defs)
        self.name_to_index: Optional[Dict[str, int]] = dict()
        self.symbols: List[str] = []
        self.positions: List[ThreeFloats] = []

    @property
    def nrows(self):
        return len(self.symbols)

    def get_units(self, kind: str, value: Union[str, Real]) -> float:
        if isinstance(value, Real):
            return float(value)
        out = self.known_units[kind].get(value.lower())
        if out is None:
            raise ValueError("Unknown {} units: {}"
                             .format(kind, value))
        return out

    def set_defs(self, defs: Union[Dict[str, float], str,
                                   List[str], None]) -> None:
        self.defs = dict()  # type: Dict[str, float]
        if defs is None:
            return

        if isinstance(defs, dict):
            self.defs.update(**defs)
            return

        if isinstance(defs, str):
            defs = _re_linesplit.split(defs.strip())

        for row in defs:
            key, val = _re_defs.split(row)
            self.defs[key] = self.get_var(val)

    def get_var(self, val: str) -> float:
        try:
            return float(val)
        except ValueError as e:
            val_out = self.defs.get(val.lstrip('+-'))
            if val_out is None:
                raise ValueError('Invalid value encountered in Z-matrix: {}'
                                 .format(val)) from e
        return val_out * (-1 if val.startswith('-') else 1)

    def get_index(self, name: str) -> int:
        """Find index for a given atom name"""
        try:
            return int(name) - 1
        except ValueError as e:
            if self.name_to_index is None or name not in self.name_to_index:
                raise ValueError('Failed to determine index for name "{}"'
                                 .format(name)) from e
        return self.name_to_index[name]

    def set_index(self, name: str) -> None:
        """Assign index to a given atom name for name -> index lookup"""
        if self.name_to_index is None:
            return

        if name in self.name_to_index:
            # "name" has been encountered before, so name_to_index is no
            # longer meaningful. Destroy the map.
            self.name_to_index = None
            return

        self.name_to_index[name] = self.nrows

    def validate_indices(self, *indices: int) -> None:
        """Raises an error if indices in a Z-matrix row are invalid."""
        if any(np.array(indices) >= self.nrows):
            raise ValueError('An invalid Z-matrix was provided! Row {} refers '
                             'to atom indices {}, at least one of which '
                             "hasn't been defined yet!"
                             .format(self.nrows, indices))

        if len(indices) != len(set(indices)):
            raise ValueError('An atom index has been used more than once a '
                             'row of the Z-matrix! Row numbers {}, '
                             'referred indices: {}'
                             .format(self.nrows, indices))

    def parse_row(self, row: str) -> Tuple[
            str, Union[_ZMatrixRow, ThreeFloats],
    ]:
        tokens = row.split()
        name = tokens[0]
        self.set_index(name)
        if len(tokens) == 1:
            assert self.nrows == 0
            return name, np.zeros(3, dtype=float)

        ind1 = self.get_index(tokens[1])
        if ind1 == -1:
            assert len(tokens) == 5
            return name, np.array(list(map(self.get_var, tokens[2:])),
                                  dtype=float)

        dist = self.dconv * self.get_var(tokens[2])

        if len(tokens) == 3:
            assert self.nrows == 1
            self.validate_indices(ind1)
            return name, np.array([dist, 0, 0], dtype=float)

        ind2 = self.get_index(tokens[3])
        a_bend = self.aconv * self.get_var(tokens[4])

        if len(tokens) == 5:
            assert self.nrows == 2
            self.validate_indices(ind1, ind2)
            return name, _ZMatrixRow(ind1, dist, ind2, a_bend, None, None)

        ind3 = self.get_index(tokens[5])
        a_dihedral = self.aconv * self.get_var(tokens[6])
        self.validate_indices(ind1, ind2, ind3)
        return name, _ZMatrixRow(ind1, dist, ind2, a_bend, ind3,
                                 a_dihedral)

    def add_atom(self, name: str, pos: ThreeFloats) -> None:
        """Sets the symbol and position of an atom."""
        self.symbols.append(
            ''.join([c for c in name if c not in digits]).capitalize()
        )
        self.positions.append(pos)

    def add_row(self, row: str) -> None:
        name, zrow = self.parse_row(row)

        if not isinstance(zrow, _ZMatrixRow):
            self.add_atom(name, zrow)
            return

        if zrow.ind3 is None:
            # This is the third atom, so only a bond distance and an angle
            # have been provided.
            pos = self.positions[zrow.ind1].copy()
            pos[0] += zrow.dist * np.cos(zrow.a_bend) * (zrow.ind2 - zrow.ind1)
            pos[1] += zrow.dist * np.sin(zrow.a_bend)
            self.add_atom(name, pos)
            return

        # ax1 is the dihedral axis, which is defined by the bond vector
        # between the two inner atoms in the dihedral, ind1 and ind2
        ax1 = self.positions[zrow.ind2] - self.positions[zrow.ind1]
        ax1 /= np.linalg.norm(ax1)

        # ax2 lies within the 1-2-3 plane, and it is perpendicular
        # to the dihedral axis
        ax2 = self.positions[zrow.ind2] - self.positions[zrow.ind3]
        ax2 -= ax1 * (ax2 @ ax1)
        ax2 /= np.linalg.norm(ax2)

        # ax3 is a vector that forms the appropriate dihedral angle, though
        # the bending angle is 90 degrees, rather than a_bend. It is formed
        # from a linear combination of ax2 and (ax2 x ax1)
        ax3 = (ax2 * np.cos(zrow.a_dihedral)
               + np.cross(ax2, ax1) * np.sin(zrow.a_dihedral))

        # The final position vector is a linear combination of ax1 and ax3.
        pos = ax1 * np.cos(zrow.a_bend) - ax3 * np.sin(zrow.a_bend)
        pos *= zrow.dist / np.linalg.norm(pos)
        pos += self.positions[zrow.ind1]
        self.add_atom(name, pos)

    def to_atoms(self) -> Atoms:
        return Atoms(self.symbols, self.positions)


def parse_zmatrix(zmat: Union[str, List[str]],
                  distance_units: Union[str, Real] = 'angstrom',
                  angle_units: Union[str, Real] = 'degrees',
                  defs: Optional[Union[Dict[str, float], str,
                                       List[str]]] = None) -> Atoms:
    """Converts a Z-matrix into an Atoms object.

    Parameters:

    zmat: Iterable or str
        The Z-matrix to be parsed. Iteration over `zmat` should yield the rows
        of the Z-matrix. If `zmat` is a str, it will be automatically split
        into a list at newlines.
    distance_units: str or float, optional
        The units of distance in the provided Z-matrix.
        Defaults to Angstrom.
    angle_units: str or float, optional
        The units for angles in the provided Z-matrix.
        Defaults to degrees.
    defs: dict or str, optional
        If `zmat` contains symbols for bond distances, bending angles, and/or
        dihedral angles instead of numeric values, then the definition of
        those symbols should be passed to this function using this keyword
        argument.
        Note: The symbol definitions are typically printed adjacent to the
        Z-matrix itself, but this function will not automatically separate
        the symbol definitions from the Z-matrix.

    Returns:

    atoms: Atoms object
    """
    zmatrix = _ZMatrixToAtoms(distance_units, angle_units, defs=defs)

    # zmat should be a list containing the rows of the z-matrix.
    # for convenience, allow block strings and split at newlines.
    if isinstance(zmat, str):
        zmat = _re_linesplit.split(zmat.strip())

    for row in zmat:
        zmatrix.add_row(row)

    return zmatrix.to_atoms()
