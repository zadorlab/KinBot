from warnings import warn

import numpy as np
from scipy.linalg import expm, logm
from ase.calculators.calculator import PropertyNotImplementedError
from ase.geometry import (find_mic, wrap_positions, get_distances_derivatives,
                          get_angles_derivatives, get_dihedrals_derivatives,
                          conditional_find_mic, get_angles, get_dihedrals)
from ase.utils.parsemath import eval_expression
from ase.stress import (full_3x3_to_voigt_6_stress,
                        voigt_6_to_full_3x3_stress)

__all__ = [
    'FixCartesian', 'FixBondLength', 'FixedMode',
    'FixAtoms', 'UnitCellFilter', 'ExpCellFilter',
    'FixScaled', 'StrainFilter', 'FixCom', 'FixedPlane', 'Filter',
    'FixConstraint', 'FixedLine', 'FixBondLengths', 'FixLinearTriatomic',
    'FixInternals', 'Hookean', 'ExternalForce', 'MirrorForce', 'MirrorTorque',
    "FixScaledParametricRelations", "FixCartesianParametricRelations"]


def dict2constraint(dct):
    if dct['name'] not in __all__:
        raise ValueError
    return globals()[dct['name']](**dct['kwargs'])


def slice2enlist(s, n):
    """Convert a slice object into a list of (new, old) tuples."""
    if isinstance(s, slice):
        return enumerate(range(*s.indices(n)))
    return enumerate(s)


def constrained_indices(atoms, only_include=None):
    """Returns a list of indices for the atoms that are constrained
    by a constraint that is applied.  By setting only_include to a
    specific type of constraint you can make it only look for that
    given constraint.
    """
    indices = []
    for constraint in atoms.constraints:
        if only_include is not None:
            if not isinstance(constraint, only_include):
                continue
        indices.extend(np.array(constraint.get_indices()))
    return np.array(np.unique(indices))


class FixConstraint:
    """Base class for classes that fix one or more atoms in some way."""

    def index_shuffle(self, atoms, ind):
        """Change the indices.

        When the ordering of the atoms in the Atoms object changes,
        this method can be called to shuffle the indices of the
        constraints.

        ind -- List or tuple of indices.

        """
        raise NotImplementedError

    def repeat(self, m, n):
        """ basic method to multiply by m, needs to know the length
        of the underlying atoms object for the assignment of
        multiplied constraints to work.
        """
        msg = ("Repeat is not compatible with your atoms' constraints."
               ' Use atoms.set_constraint() before calling repeat to '
               'remove your constraints.')
        raise NotImplementedError(msg)

    def adjust_momenta(self, atoms, momenta):
        """Adjusts momenta in identical manner to forces."""
        self.adjust_forces(atoms, momenta)

    def copy(self):
        return dict2constraint(self.todict().copy())


class IndexedConstraint(FixConstraint):
    def __init__(self, indices=None, mask=None):
        """Constrain chosen atoms.

        Parameters
        ----------
        indices : sequence of int
           Indices for those atoms that should be constrained.
        mask : sequence of bool
           One boolean per atom indicating if the atom should be
           constrained or not.
        """

        if mask is not None:
            if indices is not None:
                raise ValueError('Use only one of "indices" and "mask".')
            indices = mask
        indices = np.atleast_1d(indices)
        if np.ndim(indices) > 1:
            raise ValueError('indices has wrong amount of dimensions. '
                             f'Got {np.ndim(indices)}, expected ndim <= 1')

        if indices.dtype == bool:
            indices = np.arange(len(indices))[indices]
        elif len(indices) == 0:
            indices = np.empty(0, dtype=int)
        elif not np.issubdtype(indices.dtype, np.integer):
            raise ValueError('Indices must be integers or boolean mask, '
                             f'not dtype={indices.dtype}')

        if len(set(indices)) < len(indices):
            raise ValueError(
                'The indices array contains duplicates. '
                'Perhaps you want to specify a mask instead, but '
                'forgot the mask= keyword.')

        self.index = indices

    def index_shuffle(self, atoms, ind):
        # See docstring of superclass
        index = []

        # Resolve negative indices:
        actual_indices = set(np.arange(len(atoms))[self.index])

        for new, old in slice2enlist(ind, len(atoms)):
            if old in actual_indices:
                index.append(new)
        if len(index) == 0:
            raise IndexError('All indices in FixAtoms not part of slice')
        self.index = np.asarray(index, int)
        # XXX make immutable

    def get_indices(self):
        return self.index.copy()

    def repeat(self, m, n):
        i0 = 0
        natoms = 0
        if isinstance(m, int):
            m = (m, m, m)
        index_new = []
        for m2 in range(m[2]):
            for m1 in range(m[1]):
                for m0 in range(m[0]):
                    i1 = i0 + n
                    index_new += [i + natoms for i in self.index]
                    i0 = i1
                    natoms += n
        self.index = np.asarray(index_new, int)
        # XXX make immutable
        return self

    def delete_atoms(self, indices, natoms):
        """Removes atoms from the index array, if present.

        Required for removing atoms with existing constraint.
        """

        i = np.zeros(natoms, int) - 1
        new = np.delete(np.arange(natoms), indices)
        i[new] = np.arange(len(new))
        index = i[self.index]
        self.index = index[index >= 0]
        # XXX make immutable
        if len(self.index) == 0:
            return None
        return self


class FixAtoms(IndexedConstraint):
    """Fix chosen atoms.

    Examples
    --------
    Fix all Copper atoms:

    >>> mask = (atoms.symbols == 'Cu')
    >>> c = FixAtoms(mask=mask)
    >>> atoms.set_constraint(c)

    Fix all atoms with z-coordinate less than 1.0 Angstrom:

    >>> c = FixAtoms(mask=atoms.positions[:, 2] < 1.0)
    >>> atoms.set_constraint(c)
    """

    def get_removed_dof(self, atoms):
        return 3 * len(self.index)

    def adjust_positions(self, atoms, new):
        new[self.index] = atoms.positions[self.index]

    def adjust_forces(self, atoms, forces):
        forces[self.index] = 0.0

    def __repr__(self):
        clsname = type(self).__name__
        indices = ints2string(self.index)
        return f'{clsname}(indices={indices})'

    def todict(self):
        return {'name': 'FixAtoms',
                'kwargs': {'indices': self.index.tolist()}}


class FixCom(FixConstraint):
    """Constraint class for fixing the center of mass.

    References

    https://pubs.acs.org/doi/abs/10.1021/jp9722824

    """

    def get_removed_dof(self, atoms):
        return 3

    def adjust_positions(self, atoms, new):
        masses = atoms.get_masses()
        old_cm = atoms.get_center_of_mass()
        new_cm = np.dot(masses, new) / masses.sum()
        d = old_cm - new_cm
        new += d

    def adjust_forces(self, atoms, forces):
        m = atoms.get_masses()
        mm = np.tile(m, (3, 1)).T
        lb = np.sum(mm * forces, axis=0) / sum(m**2)
        forces -= mm * lb

    def todict(self):
        return {'name': 'FixCom',
                'kwargs': {}}


def ints2string(x, threshold=None):
    """Convert ndarray of ints to string."""
    if threshold is None or len(x) <= threshold:
        return str(x.tolist())
    return str(x[:threshold].tolist())[:-1] + ', ...]'


class FixBondLengths(FixConstraint):
    maxiter = 500

    def __init__(self, pairs, tolerance=1e-13,
                 bondlengths=None, iterations=None):
        """iterations:
                Ignored"""
        self.pairs = np.asarray(pairs)
        self.tolerance = tolerance
        self.bondlengths = bondlengths

    def get_removed_dof(self, atoms):
        return len(self.pairs)

    def adjust_positions(self, atoms, new):
        old = atoms.positions
        masses = atoms.get_masses()

        if self.bondlengths is None:
            self.bondlengths = self.initialize_bond_lengths(atoms)

        for i in range(self.maxiter):
            converged = True
            for j, ab in enumerate(self.pairs):
                a = ab[0]
                b = ab[1]
                cd = self.bondlengths[j]
                r0 = old[a] - old[b]
                d0, _ = find_mic(r0, atoms.cell, atoms.pbc)
                d1 = new[a] - new[b] - r0 + d0
                m = 1 / (1 / masses[a] + 1 / masses[b])
                x = 0.5 * (cd**2 - np.dot(d1, d1)) / np.dot(d0, d1)
                if abs(x) > self.tolerance:
                    new[a] += x * m / masses[a] * d0
                    new[b] -= x * m / masses[b] * d0
                    converged = False
            if converged:
                break
        else:
            raise RuntimeError('Did not converge')

    def adjust_momenta(self, atoms, p):
        old = atoms.positions
        masses = atoms.get_masses()

        if self.bondlengths is None:
            self.bondlengths = self.initialize_bond_lengths(atoms)

        for i in range(self.maxiter):
            converged = True
            for j, ab in enumerate(self.pairs):
                a = ab[0]
                b = ab[1]
                cd = self.bondlengths[j]
                d = old[a] - old[b]
                d, _ = find_mic(d, atoms.cell, atoms.pbc)
                dv = p[a] / masses[a] - p[b] / masses[b]
                m = 1 / (1 / masses[a] + 1 / masses[b])
                x = -np.dot(dv, d) / cd**2
                if abs(x) > self.tolerance:
                    p[a] += x * m * d
                    p[b] -= x * m * d
                    converged = False
            if converged:
                break
        else:
            raise RuntimeError('Did not converge')

    def adjust_forces(self, atoms, forces):
        self.constraint_forces = -forces
        self.adjust_momenta(atoms, forces)
        self.constraint_forces += forces

    def initialize_bond_lengths(self, atoms):
        bondlengths = np.zeros(len(self.pairs))

        for i, ab in enumerate(self.pairs):
            bondlengths[i] = atoms.get_distance(ab[0], ab[1], mic=True)

        return bondlengths

    def get_indices(self):
        return np.unique(self.pairs.ravel())

    def todict(self):
        return {'name': 'FixBondLengths',
                'kwargs': {'pairs': self.pairs.tolist(),
                           'tolerance': self.tolerance}}

    def index_shuffle(self, atoms, ind):
        """Shuffle the indices of the two atoms in this constraint"""
        map = np.zeros(len(atoms), int)
        map[ind] = 1
        n = map.sum()
        map[:] = -1
        map[ind] = range(n)
        pairs = map[self.pairs]
        self.pairs = pairs[(pairs != -1).all(1)]
        if len(self.pairs) == 0:
            raise IndexError('Constraint not part of slice')


def FixBondLength(a1, a2):
    """Fix distance between atoms with indices a1 and a2."""
    return FixBondLengths([(a1, a2)])


class FixLinearTriatomic(FixConstraint):
    """Holonomic constraints for rigid linear triatomic molecules."""

    def __init__(self, triples):
        """Apply RATTLE-type bond constraints between outer atoms n and m
           and linear vectorial constraints to the position of central
           atoms o to fix the geometry of linear triatomic molecules of the
           type:

           n--o--m

           Parameters:

           triples: list
               Indices of the atoms forming the linear molecules to constrain
               as triples. Sequence should be (n, o, m) or (m, o, n).

           When using these constraints in molecular dynamics or structure
           optimizations, atomic forces need to be redistributed within a
           triple. The function redistribute_forces_optimization implements
           the redistribution of forces for structure optimization, while
           the function redistribute_forces_md implements the redistribution
           for molecular dynamics.

           References:

           Ciccotti et al. Molecular Physics 47 (1982)
           :doi:`10.1080/00268978200100942`
        """
        self.triples = np.asarray(triples)
        if self.triples.shape[1] != 3:
            raise ValueError('"triples" has wrong size')
        self.bondlengths = None

    def get_removed_dof(self, atoms):
        return 4 * len(self.triples)

    @property
    def n_ind(self):
        return self.triples[:, 0]

    @property
    def m_ind(self):
        return self.triples[:, 2]

    @property
    def o_ind(self):
        return self.triples[:, 1]

    def initialize(self, atoms):
        masses = atoms.get_masses()
        self.mass_n, self.mass_m, self.mass_o = self.get_slices(masses)

        self.bondlengths = self.initialize_bond_lengths(atoms)
        self.bondlengths_nm = self.bondlengths.sum(axis=1)

        C1 = self.bondlengths[:, ::-1] / self.bondlengths_nm[:, None]
        C2 = (C1[:, 0] ** 2 * self.mass_o * self.mass_m +
              C1[:, 1] ** 2 * self.mass_n * self.mass_o +
              self.mass_n * self.mass_m)
        C2 = C1 / C2[:, None]
        C3 = self.mass_n * C1[:, 1] - self.mass_m * C1[:, 0]
        C3 = C2 * self.mass_o[:, None] * C3[:, None]
        C3[:, 1] *= -1
        C3 = (C3 + 1) / np.vstack((self.mass_n, self.mass_m)).T
        C4 = (C1[:, 0]**2 + C1[:, 1]**2 + 1)
        C4 = C1 / C4[:, None]

        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4

    def adjust_positions(self, atoms, new):
        old = atoms.positions
        new_n, new_m, new_o = self.get_slices(new)

        if self.bondlengths is None:
            self.initialize(atoms)

        r0 = old[self.n_ind] - old[self.m_ind]
        d0, _ = find_mic(r0, atoms.cell, atoms.pbc)
        d1 = new_n - new_m - r0 + d0
        a = np.einsum('ij,ij->i', d0, d0)
        b = np.einsum('ij,ij->i', d1, d0)
        c = np.einsum('ij,ij->i', d1, d1) - self.bondlengths_nm ** 2
        g = (b - (b**2 - a * c)**0.5) / (a * self.C3.sum(axis=1))
        g = g[:, None] * self.C3
        new_n -= g[:, 0, None] * d0
        new_m += g[:, 1, None] * d0
        if np.allclose(d0, r0):
            new_o = (self.C1[:, 0, None] * new_n
                     + self.C1[:, 1, None] * new_m)
        else:
            v1, _ = find_mic(new_n, atoms.cell, atoms.pbc)
            v2, _ = find_mic(new_m, atoms.cell, atoms.pbc)
            rb = self.C1[:, 0, None] * v1 + self.C1[:, 1, None] * v2
            new_o = wrap_positions(rb, atoms.cell, atoms.pbc)

        self.set_slices(new_n, new_m, new_o, new)

    def adjust_momenta(self, atoms, p):
        old = atoms.positions
        p_n, p_m, p_o = self.get_slices(p)

        if self.bondlengths is None:
            self.initialize(atoms)

        mass_nn = self.mass_n[:, None]
        mass_mm = self.mass_m[:, None]
        mass_oo = self.mass_o[:, None]

        d = old[self.n_ind] - old[self.m_ind]
        d, _ = find_mic(d, atoms.cell, atoms.pbc)
        dv = p_n / mass_nn - p_m / mass_mm
        k = np.einsum('ij,ij->i', dv, d) / self.bondlengths_nm ** 2
        k = self.C3 / (self.C3.sum(axis=1)[:, None]) * k[:, None]
        p_n -= k[:, 0, None] * mass_nn * d
        p_m += k[:, 1, None] * mass_mm * d
        p_o = (mass_oo * (self.C1[:, 0, None] * p_n / mass_nn +
                          self.C1[:, 1, None] * p_m / mass_mm))

        self.set_slices(p_n, p_m, p_o, p)

    def adjust_forces(self, atoms, forces):

        if self.bondlengths is None:
            self.initialize(atoms)

        A = self.C4 * np.diff(self.C1)
        A[:, 0] *= -1
        A -= 1
        B = np.diff(self.C4) / (A.sum(axis=1))[:, None]
        A /= (A.sum(axis=1))[:, None]

        self.constraint_forces = -forces
        old = atoms.positions

        fr_n, fr_m, fr_o = self.redistribute_forces_optimization(forces)

        d = old[self.n_ind] - old[self.m_ind]
        d, _ = find_mic(d, atoms.cell, atoms.pbc)
        df = fr_n - fr_m
        k = -np.einsum('ij,ij->i', df, d) / self.bondlengths_nm ** 2
        forces[self.n_ind] = fr_n + k[:, None] * d * A[:, 0, None]
        forces[self.m_ind] = fr_m - k[:, None] * d * A[:, 1, None]
        forces[self.o_ind] = fr_o + k[:, None] * d * B

        self.constraint_forces += forces

    def redistribute_forces_optimization(self, forces):
        """Redistribute forces within a triple when performing structure
        optimizations.

        The redistributed forces needs to be further adjusted using the
        appropriate Lagrange multipliers as implemented in adjust_forces."""
        forces_n, forces_m, forces_o = self.get_slices(forces)
        C1_1 = self.C1[:, 0, None]
        C1_2 = self.C1[:, 1, None]
        C4_1 = self.C4[:, 0, None]
        C4_2 = self.C4[:, 1, None]

        fr_n = ((1 - C4_1 * C1_1) * forces_n -
                C4_1 * (C1_2 * forces_m - forces_o))
        fr_m = ((1 - C4_2 * C1_2) * forces_m -
                C4_2 * (C1_1 * forces_n - forces_o))
        fr_o = ((1 - 1 / (C1_1**2 + C1_2**2 + 1)) * forces_o +
                C4_1 * forces_n + C4_2 * forces_m)

        return fr_n, fr_m, fr_o

    def redistribute_forces_md(self, atoms, forces, rand=False):
        """Redistribute forces within a triple when performing molecular
        dynamics.

        When rand=True, use the equations for random force terms, as
        used e.g. by Langevin dynamics, otherwise apply the standard
        equations for deterministic forces (see Ciccotti et al. Molecular
        Physics 47 (1982))."""
        if self.bondlengths is None:
            self.initialize(atoms)
        forces_n, forces_m, forces_o = self.get_slices(forces)
        C1_1 = self.C1[:, 0, None]
        C1_2 = self.C1[:, 1, None]
        C2_1 = self.C2[:, 0, None]
        C2_2 = self.C2[:, 1, None]
        mass_nn = self.mass_n[:, None]
        mass_mm = self.mass_m[:, None]
        mass_oo = self.mass_o[:, None]
        if rand:
            mr1 = (mass_mm / mass_nn) ** 0.5
            mr2 = (mass_oo / mass_nn) ** 0.5
            mr3 = (mass_nn / mass_mm) ** 0.5
            mr4 = (mass_oo / mass_mm) ** 0.5
        else:
            mr1 = 1.0
            mr2 = 1.0
            mr3 = 1.0
            mr4 = 1.0

        fr_n = ((1 - C1_1 * C2_1 * mass_oo * mass_mm) * forces_n -
                C2_1 * (C1_2 * mr1 * mass_oo * mass_nn * forces_m -
                        mr2 * mass_mm * mass_nn * forces_o))

        fr_m = ((1 - C1_2 * C2_2 * mass_oo * mass_nn) * forces_m -
                C2_2 * (C1_1 * mr3 * mass_oo * mass_mm * forces_n -
                        mr4 * mass_mm * mass_nn * forces_o))

        self.set_slices(fr_n, fr_m, 0.0, forces)

    def get_slices(self, a):
        a_n = a[self.n_ind]
        a_m = a[self.m_ind]
        a_o = a[self.o_ind]

        return a_n, a_m, a_o

    def set_slices(self, a_n, a_m, a_o, a):
        a[self.n_ind] = a_n
        a[self.m_ind] = a_m
        a[self.o_ind] = a_o

    def initialize_bond_lengths(self, atoms):
        bondlengths = np.zeros((len(self.triples), 2))

        for i in range(len(self.triples)):
            bondlengths[i, 0] = atoms.get_distance(self.n_ind[i],
                                                   self.o_ind[i], mic=True)
            bondlengths[i, 1] = atoms.get_distance(self.o_ind[i],
                                                   self.m_ind[i], mic=True)

        return bondlengths

    def get_indices(self):
        return np.unique(self.triples.ravel())

    def todict(self):
        return {'name': 'FixLinearTriatomic',
                'kwargs': {'triples': self.triples.tolist()}}

    def index_shuffle(self, atoms, ind):
        """Shuffle the indices of the three atoms in this constraint"""
        map = np.zeros(len(atoms), int)
        map[ind] = 1
        n = map.sum()
        map[:] = -1
        map[ind] = range(n)
        triples = map[self.triples]
        self.triples = triples[(triples != -1).all(1)]
        if len(self.triples) == 0:
            raise IndexError('Constraint not part of slice')


class FixedMode(FixConstraint):
    """Constrain atoms to move along directions orthogonal to
    a given mode only."""

    def __init__(self, mode):
        self.mode = (np.asarray(mode) / np.sqrt((mode**2).sum())).reshape(-1)

    def get_removed_dof(self, atoms):
        return len(atoms)

    def adjust_positions(self, atoms, newpositions):
        newpositions = newpositions.ravel()
        oldpositions = atoms.positions.ravel()
        step = newpositions - oldpositions
        newpositions -= self.mode * np.dot(step, self.mode)

    def adjust_forces(self, atoms, forces):
        forces = forces.ravel()
        forces -= self.mode * np.dot(forces, self.mode)

    def index_shuffle(self, atoms, ind):
        eps = 1e-12
        mode = self.mode.reshape(-1, 3)
        excluded = np.ones(len(mode), dtype=bool)
        excluded[ind] = False
        if (abs(mode[excluded]) > eps).any():
            raise IndexError('All nonzero parts of mode not in slice')
        self.mode = mode[ind].ravel()

    def get_indices(self):
        # This function will never properly work because it works on all
        # atoms and it has no idea how to tell how many atoms it is
        # attached to.  If it is being used, surely the user knows
        # everything is being constrained.
        return []

    def todict(self):
        return {'name': 'FixedMode',
                'kwargs': {'mode': self.mode.tolist()}}

    def __repr__(self):
        return 'FixedMode(%s)' % self.mode.tolist()


def _normalize(direction):
    if np.shape(direction) != (3,):
        raise ValueError("len(direction) is {len(direction)}. Has to be 3")

    direction = np.asarray(direction) / np.linalg.norm(direction)
    return direction


class FixedPlane(IndexedConstraint):
    """
    Constraint object for fixing chosen atoms to only move in a plane.

    The plane is defined by its normal vector *direction*
    """

    def __init__(self, indices, direction):
        """Constrain chosen atoms.

        Parameters
        ----------
        indices : int or list of int
            Index or indices for atoms that should be constrained
        direction : list of 3 int
            Direction of the normal vector

        Examples
        --------
        Fix all Copper atoms to only move in the yz-plane:

        >>> from ase.constraints import FixedPlane
        >>> c = FixedPlane(
        >>>     indices=[atom.index for atom in atoms if atom.symbol == 'Cu'],
        >>>     direction=[1, 0, 0],
        >>> )
        >>> atoms.set_constraint(c)

        or constrain a single atom with the index 0 to move in the xy-plane:

        >>> c = FixedPlane(indices=0, direction=[0, 0, 1])
        >>> atoms.set_constraint(c)
        """
        super().__init__(indices=indices)
        self.dir = _normalize(direction)

    def adjust_positions(self, atoms, newpositions):
        step = newpositions[self.index] - atoms.positions[self.index]
        newpositions[self.index] -= _projection(step, self.dir)

    def adjust_forces(self, atoms, forces):
        forces[self.index] -= _projection(forces[self.index], self.dir)

    def get_removed_dof(self, atoms):
        return len(self.index)

    def todict(self):
        return {
            'name': 'FixedPlane',
            'kwargs': {'indices': self.index.tolist(),
                       'direction': self.dir.tolist()}
        }

    def __repr__(self):
        return f'FixedPlane(indices={self.index}, {self.dir.tolist()})'


def _projection(vectors, direction):
    dotprods = vectors @ direction
    projection = direction[None, :] * dotprods[:, None]
    return projection


class FixedLine(IndexedConstraint):
    """
    Constrain an atom index or a list of atom indices to move on a line only.

    The line is defined by its vector *direction*
    """
    def __init__(self, indices, direction):
        """Constrain chosen atoms.

        Parameters
        ----------
        indices : int or list of int
            Index or indices for atoms that should be constrained
        direction : list of 3 int
            Direction of the vector defining the line

        Examples
        --------
        Fix all Copper atoms to only move in the x-direction:

        >>> from ase.constraints import FixedLine
        >>> c = FixedLine(
        >>>     indices=[atom.index for atom in atoms if atom.symbol == 'Cu'],
        >>>     direction=[1, 0, 0],
        >>> )
        >>> atoms.set_constraint(c)

        or constrain a single atom with the index 0 to move in the z-direction:

        >>> c = FixedLine(indices=0, direction=[0, 0, 1])
        >>> atoms.set_constraint(c)
        """
        super().__init__(indices)
        self.dir = _normalize(direction)

    def adjust_positions(self, atoms, newpositions):
        step = newpositions[self.index] - atoms.positions[self.index]
        projection = _projection(step, self.dir)
        newpositions[self.index] = atoms.positions[self.index] + projection

    def adjust_forces(self, atoms, forces):
        forces[self.index] = _projection(forces[self.index], self.dir)

    def get_removed_dof(self, atoms):
        return 2 * len(self.index)

    def __repr__(self):
        return f'FixedLine(indices={self.index}, {self.dir.tolist()})'

    def todict(self):
        return {
            'name': 'FixedLine',
            'kwargs': {'indices': self.index.tolist(),
                       'direction': self.dir.tolist()}
        }


class FixCartesian(IndexedConstraint):
    'Fix an atom index *a* in the directions of the cartesian coordinates.'

    def __init__(self, a, mask=(1, 1, 1)):
        super().__init__(indices=a)
        self.mask = ~np.asarray(mask, bool)

    def get_removed_dof(self, atoms):
        return (3 - self.mask.sum()) * len(self.index)

    def adjust_positions(self, atoms, new):
        step = new[self.index] - atoms.positions[self.index]
        step *= self.mask[None, :]
        new[self.index] = atoms.positions[self.index] + step

    def adjust_forces(self, atoms, forces):
        forces[self.index] *= self.mask[None, :]

    def __repr__(self):
        return 'FixCartesian(indices={}, mask={})'.format(
            self.index.tolist(), list(~self.mask))

    def todict(self):
        return {'name': 'FixCartesian',
                'kwargs': {'a': self.index.tolist(),
                           'mask': (~self.mask).tolist()}}


class FixScaled(IndexedConstraint):
    'Fix an atom index *a* in the directions of the unit vectors.'

    def __init__(self, a, mask=(1, 1, 1), cell=None):
        # XXX The unused cell keyword is there for compatibility
        # with old trajectory files.
        super().__init__(a)
        self.mask = np.array(mask, bool)

    def get_removed_dof(self, atoms):
        return self.mask.sum() * len(self.index)

    def adjust_positions(self, atoms, new):
        cell = atoms.cell
        scaled_old = cell.scaled_positions(atoms.positions[self.index])
        scaled_new = cell.scaled_positions(new[self.index])
        scaled_new[:, self.mask] = scaled_old[:, self.mask]
        new[self.index] = cell.cartesian_positions(scaled_new)

    def adjust_forces(self, atoms, forces):
        # Forces are covariant to the coordinate transformation,
        # use the inverse transformations
        cell = atoms.cell
        scaled_forces = cell.cartesian_positions(forces[self.index])
        scaled_forces *= -(self.mask - 1)
        forces[self.index] = cell.scaled_positions(scaled_forces)

    def todict(self):
        return {'name': 'FixScaled',
                'kwargs': {'a': self.index.tolist(),
                           'mask': self.mask.tolist()}}

    def __repr__(self):
        return 'FixScaled({}, {})'.format(self.index.tolist(), self.mask)


# TODO: Better interface might be to use dictionaries in place of very
# nested lists/tuples
class FixInternals(FixConstraint):
    """Constraint object for fixing multiple internal coordinates.

    Allows fixing bonds, angles, and dihedrals as well as linear combinations
    of bond lengths (bondcombos).
    Please provide angular units in degrees using angles_deg and
    dihedrals_deg.
    Fixing planar angles is not supported at the moment.
    """
    def __init__(self, bonds=None, angles=None, dihedrals=None,
                 angles_deg=None, dihedrals_deg=None,
                 bondcombos=None,
                 mic=False, epsilon=1.e-7):

        # deprecate public API using radians; degrees is preferred
        warn_msg = 'Please specify {} in degrees using the {} argument.'
        if angles:
            warn(FutureWarning(warn_msg.format('angles', 'angle_deg')))
            angles = np.asarray(angles)
            angles[:, 0] = angles[:, 0] / np.pi * 180
            angles = angles.tolist()
        else:
            angles = angles_deg
        if dihedrals:
            warn(FutureWarning(warn_msg.format('dihedrals', 'dihedrals_deg')))
            dihedrals = np.asarray(dihedrals)
            dihedrals[:, 0] = dihedrals[:, 0] / np.pi * 180
            dihedrals = dihedrals.tolist()
        else:
            dihedrals = dihedrals_deg

        self.bonds = bonds or []
        self.angles = angles or []
        self.dihedrals = dihedrals or []
        self.bondcombos = bondcombos or []
        self.mic = mic
        self.epsilon = epsilon

        self.n = (len(self.bonds) + len(self.angles) + len(self.dihedrals)
                  + len(self.bondcombos))

        # Initialize these at run-time:
        self.constraints = []
        self.initialized = False

    def get_removed_dof(self, atoms):
        return self.n

    def initialize(self, atoms):
        if self.initialized:
            return
        masses = np.repeat(atoms.get_masses(), 3)
        cell = None
        pbc = None
        if self.mic:
            cell = atoms.cell
            pbc = atoms.pbc
        self.constraints = []
        for data, make_constr in [(self.bonds, self.FixBondLengthAlt),
                                  (self.angles, self.FixAngle),
                                  (self.dihedrals, self.FixDihedral),
                                  (self.bondcombos, self.FixBondCombo)]:
            for datum in data:
                constr = make_constr(datum[0], datum[1], masses, cell, pbc)
                self.constraints.append(constr)
        self.initialized = True

    def shuffle_definitions(self, shuffle_dic, internal_type):
        dfns = []  # definitions
        for dfn in internal_type:  # e.g. for bond in self.bonds
            append = True
            new_dfn = [dfn[0], list(dfn[1])]
            for old in dfn[1]:
                if old in shuffle_dic:
                    new_dfn[1][dfn[1].index(old)] = shuffle_dic[old]
                else:
                    append = False
                    break
            if append:
                dfns.append(new_dfn)
        return dfns

    def shuffle_combos(self, shuffle_dic, internal_type):
        dfns = []  # definitions
        for dfn in internal_type:  # e.g. for bondcombo in self.bondcombos
            append = True
            all_indices = [idx[0:-1] for idx in dfn[1]]
            new_dfn = [dfn[0], list(dfn[1])]
            for i, indices in enumerate(all_indices):
                for old in indices:
                    if old in shuffle_dic:
                        new_dfn[1][i][indices.index(old)] = shuffle_dic[old]
                    else:
                        append = False
                        break
                if not append:
                    break
            if append:
                dfns.append(new_dfn)
        return dfns

    def index_shuffle(self, atoms, ind):
        # See docstring of superclass
        self.initialize(atoms)
        shuffle_dic = dict(slice2enlist(ind, len(atoms)))
        shuffle_dic = {old: new for new, old in shuffle_dic.items()}
        self.bonds = self.shuffle_definitions(shuffle_dic, self.bonds)
        self.angles = self.shuffle_definitions(shuffle_dic, self.angles)
        self.dihedrals = self.shuffle_definitions(shuffle_dic, self.dihedrals)
        self.bondcombos = self.shuffle_combos(shuffle_dic, self.bondcombos)
        self.initialized = False
        self.initialize(atoms)
        if len(self.constraints) == 0:
            raise IndexError('Constraint not part of slice')

    def get_indices(self):
        cons = []
        for dfn in self.bonds + self.dihedrals + self.angles:
            cons.extend(dfn[1])
        for dfn in self.bondcombos:
            for partial_dfn in dfn[1]:
                cons.extend(partial_dfn[0:-1])  # last index is the coefficient
        return list(set(cons))

    def todict(self):
        return {'name': 'FixInternals',
                'kwargs': {'bonds': self.bonds,
                           'angles_deg': self.angles,
                           'dihedrals_deg': self.dihedrals,
                           'bondcombos': self.bondcombos,
                           'mic': self.mic,
                           'epsilon': self.epsilon}}

    def adjust_positions(self, atoms, newpos):
        self.initialize(atoms)
        for constraint in self.constraints:
            constraint.setup_jacobian(atoms.positions)
        for j in range(50):
            maxerr = 0.0
            for constraint in self.constraints:
                constraint.adjust_positions(atoms.positions, newpos)
                maxerr = max(abs(constraint.sigma), maxerr)
            if maxerr < self.epsilon:
                return
        msg = 'FixInternals.adjust_positions did not converge.'
        if any([constr.targetvalue > 175. or constr.targetvalue < 5. for constr
                in self.constraints if type(constr) is self.FixAngle]):
            msg += (' This may be caused by an almost planar angle.'
                    ' Support for planar angles would require the'
                    ' implementation of ghost, i.e. dummy, atoms.'
                    ' See issue #868.')
        raise ValueError(msg)

    def adjust_forces(self, atoms, forces):
        """Project out translations and rotations and all other constraints"""
        self.initialize(atoms)
        positions = atoms.positions
        N = len(forces)
        list2_constraints = list(np.zeros((6, N, 3)))
        tx, ty, tz, rx, ry, rz = list2_constraints

        list_constraints = [r.ravel() for r in list2_constraints]

        tx[:, 0] = 1.0
        ty[:, 1] = 1.0
        tz[:, 2] = 1.0
        ff = forces.ravel()

        # Calculate the center of mass
        center = positions.sum(axis=0) / N

        rx[:, 1] = -(positions[:, 2] - center[2])
        rx[:, 2] = positions[:, 1] - center[1]
        ry[:, 0] = positions[:, 2] - center[2]
        ry[:, 2] = -(positions[:, 0] - center[0])
        rz[:, 0] = -(positions[:, 1] - center[1])
        rz[:, 1] = positions[:, 0] - center[0]

        # Normalizing transl., rotat. constraints
        for r in list2_constraints:
            r /= np.linalg.norm(r.ravel())

        # Add all angle, etc. constraint vectors
        for constraint in self.constraints:
            constraint.setup_jacobian(positions)
            constraint.adjust_forces(positions, forces)
            list_constraints.insert(0, constraint.jacobian)
        # QR DECOMPOSITION - GRAM SCHMIDT

        list_constraints = [r.ravel() for r in list_constraints]
        aa = np.column_stack(list_constraints)
        (aa, bb) = np.linalg.qr(aa)
        # Projection
        hh = []
        for i, constraint in enumerate(self.constraints):
            hh.append(aa[:, i] * np.row_stack(aa[:, i]))

        txx = aa[:, self.n] * np.row_stack(aa[:, self.n])
        tyy = aa[:, self.n + 1] * np.row_stack(aa[:, self.n + 1])
        tzz = aa[:, self.n + 2] * np.row_stack(aa[:, self.n + 2])
        rxx = aa[:, self.n + 3] * np.row_stack(aa[:, self.n + 3])
        ryy = aa[:, self.n + 4] * np.row_stack(aa[:, self.n + 4])
        rzz = aa[:, self.n + 5] * np.row_stack(aa[:, self.n + 5])
        T = txx + tyy + tzz + rxx + ryy + rzz
        for vec in hh:
            T += vec
        ff = np.dot(T, np.row_stack(ff))
        forces[:, :] -= np.dot(T, np.row_stack(ff)).reshape(-1, 3)

    def __repr__(self):
        constraints = repr(self.constraints)
        return 'FixInternals(_copy_init=%s, epsilon=%s)' % (constraints,
                                                            repr(self.epsilon))

    def __str__(self):
        return '\n'.join([repr(c) for c in self.constraints])

    # Classes for internal use in FixInternals
    class FixInternalsBase:
        """Base class for subclasses of FixInternals."""
        def __init__(self, targetvalue, indices, masses, cell, pbc):
            self.targetvalue = targetvalue  # constant target value
            self.indices = [defin[0:-1] for defin in indices]  # indices, defs
            self.coefs = np.asarray([defin[-1] for defin in indices])  # coefs
            self.masses = masses
            self.jacobian = []  # geometric Jacobian matrix, Wilson B-matrix
            self.sigma = 1.  # difference between current and target value
            self.projected_force = None  # helps optimizers scan along constr.
            self.cell = cell
            self.pbc = pbc

        def finalize_jacobian(self, pos, n_internals, n, derivs):
            """Populate jacobian with derivatives for `n_internals` defined
            internals. n = 2 (bonds), 3 (angles), 4 (dihedrals)."""
            jacobian = np.zeros((n_internals, *pos.shape))
            for i, idx in enumerate(self.indices):
                for j in range(n):
                    jacobian[i, idx[j]] = derivs[i, j]
            jacobian = jacobian.reshape((n_internals, 3 * len(pos)))
            return self.coefs @ jacobian

        def finalize_positions(self, newpos):
            jacobian = self.jacobian / self.masses
            lamda = -self.sigma / (jacobian @ self.get_jacobian(newpos))
            dnewpos = lamda * jacobian
            newpos += dnewpos.reshape(newpos.shape)

        def adjust_forces(self, positions, forces):
            self.projected_force = np.dot(self.jacobian, forces.ravel())
            self.jacobian /= np.linalg.norm(self.jacobian)

    class FixBondCombo(FixInternalsBase):
        """Constraint subobject for fixing linear combination of bond lengths
        within FixInternals.

        sum_i( coef_i * bond_length_i ) = constant
        """
        def get_jacobian(self, pos):
            bondvectors = [pos[k] - pos[h] for h, k in self.indices]
            derivs = get_distances_derivatives(bondvectors, cell=self.cell,
                                               pbc=self.pbc)
            return self.finalize_jacobian(pos, len(bondvectors), 2, derivs)

        def setup_jacobian(self, pos):
            self.jacobian = self.get_jacobian(pos)

        def adjust_positions(self, oldpos, newpos):
            bondvectors = [newpos[k] - newpos[h] for h, k in self.indices]
            (_, ), (dists, ) = conditional_find_mic([bondvectors],
                                                    cell=self.cell,
                                                    pbc=self.pbc)
            value = self.coefs @ dists
            self.sigma = value - self.targetvalue
            self.finalize_positions(newpos)

        def __repr__(self):
            return 'FixBondCombo({}, {}, {})'.format(repr(self.targetvalue),
                                                     self.indices, self.coefs)

    class FixBondLengthAlt(FixBondCombo):
        """Constraint subobject for fixing bond length within FixInternals.
        Fix distance between atoms with indices a1, a2."""
        def __init__(self, targetvalue, indices, masses, cell, pbc):
            if targetvalue <= 0.: 
                raise ZeroDivisionError('Invalid targetvalue for fixed bond')
            indices = [list(indices) + [1.]]  # bond definition with coef 1.
            super().__init__(targetvalue, indices, masses, cell=cell, pbc=pbc)

        def __repr__(self):
            return 'FixBondLengthAlt({}, {})'.format(self.targetvalue,
                                                     *self.indices)

    class FixAngle(FixInternalsBase):
        """Constraint subobject for fixing an angle within FixInternals.

        Convergence is potentially problematic for angles very close to
        0 or 180 degrees as there is a singularity in the Cartesian derivative.
        Fixing planar angles is therefore not supported at the moment.
        """
        def __init__(self, targetvalue, indices, masses, cell, pbc):
            """Fix atom movement to construct a constant angle."""
            if targetvalue <= 0. or targetvalue >= 180.: 
                raise ZeroDivisionError('Invalid targetvalue for fixed angle')
            indices = [list(indices) + [1.]]  # angle definition with coef 1.
            super().__init__(targetvalue, indices, masses, cell=cell, pbc=pbc)

        def gather_vectors(self, pos):
            v0 = [pos[h] - pos[k] for h, k, l in self.indices]
            v1 = [pos[l] - pos[k] for h, k, l in self.indices]
            return v0, v1

        def get_jacobian(self, pos):
            v0, v1 = self.gather_vectors(pos)
            derivs = get_angles_derivatives(v0, v1, cell=self.cell,
                                            pbc=self.pbc)
            return self.finalize_jacobian(pos, len(v0), 3, derivs)

        def setup_jacobian(self, pos):
            self.jacobian = self.get_jacobian(pos)

        def adjust_positions(self, oldpos, newpos):
            v0, v1 = self.gather_vectors(newpos)
            value = get_angles(v0, v1, cell=self.cell, pbc=self.pbc)
            self.sigma = value - self.targetvalue
            self.finalize_positions(newpos)

        def __repr__(self):
            return 'FixAngle({}, {})'.format(self.targetvalue, *self.indices)

    class FixDihedral(FixInternalsBase):
        """Constraint subobject for fixing a dihedral angle within FixInternals.

        A dihedral becomes undefined when at least one of the inner two angles
        becomes planar. Make sure to avoid this situation.
        """
        def __init__(self, targetvalue, indices, masses, cell, pbc):
            indices = [list(indices) + [1.]]  # dihedral def. with coef 1.
            super().__init__(targetvalue, indices, masses, cell=cell, pbc=pbc)

        def gather_vectors(self, pos):
            v0 = [pos[k] - pos[h] for h, k, l, m in self.indices]
            v1 = [pos[l] - pos[k] for h, k, l, m in self.indices]
            v2 = [pos[m] - pos[l] for h, k, l, m in self.indices]
            return v0, v1, v2

        def get_jacobian(self, pos):
            v0, v1, v2 = self.gather_vectors(pos)
            derivs = get_dihedrals_derivatives(v0, v1, v2, cell=self.cell,
                                               pbc=self.pbc)
            return self.finalize_jacobian(pos, len(v0), 4, derivs)

        def setup_jacobian(self, pos):
            self.jacobian = self.get_jacobian(pos)

        def adjust_positions(self, oldpos, newpos):
            v0, v1, v2 = self.gather_vectors(newpos)
            value = get_dihedrals(v0, v1, v2, cell=self.cell, pbc=self.pbc)
            # apply minimum dihedral difference 'convention': (diff <= 180)
            self.sigma = (value - self.targetvalue + 180) % 360 - 180
            self.finalize_positions(newpos)

        def __repr__(self):
            return 'FixDihedral({}, {})'.format(self.targetvalue, *self.indices)


class FixParametricRelations(FixConstraint):

    def __init__(
        self,
        indices,
        Jacobian,
        const_shift,
        params=None,
        eps=1e-12,
        use_cell=False,
    ):
        """Constrains the degrees of freedom to act in a reduced parameter space defined by the Jacobian

        These constraints are based off the work in: https://arxiv.org/abs/1908.01610

        The constraints linearly maps the full 3N degrees of freedom, where N is number of active
        lattice vectors/atoms onto a reduced subset of M free parameters, where M <= 3*N. The
        Jacobian matrix and constant shift vector map the full set of degrees of freedom onto the
        reduced parameter space.

        Currently the constraint is set up to handle either atomic positions or lattice vectors
        at one time, but not both. To do both simply add a two constraints for each set. This is
        done to keep the mathematics behind the operations separate.

        It would be possible to extend these constraints to allow non-linear transformations
        if functionality to update the Jacobian at each position update was included. This would
        require passing an update function evaluate it every time adjust_positions is callled.
        This is currently NOT supported, and there are no plans to implement it in the future.

        Args:
            indices (list of int): indices of the constrained atoms
                (if not None or empty then cell_indices must be None or Empty)
            Jacobian (np.ndarray(shape=(3*len(indices), len(params)))): The Jacobian describing
                the parameter space transformation
            const_shift (np.ndarray(shape=(3*len(indices)))): A vector describing the constant term
                in the transformation not accounted for in the Jacobian
            params (list of str): parameters used in the parametric representation
                if None a list is generated based on the shape of the Jacobian
            eps (float): a small number to compare the similarity of numbers and set the precision used
                to generate the constraint expressions
            use_cell (bool): if True then act on the cell object
        """
        self.indices = np.array(indices)
        self.Jacobian = np.array(Jacobian)
        self.const_shift = np.array(const_shift)

        assert self.const_shift.shape[0] == 3*len(self.indices)
        assert self.Jacobian.shape[0] == 3*len(self.indices)

        self.eps = eps
        self.use_cell = use_cell

        if params is None:
            params = []
            if self.Jacobian.shape[1] > 0:
                int_fmt_str = "{:0" + str(int(np.ceil(np.log10(self.Jacobian.shape[1])))) + "d}"
                for param_ind in range(self.Jacobian.shape[1]):
                    params.append("param_" + int_fmt_str.format(param_ind))
        else:
            assert len(params) == self.Jacobian.shape[-1]

        self.params = params

        self.Jacobian_inv = np.linalg.inv(self.Jacobian.T @ self.Jacobian) @ self.Jacobian.T

    @classmethod
    def from_expressions(cls, indices, params, expressions, eps=1e-12, use_cell=False):
        """Converts the expressions into a Jacobian Matrix/const_shift vector and constructs a FixParametricRelations constraint

        The expressions must be a list like object of size 3*N and elements must be ordered as:
        [n_0,i; n_0,j; n_0,k; n_1,i; n_1,j; .... ; n_N-1,i; n_N-1,j; n_N-1,k],
        where i, j, and k are the first, second and third component of the atomic position/lattice
        vector. Currently only linear operations are allowed to be included in the expressions so
        only terms like:
            - const * param_0
            - sqrt[const] * param_1
            - const * param_0 +/- const * param_1 +/- ... +/- const * param_M
        where const is any real number and param_0, param_1, ..., param_M are the parameters passed in
        params, are allowed.

        For example, the fractional atomic position constraints for wurtzite are:
        params = ["z1", "z2"]
        expressions = [
            "1.0/3.0", "2.0/3.0", "z1",
            "2.0/3.0", "1.0/3.0", "0.5 + z1",
            "1.0/3.0", "2.0/3.0", "z2",
            "2.0/3.0", "1.0/3.0", "0.5 + z2",
        ]

        For diamond are:
        params = []
        expressions = [
            "0.0", "0.0", "0.0",
            "0.25", "0.25", "0.25",
        ],

        and for stannite are
        params=["x4", "z4"]
        expressions = [
            "0.0", "0.0", "0.0",
            "0.0", "0.5", "0.5",
            "0.75", "0.25", "0.5",
            "0.25", "0.75", "0.5",
            "x4 + z4", "x4 + z4", "2*x4",
            "x4 - z4", "x4 - z4", "-2*x4",
             "0.0", "-1.0 * (x4 + z4)", "x4 - z4",
             "0.0", "x4 - z4", "-1.0 * (x4 + z4)",
        ]

        Args:
            indices (list of int): indices of the constrained atoms
                (if not None or empty then cell_indices must be None or Empty)
            params (list of str): parameters used in the parametric representation
            expressions (list of str): expressions used to convert from the parametric to the real space
                representation
            eps (float): a small number to compare the similarity of numbers and set the precision used
                to generate the constraint expressions
            use_cell (bool): if True then act on the cell object

        Returns:
            cls(
                indices,
                Jacobian generated from expressions,
                const_shift generated from expressions,
                params,
                eps-12,
                use_cell,
            )
        """
        Jacobian = np.zeros((3*len(indices), len(params)))
        const_shift = np.zeros(3*len(indices))

        for expr_ind, expression in enumerate(expressions):
            expression = expression.strip()

            # Convert subtraction to addition
            expression = expression.replace("-", "+(-1.0)*")
            if "+" == expression[0]:
                expression = expression[1:]
            elif "(+" == expression[:2]:
                expression = "(" + expression[2:]

            # Explicitly add leading zeros so when replacing param_1 with 0.0 param_11 does not become 0.01
            int_fmt_str = "{:0" + str(int(np.ceil(np.log10(len(params)+1)))) + "d}"

            param_dct = dict()
            param_map = dict()

            # Construct a standardized param template for A/B filling
            for param_ind, param in enumerate(params):
                param_str = "param_" + int_fmt_str.format(param_ind)
                param_map[param] = param_str
                param_dct[param_str] = 0.0

            # Replace the parameters according to the map
            # Sort by string length (long to short) to prevent cases like x11 becoming f"{param_map["x1"]}1"
            for param in sorted(params, key=lambda s: -1.0*len(s)):
                expression = expression.replace(param, param_map[param])

            # Partial linearity check
            for express_sec in expression.split("+"):
                in_sec = [param in express_sec for param in param_dct]
                n_params_in_sec = len(np.where(np.array(in_sec))[0])
                if n_params_in_sec > 1:
                    raise ValueError("The FixParametricRelations expressions must be linear.")

            const_shift[expr_ind] = float(eval_expression(expression, param_dct))

            for param_ind in range(len(params)):
                param_str = "param_" + int_fmt_str.format(param_ind)
                if param_str not in expression:
                    Jacobian[expr_ind, param_ind] = 0.0
                    continue
                param_dct[param_str] = 1.0
                test_1 = float(eval_expression(expression, param_dct))
                test_1 -= const_shift[expr_ind]
                Jacobian[expr_ind, param_ind] = test_1

                param_dct[param_str] = 2.0
                test_2 = float(eval_expression(expression, param_dct))
                test_2 -= const_shift[expr_ind]
                if abs(test_2 / test_1 - 2.0) > eps:
                    raise ValueError("The FixParametricRelations expressions must be linear.")
                param_dct[param_str] = 0.0

        args = [
            indices,
            Jacobian,
            const_shift,
            params,
            eps,
            use_cell,
        ]
        if cls is FixScaledParametricRelations:
            args = args[:-1]
        return cls(*args)

    @property
    def expressions(self):
        """Generate the expressions represented by the current self.Jacobian and self.const_shift objects"""
        expressions = []
        per = int(round(-1 * np.log10(self.eps)))
        fmt_str = "{:." + str(per + 1) + "g}"
        for index, shift_val in enumerate(self.const_shift):
            exp = ""
            if np.all(np.abs(self.Jacobian[index]) < self.eps) or np.abs(shift_val) > self.eps:
                exp += fmt_str.format(shift_val)

            param_exp = ""
            for param_index, jacob_val in enumerate(self.Jacobian[index]):
                abs_jacob_val = np.round(np.abs(jacob_val), per + 1)
                if abs_jacob_val < self.eps:
                    continue

                param = self.params[param_index]
                if param_exp or exp:
                    if jacob_val > -1.0*self.eps:
                        param_exp += " + "
                    else:
                        param_exp += " - "
                elif (not exp) and (not param_exp) and (jacob_val < -1.0*self.eps):
                    param_exp += "-"

                if np.abs(abs_jacob_val-1.0) <= self.eps:
                    param_exp += "{:s}".format(param)
                else:
                    param_exp += (fmt_str + "*{:s}").format(abs_jacob_val, param)

            exp += param_exp

            expressions.append(exp)
        return np.array(expressions).reshape((-1, 3))

    def todict(self):
        """Create a dictionary representation of the constraint"""
        return {
            "name": type(self).__name__,
            "kwargs": {
                "indices": self.indices,
                "params": self.params,
                "Jacobian": self.Jacobian,
                "const_shift": self.const_shift,
                "eps": self.eps,
                "use_cell": self.use_cell,
            }
        }

    def __repr__(self):
        """The str representation of the constraint"""
        if len(self.indices) > 1:
            indices_str = "[{:d}, ..., {:d}]".format(self.indices[0], self.indices[-1])
        else:
            indices_str = "[{:d}]".format(self.indices[0])

        if len(self.params) > 1:
            params_str = "[{:s}, ..., {:s}]".format(self.params[0], self.params[-1])
        elif len(self.params) == 1:
            params_str = "[{:s}]".format(self.params[0])
        else:
            params_str = "[]"

        return '{:s}({:s}, {:s}, ..., {:e})'.format(
            type(self).__name__,
            indices_str,
            params_str,
            self.eps
        )


class FixScaledParametricRelations(FixParametricRelations):

    def __init__(
        self,
        indices,
        Jacobian,
        const_shift,
        params=None,
        eps=1e-12,
    ):
        """The fractional coordinate version of FixParametricRelations

        All arguments are the same, but since this is for fractional coordinates use_cell is false
        """
        super(FixScaledParametricRelations, self).__init__(
            indices,
            Jacobian,
            const_shift,
            params,
            eps,
            False,
        )

    def adjust_contravariant(self, cell, vecs, B):
        """Adjust the values of a set of vectors that are contravariant with the unit transformation"""
        scaled = cell.scaled_positions(vecs).flatten()
        scaled = self.Jacobian_inv @ (scaled - B)
        scaled = ((self.Jacobian @ scaled) + B).reshape((-1, 3))

        return cell.cartesian_positions(scaled)

    def adjust_positions(self, atoms, positions):
        """Adjust positions of the atoms to match the constraints"""
        positions[self.indices] = self.adjust_contravariant(
            atoms.cell,
            positions[self.indices],
            self.const_shift,
        )
        positions[self.indices] = self.adjust_B(atoms.cell, positions[self.indices])

    def adjust_B(self, cell, positions):
        """Wraps the positions back to the unit cell and adjust B to keep track of this change"""
        fractional = cell.scaled_positions(positions)
        wrapped_fractional = (fractional % 1.0) % 1.0
        self.const_shift += np.round(wrapped_fractional - fractional).flatten()
        return cell.cartesian_positions(wrapped_fractional)

    def adjust_momenta(self, atoms, momenta):
        """Adjust momenta of the atoms to match the constraints"""
        momenta[self.indices] = self.adjust_contravariant(
            atoms.cell,
            momenta[self.indices],
            np.zeros(self.const_shift.shape),
        )

    def adjust_forces(self, atoms, forces):
        """Adjust forces of the atoms to match the constraints"""
        # Forces are coavarient to the coordinate transformation, use the inverse transformations
        cart2frac_jacob = np.zeros(2*(3*len(atoms),))
        for i_atom in range(len(atoms)):
            cart2frac_jacob[3*i_atom:3*(i_atom+1), 3*i_atom:3*(i_atom+1)] = atoms.cell.T

        jacobian = cart2frac_jacob @ self.Jacobian
        jacobian_inv = np.linalg.inv(jacobian.T @ jacobian) @ jacobian.T

        reduced_forces = jacobian.T @ forces.flatten()
        forces[self.indices] = (jacobian_inv.T @ reduced_forces).reshape(-1, 3)

    def todict(self):
        """Create a dictionary representation of the constraint"""
        dct = super(FixScaledParametricRelations, self).todict()
        del(dct["kwargs"]["use_cell"])
        return dct


class FixCartesianParametricRelations(FixParametricRelations):

    def __init__(
        self,
        indices,
        Jacobian,
        const_shift,
        params=None,
        eps=1e-12,
        use_cell=False,
    ):
        """The Cartesian coordinate version of FixParametricRelations"""
        super(FixCartesianParametricRelations, self).__init__(
            indices,
            Jacobian,
            const_shift,
            params,
            eps,
            use_cell,
        )

    def adjust_contravariant(self, vecs, B):
        """Adjust the values of a set of vectors that are contravariant with the unit transformation"""
        vecs = self.Jacobian_inv @ (vecs.flatten() - B)
        vecs = ((self.Jacobian @ vecs) + B).reshape((-1, 3))
        return vecs

    def adjust_positions(self, atoms, positions):
        """Adjust positions of the atoms to match the constraints"""
        if self.use_cell:
            return
        positions[self.indices] = self.adjust_contravariant(
            positions[self.indices],
            self.const_shift,
        )

    def adjust_momenta(self, atoms, momenta):
        """Adjust momenta of the atoms to match the constraints"""
        if self.use_cell:
            return
        momenta[self.indices] = self.adjust_contravariant(
            momenta[self.indices],
            np.zeros(self.const_shift.shape),
        )

    def adjust_forces(self, atoms, forces):
        """Adjust forces of the atoms to match the constraints"""
        if self.use_cell:
            return

        forces_reduced = self.Jacobian.T @ forces[self.indices].flatten()
        forces[self.indices] = (self.Jacobian_inv.T @ forces_reduced).reshape(-1, 3)

    def adjust_cell(self, atoms, cell):
        """Adjust the cell of the atoms to match the constraints"""
        if not self.use_cell:
            return
        cell[self.indices] = self.adjust_contravariant(
            cell[self.indices],
            np.zeros(self.const_shift.shape),
        )

    def adjust_stress(self, atoms, stress):
        """Adjust the stress of the atoms to match the constraints"""
        if not self.use_cell:
            return

        stress_3x3 = voigt_6_to_full_3x3_stress(stress)
        stress_reduced = self.Jacobian.T @ stress_3x3[self.indices].flatten()
        stress_3x3[self.indices] = (self.Jacobian_inv.T @ stress_reduced).reshape(-1, 3)

        stress[:] = full_3x3_to_voigt_6_stress(stress_3x3)


class Hookean(FixConstraint):
    """Applies a Hookean restorative force between a pair of atoms, an atom
    and a point, or an atom and a plane."""

    def __init__(self, a1, a2, k, rt=None):
        """Forces two atoms to stay close together by applying no force if
        they are below a threshold length, rt, and applying a Hookean
        restorative force when the distance between them exceeds rt. Can
        also be used to tether an atom to a fixed point in space or to a
        distance above a plane.

        a1 : int
           Index of atom 1
        a2 : one of three options
           1) index of atom 2
           2) a fixed point in cartesian space to which to tether a1
           3) a plane given as (A, B, C, D) in A x + B y + C z + D = 0.
        k : float
           Hooke's law (spring) constant to apply when distance
           exceeds threshold_length. Units of eV A^-2.
        rt : float
           The threshold length below which there is no force. The
           length is 1) between two atoms, 2) between atom and point.
           This argument is not supplied in case 3. Units of A.

        If a plane is specified, the Hooke's law force is applied if the atom
        is on the normal side of the plane. For instance, the plane with
        (A, B, C, D) = (0, 0, 1, -7) defines a plane in the xy plane with a z
        intercept of +7 and a normal vector pointing in the +z direction.
        If the atom has z > 7, then a downward force would be applied of
        k * (atom.z - 7). The same plane with the normal vector pointing in
        the -z direction would be given by (A, B, C, D) = (0, 0, -1, 7).

        References:

           Andrew A. Peterson,  Topics in Catalysis volume 57, pages40–53 (2014)
           https://link.springer.com/article/10.1007%2Fs11244-013-0161-8
        """

        if isinstance(a2, int):
            self._type = 'two atoms'
            self.indices = [a1, a2]
        elif len(a2) == 3:
            self._type = 'point'
            self.index = a1
            self.origin = np.array(a2)
        elif len(a2) == 4:
            self._type = 'plane'
            self.index = a1
            self.plane = a2
        else:
            raise RuntimeError('Unknown type for a2')
        self.threshold = rt
        self.spring = k

    def get_removed_dof(self, atoms):
        return 0

    def todict(self):
        dct = {'name': 'Hookean'}
        dct['kwargs'] = {'rt': self.threshold,
                         'k': self.spring}
        if self._type == 'two atoms':
            dct['kwargs']['a1'] = self.indices[0]
            dct['kwargs']['a2'] = self.indices[1]
        elif self._type == 'point':
            dct['kwargs']['a1'] = self.index
            dct['kwargs']['a2'] = self.origin
        elif self._type == 'plane':
            dct['kwargs']['a1'] = self.index
            dct['kwargs']['a2'] = self.plane
        else:
            raise NotImplementedError('Bad type: %s' % self._type)
        return dct

    def adjust_positions(self, atoms, newpositions):
        pass

    def adjust_momenta(self, atoms, momenta):
        pass

    def adjust_forces(self, atoms, forces):
        positions = atoms.positions
        if self._type == 'plane':
            A, B, C, D = self.plane
            x, y, z = positions[self.index]
            d = ((A * x + B * y + C * z + D) /
                 np.sqrt(A**2 + B**2 + C**2))
            if d < 0:
                return
            magnitude = self.spring * d
            direction = - np.array((A, B, C)) / np.linalg.norm((A, B, C))
            forces[self.index] += direction * magnitude
            return
        if self._type == 'two atoms':
            p1, p2 = positions[self.indices]
        elif self._type == 'point':
            p1 = positions[self.index]
            p2 = self.origin
        displace, _ = find_mic(p2 - p1, atoms.cell, atoms.pbc)
        bondlength = np.linalg.norm(displace)
        if bondlength > self.threshold:
            magnitude = self.spring * (bondlength - self.threshold)
            direction = displace / np.linalg.norm(displace)
            if self._type == 'two atoms':
                forces[self.indices[0]] += direction * magnitude
                forces[self.indices[1]] -= direction * magnitude
            else:
                forces[self.index] += direction * magnitude

    def adjust_potential_energy(self, atoms):
        """Returns the difference to the potential energy due to an active
        constraint. (That is, the quantity returned is to be added to the
        potential energy.)"""
        positions = atoms.positions
        if self._type == 'plane':
            A, B, C, D = self.plane
            x, y, z = positions[self.index]
            d = ((A * x + B * y + C * z + D) /
                 np.sqrt(A**2 + B**2 + C**2))
            if d > 0:
                return 0.5 * self.spring * d**2
            else:
                return 0.
        if self._type == 'two atoms':
            p1, p2 = positions[self.indices]
        elif self._type == 'point':
            p1 = positions[self.index]
            p2 = self.origin
        displace, _ = find_mic(p2 - p1, atoms.cell, atoms.pbc)
        bondlength = np.linalg.norm(displace)
        if bondlength > self.threshold:
            return 0.5 * self.spring * (bondlength - self.threshold)**2
        else:
            return 0.

    def get_indices(self):
        if self._type == 'two atoms':
            return self.indices
        elif self._type == 'point':
            return self.index
        elif self._type == 'plane':
            return self.index

    def index_shuffle(self, atoms, ind):
        # See docstring of superclass
        if self._type == 'two atoms':
            newa = [-1, -1]  # Signal error
            for new, old in slice2enlist(ind, len(atoms)):
                for i, a in enumerate(self.indices):
                    if old == a:
                        newa[i] = new
            if newa[0] == -1 or newa[1] == -1:
                raise IndexError('Constraint not part of slice')
            self.indices = newa
        elif (self._type == 'point') or (self._type == 'plane'):
            newa = -1   # Signal error
            for new, old in slice2enlist(ind, len(atoms)):
                if old == self.index:
                    newa = new
                    break
            if newa == -1:
                raise IndexError('Constraint not part of slice')
            self.index = newa

    def __repr__(self):
        if self._type == 'two atoms':
            return 'Hookean(%d, %d)' % tuple(self.indices)
        elif self._type == 'point':
            return 'Hookean(%d) to cartesian' % self.index
        else:
            return 'Hookean(%d) to plane' % self.index


class ExternalForce(FixConstraint):
    """Constraint object for pulling two atoms apart by an external force.

    You can combine this constraint for example with FixBondLength but make
    sure that *ExternalForce* comes first in the list if there are overlaps
    between atom1-2 and atom3-4:

    >>> con1 = ExternalForce(atom1, atom2, f_ext)
    >>> con2 = FixBondLength(atom3, atom4)
    >>> atoms.set_constraint([con1, con2])

    see ase/test/external_force.py"""

    def __init__(self, a1, a2, f_ext):
        self.indices = [a1, a2]
        self.external_force = f_ext

    def get_removed_dof(self, atoms):
        return 0

    def adjust_positions(self, atoms, new):
        pass

    def adjust_forces(self, atoms, forces):
        dist = np.subtract.reduce(atoms.positions[self.indices])
        force = self.external_force * dist / np.linalg.norm(dist)
        forces[self.indices] += (force, -force)

    def adjust_potential_energy(self, atoms):
        dist = np.subtract.reduce(atoms.positions[self.indices])
        return -np.linalg.norm(dist) * self.external_force

    def index_shuffle(self, atoms, ind):
        """Shuffle the indices of the two atoms in this constraint"""
        newa = [-1, -1]  # Signal error
        for new, old in slice2enlist(ind, len(atoms)):
            for i, a in enumerate(self.indices):
                if old == a:
                    newa[i] = new
        if newa[0] == -1 or newa[1] == -1:
            raise IndexError('Constraint not part of slice')
        self.indices = newa

    def __repr__(self):
        return 'ExternalForce(%d, %d, %f)' % (self.indices[0],
                                              self.indices[1],
                                              self.external_force)

    def todict(self):
        return {'name': 'ExternalForce',
                'kwargs': {'a1': self.indices[0], 'a2': self.indices[1],
                           'f_ext': self.external_force}}


class MirrorForce(FixConstraint):
    """Constraint object for mirroring the force between two atoms.

    This class is designed to find a transition state with the help of a
    single optimization. It can be used if the transition state belongs to a
    bond breaking reaction. First the given bond length will be fixed until
    all other degrees of freedom are optimized, then the forces of the two
    atoms will be mirrored to find the transition state. The mirror plane is
    perpendicular to the connecting line of the atoms. Transition states in
    dependence of the force can be obtained by stretching the molecule and
    fixing its total length with *FixBondLength* or by using *ExternalForce*
    during the optimization with *MirrorForce*.

    Parameters
    ----------
    a1: int
        First atom index.
    a2: int
        Second atom index.
    max_dist: float
        Upper limit of the bond length interval where the transition state
        can be found.
    min_dist: float
        Lower limit of the bond length interval where the transition state
        can be found.
    fmax: float
        Maximum force used for the optimization.

    Notes
    -----
    You can combine this constraint for example with FixBondLength but make
    sure that *MirrorForce* comes first in the list if there are overlaps
    between atom1-2 and atom3-4:

    >>> con1 = MirrorForce(atom1, atom2)
    >>> con2 = FixBondLength(atom3, atom4)
    >>> atoms.set_constraint([con1, con2])

    """

    def __init__(self, a1, a2, max_dist=2.5, min_dist=1., fmax=0.1):
        self.indices = [a1, a2]
        self.min_dist = min_dist
        self.max_dist = max_dist
        self.fmax = fmax

    def adjust_positions(self, atoms, new):
        pass

    def adjust_forces(self, atoms, forces):
        dist = np.subtract.reduce(atoms.positions[self.indices])
        d = np.linalg.norm(dist)
        if (d < self.min_dist) or (d > self.max_dist):
            # Stop structure optimization
            forces[:] *= 0
            return
        dist /= d
        df = np.subtract.reduce(forces[self.indices])
        f = df.dot(dist)
        con_saved = atoms.constraints
        try:
            con = [con for con in con_saved
                   if not isinstance(con, MirrorForce)]
            atoms.set_constraint(con)
            forces_copy = atoms.get_forces()
        finally:
            atoms.set_constraint(con_saved)
        df1 = -1 / 2. * f * dist
        forces_copy[self.indices] += (df1, -df1)
        # Check if forces would be converged if the bond with mirrored forces
        # would also be fixed
        if (forces_copy**2).sum(axis=1).max() < self.fmax**2:
            factor = 1.
        else:
            factor = 0.
        df1 = -(1 + factor) / 2. * f * dist
        forces[self.indices] += (df1, -df1)

    def index_shuffle(self, atoms, ind):
        """Shuffle the indices of the two atoms in this constraint

        """
        newa = [-1, -1]  # Signal error
        for new, old in slice2enlist(ind, len(atoms)):
            for i, a in enumerate(self.indices):
                if old == a:
                    newa[i] = new
        if newa[0] == -1 or newa[1] == -1:
            raise IndexError('Constraint not part of slice')
        self.indices = newa

    def __repr__(self):
        return 'MirrorForce(%d, %d, %f, %f, %f)' % (
            self.indices[0], self.indices[1], self.max_dist, self.min_dist,
            self.fmax)

    def todict(self):
        return {'name': 'MirrorForce',
                'kwargs': {'a1': self.indices[0], 'a2': self.indices[1],
                           'max_dist': self.max_dist,
                           'min_dist': self.min_dist, 'fmax': self.fmax}}


class MirrorTorque(FixConstraint):
    """Constraint object for mirroring the torque acting on a dihedral
    angle defined by four atoms.

    This class is designed to find a transition state with the help of a
    single optimization. It can be used if the transition state belongs to a
    cis-trans-isomerization with a change of dihedral angle. First the given
    dihedral angle will be fixed until all other degrees of freedom are
    optimized, then the torque acting on the dihedral angle will be mirrored
    to find the transition state. Transition states in
    dependence of the force can be obtained by stretching the molecule and
    fixing its total length with *FixBondLength* or by using *ExternalForce*
    during the optimization with *MirrorTorque*.

    This constraint can be used to find
    transition states of cis-trans-isomerization.

    a1    a4
    |      |
    a2 __ a3

    Parameters
    ----------
    a1: int
        First atom index.
    a2: int
        Second atom index.
    a3: int
        Third atom index.
    a4: int
        Fourth atom index.
    max_angle: float
        Upper limit of the dihedral angle interval where the transition state
        can be found.
    min_angle: float
        Lower limit of the dihedral angle interval where the transition state
        can be found.
    fmax: float
        Maximum force used for the optimization.

    Notes
    -----
    You can combine this constraint for example with FixBondLength but make
    sure that *MirrorTorque* comes first in the list if there are overlaps
    between atom1-4 and atom5-6:

    >>> con1 = MirrorTorque(atom1, atom2, atom3, atom4)
    >>> con2 = FixBondLength(atom5, atom6)
    >>> atoms.set_constraint([con1, con2])

    """

    def __init__(self, a1, a2, a3, a4, max_angle=2 * np.pi, min_angle=0.,
                 fmax=0.1):
        self.indices = [a1, a2, a3, a4]
        self.min_angle = min_angle
        self.max_angle = max_angle
        self.fmax = fmax

    def adjust_positions(self, atoms, new):
        pass

    def adjust_forces(self, atoms, forces):
        angle = atoms.get_dihedral(self.indices[0], self.indices[1],
                                   self.indices[2], self.indices[3])
        angle *= np.pi / 180.
        if (angle < self.min_angle) or (angle > self.max_angle):
            # Stop structure optimization
            forces[:] *= 0
            return
        p = atoms.positions[self.indices]
        f = forces[self.indices]

        f0 = (f[1] + f[2]) / 2.
        ff = f - f0
        p0 = (p[2] + p[1]) / 2.
        m0 = np.cross(p[1] - p0, ff[1]) / (p[1] - p0).dot(p[1] - p0)
        fff = ff - np.cross(m0, p - p0)
        d1 = np.cross(np.cross(p[1] - p0, p[0] - p[1]), p[1] - p0) / \
            (p[1] - p0).dot(p[1] - p0)
        d2 = np.cross(np.cross(p[2] - p0, p[3] - p[2]), p[2] - p0) / \
            (p[2] - p0).dot(p[2] - p0)
        omegap1 = (np.cross(d1, fff[0]) / d1.dot(d1)).dot(p[1] - p0) / \
            np.linalg.norm(p[1] - p0)
        omegap2 = (np.cross(d2, fff[3]) / d2.dot(d2)).dot(p[2] - p0) / \
            np.linalg.norm(p[2] - p0)
        omegap = omegap1 + omegap2
        con_saved = atoms.constraints
        try:
            con = [con for con in con_saved
                   if not isinstance(con, MirrorTorque)]
            atoms.set_constraint(con)
            forces_copy = atoms.get_forces()
        finally:
            atoms.set_constraint(con_saved)
        df1 = -1 / 2. * omegap * np.cross(p[1] - p0, d1) / \
            np.linalg.norm(p[1] - p0)
        df2 = -1 / 2. * omegap * np.cross(p[2] - p0, d2) / \
            np.linalg.norm(p[2] - p0)
        forces_copy[self.indices] += (df1, [0., 0., 0.], [0., 0., 0.], df2)
        # Check if forces would be converged if the dihedral angle with
        # mirrored torque would also be fixed
        if (forces_copy**2).sum(axis=1).max() < self.fmax**2:
            factor = 1.
        else:
            factor = 0.
        df1 = -(1 + factor) / 2. * omegap * np.cross(p[1] - p0, d1) / \
            np.linalg.norm(p[1] - p0)
        df2 = -(1 + factor) / 2. * omegap * np.cross(p[2] - p0, d2) / \
            np.linalg.norm(p[2] - p0)
        forces[self.indices] += (df1, [0., 0., 0.], [0., 0., 0.], df2)

    def index_shuffle(self, atoms, ind):
        # See docstring of superclass
        indices = []
        for new, old in slice2enlist(ind, len(atoms)):
            if old in self.indices:
                indices.append(new)
        if len(indices) == 0:
            raise IndexError('All indices in MirrorTorque not part of slice')
        self.indices = np.asarray(indices, int)

    def __repr__(self):
        return 'MirrorTorque(%d, %d, %d, %d, %f, %f, %f)' % (
            self.indices[0], self.indices[1], self.indices[2],
            self.indices[3], self.max_angle, self.min_angle, self.fmax)

    def todict(self):
        return {'name': 'MirrorTorque',
                'kwargs': {'a1': self.indices[0], 'a2': self.indices[1],
                           'a3': self.indices[2], 'a4': self.indices[3],
                           'max_angle': self.max_angle,
                           'min_angle': self.min_angle, 'fmax': self.fmax}}


class Filter:
    """Subset filter class."""

    def __init__(self, atoms, indices=None, mask=None):
        """Filter atoms.

        This filter can be used to hide degrees of freedom in an Atoms
        object.

        Parameters
        ----------
        indices : list of int
           Indices for those atoms that should remain visible.
        mask : list of bool
           One boolean per atom indicating if the atom should remain
           visible or not.

        If a Trajectory tries to save this object, it will instead
        save the underlying Atoms object.  To prevent this, override
        the iterimages method.
        """

        self.atoms = atoms
        self.constraints = []
        # Make self.info a reference to the underlying atoms' info dictionary.
        self.info = self.atoms.info

        if indices is None and mask is None:
            raise ValueError('Use "indices" or "mask".')
        if indices is not None and mask is not None:
            raise ValueError('Use only one of "indices" and "mask".')

        if mask is not None:
            self.index = np.asarray(mask, bool)
            self.n = self.index.sum()
        else:
            self.index = np.asarray(indices, int)
            self.n = len(self.index)

    def iterimages(self):
        # Present the real atoms object to Trajectory and friends
        return self.atoms.iterimages()

    def get_cell(self):
        """Returns the computational cell.

        The computational cell is the same as for the original system.
        """
        return self.atoms.get_cell()

    def get_pbc(self):
        """Returns the periodic boundary conditions.

        The boundary conditions are the same as for the original system.
        """
        return self.atoms.get_pbc()

    def get_positions(self):
        'Return the positions of the visible atoms.'
        return self.atoms.get_positions()[self.index]

    def set_positions(self, positions, **kwargs):
        'Set the positions of the visible atoms.'
        pos = self.atoms.get_positions()
        pos[self.index] = positions
        self.atoms.set_positions(pos, **kwargs)

    positions = property(get_positions, set_positions,
                         doc='Positions of the atoms')

    def get_momenta(self):
        'Return the momenta of the visible atoms.'
        return self.atoms.get_momenta()[self.index]

    def set_momenta(self, momenta, **kwargs):
        'Set the momenta of the visible atoms.'
        mom = self.atoms.get_momenta()
        mom[self.index] = momenta
        self.atoms.set_momenta(mom, **kwargs)

    def get_atomic_numbers(self):
        'Return the atomic numbers of the visible atoms.'
        return self.atoms.get_atomic_numbers()[self.index]

    def set_atomic_numbers(self, atomic_numbers):
        'Set the atomic numbers of the visible atoms.'
        z = self.atoms.get_atomic_numbers()
        z[self.index] = atomic_numbers
        self.atoms.set_atomic_numbers(z)

    def get_tags(self):
        'Return the tags of the visible atoms.'
        return self.atoms.get_tags()[self.index]

    def set_tags(self, tags):
        'Set the tags of the visible atoms.'
        tg = self.atoms.get_tags()
        tg[self.index] = tags
        self.atoms.set_tags(tg)

    def get_forces(self, *args, **kwargs):
        return self.atoms.get_forces(*args, **kwargs)[self.index]

    def get_stress(self, *args, **kwargs):
        return self.atoms.get_stress(*args, **kwargs)

    def get_stresses(self, *args, **kwargs):
        return self.atoms.get_stresses(*args, **kwargs)[self.index]

    def get_masses(self):
        return self.atoms.get_masses()[self.index]

    def get_potential_energy(self, **kwargs):
        """Calculate potential energy.

        Returns the potential energy of the full system.
        """
        return self.atoms.get_potential_energy(**kwargs)

    def get_chemical_symbols(self):
        return self.atoms.get_chemical_symbols()

    def get_initial_magnetic_moments(self):
        return self.atoms.get_initial_magnetic_moments()

    def get_calculator(self):
        """Returns the calculator.

        WARNING: The calculator is unaware of this filter, and sees a
        different number of atoms.
        """
        return self.atoms.calc

    @property
    def calc(self):
        return self.atoms.calc

    def get_celldisp(self):
        return self.atoms.get_celldisp()

    def has(self, name):
        'Check for existence of array.'
        return self.atoms.has(name)

    def __len__(self):
        'Return the number of movable atoms.'
        return self.n

    def __getitem__(self, i):
        'Return an atom.'
        return self.atoms[self.index[i]]


class StrainFilter(Filter):
    """Modify the supercell while keeping the scaled positions fixed.

    Presents the strain of the supercell as the generalized positions,
    and the global stress tensor (times the volume) as the generalized
    force.

    This filter can be used to relax the unit cell until the stress is
    zero.  If MDMin is used for this, the timestep (dt) to be used
    depends on the system size. 0.01/x where x is a typical dimension
    seems like a good choice.

    The stress and strain are presented as 6-vectors, the order of the
    components follow the standard engingeering practice: xx, yy, zz,
    yz, xz, xy.

    """

    def __init__(self, atoms, mask=None, include_ideal_gas=False):
        """Create a filter applying a homogeneous strain to a list of atoms.

        The first argument, atoms, is the atoms object.

        The optional second argument, mask, is a list of six booleans,
        indicating which of the six independent components of the
        strain that are allowed to become non-zero.  It defaults to
        [1,1,1,1,1,1].

        """

        self.strain = np.zeros(6)
        self.include_ideal_gas = include_ideal_gas

        if mask is None:
            mask = np.ones(6)
        else:
            mask = np.array(mask)

        Filter.__init__(self, atoms, mask=mask)
        self.mask = mask
        self.origcell = atoms.get_cell()

    def get_positions(self):
        return self.strain.reshape((2, 3)).copy()

    def set_positions(self, new):
        new = new.ravel() * self.mask
        eps = np.array([[1.0 + new[0], 0.5 * new[5], 0.5 * new[4]],
                        [0.5 * new[5], 1.0 + new[1], 0.5 * new[3]],
                        [0.5 * new[4], 0.5 * new[3], 1.0 + new[2]]])

        self.atoms.set_cell(np.dot(self.origcell, eps), scale_atoms=True)
        self.strain[:] = new

    def get_forces(self, **kwargs):
        stress = self.atoms.get_stress(include_ideal_gas=self.include_ideal_gas)
        return -self.atoms.get_volume() * (stress * self.mask).reshape((2, 3))

    def has(self, x):
        return self.atoms.has(x)

    def __len__(self):
        return 2


class UnitCellFilter(Filter):
    """Modify the supercell and the atom positions. """
    def __init__(self, atoms, mask=None,
                 cell_factor=None,
                 hydrostatic_strain=False,
                 constant_volume=False,
                 scalar_pressure=0.0):
        """Create a filter that returns the atomic forces and unit cell
        stresses together, so they can simultaneously be minimized.

        The first argument, atoms, is the atoms object. The optional second
        argument, mask, is a list of booleans, indicating which of the six
        independent components of the strain are relaxed.

        - True = relax to zero
        - False = fixed, ignore this component

        Degrees of freedom are the positions in the original undeformed cell,
        plus the deformation tensor (extra 3 "atoms"). This gives forces
        consistent with numerical derivatives of the potential energy
        with respect to the cell degreees of freedom.

        For full details see:
            E. B. Tadmor, G. S. Smith, N. Bernstein, and E. Kaxiras,
            Phys. Rev. B 59, 235 (1999)

        You can still use constraints on the atoms, e.g. FixAtoms, to control
        the relaxation of the atoms.

        >>> # this should be equivalent to the StrainFilter
        >>> atoms = Atoms(...)
        >>> atoms.set_constraint(FixAtoms(mask=[True for atom in atoms]))
        >>> ucf = UnitCellFilter(atoms)

        You should not attach this UnitCellFilter object to a
        trajectory. Instead, create a trajectory for the atoms, and
        attach it to an optimizer like this:

        >>> atoms = Atoms(...)
        >>> ucf = UnitCellFilter(atoms)
        >>> qn = QuasiNewton(ucf)
        >>> traj = Trajectory('TiO2.traj', 'w', atoms)
        >>> qn.attach(traj)
        >>> qn.run(fmax=0.05)

        Helpful conversion table:

        - 0.05 eV/A^3   = 8 GPA
        - 0.003 eV/A^3  = 0.48 GPa
        - 0.0006 eV/A^3 = 0.096 GPa
        - 0.0003 eV/A^3 = 0.048 GPa
        - 0.0001 eV/A^3 = 0.02 GPa

        Additional optional arguments:

        cell_factor: float (default float(len(atoms)))
            Factor by which deformation gradient is multiplied to put
            it on the same scale as the positions when assembling
            the combined position/cell vector. The stress contribution to
            the forces is scaled down by the same factor. This can be thought
            of as a very simple preconditioners. Default is number of atoms
            which gives approximately the correct scaling.

        hydrostatic_strain: bool (default False)
            Constrain the cell by only allowing hydrostatic deformation.
            The virial tensor is replaced by np.diag([np.trace(virial)]*3).

        constant_volume: bool (default False)
            Project out the diagonal elements of the virial tensor to allow
            relaxations at constant volume, e.g. for mapping out an
            energy-volume curve. Note: this only approximately conserves
            the volume and breaks energy/force consistency so can only be
            used with optimizers that do require do a line minimisation
            (e.g. FIRE).

        scalar_pressure: float (default 0.0)
            Applied pressure to use for enthalpy pV term. As above, this
            breaks energy/force consistency.
        """

        Filter.__init__(self, atoms, indices=range(len(atoms)))
        self.atoms = atoms
        self.orig_cell = atoms.get_cell()
        self.stress = None

        if mask is None:
            mask = np.ones(6)
        mask = np.asarray(mask)
        if mask.shape == (6,):
            self.mask = voigt_6_to_full_3x3_stress(mask)
        elif mask.shape == (3, 3):
            self.mask = mask
        else:
            raise ValueError('shape of mask should be (3,3) or (6,)')

        if cell_factor is None:
            cell_factor = float(len(atoms))
        self.hydrostatic_strain = hydrostatic_strain
        self.constant_volume = constant_volume
        self.scalar_pressure = scalar_pressure
        self.cell_factor = cell_factor
        self.copy = self.atoms.copy
        self.arrays = self.atoms.arrays

    def deform_grad(self):
        return np.linalg.solve(self.orig_cell, self.atoms.cell).T

    def get_positions(self):
        """
        this returns an array with shape (natoms + 3,3).

        the first natoms rows are the positions of the atoms, the last
        three rows are the deformation tensor associated with the unit cell,
        scaled by self.cell_factor.
        """

        cur_deform_grad = self.deform_grad()
        natoms = len(self.atoms)
        pos = np.zeros((natoms + 3, 3))
        # UnitCellFilter's positions are the self.atoms.positions but without
        # the applied deformation gradient
        pos[:natoms] = np.linalg.solve(cur_deform_grad,
                                       self.atoms.positions.T).T
        # UnitCellFilter's cell DOFs are the deformation gradient times a
        # scaling factor
        pos[natoms:] = self.cell_factor * cur_deform_grad
        return pos

    def set_positions(self, new, **kwargs):
        """
        new is an array with shape (natoms+3,3).

        the first natoms rows are the positions of the atoms, the last
        three rows are the deformation tensor used to change the cell shape.

        the new cell is first set from original cell transformed by the new
        deformation gradient, then the positions are set with respect to the
        current cell by transforming them with the same deformation gradient
        """

        natoms = len(self.atoms)
        new_atom_positions = new[:natoms]
        new_deform_grad = new[natoms:] / self.cell_factor
        # Set the new cell from the original cell and the new
        # deformation gradient.  Both current and final structures should
        # preserve symmetry, so if set_cell() calls FixSymmetry.adjust_cell(),
        # it should be OK
        self.atoms.set_cell(self.orig_cell @ new_deform_grad.T,
                            scale_atoms=True)
        # Set the positions from the ones passed in (which are without the
        # deformation gradient applied) and the new deformation gradient.
        # This should also preserve symmetry, so if set_positions() calls
        # FixSymmetyr.adjust_positions(), it should be OK
        self.atoms.set_positions(new_atom_positions @ new_deform_grad.T,
                                 **kwargs)

    def get_potential_energy(self, force_consistent=True):
        """
        returns potential energy including enthalpy PV term.
        """
        atoms_energy = self.atoms.get_potential_energy(
            force_consistent=force_consistent)
        return atoms_energy + self.scalar_pressure * self.atoms.get_volume()

    def get_forces(self, **kwargs):
        """
        returns an array with shape (natoms+3,3) of the atomic forces
        and unit cell stresses.

        the first natoms rows are the forces on the atoms, the last
        three rows are the forces on the unit cell, which are
        computed from the stress tensor.
        """

        stress = self.atoms.get_stress(**kwargs)
        atoms_forces = self.atoms.get_forces(**kwargs)

        volume = self.atoms.get_volume()
        virial = -volume * (voigt_6_to_full_3x3_stress(stress) +
                            np.diag([self.scalar_pressure] * 3))
        cur_deform_grad = self.deform_grad()
        atoms_forces = atoms_forces @ cur_deform_grad
        virial = np.linalg.solve(cur_deform_grad, virial.T).T

        if self.hydrostatic_strain:
            vtr = virial.trace()
            virial = np.diag([vtr / 3.0, vtr / 3.0, vtr / 3.0])

        # Zero out components corresponding to fixed lattice elements
        if (self.mask != 1.0).any():
            virial *= self.mask

        if self.constant_volume:
            vtr = virial.trace()
            np.fill_diagonal(virial, np.diag(virial) - vtr / 3.0)

        natoms = len(self.atoms)
        forces = np.zeros((natoms + 3, 3))
        forces[:natoms] = atoms_forces
        forces[natoms:] = virial / self.cell_factor

        self.stress = -full_3x3_to_voigt_6_stress(virial) / volume
        return forces

    def get_stress(self):
        raise PropertyNotImplementedError

    def has(self, x):
        return self.atoms.has(x)

    def __len__(self):
        return (len(self.atoms) + 3)


class ExpCellFilter(UnitCellFilter):
    """Modify the supercell and the atom positions."""
    def __init__(self, atoms, mask=None,
                 cell_factor=None,
                 hydrostatic_strain=False,
                 constant_volume=False,
                 scalar_pressure=0.0):
        r"""Create a filter that returns the atomic forces and unit cell
        stresses together, so they can simultaneously be minimized.

        The first argument, atoms, is the atoms object. The optional second
        argument, mask, is a list of booleans, indicating which of the six
        independent components of the strain are relaxed.

        - True = relax to zero
        - False = fixed, ignore this component

        Degrees of freedom are the positions in the original undeformed cell,
        plus the log of the deformation tensor (extra 3 "atoms"). This gives
        forces consistent with numerical derivatives of the potential energy
        with respect to the cell degrees of freedom.

        For full details see:
            E. B. Tadmor, G. S. Smith, N. Bernstein, and E. Kaxiras,
            Phys. Rev. B 59, 235 (1999)

        You can still use constraints on the atoms, e.g. FixAtoms, to control
        the relaxation of the atoms.

        >>> # this should be equivalent to the StrainFilter
        >>> atoms = Atoms(...)
        >>> atoms.set_constraint(FixAtoms(mask=[True for atom in atoms]))
        >>> ecf = ExpCellFilter(atoms)

        You should not attach this ExpCellFilter object to a
        trajectory. Instead, create a trajectory for the atoms, and
        attach it to an optimizer like this:

        >>> atoms = Atoms(...)
        >>> ecf = ExpCellFilter(atoms)
        >>> qn = QuasiNewton(ecf)
        >>> traj = Trajectory('TiO2.traj', 'w', atoms)
        >>> qn.attach(traj)
        >>> qn.run(fmax=0.05)

        Helpful conversion table:

        - 0.05 eV/A^3   = 8 GPA
        - 0.003 eV/A^3  = 0.48 GPa
        - 0.0006 eV/A^3 = 0.096 GPa
        - 0.0003 eV/A^3 = 0.048 GPa
        - 0.0001 eV/A^3 = 0.02 GPa

        Additional optional arguments:

        cell_factor: (DEPRECATED)
            Retained for backwards compatibility, but no longer used.

        hydrostatic_strain: bool (default False)
            Constrain the cell by only allowing hydrostatic deformation.
            The virial tensor is replaced by np.diag([np.trace(virial)]*3).

        constant_volume: bool (default False)
            Project out the diagonal elements of the virial tensor to allow
            relaxations at constant volume, e.g. for mapping out an
            energy-volume curve.

        scalar_pressure: float (default 0.0)
            Applied pressure to use for enthalpy pV term. As above, this
            breaks energy/force consistency.

        Implementation details:

        The implementation is based on that of Christoph Ortner in JuLIP.jl:
        https://github.com/libAtoms/JuLIP.jl/blob/expcell/src/Constraints.jl#L244

        We decompose the deformation gradient as

            F = exp(U) F0
            x =  F * F0^{-1} z  = exp(U) z

        If we write the energy as a function of U we can transform the
        stress associated with a perturbation V into a derivative using a
        linear map V -> L(U, V).

        \phi( exp(U+tV) (z+tv) ) ~ \phi'(x) . (exp(U) v) + \phi'(x) .
                                                    ( L(U, V) exp(-U) exp(U) z )

        where

               \nabla E(U) : V  =  [S exp(-U)'] : L(U,V)
                                =  L'(U, S exp(-U)') : V
                                =  L(U', S exp(-U)') : V
                                =  L(U, S exp(-U)) : V     (provided U = U')

        where the : operator represents double contraction,
        i.e. A:B = trace(A'B), and

          F = deformation tensor - 3x3 matrix
          F0 = reference deformation tensor - 3x3 matrix, np.eye(3) here
          U = cell degrees of freedom used here - 3x3 matrix
          V = perturbation to cell DoFs - 3x3 matrix
          v = perturbation to position DoFs
          x = atomic positions in deformed cell
          z = atomic positions in original cell
          \phi = potential energy
          S = stress tensor [3x3 matrix]
          L(U, V) = directional derivative of exp at U in direction V, i.e
          d/dt exp(U + t V)|_{t=0} = L(U, V)

        This means we can write

          d/dt E(U + t V)|_{t=0} = L(U, S exp (-U)) : V

        and therefore the contribution to the gradient of the energy is

          \nabla E(U) / \nabla U_ij =  [L(U, S exp(-U))]_ij
        """

        Filter.__init__(self, atoms, indices=range(len(atoms)))
        UnitCellFilter.__init__(self, atoms, mask, cell_factor,
                                hydrostatic_strain,
                                constant_volume, scalar_pressure)
        if cell_factor is not None:
            warn("cell_factor is deprecated")
        self.cell_factor = 1.0

    def get_positions(self):
        pos = UnitCellFilter.get_positions(self)
        natoms = len(self.atoms)
        pos[natoms:] = logm(self.deform_grad())
        return pos

    def set_positions(self, new, **kwargs):
        natoms = len(self.atoms)
        new2 = new.copy()
        new2[natoms:] = expm(new[natoms:])
        UnitCellFilter.set_positions(self, new2, **kwargs)

    def get_forces(self, **kwargs):
        forces = UnitCellFilter.get_forces(self, **kwargs)

        # forces on atoms are same as UnitCellFilter, we just
        # need to modify the stress contribution
        stress = self.atoms.get_stress(**kwargs)
        volume = self.atoms.get_volume()
        virial = -volume * (voigt_6_to_full_3x3_stress(stress) +
                            np.diag([self.scalar_pressure] * 3))

        cur_deform_grad = self.deform_grad()
        cur_deform_grad_log = logm(cur_deform_grad)

        if self.hydrostatic_strain:
            vtr = virial.trace()
            virial = np.diag([vtr / 3.0, vtr / 3.0, vtr / 3.0])

        # Zero out components corresponding to fixed lattice elements
        if (self.mask != 1.0).any():
            virial *= self.mask

        deform_grad_log_force_naive = virial.copy()
        Y = np.zeros((6, 6))
        Y[0:3, 0:3] = cur_deform_grad_log
        Y[3:6, 3:6] = cur_deform_grad_log
        Y[0:3, 3:6] = - virial @ expm(-cur_deform_grad_log)
        deform_grad_log_force = -expm(Y)[0:3, 3:6]
        for (i1, i2) in [(0, 1), (0, 2), (1, 2)]:
            ff = 0.5 * (deform_grad_log_force[i1, i2] +
                        deform_grad_log_force[i2, i1])
            deform_grad_log_force[i1, i2] = ff
            deform_grad_log_force[i2, i1] = ff

        # check for reasonable alignment between naive and
        # exact search directions
        all_are_equal = np.all(np.isclose(deform_grad_log_force,
                                          deform_grad_log_force_naive))
        if all_are_equal or \
            (np.sum(deform_grad_log_force * deform_grad_log_force_naive) /
             np.sqrt(np.sum(deform_grad_log_force**2) *
                     np.sum(deform_grad_log_force_naive**2)) > 0.8):
            deform_grad_log_force = deform_grad_log_force_naive

        # Cauchy stress used for convergence testing
        convergence_crit_stress = -(virial / volume)
        if self.constant_volume:
            # apply constraint to force
            dglf_trace = deform_grad_log_force.trace()
            np.fill_diagonal(deform_grad_log_force,
                             np.diag(deform_grad_log_force) - dglf_trace / 3.0)
            # apply constraint to Cauchy stress used for convergence testing
            ccs_trace = convergence_crit_stress.trace()
            np.fill_diagonal(convergence_crit_stress,
                             np.diag(convergence_crit_stress) - ccs_trace / 3.0)

        # pack gradients into vector
        natoms = len(self.atoms)
        forces[natoms:] = deform_grad_log_force
        self.stress = full_3x3_to_voigt_6_stress(convergence_crit_stress)
        return forces
