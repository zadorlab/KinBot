"""Local motion diagnostics; never a pathway or state-counting decision.

Compare an observed Cartesian displacement with physical internal-rotation
tangents in the mass metric after removing external translations/rotations.
The residual includes bends, stretches, and inversion; it is not automatically
classified as an in-plane bend. A tangent projection says nothing about which
finite path has the lowest barrier or whether a partition function covers it.
"""
from numbers import Integral

import networkx as nx
import numpy as np

from kinbot.constants import exact_mass


def _basis(rows, width):
    rows = np.asarray(rows, float).reshape(-1, width)
    if not len(rows):
        return np.empty((0, width))
    norms = np.linalg.norm(rows, axis=1)
    rows = rows[norms > 1.e-12] / norms[norms > 1.e-12, None]
    if not len(rows):
        return np.empty((0, width))
    _, singular, vectors = np.linalg.svd(rows, full_matrices=False)
    return vectors[singular > 1.e-10 * singular[0]]


def tangent_subspaces(geometry, masses, bond, axes):
    """Return external and torsional orthonormal rows and individual tangents.

    Axes are zero based, must be actual graph bridges, and do not themselves
    assert that any scan or chemical interconversion exists. Rotating either
    side gives the same tangent after removing overall molecular rotation.
    """
    xyz = np.asarray(geometry, float)
    masses = np.asarray(masses, float)
    n = len(xyz)
    bond = np.asarray(bond)
    if (n == 0 or xyz.shape != (n, 3) or masses.shape != (n,)
            or bond.shape != (n, n) or not np.isfinite(xyz).all()
            or not np.isfinite(masses).all() or np.any(masses <= 0.)
            or not np.isfinite(bond).all() or not np.array_equal(bond, bond.T)):
        raise ValueError('Finite geometry, positive masses, and a symmetric graph are required.')
    xyz = xyz - np.average(xyz, axis=0, weights=masses)
    sqrtmass = np.sqrt(masses)[:, None]
    external = [(np.tile(v, (n, 1))*sqrtmass).ravel() for v in np.eye(3)]
    external += [(np.cross(xyz, v)*sqrtmass).ravel() for v in np.eye(3)]
    external = _basis(external, 3*n)
    graph = nx.from_numpy_array(bond > 0)
    bridges = {tuple(sorted(edge)) for edge in nx.bridges(graph)}
    tangents, records, seen = [], [], set()
    for pair in axes:
        if (len(pair) != 2 or any(not isinstance(i, Integral) or isinstance(i, (bool, np.bool_))
                                  for i in pair)):
            raise ValueError('Axes must contain two zero-based integer atom indices.')
        axis = tuple(sorted(map(int, pair)))
        if axis in seen:
            raise ValueError('Duplicate torsional axis.')
        seen.add(axis)
        if axis not in bridges:
            raise ValueError('An internal-rotation axis must be a bond that splits the graph.')
        cut = graph.copy()
        cut.remove_edge(*axis)
        top = sorted(nx.node_connected_component(cut, axis[1]))
        unit = xyz[axis[1]]-xyz[axis[0]]
        norm = np.linalg.norm(unit)
        if norm < 1.e-10:
            raise ValueError('Zero-length torsional axis.')
        unit /= norm
        vector = np.zeros_like(xyz)
        vector[top] = np.cross(xyz[top]-xyz[axis[0]], unit)
        vector = (vector*sqrtmass).ravel()
        vector -= external.T @ (external @ vector)
        norm = np.linalg.norm(vector)
        # E.g. rotating a terminal atom on the axis has no internal motion.
        tangents.append(vector/norm if norm > 1.e-10 else np.zeros(3*n))
        records.append(dict(axis=list(axis), top=top, nonzero=bool(norm > 1.e-10)))
    return external, _basis(tangents, 3*n), tangents, records


def characterize_modes(species, cartesian_modes, *, masses=None):
    """Report squared local projections for supplied Cartesian displacements.

    Gaussian printed modes can be used for character diagnostics, not Hessian
    reconstruction or quantitative replacement of a vibrational spectrum.
    Isotope labels require explicit corresponding masses.
    """
    if masses is None:
        if any(getattr(species, 'isotopes', ())):
            raise ValueError('Supply isotope-specific masses for motion diagnostics.')
        masses = [exact_mass[str(atom)] for atom in species.atom]
    modes = np.asarray(cartesian_modes, float)
    if modes.ndim != 3 or modes.shape[1:] != (len(species.atom), 3) or not np.isfinite(modes).all():
        raise ValueError('Finite Cartesian modes with shape (modes, atoms, 3) are required.')
    axes = [list(map(int, rotor[1:3])) for rotor in species.dihed]
    external, torsions, individual, axes = tangent_subspaces(
        species.geom, masses, species.bond, axes)
    mw = np.sqrt(np.repeat(masses, 3))
    records = []
    for mode in modes:
        original = mode.ravel()*mw
        norm = np.linalg.norm(original)
        vector = original - external.T @ (external @ original)
        internal_norm = np.linalg.norm(vector)
        if norm == 0. or internal_norm < 1.e-8*norm:
            raise ValueError('A supplied displacement has no resolvable internal motion.')
        vector /= internal_norm
        fraction = float(np.clip(np.sum((torsions @ vector)**2), 0., 1.))
        records.append(dict(torsional_fraction=fraction,
            residual_fraction=1.-fraction,
            individual_axis_squared_overlaps=[float(np.clip(np.dot(v, vector)**2, 0., 1.))
                                              for v in individual],
            external_fraction=float(np.clip(1.-(internal_norm/norm)**2, 0., 1.))))
    return dict(schema=1, method='local_mass_weighted_tangent_projection',
        definition='squared norm of projection onto an orthonormal torsional subspace',
        individual_overlaps_are_additive=False,
        external_rank=len(external), torsional_rank=len(torsions),
        masses_amu=list(map(float, masses)), axes=axes, modes=records,
        lowest_path_established=False, alternative_paths_excluded=False,
        anharmonic_partition_function_established=False,
        population_merge_authorized=False)
