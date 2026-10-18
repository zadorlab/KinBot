"""Read native orbital cubes and obtain local lobe directions for VRC pivots."""

from pathlib import Path

import numpy as np
from ase.io.cube import read_cube


def read_orbital_cube(path):
    """Read one scalar orbital from either Gaussian's or Molpro's cube layout.

    ASE handles the optional orbital-label row, optional fifth header field,
    arbitrary value wrapping, and Bohr-to-Angstrom conversion.  A density or
    orbital's identity must still be established by its producing calculation.
    """
    path = Path(path)
    try:
        with path.open() as stream:
            cube = read_cube(stream)
    except (ValueError, IndexError) as error:
        raise ValueError(f'{path}: malformed or incomplete orbital cube: {error}') from error
    if len(cube['datas']) != 1 or len(cube['labels']) > 1:
        raise ValueError(f'{path}: VRC pivots require exactly one orbital.')
    arrays = (cube['data'], cube['origin'], cube['spacing'], cube['atoms'].positions)
    if (not len(cube['atoms']) or any(not np.isfinite(a).all() for a in arrays)
            or any(n < 2 for n in cube['data'].shape)
            or abs(np.linalg.det(cube['spacing'])) < 1.e-15):
        raise ValueError(f'{path}: invalid orbital cube geometry or grid.')
    return cube


def lobe_directions(cube, center, *, radius=0.8, relative_lobe=0.8):
    """Find significant opposite-sign extrema within ``radius`` Angstrom.

    Use the full voxel matrix, including rotated/skew grids. The stronger
    lobe is always retained; the other is retained at the historic 80% cutoff.
    Choosing by absolute amplitude makes the result invariant to orbital sign.
    """
    center = np.asarray(center, dtype=float)
    spacing = cube['spacing']
    inverse = np.linalg.inv(spacing)
    fractional = (center - cube['origin']) @ inverse
    shape = np.asarray(cube['data'].shape)
    if np.any(fractional < 0) or np.any(fractional > shape - 1):
        raise ValueError('Reactive atom lies outside the orbital cube grid.')
    extent = radius * np.linalg.norm(inverse, axis=0)
    low = np.maximum(0, np.floor(fractional - extent).astype(int))
    high = np.minimum(shape, np.ceil(fractional + extent).astype(int) + 1)
    indices = np.indices(high - low).reshape(3, -1).T + low
    vectors = cube['origin'] + indices @ spacing - center
    distances = np.linalg.norm(vectors, axis=1)
    mask = (distances < radius) & (distances > 1.e-8)
    if not mask.any():
        raise ValueError('Orbital cube has no grid points near the reactive atom.')
    indices, vectors, distances = indices[mask], vectors[mask], distances[mask]
    values = cube['data'][tuple(indices.T)]
    peaks = []
    for sign in (1, -1):
        strength = sign * values
        candidates = np.flatnonzero(strength == strength.max())
        best = candidates[np.argmin(distances[candidates])]
        if strength[best] > 0:
            peaks.append((strength[best], vectors[best] / distances[best]))
    if not peaks:
        raise ValueError('Orbital has zero amplitude near the reactive atom.')
    peaks.sort(key=lambda item: -item[0])
    directions = [peaks[0][1]]
    if (len(peaks) == 2 and peaks[1][0] > relative_lobe * peaks[0][0]
            and np.linalg.norm(peaks[1][1] - peaks[0][1]) > 0.2):
        directions.append(peaks[1][1])
    return np.asarray(sorted(directions, key=lambda v: tuple(v)), dtype=np.float32)
