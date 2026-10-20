"""Read native orbital cubes and obtain local lobe directions for VRC pivots."""

from pathlib import Path

import numpy as np
import re
from ase import Atoms
from ase.data import chemical_symbols
from ase.units import Bohr


def _number(value):
    value = value.replace('D', 'E').replace('d', 'e')
    # Fortran may omit E when a three-digit exponent fills the field.
    if re.fullmatch(r'[+-]?\d*\.\d+[+-]\d{3}', value):
        value = re.sub(r'([+-]\d{3})$', r'E\1', value)
    return float(value)


def _fields(line, count):
    fields = line.split()
    if len(fields) == count:
        return fields
    # Native Molpro headers use I5 followed by F12.6 fields. Negative values
    # can fill adjacent fields completely, without a separating blank.
    if len(line.rstrip('\n')) == 5 + 12 * (count - 1):
        return [line[:5]] + [line[i:i+12] for i in range(5, 5+12*(count-1), 12)]
    return fields


def read_orbital_cube(path):
    """Strict single-orbital Gaussian/Molpro formatted cube reader.

    Molpro documents Bohr coordinates and one scalar file per orbital.
    Gaussian's negative atom count introduces a counted orbital-label list.
    Check that list at EOF rather than looping indefinitely on a short file.
    """
    path = Path(path)
    try:
        with path.open() as stream:
            title, description = stream.readline(), stream.readline()
            header = _fields(stream.readline(), 4)
            if len(header) not in (4, 5):
                raise ValueError('invalid origin header')
            signed_atoms = int(header[0])
            nval = int(header[4]) if len(header) == 5 else 1
            if signed_atoms == 0 or nval != 1:
                raise ValueError('exactly one orbital and a nonempty molecule required')
            number = _number
            origin = np.array([number(v) for v in header[1:4]])
            shape, spacing = [], []
            for _ in range(3):
                row = _fields(stream.readline(), 4)
                if len(row) != 4 or int(row[0]) < 2:
                    raise ValueError('invalid grid; Molpro/Gaussian Bohr grid required')
                shape.append(int(row[0]))
                spacing.append([number(v) for v in row[1:]])
            numbers, positions = [], []
            for _ in range(abs(signed_atoms)):
                row = _fields(stream.readline(), 5)
                if len(row) != 5:
                    raise ValueError('incomplete atom row')
                numbers.append(int(row[0]))
                positions.append([number(v) for v in row[2:]])
            labels = []
            if signed_atoms < 0:
                row = stream.readline().split()
                if len(row) < 1 or int(row[0]) != 1:
                    raise ValueError('exactly one labeled orbital required')
                labels = row[1:]
                if not labels:
                    labels = stream.readline().split()
                if len(labels) != 1:
                    raise ValueError('incomplete or ambiguous orbital labels')
                labels = [int(labels[0])]
            values = []
            for line in stream:
                try:
                    row = [number(v) for v in line.split()]
                except ValueError:
                    raw = line.rstrip('\n')
                    if len(raw) % 13:
                        raise ValueError('invalid scalar field width')
                    row = [number(raw[i:i+13]) for i in range(0, len(raw), 13)]
                values.extend(row)
            values = np.asarray(values)
        if values.size != int(np.prod(shape)):
            raise ValueError('incomplete or excess scalar data')
        axes = [0, 1, 2]
        if 'OUTER LOOP' in description.upper():
            names = re.findall(r'(?:OUTER|MIDDLE|INNER) LOOP:\s*([XYZ])', description.upper())
            if len(names) != 3 or len(set(names)) != 3:
                raise ValueError('invalid loop order')
            axes = ['XYZ'.index(name) for name in names]
        data = values.reshape([shape[a] for a in axes]).transpose(np.argsort(axes))
        origin *= Bohr
        spacing = np.asarray(spacing) * Bohr
        atoms = Atoms(numbers=numbers, positions=np.asarray(positions) * Bohr)
        arrays = (data, origin, spacing, atoms.positions)
        if (any(not np.isfinite(a).all() for a in arrays)
                or any(not 1 <= n < len(chemical_symbols) for n in numbers)
                or abs(np.linalg.det(spacing)) < 1.e-15):
            raise ValueError('invalid orbital cube geometry or grid')
        return {'data': data, 'datas': data[None], 'origin': origin,
                'spacing': spacing, 'atoms': atoms, 'labels': labels,
                'title': title.strip(), 'description': description.strip()}
    except (ValueError, IndexError, OverflowError) as error:
        raise ValueError(f'{path}: malformed or incomplete orbital cube: {error}') from error


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
