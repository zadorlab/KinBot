"""Gaussian formatted-checkpoint Hessians in the selected geometry frame."""
from pathlib import Path
import re

import numpy as np
from ase.data import atomic_numbers
from ase.units import Bohr


def _array(lines, label):
    headers = [i for i, line in enumerate(lines) if line.startswith(label)]
    if len(headers) != 1:
        raise ValueError(f'Gaussian checkpoint needs one {label} array.')
    index = headers[0]
    match = re.search(r'\bN=\s*(\d+)\s*$', lines[index])
    if not match:
        raise ValueError(f'Gaussian checkpoint has invalid {label} header.')
    count, values = int(match[1]), []
    for line in lines[index + 1:]:
        if len(values) >= count:
            break
        values.extend(float(v.replace('D', 'E')) for v in line.split())
    if len(values) != count or not np.isfinite(values).all():
        raise ValueError(f'Gaussian checkpoint has incomplete/nonfinite {label}.')
    return np.asarray(values)


def read_formatted_hessian(path, symbols, positions):
    """Verify atom order/geometry and rotate every Cartesian Hessian block.

    Gaussian may use a standard orientation even for a fixed-geometry ASE
    job. The checkpoint's coordinates and force constants share that frame;
    the database geometry need not. No atom permutation or reflection is used.
    """
    lines = Path(path).read_text().splitlines()
    numbers = _array(lines, 'Atomic numbers')
    if not np.array_equal(numbers, [atomic_numbers[str(s)] for s in symbols]):
        raise ValueError('Gaussian checkpoint atom order differs from selected geometry.')
    n = len(symbols)
    native = _array(lines, 'Current cartesian coordinates').reshape(n, 3) * Bohr
    target = np.asarray(positions, dtype=float)
    if target.shape != native.shape or not np.isfinite(target).all():
        raise ValueError('Invalid selected Hessian geometry.')
    x, y = native - native.mean(0), target - target.mean(0)
    u, _, vt = np.linalg.svd(x.T @ y)
    if np.linalg.det(u @ vt) < 0:
        u[:, -1] *= -1
    rotation = u @ vt
    if np.max(np.linalg.norm(x @ rotation - y, axis=1)) > 1.e-5:
        raise ValueError('Gaussian checkpoint geometry differs from selected calculation.')
    values = _array(lines, 'Cartesian Force Constants')
    hessian = np.zeros((3*n, 3*n))
    if len(values) != 3*n*(3*n+1)//2:
        raise ValueError('Gaussian checkpoint force-constant count differs.')
    hessian[np.tril_indices(3*n)] = values
    hessian += np.tril(hessian, -1).T
    transform = np.kron(np.eye(n), rotation)
    return transform.T @ hessian @ transform
