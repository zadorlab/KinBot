"""Audit explicitly declared geometry symmetries for orbital reconnaissance.

Operations and atom permutations are input data. This module does not infer a
reaction path or prescribe a molecular point group. Electronic symmetry and
orbital/state continuity require separate checks.
"""
import numpy as np


def operation_arrays(symbols, operation):
    """Validate a Cartesian column-vector operation and its atom mapping."""
    matrix = np.asarray(operation['matrix'], dtype=float)
    origin = np.asarray(operation['origin_angstrom'], dtype=float)
    permutation = operation['permutation_0based']
    if (matrix.shape != (3, 3) or origin.shape != (3,)
            or not np.isfinite(matrix).all() or not np.isfinite(origin).all()
            or not np.allclose(matrix.T @ matrix, np.eye(3), atol=1.e-10, rtol=0)
            or any(type(i) is not int for i in permutation)
            or sorted(permutation) != list(range(len(symbols)))
            or any(symbols[i] != symbols[j] for i, j in enumerate(permutation))):
        raise ValueError('Invalid symmetry operation or element-preserving atom permutation')
    return matrix, origin, permutation


def audit_symmetry(symbols, coordinates, operations, tolerance_angstrom):
    xyz = np.asarray(coordinates, dtype=float)
    if (xyz.shape != (len(symbols), 3) or not np.isfinite(xyz).all()
            or not np.isfinite(tolerance_angstrom) or tolerance_angstrom <= 0
            or not operations):
        raise ValueError('Finite coordinates, operations and a positive tolerance are required')
    rows = []
    for operation in operations:
        matrix, origin, permutation = operation_arrays(symbols, operation)
        mapped = (xyz - origin) @ matrix.T + origin
        residual = float(np.linalg.norm(mapped - xyz[permutation], axis=1).max())
        rows.append({'name': operation['name'], 'maximum_residual_angstrom': residual})
    maximum = max(row['maximum_residual_angstrom'] for row in rows)
    return {'passed': maximum <= tolerance_angstrom,
            'maximum_residual_angstrom': maximum,
            'tolerance_angstrom': tolerance_angstrom, 'operations': rows,
            'electronic_symmetry_verified': False, 'orbital_continuity_verified': False}


def project_symmetry(symbols, coordinates, operations, *, maximum_displacement_angstrom,
                     tolerance_angstrom=1.e-10):
    """Project onto declared constraints, rejecting large or inconsistent changes.

The caller must save both geometries and the displacement. Projection is a
geometry construction, not an energy optimization. No existing evidence is
modified. Iterative averaging also accepts a set of generators rather than a
complete group, provided the requested constraints have a common solution.
"""
    original = np.asarray(coordinates, dtype=float)
    audit_symmetry(symbols, original, operations, tolerance_angstrom)
    if not np.isfinite(maximum_displacement_angstrom) or maximum_displacement_angstrom <= 0:
        raise ValueError('A positive maximum symmetry projection displacement is required')
    xyz = original.copy()
    for _ in range(1000):
        transformed = [xyz]
        for operation in operations:
            matrix, origin, permutation = operation_arrays(symbols, operation)
            transformed.append(((xyz-origin) @ matrix.T + origin)[np.argsort(permutation)])
        xyz = np.mean(transformed, axis=0)
        displacement = float(np.linalg.norm(xyz-original, axis=1).max())
        if displacement > maximum_displacement_angstrom:
            raise ValueError('Requested symmetry projection exceeds the allowed displacement')
        report = audit_symmetry(symbols, xyz, operations, tolerance_angstrom)
        if report['passed']:
            return xyz, {**report, 'maximum_displacement_angstrom': displacement,
                         'geometry_optimization': False}
    raise ValueError('Declared symmetry constraints did not converge to a common geometry')
