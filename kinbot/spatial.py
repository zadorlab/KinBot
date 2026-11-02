"""Common Cartesian transformations; callers decide chemical equivalence.

Atom mapping, allowed reflections, conformational motion and numerical
tolerances belong to the scientific caller. A vector is transformed with its
geometry, never independently fitted or translated to improve its overlap.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Alignment:
    rotation: np.ndarray
    translation: np.ndarray
    errors: np.ndarray

    @property
    def rmsd(self):
        return float(np.sqrt(np.mean(self.errors**2)))


def proper_alignment(source, target):
    """Map fixed-order row coordinates onto target with determinant +1."""
    source, target = np.asarray(source, float), np.asarray(target, float)
    if (source.ndim != 2 or source.shape[1:] != (3,) or not len(source)
            or source.shape != target.shape
            or not np.isfinite(source).all() or not np.isfinite(target).all()):
        raise ValueError('Alignment requires finite, matching N by 3 coordinates.')
    x, y = source-source.mean(0), target-target.mean(0)
    u, _, vt = np.linalg.svd(x.T @ y)
    sign = np.eye(3)
    sign[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ sign @ vt
    translation = target.mean(0)-source.mean(0) @ rotation
    errors = np.linalg.norm(x @ rotation-y, axis=1)
    return Alignment(rotation, translation, errors)


def mode_overlap(left, right, *, masses=None, signed=False):
    """Normalized displacement overlap in an already common atomic frame.

    Masses give the vibrational inner product for Cartesian displacements.
    Omit masses only for vectors already mass-weighted. For comparisons of
    two independent eigenvectors their arbitrary signs are ignored. Mirror
    parity uses signed=True: the mode and its own transform share one sign.
    """
    a, b = np.asarray(left, float), np.asarray(right, float)
    if (a.shape != b.shape or a.ndim != 2 or a.shape[1:] != (3,)
            or not np.isfinite(a).all() or not np.isfinite(b).all()):
        raise ValueError('Mode overlap requires finite, matching N by 3 vectors.')
    if masses is not None:
        m = np.asarray(masses, float)
        if m.shape != (len(a),) or not np.isfinite(m).all() or np.any(m <= 0):
            raise ValueError('Mode overlap requires one positive mass per atom.')
        a, b = a*np.sqrt(m[:, None]), b*np.sqrt(m[:, None])
    norm = np.linalg.norm(a)*np.linalg.norm(b)
    if not norm:
        raise ValueError('A zero displacement cannot establish mode overlap.')
    overlap = float(np.clip(np.sum(a*b)/norm, -1., 1.))
    return overlap if signed else abs(overlap)


def aligned_mode_overlap(left_geometry, right_geometry, left_mode, right_mode, *, masses):
    """Compare fixed-order modes after aligning right geometry to left."""
    fit = proper_alignment(right_geometry, left_geometry)
    return mode_overlap(left_mode, np.asarray(right_mode) @ fit.rotation, masses=masses)
