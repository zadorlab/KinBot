"""Versioned rigid-rotor counting in KinBot's cached single-state templates.

The geometry correction constrains external rotations of a whole saddle or
disconnected assembly. It does not infer missing optical or rotor coverage.
"""
import re
import numpy as np

MARKER = '! kinbot_whole_geometry_sigma_v1'


def mark_current(text):
    """New writers have already applied the whole-geometry constraint."""
    return text if MARKER in text else MARKER + '\n' + text


def upgrade_cached_rigid_rotor(text, charge, multiplicity, *, saddle=False):
    """Update a legacy single-RRHO worker fragment in memory, at most once.

    Retain its existing optical division. Migration is limited to the legacy
    one/two optical convention and a reconstructible graph sigma. Multiple
    RRHO members, HIR models and custom divisors require their original
    member/domain evidence and are not inferred from a serialized geometry.
    """
    if MARKER in text or len(re.findall(r'^\s*RRHO\s*$', text, re.M)) != 1:
        return text
    if re.search(r'Rotor\s+Hindered|Core\s+(?:PhaseSpaceTheory|Rotd)\b', text):
        return text
    geom = re.search(r'Geometry\[angstrom\]\s+(\d+)\s*\n', text)
    factor = re.search(r'(SymmetryFactor\s+)([-+\d.eE]+)', text)
    if geom is None or factor is None:
        return text
    rows = [row.split() for row in text[geom.end():].splitlines()[:int(geom[1])]]
    if any(len(row) != 4 for row in rows):
        return text
    from kinbot.stationary_pt import StationaryPoint
    from kinbot import symmetry
    point = StationaryPoint('cached_rrho', charge, multiplicity,
        atom=[r[0] for r in rows], geom=np.array([r[1:] for r in rows], float),
        wellorts=int(saddle))
    point.characterize()
    symmetry.calculate_symmetry(point)
    proof = getattr(point, 'rigid_rotation_evidence', {})
    old = proof.get('previous_graph_sigma')
    if old is None or point.sigma_ext == old:
        return text
    value = float(factor[2])
    if value <= 0 or not any(abs(old/value-n) < 1.e-8 for n in (1, 2)):
        return text
    corrected = value * point.sigma_ext / old
    text = text[:factor.start(2)] + f'{corrected:.12g}' + text[factor.end(2):]
    return (MARKER + '\n! External rotational sigma corrected from '
        f'{old:g} to {point.sigma_ext:g} using proper rotations of the whole '
        'serialized geometry; existing optical division retained.\n' + text)
