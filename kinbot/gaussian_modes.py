"""Read completed native Gaussian normal modes in a selected geometry frame."""
from pathlib import Path
import hashlib
import re
import numpy as np
from ase.data import atomic_numbers


def read_normal_modes(path, symbols, positions):
    """Require the same ordered atoms and geometry; retain native provenance.

    Printed Gaussian displacement vectors have limited precision. They are
    suitable for a well-separated symmetry-parity test, not a new Hessian or
    a frequency refinement. No QC or checkpoint conversion is performed.
    """
    path = Path(path)
    data = path.read_bytes()
    lines = data.decode(errors='replace').splitlines()
    geometry = mode_geometry = None
    frequencies, modes = [], []
    completed = False
    for i, line in enumerate(lines):
        if 'Harmonic frequencies (cm**-1)' in line:
            frequencies, modes = [], []
            completed = False
            mode_geometry = None
        if 'orientation:' in line and ('Standard' in line or 'Input' in line):
            rows = []
            for row in lines[i+5:]:
                words = row.split()
                if len(words) != 6 or not words[0].isdigit():
                    break
                rows.append(words)
            if rows:
                geometry = (np.array([int(row[1]) for row in rows]),
                            np.array([row[3:6] for row in rows], float))
        if 'Frequencies --' in line:
            if geometry is None:
                raise ValueError('Normal modes have no preceding native orientation.')
            values = [float(v) for v in line.split('--', 1)[1].split()]
            start = next((j for j in range(i+1, min(i+10, len(lines)))
                          if re.match(r'\s*Atom\s+AN\s+X', lines[j])), None)
            if start is None:
                raise ValueError('Native normal-mode vectors are missing.')
            vectors = []
            for row in lines[start+1:start+1+len(symbols)]:
                words = row.split()
                if len(words) != 2+3*len(values):
                    raise ValueError('Incomplete native displacement table.')
                atom = len(vectors)
                if (int(words[0]) != atom+1
                        or int(words[1]) != atomic_numbers[str(symbols[atom])]):
                    raise ValueError('Native displacement atom order differs from the selected calculation.')
                vectors.append([float(v) for v in words[2:]])
            if len(vectors) != len(symbols):
                raise ValueError('Incomplete native displacement table.')
            if mode_geometry is not None and (
                    not np.array_equal(mode_geometry[0], geometry[0])
                    or not np.array_equal(mode_geometry[1], geometry[1])):
                raise ValueError('Native displacement tables mix different geometries.')
            frequencies.extend(values)
            mode_geometry = geometry
            modes.extend(np.array(vectors).reshape(len(symbols), len(values), 3).transpose(1, 0, 2))
            completed = False
        if 'Normal termination of Gaussian' in line:
            completed = True
        if 'Error termination' in line:
            completed = False
    if not completed or not modes:
        raise ValueError('No normally completed native frequency calculation.')
    if len(modes) not in (3*len(symbols)-6, 3*len(symbols)-5):
        raise ValueError('Incomplete native normal-mode spectrum.')
    if not np.isfinite(frequencies).all() or not np.isfinite(modes).all():
        raise ValueError('Nonfinite native normal modes.')
    numbers, native = mode_geometry
    if not np.array_equal(numbers, [atomic_numbers[str(s)] for s in symbols]):
        raise ValueError('Native mode atom order differs from the selected calculation.')
    target = np.asarray(positions, float)
    from kinbot.spatial import proper_alignment
    fit = proper_alignment(native, target)
    if np.max(fit.errors) > 1.e-5:
        raise ValueError('Native mode geometry differs from the selected calculation.')
    return dict(frequencies_cm1=frequencies,
                modes=(np.asarray(modes) @ fit.rotation).tolist(),
                source=str(path), sha256=hashlib.sha256(data).hexdigest(),
                vector_precision='native printed displacement vectors')
