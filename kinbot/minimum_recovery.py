"""Bounded descent from a frequency-rejected Gaussian minimum candidate.

This changes geometry, never the electronic method, multiplicity, or frequency
acceptance criterion. Each rejected point and both trial calculations remain
available for review. A new optimization and Hessian are required afterward.
"""
from pathlib import Path
import json
import shutil
import tempfile

import numpy as np


RECOVERY_DEFAULTS = {'minimum_recovery_attempts': 2,
                     'minimum_recovery_displacement': 0.02,
                     'imagfreq_threshold': 50.}


def worker_sella_kwargs(parameters):
    """Carry controls through the existing template field, not new fields.

    Live controllers may retain the previous qc.py formatter after git pull.
    Templates also accept its plain Sella options and use defaults. The private
    control entry is removed before constructing Sella.
    """
    options = dict(parameters['sella_kwargs'])
    options['_kinbot_worker_controls'] = {
        key: parameters.get(key, default) for key, default in RECOVERY_DEFAULTS.items()}
    return options


def pop_worker_controls(sella_kwargs):
    controls = dict(RECOVERY_DEFAULTS)
    controls.update(sella_kwargs.pop('_kinbot_worker_controls', {}))
    return controls


def aligned_mode(source, target, mode):
    """Rotate a native normal mode into the optimizer's Cartesian frame.

    Atom order must already agree. Use a proper rotation, and reject a geometry
    mismatch rather than applying a mode from another point or a reflection.
    """
    source, target, mode = (np.asarray(x, dtype=float) for x in (source, target, mode))
    if (source.shape != target.shape or mode.shape != target.shape
            or target.ndim != 2 or target.shape[1] != 3
            or not all(np.isfinite(x).all() for x in (source, target, mode))):
        raise ValueError('Invalid geometry or normal mode for minimum recovery.')
    x, y = source - source.mean(axis=0), target - target.mean(axis=0)
    u, _, vt = np.linalg.svd(x.T @ y)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    if np.max(np.linalg.norm(x @ rotation - y, axis=1)) > 1.e-4:
        raise ValueError('Frequency geometry differs from the optimization geometry.')
    norm = np.linalg.norm(mode)
    if norm < 1.e-10:
        raise ValueError('Missing or zero imaginary normal mode.')
    return mode @ rotation / norm


def gaussian_minimum_restart(mol, parameters, label, attempt, displacement):
    """Try both mode signs with independent checkpoint copies; select descent.

    ``displacement`` is the Euclidean length of the full Cartesian displacement
    in Angstrom. The caller controls the finite attempt budget and must repeat
    optimization/frequencies. Any native error propagates instead of being
    disguised as a successful recovery.
    """
    from kinbot.ase_modules.calculators.gaussian import Gaussian
    from kinbot.ase_modules.io.formats import read
    from kinbot.reader_gauss import read_freq, read_imag_mode

    native = read(label + '.log', format='gaussian-out')
    if list(native.symbols) != list(mol.symbols):
        raise ValueError('Frequency atom order differs from the optimization.')
    frequencies = np.asarray(read_freq(label + '.log', mol.symbols), dtype=float)
    if (not frequencies.size or not np.isfinite(frequencies).all()
            or frequencies[0] >= 0 or np.argmin(frequencies) != 0):
        raise ValueError('No lowest imaginary mode available for minimum recovery.')
    direction = aligned_mode(native.positions, mol.positions,
                             read_imag_mode(label, len(mol)))
    energy = float(native.get_potential_energy())
    if not np.isfinite(energy):
        raise ValueError('Nonfinite reference energy for minimum recovery.')
    checkpoint = Path(label + '.chk')
    if not checkpoint.is_file():
        raise FileNotFoundError(checkpoint)
    archive_root = Path(label + '_minimum_recovery')
    archive_root.mkdir(parents=True, exist_ok=True)
    # A resubmitted worker must not overwrite an earlier recovery's evidence.
    archive = Path(tempfile.mkdtemp(prefix=f'attempt_{attempt}_', dir=archive_root))
    for suffix in ('.com', '.log', '.chk', '.traj', '_sella.log'):
        source = Path(label + suffix)
        if source.is_file():
            shutil.copy2(source, archive / source.name)
    options = {key: value for key, value in parameters.items()
               if key.lower() not in ('label', 'chk', 'oldchk', 'guess', 'force',
                                      'freq', 'opt', 'symm', 'symmetry')}
    options.update(Symm='None', guess='Read')
    candidates = []
    report = {'attempt': attempt, 'frequency_cm-1': float(frequencies[0]),
              'archive': str(archive),
              'displacement_angstrom': displacement, 'reference_energy_eV': energy,
              'trials': [], 'status': 'running'}
    report_path = archive / 'recovery.json'
    try:
        for sign, name in ((1, 'plus'), (-1, 'minus')):
            directory = archive / name
            directory.mkdir()
            shutil.copy2(checkpoint, directory / 'probe.chk')
            trial = mol.copy()
            trial.positions += sign * displacement * direction
            trial.calc = Gaussian(label=str(directory / 'probe'), chk='probe.chk', **options)
            forces = trial.get_forces()
            trial_energy = float(trial.get_potential_energy())
            if not np.isfinite(trial_energy) or not np.isfinite(forces).all():
                raise ValueError('Nonfinite minimum-recovery trial result.')
            report['trials'].append({'name': name, 'energy_eV': trial_energy,
                                    'geometry_angstrom': trial.positions.tolist()})
            candidates.append((trial_energy, name, trial.positions.copy()))
        selected_energy, selected, positions = min(candidates, key=lambda item: item[0])
        # Do not restart on a numerical tie or move uphill. This small guard
        # is an energy-decrease test, not a frequency acceptance threshold.
        if selected_energy >= energy - 1.e-6:
            raise ValueError('Neither imaginary-mode displacement lowers the energy '
                             'by more than 1e-6 eV; minimum recovery stopped.')
        report.update(status='descent_selected', selected=selected,
                      energy_change_eV=selected_energy-energy)
        shutil.copy2(archive / selected / 'probe.chk', checkpoint)
        mol.positions = positions
        mol.calc = Gaussian(**parameters)
        return report
    except Exception as error:
        report.update(status='failed', error=f'{type(error).__name__}: {error}')
        raise
    finally:
        report_path.write_text(json.dumps(report, indent=2) + '\n')
