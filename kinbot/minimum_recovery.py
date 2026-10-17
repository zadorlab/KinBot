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

# ANL geometries must be accepted before their dependent energy tasks start.
# Its independent native harmonic parser retains its stricter 1 cm^-1 gate.
ANL_RECOVERY_DEFAULTS = {'attempts': 2, 'displacement_angstrom': 0.02,
                         'difference_step_angstrom': 0.01,
                         'imaginary_threshold_cm_inverse': 1.0}


def validate_anl_recovery(value):
    """Complete and validate the public, molecule-independent ANL controls."""
    if not isinstance(value, dict) or set(value) - ANL_RECOVERY_DEFAULTS.keys():
        raise ValueError('Invalid anl_minimum_recovery controls.')
    options = {**ANL_RECOVERY_DEFAULTS, **value}
    attempts = options['attempts']
    if isinstance(attempts, bool) or not isinstance(attempts, int) or attempts < 0:
        raise ValueError('ANL minimum recovery attempts must be a nonnegative integer.')
    for name, number in options.items():
        if name == 'attempts':
            continue
        if isinstance(number, bool) or not isinstance(number, (int, float)) or not np.isfinite(number) or number <= 0:
            raise ValueError(f'ANL minimum recovery {name} must be positive and finite.')
    return options


def projected_force_modes(atoms, step, report):
    """Central force differences, with rigid translations/rotations removed.

    Work on copies so a failed probe never changes the accepted candidate.
    Store every energy/force/geometry in the report as each probe finishes.
    The original calculator retains its native files, including Molpro's
    numerical-gradient evaluations. This is a real additional curvature
    calculation, not a decision based on the optimizer's approximate Hessian.
    """
    from ase import units
    origin = atoms.positions.copy()
    size = origin.size
    hessian = np.empty((size, size))
    probes = []
    data = {'geometry_angstrom': origin.tolist(), 'step_angstrom': step,
            'probes': probes, 'status': 'running'}
    try:
        for index in range(size):
            forces = []
            for sign in (-1, 1):
                trial = atoms.copy()
                trial.calc = atoms.calc
                trial.positions.flat[index] += sign * step
                force = np.asarray(trial.get_forces(), dtype=float)
                energy = float(trial.get_potential_energy())
                if force.shape != origin.shape or not np.isfinite(force).all() or not np.isfinite(energy):
                    raise ValueError('Nonfinite ANL curvature probe.')
                forces.append(force.ravel())
                probes.append({'coordinate': index, 'sign': sign,
                               'energy_eV': energy, 'forces_eV_angstrom': force.tolist(),
                               'geometry_angstrom': trial.positions.tolist()})
                Path(report).write_text(json.dumps(data, indent=2) + '\n')
            hessian[:, index] = -(forces[1] - forces[0]) / (2. * step)
        hessian = (hessian + hessian.T) / 2.
        masses = atoms.get_masses()
        weights = np.repeat(np.sqrt(masses), 3)
        xyz = origin - np.average(origin, axis=0, weights=masses)
        rigid = [np.tile(axis, len(atoms)) * weights for axis in np.eye(3)]
        rigid += [np.cross(xyz, axis).ravel() * weights for axis in np.eye(3)]
        u, singular, _ = np.linalg.svd(np.array(rigid).T, full_matrices=True)
        rank = int(np.count_nonzero(singular > singular[0] * 1.e-10))
        vibration = u[:, rank:]
        eigenvalues, eigenvectors = np.linalg.eigh(
            vibration.T @ (hessian / np.outer(weights, weights)) @ vibration)
        conversion = units._hbar * units.m / np.sqrt(units._e * units._amu) / units.invcm
        frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * conversion
        modes = ((vibration @ eigenvectors).T / weights).reshape(-1, len(atoms), 3)
        modes /= np.linalg.norm(modes, axis=(1, 2))[:, None, None]
        data.update(status='complete', hessian_eV_angstrom2=hessian.tolist(),
                    rigid_rank=rank, frequencies_cm_inverse=frequencies.tolist(),
                    modes_cartesian=modes.tolist())
        return frequencies, modes
    except Exception as error:
        data.update(status='failed', error=f'{type(error).__name__}: {error}')
        raise
    finally:
        Path(report).write_text(json.dumps(data, indent=2) + '\n')


def optimize_anl_minimum(atoms, optimize, directory, controls):
    """Optimize/check/recover before releasing geometry to dependent tasks.

    A zero attempt budget disables this extra curvature check; the later native
    harmonic task remains mandatory. Call only for unconstrained minima, never
    TS or VRC constrained paths. Each accepted displacement must lower energy.
    """
    from ase.io import write
    options = validate_anl_recovery(controls)
    directory = Path(directory)
    generated = []
    report = {'controls': options, 'attempts': [], 'status': 'running'}
    summary = directory / 'minimum_recovery.json'
    generated.append(summary.name)
    try:
        for attempt in range(options['attempts'] + 1):
            if not optimize():
                raise RuntimeError('ASE geometry optimization did not converge.')
            if not options['attempts']:
                report['status'] = 'disabled'
                return generated
            point = directory / f'minimum_candidate_{attempt}.xyz'
            write(point, atoms)
            generated.append(point.name)
            for name in ('optimization.log', 'optimization.traj'):
                source = directory / name
                if source.is_file():
                    archived = directory / f'minimum_candidate_{attempt}_{name}'
                    shutil.copy2(source, archived)
                    generated.append(archived.name)
            spectrum = directory / f'minimum_curvature_{attempt}.json'
            generated.append(spectrum.name)
            frequencies, modes = projected_force_modes(
                atoms, options['difference_step_angstrom'], spectrum)
            record = {'candidate': point.name, 'curvature': spectrum.name,
                      'lowest_frequency_cm_inverse': float(frequencies[0])}
            report['attempts'].append(record)
            if frequencies[0] >= -options['imaginary_threshold_cm_inverse']:
                report['status'] = 'minimum_curvature_accepted'
                return generated
            if attempt == options['attempts']:
                raise RuntimeError('ANL minimum retains imaginary curvature after bounded recovery.')
            reference = float(atoms.get_potential_energy())
            if not np.isfinite(reference):
                raise ValueError('Nonfinite ANL minimum reference energy.')
            trials = []
            distance = options['displacement_angstrom'] * (attempt + 1)
            record.update(reference_energy_eV=reference,
                          displacement_angstrom=distance, trials=[])
            for sign in (-1, 1):
                trial = atoms.copy()
                trial.calc = atoms.calc
                trial.positions += sign * distance * modes[0]
                energy = float(trial.get_potential_energy())
                if not np.isfinite(energy):
                    raise ValueError('Nonfinite ANL minimum displacement energy.')
                trials.append((energy, trial.positions.copy()))
                record['trials'].append({'sign': sign, 'energy_eV': energy,
                                        'geometry_angstrom': trial.positions.tolist()})
            energy, positions = min(trials, key=lambda item: item[0])
            if energy >= reference - 1.e-6:
                raise RuntimeError('Neither ANL imaginary-mode displacement lowers energy by 1e-6 eV.')
            record['selected_energy_eV'] = energy
            atoms.positions = positions
            with (directory / 'optimization.log').open('a') as stream:
                stream.write(f'Minimum recovery {attempt + 1}: displaced {distance:g} A '
                             f'along mode {frequencies[0]:.3f} cm^-1; reoptimizing.\n')
            summary.write_text(json.dumps(report, indent=2) + '\n')
        raise AssertionError('Unreachable minimum recovery state.')
    except Exception as error:
        report.update(status='failed', error=f'{type(error).__name__}: {error}')
        raise
    finally:
        summary.write_text(json.dumps(report, indent=2) + '\n')


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
