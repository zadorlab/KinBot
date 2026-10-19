import os
import copy
import logging
from contextlib import contextmanager
from pathlib import Path
import shutil
import subprocess
import tempfile
from collections.abc import Iterable
from shutil import which
from typing import Dict, Optional

from ase.calculators.calculator import FileIOCalculator
from kinbot.ase_modules.io.formats import read, write


def qc_scf_failure(path):
    """Recognize a native QC-SCF abort, not syntax, I/O or geometry errors."""
    try:
        with Path(path).open('rb') as stream:
            stream.seek(max(0, Path(path).stat().st_size - 65536))
            text = stream.read().decode(errors='replace').lower()
    except OSError:
        return False
    return ('error termination' in text and 'l508.exe' in text and
            ('no lower point found' in text or 'convergence failure' in text))


def scf_recovery_parameters(parameters):
    """Change the SCF algorithm only; retain the requested density threshold."""
    from kinbot.gaussian_input import route_options
    result = dict(parameters)
    key = next((k for k in result if k.lower() == 'scf'), 'scf')
    options = route_options(result.pop(key, ''))
    replace = {'qc', 'xqc', 'yqc', 'sd', 'ssd', 'dm', 'diis', 'nodiis',
               'cdiis', 'ediis', 'fermi', 'nofermi', 'damp', 'nodamp',
               'varacc', 'novaracc', 'varint', 'novarint', 'incfock', 'noincfock',
               'maxconventionalcycles', 'maxconventional', 'maxcycle', 'maxcycles'}
    kept = [opt for opt in options if opt.split('=', 1)[0].lower() not in replace]
    cycles = 256
    for opt in options:
        if opt.split('=', 1)[0].lower() in {'maxcycle', 'maxcycles'}:
            cycles = max(cycles, int(opt.split('=', 1)[1]))
    result['scf'] = ','.join(kept + ['CDIIS', 'NoVarAcc', 'NoIncFock', f'MaxCycle={cycles}'])
    return result


def fixed_frequency_parameters(parameters):
    """Return Gaussian options for a fixed-geometry Hessian calculation.

    External Cartesian optimizers use NoSymm so gradients retain a stable
    frame.  Gaussian's vibrational analysis still needs its symmetry
    machinery, especially for linear rotors, so use a C1 point-group cap for
    the separate frequency calculation.
    """
    result = dict(parameters)
    result.pop('force', None)
    result.pop('opt', None)
    for key in tuple(result):
        if key.lower() in ('symm', 'symmetry'):
            result.pop(key)
    result['symmetry'] = 'PG=C1'
    result['freq'] = ''
    return result


@contextmanager
def gaussian_scratch_environment():
    """Compatibility context using the common verified-local scratch policy."""
    from kinbot.anl.runtime import _check_scratch_io
    from kinbot.qc_scratch import allocate
    original = os.environ.get('GAUSS_SCRDIR')
    scratch, provenance = allocate(dict(os.environ), 'gaussian', 4096,
                                   probe=_check_scratch_io)
    os.environ['GAUSS_SCRDIR'] = str(scratch)
    try:
        yield scratch, provenance['scratch_source']
    finally:
        if original is None:
            os.environ.pop('GAUSS_SCRDIR', None)
        else:
            os.environ['GAUSS_SCRDIR'] = original
        shutil.rmtree(scratch, ignore_errors=True)


class GaussianDynamics:
    calctype = 'optimizer'
    delete = ['force']
    keyword: Optional[str] = None
    special_keywords: Dict[str, str] = dict()

    def __init__(self, atoms, calc=None):
        self.atoms = atoms
        if calc is not None:
            self.calc = calc
        else:
            if self.atoms.calc is None:
                raise ValueError("{} requires a valid Gaussian calculator "
                                 "object!".format(self.__class__.__name__))

            self.calc = self.atoms.calc

    def todict(self):
        return {'type': self.calctype,
                'optimizer': self.__class__.__name__}

    def delete_keywords(self, kwargs):
        """removes list of keywords (delete) from kwargs"""
        for d in self.delete:
            kwargs.pop(d, None)

    def set_keywords(self, kwargs):
        args = kwargs.pop(self.keyword, [])
        if isinstance(args, str):
            args = [args]
        elif isinstance(args, Iterable):
            args = list(args)

        for key, template in self.special_keywords.items():
            if key in kwargs:
                val = kwargs.pop(key)
                args.append(template.format(val))

        kwargs[self.keyword] = args

    def run(self, **kwargs):
        calc_old = self.atoms.calc
        params_old = copy.deepcopy(self.calc.parameters)

        self.delete_keywords(kwargs)
        self.delete_keywords(self.calc.parameters)
        self.set_keywords(kwargs)

        self.calc.set(**kwargs)
        self.atoms.calc = self.calc

        try:
            self.atoms.get_potential_energy()
        except OSError:
            converged = False
        else:
            converged = True

        atoms = read(self.calc.label + '.log')
        self.atoms.cell = atoms.cell
        self.atoms.positions = atoms.positions

        self.calc.parameters = params_old
        self.calc.reset()
        if calc_old is not None:
            self.atoms.calc = calc_old

        return converged


class GaussianOptimizer(GaussianDynamics):
    keyword = 'opt'
    special_keywords = {
        'fmax': '{}',
        'steps': 'maxcycle={}',
    }


class GaussianIRC(GaussianDynamics):
    keyword = 'irc'
    special_keywords = {
        'direction': '{}',
        'steps': 'maxpoints={}',
    }


class Gaussian(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'dipole']
    command = 'GAUSSIAN < PREFIX.com > PREFIX.log'
    discard_results_on_any_change = True

    def __init__(self, *args, label='Gaussian', scf_recovery=True, **kwargs):
        if scf_recovery not in (True, False, 0, 1):
            raise ValueError('scf_recovery must be true or false.')
        self.scf_recovery = bool(scf_recovery)
        FileIOCalculator.__init__(self, *args, label=label, **kwargs)

    def _initialize_profile(self, command):
        # KinBot resolves the Gaussian executable dynamically in calculate();
        # skip ASE's profile/config-file system entirely.
        return None

    def execute(self):
        command = self.command.replace('PREFIX', self.prefix)
        directory = getattr(self, 'directory', '.')
        from kinbot.anl.runtime import qc_runtime_environment, cleanup_qc_runtime
        from kinbot.qc_scratch import workspace
        env, runtime = qc_runtime_environment('g16', 'gaussian', work_directory=directory)
        self.last_scratch = runtime
        try:
            with workspace(directory, runtime, backend='gaussian',
                           inputs=[self.prefix + '.com'], label=self.prefix) as local:
                proc = subprocess.Popen(command, shell=True, cwd=local, env=env)
                errorcode = proc.wait()
        finally:
            cleanup_qc_runtime(runtime)
        if errorcode:
            raise RuntimeError(
                f'Gaussian exited with error code {errorcode} '
                f'(command: {command})')

    def calculate(self, *args, **kwargs):
        gaussians = ('g16', 'g09', 'g03')
        if 'GAUSSIAN' in self.command:
            for gau in gaussians:
                if which(gau):
                    self.command = self.command.replace('GAUSSIAN', gau)
                    break
            else:
                raise ValueError('Missing Gaussian executable {}'
                                 .format(gaussians))

        try:
            FileIOCalculator.calculate(self, *args, **kwargs)
        except RuntimeError:
            output = Path(self.label + '.log')
            if not self.scf_recovery or not qc_scf_failure(output):
                raise
            # Retry the same geometry and requested electronic state once.
            # No checkpoint/guess replacement or relaxed convergence criterion.
            from kinbot.execution import atomic_json, file_sha256
            archive_root = output.parent / 'restart_gaussian_scf'
            archive_root.mkdir(exist_ok=True)
            archive = Path(tempfile.mkdtemp(prefix=output.stem + '.', dir=archive_root))
            preserved = {}
            paths = [Path(self.label + suffix) for suffix in ('.com', '.log')]
            checkpoint = self.parameters.get('chk')
            if checkpoint:
                path = Path(checkpoint)
                if path.suffix.lower() != '.chk':
                    path = Path(str(path) + '.chk')
                if not path.is_absolute():
                    path = Path(self.directory) / path
                paths.append(path)
            for path in paths:
                if path.is_file():
                    target = archive / path.name
                    shutil.copy2(path, target)
                    preserved[str(path)] = {'archive': str(target), 'sha256': file_sha256(target)}
            original = copy.deepcopy(self.parameters)
            replacement = scf_recovery_parameters(original)
            record = {'status': 'retrying', 'reason': 'Gaussian QC-SCF failure in Link 508',
                      'original_scf': next((v for k, v in original.items() if k.lower() == 'scf'), ''),
                      'retry_scf': replacement['scf'], 'preserved_files': preserved}
            report = Path(self.label + '.scf_recovery.json')
            atomic_json(report, record)
            logging.getLogger('KinBot').warning(
                '%s: QC-SCF failed; one CDIIS retry at unchanged theory and convergence. Evidence: %s',
                self.label, report)
            self.parameters.update(replacement)
            for key in tuple(self.parameters):
                if key.lower() == 'scf' and key != 'scf':
                    del self.parameters[key]
            try:
                FileIOCalculator.calculate(self, *args, **kwargs)
            except Exception as error:
                record.update(status='failed', error=f'{type(error).__name__}: {error}')
                raise
            else:
                record.update(status='complete', output_sha256=file_sha256(output))
            finally:
                self.parameters = original
                atomic_json(report, record)

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write(self.label + '.com', atoms, properties=properties,
              format='gaussian-in', parallel=False, **self.parameters)

    def read_results(self):
        # A rejected native SCF must not leave a previous geometry's cached
        # energy/forces available to the caller.
        self.results = {}
        self.calc = None
        output = read(self.label + '.log', format='gaussian-out')
        self.calc = output.calc
        self.results = output.calc.results

    # Method(s) defined in the old calculator, added here for
    # backwards compatibility
    def clean(self):
        for suffix in ['.com', '.chk', '.log']:
            try:
                os.remove(os.path.join(self.directory, self.label + suffix))
            except OSError:
                pass

    def get_version(self):
        raise NotImplementedError  # not sure how to do this yet
