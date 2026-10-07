import os
import copy
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


@contextmanager
def gaussian_scratch_environment():
    """Give one Gaussian invocation a writable, collision-free scratch.

    Site Gaussian profiles sometimes export a node path that is absent on a
    different partition.  Treat environment variables as candidate roots,
    validate them on the compute node, and fall back through scheduler/site
    scratch and finally the user's cache or working directory.  Restore the
    caller's environment after Gaussian exits.
    """
    original = os.environ.get('GAUSS_SCRDIR')
    candidates = [
        ('GAUSS_SCRDIR', original),
        ('SLURM_TMPDIR', os.environ.get('SLURM_TMPDIR')),
        ('SCRATCH', os.environ.get('SCRATCH')),
        ('TMPDIR', os.environ.get('TMPDIR')),
        ('HOME', str(Path.home() / '.cache' / 'kinbot' / 'gaussian')),
        ('PWD', os.getcwd()),
    ]
    attempted = []
    scratch = None
    source = None
    for name, value in candidates:
        if not value:
            continue
        root = Path(value).expanduser()
        if root in attempted:
            continue
        attempted.append(root)
        try:
            root.mkdir(parents=True, exist_ok=True)
            if not root.is_dir() or not os.access(root, os.W_OK | os.X_OK):
                continue
            scratch = Path(tempfile.mkdtemp(
                prefix='kinbot-gaussian-', dir=root))
        except OSError:
            continue
        source = name
        break
    if scratch is None:
        paths = ', '.join(str(path) for path in attempted)
        raise RuntimeError('No writable Gaussian scratch root was found. '
                           f'Tried: {paths}')
    os.environ['GAUSS_SCRDIR'] = str(scratch)
    try:
        yield scratch, source
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

    def __init__(self, *args, label='Gaussian', **kwargs):
        FileIOCalculator.__init__(self, *args, label=label, **kwargs)

    def _initialize_profile(self, command):
        # KinBot resolves the Gaussian executable dynamically in calculate();
        # skip ASE's profile/config-file system entirely.
        return None

    def execute(self):
        command = self.command.replace('PREFIX', self.prefix)
        directory = getattr(self, 'directory', '.')
        with gaussian_scratch_environment() as (scratch, source):
            self.last_scratch = {
                'directory': str(scratch), 'source': source,
            }
            proc = subprocess.Popen(command, shell=True, cwd=directory)
            errorcode = proc.wait()
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

        FileIOCalculator.calculate(self, *args, **kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write(self.label + '.com', atoms, properties=properties,
              format='gaussian-in', parallel=False, **self.parameters)

    def read_results(self):
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
