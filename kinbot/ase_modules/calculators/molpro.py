"""ASE energy/force calculator for a Molpro CCSD(T) geometry step.

Molpro evaluates one geometry per invocation. ASE/Sella, never Molpro OPTG,
chooses the next geometry. See docs/composite_qc_validation.md for the
documented FORCE and PUT,XYZGRAD formats and the offsite acceptance gate.
"""

from pathlib import Path
import re
import shlex
import subprocess

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree
import numpy as np


_ENERGY = re.compile(
    r'^\s*SETTING\s+KB_GEOM_ENERGY\s*=\s*'
    r'([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)\s+AU\b',
    re.IGNORECASE | re.MULTILINE,
)
_SAFE_BASIS = re.compile(r'[A-Za-z0-9_+./*()-]+\Z')


def render_input(atoms, *, basis, gradient_name):
    """Render a neutral-singlet conventional CCSD(T) energy and gradient."""
    if not basis or not _SAFE_BASIS.fullmatch(basis):
        raise ValueError('Molpro basis must be one safe basis identifier.')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*\.xyz', gradient_name):
        raise ValueError('Molpro gradient filename must be a safe XYZ basename.')
    xyz = '\n'.join(
        f'{symbol} {x:.12f} {y:.12f} {z:.12f}'
        for symbol, (x, y, z) in zip(atoms.get_chemical_symbols(), atoms.positions)
    )
    return (f'***,KinBot ASE/Sella gradient\n'
            f'symmetry,nosym\norient,noorient\ngeomtyp=xyz\n'
            f'geometry={{\n{len(atoms)}\nKinBot ASE geometry (angstrom)\n{xyz}\n}}\n'
            f'basis={basis}\nset,charge=0\nset,spin=0\n'
            f'gthresh,energy=1.d-9\nrhf\nccsd(t)\n'
            f'kb_geom_energy=energy\n'
            f'forces,numerical,variable=kb_geom_energy,startcmd=rhf\n'
            f'put,xyzgrad,{gradient_name}\n')


def parse_xyzgrad(path, atoms):
    """Read Molpro PUT,XYZGRAD forces (-eV/Angstrom), matching atom positions.

    Molpro can reorder centres by element. Matching coordinates as well as
    symbols prevents a silent force permutation; ambiguous matches fail.
    """
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError('Molpro XYZGRAD is incomplete.')
    try:
        count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError('Molpro XYZGRAD lacks its atom count.') from exc
    if count != len(atoms) or len(lines) < count + 2:
        raise ValueError('Molpro XYZGRAD atom count disagrees with ASE geometry.')
    rows = []
    for line in lines[2:2 + count]:
        fields = line.split()
        if len(fields) != 8:
            raise ValueError('Molpro XYZGRAD needs symbol, XYZ, charge, and three forces.')
        atom = re.fullmatch(r'([A-Za-z]{1,2})\d*', fields[0])
        if atom is None:
            raise ValueError('Molpro XYZGRAD contains an invalid atom label.')
        try:
            values = np.array([float(value.replace('D', 'E').replace('d', 'e'))
                               for value in fields[1:]], dtype=float)
        except ValueError as exc:
            raise ValueError('Molpro XYZGRAD contains a nonnumeric value.') from exc
        if not np.isfinite(values).all():
            raise ValueError('Molpro XYZGRAD contains a nonfinite value.')
        rows.append((atom.group(1).capitalize(), values[:3], values[4:]))
    forces = np.empty((count, 3), dtype=float)
    unused = set(range(count))
    for index, (symbol, position) in enumerate(
            zip(atoms.get_chemical_symbols(), atoms.positions)):
        matches = [row for row in unused if rows[row][0] == symbol
                   and np.max(np.abs(rows[row][1] - position)) < 1e-4]
        if len(matches) != 1:
            raise ValueError('Molpro XYZGRAD atom order/coordinates cannot be mapped uniquely.')
        row = matches[0]
        forces[index] = rows[row][2]
        unused.remove(row)
    return forces


def parse_output(path):
    """Read the explicitly saved undisplaced energy and normal termination."""
    output = Path(path).read_text(errors='replace')
    if 'Molpro calculation terminated' not in output or 'ERROR EXIT' in output:
        raise ValueError('Molpro did not terminate normally.')
    # Molpro 2024.1 prints SETTING for a user assignment, as in the supplied
    # CH4/H2O2 outputs. Use the first such line: numerical FORCE then repeats
    # the energy procedure for displaced geometries in the same invocation.
    match = _ENERGY.search(output)
    if match is None:
        raise ValueError('Undisplaced KB_GEOM_ENERGY missing from Molpro output.')
    energy_hartree = float(match.group(1).replace('D', 'E').replace('d', 'e'))
    if not np.isfinite(energy_hartree):
        raise ValueError('Molpro energy is nonfinite.')
    return energy_hartree * Hartree


class Molpro(Calculator):
    """Run one Molpro force evaluation for each new ASE geometry."""

    implemented_properties = ['energy', 'forces']

    def __init__(self, *, directory, label, method, basis, command='molpro',
                 charge=0, mult=1, nproc=1, stack_mw=256):
        super().__init__()
        if method.upper() != 'CCSD(T)':
            raise ValueError('Molpro ASE geometry currently supports conventional CCSD(T).')
        if charge != 0 or mult != 1:
            raise ValueError('Molpro ASE geometry currently supports a neutral singlet only.')
        if not isinstance(nproc, int) or nproc < 1:
            raise ValueError('Molpro nproc must be positive.')
        if not isinstance(stack_mw, int) or stack_mw < 32:
            raise ValueError('Molpro stack_mw must be at least 32 MW.')
        if not isinstance(label, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', label):
            raise ValueError('Molpro label must be a safe basename.')
        self.work_directory = Path(directory)
        self.work_directory.mkdir(parents=True, exist_ok=True)
        self.step_label = label
        self.basis = basis
        self.command = shlex.split(command)
        if not self.command:
            raise ValueError('Molpro command is empty.')
        self.nproc = nproc
        self.stack_mw = stack_mw
        self.evaluations = 0
        self.generated_files = []

    def calculate(self, atoms=None, properties=('energy',),
                  system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.evaluations += 1
        stem = f'{self.step_label}_step_{self.evaluations:04d}'
        input_name = f'{stem}.inp'
        output_name = f'{stem}.out'
        gradient_name = f'{stem}.xyz'
        stdout_name = f'{stem}.stdout'
        stderr_name = f'{stem}.stderr'
        (self.work_directory / input_name).write_text(
            render_input(atoms, basis=self.basis, gradient_name=gradient_name))
        self.generated_files.append(input_name)
        command = [*self.command, '-g', '-n', str(self.nproc), '-m',
                   str(self.stack_mw), input_name]
        with (self.work_directory / stdout_name).open('w') as stdout, \
                (self.work_directory / stderr_name).open('w') as stderr:
            result = subprocess.run(command, cwd=self.work_directory, stdout=stdout,
                                    stderr=stderr, check=False)
        self.generated_files.extend([stdout_name, stderr_name])
        if result.returncode:
            raise RuntimeError(f'Molpro force evaluation {stem} exited with '
                               f'status {result.returncode}.')
        log_name = f'{stem}.log'
        for name in (output_name, log_name, gradient_name):
            path = self.work_directory / name
            if not path.is_file() or not path.stat().st_size:
                raise RuntimeError(f'Molpro force evaluation {stem} lacks nonempty {name}.')
        self.generated_files.extend([output_name, log_name, gradient_name])
        xml_name = f'{stem}.xml'
        if (self.work_directory / xml_name).is_file():
            self.generated_files.append(xml_name)
        self.results = {'energy': parse_output(self.work_directory / output_name),
                        'forces': parse_xyzgrad(self.work_directory / gradient_name,
                                                atoms)}
