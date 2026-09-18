"""ASE energy/force calculator for a Molpro CCSD(T) geometry step.

Molpro evaluates one geometry per invocation. ASE/Sella, never Molpro OPTG,
chooses the next geometry. See docs/composite_qc_validation.md for the
documented numerical FORCE and PUT,XYZ formats and the offsite acceptance gate.
"""

from pathlib import Path
import re
import shlex
import subprocess

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr, Hartree
import numpy as np


_ENERGY = re.compile(
    r'^\s*SETTING\s+KB_GEOM_ENERGY\s*=\s*'
    r'([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)\s+AU\b',
    re.IGNORECASE | re.MULTILINE,
)
_SAFE_BASIS = re.compile(r'[A-Za-z0-9_+./*()-]+\Z')
_PROCESS_COUNT = re.compile(r'Distribution of processes:\s*nprocs\(total\)=\s*(\d+)',
                            re.IGNORECASE)


def check_process_count(output_path, allocated_ranks):
    """Reject a Molpro launcher that starts more ranks than Slurm requested."""
    output = Path(output_path).read_text(errors='replace')
    counts = {int(match) for match in _PROCESS_COUNT.findall(output)}
    if len(counts) > 1:
        raise RuntimeError('Molpro output reports inconsistent MPI process counts.')
    count = next(iter(counts), None)
    if count is not None and count > allocated_ranks:
        raise RuntimeError(f'Molpro launched {count} MPI processes '
                           f'for an allocation of {allocated_ranks} ranks.')
    return count


def render_input(atoms, *, basis, geometry_name):
    """Render a neutral-singlet conventional CCSD(T) energy and gradient."""
    if not basis or not _SAFE_BASIS.fullmatch(basis):
        raise ValueError('Molpro basis must be one safe basis identifier.')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*\.xyz', geometry_name):
        raise ValueError('Molpro geometry filename must be a safe XYZ basename.')
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
            f'put,xyz,{geometry_name}\n')


def _xyz_atom_order(path, atoms):
    """Map Molpro's printed atom order to ASE's order using its XYZ geometry."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError('Molpro XYZ geometry is incomplete.')
    try:
        count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError('Molpro XYZ geometry lacks its atom count.') from exc
    if count != len(atoms) or len(lines) < count + 2:
        raise ValueError('Molpro XYZ atom count disagrees with ASE geometry.')
    rows = []
    for line in lines[2:2 + count]:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError('Molpro XYZ needs a symbol and three coordinates.')
        atom = re.fullmatch(r'([A-Za-z]{1,2})\d*', fields[0])
        if atom is None:
            raise ValueError('Molpro XYZ contains an invalid atom label.')
        try:
            values = np.array([float(value.replace('D', 'E').replace('d', 'e'))
                               for value in fields[1:]], dtype=float)
        except ValueError as exc:
            raise ValueError('Molpro XYZ contains a nonnumeric coordinate.') from exc
        if not np.isfinite(values).all():
            raise ValueError('Molpro XYZ contains a nonfinite coordinate.')
        rows.append((atom.group(1).capitalize(), values))
    order = []
    unused = set(range(count))
    for symbol, position in zip(atoms.get_chemical_symbols(), atoms.positions):
        matches = [row for row in unused if rows[row][0] == symbol
                   and np.max(np.abs(rows[row][1] - position)) < 1e-4]
        if len(matches) != 1:
            raise ValueError('Molpro XYZ atom order/coordinates cannot be mapped uniquely.')
        row = matches[0]
        order.append(row)
        unused.remove(row)
    return order


def parse_numerical_gradient(output_path, geometry_path, atoms):
    """Read Molpro's numerical dE/dx table (Hartree/Bohr) as ASE forces.

    Molpro 2024 prints the numerical gradient in .out but PUT,XYZGRAD fails
    after FORCE,NUMERICAL. PUT,XYZ supplies the printed atom order instead.
    """
    output = Path(output_path).read_text(errors='replace').splitlines()
    if ('Molpro calculation terminated' not in '\n'.join(output)
            or any('ERROR EXIT' in line for line in output)):
        raise ValueError('Molpro did not terminate normally.')
    headings = [index for index, line in enumerate(output)
                if re.fullmatch(r'\s*Numerical gradient for KB_GEOM_ENERGY\s*',
                                line, re.IGNORECASE)]
    if len(headings) != 1:
        raise ValueError('Expected exactly one KB_GEOM_ENERGY numerical gradient table.')
    start = headings[0]
    headers = [index for index in range(start + 1, min(start + 12, len(output)))
               if re.search(r'\bAtom\s+dE/dx\s+dE/dy\s+dE/dz\b',
                            output[index], re.IGNORECASE)]
    if len(headers) != 1:
        raise ValueError('Molpro numerical gradient column header is missing.')
    gradients = np.empty((len(atoms), 3), dtype=float)
    for atom_index in range(len(atoms)):
        line_index = headers[0] + 1 + atom_index
        fields = output[line_index].split() if line_index < len(output) else []
        if len(fields) != 7 or fields[0] != str(atom_index + 1):
            raise ValueError('Molpro numerical gradient rows are incomplete or unordered.')
        try:
            values = [float(value.replace('D', 'E').replace('d', 'e'))
                      for value in fields[1:]]
        except ValueError as exc:
            raise ValueError('Molpro numerical gradient has a nonnumeric value.') from exc
        if not np.isfinite(values).all():
            raise ValueError('Molpro numerical gradient has a nonfinite value.')
        gradients[atom_index] = values[:3]
    order = _xyz_atom_order(geometry_path, atoms)
    return -gradients[order] * Hartree / Bohr


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
        geometry_name = f'{stem}.xyz'
        stdout_name = f'{stem}.stdout'
        stderr_name = f'{stem}.stderr'
        (self.work_directory / input_name).write_text(
            render_input(atoms, basis=self.basis, geometry_name=geometry_name))
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
        for name in (output_name, log_name, geometry_name):
            path = self.work_directory / name
            if not path.is_file() or not path.stat().st_size:
                raise RuntimeError(f'Molpro force evaluation {stem} lacks nonempty {name}.')
        check_process_count(self.work_directory / output_name, self.nproc)
        self.generated_files.extend([output_name, log_name, geometry_name])
        xml_name = f'{stem}.xml'
        if (self.work_directory / xml_name).is_file():
            self.generated_files.append(xml_name)
        self.results = {'energy': parse_output(self.work_directory / output_name),
                        'forces': parse_numerical_gradient(
                            self.work_directory / output_name,
                            self.work_directory / geometry_name, atoms)}
