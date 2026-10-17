"""ASE energy/force calculator for Molpro geometry steps.

Molpro evaluates one geometry per invocation. ASE/Sella, never Molpro OPTG,
chooses the next geometry. Each step preserves its native input and output.
Conventional CCSD(T) uses numerical gradients.  A VRC relaxed path can use
an analytic RS2-CASPT2 gradient or a numerical gradient of an explicitly
selected CASPT2-F12 or Davidson-corrected MRCI energy.
"""

from pathlib import Path
from copy import deepcopy
import re
import shlex
import subprocess

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr, Hartree
import numpy as np

from kinbot.anl.runtime import cleanup_qc_runtime, qc_runtime_environment
from kinbot.anl.results import validate_molpro_triples
from kinbot.molpro import (
    _vrc_multireference_method, stage_vrc_state_model,
    molpro_wfu_name,
    validate_molpro_convergence,
    validate_vrc_state_model,
    vrc_multireference_energy_contract,
)
from kinbot.theory import fortran_threshold
from kinbot.execution import file_sha256


def _mr_gradient_variable(contract, state_model, mult):
    selection = (state_model or {}).get('energy_selection', {})
    if selection.get(mult, selection.get(str(mult))) == 'lowest':
        raise ValueError('Lowest-manifold selection is supported for energies only; '
                         'a geometry gradient requires an explicit root.')
    if ((state_model or {}).get('correlation_options', {}).get('state_policy') == 'simultaneous'
            and contract['method'].startswith('MRCI')):
        roots = state_model['selected_roots']
        return f'ENERGD({roots.get(mult, roots.get(str(mult)))})'
    return 'ENERGD' if contract['method'] == 'MRCI+Q' else contract['energy_variable']


_ENERGY = re.compile(
    r'^\s*SETTING\s+KB_GEOM_ENERGY\s*=\s*'
    r'([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)\s+(?:AU|HARTREE)\b',
    re.IGNORECASE | re.MULTILINE,
)
_SAFE_BASIS = re.compile(r'[A-Za-z0-9_+./*()-]+\Z')
_PROCESS_COUNT = re.compile(r'Distribution of processes:\s*nprocs\(total\)=\s*(\d+)',
                            re.IGNORECASE)
_CASPT2 = re.compile(r'caspt2\((\d+),(\d+)\)\Z', re.IGNORECASE)


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


def _caspt2_space(method, electron_count):
    match = _CASPT2.fullmatch(str(method).strip())
    if match is None:
        return None
    active_electrons, active_orbitals = map(int, match.groups())
    if active_electrons < 1 or active_orbitals < 1:
        raise ValueError('CASPT2 active electrons and orbitals must be positive.')
    inactive_electrons = electron_count - active_electrons
    if inactive_electrons < 0 or inactive_electrons % 2:
        raise ValueError(
            f'{method} is incompatible with {electron_count} electrons: the '
            'inactive electron count must be nonnegative and even.')
    closed = inactive_electrons // 2
    return active_electrons, active_orbitals, closed, closed + active_orbitals


def render_input(atoms, *, basis, geometry_name, gradient_name=None, charge=0,
                 mult=1, method='CCSD(T)', effective_method=None,
                 caspt2_shift=0.1, caspt2_ipea=0.25,
                 state_model=None,
                 orbital_file=None,
                 gradient_mode='auto', scf_convergence=1.e-10,
                 correlation_convergence=1.e-8):
    """Render a CCSD(T), CASPT2, or Davidson-corrected MRCI gradient.

    Molpro's RHF step supplies RHF or ROHF orbitals according to multiplicity.
    Bare ``UCCSD(T)`` then uses the restricted-reference unrestricted-CC
    formalism.  ``UHF_UCCSD=1`` must not be added: it forces Molpro's separate
    UHF-UCC engine and produced a zero triples correction for closed-shell
    ethane in Molpro 2024.1.
    """
    if not basis or not _SAFE_BASIS.fullmatch(basis):
        raise ValueError('Molpro basis must be one safe basis identifier.')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*\.xyz', geometry_name):
        raise ValueError('Molpro geometry filename must be a safe XYZ basename.')
    if isinstance(charge, bool) or not isinstance(charge, int):
        raise ValueError('Molpro charge must be an integer.')
    if isinstance(mult, bool) or not isinstance(mult, int) or mult < 1:
        raise ValueError('Molpro multiplicity must be a positive integer.')
    requested_method = str(method).strip()
    electron_count = int(sum(atoms.numbers)) - charge
    contract = vrc_multireference_energy_contract(requested_method)
    if contract is not None:
        # Never differentiate the final root while reporting MIN(all roots).
        _mr_gradient_variable(contract, state_model, mult)
    if gradient_mode not in ('auto', 'analytic', 'numerical'):
        raise ValueError('Molpro gradient_mode must be auto, analytic, or '
                         'numerical.')
    analytic_available = (
        contract is not None and contract['method'] == 'CASPT2'
        and not contract['explicitly_correlated'])
    if gradient_mode == 'analytic' and not analytic_available:
        raise ValueError(
            f'{requested_method} has no supported analytic Molpro gradient.')
    analytic_caspt2 = analytic_available and gradient_mode != 'numerical'
    if contract is None:
        cc_command = {None: 'uccsd(t)', 'HF': '', 'CCSD': 'uccsd'}.get(
            effective_method)
        if requested_method.upper() != 'CCSD(T)' or cc_command is None:
            raise ValueError('Molpro ASE geometry method is unsupported.')
    else:
        if effective_method is not None:
            raise ValueError(
                'Multireference geometry does not accept effective_method.')
        # Validate the active space before a job can be submitted.  The
        # renderer below performs the same parity check while constructing
        # the exact Molpro program and selected energy variable.
        _caspt2_space(
            f'caspt2({contract["active_electrons"]},'
            f'{contract["active_orbitals"]})', electron_count)
        validate_vrc_state_model(
            state_model, active_electrons=contract['active_electrons'],
            active_orbitals=contract['active_orbitals'],
            nelectron=electron_count, atom_count=len(atoms))
        if (analytic_caspt2 and
                (gradient_name is None or not re.fullmatch(
                    r'[A-Za-z0-9][A-Za-z0-9_.-]*\.xyzgrad',
                    gradient_name))):
            raise ValueError(
                'Analytic CASPT2 geometry requires a safe XYZGRAD filename.')
        if (isinstance(caspt2_shift, bool)
                or not isinstance(caspt2_shift, (int, float))
                or not np.isfinite(caspt2_shift) or caspt2_shift < 0):
            raise ValueError('CASPT2 level shift must be a finite nonnegative number.')
        if (isinstance(caspt2_ipea, bool)
                or not isinstance(caspt2_ipea, (int, float))
                or not np.isfinite(caspt2_ipea) or caspt2_ipea < 0):
            raise ValueError('CASPT2 IPEA shift must be a finite nonnegative number.')
    scf_threshold = fortran_threshold(
        scf_convergence, 'Molpro SCF convergence')
    correlation_threshold = fortran_threshold(
        correlation_convergence, 'Molpro correlation convergence')
    xyz = '\n'.join(
        f'{symbol} {x:.12f} {y:.12f} {z:.12f}'
        for symbol, (x, y, z) in zip(atoms.get_chemical_symbols(), atoms.positions)
    )
    header = (f'***,KinBot ASE/Sella gradient\n'
            f'symmetry,nosym\norient,noorient\ngeomtyp=xyz\n'
            f'geometry={{\n{len(atoms)}\nKinBot ASE geometry (angstrom)\n{xyz}\n}}\n'
            f'basis={basis}\nset,charge={charge}\nset,spin={mult - 1}\n'
            f'gthresh,energy={scf_threshold},coeff={correlation_threshold},'
            f'civec={correlation_threshold}\n')
    selection = ((state_model or {}).get('orbital_selection')
                 if isinstance(state_model, dict) else None)
    if (orbital_file is not None
            and not (isinstance(selection, dict)
                     and selection.get('mode') == 'molpro_record')):
        if (len(orbital_file) > 31 or
                not re.fullmatch(r'[a-z0-9][a-z0-9_.-]*\.wfu', orbital_file)):
            raise ValueError('Molpro orbital checkpoint must be a short lowercase WFU basename.')
        header += f'file,2,{orbital_file}\n'
    header += 'gprint,orbitals,orben,civector\n'
    if contract is not None:
        if analytic_caspt2 and state_model is None:
            active_electrons, active_orbitals, closed, occupied = \
                _caspt2_space(requested_method, electron_count)
            method_block = (
                '{multi,\n'
                f'occ,{occupied}\nclosed,{closed}\n'
                f'wf,{electron_count},1,{mult - 1}\nmaxit,50;}}\n\n'
                f'{{rs2,shift={float(caspt2_shift):g},'
                f'ipea={float(caspt2_ipea):g}}}\n'
                '! KinBot energy contract: CASPT2/RS2 = ENERGY\n')
            selected = 'energy'
        else:
            method_block, selected = _vrc_multireference_method(
                requested_method, electron_count, spin=mult - 1,
                shift=caspt2_shift, ipea=caspt2_ipea,
                state_model=state_model)
        body = method_block + f'kb_geom_energy={selected}\n'
        if analytic_caspt2:
            body += f'forces\nput,xyzgrad,{gradient_name}\n'
        else:
            # Differentiate the selected native variable directly.  For
            # MRCI this is ENERGD(1), never ENERGY, so the geometry path
            # includes the relaxed-reference Davidson correction.
            # Match the selected calculated-state slot, including explicitly
            # requested simultaneous MRCI roots. F12 has no CCSD A/B variants here.
            gradient_variable = _mr_gradient_variable(contract, state_model, mult).casefold()
            body += (f'forces,numerical,variable={gradient_variable},'
                     'startcmd=multi\n')
        body += f'put,xyz,{geometry_name}\n'
        return header + body
    return (header + f'rhf\n{cc_command}\n'
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


def parse_numerical_gradient(output_path, geometry_path, atoms,
                             variable='KB_GEOM_ENERGY'):
    """Read Molpro's numerical dE/dx table (Hartree/Bohr) as ASE forces.

    Molpro 2024 prints the numerical gradient in .out but PUT,XYZGRAD fails
    after FORCE,NUMERICAL. PUT,XYZ supplies the printed atom order instead.
    """
    output = Path(output_path).read_text(errors='replace').splitlines()
    if ('Molpro calculation terminated' not in '\n'.join(output)
            or any('ERROR EXIT' in line for line in output)):
        raise ValueError('Molpro did not terminate normally.')
    headings = [index for index, line in enumerate(output)
                if re.fullmatch(
                    rf'\s*Numerical gradient for {re.escape(variable)}\s*',
                    line, re.IGNORECASE)]
    if len(headings) != 1:
        raise ValueError(
            f'Expected exactly one {variable} numerical gradient table.')
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


def parse_xyzgrad(path, atoms):
    """Read Molpro PUT,XYZGRAD forces in eV/angstrom in ASE atom order."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < len(atoms) + 2:
        raise ValueError('Molpro XYZGRAD file is incomplete.')
    try:
        count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError('Molpro XYZGRAD lacks its atom count.') from exc
    if count != len(atoms):
        raise ValueError('Molpro XYZGRAD atom count disagrees with geometry.')
    rows = []
    for line in lines[2:2 + count]:
        fields = line.split()
        if len(fields) != 8:
            raise ValueError('Molpro XYZGRAD row must contain a symbol and seven values.')
        atom = re.fullmatch(r'([A-Za-z]{1,2})\d*', fields[0])
        if atom is None:
            raise ValueError('Molpro XYZGRAD contains an invalid atom label.')
        try:
            values = np.asarray([
                float(value.replace('D', 'E').replace('d', 'e'))
                for value in fields[1:]], dtype=float)
        except ValueError as exc:
            raise ValueError('Molpro XYZGRAD contains a nonnumeric value.') from exc
        if not np.isfinite(values).all():
            raise ValueError('Molpro XYZGRAD contains a nonfinite value.')
        rows.append((atom.group(1).capitalize(), values[:3], values[4:]))
    order, unused = [], set(range(count))
    for symbol, position in zip(atoms.get_chemical_symbols(), atoms.positions):
        matches = [row for row in unused if rows[row][0] == symbol
                   and np.max(np.abs(rows[row][1] - position)) < 1e-4]
        if len(matches) != 1:
            raise ValueError('Molpro XYZGRAD atom order cannot be mapped uniquely.')
        order.append(matches[0])
        unused.remove(matches[0])
    # Molpro documents the XYZGRAD vector in -eV/angstrom, i.e. force units.
    return np.asarray([rows[index][2] for index in order], dtype=float)


def parse_output(path, *, method='CCSD(T)', effective_method=None):
    """Read the explicitly saved undisplaced energy and normal termination."""
    output = Path(path).read_text(errors='replace')
    if 'Molpro calculation terminated' not in output or 'ERROR EXIT' in output:
        raise ValueError('Molpro did not terminate normally.')
    contract = vrc_multireference_energy_contract(method)
    validate_molpro_convergence(output, allow_scf_guess_failure=contract is not None)
    if contract is not None:
        f12 = contract['explicitly_correlated']
        if contract['method'].startswith('CASPT2'):
            marker = r'RS2-F12\b' if f12 else r'RS2(?=\s|\()'
        else:
            marker = (r'(?:CI|MRCI)-F12\b' if f12
                      else r'(?:CI|MRCI)(?=\s|\()')
        if not re.search(rf'PROGRAM\s*\*\s*{marker}', output,
                         re.IGNORECASE):
            raise ValueError(
                f'Molpro output lacks the requested {contract["method"]} '
                'calculation.')
    elif effective_method is None:
        validate_molpro_triples(output, scaled=False)
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
                 charge=0, mult=1, nproc=1, stack_mw=256,
                 scratch_min_mb=4096, effective_method=None,
                 caspt2_shift=0.1, caspt2_ipea=0.25,
                 state_model=None,
                 gradient_mode='auto', scf_convergence=1.e-10,
                 correlation_convergence=1.e-8):
        super().__init__()
        if (method.upper() != 'CCSD(T)'
                and vrc_multireference_energy_contract(method) is None):
            raise ValueError(
                'Molpro ASE geometry supports CCSD(T), CASPT2 variants, '
                'MRCI+Q, and MRCI-F12+Q active-space methods.')
        if effective_method not in (None, 'HF', 'CCSD'):
            raise ValueError('Molpro effective geometry method is unsupported.')
        if isinstance(charge, bool) or not isinstance(charge, int):
            raise ValueError('Molpro charge must be an integer.')
        if isinstance(mult, bool) or not isinstance(mult, int) or mult < 1:
            raise ValueError('Molpro multiplicity must be a positive integer.')
        if not isinstance(nproc, int) or nproc < 1:
            raise ValueError('Molpro nproc must be positive.')
        if not isinstance(stack_mw, int) or stack_mw < 32:
            raise ValueError('Molpro stack_mw must be at least 32 MW.')
        if (isinstance(scratch_min_mb, bool)
                or not isinstance(scratch_min_mb, int) or scratch_min_mb < 1):
            raise ValueError('Molpro scratch_min_mb must be positive.')
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
        self.scratch_min_mb = scratch_min_mb
        self.charge = charge
        self.mult = mult
        self.effective_method = effective_method
        self.method = method
        self.caspt2_shift = caspt2_shift
        self.caspt2_ipea = caspt2_ipea
        self.state_model = validate_vrc_state_model(state_model)
        self._orbital_checkpoint_model = deepcopy(self.state_model)
        # Validate at construction so a bad public input never reaches the
        # scheduler, including conventional CC jobs that do not use CASPT2.
        fortran_threshold(scf_convergence, 'Molpro SCF convergence')
        fortran_threshold(correlation_convergence,
                          'Molpro correlation convergence')
        self.scf_convergence = scf_convergence
        self.correlation_convergence = correlation_convergence
        if gradient_mode not in ('auto', 'analytic', 'numerical'):
            raise ValueError('Molpro gradient_mode must be auto, analytic, '
                             'or numerical.')
        self.gradient_mode = gradient_mode
        contract = vrc_multireference_energy_contract(method)
        if contract is not None:
            _mr_gradient_variable(contract, self.state_model, self.mult)
        analytic_available = (
            contract is not None and contract['method'] == 'CASPT2'
            and not contract['explicitly_correlated'])
        if gradient_mode == 'analytic' and not analytic_available:
            raise ValueError(
                f'{method} has no supported analytic Molpro gradient.')
        self.analytic_caspt2 = (
            analytic_available and gradient_mode != 'numerical')
        self.gradient_variable = (
            _mr_gradient_variable(contract, self.state_model, self.mult) if contract is not None
            and not self.analytic_caspt2 else 'KB_GEOM_ENERGY')
        self.evaluations = 0
        self.generated_files = []

    @property
    def orbital_checkpoint_model(self):
        """Return the last hash-bound active-orbital model for provenance."""
        return deepcopy(self._orbital_checkpoint_model)

    def calculate(self, atoms=None, properties=('energy',),
                  system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.evaluations += 1
        stem = f'{self.step_label}_step_{self.evaluations:04d}'
        input_name = f'{stem}.inp'
        output_name = f'{stem}.out'
        geometry_name = f'{stem}.xyz'
        gradient_name = f'{stem}.xyzgrad' if self.analytic_caspt2 else None
        stdout_name = f'{stem}.stdout'
        stderr_name = f'{stem}.stderr'
        state_model = stage_vrc_state_model(
            self._orbital_checkpoint_model, self.work_directory, stem)
        selection = ((state_model or {}).get('orbital_selection')
                     if isinstance(state_model, dict) else None)
        if (isinstance(selection, dict)
                and selection.get('mode') == 'molpro_record'):
            orbital_name = Path(selection['path']).name
        else:
            orbital_name = molpro_wfu_name(f'{stem}.orbitals')
        (self.work_directory / input_name).write_text(
            render_input(atoms, basis=self.basis, geometry_name=geometry_name,
                         gradient_name=gradient_name, method=self.method,
                         charge=self.charge, mult=self.mult,
                         effective_method=self.effective_method,
                         caspt2_shift=self.caspt2_shift,
                         caspt2_ipea=self.caspt2_ipea,
                         state_model=state_model,
                         orbital_file=orbital_name,
                         gradient_mode=self.gradient_mode,
                         scf_convergence=self.scf_convergence,
                         correlation_convergence=
                             self.correlation_convergence))
        self.generated_files.append(input_name)
        command = [*self.command, '-W', str(self.work_directory.resolve()),
                   '-g', '-n', str(self.nproc), '-m',
                   str(self.stack_mw), input_name]
        child_env, runtime = qc_runtime_environment(
            command[0], 'molpro', work_directory=self.work_directory,
            scratch_min_mb=self.scratch_min_mb)
        try:
            with (self.work_directory / stdout_name).open('w') as stdout, \
                    (self.work_directory / stderr_name).open('w') as stderr:
                result = subprocess.run(command, cwd=self.work_directory, stdout=stdout,
                                        stderr=stderr, check=False, env=child_env)
        finally:
            cleanup_qc_runtime(runtime)
        self.generated_files.extend([stdout_name, stderr_name])
        if result.returncode:
            raise RuntimeError(f'Molpro force evaluation {stem} exited with '
                               f'status {result.returncode}.')
        log_name = f'{stem}.log'
        required = [output_name, log_name, geometry_name]
        if vrc_multireference_energy_contract(self.method) is not None:
            required.append(orbital_name)
        if gradient_name is not None:
            required.append(gradient_name)
        for name in required:
            path = self.work_directory / name
            if not path.is_file() or not path.stat().st_size:
                raise RuntimeError(f'Molpro force evaluation {stem} lacks nonempty {name}.')
        check_process_count(self.work_directory / output_name, self.nproc)
        self.generated_files.extend(required)
        if vrc_multireference_energy_contract(self.method) is not None:
            checkpoint = self.work_directory / orbital_name
            if not checkpoint.is_file() or not checkpoint.stat().st_size:
                raise RuntimeError(
                    f'Molpro force evaluation {stem} did not save its '
                    'active-orbital checkpoint.')
            updated = deepcopy(state_model) if state_model is not None else {
                'state_average': [{
                    'multiplicity': self.mult, 'roots': 1, 'weights': [1.],
                }],
                'selected_roots': {self.mult: 1},
                'canonical_orbitals': True,
            }
            updated['orbital_selection'] = {
                'mode': 'molpro_record',
                'path': str(checkpoint.resolve()),
                'record': '2140.2',
                'sha256': file_sha256(checkpoint),
            }
            self._orbital_checkpoint_model = validate_vrc_state_model(updated)
        xml_name = f'{stem}.xml'
        if (self.work_directory / xml_name).is_file():
            self.generated_files.append(xml_name)
        self.results = {'energy': parse_output(
                            self.work_directory / output_name,
                            method=self.method,
                            effective_method=self.effective_method),
                        'forces': (parse_xyzgrad(
                            self.work_directory / gradient_name, atoms)
                            if self.analytic_caspt2 else parse_numerical_gradient(
                                self.work_directory / output_name,
                                self.work_directory / geometry_name, atoms,
                                variable=self.gradient_variable))}
