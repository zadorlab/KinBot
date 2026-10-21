"""Molpro sampling bridge for the pinned ROTD_py worker interface.

ROTD_py owns sampling, scheduling, restart files, and flux corrections. This
bridge supplies its Molpro energy evaluations with KinBot's reviewed orbital
model and CBS algebra, without modifying the pinned external package.
"""

from copy import deepcopy
from pathlib import Path
import re
import subprocess

import numpy as np
from ase.units import Hartree
from rotd_py.molpro.molpro import Molpro as LegacyMolpro
from rotd_py.sample import sample as legacy_sample
from rotd_py.sample.multi_sample import MultiSample

from kinbot.anl.extrapolation import two_point_cbs
from kinbot.ase_modules.calculators.molpro import energy_assignment, parse_output
from kinbot.molpro import (_vrc_multireference_method, stage_vrc_state_model,
                           initialize_molpro_restart)
from kinbot.theory import fortran_threshold
from kinbot.vrc_correction import sampling_settings, _basis_cardinal
from kinbot.anl.runtime import qc_runtime_environment, cleanup_qc_runtime


class LocalMolpro(LegacyMolpro):
    """Keep even legacy RotdPy sampling on the common local-scratch path."""

    def run(self):
        MolproSample._run_basis_job(self, self)


class MolproSample(LocalMolpro):
    """Evaluate each requested basis independently, then combine in hartree."""

    def create_input(self):
        model = deepcopy(self.calc['kinbot_vrc'])
        symbols = list(self.mol.get_chemical_symbols())
        order = model['sampling_atom_order']
        if (sorted(order) != list(range(len(symbols)))
                or [model['parent_symbols'][index] for index in order] != symbols):
            raise ValueError('ROTD sample atom order differs from its orbital model.')
        # WFU coefficients and atom-labelled AVAS centers refer to the parent
        # ordering, whereas ROTD concatenates the two optimized fragments.
        atoms = self.mol[np.argsort(order)]
        electron_count = int(sum(atoms.numbers)) - model['charge']
        mult = model['multiplicity']
        if ((electron_count - mult + 1) % 2 or mult - 1 > electron_count):
            raise ValueError('ROTD sample has an incompatible molecular spin.')
        level = sampling_settings({
            'vrc_tst_sample_method': model['method'],
            'vrc_tst_sample_basis': model['basis'],
            'vrc_tst_sample_extrapolation_power': model['extrapolation_power'],
        })
        bases = level['basis'] if isinstance(level['basis'], list) else [level['basis']]
        self.kinbot_level = level
        self.basis_jobs = []
        rows = '\n'.join(
            f'{symbol} {x:.12f} {y:.12f} {z:.12f}'
            for symbol, (x, y, z) in zip(atoms.get_chemical_symbols(), atoms.positions))
        scf = fortran_threshold(model['scf_convergence'], 'VRC SCF convergence')
        correlation = fortran_threshold(model['correlation_convergence'],
                                         'VRC correlation convergence')
        for index, basis in enumerate(bases):
            if not re.fullmatch(r'[A-Za-z0-9_+./*()-]+', basis):
                raise ValueError('ROTD sampling basis must be a safe identifier.')
            name = f'{self.name}_basis_{index}'
            selection = stage_vrc_state_model(
                model['state_models_by_basis'][basis], Path.cwd(), name)
            rendered = _vrc_multireference_method(
                level['method'], electron_count, spin=mult - 1,
                shift=model['shift'], ipea=model['ipea'], state_model=selection)
            if rendered is None:
                raise ValueError('The profiled ROTD sampling bridge requires '
                                 'CASPT2 or Davidson-corrected MRCI.')
            method, selected_energy = rendered
            child = LegacyMolpro(name, atoms, self.calc.copy())
            Path(f'{name}.inp').write_text(initialize_molpro_restart(
                f'***,KinBot ROTD sampling\nmemory,{self.memory_in_MW},m\n'
                'angstrom\norient,noorient\nsymmetry,nosym\n'
                f'set,charge={model["charge"]}\nset,spin={mult - 1}\n'
                f'geometry={{\n{len(atoms)}\n{name}\n{rows}\n}}\n'
                + f'basis={basis}\n'
                +
                f'gthresh,energy={scf},coeff={correlation},civec={correlation}\n'
                + method + energy_assignment(selected_energy) +
                f'MYENERGY={selected_energy}\n'))
            self.basis_jobs.append(child)

    def run(self):
        for child in self.basis_jobs:
            self._run_basis_job(child)

    def _run_basis_job(self, child):
        env, runtime = qc_runtime_environment('molpro', 'molpro', work_directory=Path.cwd())
        command = ['molpro', '-d', runtime['scratch_directory'], '-W', str(Path.cwd()), '-t', '1',
                   '-n', str(child.procs), f'{child.name}.inp']
        try:
            from kinbot.qc_scratch import workspace, molpro_command
            with workspace(Path.cwd(), runtime, backend='molpro',
                           inputs=[f'{child.name}.inp'], label=child.name) as local:
                command = molpro_command(command, runtime, local)
                with (local / f'{child.name}.stdout').open('w') as stdout, \
                        (local / f'{child.name}.stderr').open('w') as stderr:
                    result = subprocess.run(command, cwd=local, env=env, stdout=stdout,
                                            stderr=stderr, check=False)
            from kinbot.molpro_diagnostics import write_native_report
            write_native_report(Path(f'{child.name}.out'))
            if result.returncode:
                raise RuntimeError(f'Molpro sample {child.name} exited with status {result.returncode}.')
        finally:
            cleanup_qc_runtime(runtime)

    def read_energy(self):
        energies = []
        for child in self.basis_jobs:
            energies.append(parse_output(
                Path(f'{child.name}.out'), method=self.kinbot_level['method']) / Hartree)
        if len(energies) == 2:
            value = two_point_cbs(
                *energies, upper_cardinal=_basis_cardinal(self.kinbot_level['basis'][1]),
                power=self.kinbot_level['extrapolation_power'])
        else:
            value = energies[0]
        return value * Hartree


class KinBotMultiSample(MultiSample):
    """Preserve ROTD's energy caching and corrections with a profiled adapter.

    Each serialized ROTD worker evaluates samples synchronously in its own
    process. Swap only its Molpro adapter for this call and restore it even
    when a native calculation fails; legacy inputs also use local QC storage.
    """

    def get_energies(self, calculator, face_id=0, flux_id=0):
        adapter = legacy_sample.Molpro
        try:
            legacy_sample.Molpro = MolproSample if 'kinbot_vrc' in calculator else LocalMolpro
            return super().get_energies(calculator, face_id, flux_id)
        finally:
            legacy_sample.Molpro = adapter
