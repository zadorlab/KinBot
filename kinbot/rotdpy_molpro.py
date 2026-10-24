"""Molpro sampling bridge for the pinned ROTD_py worker interface.

ROTD_py owns sampling, scheduling, restart files, and flux corrections. This
bridge supplies its Molpro energy evaluations with KinBot's reviewed orbital
model and CBS algebra, without modifying the pinned external package.
"""

from copy import deepcopy
from pathlib import Path
import json
import re
import shutil
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
from kinbot.execution import atomic_json, file_sha256


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
            header = (
                f'***,KinBot ROTD sampling\nmemory,{self.memory_in_MW},m\n'
                'angstrom\norient,noorient\nsymmetry,nosym\n'
                f'set,charge={model["charge"]}\nset,spin={mult - 1}\n'
                f'geometry={{\n{len(atoms)}\n{name}\n{rows}\n}}\n'
                + f'basis={basis}\n'
                + f'gthresh,energy={scf},coeff={correlation},civec={correlation}\n')
            child.kinbot_retry = dict(header=header, model=deepcopy(model['state_models_by_basis'][basis]),
                method=level['method'], electrons=electron_count, mult=mult,
                shift=model['shift'], ipea=model['ipea'], scf=model['scf_convergence'])
            Path(f'{name}.inp').write_text(initialize_molpro_restart(header
                + method + energy_assignment(selected_energy) + f'MYENERGY={selected_energy}\n'))
            self.basis_jobs.append(child)

    def run(self):
        for child in self.basis_jobs:
            self._run_with_solver_recovery(child)

    def _run_with_solver_recovery(self, child):
        """Retry this exact geometry, never substitute a fresh random sample."""
        from kinbot.anl.molpro_solver_recovery import solver_options
        path = Path(child.name+'.inp')
        original = file_sha256(path)
        archive = Path('.kinbot_sample_attempts')/child.name/original
        journal = archive/'recovery.json'
        record = json.loads(journal.read_text()) if journal.exists() else dict(attempts=[])
        if record.get('status') == 'exhausted':
            raise RuntimeError(f'Molpro sample {child.name}: solver retries exhausted for this geometry.')
        model = deepcopy(record.get('next_model', child.kinbot_retry['model']))
        if record.get('next_model'):
            self._retry_input(child, model, len(record['attempts']))
        while True:
            try:
                self._run_basis_job(child)
            except RuntimeError as error:
                # Only a native nonzero return is eligible, not parse, setup,
                # scratch, subprocess or validation errors.
                if (not getattr(error, 'native_output_current', False)
                        or str(error) != f'Molpro sample {child.name} exited with status {getattr(error, "native_returncode", None)}.'):
                    raise
                output = Path(child.name+'.out')
                text = output.read_text(errors='replace') if output.exists() else ''
                next_model = solver_options(model, text, len(record['attempts']))
                if not next_model and not record['attempts']:
                    raise
                dest = archive/str(len(record['attempts']))
                dest.mkdir(parents=True, exist_ok=True)
                hashes = {}
                sources = list(Path.cwd().glob(child.name+'.*'))
                restart = re.search(r'(?m)^file,2,([^,\n]+)', path.read_text())
                if restart and Path(restart[1]).name == restart[1]:
                    sources.append(Path.cwd()/restart[1])
                for source in sources:
                    if source.is_file():
                        digest = file_sha256(source)
                        saved = dest/source.name
                        if saved.exists() and file_sha256(saved) != digest:
                            raise RuntimeError('Sample recovery archive changed; inspect before retry.')
                        if not saved.exists():
                            shutil.copy2(source, saved)
                        hashes[source.name] = digest
                record['attempts'].append(dict(files=hashes, model=model,
                    returncode=error.native_returncode, archive=str(dest.resolve())))
                record.update(status='retry_ready' if next_model else 'exhausted', next_model=next_model)
                atomic_json(journal, record)
                if next_model is None:
                    raise
                model = next_model
                self._retry_input(child, model, len(record['attempts']))
            else:
                if record['attempts']:
                    record.update(status='native_executed', output_sha256=file_sha256(Path(child.name+'.out')))
                    atomic_json(journal, record)
                return

    @staticmethod
    def _retry_input(child, model, attempt):
        settings = child.kinbot_retry
        selection = stage_vrc_state_model(model, Path.cwd(), child.name)
        method, energy = _vrc_multireference_method(settings['method'], settings['electrons'],
            spin=settings['mult']-1, shift=settings['shift'], ipea=settings['ipea'], state_model=selection)
        header = settings['header']
        if attempt >= 2:
            scf = fortran_threshold(min(settings['scf'], 1.e-10), 'VRC SCF convergence')
            header = re.sub(r'gthresh,energy=[^,\n]+', 'gthresh,energy='+scf, header)
        Path(child.name+'.inp').write_text(initialize_molpro_restart(header+method
            + energy_assignment(energy)+f'MYENERGY={energy}\n'))

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
                # A launcher failure may leave an old output in the shared
                # worker directory. Only this fresh scratch attempt is evidence.
                output = local/f'{child.name}.out'
                has_output = output.is_file() and output.stat().st_size > 0
            from kinbot.molpro_diagnostics import write_native_report
            if has_output:
                write_native_report(Path(f'{child.name}.out'))
            if result.returncode:
                error = RuntimeError(f'Molpro sample {child.name} exited with status {result.returncode}.')
                error.native_returncode = result.returncode
                error.native_output_current = has_output
                raise error
            if not has_output:
                raise RuntimeError(f'Molpro sample {child.name} produced no native output in its scratch attempt.')
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
