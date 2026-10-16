"""Persistent L0/L1/L2 routing for KinBot's existing QC interface.

Preliminary searches use L0; discovery and ordinary conformers use L1.
Refinement uses L2 when enabled, and rotors follow that selected level.  The router delegates to the mature legacy job
writers and records every job's backend so result polling remains correct
after restarting KinBot.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import os
import re
from pathlib import Path

from kinbot.theory import resolve_profiles


class MissingL1Result(FileNotFoundError):
    """A reuse-only continuation needs an L1 result that cannot be reused."""


_ROUTE_FILE = Path('.kinbot_theory_jobs.json')
_QC_NAMES = {
    'fairchem': 'fc',
    'gaussian': 'gauss',
    'qchem': 'qchem',
    'nwchem': 'nwchem',
    'orca': 'orca',
}
_DEFAULT_COMMANDS = {
    'gaussian': 'g16',
    'qchem': 'qchem',
    'nwchem': 'nwchem',
    'orca': 'orca',
}


def _rotor_job(job):
    return job.startswith('hir/') or '/hir_' in job


def _route_file(par):
    filename = par.get('profiled_route_file', str(_ROUTE_FILE))
    if not isinstance(filename, str) or not filename.strip():
        raise ValueError('profiled_route_file must be a nonempty path.')
    return Path(filename).expanduser()


def validate_theory_restart(par, directory='.', *, profiles=None):
    """Check recorded theory without opening calculators or changing results."""
    from kinbot.run_format import MARKER, _existing_results, ensure_current_run

    directory = Path(directory)
    marker = directory / MARKER
    recorded = ensure_current_run(directory) if marker.exists() else {}
    if profiles is None:
        profiles = resolve_profiles(par)
    route_file = directory / _route_file(par)
    if not profiles:
        paths = {directory / _ROUTE_FILE, route_file}
        if recorded.get('profiled_theory') is True:
            paths.add(marker)
        for path in paths:
            if path.exists():
                raise ValueError(f'Theory profiles were removed, but {path} still '
                                 'records the previous use of profiles. Restore '
                                 'the profiles or start in a fresh directory; keep '
                                 'the theory record with the existing results.')
        return {}
    if not route_file.exists():
        existing = _existing_results(directory, ignore_startup=True)
        if existing is not None or recorded.get('profiled_theory') is True:
            raise ValueError(f'Cannot verify the theory used for {existing or directory}: '
                             f'{route_file} is missing. Start the profiled '
                             'calculation in a fresh directory; keep the '
                             'existing calculation files unchanged.')
        return {}
    try:
        payload = json.loads(route_file.read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f'Cannot read {route_file}.') from exc
    if (not isinstance(payload, dict) or payload.get('schema') != 2
            or payload.get('profile_sha256') != _fingerprint(profiles, par)
            or not isinstance(payload.get('jobs'), dict)):
        raise ValueError('Theory profiles or calculator settings changed, or the '
                         'saved theory record lacks required settings; '
                         'start a new calculation in a fresh directory. '
                         f'Keep {route_file} with the existing results; '
                         'removing it does not remove calculations made with the old profiles.')
    routes = payload['jobs']
    for job, record in routes.items():
        if (not isinstance(job, str) or not isinstance(record, dict)
                or record.get('level') not in profiles
                or (record.get('job_id') is not None
                    and not isinstance(record.get('job_id'), str))):
            raise ValueError(f'Invalid profiled job route for {job!r}.')
        # Rotor job names do not distinguish L1 from L2. Check before either
        # calculation workers or PES can read scans from the previous level.
        if _rotor_job(job):
            level = 'l2' if par['high_level'] else 'l1'
            if record['level'] != level:
                raise ValueError(f'{job}: saved rotor scan uses '
                                 f'{record["level"].upper()}, but high_level '
                                 f'requests {level.upper()}. Run the new level '
                                 'in a fresh directory; existing scans were preserved.')
    return routes


def _profile_payload(profile, par):
    payload = {
        'calculator': profile.calculator,
        'method': profile.method,
        'basis': profile.basis,
        'calculator_kwargs': dict(profile.calculator_kwargs),
        'optimizer': profile.optimizer,
        'frequency_mode': profile.frequency_mode,
        'model_path': profile.model_path,
        'task_name': profile.task_name,
    }
    # These method choices still come from the ordinary input, even when
    # calculators use profiles. Omitting them can relabel a saved energy.
    keys = []
    if profile.calculator == 'orca':
        keys += ['orcasimpleinput', 'orcablocks']
    if profile.calculator == 'gaussian':
        keys += ['integral', 'guessmix']
    if profile.calculator == 'gaussian' and profile.name == 'l2':
        keys += ['vrc_tst_scan_method', 'vrc_tst_scan_basis']
    if profile.calculator in ('gaussian', 'qchem'):
        keys += ['scan_method', 'scan_basis',
                 'barrierless_saddle_method', 'barrierless_saddle_basis',
                 'barrierless_saddle_method_high', 'barrierless_saddle_basis_high']
    if keys:
        payload['calculation_settings'] = {key: par[key] for key in keys}
    return payload


def _fingerprint(profiles, par):
    payload = {level: _profile_payload(profile, par)
               for level, profile in sorted(profiles.items())}
    return hashlib.sha256(json.dumps(payload, sort_keys=True,
                                     separators=(',', ':')).encode()).hexdigest()


def _backend_parameters(par, profile):
    """Translate one resolved profile into the existing backend settings."""
    values = deepcopy(par)
    for key, empty in (
            ('profiled_theory', 0), ('l0', ''),
            ('l0_profile', {}), ('l1', ''),
            ('l2', ''), ('l1_profile', {}), ('l2_profile', {})):
        values[key] = empty
    calculator = profile.calculator.lower()
    try:
        values['qc'] = _QC_NAMES[calculator]
    except KeyError as exc:
        raise ValueError(f'Profiled QC backend {calculator!r} is not supported '
                         'by the KinBot reaction workflow.') from exc
    values['method'] = profile.method
    values['basis'] = profile.basis
    values['high_level_method'] = profile.method
    values['high_level_basis'] = profile.basis
    values['calc_kwargs'] = dict(profile.calculator_kwargs)
    values['use_sella'] = profile.optimizer in ('sella', 'bfgs')
    values['frequency_mode'] = profile.frequency_mode
    if calculator == 'fairchem':
        values['fc_model_path'] = profile.model_path
        values['fc_task_name'] = profile.task_name
        values['fc_device'] = profile.device
        values['qc_command'] = ''
    else:
        values['qc_command'] = (profile.command
                                or _DEFAULT_COMMANDS.get(calculator, ''))
    # Generic Sella templates instantiate calculators directly. Supply the
    # executable to the instance so model and command follow the same profile.
    if calculator == 'gaussian':
        values['calc_kwargs']['command'] = values['qc_command'] + ' < PREFIX.com > PREFIX.log'
    elif calculator == 'qchem':
        values['calc_kwargs']['command'] = values['qc_command']
    return values


class ProfiledQuantumChemistry:
    """Route the existing QuantumChemistry surface by calculation level."""

    def __init__(self, par):
        from kinbot.qc import QuantumChemistry
        from kinbot.run_format import ensure_current_run, record_profiled_run

        self.par = par
        self.route_file = _route_file(par)
        reuse_only = par.get('profiled_l1_reuse_only', 0)
        if reuse_only not in (0, 1, False, True):
            raise ValueError('profiled_l1_reuse_only must be 0 or 1.')
        self.l1_reuse_only = bool(reuse_only)
        self.profiles = resolve_profiles(par)
        if set(self.profiles) != {'l0', 'l1', 'l2'}:
            raise ValueError('Profiled theory requires resolved L0, L1 and L2 profiles.')
        self._fingerprint = _fingerprint(self.profiles, par)
        self._routes = validate_theory_restart(par, profiles=self.profiles)
        ensure_current_run(create=True)
        # Save the methods before marking profile use or opening calculators,
        # so a failed initialization still has the record needed for a restart.
        execution = {key: {'command': value.command, 'device': value.device}
                     for key, value in self.profiles.items()}
        if (not self.route_file.exists()
                or json.loads(self.route_file.read_text()).get('execution') != execution):
            self._save_routes()
        record_profiled_run()
        self.backends = {
            level: QuantumChemistry(_backend_parameters(par, profile), profile_backend=True)
            for level, profile in self.profiles.items()
        }
        # L3 and AIE retain their ordinary global settings.
        self.legacy = QuantumChemistry(par, profile_backend=True)
        for level, backend in self.backends.items():
            backend.theory_level = level
            backend.optimizer = self.profiles[level].optimizer
            if backend.use_sella or backend.qc == 'fc':
                backend.par['sella_kwargs'].setdefault('internal', True)
                if backend.par['calcall_ts']:
                    backend.par['sella_kwargs'].setdefault('gamma', 0.)
        self.l0 = self.backends['l0']
        self.l1 = self.backends['l1']
        self.l2 = self.backends['l2']
        # Existing reaction code inspects these attributes directly. Reaction
        # searches are L1 by definition; explicit high-level calls are routed
        # by the methods below.
        self.qc = self.l1.qc
        self.use_sella = self.l1.use_sella
        self.ppn = self.l1.ppn
        self.db = self.l1.db
        self.job_ids = {}
        for level, backend in self.backends.items():
            original = backend.submit_qc
            original_check = backend.check_qc

            def check(job, *, _level=level, _check=original_check):
                # An existing rotor name does not identify its calculation
                # level. Refuse it before a backend can read or ingest results.
                if _rotor_job(job) and job not in self._routes:
                    saved = self.db.count(name=job) or any(
                        Path(job + suffix).exists()
                        for suffix in ('.py', '.pkl', '.log', '.out', '_sella.log'))
                    if saved:
                        raise ValueError(f'{job}: saved rotor has no recorded theory level. '
                                         'Use a fresh directory; existing files were preserved.')
                if job in self._routes and self._routes[job]['level'] != _level:
                    raise ValueError(f'{job}: requested {_level.upper()} differs from '
                                     'its recorded theory level.')
                return _check(job)

            backend.check_qc = check

            def tracked(job, nproc, singlejob=1, jobtype=None, *,
                        _level=level, _submit=original):
                if _level == 'l1' and self.l1_reuse_only:
                    status = self.backends[_level].check_qc(job)
                    if status == 0:
                        raise MissingL1Result(
                            f'{job}: completed L1 result is required by the '
                            'reuse-only continuation; no replacement job was '
                            'submitted.')
                    return 0
                self._record(job, _level)
                if re.search(r'_IRC_[FR]$', job):
                    self._record(job + '_prod', _level)
                result = _submit(job, nproc, singlejob=singlejob,
                                 jobtype=jobtype)
                job_id = self.backends[_level].job_ids.get(job)
                self._record(job, _level, job_id)
                return result

            backend.submit_qc = tracked
        self._restore_job_ids()

    def __getattr__(self, name):
        # Missing construction or serialization attributes must raise
        # AttributeError, not recursively request l1 before it exists.
        backend = self.__dict__.get('l1')
        if backend is None or name.startswith('__'):
            raise AttributeError(name)
        return getattr(backend, name)

    def _save_routes(self):
        payload = {'schema': 2, 'profile_sha256': self._fingerprint,
                   'profiles': {key: _profile_payload(value, self.par)
                                for key, value in self.profiles.items()},
                   'execution': {key: {'command': value.command, 'device': value.device}
                                 for key, value in self.profiles.items()},
                   'jobs': self._routes}
        self.route_file.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.route_file.with_suffix(self.route_file.suffix + '.tmp')
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
        os.replace(temporary, self.route_file)

    def _record(self, job, level, job_id=None):
        previous = self._routes.get(job, {})
        if previous and previous['level'] != level:
            raise ValueError(f'{job}: cannot replace recorded {previous["level"].upper()} '
                             f'with {level.upper()}.')
        record = {'level': level,
                  'job_id': str(job_id) if job_id is not None
                  else previous.get('job_id')}
        if job_id is not None:
            self.job_ids[job] = str(job_id)
        if previous != record:
            self._routes[job] = record
            self._save_routes()

    def _clear_job_id(self, job):
        self.job_ids.pop(job, None)
        record = self._routes.get(job)
        if record and record.get('job_id') is not None:
            record['job_id'] = None
            self._save_routes()

    def _restore_job_ids(self):
        for job, record in self._routes.items():
            job_id = record.get('job_id')
            if job_id:
                self.backends[record['level']].job_ids[job] = job_id
                self.job_ids[job] = job_id

    def _level_for(self, job):
        record = self._routes.get(job)
        if record:
            return record['level']
        # New ordinary names identify their level before submission. Rotor
        # names do not, so results require the record saved before launch.
        if '_semi_emp_' in job or (job.startswith('conf/')
                and re.search(r'_r\d+_\d+$', job)):
            return 'l0'
        if _rotor_job(job):
            raise ValueError(f'{job}: rotor result has no recorded theory level; '
                             'do not infer it from high_level.')
        if job.endswith('_high') \
                or '_high_freq_recovery_' in job:
            return 'l2'
        return 'l1'

    def for_level(self, level):
        return self.backends[level]

    def for_job(self, job):
        return self.backends[self._level_for(job)]

    _backend_for = for_job

    def qc_ring_conf(self, *args, **kwargs):
        if self.l0.qc == 'qchem':
            raise ValueError('Q-Chem cannot prepare ring conformers in KinBot. '
                             'Select Gaussian or FairChem in l0_profile for '
                             'ring preparation; Q-Chem can still be used at L1/L2.')
        return self.l0.qc_ring_conf(*args, **kwargs)

    def qc_conf(self, species, geom, index, semi_emp=0):
        backend = self.l0 if semi_emp else self.l1
        return backend.qc_conf(species, geom, index, semi_emp=semi_emp)

    def qc_aie(self, *args, **kwargs):
        return self.legacy.qc_aie(*args, **kwargs)

    def get_qc_arguments(self, *args, **kwargs):
        if kwargs.get('L3') or (len(args) > 14 and args[14]):
            return self.legacy.get_qc_arguments(*args, **kwargs)
        high_level = (args[9] if len(args) > 9 else kwargs.get('high_level'))
        level = 'l2' if high_level else 'l1'
        return self.backends[level].get_qc_arguments(*args, **kwargs)

    def submit_qc(self, job, nproc, singlejob=1, jobtype=None):
        return self.l1.submit_qc(job, nproc, singlejob=singlejob,
                                 jobtype=jobtype)

    def qc_opt(self, *args, **kwargs):
        high_level = (args[2] if len(args) > 2 else kwargs.get('high_level'))
        level = 'l2' if high_level else 'l1'
        return self.backends[level].qc_opt(*args, **kwargs)

    def qc_opt_ts(self, *args, **kwargs):
        high_level = (args[2] if len(args) > 2 else kwargs.get('high_level'))
        level = 'l2' if high_level else 'l1'
        return self.backends[level].qc_opt_ts(*args, **kwargs)

    def qc_hir(self, *args, **kwargs):
        level = 'l2' if self.par['high_level'] else 'l1'
        return self.backends[level].qc_hir(*args, **kwargs)

    def qc_vts_frag(self, frag):
        """Run VRC fragment preparation on its Gaussian theory context."""
        if self.l2.qc != 'gauss':
            raise ValueError('The current VRC implementation requires a '
                             'Gaussian L2 backend.')
        job = self.l2.qc_vts_frag(frag)
        self._record(job, 'l2', self.l2.job_ids.get(job))
        return job

    def qc_vts(self, *args, **kwargs):
        """Run the constrained VRC geometry point with Gaussian, not UMA."""
        if self.l2.qc != 'gauss':
            raise ValueError('The current VRC implementation requires a '
                             'Gaussian L2 backend.')
        job = self.l2.qc_vts(*args, **kwargs)
        self._record(job, 'l2', self.l2.job_ids.get(job))
        return job

    def qc_freq(self, species, source_job, high_level=0):
        level = self._level_for(source_job)
        result = self.backends[level].qc_freq(
            species, source_job, high_level=int(level == 'l2'))
        self._record(result, level,
                     self.backends[level].job_ids.get(result))
        return result

    def read_imag_mode(self, job, natom):
        return self.for_job(job).read_imag_mode(job, natom)

    def check_qc(self, job):
        return self._backend_for(job).check_qc(job)

    def is_in_database(self, job):
        return self._backend_for(job).is_in_database(job)

    def get_qc_geom(self, job, *args, **kwargs):
        return self._backend_for(job).get_qc_geom(job, *args, **kwargs)

    def get_qc_freq(self, job, *args, **kwargs):
        return self._backend_for(job).get_qc_freq(job, *args, **kwargs)

    def get_qc_energy(self, job, *args, **kwargs):
        return self._backend_for(job).get_qc_energy(job, *args, **kwargs)

    def get_qc_zpe(self, job, *args, **kwargs):
        return self._backend_for(job).get_qc_zpe(job, *args, **kwargs)

    def read_qc_hess(self, job, *args, **kwargs):
        return self._backend_for(job).read_qc_hess(job, *args, **kwargs)

    def hessian_is_massweighted(self, job=None):
        level = self._level_for(job) if job else ('l2' if self.par.get('high_level') else 'l1')
        return self.backends[level].hessian_is_massweighted()

    def invalidate_qc(self, job):
        result = self._backend_for(job).invalidate_qc(job)
        self._clear_job_id(job)
        return result

    def publish_result(self, source, target, chemical_context=None):
        level = self._level_for(source.name)
        self._record(target, level)
        result = self.backends[level].publish_result(
            source, target, chemical_context=chemical_context)
        self._clear_job_id(target)
        return result
