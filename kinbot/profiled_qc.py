"""Persistent L1/L2 routing for KinBot's existing QC interface.

Reaction discovery and conformer sampling use L1. Accepted stationary points
and hindered rotors use L2.  The router delegates to the mature legacy job
writers and records every job's backend so result polling remains correct
after restarting KinBot.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path

from kinbot.theory import resolve_profiles


_ROUTE_FILE = Path('.kinbot_theory_jobs.json')
_QC_NAMES = {
    'fairchem': 'fc',
    'gaussian': 'gauss',
    'qchem': 'qchem',
    'nwchem': 'nwchem',
    'orca': 'orca',
    'nn_pes': 'nn_pes',
}
_DEFAULT_COMMANDS = {
    'gaussian': 'g16',
    'qchem': 'qchem',
    'nwchem': 'nwchem',
    'orca': 'orca',
}


def _profile_payload(profile):
    return {
        'calculator': profile.calculator,
        'method': profile.method,
        'basis': profile.basis,
        'command': profile.command,
        'calculator_kwargs': dict(profile.calculator_kwargs),
        'optimizer': profile.optimizer,
        'frequency_mode': profile.frequency_mode,
        'model_path': profile.model_path,
        'task_name': profile.task_name,
        'device': profile.device,
    }


def _fingerprint(profiles):
    payload = {level: _profile_payload(profile)
               for level, profile in sorted(profiles.items())}
    return hashlib.sha256(json.dumps(payload, sort_keys=True,
                                     separators=(',', ':')).encode()).hexdigest()


def _backend_parameters(par, profile):
    """Translate one resolved profile into the existing backend settings."""
    values = deepcopy(par)
    for key, empty in (
            ('profiled_theory', 0), ('theory_preset', ''), ('l1', ''),
            ('l2', ''), ('l1_profile', {}), ('l2_profile', {}),
            ('composite_method', ''), ('l3_overrides', {}),
            ('l3_resource_overrides', {})):
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
    values['use_sella'] = profile.optimizer == 'sella'
    if calculator == 'fairchem':
        values['fc_model_path'] = profile.model_path
        values['fc_task_name'] = profile.task_name
        values['fc_device'] = profile.device
        values['qc_command'] = ''
    else:
        values['qc_command'] = (profile.command
                                or _DEFAULT_COMMANDS.get(calculator, ''))
    return values


class ProfiledQuantumChemistry:
    """Route the existing QuantumChemistry surface by calculation level."""

    def __init__(self, par):
        from kinbot.qc import QuantumChemistry

        self.par = par
        self.profiles = resolve_profiles(par)
        if set(self.profiles) != {'l1', 'l2'}:
            raise ValueError('Profiled theory requires resolved L1 and L2 profiles.')
        self._fingerprint = _fingerprint(self.profiles)
        self.backends = {
            level: QuantumChemistry(_backend_parameters(par, profile))
            for level, profile in self.profiles.items()
        }
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
        self._routes = self._load_routes()
        for level, backend in self.backends.items():
            original = backend.submit_qc

            def tracked(job, nproc, singlejob=1, jobtype=None, *,
                        _level=level, _submit=original):
                result = _submit(job, nproc, singlejob=singlejob,
                                 jobtype=jobtype)
                job_id = self.backends[_level].job_ids.get(job)
                self._record(job, _level, job_id)
                return result

            backend.submit_qc = tracked
        self._restore_job_ids()

    def __getattr__(self, name):
        # Unclassified operations (reaction scans, IRC, VTS, conformers) are L1.
        return getattr(self.l1, name)

    def _load_routes(self):
        if not _ROUTE_FILE.exists():
            return {}
        try:
            payload = json.loads(_ROUTE_FILE.read_text())
        except (OSError, ValueError) as exc:
            raise ValueError(f'Cannot read {_ROUTE_FILE}.') from exc
        if (payload.get('schema') != 1
                or payload.get('profile_sha256') != self._fingerprint
                or not isinstance(payload.get('jobs'), dict)):
            raise ValueError('Theory profiles changed since jobs were submitted; '
                             f'archive or remove {_ROUTE_FILE} before a new run.')
        routes = payload['jobs']
        for job, record in routes.items():
            if (not isinstance(job, str) or not isinstance(record, dict)
                    or record.get('level') not in self.backends
                    or (record.get('job_id') is not None
                        and not isinstance(record.get('job_id'), str))):
                raise ValueError(f'Invalid profiled job route for {job!r}.')
        return routes

    def _save_routes(self):
        payload = {'schema': 1, 'profile_sha256': self._fingerprint,
                   'profiles': {key: _profile_payload(value)
                                for key, value in self.profiles.items()},
                   'jobs': self._routes}
        temporary = _ROUTE_FILE.with_suffix(_ROUTE_FILE.suffix + '.tmp')
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
        os.replace(temporary, _ROUTE_FILE)

    def _record(self, job, level, job_id=None):
        previous = self._routes.get(job, {})
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
        # This also lets a newly introduced router resume conventional KinBot
        # names produced by the same two-level input before a route file existed.
        if '/hir_' in job or job.startswith('hir/') or job.endswith('_high') \
                or '_high_freq_recovery_' in job:
            return 'l2'
        return 'l1'

    def _backend_for(self, job):
        return self.backends[self._level_for(job)]

    def get_qc_arguments(self, *args, **kwargs):
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
        # Hindered rotors are evaluated on the accepted L2 surface.
        return self.l2.qc_hir(*args, **kwargs)

    def qc_freq(self, species, source_job, high_level=0):
        level = 'l2' if high_level else self._level_for(source_job)
        result = self.backends[level].qc_freq(
            species, source_job, high_level=int(level == 'l2'))
        self._record(result, level,
                     self.backends[level].job_ids.get(result))
        return result

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

    def hessian_is_massweighted(self):
        level = 'l2' if self.par.get('high_level') else 'l1'
        return self.backends[level].hessian_is_massweighted()

    def invalidate_qc(self, job):
        result = self._backend_for(job).invalidate_qc(job)
        self._clear_job_id(job)
        return result

    def publish_result(self, source, target):
        level = self._level_for(source.name)
        result = self.backends[level].publish_result(source, target)
        self._record(target, level)
        self._clear_job_id(target)
        return result
