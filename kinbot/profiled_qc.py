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
from pathlib import Path
import re

from kinbot.execution import atomic_json
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


def relaxed_scf_routes(payload, par):
    """Plan an explicit, idle-run migration; ordinary restarts remain strict.

    Previously accepted results may satisfy a looser SCF density threshold.
    No other profile change authorizes reuse. The caller must verify there
    are no live routed jobs and archive the original ledger before writing.
    """
    old = payload.get('profiles')
    digest = hashlib.sha256(json.dumps(old, sort_keys=True,
                                      separators=(',', ':')).encode()).hexdigest()
    if (payload.get('schema') != 1 or digest != payload.get('profile_sha256')
            or not isinstance(old, dict) or not isinstance(payload.get('jobs'), dict)):
        raise ValueError('Cannot migrate an invalid theory route ledger.')
    new_profiles = resolve_profiles(par)
    new = {key: _profile_payload(profile) for key, profile in new_profiles.items()}
    if old.keys() != new.keys():
        raise ValueError('Theory profile levels changed; retained results cannot be reused.')
    changes = {}
    for level in old:
        if old[level] == new[level]:
            continue
        before, after = deepcopy(old[level]), deepcopy(new[level])
        if before.get('calculator') != 'gaussian' or after.get('calculator') != 'gaussian':
            raise ValueError(f'{level}: only Gaussian SCF threshold relaxation can be migrated.')
        def split_scf(profile):
            from kinbot.gaussian_input import route_options
            kwargs = profile.get('calculator_kwargs', {})
            keys = [k for k in kwargs if k.casefold() == 'scf']
            if len(keys) != 1:
                raise ValueError(f'{level}: an explicit SCF convergence threshold is required.')
            options = route_options(kwargs.pop(keys[0]))
            convergence = [re.fullmatch(r'conver(?:gence)?\s*=\s*(\d+)', s, re.I)
                           for s in options]
            values = [int(m.group(1)) for m in convergence if m]
            if len(values) != 1:
                raise ValueError(f'{level}: one explicit SCF convergence threshold is required.')
            return values[0], sorted(s.casefold() for s, m in zip(options, convergence) if not m)
        previous, old_options = split_scf(before)
        requested, new_options = split_scf(after)
        if before != after or old_options != new_options or requested > previous:
            raise ValueError(f'{level}: theory changed beyond an SCF threshold relaxation.')
        changes[level] = {'previous_exponent': previous, 'requested_exponent': requested}
    result = deepcopy(payload)
    fingerprint = _fingerprint(new_profiles)
    for job, record in result['jobs'].items():
        if (not isinstance(job, str) or not isinstance(record, dict)
                or record.get('level') not in new
                or record.get('job_id') is not None and not isinstance(record['job_id'], str)
                or record.get('generation_tagged', False) not in (False, True)):
            raise ValueError(f'Invalid profiled job route for {job!r}.')
        # The archived ledger retains historical scheduler IDs. Only active
        # dispatcher jobs may remain in an explicitly migrated idle run.
        record.setdefault('source_profile_sha256', payload['profile_sha256'])
        record.update(job_id=None, generation_tagged=False)
    result.update(profile_sha256=fingerprint, profiles=new)
    if changes:
        result.setdefault('profile_migrations', []).append({
            'source_profile_sha256': payload['profile_sha256'],
            'target_profile_sha256': fingerprint, 'scf_changes': changes})
    return result


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
    # Backend-local marker: an explicit level profile also owns preliminary
    # constrained TS/ring searches; legacy AM1 shortcuts must not replace it.
    values['_profiled_backend'] = True
    values['high_level_method'] = profile.method
    values['high_level_basis'] = profile.basis
    values['calc_kwargs'] = dict(profile.calculator_kwargs)
    values['use_sella'] = profile.optimizer == 'sella'
    values['frequency_mode'] = profile.frequency_mode
    if calculator == 'fairchem':
        values['fc_model_path'] = profile.model_path
        values['fc_task_name'] = profile.task_name
        values['fc_device'] = profile.device
        values['qc_command'] = ''
    else:
        values['qc_command'] = (profile.command
                                or _DEFAULT_COMMANDS.get(calculator, ''))
    return values


def l1_runtime_parameters(qc, fallback):
    """Return the concrete settings used by an L1 worker.

    Most job writers are methods on the routed backend and therefore see its
    parameters directly.  Reaction-path and IRC writers live outside that
    class, so they must resolve the L1 backend explicitly instead of reading
    the profile container's legacy top-level placeholders.
    """
    backend = getattr(qc, 'l1', qc)
    return getattr(backend, 'par', fallback)


class ProfiledQuantumChemistry:
    """Route the existing QuantumChemistry surface by calculation level."""

    def __init__(self, par):
        from kinbot.qc import QuantumChemistry

        self.par = par
        route_file = par.get('profiled_route_file', str(_ROUTE_FILE))
        if not isinstance(route_file, str) or not route_file.strip():
            raise ValueError('profiled_route_file must be a nonempty path.')
        self.route_file = Path(route_file).expanduser()
        reuse_only = par.get('profiled_l1_reuse_only', 0)
        if reuse_only not in (0, 1, False, True):
            raise ValueError('profiled_l1_reuse_only must be 0 or 1.')
        self.l1_reuse_only = bool(reuse_only)
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
                if _level == 'l1' and self.l1_reuse_only:
                    status = self.backends[_level].check_qc(job)
                    if status == 0:
                        raise FileNotFoundError(
                            f'{job}: completed L1 result is required by the '
                            'reuse-only continuation; no replacement job was '
                            'submitted.')
                    return 0
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

    def _get_queue_setting(self, name):
        return getattr(self.l1, name)

    def _set_queue_setting(self, name, value):
        # Queue admission belongs to the workflow, rather than a theory
        # level.  Propagate temporary branch reservations to both routed
        # backends instead of shadowing the L1 attribute on this wrapper.
        for backend in getattr(self, 'backends', {}).values():
            setattr(backend, name, value)

    @property
    def queue_job_limit(self):
        return self._get_queue_setting('queue_job_limit')

    @queue_job_limit.setter
    def queue_job_limit(self, value):
        self._set_queue_setting('queue_job_limit', value)

    @property
    def queue_job_root(self):
        return self._get_queue_setting('queue_job_root')

    @queue_job_root.setter
    def queue_job_root(self, value):
        self._set_queue_setting('queue_job_root', value)

    @property
    def queue_job_root_recursive(self):
        return self._get_queue_setting('queue_job_root_recursive')

    @queue_job_root_recursive.setter
    def queue_job_root_recursive(self, value):
        self._set_queue_setting('queue_job_root_recursive', value)

    def _load_routes(self):
        self._profile_migrations = []
        if not self.route_file.exists():
            return {}
        try:
            payload = json.loads(self.route_file.read_text())
        except (OSError, ValueError) as exc:
            raise ValueError(f'Cannot read {self.route_file}.') from exc
        if (payload.get('schema') != 1
                or payload.get('profile_sha256') != self._fingerprint
                or not isinstance(payload.get('jobs'), dict)):
            raise ValueError('Theory profiles changed since jobs were submitted; '
                             f'archive or remove {self.route_file} before a new run.')
        routes = payload['jobs']
        self._profile_migrations = payload.get('profile_migrations', [])
        for job, record in routes.items():
            if (not isinstance(job, str) or not isinstance(record, dict)
                    or record.get('level') not in self.backends
                    or (record.get('job_id') is not None
                        and not isinstance(record.get('job_id'), str))
                    or record.get('generation_tagged', False)
                    not in (False, True)):
                raise ValueError(f'Invalid profiled job route for {job!r}.')
        return routes

    def _save_routes(self):
        payload = {'schema': 1, 'profile_sha256': self._fingerprint,
                   'profiles': {key: _profile_payload(value)
                                for key, value in self.profiles.items()},
                   'jobs': self._routes}
        if self._profile_migrations:
            payload['profile_migrations'] = self._profile_migrations
        self.route_file.parent.mkdir(parents=True, exist_ok=True)
        atomic_json(self.route_file, payload)

    def _record(self, job, level, job_id=None):
        previous = self._routes.get(job, {})
        record = {
            'level': level,
            'job_id': (str(job_id) if job_id is not None
                       else previous.get('job_id')),
            'generation_tagged': (
                bool(self.backends[level].job_generation_tagged.get(job))
                if job_id is not None
                else previous.get('generation_tagged', False)),
        }
        if 'source_profile_sha256' in previous:
            record['source_profile_sha256'] = (self._fingerprint if job_id is not None
                                               else previous['source_profile_sha256'])
        if job_id is not None:
            self.job_ids[job] = str(job_id)
        if previous != record:
            self._routes[job] = record
            self._save_routes()

    def _clear_job_id(self, job):
        self.job_ids.pop(job, None)
        record = self._routes.get(job)
        if record and record.get('job_id') is not None:
            backend = self.backends[record['level']]
            backend.job_ids.pop(job, None)
            backend.job_generation_tagged.pop(job, None)
            record['job_id'] = None
            record['generation_tagged'] = False
            self._save_routes()

    def _restore_job_ids(self):
        for job, record in self._routes.items():
            job_id = record.get('job_id')
            if job_id:
                backend = self.backends[record['level']]
                backend.job_ids[job] = job_id
                backend.job_generation_tagged[job] = bool(
                    record.get('generation_tagged', False))
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

    def read_imag_mode(self, job, natom):
        """Use the artifact contract of the backend that ran this job."""
        return self._backend_for(job).read_imag_mode(job, natom)

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
