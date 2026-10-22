from kinbot.species_routing import routing_name
import getpass
import os
import sys
import fcntl
import subprocess
import logging
import re
import shlex
import time
import json
import threading
from datetime import datetime
from pathlib import Path
import copy
import pickle
from shutil import copyfile

import numpy as np
from ase.db import connect
from kinbot.run_format import ensure_current_run
from kinbot.stereo_identity import optical_scope
from ase import Atoms

from kinbot.utils import plain_symbols, plain_geometry
from kinbot import kb_path
from kinbot import constants
from kinbot import geometry
from kinbot import exceptions
from kinbot.execution import scheduler_unavailable as _scheduler_unavailable
from kinbot.execution import slurm_job_record, SLURM_TERMINAL_STATES as _SLURM_TERMINAL_STATES
from kinbot.theory import profiled_requested
from kinbot.minimum_recovery import worker_sella_kwargs

logger = logging.getLogger('KinBot')


_SLURM_PARTITION = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')
_SLURM_TIME = re.compile(r'(?:(\d+)-)?(\d{1,2}):(\d{2}):(\d{2})\Z')
_SLURM_WALLTIME_CACHE = {}
_SLURM_NODE_STATE = '.kinbot_slurm_nodes.json'
_SCHEDULER_QUERY_TIMEOUT = 15
_SCHEDULER_QUERY_FAILURE_TIMEOUT = 60
_SCHEDULER_SUBMIT_TIMEOUT = 30
_SUBMISSION_LOCK = threading.RLock()
_QC_SITE_BACKENDS = {
    'gauss': ('gaussian', 'g16'),
    'nwchem': ('nwchem', 'nwchem'),
    'orca': ('orca', 'orca'),
    'qchem': ('qchem', 'qchem'),
}

def exclusive_slurm_script(script, enabled=True):
    """Enforce whole-node scheduling, including user-supplied headers.

    Reject conflicting sharing requests rather than silently weakening them.
    Insert before shell commands, where Slurm actually reads directives.
    """
    if enabled not in (True, False, 0, 1):
        raise ValueError('queue_exclusive must be 0 or 1')
    if not enabled:
        return script
    for line in script.splitlines():
        if re.match(r'\s*#SBATCH\s+', line):
            options = shlex.split(line.strip()[len('#SBATCH'):])
            if any(value in ('-s', '--share', '--oversubscribe') or
                   value.startswith('--exclusive=') for value in options):
                raise ValueError('Exclusive QC scheduling conflicts with template sharing options')
    lines = script.splitlines(keepends=True)
    # A first directive also works with headers having no scheduler options.
    lines.insert(1 if lines and lines[0].startswith('#!') else 0,
                 '#SBATCH --exclusive\n')
    return ''.join(lines)


def _worker_exit_status(path, expected_job_id=None):
    """Read one atomic worker-completion marker.

    New scheduler wrappers write ``<scheduler job id> <exit code>`` so a
    logical KinBot job name can be reused without accepting the preceding
    submission's files. One-field markers remain readable for calculations
    created by older KinBot revisions.

    Return ``(exit_code, generation_is_explicit)``. Return ``None`` when the
    marker is absent or explicitly belongs to a different scheduler job.
    Raise ``ValueError`` for malformed content.
    """
    try:
        fields = open(path).read().strip().split()
    except FileNotFoundError:
        return None
    if len(fields) == 1:
        return int(fields[0]), False
    if len(fields) != 2:
        raise ValueError('worker status must contain an exit code, optionally '
                         'preceded by its scheduler job id')
    marker_job_id, raw_code = fields
    marker_job_id = marker_job_id.split('.', 1)[0]
    if (expected_job_id is not None
            and marker_job_id != str(expected_job_id).split('.', 1)[0]):
        return None
    return int(raw_code), True


def _pbs_job_active(job_id):
    """Return PBS activity, or ``None`` when scheduler status is unavailable."""
    try:
        result = subprocess.run(
            ['qstat', '-f', str(job_id)], capture_output=True, text=True,
            check=False, timeout=_SCHEDULER_QUERY_TIMEOUT)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return result.returncode == 0


def _vrc_scan_requires_frozen_handoff(job):
    """Return whether one worker must publish a paired frozen VRC point."""
    return bool(re.fullmatch(
        r'vrctst/.+_vts_pt(?:\d+|_asymptote)', str(job)))


def _slurm_node_names(values):
    """Validate explicit single-node names used in ``--exclude``."""
    if values in (None, ''):
        return set()
    if not isinstance(values, (list, tuple, set)):
        raise ValueError('slurm_exclude_nodes must be a list of node names.')
    names = {str(value).strip() for value in values}
    if '' in names or any(not _SLURM_PARTITION.fullmatch(name)
                          for name in names):
        raise ValueError('slurm_exclude_nodes contains an invalid node name.')
    return names


def _load_slurm_node_state(parameters=None):
    """Read run-local node exclusions, preserving explicit site settings."""
    names = _slurm_node_names((parameters or {}).get('slurm_exclude_nodes', []))
    inherited = [item.strip() for item in
                 os.environ.get('SBATCH_EXCLUDE', '').split(',')
                 if item.strip()]
    names.update(_slurm_node_names(inherited))
    try:
        payload = json.loads(open(_SLURM_NODE_STATE).read())
    except FileNotFoundError:
        return names
    except (OSError, ValueError, TypeError) as error:
        logger.warning('Ignoring invalid %s: %s', _SLURM_NODE_STATE, error)
        return names
    recorded = payload.get('excluded_nodes', []) if isinstance(payload, dict) else []
    try:
        return names | _slurm_node_names(recorded)
    except ValueError as error:
        logger.warning('Ignoring invalid node exclusions in %s: %s',
                       _SLURM_NODE_STATE, error)
        return names


def _slurm_exclude_directive(parameters=None, nodes=None):
    """Render explicitly configured and run-local node exclusions."""
    selected = (_load_slurm_node_state(parameters)
                if nodes is None else _slurm_node_names(nodes))
    return ('#SBATCH --exclude=' + ','.join(sorted(selected))
            if selected else '')


def _export_slurm_excluded_nodes(nodes):
    """Propagate run exclusions to nested Slurm submitters.

    Automatic ANL dispatchers and ROTD_py launch their own ``sbatch``
    processes after the ordinary KinBot search.  Slurm's ``SBATCH_EXCLUDE``
    environment setting lets those subprocesses inherit the same site
    recovery decision without adding a cluster-specific node to molecular
    input or to generated workflow schemas.
    """
    selected = sorted(_slurm_node_names(nodes))
    if not selected:
        return
    inherited = [item.strip() for item in
                 os.environ.get('SBATCH_EXCLUDE', '').split(',')
                 if item.strip()]
    for node in selected:
        if node not in inherited:
            inherited.append(node)
    os.environ['SBATCH_EXCLUDE'] = ','.join(inherited)


def _save_slurm_node_failure(owner, node, job, job_id, record):
    """Exclude a node after Slurm kills a batch before its script can run."""
    if not node or node in ('(null)', 'None') or not _SLURM_PARTITION.fullmatch(node):
        return
    names = set(getattr(owner, '_slurm_excluded_nodes', set()))
    names.add(node)
    from kinbot.run_policy import exclude_failed_node
    exclude_failed_node(Path.cwd(), node)
    owner._slurm_excluded_nodes = names
    _export_slurm_excluded_nodes(names)
    payload = {
        'schema': 1,
        'excluded_nodes': sorted(names),
        'last_failure': {
            'job': job,
            'job_id': str(job_id),
            'node': node,
            'state': record.get('JobState'),
            'reason': record.get('Reason'),
            'exit_code': record.get('ExitCode'),
        },
    }
    temporary = f'{_SLURM_NODE_STATE}.tmp.{os.getpid()}'
    try:
        with open(temporary, 'w') as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write('\n')
        os.replace(temporary, _SLURM_NODE_STATE)
    except OSError as error:
        logger.warning('Could not record failed Slurm node %s: %s', node, error)
        try:
            os.unlink(temporary)
        except OSError:
            pass


def _slurm_job_record(job_id):
    """Return the retained one-line Slurm record for a completed job."""
    return slurm_job_record(job_id, timeout=_SCHEDULER_QUERY_TIMEOUT)


def _slurm_job_active(job_id):
    """Return whether a job is in Slurm, or ``None`` if Slurm is unavailable.

    A targeted, headerless query keeps the response small on busy clusters.
    Most importantly, the timeout prevents a stalled Slurm controller from
    freezing the KinBot controller. An unknown result is deliberately
    distinct from an absent job; callers must not interpret a scheduler
    outage as successful completion.
    """
    if job_id in (None, '', '-1'):
        return False
    try:
        result = subprocess.run(
            ['squeue', '-h', '-j', str(job_id), '-o', '%i'],
            capture_output=True, text=True, check=False,
            timeout=_SCHEDULER_QUERY_TIMEOUT)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if result.returncode:
        if re.search(r'\binvalid job id specified\b', result.stderr,
                     re.IGNORECASE):
            return False
        return None
    return any(line.split()[0] == str(job_id)
               for line in result.stdout.splitlines() if line.split())


def _slurm_seconds(value):
    """Return seconds for a finite Slurm maximum-time display."""
    if value.lower() in ('infinite', 'unlimited'):
        return None
    match = _SLURM_TIME.fullmatch(value)
    if not match:
        # ``sinfo %l`` may omit hours for short limits.
        fields = value.split(':')
        if len(fields) == 2 and all(field.isdigit() for field in fields):
            return int(fields[0]) * 60 + int(fields[1])
        return None
    days, hours, minutes, seconds = (
        int(field or 0) for field in match.groups())
    return ((days * 24 + hours) * 60 + minutes) * 60 + seconds


def _slurm_time(seconds):
    days, remainder = divmod(int(seconds), 24 * 60 * 60)
    hours, remainder = divmod(remainder, 60 * 60)
    minutes, seconds = divmod(remainder, 60)
    clock = f'{hours:02d}:{minutes:02d}:{seconds:02d}'
    return f'{days}-{clock}' if days else clock


def _ordinary_slurm_walltime(par, queue_name=None):
    """Resolve the explicit wall time used by bundled Slurm templates."""
    requested = str(par.get('queue_walltime', 'partition_max')).strip()
    if not requested or requested.lower() in ('default', 'none'):
        return ''
    if requested != 'partition_max':
        if _slurm_seconds(requested) is None:
            raise ValueError('queue_walltime must be partition_max or a finite '
                             'Slurm D-HH:MM:SS/HH:MM:SS value.')
        return f'#SBATCH --time={requested}'
    partition = str(queue_name or par.get('queue_name', '')).strip()
    if not _SLURM_PARTITION.fullmatch(partition):
        raise ValueError('queue_name is not a valid Slurm partition name.')
    cached = _SLURM_WALLTIME_CACHE.get(partition)
    if cached and time.monotonic() - cached[0] < 60:
        return cached[1]
    try:
        result = subprocess.run(
            ['sinfo', '-h', '-p', partition, '-o', '%l'],
            capture_output=True, text=True, check=False, timeout=15)
    except (OSError, subprocess.TimeoutExpired):
        directive = ''
    else:
        displays = [line.strip() for line in result.stdout.splitlines()
                    if line.strip()]
        limits = [_slurm_seconds(value) for value in displays]
        finite = [limit for limit in limits if limit is not None]
        valid = all(limit is not None or value.lower() in ('infinite', 'unlimited')
                    for value, limit in zip(displays, limits))
        if result.returncode == 0 and valid and finite:
            directive = f'#SBATCH --time={_slurm_time(min(finite))}'
        elif (result.returncode == 0 and displays
              and all(value.lower() in ('infinite', 'unlimited')
                      for value in displays)):
            directive = '#SBATCH --time=0'
        else:
            directive = ''
    if not directive:
        from kinbot.anl.dispatch import SchedulerUnavailableError
        raise SchedulerUnavailableError('Could not discover the maximum walltime for '
                                        f'partition {partition}; no job was submitted.')
    _SLURM_WALLTIME_CACHE[partition] = (time.monotonic(), directive)
    return directive


class QuantumChemistry:
    '''
    This class provides the link between KinBot and the qc code
    It choses the write options, submits the jobs and checks
    the jobs for success or failure
    '''

    def __new__(cls, par=None):
        """Select the opt-in L1/L2 router without changing legacy inputs."""
        if cls is QuantumChemistry and par is not None and profiled_requested(par):
            from kinbot.profiled_qc import ProfiledQuantumChemistry
            return ProfiledQuantumChemistry(par)
        return super().__new__(cls)

    def __init__(self, par):
        ensure_current_run(create=True)
        self.par = par
        self.qc = par['qc'].lower()
        self.method = par['method']
        self.basis = par['basis']
        self.scan_method = par['scan_method']
        self.scan_basis = par['scan_basis']
        self.bls_method = par['barrierless_saddle_method']
        self.bls_basis = par['barrierless_saddle_basis']
        self.high_level_method = par['high_level_method']
        self.high_level_basis = par['high_level_basis']
        self.bls_high_level_method = par['barrierless_saddle_method_high']
        self.bls_high_level_basis = par['barrierless_saddle_basis_high']
        self.vts_method = par['vrc_tst_scan_method']
        self.vts_basis = par['vrc_tst_scan_basis']
        self.single_point_template = par['single_point_template']
        self.single_point_key = par['single_point_key']
        self.integral = par['integral']
        self.orcablocks = par['orcablocks']
        self.orcasimpleinput = par['orcasimpleinput']
        self.opt = par['opt']
        self.ppn = par['ppn']
        self.queuing = par['queuing']
        self.queue_name = par['queue_name']
        self.slurm_feature = par['slurm_feature']
        self.zf = par['zf']
        self.db = connect('kinbot.db')
        self.job_ids = {}
        # True means this controller submitted the job with a wrapper that
        # writes ``<scheduler id> <exit code>``. Profiled restarts may restore
        # older jobs whose one-field marker must retain legacy semantics.
        self.job_generation_tagged = {}
        self._slurm_excluded_nodes = _load_slurm_node_state(par)
        if self.queuing == 'slurm':
            _export_slurm_excluded_nodes(self._slurm_excluded_nodes)
        self.irc_maxpoints = par['irc_maxpoints']
        self.irc_stepsize = par['irc_stepsize']
        self.irc_maxcycle = par['irc_maxcycle']
        self.qc_command = par['qc_command']
        if par['slurm_feature'] == '':
            self.slurm_feature = ''
        else:
            self.slurm_feature = '#SBATCH -C ' + par['slurm_feature']
        self.queue_job_limit = par['queue_job_limit']
        queue_root = os.environ.get('KINBOT_QUEUE_ROOT', '').strip()
        self.queue_job_root = (os.path.realpath(os.path.expanduser(queue_root))
                               if queue_root else None)
        self.queue_job_root_recursive = (
            os.environ.get('KINBOT_QUEUE_ROOT_RECURSIVE', '') == '1')
        # An empty username is the portable input default.  Resolve it at
        # runtime so scheduler admission control cannot silently query an
        # empty user and submit past ``queue_job_limit``.
        self.username = str(par.get('username', '') or '').strip()
        if not self.username:
            self.username = (os.environ.get('USER', '').strip()
                             or os.environ.get('LOGNAME', '').strip()
                             or getpass.getuser().strip())
        if self.queuing in ('slurm', 'pbs') and not self.username:
            raise ValueError(
                'Cannot determine the scheduler username; set username or '
                'the USER environment variable.')
        self.use_sella = par['use_sella']
        self.frequency_mode = par.get('frequency_mode', 'auto')
        if not self.use_sella and self.qc.lower() == 'nn_pes':
            logger.warning('NNPES needs Sella optimizer. Turning "use_sella" on.')
            self.use_sella = True
        if not self.use_sella and self.qc.lower() == 'orca':
            logger.warning('ORCA not supported without Sella. Turning "use_sella" on.')
            self.use_sella = True
        if self.use_sella and self.qc.lower() == 'fc':
            self.use_sella = False  # Sella template files are separate from fairchem

    def _prepare_stereochemistry(self, species):
        """Keep the requested well identity before an optimization changes it."""
        if not getattr(species, 'wellorts', 0):
            species.optical_population = getattr(self, 'par', {}).get('optical_population', 'specified')
            optical_scope(species, species.optical_population)

    def _worker_site_setup(self):
        """Create the portable environment sourced by ordinary Slurm jobs."""
        if self.queuing != 'slurm' or self.qc not in _QC_SITE_BACKENDS:
            return ''
        backend, default_program = _QC_SITE_BACKENDS[self.qc]
        try:
            command = shlex.split(str(getattr(self, 'qc_command', '') or ''))
        except ValueError as error:
            raise ValueError('qc_command is not valid shell syntax.') from error
        program = command[0] if command else default_program
        cache = getattr(self, '_worker_site_setup_cache', {})
        key = (backend, program)
        path = cache.get(key)
        if path is None:
            from kinbot.anl.site import _executable, render_site_setup
            if _executable(program, backend) is None:
                raise FileNotFoundError(
                    f'Required {backend} executable {program!r} is unavailable; '
                    'configure its command, installation root, or environment '
                    'module before submitting this calculation.')
            path = Path(f'.kinbot_site_setup_{backend}.sh').resolve()
            content = render_site_setup({backend: {program}})
            with _SUBMISSION_LOCK:
                temporary = path.with_name(path.name + f'.tmp.{os.getpid()}')
                temporary.write_text(content)
                os.chmod(temporary, 0o700)
                os.replace(temporary, path)
            cache[key] = path
            self._worker_site_setup_cache = cache
        return (f'export KINBOT_BACKEND={shlex.quote(backend)}\n'
                f'source {shlex.quote(str(path))} || exit $?\n')

    def get_qc_arguments(self, job, mult, charge, nel, ts=0, step=0, max_step=0,
                         irc=None, scan=0, high_level=0, hir=0,
                         start_from_geom=0, rigid=0, aie=0, L3=False, vts=0):
        '''
        Method to get the argument to pass to ase, which are then passed to the qc codes.
        Job: name of the job
        mult: multiplicity of the job
        charge: charge of the job
        ts: 1 for transition state searches, 0 for wells and products
        step: Step number in the transition state search
        max_step: total number of necessary steps to get to the final transition state structure
            This is different for every reaction family
        irc: direction of the irc, None if this is not an irc job
        scan: this calculation is part of a scan of a bond length
        '''

        if self.qc == 'gauss' or (self.qc == 'nn_pes' and step < max_step):
            # arguments for Gaussian
            kwargs = {
                'method': self.method,
                'basis': self.basis,
                'nprocshared': min(nel, self.ppn),
                'mem': '700MW',
                # ASE splits labels containing a directory into
                # ``directory`` and ``prefix`` and executes Gaussian from
                # that directory.  A checkpoint containing the original
                # path would therefore become ``hir/hir/...`` (or the
                # equivalent for any other job subdirectory).
                'chk': os.path.basename(job),
                'label': job,
                'Symm': 'None',
                'mult': mult,
                'charge': charge,
                'scf': 'xqc',
                'pop': 'None',
            }
            if self.par['guessmix'] == 1 or \
                'barrierless_saddle' in job or \
                'bls' in job or \
                (mult == 1 and 'R_Addition_MultipleBond' in job):
                kwargs['guess'] = 'Mix,Always'
            if ts:
                # arguments for transition state searches
                if not self.par.get('_profiled_backend'):
                    kwargs['method'] = 'am1'
                    kwargs['basis'] = ''

                if step == 0:
                    if not self.par['bimol']:
                        kwargs['opt'] = 'ModRedun,Loose,CalcFC'
                    else:
                        kwargs['opt'] = 'ModRedun,Loose,CalcFC'
                        kwargs['method'] = self.method
                        kwargs['basis'] = self.basis
                elif step < max_step:
                    kwargs['opt'] = 'ModRedun,Loose,CalcFC'
                    kwargs['guess'] = 'Read'
                    if self.par['guessmix'] == 1 or \
                        'barrierless_saddle' in job or \
                        (mult == 1 and 'R_Addition_MultipleBond' in job):
                        kwargs['guess'] = 'Read,Mix,Always'
                    if self.par['bimol']:
                        kwargs['method'] = self.method
                        kwargs['basis'] = self.basis
                        kwargs['opt'] = 'QST3,AddRedundant'
                        kwargs.pop('guess', None)
                else:
                    kwargs['method'] = self.method
                    kwargs['basis'] = self.basis
                    if self.par['calcall_ts'] == 1:
                        kwargs['opt'] = 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999'
                        # not sending the frequency calculation for CalcAll
                    else:
                        kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,MaxCycle=999'
                        kwargs['freq'] = 'freq'
            elif vts:
                kwargs['method'] = self.vts_method
                kwargs['basis'] = self.vts_basis
                kwargs['opt'] = 'ModRedun,CalcFC,MaxCycles=999'
                kwargs['guess'] = 'Mix, Always'
                kwargs.pop('freq', None)
                kwargs.pop('chk', None)
            else:
                kwargs['freq'] = 'freq'  # ? 
            if (scan or 'R_Addition_MultipleBond' in job) and not vts:
                kwargs['method'] = self.scan_method 
                kwargs['basis'] = self.scan_basis
            if 'barrierless_saddle' in job or 'bls' in job:
                kwargs['method'] = self.bls_method
                kwargs['basis'] = self.bls_basis
            if irc is not None:
                # arguments for the irc calculations
                if start_from_geom == 0:
                    kwargs['geom'] = 'AllCheck,NoKeepConstants'
                    if self.par['guessmix'] == 1 or \
                        'barrierless_saddle' in job or \
                        (mult == 1 and 'R_Addition_MultipleBond' in job):
                        kwargs['guess'] = 'Read,Mix'  # Always is illegal here
                    else:
                        kwargs['guess'] = 'Read'
                    kwargs['irc'] = 'RCFC,{},MaxPoints={},StepSize={},MaxCycle={}'.format(irc, self.irc_maxpoints, self.irc_stepsize, self.irc_maxcycle)
                else:
                    kwargs['irc'] = 'RCFC,CalcFC,{},MaxPoints={},StepSize={},MaxCycle={}'.format(irc, self.irc_maxpoints, self.irc_stepsize, self.irc_maxcycle)
                kwargs.pop('freq', None)
            if high_level:
                kwargs['method'] = self.high_level_method
                kwargs['basis'] = self.high_level_basis
                if len(self.opt) > 0:
                    if self.par['calcall_ts'] == 1 and ts == 1:
                        kwargs['opt'] = 'NoFreeze,TS,CalcAll,NoEigentest,' \
                                        'MaxCycle=999,{}'.format(self.opt)  
                    else:
                        kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,' \
                                        'MaxCycle=999,{}'.format(self.opt)  # to overwrite possible CalcAll
                        kwargs['freq'] = 'freq'
                else:
                    if self.par['calcall_ts'] == 1 and ts == 1:
                        kwargs['opt'] = 'NoFreeze,TS,CalcAll,NoEigentest,' \
                                        'MaxCycle=999'  
                    else:
                        kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,' \
                                        'MaxCycle=999'  # to overwrite possible CalcAll
                        kwargs['freq'] = 'freq'
                if len(self.integral) > 0:
                    kwargs['integral'] = self.integral
                if 'barrierless_saddle' in job:  # completely overwrite normal settings
                    kwargs['method'] = self.bls_high_level_method
                    kwargs['basis'] = self.bls_high_level_basis
                    kwargs['opt'] = 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999'  # to overwrite possible CalcAll
                    kwargs.pop('freq', None)
            if hir:
                kwargs['opt'] = 'ModRedun,CalcFC'
                #if (not ts) or (ts and (not self.par['calcall_ts'])):
                kwargs.pop('freq', None)
                if ts:
                    if self.par['hir_maxcycle'] is None:
                        kwargs['opt'] = 'ModRedun,CalcFC,TS,NoEigentest'
                    else:
                        kwargs['opt'] = 'ModRedun,CalcFC,TS,NoEigentest,MaxCycle={}'.format(self.par['hir_maxcycle'])
                if rigid == 1:
                    kwargs.pop('freq', None)
                    kwargs.pop('opt', None)
            if 'hom_sci' in job and not vts:
                kwargs.pop('opt', None)
            if aie:
                kwargs = {
                    'method': 'cbs-qb3',
                    'basis': None,
                    'nprocshared': min(nel, self.ppn),
                    'mem': '700MW',
                    'label': job,
                    'Symm': 'None',
                    'mult': mult,
                    'charge': charge,
                    'scf': 'xqc'
                }

            return kwargs

        elif self.qc == 'nwchem':
            # arguments for NWChem
            odft = mult > 1
            kwargs = {
                'xc': self.method,
                'basis': self.basis,
                'scratch_dir': 'scratch',
                'permanent_dir': './perm',
                'label': job,
                'mult': mult,
                'charge': charge,
                'odft': odft,
                'task': 'optimize',
                'driver': 'tight',
                'maxiter': 100,
                'clear': 1,
                'inhess': 1
            }
            if ts:
                # arguments for transition state searches
                if step == max_step:
                    kwargs['task'] = 'saddle'
                else:
                    kwargs['driver'] = 'default'
            if irc is not None:
                # arguments for the irc calculations
                irc_kwargs = {
                    'task': 'mepgs',
                    'mepgs': 1,
                    'maxmep': 30,
                    'maxiter': 20,
                    'inhess': 2,
                    'xyz': 1,
                    'evib': 0.0005,
                    'stride': 0.1,
                    'direction': irc
                }
                kwargs.update(irc_kwargs)
            return kwargs

        elif self.qc == 'qchem':
            # arguments for QChem
            kwargs = {
                'label': job,
                'jobtype': 'opt',
                'vibman_print': '4',
                'nt': min(nel, self.ppn),
                'method': self.method,
                'basis': self.basis,
                'unrestricted': 'True',
                'scf_algorithm': 'diis_gdm',
                'max_scf_cycles': '100',
                'geom_opt_max_cycles': '500',
                'multiplicity': mult,
                'charge': charge

            }
            if ts:
                if 0 < step < max_step:
                    kwargs['scf_guess'] = 'read'
                if step >= max_step:
                    kwargs['jobtype'] = 'ts'
                elif not scan and 'R_Addition_MultipleBond' not in job:
                        kwargs['method'] = 'b3lyp'
                        kwargs['basis'] = 'sto-3g'
                        kwargs['scf_convergence'] = '4'
                        kwargs['geom_opt_tol_gradient'] = '1500'
                        kwargs['geom_opt_tol_displacement'] = '6000'
                        kwargs['geom_opt_tol_energy'] = '500'

            if scan or 'R_Addition_MultipleBond' in job:
                kwargs['method'] = self.scan_method
                kwargs['basis'] = self.scan_basis
            if 'barrierless_saddle' in job or 'bls' in job:
                kwargs['method'] = self.bls_method
                kwargs['basis'] = self.bls_basis
            if vts:
                kwargs['method'] = self.vts_method
                kwargs['basis'] = self.vts_basis
            if high_level:
                kwargs['method'] = self.high_level_method
                kwargs['basis'] = self.high_level_basis
                kwargs['vibman_print'] = '4'
            if irc is not None:
                kwargs['jobtype'] = 'freq'
                kwargs['xc_grid'] = '3'
                kwargs['vibman_print'] = '4'
                kwargs['rpath_max_stepsize'] = str(self.irc_stepsize * 10)
                kwargs['rpath_max_cycles'] = str(self.irc_maxpoints)
                if irc == 'forward':
                    kwargs['rpath_direction'] = '-1'
                else:
                    kwargs['rpath_direction'] = '1'
        elif self.qc == 'orca':
            kwargs = {
                'label': job,
                'orcasimpleinput': f'EnGrad TightSCF {self.method} {self.basis} {self.orcasimpleinput}',
                'orcablocks': f'{self.orcablocks} %pal nprocs {min(nel, self.ppn)} end',
                'charge': charge,
                'mult': mult,
                'task': 'gradient',
                'profile': self.qc_command,
                'directory': os.path.basename(job)
            }
            if high_level:
                kwargs['orcasimpleinput'] = f'EnGrad TightSCF {self.high_level_method} {self.high_level_basis} {self.orcasimpleinput}'
                kwargs['orcablocks'] = f'{self.orcablocks} %pal nprocs {min(nel, self.ppn)} end'

        elif self.qc == 'nn_pes':
            if self.par['nn_model']:
                kwargs = {'fname': self.par['nn_model']}
            else:
                kwargs = {}
        elif self.qc == 'fc':
            kwargs = {
                    'mult': mult,
                    'charge': charge
            }

        return kwargs

    def qc_hir(self, species, geom, rot_index, ang_index, fix, rigid,
               recovery_attempt=0, steps=None):
        '''Creates a constrained geometry optimization input and runs it.
        wellorts: 0 for wells and 1 for saddle points
        rot_index: index of the rotor in the molecule
        ang_index: index for the current size of the angle
        fix: four atoms of the dihedral that is currently fixed
        '''
        from kinbot.geometry import calc_dihedral
        if species.wellorts:
            job = 'hir/' + species.name + '_hir_' + str(rot_index) + '_' + str(ang_index).zfill(2)
        else:
            job = 'hir/' + routing_name(species) + '_hir_' + str(rot_index) + '_' + str(ang_index).zfill(2)

        if (not isinstance(recovery_attempt, int)
                or isinstance(recovery_attempt, bool)
                or recovery_attempt < 0):
            raise ValueError('recovery_attempt must be a nonnegative integer.')
        if recovery_attempt:
            job += f'_recovery_{recovery_attempt}'
        if steps is None:
            steps = self.par['sella_steps']
        if (not isinstance(steps, int) or isinstance(steps, bool) or steps < 1):
            raise ValueError('HIR optimization steps must be a positive integer.')

        self._prepare_stereochemistry(species)
        try:
            self._check_hir_definition(job, fix[0], bool(self.db.count(name=job)),
                                       geometry=geom, atoms=species.atom)
        except ValueError as error:
            # The previous scan is evidence, not a result for this new input.
            # A stable suffix permits resuming the replacement without either
            # overwriting the old calculation or repeatedly submitting it.
            import hashlib
            import json
            signature = json.dumps([plain_symbols(species.atom), plain_geometry(geom),
                                    list(map(int, fix[0])), bool(rigid)], sort_keys=True)
            job += '_recovery_' + hashlib.sha256(signature.encode()).hexdigest()[:16]
            logger.warning('%s; using fresh HIR calculation %s.', error, job)
            self._check_hir_definition(job, fix[0], bool(self.db.count(name=job)),
                                       geometry=geom, atoms=species.atom)
        if self.check_qc(job) in ('normal', 'error', 'running'):
            return job

        kwargs = self.get_qc_arguments(job, species.mult, species.charge, species.nel,
                                       ts=species.wellorts, step=1, max_step=1,
                                       high_level=1, hir=1, rigid=rigid)
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        if self.qc == 'gauss': 
            code = 'gaussian'
            Code = 'Gaussian'
            if not self.use_sella:
                kwargs['addsec'] = f'{" ".join(str(f) for f in fix[0])} F'
                kwargs.pop('chk', None)
        elif self.qc == 'qchem':
            code = 'qchem'
            Code = 'QChem'
            if not self.use_sella:
                dihedral = calc_dihedral(geom[fix[0][0]-1], geom[fix[0][1]-1], 
                                         geom[fix[0][2]-1], geom[fix[0][3]-1])[0]
                kwargs['addsec'] = '\n$opt\nCONSTRAINT\ntors ' \
                                   f'{" ".join(str(f) for f in fix[0])} ' \
                                   f'{dihedral}\nENDCONSTRAINT\n$end'
        elif self.qc == 'nwchem':
            code = 'nwchem'
            Code = 'NWChem'
        elif self.qc == 'orca':
            code = 'orca'
            Code = 'ORCA'
        elif self.qc == 'nn_pes':
            code = 'nn_pes'
            Code = 'Nn_surr'
        elif self.qc == 'fc':
            code = 'fairchem'
            Code = 'Fairchem'
        else:
            raise ValueError(f'Unexpected value for qc parameter: {self.qc}')
        if self.use_sella:
            # Gaussian/Sella reuses the first force evaluation through
            # Guess=Read. Keep the per-scan checkpoint so later force calls
            # have a wavefunction to read.
            kwargs.pop('opt', None)
            kwargs.pop('freq', None)
            template_file = f'{kb_path}/tpl/ase_sella_hir.tpl.py'
        else:
            template_file = f'{kb_path}/tpl/ase_{self.qc}_hir.tpl.py'
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=plain_symbols(species.atom),
                                   geom=plain_geometry(geom),
                                   ppn=min(species.nel, self.ppn),
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   fix=fix[0],
                                   code=code,  # Sella
                                   Code=Code,  # Sella
                                   order=species.wellorts,  # Sella
                                   sella_kwargs=self.par['sella_kwargs'],  # Sella
                                   fmax=self.par['sella_fmax'],
                                   steps=steps,
                                   fc_model_path=self.par['fc_model_path'],
                                   fc_task_name=self.par['fc_task_name'],
                                   fc_device=self.par['fc_device'],
                                   )

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job, min(species.nel, self.ppn))

        return job

    @staticmethod
    def _check_hir_definition(job, dihedral, has_result=False, geometry=None, atoms=None):
        """Do not relabel a saved scan after changing its torsion definition.

        Read literal input only, never execute a worker. This deliberately
        checks the defining atoms and, when requested, the literal initial
        geometry. Rigid rotations/translations of that input are equivalent.
        """
        import ast
        import re
        from pathlib import Path
        path = Path(job + '.py')
        if not path.exists():
            if not has_result:
                return
            raise ValueError(f'{job}: a cached HIR result has no saved scan input. '
                             'Run fresh scans in a new directory; existing results were preserved.')
        saved = None
        saved_geometry, saved_atoms = None, None
        try:
            for node in ast.parse(path.read_text()).body:
                if not isinstance(node, ast.Assign):
                    continue
                names = [n.id for n in node.targets if isinstance(n, ast.Name)]
                if 'base_0_fix' in names and isinstance(node.value, ast.ListComp):
                    saved = ast.literal_eval(node.value.generators[0].iter)
                elif isinstance(node.value, ast.Call) and getattr(node.value.func, 'id', None) == 'Atoms':
                    values = {key.arg: ast.literal_eval(key.value) for key in node.value.keywords
                              if key.arg in ('symbols', 'positions')}
                    saved_geometry, saved_atoms = values.get('positions'), values.get('symbols')
                elif 'kwargs' in names:
                    section = ast.literal_eval(node.value).get('addsec', '')
                    match = re.search(r'(?:^|tors\s+)(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', section)
                    if match:
                        saved = list(map(int, match.groups()))
        except (ValueError, SyntaxError, TypeError, AttributeError):
            pass
        if saved is None or list(saved) not in (list(dihedral), list(dihedral)[::-1]):
            raise ValueError(f'{job}: saved HIR dihedral {saved} differs from requested {dihedral}. '
                             'Run fresh scans in a new directory; existing files were preserved.')
        if geometry is not None:
            from kinbot.molecular_symmetry import proper_rmsd
            source = np.asarray(saved_geometry, dtype=float)
            target = np.asarray(geometry, dtype=float)
            if (saved_atoms != list(map(str, atoms)) or source.shape != target.shape
                    or not np.all(np.isfinite(source))
                    or proper_rmsd(source, target) > 1.e-6):
                raise ValueError(f'{job}: saved HIR input geometry or atom order differs '
                                 'from the requested scan; existing files were preserved.')

    def qc_ring_conf(self, species, geom, fix, change, conf_idx, dih_idx):
        '''
        Creates a constrained geometry optimization input for the
        conformational search of cyclic structures and runs it.
        Makes use of Sella.

        qc: 'gauss' or 'nwchem' or 'qchem'
        wellorts: 0 for wells and 1 for saddle points
        conf_idx: index of ring conformer in the conformer search
        dih_idx: index of the dih for this conformer
        '''
        if species.wellorts:
            job = 'conf/' + species.name + '_r' + str(conf_idx).zfill(self.zf) \
                  + '_' + str(dih_idx).zfill(self.zf)
        else:
            job = 'conf/' + routing_name(species) + '_r' \
                  + str(conf_idx).zfill(self.zf) + '_' \
                  + str(dih_idx).zfill(self.zf)

        self._prepare_stereochemistry(species)

        kwargs = self.get_qc_arguments(job, species.mult, species.charge, species.nel,
                                       ts=species.wellorts, step=1, max_step=1, 
                                       hir=1)
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        kwargs.pop('opt', None)
        kwargs.pop('chk', None)
        if not self.par.get('_profiled_backend'):
            kwargs['method'] = 'am1'
            kwargs['basis'] = ''

        if self.qc == 'gauss':
            code = 'gaussian'
            Code = 'Gaussian'
        elif self.qc == 'qchem':
            code = 'qchem'
            Code = 'QChem'
        elif self.qc == 'nwchem':
            code = 'nwchem'
            Code = 'NWChem'
        elif self.qc == 'orca':
            code = 'orca'
            Code = 'ORCA'
        elif self.qc == 'nn_pes':
            kwargs.pop('method', None)
            kwargs.pop('basis', None)
            code = 'nn_pes'
            Code = 'Nn_surr'
        elif self.qc == 'fc':
            code = 'fairchem'
            Code = 'Fairchem'
        else:
            raise ValueError(f'Unexpected value for qc parameter: {self.qc}')
        if self.use_sella:
            template_file = f'{kb_path}/tpl/ase_sella_ring_conf.tpl.py'
        else:
            template_file = f'{kb_path}/tpl/ase_{self.qc}_ring_conf.tpl.py'
        skw = self.par['sella_kwargs'].copy()
        skw['internal'] = False
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=plain_symbols(species.atom),
                                   geom=plain_geometry(geom),
                                   fix=fix,
                                   change=change,
                                   ppn=min(species.nel, self.ppn),
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   code=code,  # Sella
                                   Code=Code,  # Sella
                                   sella_kwargs=skw, # internal turned off, constraints used instead
                                   fmax=self.par['sella_fmax'],
                                   steps=self.par['sella_steps'],
                                   fc_model_path=self.par['fc_model_path'],
                                   fc_task_name=self.par['fc_task_name'],
                                   fc_device=self.par['fc_device'],
                                   )

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job, min(species.nel, self.ppn))
        return 0

    def qc_conf(self, species, geom, index, semi_emp=0):
        '''
        Creates a geometry optimization input for the conformational search and runs it.
        qc: 'gauss' or 'nwchem' or 'qchem'
        wellorts: 0 for wells and 1 for saddle points
        index: >=0 for sampling, each job will get numbered with index
        '''
        add = ''
        if semi_emp:
            add = 'semi_emp_'
        if species.wellorts:
            job = 'conf/' + species.name + '_' + add + str(index).zfill(self.zf)
        else:
            job = 'conf/' + routing_name(species) + '_' + add \
                  + str(index).zfill(self.zf)

        self._prepare_stereochemistry(species)

        if species.wellorts:
            kwargs = self.get_qc_arguments(job, species.mult, species.charge, species.nel,
                                           ts=1, step=1, max_step=1)
        else:
            kwargs = self.get_qc_arguments(job, species.mult, species.charge, species.nel)
            if self.qc == 'gauss' and not self.use_sella:
                if self.par['opt'].casefold() == 'Tight'.casefold(): 
                    kwargs['opt'] = 'CalcFC, Tight'
                else:
                    kwargs['opt'] = 'CalcFC'

        if self.qc == 'gauss':
            code = 'gaussian'
            Code = 'Gaussian'
            # Native Gaussian conformers retain selected Hessian checkpoints
            # for the MC optical check. Sella also needs a checkpoint for its
            # Guess=Read force calls, independently of the MC setting.
            keep_checkpoint = (self.par.get('multi_conf_tst')
                               and not self.use_sella and not semi_emp)
            if self.use_sella:
                # The external optimizer uses Guess=Read after its first
                # force call. Gaussian runs inside the label's directory.
                kwargs['chk'] = os.path.basename(job)
            elif not keep_checkpoint:
                kwargs.pop('chk', None)
            else:
                # Native Gaussian templates execute from the run root, where
                # the calculation label is also the checkpoint path.
                kwargs['chk'] = job
        elif self.qc == 'qchem':
            code = 'qchem'
            Code = 'QChem'
        elif self.qc == 'nwchem':
            code = 'nwchem'
            Code = 'NWChem'
        elif self.qc == 'orca':
            code = 'orca'
            Code = 'ORCA'
        elif self.qc == 'nn_pes':
            code = 'nn_pes'
            Code = 'Nn_surr'
        elif self.qc == 'fc':
            code = 'fairchem'
            Code = 'Fairchem'
        else:
            raise ValueError(f'Unexpected value for qc parameter: {self.qc}')
        
        if semi_emp:
            kwargs['method'] = self.par['semi_emp_method']
            kwargs['basis'] = ''
        if species.natom < 3:
            kwargs.pop('Symm', None)
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        if self.use_sella:
            kwargs.pop('opt', None)
            kwargs.pop('freq', None)
            template_file = f'{kb_path}/tpl/ase_sella_opt_well.tpl.py'
        else:        
            template_file = f'{kb_path}/tpl/ase_{self.qc}_opt_well.tpl.py'
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   checkpoint=os.path.basename(job),
                                   kwargs=kwargs,
                                   atom=plain_symbols(species.atom),
                                   geom=plain_geometry(geom),
                                   ppn=min(species.nel, self.ppn),
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   code=code,    # Sella
                                   Code=Code,    # Sella
                                   order=species.wellorts,      # Sella
                                   sella_kwargs=(worker_sella_kwargs(self.par) if self.use_sella
                                                 else self.par['sella_kwargs']),  # Sella
                                   fmax=self.par['sella_fmax'],
                                   steps=self.par['sella_steps'],
                                   frequency_mode=self.frequency_mode,
                                   fc_model_path=self.par['fc_model_path'],
                                   fc_task_name=self.par['fc_task_name'],
                                   fc_device=self.par['fc_device'],
                                   )
        
        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job, min(species.nel, self.ppn))

        return 0

    def qc_freq(self, species, source_job, high_level=0):
        """Recover native frequencies/Hessian at the selected, fixed geometry.

        Native jobs can lack a checkpoint or printed Hessian. Keep their
        records intact and use a separate, same-level frequency calculation.
        """
        if ((self.use_sella and self.frequency_mode != 'native_hessian')
                or self.qc not in ('gauss', 'qchem')):
            raise ValueError(f'Missing stored Hessian for {source_job} ({self.qc}).')
        rows = list(self.db.select(name=source_job))
        source_id = rows[-1].id
        # Include the immutable source observation to avoid reusing a recovery
        # after the same optimization name has been calculated again.
        job = f'{source_job}_freq_recovery_{source_id}'
        kwargs = self.get_qc_arguments(job, species.mult, species.charge,
            species.nel, ts=species.wellorts, step=1, max_step=1,
            high_level=high_level)
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        kwargs.pop('opt', None)
        if self.qc == 'gauss':
            # Recovery has no predecessor checkpoint under its new job name.
            # Preserve electronic-state choices (e.g. broken-symmetry Mix).
            guess = kwargs.get('guess')
            if isinstance(guess, str):
                options = [option.strip() for option in guess.strip('()').split(',')
                           if option.strip().lower() != 'read']
                if any(options):
                    kwargs['guess'] = ','.join(options)
                else:
                    kwargs.pop('guess', None)
            kwargs.update(freq='', chk=os.path.basename(job), label=job)
            template_file = f'{kb_path}/tpl/ase_gauss_opt_well.tpl.py'
        else:
            kwargs.update(jobtype='freq', vibman_print='4', xc_grid='3',
                          label=job + '_freq')
            template_file = f'{kb_path}/tpl/ase_qchem_freq.tpl.py'
        with open(template_file) as handle:
            template = handle.read().format(label=job, kwargs=kwargs,
                atom=plain_symbols(species.atom), geom=plain_geometry(species.geom),
                ppn=min(species.nel, self.ppn), qc_command=self.qc_command,
                working_dir=os.getcwd())
        with open(job + '.py', 'w') as handle:
            handle.write(template)
        self.submit_qc(job, min(species.nel, self.ppn))
        return job

    def qc_aie(self, species, geom, ext):
        '''
        Sets up AIE calculations. Currently results are not digested by KinBot.
        Only for well, not for saddles!
        qc: 'gauss' or 'nwchem' or 'qchem'
        index: the index of the conformer as in conf search
        '''
        self._prepare_stereochemistry(species)
        job0 = f'aie/{routing_name(species)}_AIE0_{ext}'  # neutral
        job1 = f'aie/{routing_name(species)}_AIE1_{ext}'  # cation
        kwargs0 = self.get_qc_arguments(job0, species.mult, species.charge, species.nel, aie=1)
        if species.mult == 1: m1 = 2
        elif species.mult == 2: m1 = 1
        elif species.mult == 3: m1 = 2
        kwargs1 = self.get_qc_arguments(job1, m1, species.charge + 1, species.nel - 1, aie=1)
        template_file = f'{kb_path}/tpl/ase_{self.qc}_opt_well.tpl.py'
        template = open(template_file, 'r').read()
        t0 = template.format(label=job0,
                             kwargs=kwargs0,
                             atom=plain_symbols(species.atom),
                             geom=plain_geometry(geom),
                             ppn=min(species.nel, self.ppn),
                             qc_command=self.qc_command,
                             working_dir=os.getcwd())
        t1 = template.format(label=job1,
                             kwargs=kwargs1,
                             atom=plain_symbols(species.atom),
                             geom=plain_geometry(geom),
                             ppn=min(species.nel-1, self.ppn),
                             qc_command=self.qc_command,
                             working_dir=os.getcwd())
        
        for job, template, electrons in ((job0, t0, species.nel),
                                         (job1, t1, species.nel - 1)):
            with open(f'{job}.py', 'w') as f:
                f.write(template)
            self.submit_qc(job, min(electrons, self.ppn))

        return 0

    def qc_opt(self, species, geom, high_level=0, mp2=0, bls=0, ext=None, 
               fdir=None, do_vdW=False):
        '''
        Creates a geometry optimization input and runs it.
        '''
        #if mp2 and (self.qc == 'fc' or self.qc == 'orca'):
        #    return 0
        if do_vdW:
            if high_level:
                job = f'{species.name}_high'
        elif ext is None:
            job = routing_name(species) + '_well'
            if high_level:
                job = routing_name(species) + '_well_high'
            if mp2:
                job = routing_name(species) + '_well_mp2'
            if bls:
                job = routing_name(species) + '_well_bls'

        else:
            job = routing_name(species) + ext

        if fdir is not None:
            job = f'{fdir}/{job}'

        if not do_vdW:
            self._prepare_stereochemistry(species)

        if not do_vdW and not getattr(species, 'wellorts', 0):
            from kinbot.stereo_routing import guard_well_job
            guard_well_job(self, species, geom, job)

        # TODO: Code exceptions into their own function/py script that opt can call.
        # TODO: Fix symmetry numbers for calcs as well if needed
        mult = exceptions.get_multiplicity(species.chemid, species.mult)

        kwargs = self.get_qc_arguments(job, mult, species.charge, species.nel,
                                       high_level=high_level)
        if self.qc != 'qchem' and 'opt' not in kwargs:
            kwargs['opt'] = ''

        if self.qc == 'gauss':
            code = 'gaussian'
            Code = 'Gaussian'
            if self.par['opt'].casefold() == 'Tight'.casefold(): 
                kwargs['opt'] = 'CalcFC, Tight'
            else:
                kwargs['opt'] = 'CalcFC'
        elif self.qc == 'qchem':
            code = 'qchem'
            Code = 'QChem'
        elif self.qc == 'nwchem':
            code = 'nwchem'
            Code = 'NWChem'   
        elif self.qc == 'orca':
            code = 'orca'
            Code = 'ORCA'
        elif self.qc == 'nn_pes':
            code = 'nn_pes'
            Code = 'Nn_surr'
        elif self.qc == 'fc':
            code = 'fairchem'
            Code = 'Fairchem'
        else:
            raise ValueError(f'Unexpected value for qc parameter: {self.qc}')
        
        if mp2:
            kwargs['method'] = self.scan_method
            kwargs['basis'] = self.scan_basis
        if high_level and self.qc == 'gauss' and self.opt:
            kwargs['opt'] = 'CalcFC, {}'.format(self.opt)
        if species.natom < 3:
            kwargs.pop('Symm', None)
        # the integral is set in the get_qc_arguments parts, bad design
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        if self.use_sella:
            kwargs.pop('opt', None)
            kwargs.pop('freq', None)
            template_file = f'{kb_path}/tpl/ase_sella_opt_well.tpl.py'
        else:
            template_file = f'{kb_path}/tpl/ase_{self.qc}_opt_well.tpl.py'
        
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   checkpoint=os.path.basename(job),
                                   kwargs=kwargs,
                                   atom=plain_symbols(species.atom),
                                   geom=plain_geometry(geom),
                                   ppn=min(species.nel, self.ppn),
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   code=code,    # Sella
                                   Code=Code,    # Sella
                                   order=0,      # Sella
                                   sella_kwargs=(worker_sella_kwargs(self.par) if self.use_sella
                                                 else self.par['sella_kwargs']), # Sella
                                   fmax=self.par['sella_fmax'],
                                   steps=self.par['sella_steps'],
                                   frequency_mode=self.frequency_mode,
                                   charge=species.charge,
                                   spin=species.mult,
                                   fc_model_path=self.par['fc_model_path'],
                                   fc_task_name=self.par['fc_task_name'],
                                   fc_device=self.par['fc_device'],
                                   )

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job, min(species.nel, self.ppn))
        return 0

    def qc_opt_ts(self, species, geom, high_level=0, ext=None, fdir=None):
        '''Creates a ts optimization input and runs it
        '''

        job = str(species.name)
        if ext is None:
            if high_level:
                job += '_high'
        else:
            job += ext

        if fdir is not None:
            job = f'{fdir}/{job}'

        kwargs = self.get_qc_arguments(
            job, species.mult, species.charge, species.nel, ts=1,
            step=1, max_step=1, high_level=high_level)
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        if self.qc == 'gauss':
            code = 'gaussian'
            Code = 'Gaussian'
        elif self.qc == 'qchem':
            code = 'qchem'
            Code = 'QChem'
        elif self.qc == 'nwchem':
            code = 'nwchem'
            Code = 'NWChem'
        elif self.qc == 'orca':
            code = 'orca'
            Code = 'ORCA'
        elif self.qc == 'nn_pes':
            code = 'nn_pes'
            Code = 'Nn_surr'
        elif self.qc == 'fc':
            code = 'fairchem'
            Code = 'Fairchem'
        else:
            raise ValueError(f'Unrecognized qc option: {self.qc}')
        
        if self.use_sella:
            kwargs.pop('addsec', None)
            kwargs.pop('opt', None)
            kwargs.pop('freq', None)
            template_file = f'{kb_path}/tpl/ase_sella_ts_end.tpl.py'
        else:
            template_file = f'{kb_path}/tpl/ase_{self.qc}_ts_end.tpl.py'
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=plain_symbols(species.atom),
                                   geom=plain_geometry(geom),
                                   ppn=min(species.nel, self.ppn),
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   code=code, # Sella
                                   Code=Code, # Sella
                                   sella_kwargs=(worker_sella_kwargs(self.par) if self.use_sella
                                                 else self.par['sella_kwargs']), # Sella
                                   fmax=self.par['sella_fmax'],
                                   steps=self.par['sella_steps'],
                                   frequency_mode=self.frequency_mode,
                                   fc_model_path=self.par['fc_model_path'],
                                   fc_task_name=self.par['fc_task_name'],
                                   fc_device=self.par['fc_device'],
                                   )

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job, min(species.nel, self.ppn))

        return 0

    def qc_vts_frag(self, frag):
        '''
        Runs VTS fragment optimization
        '''

        job = f'vrctst/{routing_name(frag)}_vts'
        self._prepare_stereochemistry(frag)
        mult = exceptions.get_multiplicity(frag.chemid, frag.mult)
        kwargs = self.get_qc_arguments(job, mult, frag.charge, frag.nel, vts=1)
        # Fragment calculations provide orbitals for HOMO-oriented pivot
        # points.  Keep that auxiliary role separate from the relaxed VRC
        # scan level, which may be a Molpro multireference method.
        kwargs['method'] = self.par.get(
            'vrc_tst_pivot_method',
            getattr(self, 'vts_method', kwargs.get('method')))
        kwargs['basis'] = self.par.get(
            'vrc_tst_pivot_basis',
            getattr(self, 'vts_basis', kwargs.get('basis')))
        # Add chk to fragments only
        kwargs['chk'] = f'{routing_name(frag)}_vts'

        if self.qc != 'gauss':
            raise ValueError(f'Only implemeted for Gaussian. Instead I got: {self.qc}')

        if frag.natom < 3:
            kwargs.pop('Symm', None)
        vrc_kwargs = self.par.get('vrc_tst_calc_kwargs')
        if vrc_kwargs is None:
            vrc_kwargs = self.par.get('calc_kwargs', {})
        if vrc_kwargs:
            kwargs = self.merge_kwargs(kwargs, vrc_kwargs)

        template_file = f'{kb_path}/tpl/ase_{self.qc}_opt_well.tpl.py'
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=plain_symbols(frag.atom),
                                   geom=plain_geometry(frag.geom),
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   )  # TODO add sella keywords

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job, min(frag.nel, self.ppn))
        
        return job 

    def qc_vts(self, reac, geom, step, equiv, asymptote, step0_geom):
        '''
        Creates a geometry optimization along a scan and runs it.
        reac: full reaction object
        '''

        if not asymptote:
            job = f'vrctst/{reac.instance_name}_vts_pt{str(step).zfill(2)}'
        else:
            job = f'vrctst/{reac.instance_name}_vts_pt_asymptote'
        self._prepare_stereochemistry(reac.species)
        mult = exceptions.get_multiplicity(reac.species.chemid, reac.species.mult)
        kwargs = self.get_qc_arguments(job, mult, reac.species.charge, reac.species.nel, vts=1)

        if self.qc == 'gauss' and not self.par['vrc_tst_scan_sella']:
            kwargs['addsec'] = f'{reac.scan_coo[0]+1} {reac.scan_coo[1]+1} F\n'
        
        if reac.species.natom < 3:
            kwargs.pop('Symm', None)
        vrc_kwargs = self.par.get('vrc_tst_calc_kwargs')
        if vrc_kwargs is None:
            vrc_kwargs = self.par.get('calc_kwargs', {})
        if vrc_kwargs:
            kwargs = self.merge_kwargs(kwargs, vrc_kwargs)
        
        if self.par['vrc_tst_scan_sella']:  # TODO sella globally
            kwargs.pop('opt', None)
            if self.qc == 'gauss':
                code = 'gaussian'
                Code = 'Gaussian'
            else:
                raise ValueError(f'Currently only Gaussian is supported: {self.qc}')
            template_file = f'{kb_path}/tpl/ase_sella_vts.tpl.py'
            template = open(template_file, 'r').read()
            template = template.format(label=job,
                                       scan_coo=reac.scan_coo,
                                       bonds=reac.irc_prod.bondlist,
                                       kwargs=kwargs,
                                       atom=plain_symbols(reac.species.atom),
                                       init_geom=plain_geometry(geom),
                                       qc_command=self.qc_command,
                                       working_dir=os.getcwd(),
                                       vts_ang_dev=self.par['vts_ang_dev'],
                                       scan_deviation=self.par['vrc_tst_scan_deviation'],
                                       frag_maps=list([list(rm) for rm in reac.maps]),
                                       froz_A_geom=plain_geometry(reac.products[0].geom),
                                       froz_A_atom=reac.products[0].atom,
                                       froz_B_geom=plain_geometry(reac.products[1].geom),
                                       froz_B_atom=reac.products[1].atom,
                                       code=code,
                                       Code=Code,
                                       sella_kwargs=self.par['sella_kwargs'],  # Sella
                                       equiv=equiv,
                                       asymptote=asymptote,
                                       frag_bonds_0=list([list(b) for b in reac.products[0].bond01]),
                                       frag_bonds_1=list([list(b) for b in reac.products[1].bond01]),
                                       step0_geom=step0_geom
                                       )

        else:
            template_file = f'{kb_path}/tpl/ase_{self.qc}_vts.tpl.py'
            template = open(template_file, 'r').read()
            template = template.format(label=job,
                                       scan_coo=reac.scan_coo,
                                       bonds=reac.irc_prod.bondlist,
                                       kwargs=kwargs,
                                       atom=plain_symbols(reac.species.atom),
                                       init_geom=plain_geometry(geom),
                                       qc_command=self.qc_command,
                                       working_dir=os.getcwd(),
                                       scan_deviation=self.par['vrc_tst_scan_deviation'],
                                       frag_maps=list([list(rm) for rm in reac.maps]),
                                       froz_A_geom=plain_geometry(reac.products[0].geom),
                                       froz_A_atom=reac.products[0].atom,
                                       froz_B_geom=plain_geometry(reac.products[1].geom),
                                       froz_B_atom=reac.products[1].atom,
                                       )

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job, min(reac.species.nel, self.ppn))
        return job 

    def invalidate_qc(self, job):
        """Retire a completed result before a HIR restart, for every backend.

        Preserve database observations and archive output/checkpoint files.
        A database marker also invalidates backends without a completion log.
        Queued/running jobs must never be invalidated.
        """
        if self.check_qc(job) == 'running':
            raise ValueError(f'Cannot invalidate a running calculation: {job}')
        rows = list(self.db.select(name=job))
        if not rows:
            return
        previous = rows[-1]
        marker = self.db.write(previous.toatoms(), name=job, data={
            'status': 0, 'reason': 'HIR restart', 'supersedes_row_id': previous.id})
        self._archive_outputs(job, marker)
        self.job_ids.pop(job, None)
        getattr(self, 'job_generation_tagged', {}).pop(job, None)

    @staticmethod
    def _archive_outputs(job, marker):
        for suffix in ('.log', '.out', '_sella.log', '.chk', '.fchk', '_freq.out'):
            path = job + suffix
            if os.path.isfile(path):
                os.replace(path, f'{path}.restart_{marker}')

    def publish_result(self, source, target):
        """Copy an accepted result and its native outputs to a conventional name.

        The normal database row is written last. An interrupted copy leaves an
        incomplete result, rather than pairing new properties with old files.
        Original calculations and earlier conventional outputs are preserved.
        """
        if source.name == target:
            return source.id
        rows = list(self.db.select(name=target))
        if (rows and rows[-1].data.get('status') == 'normal'
                and rows[-1].data.get('copied_from_job') == source.name
                and rows[-1].data.get('copied_from_row_id') == source.id):
            return rows[-1].id
        if self.check_qc(target) == 'running':
            raise ValueError(f'Cannot publish over a running calculation: {target}')
        os.makedirs(os.path.dirname(target) or '.', exist_ok=True)
        marker = self.db.write(source.toatoms(), name=target, data={
            'status': 0, 'reason': 'Publishing accepted optimization result'})
        self._archive_outputs(target, marker)
        for suffix in ('.log', '.out', '_sella.log', '.chk', '.fchk', '_freq.out'):
            if os.path.isfile(source.name + suffix):
                copyfile(source.name + suffix, target + suffix)
        data = dict(source.data)
        data.update(copied_from_job=source.name, copied_from_row_id=source.id)
        keys = dict(source.key_value_pairs, name=target)
        result = self.db.write(source.toatoms(), key_value_pairs=keys, data=data)
        self.job_ids.pop(target, None)
        getattr(self, 'job_generation_tagged', {}).pop(target, None)
        return result

    def submit_qc(self, job, nproc, singlejob=1, jobtype=None):
        """Submit under an optional cross-process queue admission lock.

        Automatic CBH references are ordinary KinBot controllers in separate
        directories.  They share one lock so their queue-count observation
        and ``sbatch`` call are atomic with respect to one another.  Existing
        standalone inputs do not set the environment variable and preserve
        their historical behavior.
        """
        with _SUBMISSION_LOCK:
            lock_path = os.environ.get('KINBOT_QUEUE_LOCK', '').strip()
            if not lock_path:
                return QuantumChemistry._submit_qc(
                    self, job, nproc, singlejob, jobtype)
            lock_path = os.path.abspath(os.path.expanduser(lock_path))
            os.makedirs(os.path.dirname(lock_path), exist_ok=True)
            with open(lock_path, 'a') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                try:
                    return QuantumChemistry._submit_qc(
                        self, job, nproc, singlejob, jobtype)
                finally:
                    fcntl.flock(lock, fcntl.LOCK_UN)

    def _submit_qc(self, job, nproc, singlejob=1, jobtype=None):
        '''Submit a job to the queue, unless the job:
            * has finished with Normal termination
            * has finished with Error termination
            * is currently running
        However, if the optional parameter singlejob is set to zero, then
        the job is run only if it has finished earlier with normal termination.
        This is for continuations, when the continuing jobs overwrite each other.
        If the number of jobs in the queue is larger than the user-set limit,
        KinBot will park here until resources are freed up.
        '''
        # if the logfile already exists, copy it with another name

        logger.debug(f'checking {job}')
        check = self.check_qc(job)
        if singlejob == 1:
            if check != 0:
                return 0
        else:
            if check == 'running':
                return 0

        from kinbot import run_policy
        shared_pool = self.queuing == 'slurm' and run_policy.find() is not None
        if self.queue_job_limit > 0 and not shared_pool:
            self.limit_jobs()

        # Keep the completed generation's handoff while waiting for admission
        # or a scheduler response. Deleting it here loses completed work if a
        # controller is replaced before the next allocation is submitted.
        # The batch body clears these files only when the new worker starts;
        # generation-tagged status checks fence the old marker in the meantime.

        try:
            if jobtype == 'am1' and self.par['q_temp_am1']:
                template_head_file = self.par['q_temp_am1']
            elif job.endswith('mp2') and self.par['q_temp_mp2']:
                template_head_file = self.par['q_temp_mp2']
            elif job.endswith('high') and self.par['q_temp_hi']:
                template_head_file = self.par['q_temp_hi']
            elif self.par['queue_template'] == '':
                template_head_file = f'{kb_path}/tpl/{self.queuing}.tpl'
            else:
                template_head_file = self.par['queue_template']
        except OSError:
            logger.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
            logger.error('Or no file is found at {}.'.format(self.par['queue_template']))
            logger.error('Exiting')
            sys.exit()

        if self.queuing == 'local':
            err_msg = f'Job {job} is missing in the database or the output ' \
                      'file is not present. Unable to run calculations when ' \
                      'queuing is \'local\'. To ignore calculations that have ' \
                      'not finished set "error_missing_local" to False.'
            if self.par['error_missing_local']:
                logger.error(err_msg)
                raise FileNotFoundError(err_msg)
            else:
                logger.warning(err_msg)
                return -1

        template_file = f'{kb_path}/tpl/{self.queuing}_python.tpl'
        if self.qc == 'orca':
            python_file = f'{os.path.basename(job)}.py'
        else:
            python_file = f'{job}.py'
        site_setup = (self._worker_site_setup()
                      if hasattr(self, '_worker_site_setup') else '')
        job_template = (open(template_head_file, 'r').read() + site_setup
                        + open(template_file, 'r').read())
        excluded_nodes = getattr(self, '_slurm_excluded_nodes', None)
        if excluded_nodes is None:
            excluded_nodes = _load_slurm_node_state(getattr(self, 'par', {}))
        excluded_nodes = sorted(excluded_nodes)
        slurm_exclude_directive = _slurm_exclude_directive(
            getattr(self, 'par', {}), excluded_nodes)
        if (self.queuing == 'slurm' and slurm_exclude_directive
                and '{slurm_exclude_directive}' not in job_template):
            # Preserve custom templates while keeping scheduler directives in
            # the header, before the first executable shell statement.
            lines = job_template.splitlines(keepends=True)
            indices = [index for index, line in enumerate(lines)
                       if line.startswith('#SBATCH')]
            insert_at = max(indices) + 1 if indices else 1
            lines.insert(insert_at, slurm_exclude_directive + '\n')
            job_template = ''.join(lines)
        if self.qc == 'orca': 
            if 'hir' not in job: 
                job_template += '\ncp -r $SCRATCH_DIR/* $SLURM_SUBMIT_DIR/.\ncd /scratch/$USER\nrm -rf $SCRATCH_DIR'
            else:
                job_template += '\ncp -r $SCRATCH_DIR/* $SLURM_SUBMIT_DIR/hir/.\ncd /scratch/$USER\nrm -rf $SCRATCH_DIR'

        if self.queuing == 'pbs':
            job_template = job_template.format(name=job, ppn=max(1, nproc), queue_name=self.queue_name,
                                               errdir='perm', python_file=python_file, arguments='',
                                               python_executable=shlex.quote(sys.executable))
        elif self.queuing == 'slurm':
            # Existing custom templates keep complete control of scheduler
            # directives. Bundled templates and opted-in custom templates
            # expose this placeholder.
            walltime_directive = (_ordinary_slurm_walltime(
                self.par, self.queue_name)
                if '{walltime_directive}' in job_template else '')
            job_template = job_template.format(name=job, ppn=max(1, nproc), queue_name=self.queue_name, errdir='perm',
                                               slurm_feature=self.slurm_feature, python_file=python_file, arguments='',
                                               python_executable=shlex.quote(sys.executable),
                                               slurm_exclude_directive=slurm_exclude_directive,
                                               walltime_directive=walltime_directive)
        else:
            logger.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
            logger.error('Exiting')
            sys.exit()

        qu_file = '{}{}'.format(job, constants.qext[self.queuing])
        if self.queuing == 'slurm':
            from kinbot.slurm_retry import wrap_script
            job_template = exclusive_slurm_script(
                job_template, self.par.get('queue_exclusive', 1))
            job_template = wrap_script(job_template, identity=job,
                prefixes=[job, os.path.join('perm', job)],
                exit_marker=f'{job}.exitcode', qc_job=job)
        with open(qu_file, 'w') as f_out_qu:
            f_out_qu.write(job_template)

        live_submit = None
        from kinbot import live_control
        if self.queuing == 'slurm' and live_control.find():
            from kinbot.qc_live import prepare_resources
            original_worker = Path(job+'.py').read_text()
            def live_submit():
                Path(qu_file).write_text(prepare_resources(self,job,job_template,original_worker))

        if self.queuing == 'slurm':
            # Slurm signal 53 occurs before the script starts when it cannot
            # open stdout/stderr. Catch missing directories, permissions,
            # quota, and inode exhaustion before consuming a scheduler job.
            output_directory = os.path.join('perm', os.path.dirname(job))
            os.makedirs(output_directory, exist_ok=True)
            probe = os.path.join(
                output_directory, f'.kinbot-write-probe-{os.getpid()}')
            try:
                with open(probe, 'w') as stream:
                    stream.write('KinBot scheduler output probe\n')
                    stream.flush()
                    os.fsync(stream.fileno())
            except OSError as error:
                raise RuntimeError(
                    f'Cannot write Slurm output directory {output_directory}: '
                    f'{error}. Check storage quota, inode quota, and permissions.') from error
            finally:
                try:
                    os.unlink(probe)
                except OSError:
                    pass

        command = [constants.qsubmit[self.queuing],
                   job + constants.qext[self.queuing]]
        failure_started = None
        warned = False
        while True:
            try:
                if shared_pool:
                    response = run_policy.submit(
                        command, wait=True, capture_output=True, check=False,
                        timeout=_SCHEDULER_SUBMIT_TIMEOUT,
                        **({'before_submit':live_submit} if live_submit else {}))
                    raw_out, raw_err = response.stdout, response.stderr
                    returncode = response.returncode
                else:
                    if live_submit:
                        live_submit()
                    process = subprocess.Popen(
                        command, shell=False, stdout=subprocess.PIPE,
                        stdin=subprocess.PIPE, stderr=subprocess.PIPE)
                    try:
                        raw_out, raw_err = process.communicate(
                            timeout=_SCHEDULER_SUBMIT_TIMEOUT)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.communicate()
                        raise
                    returncode = getattr(process, 'returncode', 0)
            except subprocess.TimeoutExpired as error:
                raise RuntimeError(
                    f'Scheduler submission timed out for {job}; check '
                    'the queue before retrying because the submission '
                    'outcome is unknown.') from error
            except FileNotFoundError as error:
                raise RuntimeError(
                    f'Scheduler submission executable {command[0]!r} is '
                    'unavailable.') from error
            except OSError as error:
                out, err = '', str(error)
                returncode = 'launch_error'
                failed = True
            else:
                out = raw_out.decode(errors='replace')
                err = raw_err.decode(errors='replace')
                failed = isinstance(returncode, int) and returncode != 0
            detail = str(err).strip() or str(out).strip()
            if not failed:
                break
            if not _scheduler_unavailable(detail):
                raise RuntimeError(
                    f'Scheduler rejected {job} with status {returncode}: '
                    f'{detail or "no diagnostic"}')
            now = time.monotonic()
            if failure_started is None:
                failure_started = now
            if not warned:
                logger.warning(
                    'Scheduler submission is unavailable for %s; no job was '
                    'created and submission will be retried: %s', job, detail)
                warned = True
            automatic = bool(getattr(self, 'queue_job_root', None))
            if (not automatic and
                    now - failure_started >= _SCHEDULER_QUERY_FAILURE_TIMEOUT):
                raise RuntimeError(
                    f'Scheduler submission remained unavailable for 60 '
                    f'seconds; no {job} job was submitted.')
            time.sleep(5)
        if self.queuing == 'pbs':
            pid = out.split('\n', 1)[0].split('.', 1)[0].strip()
        elif self.queuing == 'slurm':
            match = re.search(r'\bSubmitted batch job\s+(\d+)\b', out)
            if match is None:
                raise RuntimeError(
                    f'Cannot parse scheduler job id for {job}. stdout: '
                    f'{out!r}; stderr: {err!r}')
            pid = match.group(1)
        if not pid or not str(pid).isdigit():
            raise RuntimeError(
                f'Cannot parse scheduler job id for {job}. stdout: {out!r}; '
                f'stderr: {err!r}')
        self.job_ids[job] = pid
        live_control.observe_qc(Path.cwd()/job,pid,'running')
        generations = getattr(self, 'job_generation_tagged', None)
        if generations is None:
            generations = {}
            self.job_generation_tagged = generations
        generations[job] = True

        now = datetime.now()
        logger.debug(f'SUBMITTED {job} on {now.ctime()}')
        return 1  # important to keep it 1, this is the natural counter of jobs submitted

    def get_qc_geom(self,
                    job,
                    natom,
                    wait=0,
                    allow_error=0,
                    previous=0,
                    reorder: bool = False):
        '''
        Get the geometry from the ase database file.
        Returns it, with the following conditions about the status of the job.
        If wait = 0, return an (1, empty array) when the job is still running (instead of the geometry).
        If wait = 1, wait for the job to finish.
        If wait = 2, return the last geometry while the job is still running.
            This option is to monitor the formation of bimolecular products.
        If allow_error = 0, do not read the final geometry if the status is not 'normal'
        if allow_error = 1, read the geometry even though there is an error in the output file
            This option is to read the final IRC geometry when it did not converge
        If previous = 0, read the last geometry, this is the normal behaviour
        if previous = 1, read the geometry before the last one, this is needed in scan types so
            that the max energy point is taken, not the one after that
        '''
        geom = np.zeros((natom, 3))
        atoms = np.full(natom, 'H')
        check = self.check_qc(job)
        if check == 'error' and not allow_error:
            if reorder:
                return -1, geom, atoms
            else:
                return -1, geom
        status = 0
        while 1:
            check = self.check_qc(job)
            if check == 'running':
                if wait == 1:
                    time.sleep(1)
                elif wait == 2:
                    status = 2
                    break
                else:
                    if reorder:
                        return 1, geom, atoms
                    else:
                        return 1, geom
            else:
                break
        if check != 'normal':
            if not allow_error:
                if wait != 2:
                    if reorder:
                        return -1, geom, atoms
                    else:
                        return -1, geom

        # open the database
        finite_entries = []
        for row in self.db.select(name=job):
            mol = row.toatoms()
            candidate = np.asarray(mol.positions, dtype=float)
            if (candidate.shape == (natom, 3)
                    and np.isfinite(candidate).all()):
                finite_entries.append((candidate, mol.symbols))
            else:
                logger.warning('%s: ignoring a database row with a nonfinite '
                               'or malformed geometry.', job)

        if finite_entries and previous == 0:
            geom, atoms = finite_entries[-1]
            if reorder:
                return status, geom, atoms
            else:
                return status, geom
        elif len(finite_entries) >= 2 and previous == 1:
            return status, finite_entries[-2][0]
        else:
            if reorder:
                return -1, geom, atoms
            else:
                return -1, geom

    def get_qc_freq(self, job, natom, wait=0, allow_error=0):
        '''
        Get the frequencies from the ase database file
        If wait is set to 1, it will wait for the job to finish.
        '''

        check = self.check_qc(job)
        if check == 'error':
            return -1, [0]
        while 1:
            check = self.check_qc(job)
            if check == 'running':
                if wait == 1:
                    time.sleep(1)
                else:
                    return 1, []
            else:
                break

        if check != 'normal':
            return -1, [0]

        freq = []

        # open the database
        rows = self.db.select(name=job)

        # take the last entry
        for row in rows:
            if hasattr(row, 'data'):
                if not row.data.get('frequencies') is None:
                    freq = list(row.data.get('frequencies'))

        if len(freq) == 0 and natom > 1:
            return -1, freq

        return 0, freq

    def get_qc_energy(self, job, wait=0):
        '''
        Read the last energy from the db for job.
        For Gaussian currently works for DFT and HF only.
        For NWChem it works for optimization jobs, using the @ handle.
        If wait is set to 1, it will wait for the job to finish, otherwise
        just reads the last energy in the file.
        Returns the error code and the energy
        Error codes:
        -1: error
         0: normal and done
         1: running
        '''

        check = self.check_qc(job)
        if check == 'error':
            return -1, 0.
        while 1:
            check = self.check_qc(job)
            if check == 'running':
                if wait == 1:
                    time.sleep(1)
                else:
                    return 1, 0.
            else:
                break

        # Get last entry
        *_, last_row = self.db.select(name=job)
        if hasattr(last_row, 'data') and 'energy' in last_row.data:
            energy = last_row.data.get('energy')
        else:
            logger.warning(f'No energy found in the database for {job}. '
                           'This will lead to erroneous energies.')
            energy = 0

        energy *= constants.EVtoHARTREE

        return 0, energy

    def get_qc_zpe(self, job, wait=1):
        '''
        Read the zero point energy from the db for job.
        If wait is set to 1 (default), it will wait for the job to finish.
        '''

        check = self.check_qc(job)
        if check == 'error':
            return -1, 0.
        while 1:
            check = self.check_qc(job)
            if check == 'running':
                if wait == 1:
                    time.sleep(1)
                else:
                    return 0, 0.
            else:
                break
        zpe = 0.0  # set as default
        # open the database
        *_, last_row = self.db.select(name=job)
        # take the last entry
        
        if hasattr(last_row, 'data'):
            zpe = last_row.data.get('zpe')
        else:
            zpe = 0.0
            logger.warning('{} has no zpe in database. ZPE SET TO 0.0'.format(job))

        if zpe is None:
            zpe = 0.00

        return 0, zpe

    def hessian_is_massweighted(self):
        """Native QChem Hessians are weighted; ASE/Sella Hessians are not."""
        return (self.qc == 'qchem'
                and (not self.use_sella
                     or getattr(self, 'frequency_mode', 'auto')
                     == 'native_hessian'))

    def read_qc_hess(self, job, natom):
        '''
        Read the hessian of a gaussian chk file
        '''

        check = self.check_qc(job)
        if check != 'normal':
            return []

        rows = list(self.db.select(name=job))
        if (not (self.qc == 'qchem' and not self.use_sella)
                and rows and rows[-1].data.get('hess') is not None):
            return np.asarray(rows[-1].data['hess'])

        if (self.use_sella and self.frequency_mode != 'native_hessian') \
                or self.qc == 'fc':
            #db = connect('kinbot.db')
            for row in self.db.select(name=job):
                last_row = row
            hess = last_row.data['hess']
        elif self.qc == 'gauss':
            from kinbot.gaussian_hessian import read_formatted_hessian
            if not rows or len(rows[-1].symbols) != natom:
                raise ValueError('Gaussian Hessian has no matching selected geometry row.')
            fchk = str(job) + '.fchk'
            checkpoint = str(job) + '.chk'
            if (not os.path.exists(fchk) or (os.path.exists(checkpoint)
                    and os.path.getmtime(checkpoint) > os.path.getmtime(fchk))):
                if not os.path.exists(checkpoint):
                    return []
                subprocess.run(['formchk', checkpoint, fchk],
                               check=True, stdout=subprocess.DEVNULL)
            hess = read_formatted_hessian(fchk, rows[-1].symbols, rows[-1].positions)
        elif self.qc == 'qchem':
            hess = []
            row = 0
            do_read = False
            if natom == 1:
                return np.zeros([3, 3])
            if not os.path.exists(job + '_freq.out'):
                return []
            with open(job + '_freq.out') as f:
                for line in f:
                    if not do_read and line.startswith(' Mass-Weighted Hessian Matrix'):
                        do_read = True
                    elif do_read and 'Translations and Rotations' in line:
                        do_read = False
                    elif do_read:
                        if len(line) > 5:
                            if len(hess) < 3 * natom:
                                hess.append([])
                            hess[row].extend([float(val) for val in line.split()])
                            row += 1
                        else:
                            row = 0
                    else:
                        continue
            hess = np.array(hess)
        else:
            raise NotImplementedError()
        return hess

    def is_in_database(self, job):
        '''
        Checks if the current job is in the database:
        '''
        # open the database
        logger.debug(f'Looking for {job} in the db.')
        rows = self.db.select(name=job)

        mol = None

        # take the last entry
        for row in rows:
            mol = row.toatoms()

        if mol is None:
            return 0

        return 1

    def read_imag_mode(self, job, natom):
        """Read a TS mode from the artifact produced by this backend.

        Sella force-difference Hessians publish ``*_vib.xyz``. A Sella
        geometry optimization followed by a native Gaussian or Q-Chem
        Hessian publishes the program's native frequency output instead.
        Treating every Sella-optimized job as an ASE vibration calculation
        therefore asks for an artifact that native-Hessian workers neither
        need nor create.
        """
        if (self.qc == 'fc'
                or (self.use_sella
                    and self.frequency_mode != 'native_hessian')):
            from kinbot import reader_sella
            return reader_sella.read_imag_mode(job, natom)
        if self.qc == 'gauss':
            from kinbot import reader_gauss
            return reader_gauss.read_imag_mode(job, natom)
        if self.qc == 'qchem':
            from kinbot import reader_qchem
            return reader_qchem.read_imag_mode(job, natom)
        raise ValueError(f'Imaginary-mode reader is unavailable for {self.qc}.')

    def ingest_pkl(self, job, max_wait=60):
        '''
        Loads the {job}.pkl handoff written by a finished compute job.
        Returns the unpickled dict, or None if the file is not there.
        The file can disappear at any moment (the job's own start-up removes
        stale pkls, another driver may ingest it, NFS may lag), so a missing
        file is never an error here. A partially written file is retried
        for up to max_wait seconds, then skipped with a warning.
        '''
        for _ in range(max_wait):
            try:
                with open(f'{job}.pkl', 'rb') as f:
                    return pickle.load(f)
            except FileNotFoundError:
                return None
            except (EOFError, pickle.UnpicklingError):
                time.sleep(1)  # still being written
        logger.warning(f'{job}.pkl could not be read after {max_wait} s, '
                       'skipping it.')
        return None

    def check_qc(self, job):
        """Reconcile the current generation and optionally recover a timeout."""
        from kinbot.pes_phase import checkpoint
        checkpoint()
        from kinbot import qc_timeout, timeout_retry
        if self.queuing == 'slurm':
            waiting = qc_timeout.resume(self, job)
            if waiting is not None:
                return waiting
        status = QuantumChemistry._check_qc_with_grace(self, job)
        if self.queuing == 'slurm' and job in self.job_ids:
            if status == 'running':
                timeout_retry.observe(Path.cwd(), self.job_ids[job])
            elif status == 'error':
                replacement = qc_timeout.recover(self, job)
                if replacement is not None:
                    return replacement
        if self.queuing == 'slurm':
            from kinbot.qc_live import process
            status = process(self, job, status)
        return status

    def _check_qc_with_grace(self, job):
        """Read job status, allowing 60 s for a submitted job's final files."""
        status = QuantumChemistry._check_qc(self, job)
        if self.queuing not in ('slurm', 'pbs') or job not in self.job_ids:
            return status

        # Each submission gets one deadline. None marks a completed or expired
        # wait, so later checks cannot start another minute for the same job.
        if not hasattr(self, '_completion_deadlines'):
            self._completion_deadlines = {}
        key = (job, self.job_ids[job])
        if status == 'running':
            self._completion_deadlines.pop(key, None)
            return status
        if status in ('normal', 'error'):
            self._completion_deadlines[key] = None
            return status
        if status != 0:
            return status

        record = (_slurm_job_record(self.job_ids[job])
                  if self.queuing == 'slurm' else {})
        scheduler_state = record.get('JobState', '').split('+', 1)[0]
        if scheduler_state in _SLURM_TERMINAL_STATES - {'COMPLETED'}:
            exit_code = record.get('ExitCode', 'unknown')
            reason = record.get('Reason', 'unknown')
            node = record.get('NodeList', '')
            if exit_code.endswith(':53'):
                _save_slurm_node_failure(
                    self, node, job, self.job_ids[job], record)
                logger.error(
                    '%s: Slurm could not initialize batch stdout/stderr on %s '
                    '(job %s, state %s, reason %s, exit %s). The batch script '
                    'never started. Check storage quota and compute-node '
                    'filesystem access; %s is excluded from later submissions '
                    'in this run.', job, node or 'an unknown node',
                    self.job_ids[job], scheduler_state, reason, exit_code,
                    node or 'the failed node')
            else:
                logger.error(
                    '%s: Slurm job %s ended before publishing a calculation '
                    'result (state %s, reason %s, exit %s, node %s).',
                    job, self.job_ids[job], scheduler_state, reason,
                    exit_code, node or 'unknown')
            return 'error'

        now = time.monotonic()
        if key not in self._completion_deadlines:
            self._completion_deadlines[key] = now + 60.
            logger.info(f'{job}: job {key[1]} has left the queue; waiting up to '
                        '60 s for its complete result.')
        deadline = self._completion_deadlines[key]
        if deadline is not None:
            if now < deadline:
                return 'running'
            self._completion_deadlines[key] = None
            logger.warning(f'{job}: complete result is still unavailable '
                           '60 s after the job left the queue.')
        # A worker that exited without publishing a normal/error calculation
        # record is a failed job.  Returning zero here leaves Optimize polling
        # forever and hides the scheduler failure from restart logic.
        exit_file = f'{job}.exitcode'
        try:
            worker_exit = _worker_exit_status(exit_file, key[1])
        except (OSError, ValueError):
            logger.error(f'{job}: invalid worker status file {exit_file}.')
        else:
            if worker_exit is None:
                logger.error(f'{job}: scheduler job {key[1]} left the queue '
                             f'without a matching {exit_file} or a calculation '
                             'result. Inspect its scheduler stderr; the job '
                             'cannot be treated as complete.')
            elif worker_exit[0]:
                logger.error(f'{job}: worker exited with status '
                             f'{worker_exit[0]}; inspect the scheduler stderr.')
            else:
                logger.error(f'{job}: worker exited successfully but did not '
                             'publish a calculation result.')
        return 'error'

    def _database_qc_status(self, job, worker_finished=False):
        """Return a complete result status, or ``None`` if incomplete.

        The scheduler wrapper writes ``.exitcode`` only after the worker has
        returned.  A zero exit status therefore makes the last database row a
        finalized handoff even when an older/reused calculation lacks the
        backend's legacy trailing ``done`` marker.  Before that boundary the
        marker remains necessary so a partly visible network-filesystem write
        cannot be accepted as complete.
        """
        logger.debug('Checking for job %s in db', job)
        if not self.is_in_database(job):
            logger.debug('job %s is not in database', job)
            return None

        logger.debug('%s is in db', job)
        if (getattr(self, 'use_sella', False)
                and self.qc != 'nn_pes'
                and '_freq_recovery_' not in job
                and not job.startswith('vrctst/')):
            log_file = job + '_sella.log'
        elif self.qc == 'gauss':
            log_file = job + '.log'
        elif self.qc in ('nwchem', 'qchem'):
            log_file = job + '.out'
        elif self.qc in ('orca', 'fc'):
            log_file = job + '_sella.log'
        elif self.qc == 'nn_pes':
            log_file = None
        else:
            raise ValueError('Unknown code')

        if log_file is not None and not worker_finished:
            if not os.path.exists(log_file):
                if (self.queuing == 'local'
                        and not self.par['error_missing_local']):
                    return None
                # A database row and its filesystem artifacts can become
                # visible at slightly different times on network storage.
                time.sleep(1)
                if not os.path.exists(log_file):
                    logger.debug('Log file %s is not present.', log_file)
                    return None
            with open(log_file) as stream:
                lines = stream.readlines()
            if not lines:
                logger.debug('Log file %s is present, but it is empty.',
                             log_file)
                return None
            if 'done' not in lines[-1]:
                logger.debug('Log file %s is present, but has no "done" '
                             'stamp.', log_file)
                return None

        data = None
        last_row = None
        for row in self.db.select(name=job):
            if hasattr(row, 'data'):
                data = row.data
                last_row = row
        if data is None or data.get('status') not in ('normal', 'error'):
            logger.debug('Complete data is not in database for %s.', job)
            return None
        if data['status'] == 'normal':
            positions = np.asarray(last_row.toatoms().positions, dtype=float)
            energy = data.get('energy')
            try:
                finite_energy = (energy is None or np.isfinite(float(energy)))
            except (TypeError, ValueError):
                finite_energy = False
            if (positions.ndim != 2 or positions.shape[1:] != (3,)
                    or not np.isfinite(positions).all()
                    or not finite_energy):
                logger.error('%s: rejecting a nominally normal database '
                             'result with nonfinite geometry or energy.', job)
                return 'error'
        # A legacy Gaussian VRC worker publishes two database records.  The
        # relaxed record used to become visible before the frozen single
        # point was written, so the controller could advance directly into
        # correction-energy assembly and fail on the missing ``_fr`` row.
        # Treat both records as one atomic worker handoff.  When the worker
        # has exited, ``worker_finished`` also lets historical runs without
        # trailing log stamps be recovered from their completed database
        # rows, consistently with ordinary KinBot jobs.
        if _vrc_scan_requires_frozen_handoff(job):
            frozen = QuantumChemistry._database_qc_status(
                self, job + '_fr', worker_finished=worker_finished)
            if frozen is None:
                logger.debug(
                    'Waiting for paired frozen VRC result %s_fr.', job)
                return None
            if frozen != 'normal':
                logger.error('%s_fr did not produce an accepted frozen VRC '
                             'geometry.', job)
                return 'error'
        logger.debug('Returning status %s', data['status'])
        return data['status']

    def _check_qc(self, job):
        '''
        Checks the status of the qc job.
        Possible returns:
        running - the job is either running or is in the queue
        status - this can be normal or error, read from the database.
                 Only happens if log file had a done stamp and it was in the db.
        0 - job is not in the db or log file is not there with a done stamp or both.
            ==> this one resets the step number to 0
        '''

        # A reaction path deliberately reuses one logical job name for its
        # successive constrained optimizations.  A replacement submission
        # can therefore coexist briefly with the previous step's exit marker
        # or database row on a shared filesystem.  Check the scheduler-owned
        # generation first: while that exact Slurm job is pending or running,
        # historical artifacts must not advance the path or trigger another
        # submission of the same script.
        current_job_id = (self.job_ids.get(job)
                          if self.queuing in ('slurm', 'pbs') else None)
        receipt_finished = False
        generation_tagged = getattr(
            self, 'job_generation_tagged', {}).get(job, True)
        if current_job_id and generation_tagged:
            # The worker publishes its generation-tagged marker atomically and
            # only after it has finished publishing its calculation artifacts.
            # Prefer that positive handoff over another scheduler round trip.
            # A one-field legacy marker is deliberately not authoritative here:
            # it cannot distinguish a replacement constrained-path submission
            # from the preceding generation that used the same logical name.
            try:
                current_exit = _worker_exit_status(
                    f'{job}.exitcode', current_job_id)
            except (OSError, ValueError):
                logger.error('Invalid worker status file %s.exitcode.', job)
                return 'error'
            current_finished = bool(current_exit and current_exit[1])
            if current_finished and current_exit[0]:
                logger.error('%s exited with worker status %s.',
                             job, current_exit[0])
                return 'error'

            active = (False if current_finished else
                      (_slurm_job_active(current_job_id)
                       if self.queuing == 'slurm'
                       else _pbs_job_active(current_job_id)))
            if active is None:
                warned = getattr(self, '_scheduler_query_warnings', set())
                key = (self.queuing, str(current_job_id))
                if key not in warned:
                    scheduler_name = ('Slurm' if self.queuing == 'slurm'
                                      else 'PBS')
                    logger.warning(
                        '%s: %s status query is unavailable; retaining the '
                        'submitted job state.', job, scheduler_name)
                    warned.add(key)
                    self._scheduler_query_warnings = warned
                return 'running'
            if active:
                logger.debug('Current %s generation of %s is active.',
                             self.queuing.upper(), job)
                return 'running'

            # The scheduler has released this generation. Do not ingest a pkl
            # or database row until its matching worker marker is visible;
            # otherwise an older constrained-path step can be mistaken for
            # the replacement on a lagging shared filesystem.
            if not current_finished:
                # Older controllers could delete the completed generation's
                # exit marker while waiting to submit its next constrained
                # step. The allocation-specific wrapper receipt is a second
                # completion record, but only after Slurm confirms inactivity.
                if self.queuing == 'slurm' and not os.path.exists(f'{job}.exitcode'):
                    from kinbot.slurm_retry import completion_receipt
                    receipt = completion_receipt(job, current_job_id)
                    if receipt is not None:
                        if receipt[0]:
                            return 'error'
                        receipt_finished = True
                if not receipt_finished:
                    return 0

        loaded_data = self.ingest_pkl(job)
        if loaded_data is not None:
            mol = Atoms(symbols=loaded_data['sym'], positions=loaded_data['pos'])
            name = loaded_data['name']
            data = loaded_data['data']
            count0 = sum(1 for _ in self.db.select(name=name))
            self.db.write(mol, name=name, data=data)
            # wait for the row to become visible, but do not hang forever
            for _ in range(60):
                count1 = sum(1 for _ in self.db.select(name=name))
                if count1 == count0 + 1:
                    break
                time.sleep(1)
            else:
                logger.warning(f'Row for {name} was not visible in kinbot.db '
                               '60 s after writing it.')
            # the pkl is a one-way handoff from the job, not needed after this
            # it may already be gone: the job's own start-up removes stale pkls,
            # a previous driver may have ingested it, or the filesystem lags
            try:
                os.remove(f'{job}.pkl')
            except FileNotFoundError:
                pass

        logger.debug('Checking job {}'.format(job))
        exit_file = f'{job}.exitcode'
        worker_finished = receipt_finished
        if not receipt_finished and os.path.isfile(exit_file):
            try:
                worker_exit = _worker_exit_status(exit_file, current_job_id)
            except (OSError, ValueError):
                logger.error(f'Invalid worker status file {exit_file}.')
                return 'error'
            if worker_exit is None:
                return 0
            if worker_exit[0]:
                logger.error(f'{job} exited with worker status '
                             f'{worker_exit[0]}.')
                return 'error'
            worker_finished = True

        database_status = QuantumChemistry._database_qc_status(
            self, job, worker_finished=worker_finished)
        if database_status is not None:
            if receipt_finished and database_status == 'normal':
                notices = getattr(self, '_receipt_completion_notices', set())
                key = (job, str(current_job_id))
                if key not in notices:
                    logger.info('%s: recovered missing completion handoff for Slurm job %s '
                                'from its completed wrapper receipt and normal QC result.',
                                job, current_job_id)
                    notices.add(key)
                    self._receipt_completion_notices = notices
            return database_status

        # The batch wrapper publishes the exit file only after the worker has
        # returned. Once it says zero, read the native result directly: a
        # scheduler query cannot add useful information and may itself stall.
        if not worker_finished:
            if self.queuing == 'pbs':
                job_id = self.job_ids.get(job)
                if job_id:
                    try:
                        result = subprocess.run(
                            ['qstat', '-f', str(job_id)], capture_output=True,
                            text=True, check=False,
                            timeout=_SCHEDULER_QUERY_TIMEOUT)
                    except (OSError, subprocess.TimeoutExpired):
                        result = None
                    if result is None:
                        warned = getattr(
                            self, '_scheduler_query_warnings', set())
                        key = ('pbs', str(job_id))
                        if key not in warned:
                            logger.warning(
                                '%s: PBS status query is unavailable; '
                                'retaining the submitted job state.', job)
                            warned.add(key)
                            self._scheduler_query_warnings = warned
                        return 'running'
                    if result.returncode == 0:
                        logger.debug('Job is running')
                        return 'running'
            elif self.queuing == 'slurm':
                job_id = self.job_ids.get(job)
                if job_id:
                    active = _slurm_job_active(job_id)
                    if active is None:
                        warned = getattr(
                            self, '_scheduler_query_warnings', set())
                        key = ('slurm', str(job_id))
                        if key not in warned:
                            logger.warning(
                                '%s: Slurm status query is unavailable; '
                                'retaining the submitted job state.', job)
                            warned.add(key)
                            self._scheduler_query_warnings = warned
                        return 'running'
                    if active:
                        logger.debug('Job is running')
                        return 'running'
            elif self.queuing == 'local':
                pass
            else:
                logger.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
                logger.error('Exiting')
                sys.exit()

        if worker_finished and current_job_id:
            # Let check_qc apply its bounded shared-filesystem grace period.
            return 0
        if worker_finished:
            # No scheduler transition can repair a worker that has already
            # exited without publishing a database result.  Returning zero
            # here made Optimize poll such jobs forever, especially after a
            # restart where the historical Slurm id is intentionally absent.
            logger.error('%s: worker exited successfully but did not publish '
                         'a calculation result.', job)
            return 'error'
        if self.queuing == 'local' and not self.par['error_missing_local']:
            logger.debug('local qu error')
            return 'error'
        return 0

    def limit_jobs(self):
        '''
        Check how many jobs are in the queue from the user, and if larger than the limit,
        then wait for resources to free up.
        '''

        warned = False
        failure_started = None
        while 1:
            if self.queuing == 'slurm':
                command = ['squeue', '-h', '-u', '{}'.format(self.username)]
                if getattr(self, 'queue_job_root', None):
                    command.extend(('-o', '%Z'))
            elif self.queuing == 'pbs':
                command = ['qselect', '-u', '{}'.format(self.username)]
            elif self.queuing == 'local':
                command = command = ['echo', '']
            try:
                result = subprocess.run(
                    command, capture_output=True, text=True, check=False,
                    timeout=_SCHEDULER_QUERY_TIMEOUT)
            except (OSError, subprocess.TimeoutExpired) as error:
                now = time.monotonic()
                if failure_started is None:
                    failure_started = now
                if not warned:
                    logger.warning('Scheduler queue query is unavailable; '
                                   'waiting before job submission: %s', error)
                    warned = True
                if (not getattr(self, 'queue_job_root', None)
                        and now - failure_started
                        >= _SCHEDULER_QUERY_FAILURE_TIMEOUT):
                    raise RuntimeError(
                        'Scheduler queue query remained unavailable for 60 '
                        'seconds; no job was submitted.') from error
                time.sleep(5)
                continue
            if result.returncode:
                now = time.monotonic()
                if failure_started is None:
                    failure_started = now
                detail = result.stderr.strip() or str(result.returncode)
                if not warned:
                    logger.warning('Scheduler queue query failed; waiting '
                                   'before job submission: %s',
                                   detail)
                    warned = True
                automatic_outage = (
                    bool(getattr(self, 'queue_job_root', None))
                    and _scheduler_unavailable(detail))
                if (not automatic_outage and now - failure_started
                        >= _SCHEDULER_QUERY_FAILURE_TIMEOUT):
                    raise RuntimeError(
                        'Scheduler queue query failed continuously for 60 '
                        'seconds; no job was submitted. Last error: ' +
                        detail)
                time.sleep(5)
                continue
            warned = False
            failure_started = None
            jobs = [line for line in result.stdout.splitlines() if line.strip()]
            queue_root = getattr(self, 'queue_job_root', None)
            if self.queuing == 'slurm' and queue_root:
                recursive = getattr(
                    self, 'queue_job_root_recursive', False)
                selected = []
                for value in jobs:
                    try:
                        workdir = os.path.realpath(value.strip())
                        common = os.path.commonpath((queue_root, workdir))
                    except (OSError, ValueError):
                        continue
                    if (workdir == queue_root
                            or (recursive and common == queue_root)):
                        selected.append(value)
                jobs = selected
            if len(jobs) < self.queue_job_limit:
                return 0
            time.sleep(1)

    def add_dummy(self, spatom, geom, spbond):
        '''
        Add a dummy atom for each close to linear angle.
        Obsolete.
        '''
        atom = copy.deepcopy(spatom)
        dummy = geometry.is_linear(geom, spbond)
        if len(dummy) > 0:
            for d in dummy:
                atom = np.append(atom, ['X'])
                geom = np.concatenate((geom, [d]), axis=0)
        dummy = [d.tolist() for d in dummy]
        return atom, geom, dummy

    def merge_kwargs(self, kwargs, user_kwargs=None):
        if user_kwargs is None:
            user_kwargs = self.par['calc_kwargs']
        if self.qc == 'gauss':
            from kinbot.gaussian_input import merge_gaussian_kwargs
            return merge_gaussian_kwargs(kwargs, user_kwargs)
        for key, uvalue in user_kwargs.items():
            if key in kwargs:
                user_val_dict, val_dict = {}, {}
                for field in uvalue.split(','):
                    if len(field.split('=')) == 1:
                        user_val_dict[field] = None
                    elif len(field.split('=')) == 2:
                        user_val_dict[field.split('=')[0]] = field.split('=')[1]
                    else:
                        raise ValueError(f'Unexpected value: {uvalue}')
                for field in kwargs[key].split(','):
                    if len(field.split('=')) == 1:
                        val_dict[field] = None
                    elif len(field.split('=')) == 2:
                        val_dict[field.split('=')[0]] = field.split('=')[1]
                    else:
                        raise ValueError(f'Unexpected value: {uvalue}')
                val_dict.update(user_val_dict)
                kwargs[key] = ','.join([f'{k}={v}' if v is not None else k for k, v in val_dict.items() ])
                
            else:
                kwargs[key] = uvalue
        return kwargs
