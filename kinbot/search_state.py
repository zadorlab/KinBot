"""Durable reaction-search progress, independent of the last QC exit status.

A successful constrained QC point is not a successful reaction search. Keep
controller decisions and path steps separately; rebuild higher-stage Python
objects from QC records without reopening terminal searches.
"""
from datetime import datetime, timezone
import ast
import hashlib
import json
from pathlib import Path

from kinbot.execution import atomic_json, file_sha256


NAME = '.kinbot_search_state.json'


def _read(path):
    return json.loads(path.read_text()) if path.is_file() else {}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                    separators=(',', ':')).encode()).hexdigest()


def contract(inp):
    inp = Path(inp).resolve()
    if not inp.is_file():
        raise FileNotFoundError(f'Reaction-search input is missing: {inp}')
    routes = _read(inp.parent/'.kinbot_theory_jobs.json')
    return dict(input_sha256=_digest(_read(inp)),
                profile_sha256=routes.get('profile_sha256'))


def validate_input(inp):
    """Check restart compatibility before starting or rewriting a worker."""
    inp = Path(inp).resolve()
    record = _read(inp.parent/NAME)
    if record and (record.get('schema') != 1 or record.get('contract') != contract(inp)):
        raise RuntimeError(f'{inp}: reaction-search checkpoint inputs/theory changed. '
                           'Restore the recorded input or use an explicit search revision; '
                           'no search was restarted.')


def input_restore_plan(inp, logs=()):
    """Recover only a regenerated embedding, proven by the full input hash.

    Old checkpoints did not retain input snapshots. A structure printed by the
    old worker is usable only when substituting it reproduces the *entire*
    recorded input. No theory, model, policy or search decision is rebased.
    This function is read-only; refresh applies the plan after owners stop.
    """
    inp = Path(inp).resolve()
    record = _read(inp.parent/NAME)
    current = contract(inp)
    if not record or record.get('contract') == current:
        validate_input(inp)
        return None
    expected = record.get('contract', {})
    if (record.get('schema') != 1
            or expected.get('profile_sha256') != current['profile_sha256']):
        validate_input(inp)
    values = _read(inp)
    candidates = []
    if 'input_snapshot' in record:
        candidates.append(('checkpoint input snapshot', record['input_snapshot']))
    for log in logs:
        log = Path(log)
        if not log.is_file():
            continue
        # Parameters precede native progress; never read a whole QC log or
        # parse/evaluate arbitrary Python from retained evidence.
        with log.open(errors='replace') as stream:
            reading = False
            for line in stream.read(2_000_000).splitlines():
                if 'Input parameters' in line:
                    reading = True
                elif reading and 'Starting KinBot' in line:
                    reading = False
                elif reading and line.strip().startswith('structure: '):
                    try:
                        structure = ast.literal_eval(line.strip().partition(': ')[2])
                    except (ValueError, SyntaxError):
                        continue
                    candidates.append((str(log), dict(values, structure=structure)))
    for source, candidate in candidates:
        if (isinstance(candidate, dict)
                and _digest(candidate) == expected.get('input_sha256')
                and {k for k in candidate.keys() | values.keys()
                     if candidate.get(k) != values.get(k)} == {'structure'}):
            return dict(before=current, after=expected, input=candidate, source=source,
                        reason='Restore original PES embedding with exact checkpoint input hash')
    validate_input(inp)  # No verified recovery: preserve the strict mismatch.


def restore_input(inp, plan, archive):
    """Apply a verified plan under refresh's stopped-controller locks."""
    inp, archive = Path(inp), Path(archive)
    current = contract(inp)
    record = _read(inp.parent/NAME)
    if record.get('contract') != plan['after'] or _digest(plan['input']) != plan['after']['input_sha256']:
        raise RuntimeError(f'{inp}: input recovery proof changed; nothing restored.')
    if current == plan['after']:
        return  # An interrupted handoff already restored this input.
    if current != plan['before']:
        raise RuntimeError(f'{inp}: input changed since recovery planning; nothing restored.')
    archive.parent.mkdir(parents=True, exist_ok=True)
    if archive.exists() and archive.read_bytes() != inp.read_bytes():
        raise RuntimeError(f'{archive}: conflicting recovery archive; retained.')
    if not archive.exists():
        archive.write_bytes(inp.read_bytes())
    atomic_json(inp, plan['input'])
    validate_input(inp)


def capture_legacy(inp, *, publish=True):
    """Import the old monitor once, with controllers stopped/frozen.

    Admission receipts distinguish a submission newer than the monitor from
    an input merely rendered while waiting for pool capacity. No QC is edited.
    """
    inp = Path(inp).resolve()
    path, monitor = inp.parent/NAME, inp.parent/'kinbot_monitor.out'
    if path.exists() or not monitor.is_file():
        return
    routes = _read(inp.parent/'.kinbot_theory_jobs.json').get('jobs', {})
    from kinbot.run_policy import find
    policy = find(inp.parent)
    entries = {}
    for line in monitor.read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        if len(fields) < 3:
            raise RuntimeError(f'{monitor}: incomplete search monitor; retained.')
        status, step = map(int, fields[:2])
        name = fields[2]
        if name in entries or status not in (-999, -1, 0, 1, 2, 3, 4, 5, 6):
            raise RuntimeError(f'{monitor}: invalid/duplicate search record for {name}.')
        job = routes.get(name, {}).get('job_id')
        entry = dict(status=status, step=step, job_id=job, legacy=True)
        if status == 0 and job:
            receipt = (_read(policy.parent/'.kinbot_policy_submissions'/(str(job)+'.json'))
                       if policy else {})
            created = receipt.get('created_at')
            if created and Path(receipt.get('directory', '')).resolve() == inp.parent:
                command = receipt.get('command', [])
                # Timeout/I/O retries resubmit the same point using another
                # batch script; they do not advance the reaction coordinate.
                path_submission = bool(command and Path(command[-1]).name == name+'.sbatch')
                entry['submission_after_monitor'] = (
                    path_submission and
                    datetime.fromisoformat(created).timestamp() > monitor.stat().st_mtime)
            else:
                # A frozen old controller can be one step ahead of its monitor.
                # Do not silently guess when no submission provenance exists.
                entry['unverified_step'] = True
        entries[name] = entry
    record = dict(schema=1, contract=contract(inp), input_snapshot=_read(inp), reactions=entries,
                  imported_monitor_sha256=file_sha256(monitor))
    if publish:
        atomic_json(path, record)
    return record


def rebase(inp, reason, *, reassess=False):
    """An explicit input migration may revalidate retained accepted results."""
    inp = Path(inp).resolve()
    path = inp.parent/NAME
    if not path.exists():
        return
    record = _read(path)
    record.setdefault('revisions', []).append(dict(
        reason=reason, contract=record['contract'], reassess=bool(reassess)))
    record['contract'] = contract(inp)
    record['input_snapshot'] = _read(inp)
    if reassess:
        # Terminal failures remain closed unless the bounded revisit control
        # names them. Accepted searches reconstruct and revalidate under policy.
        for entry in record['reactions'].values():
            if entry['status'] == -1:
                entry['reassess'] = True
    atomic_json(path, record)


def _backends(qc):
    return list(qc.backends.values()) if hasattr(qc, 'backends') else [qc]


def _job(qc, name):
    backend = qc._backend_for(name) if hasattr(qc, '_backend_for') else qc
    return getattr(backend, 'job_ids', {}).get(name)


class SearchState:
    def __init__(self, generator):
        self.generator = generator
        self.path = Path.cwd()/NAME
        self.source = contract(generator.inp)
        # Direct, single-owner restarts also import the pre-checkpoint monitor.
        capture_legacy(generator.inp)
        self.record = _read(self.path) or dict(schema=1, contract=self.source, reactions={})
        validate_input(generator.inp)
        self.record['input_snapshot'] = _read(Path(generator.inp))
        self.entries = self.record['reactions']
        self.replaying = set()
        self.scan_seen = {}
        self.restore()

    def write(self):
        atomic_json(self.path, self.record)

    def signature(self, obj):
        return dict(family=obj.family_name, instance=[int(x) for x in obj.instance],
                    max_step=int(obj.max_step), scan=bool(obj.scan))

    def restore(self):
        species, qc = self.generator.species, self.generator.qc
        for i, obj in enumerate(species.reac_obj):
            obj.search_checkpoint = self
            obj.search_checkpoint_index = i
            entry = self.entries.get(obj.instance_name)
            if entry is None:
                self.save(i)
                continue
            signature = self.signature(obj)
            if entry.get('signature', signature) != signature:
                raise RuntimeError(f'{obj.instance_name}: reaction definition changed; '
                                   'an explicit search revision is required.')
            entry['signature'] = signature
            status = entry['status']
            if status == -999:
                species.reac_ts_done[i] = -999
                species.reac_step[i] = entry['step']
                if entry.get('stereochemical_rejection'):
                    obj.stereochemical_rejection = entry['stereochemical_rejection']
            elif status == 0:
                step = entry['step']
                if entry.get('unverified_step'):
                    raise RuntimeError(f'{obj.instance_name}: old constrained-search step '
                                       'has no admission receipt; refusing to repeat or skip it.')
                if entry.get('legacy') and obj.scan and step > 1:
                    raise RuntimeError(f'{obj.instance_name}: legacy scan lacks its energy '
                                       'history; review is required before resuming.')
                if entry.pop('submission_after_monitor', False):
                    # Old non-scan families only skip the initial torsion steps.
                    if step == 0 and obj.skip and len(obj.instance) < 4:
                        step = 12
                    step += 1
                pending = entry.get('pending_submission')
                if pending:
                    if pending.get('submitted_job_id'):
                        step = pending['next_step']
                    else:
                        # A changed allocation alone can be an infrastructure
                        # retry of the preceding point, not a new path step.
                        step = pending['step']
                    entry.pop('pending_submission')
                species.reac_step[i] = step
                species.reac_scan_energy[i] = list(entry.get('scan_energy', []))
                self.scan_seen[obj.instance_name] = entry.get('scan_energy_job')
                entry['step'] = step
            else:
                # Stages 1+ need IRC/products/Optimize objects, not another TS
                # search. Homolytic scission constructs its products in stage 0.
                species.reac_ts_done[i] = 0 if species.reac_type[i] == 'hom_sci' else 1
                species.reac_step[i] = entry['step']
                if status == -1 and not entry.pop('reassess', False):
                    self.replaying.add(obj.instance_name)
            entry.pop('legacy', None)
        self.write()

    def enter(self, index):
        obj = self.generator.species.reac_obj[index]
        for backend in _backends(self.generator.qc):
            backend._search_replay = obj.instance_name if obj.instance_name in self.replaying else None

    def leave(self):
        for backend in _backends(self.generator.qc):
            backend._search_replay = None

    def adopt_revision(self):
        self.record = _read(self.path)
        self.source = self.record['contract']
        self.entries = self.record['reactions']

    def save(self, index):
        species = self.generator.species
        obj = species.reac_obj[index]
        name, status = obj.instance_name, int(species.reac_ts_done[index])
        previous = self.entries.get(name, {})
        # A second refresh during reconstruction must still see the terminal
        # result, rather than converting it to an unfinished search.
        if name in self.replaying and status >= 0:
            return
        entry = dict(status=status, step=int(species.reac_step[index]),
                     signature=self.signature(obj),
                     job_id=_job(self.generator.qc, name),
                     scan_energy=list(species.reac_scan_energy[index]),
                     scan_energy_job=self.scan_seen.get(name))
        if getattr(obj, 'stereochemical_rejection', None):
            entry['stereochemical_rejection'] = obj.stereochemical_rejection
        if previous.get('pending_submission'):
            entry['pending_submission'] = previous['pending_submission']
        if previous.get('history'):
            entry['history'] = previous['history']
        if entry != previous:
            self.entries[name] = entry
            self.write()
        if status < 0:
            self.replaying.discard(name)

    def scan_energy(self, index, energy):
        species = self.generator.species
        name = species.reac_obj[index].instance_name
        # Include the path step for non-profiled/reuse-only backends with no ID.
        key = [_job(self.generator.qc, name), int(species.reac_step[index])]
        if self.scan_seen.get(name) != key:
            species.reac_scan_energy[index].append(float(energy))
            self.scan_seen[name] = key
            self.save(index)

    def prepare_submission(self, obj, step):
        self.save(obj.search_checkpoint_index)
        entry = self.entries[obj.instance_name]
        entry['pending_submission'] = dict(step=int(step), next_step=int(step)+1,
                                          previous_job_id=_job(self.generator.qc, obj.instance_name))
        self.write()
        for backend in _backends(self.generator.qc):
            backend._reaction_path_checkpoint = (self, obj.instance_name)

    def submitted(self, name, job):
        entry = self.entries[name]
        pending = entry['pending_submission']
        pending['submitted_job_id'] = str(job)
        entry.update(step=pending['next_step'], job_id=str(job))
        self.write()

    def finish_submission(self, obj, step):
        entry = self.entries[obj.instance_name]
        entry.pop('pending_submission', None)
        entry['step'] = int(step)
        self.write()
        for backend in _backends(self.generator.qc):
            backend._reaction_path_checkpoint = None

    def reset(self, index, reason):
        """Only owning live controls call this after archiving the old attempt."""
        obj = self.generator.species.reac_obj[index]
        old = self.entries.pop(obj.instance_name, {})
        history = old.pop('history', [])
        self.replaying.discard(obj.instance_name)
        self.scan_seen.pop(obj.instance_name, None)
        self.generator.species.reac_scan_energy[index] = []
        obj.search_checkpoint = self
        obj.search_checkpoint_index = index
        self.save(index)
        self.entries[obj.instance_name]['history'] = history + [dict(
            reason=reason, previous=old, reset_at=datetime.now(timezone.utc).isoformat())]
        self.write()


def submitted(qc, job, job_id):
    pending = getattr(qc, '_reaction_path_checkpoint', None)
    if pending and pending[1] == job:
        pending[0].submitted(job, job_id)


def reset(generator, index, reason):
    checkpoint = getattr(generator, 'search_checkpoint', None)
    if checkpoint:
        checkpoint.reset(index, reason)
