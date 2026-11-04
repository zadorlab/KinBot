"""Share identical new Gaussian/Sella HIR and conformer attempts within a pool.

This is calculation reuse, not species/TS equivalence or population merging.
An immutable native bundle is copied to each consumer before its DB row appears.
Existing local attempts are never replaced and Slurm ownership is never aliased.
"""
import ast
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
from functools import lru_cache
import hashlib
from importlib.metadata import version, PackageNotFoundError
import json
import logging
import os
from pathlib import Path
import pickle
import shutil
import sqlite3
import sys
from types import SimpleNamespace
import tempfile
import uuid

import numpy as np

from kinbot.execution import atomic_json, file_sha256
from kinbot.spatial import proper_alignment

logger = logging.getLogger('KinBot')
NAME = '.kinbot_qc_reuse'
# Coordinate serialization noise, not a tolerance for different conformers.
GEOMETRY_TOLERANCE = 1.e-6
SUFFIXES = ('.log', '.com', '.chk', '.fchk', '_sella.log', '.traj', '.xyz',
            '_vib.xyz', '.result.json')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


@lru_cache(maxsize=1)
def runtime():
    h = hashlib.sha256()
    base = Path(__file__).parent
    for path in sorted(base.rglob('*.py')):
        h.update(str(path.relative_to(base)).encode())
        h.update(path.read_bytes())
    packages = {}
    for name in ('ase', 'numpy', 'sella'):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = None
    return dict(kinbot_sha256=h.hexdigest(), packages=packages,
                python=list(sys.version_info[:3]))


def definition(qc, job, *, path=None):
    """Parse literals without executing a worker; unsupported forms opt out."""
    if (getattr(qc, 'qc', None) != 'gauss' or not getattr(qc, 'use_sella', False)
            or not Path(job).parts or Path(job).parts[0] not in ('hir', 'conf')):
        return None
    try:
        raw = Path(path or job+'.py').read_text()
        tree = ast.parse(raw)
        atoms = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                 and isinstance(n.func, ast.Name) and n.func.id == 'Atoms']
        if len(atoms) != 1:
            return None
        values = {k.arg: ast.literal_eval(k.value) for k in atoms[0].keywords}
        symbols, positions = values['symbols'], np.asarray(values['positions'], float)
        if (positions.shape != (len(symbols), 3) or not len(symbols)
                or not np.isfinite(positions).all()):
            return None
        assigns = [n for n in tree.body if isinstance(n, ast.Assign)
                   and any(isinstance(t, ast.Name) and t.id == 'kwargs' for t in n.targets)]
        if len(assigns) != 1:
            return None
        kwargs = ast.literal_eval(assigns[0].value)
        # Fields, embedding, external basis/ECP files and supplied guesses
        # need additional frame/dependency evidence. They run normally.
        allowed = {'method', 'basis', 'charge', 'mult', 'Symm', 'symm', 'pop',
                   'scf', 'integral', 'EmpiricalDispersion', 'force',
                   'label', 'chk', 'mem', 'nprocshared'}
        if set(kwargs) - allowed or str(kwargs.get('basis', '')).lower() in ('gen', 'genecp', 'chkbasis'):
            return None
        # External starting wavefunctions and arbitrary user scripts require
        # an explicit dependency contract; they are not inferred from paths.
        if any(str(k).lower() in ('guess', 'oldchk', 'restart') for k in kwargs):
            return None
        if kwargs.get('chk') not in (None, job, Path(job).name):
            return None
        if any(not isinstance(x, str) for x in symbols):
            return None
        if not any(isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                   and n.func.id == 'run_optimizer' for n in ast.walk(tree)):
            return None
        reduced = {k: v for k, v in sorted(kwargs.items())
                   if k not in ('label', 'chk', 'mem', 'nprocshared')}
        assigns[0].value = ast.parse(repr(reduced), mode='eval').body
        for k in atoms[0].keywords:
            if k.arg == 'positions':
                k.value = ast.Constant('GEOMETRY')
        class Labels(ast.NodeTransformer):
            def visit_Constant(self, node):
                if isinstance(node.value, str):
                    for label, replacement in ((job, 'JOB'), (Path(job).name, 'BASE')):
                        if node.value == label or node.value.startswith((label+'.', label+'_')):
                            return ast.Constant(replacement+node.value[len(label):])
                return node
        tree = Labels().visit(tree)
        contract = dict(program=ast.dump(tree, include_attributes=False),
            command=getattr(qc, 'qc_command', None), runtime=runtime(),
            # Site setup/environment and user batch wrappers can change QC.
            environment={k: os.environ.get(k) for k in
                         ('g16root', 'GAUSS_EXEDIR', 'GAUSS_ARCHDIR', 'LOADEDMODULES',
                          'PATH', 'LD_LIBRARY_PATH', 'KINBOT_BACKEND')},
            queue_template=getattr(qc, 'par', {}).get('queue_template', ''))
        template = contract['queue_template']
        if template:
            contract['queue_template'] = file_sha256(template)
        return dict(contract_sha256=digest(contract), symbols=symbols,
                    positions=positions.tolist(), worker_sha256=hashlib.sha256(raw.encode()).hexdigest(),
                    program_sha256=digest(contract['program']))
    except (OSError, SyntaxError, ValueError, TypeError, KeyError):
        return None


def equivalent(a, b):
    if a['contract_sha256'] != b['contract_sha256'] or a['symbols'] != b['symbols']:
        return False
    fit = proper_alignment(a['positions'], b['positions'])
    return max(fit.errors) <= GEOMETRY_TOLERANCE


def _same_program(saved, job, path=None):
    current = definition(SimpleNamespace(qc="gauss", use_sella=True), job, path=path)
    return _program_matches(saved, current)


def _program_matches(saved, current):
    return (current is not None and current["program_sha256"] == saved["program_sha256"]
            and current["symbols"] == saved["symbols"]
            and np.max(proper_alignment(current["positions"], saved["positions"]).errors)
            <= GEOMETRY_TOLERANCE)


def _marker(job, role):
    return Path(job+f'.qc_reuse_{role}.json')


def _read(path):
    return json.loads(Path(path).read_text())


def _inside(root, name):
    path = (root/name).resolve()
    if not path.is_relative_to(root):
        raise ValueError('QC reuse path escapes the shared pool.')
    return path


@contextmanager
def _locked(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        yield


def _record(root, token):
    if not isinstance(token, str) or len(token) != 32 or any(c not in '0123456789abcdef' for c in token):
        raise ValueError('Invalid QC reuse record identity.')
    return root/NAME/'records'/(token+'.json')


def _scope():
    from kinbot.run_policy import find
    policy = find()
    return policy.parent.resolve() if policy else None


def submit(qc, job, singlejob, launch):
    """Atomically choose one new producer; followers never own its Slurm ID."""
    if (singlejob != 1 or getattr(qc, 'queuing', None) != 'slurm'
            or getattr(qc, 'qc', None) != 'gauss' or not getattr(qc, 'use_sella', False)
            or not Path(job).parts or Path(job).parts[0] not in ('hir', 'conf')):
        return launch()
    root = _scope()
    if root is None or getattr(qc, '_search_replay', None):
        return launch()
    # Keep already-submitted, completed, failed and partially published local
    # attempts authoritative. This also makes activation on a live run safe.
    consumer = _marker(job, 'consumer')
    if consumer.exists():
        status = poll(qc, job)
        if status is not None:
            return 0
    if (job in getattr(qc, 'job_ids', {}) or qc.is_in_database(job)
            or any(Path(job+s).exists() for s in (*SUFFIXES, '.pkl', '.exitcode', '.sbatch'))):
        return launch()
    wanted = definition(qc, job)
    if wanted is None:
        return launch()
    directory = Path.cwd().resolve()
    source = str((directory/job).relative_to(root))
    bucket = root/NAME/'index'/(wanted['contract_sha256']+'.json')
    token = None
    with _locked(bucket.with_suffix('.lock')):
        candidates = _read(bucket) if bucket.exists() else []
        for ident in candidates:
            entry = _read(_record(root, ident))
            if (entry['source'] == source and entry['state'] in ('reserved', 'submitted')
                    and equivalent(entry['definition'], wanted)):
                # The owning driver resumes its existing reservation. Its
                # ordinary scheduler ledger remains responsible for adoption.
                token = ident
                break
            if entry['state'] not in ('reserved', 'submitted', 'complete'):
                continue
            if equivalent(entry['definition'], wanted):
                atomic_json(consumer, dict(schema=1, root=str(root), token=ident,
                    definition=wanted, state='waiting'))
                logger.info('%s: waiting for equivalent QC %s; no duplicate submission.', job, entry['source'])
                return 0
        if token is None:
            from kinbot.pes_phase import process_identity
            token = uuid.uuid4().hex
            atomic_json(_record(root, token), dict(schema=1, token=token, source=source,
                directory=str(directory.relative_to(root)), job=job, definition=wanted,
                state='reserved', owner=process_identity(os.getpid()),
                created_at=datetime.now(timezone.utc).isoformat()))
            candidates.append(token)
            atomic_json(bucket, candidates)
        from kinbot.pes_phase import process_identity
        path = _record(root, token)
        with _locked(path.with_suffix('.lock')):
            entry = _read(path)
            entry['owner'] = process_identity(os.getpid())
            atomic_json(path, entry)
        atomic_json(_marker(job, 'producer'), dict(root=str(root), token=token))
    # Never hold the reuse lock while waiting for resource admission.
    result = launch()
    path = _record(root, token)
    with _locked(path.with_suffix('.lock')):
        entry = _read(path)
        if qc.job_ids.get(job) and entry['state'] not in ('complete', 'failed', 'unshareable'):
            entry.update(state='submitted', job_id=str(qc.job_ids[job]))
            atomic_json(path, entry)
    return result


def _source_payload(entry, root):
    """Read a completed generation without consuming its worker handoff."""
    source = _inside(root, entry['source'])
    fields = Path(str(source)+'.exitcode').read_text().split()
    if fields != [str(entry.get('job_id')), '0']:
        return None
    log = Path(str(source)+'_sella.log').read_text()
    if not log.rstrip().endswith('done'):
        return None
    # This registry only enrolls new jobs with no pre-existing database rows.
    # Thus a row cannot be mistaken for an older generation of the job name.
    pkl = Path(str(source)+'.pkl')
    try:
        raw = pkl.read_bytes()
    except FileNotFoundError:
        from ase.db.sqlite import SQLite3Database
        database = _inside(root, entry['directory'])/'kinbot.db'
        with sqlite3.connect(database.as_uri()+'?mode=ro', uri=True) as connection:
            db = SQLite3Database(str(database), use_lock_file=False)
            db.connection = connection
            db.change_count = 0
            rows = list(db.select(name=entry['job']))
        if not rows:
            return None
        row = rows[-1]
        return row.toatoms(), dict(row.data), row.id
    else:
        from ase import Atoms
        payload = pickle.loads(raw)
        if payload['name'] != entry['job']:
            raise ValueError('Shared QC handoff names a different job.')
        return Atoms(payload['sym'], positions=payload['pos']), payload['data'], None


def _finish(root, token, *, failed=False):
    """Freeze native evidence before another controller can reuse a result."""
    path = _record(root, token)
    with _locked(path.with_suffix('.lock')):
        entry = _read(path)
        if entry['state'] in ('complete', 'failed', 'unshareable'):
            return entry
        if failed:
            entry['state'] = 'failed'
            atomic_json(path, entry)
            return entry
        source = _inside(root, entry['source'])
        # A changed allocation must be adopted by its owning driver first.
        # Its ordinary timeout/manual recovery remains authoritative.
        retry_token = hashlib.sha256(entry['job'].encode()).hexdigest()[:24]
        retry_record = _inside(root, entry['directory'])/'.kinbot_timeout_retries'/retry_token/'current.json'
        if retry_record.exists():
            entry.update(state='unshareable', reason='Producer entered an independent retry.')
            atomic_json(path, entry)
            return entry
        try:
            payload = _source_payload(entry, root)
        except (FileNotFoundError, sqlite3.OperationalError, EOFError, pickle.UnpicklingError):
            return entry
        if payload is None:
            return entry
        # A staged worker can be rewritten while its controller resumes.
        # Do not release a pending alias merely because that file is transient.
        current = definition(SimpleNamespace(qc='gauss', use_sella=True),
                             entry['job'], path=Path(str(source)+'.py'))
        if current is None:
            return entry
        if not _program_matches(entry['definition'], current):
            entry.update(state='unshareable', reason='Producer calculation definition changed.')
            atomic_json(path, entry)
            return entry
        atoms, data, row_id = payload
        if (data.get('status') != 'normal' or atoms.get_chemical_symbols() != entry['definition']['symbols']
                or not np.isfinite(atoms.positions).all()
                or data.get('energy') is None or not np.isfinite(data['energy'])):
            return entry
        try:
            native = Path(str(source)+'.log').read_text(errors='replace')
        except FileNotFoundError:
            return entry
        if 'Normal termination of Gaussian' not in native[-6000:]:
            return entry
        if entry['job'].startswith('conf/') and len(atoms)>1:
            if (data.get('zpe') is None or not np.isfinite(data['zpe'])
                    or data.get('frequencies') is None or not len(data['frequencies'])
                    or not np.isfinite(data['frequencies']).all()):
                return entry
        # Original code and atoms stay in the bundle. No synthetic native
        # output is produced and no Hessian/vector is rotated independently.
        from ase.io.jsonio import encode
        bundles = root/NAME/'bundles'
        bundles.mkdir(parents=True, exist_ok=True)
        temporary = Path(tempfile.mkdtemp(prefix='.building-', dir=bundles))
        hashes = {}
        try:
            for suffix in (*SUFFIXES, '.py'):
                original = Path(str(source)+suffix)
                if not original.exists():
                    continue
                before = file_sha256(original)
                shutil.copyfile(original, temporary/('native'+suffix))
                if file_sha256(original) != before or file_sha256(temporary/('native'+suffix)) != before:
                    raise RuntimeError('QC reuse source changed during capture: '+str(original))
                hashes[suffix] = before
            # Check the whole capture again, including the generation fence.
            # A concurrent explicit retirement cannot mix attempts in a bundle.
            for suffix, sha in hashes.items():
                if file_sha256(Path(str(source)+suffix)) != sha:
                    raise RuntimeError('Shared source was revised during capture.')
            if Path(str(source)+'.exitcode').read_text().split() != [entry['job_id'], '0']:
                raise RuntimeError('Shared source generation changed during capture.')
            atomic_json(temporary/'record.json', dict(atoms=json.loads(encode(atoms.todict())),
                data=json.loads(encode(dict(data))), source_row_id=row_id,
                source=entry['source'], job_id=entry['job_id']))
            hashes['record'] = file_sha256(temporary/'record.json')
            atomic_json(temporary/'manifest.json', hashes)
            bundle = bundles/uuid.uuid4().hex
            temporary.rename(bundle)
        except BaseException:
            shutil.rmtree(temporary)
            raise
        entry.update(state='complete', bundle=str(bundle.relative_to(root)), hashes=hashes)
        atomic_json(path, entry)
        return entry


def observe(qc, job, status):
    marker = _marker(job, 'producer')
    if not marker.exists():
        return
    saved = _read(marker)
    root = Path(saved['root']).resolve()
    path = _record(root, saved['token'])
    with _locked(path.with_suffix('.lock')):
        entry = _read(path)
        if entry['state'] in ('complete', 'failed', 'unshareable'):
            return
        ident = getattr(qc, 'job_ids', {}).get(job)
        if ident and entry.get('job_id') not in (None, str(ident)):
            entry.update(state='unshareable', reason='Producer allocation replaced; original work retained.')
            atomic_json(path, entry)
            return
        if ident:
            entry.update(job_id=str(ident), state='submitted')
            atomic_json(path, entry)
    if status in ('normal', 'error'):
        _finish(root, saved['token'], failed=status == 'error')


def _adopt(qc, job, receipt, entry, root):
    from ase import Atoms
    from ase.io.jsonio import decode
    bundle = _inside(root, entry['bundle'])
    hashes = entry['hashes']
    if file_sha256(bundle/'record.json') != hashes['record']:
        raise RuntimeError('Shared QC result record changed.')
    payload = _read(bundle/'record.json')
    rows = list(qc.db.select(name=job))
    if rows and rows[-1].data.get('shared_qc_token') != entry['token']:
        raise RuntimeError(f'{job}: existing result remains authoritative; refusing shared overwrite.')
    marker = _marker(job, 'consumer')
    previous = receipt['state']
    receipt['state'] = 'copying'
    atomic_json(marker, receipt)
    for suffix, sha in hashes.items():
        if suffix == 'record':
            continue
        source, target = bundle/('native'+suffix), Path(job+suffix)
        if file_sha256(source) != sha:
            raise RuntimeError('Shared native artifact changed: '+str(source))
        if suffix == '.py':
            continue  # Keep the consumer's original request script.
        if target.exists():
            if previous not in ('copying', 'reused') or file_sha256(target) != sha:
                raise RuntimeError(f'{target}: existing artifact remains authoritative.')
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_name(target.name+'.reuse-'+uuid.uuid4().hex)
        shutil.copyfile(source, temporary)
        if file_sha256(temporary) != sha:
            raise RuntimeError('Shared QC copy verification failed.')
        os.replace(temporary, target)
    if not rows:
        atoms = Atoms.fromdict(decode(json.dumps(payload['atoms'])))
        data = decode(json.dumps(payload['data']))
        data.update(shared_qc_token=entry['token'], shared_qc_source=str(_inside(root, entry['source'])),
                    shared_qc_job_id=entry['job_id'], shared_qc_source_row_id=payload['source_row_id'])
        qc.db.write(atoms, name=job, data=data)
    receipt.update(state='reused', source=entry['source'], source_job_id=entry['job_id'],
                   bundle=entry['bundle'])
    atomic_json(marker, receipt)
    logger.info('%s: adopted complete native QC from %s (Slurm %s).', job, entry['source'], entry['job_id'])
    return 'normal'


def poll(qc, job):
    marker = _marker(job, 'consumer')
    if not marker.exists() or job in getattr(qc, 'job_ids', {}):
        return None
    receipt = _read(marker)
    if receipt['state'] == 'reused':
        return None  # Usual local native/database validation remains active.
    if not _same_program(receipt['definition'], job):
        raise RuntimeError(f'{job}: shared request definition changed; preserve its receipt and '
                           'revise through the owning driver before reusing QC.')
    root = Path(receipt['root']).resolve()
    path = _record(root, receipt['token'])
    entry = _read(path)
    if entry.get('job_id') is None:
        # Submission may have finished before the producer recorded its ID.
        # Read the durable shared admission ledger; never infer an absent job
        # merely from a dead controller or a temporarily empty queue listing.
        admissions = root/'.kinbot_policy_admissions.json'
        claims = list(_read(admissions).values()) if admissions.exists() else []
        for p in (root/'.kinbot_policy_submissions').glob('*.json'):
            claims.append(_read(p))
        matched = [c for c in claims if c.get('directory') == str(_inside(root, entry['directory']))
                   and c.get('command') == ['sbatch', entry['job']+'.sbatch']
                   and c.get('created_at', '') >= entry['created_at']]
        ids = {str(c['job_id']) for c in matched if c.get('job_id')}
        if len(ids) == 1:
            with _locked(path.with_suffix('.lock')):
                entry = _read(path)
                entry.update(job_id=ids.pop(), state='submitted')
                atomic_json(path, entry)
        # An uncertain admission or dead owner is not permission to submit
        # another copy. A normal controller refresh resumes the producer.
        # Keep the consumer pending across the staging/admission window.
    entry = _finish(root, receipt['token'])
    if entry['state'] == 'complete':
        return _adopt(qc, job, receipt, entry, root)
    if entry['state'] == 'unshareable':
        # Cache ineligibility is not scientific rejection. Preserve the
        # association and let the consumer's own normal driver proceed.
        marker.rename(marker.with_name(marker.name+'.released-'+uuid.uuid4().hex))
        return None
    if entry['state'] == 'failed':
        receipt.update(state='source_failed', source=entry['source'], source_job_id=entry.get('job_id'))
        atomic_json(marker, receipt)
        logger.warning('%s: equivalent shared attempt %s was unsuccessful; no duplicate retry.', job, entry['source'])
        return 'error'
    return 'running'
