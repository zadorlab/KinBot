"""Share equivalent new Gaussian/Sella conformer, minimum L2 and HIR attempts.

This is calculation reuse, not species/TS equivalence or population merging.
An immutable native bundle is copied to each consumer before its DB row appears.
Existing local attempts are never replaced and Slurm ownership is never aliased.
"""
import ast
from itertools import islice
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


def definition(qc, job, *, path=None, minimum_source=False, legacy=False):
    """Parse literals without executing a worker; unsupported forms opt out."""
    if (getattr(qc, 'qc', None) != 'gauss' or not getattr(qc, 'use_sella', False)
            or not Path(job).parts or (not minimum_source
                and Path(job).parts[0] not in ('hir', 'conf')
                and not job.endswith('_well_high'))):
        return None
    try:
        raw = Path(path or job+'.py').read_text()
        tree = ast.parse(raw)
        atoms = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                 and isinstance(n.func, ast.Name) and n.func.id == 'Atoms']
        if len(atoms) != 1:
            return None
        values = {k.arg: ast.literal_eval(k.value) for k in atoms[0].keywords}
        if set(values) != {'symbols', 'positions'}:
            return None  # External masses/charges need an explicit contract.
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
        constraints = []
        kind = 'minimum'
        order = [ast.literal_eval(n.value) for n in ast.walk(tree)
                 if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name)
                 and t.id == 'order' for t in n.targets)]
        for n in ast.walk(tree):
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == 'Sella':
                for k in n.keywords:
                    if k.arg == 'order' and isinstance(k.value, ast.Constant):
                        order.append(ast.literal_eval(k.value))
            if (isinstance(n, ast.Assign) and any(isinstance(t, ast.Name)
                    and t.id == 'base_0_fix' for t in n.targets)):
                if not isinstance(n.value, ast.ListComp) or len(n.value.generators) != 1:
                    return None
                indices = ast.literal_eval(n.value.generators[0].iter)
                if len(indices) != 4 or len(set(indices)) != 4 or any(
                        type(i) is not int or i < 1 or i > len(symbols) for i in indices):
                    return None
                constraints = [i-1 for i in indices]
                if not legacy:
                    n.value.generators[0].iter = ast.Constant('DIHEDRAL_ATOMS')
                kind = 'hir'
        if not order or any(o != 0 for o in order):
            kind = 'transition_state'
        if minimum_source and kind != 'minimum':
            return None
        if job.endswith('_well_high') and kind != 'minimum':
            return None
        reduced = {k: v for k, v in sorted(kwargs.items())
                   if k not in ('label', 'chk', 'mem', 'nprocshared')}
        assigns[0].value = ast.parse(repr(reduced), mode='eval').body
        for k in atoms[0].keywords:
            if k.arg == 'positions':
                k.value = ast.Constant('GEOMETRY')
            elif k.arg == 'symbols' and not legacy:
                k.value = ast.Constant('ORDERED_ATOMS')
        class Labels(ast.NodeTransformer):
            def visit_List(self, node):
                # Native frequency readers also receive the ordered symbols.
                # Normalize every literal copy of that same atom list.
                if not legacy:
                    try:
                        if ast.literal_eval(node) == symbols:
                            return ast.Constant('ORDERED_ATOMS')
                    except (ValueError, TypeError):
                        pass
                return self.generic_visit(node)

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
        if not legacy:
            contract['atoms'] = sorted(symbols)
        if template:
            contract['queue_template'] = file_sha256(template)
        result = dict(contract_sha256=digest(contract), symbols=symbols,
                    positions=positions.tolist(), worker_sha256=hashlib.sha256(raw.encode()).hexdigest(),
                    program_sha256=digest(contract['program']), kind=kind,
                    constraint_atoms=constraints, charge=kwargs.get('charge'),
                    multiplicity=kwargs.get('mult'))
        if legacy:
            result.pop('kind')
            result.pop('constraint_atoms')
        proof = Path(job+'.minimum_reuse.json')
        if path is None and proof.is_file():
            saved = _read(proof)
            if saved.get('worker_sha256') == result['worker_sha256']:
                result['minimum'] = saved['minimum']
        return result
    except (OSError, SyntaxError, ValueError, TypeError, KeyError):
        return None


def matching(a, b):
    """Return b indices in a order; proper rotations only, constraints mapped.

    Exact requests can share without a minimum proof. Relaxing coordinate
    serialization tolerance requires two independently verified minima AND
    the same constrained geometry under the very same atom permutation.
    """
    if a['contract_sha256'] != b['contract_sha256'] or sorted(a['symbols']) != sorted(b['symbols']):
        return None
    from kinbot.minimum_reuse import comparison, LIMITS
    tolerance = (LIMITS['max_distance_angstrom'] if a.get('minimum') and b.get('minimum')
                 else GEOMETRY_TOLERANCE)
    x, y = np.asarray(a['positions']), np.asarray(b['positions'])
    def check(p):
        if a['symbols'] != [b['symbols'][i] for i in p]:
            return None
        ca, cb = a.get('constraint_atoms', []), b.get('constraint_atoms', [])
        if [p[i] for i in ca] not in (cb, cb[::-1]):
            return None
        error = float(max(proper_alignment(x, y[p]).errors))
        proof = comparison(a.get('minimum'), b.get('minimum'), p)
        if error > GEOMETRY_TOLERANCE and (error > tolerance or proof is None):
            return None
        return dict(permutation=p, max_distance_angstrom=error, minimum_comparison=proof)
    identity = check(list(range(len(x))))
    if identity:
        return identity
    # TS modes/path definitions need their own permutation contract. Existing
    # same-order sharing remains available; this extension is for minima.
    if a.get('kind') == 'transition_state' or b.get('kind') == 'transition_state':
        return None
    import networkx as nx
    def graph(symbols, positions):
        g = nx.Graph()
        g.add_nodes_from((i, dict(element=s)) for i, s in enumerate(symbols))
        g.add_edges_from((i,j,dict(distance=float(np.linalg.norm(positions[i]-positions[j]))))
                         for i in range(len(positions)) for j in range(i))
        return g
    matcher = nx.algorithms.isomorphism.GraphMatcher(graph(a['symbols'], x), graph(b['symbols'], y),
        node_match=lambda u,v: u['element']==v['element'],
        edge_match=lambda u,v: abs(u['distance']-v['distance']) <= 2*tolerance+1.e-12)
    for m in islice(matcher.isomorphisms_iter(), 10000):
        result = check([m[i] for i in range(len(x))])
        if result:
            return result
    return None  # A bounded mapping search cannot authorize a false match.


def equivalent(a, b):
    return matching(a, b) is not None


def _same_program(saved, job, path=None):
    current = definition(SimpleNamespace(qc="gauss", use_sella=True), job, path=path,
                         legacy='kind' not in saved)
    return _program_matches(saved, current)


def _program_matches(saved, current):
    return (current is not None and current["program_sha256"] == saved["program_sha256"]
            and current["symbols"] == saved["symbols"]
            and current.get('constraint_atoms', []) == saved.get('constraint_atoms', [])
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
            or not Path(job).parts or (Path(job).parts[0] not in ('hir', 'conf')
                                      and not job.endswith('_well_high'))):
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
            match = matching(wanted, entry['definition'])
            if match is not None:
                atomic_json(consumer, dict(schema=1, root=str(root), token=ident,
                    definition=wanted, mapping=match, state='waiting'))
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
                             entry['job'], path=Path(str(source)+'.py'),
                             legacy='kind' not in entry['definition'])
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
        if (entry['definition'].get('kind') == 'minimum' or entry['job'].startswith('conf/')) and len(atoms)>1:
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
    verified = matching(receipt['definition'], entry['definition'])
    if (verified is None or receipt.get('mapping', {}).get('permutation', verified['permutation'])
            != verified['permutation']):
        raise RuntimeError('Shared QC atom mapping no longer matches its captured request evidence.')
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
        permutation = verified['permutation']
        if permutation != list(range(len(atoms))):
            # Use the producer's Cartesian frame, with the consumer's atom
            # order. Geometry, masses and every Hessian block move together.
            # Native files remain byte-for-byte original; readers below use
            # the native frame descriptor instead of pretending otherwise.
            atoms = atoms[permutation]
            if data.get('hess') is not None:
                data['hess'] = permute_hessian(data['hess'], permutation)
            data['shared_qc_native'] = dict(root=str(root), bundle=entry['bundle'],
                hashes=entry['hashes'], permutation=permutation, job=entry['job'])
        data.update(shared_qc_token=entry['token'], shared_qc_source=str(_inside(root, entry['source'])),
                    shared_qc_job_id=entry['job_id'], shared_qc_source_row_id=payload['source_row_id'])
        qc.db.write(atoms, name=job, data=data)
    receipt.update(state='reused', source=entry['source'], source_job_id=entry['job_id'],
                   bundle=entry['bundle'])
    atomic_json(marker, receipt)
    logger.info('%s: adopted complete native QC from %s (Slurm %s).', job, entry['source'], entry['job_id'])
    return 'normal'


def permute_hessian(hessian, permutation):
    """H_target = P H_native P^T, with each atomic Cartesian block intact."""
    n = len(permutation)
    h = np.asarray(hessian, float)
    if h.shape != (3*n, 3*n) or not np.isfinite(h).all():
        raise ValueError('Shared minimum has an invalid Cartesian Hessian.')
    indices = (3*np.asarray(permutation)[:, None]+np.arange(3)).ravel()
    return h[np.ix_(indices, indices)]


def native_record(row):
    """Original atoms/frame for a permuted reused row, including published aliases."""
    saved = row.data.get('shared_qc_native')
    if not saved:
        return None
    from ase import Atoms
    from ase.io.jsonio import decode
    root = Path(saved['root']).resolve()
    bundle = _inside(root, saved['bundle'])
    if file_sha256(bundle/'record.json') != saved['hashes']['record']:
        raise ValueError('Shared native frame record changed.')
    record = _read(bundle/'record.json')
    atoms = Atoms.fromdict(decode(json.dumps(record['atoms'])))
    p = list(saved['permutation'])
    if (sorted(p) != list(range(len(atoms))) or list(atoms[p].symbols) != list(row.symbols)
            or not np.array_equal(atoms[p].positions, row.positions)):
        raise ValueError('Shared native frame does not match the selected result.')
    return bundle, dict(symbols=list(atoms.symbols), positions=atoms.positions,
                        job=saved['job']), p


def mapped_hessian(qc, job, row):
    """Read native force constants in their actual order, then permute once."""
    shared = native_record(row)
    if shared is None:
        return None
    bundle, native, p = shared
    from kinbot.gaussian_hessian import read_formatted_hessian
    import subprocess
    fchk = bundle/'native.fchk'
    if not fchk.exists():
        # Conversion writes only a consumer-local derivative of the immutable
        # original checkpoint. It never edits the bundle or launches QC.
        chk = bundle/'native.chk'
        if not chk.exists():
            return []
        if file_sha256(chk) != row.data['shared_qc_native']['hashes']['.chk']:
            raise ValueError('Shared native checkpoint changed.')
        fchk = Path(job+'.shared_native.fchk')
        if not fchk.exists():
            temp = fchk.with_name(fchk.name+'.'+uuid.uuid4().hex)
            subprocess.run(['formchk', str(chk), str(temp)], check=True, stdout=subprocess.DEVNULL)
            os.replace(temp, fchk)
    else:
        if file_sha256(fchk) != row.data['shared_qc_native']['hashes']['.fchk']:
            raise ValueError('Shared formatted checkpoint changed.')
    return permute_hessian(read_formatted_hessian(fchk, native['symbols'], native['positions']), p)


def mapped_modes(row, *, geometry=None, symbols=None):
    """Map native displacements with their geometry; no independent mode fit."""
    shared = native_record(row)
    if shared is None:
        return None
    bundle, native, p = shared
    path = bundle/'native.log'
    if file_sha256(path) != row.data['shared_qc_native']['hashes']['.log']:
        raise ValueError('Shared native frequency output changed.')
    from kinbot.gaussian_modes import read_normal_modes
    result = read_normal_modes(path, native['symbols'], native['positions'])
    modes = np.asarray(result['modes'])[:, p, :]
    if geometry is not None:
        if list(map(str, symbols)) != list(row.symbols):
            raise ValueError('Selected mode atom order differs from the reused result.')
        fit = proper_alignment(row.positions, geometry)
        if max(fit.errors) > 1.e-5:
            raise ValueError('Selected mode geometry differs from the reused result.')
        modes = modes @ fit.rotation
    result['modes'] = modes.tolist()
    return result


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
