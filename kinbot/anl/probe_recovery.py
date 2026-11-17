"""Finite repairs of failed native diagnostic probes, never PES acceptance.

Return new tasks/graphs with hashed native provenance. The owner applies them;
failed attempts remain intact. Each recipe has one trial and explicit ancestry.
"""
from copy import deepcopy
import json
from pathlib import Path
import re
from statistics import median

from kinbot.anl import checkpoint_inputs, cc_iteration_recovery, path_review
from kinbot.anl.tasks import molpro_task


def dvv_step_repair(graph, ident):
    """Repair only a negative-time DVV rejection before any propagated point."""
    graph = Path(graph).resolve()
    _, task, outcome, proof = path_review.source(graph, ident, [ident+'.log'])
    metadata = task.get('native_path_probe', {})
    if (task.get('backend') != 'gaussian' or metadata.get('recipe') != 'dvv'
            or task.get('native_probe_repair') or outcome.get('status') != 'failed'
            or outcome.get('failure_phase') != 'execution'):
        raise ValueError('DVV repair requires a failed, unrepaired native DVV probe.')
    text = task['input_template']
    route = re.search(r'IRC=\(([^()]*)\)', text, re.I)
    if not route or '--Link1--' in text: raise ValueError('Expected one raw DVV path.')
    options = route[1].lower().split(',')
    required = {'rcfc','gradientonly','dvv','massweighted','verytight','report=cartesians','maxcycle=100'}
    steps = re.findall(r'\bstepsize=(-[1-9]\d*)\b', route[1], re.I)
    if (len(options) != 10 or not required <= set(options) or len(steps) != 1
            or ('forward' in options) == ('reverse' in options)
            or not any(re.fullmatch(r'maxpoints=[1-9]\d*', s) for s in options)):
        raise ValueError('DVV input is not the negative-time failure class.')
    directory = graph/'tasks'/ident
    native = (directory/(ident+'.log')).read_text(errors='replace')
    reported = re.search(r'RxnCyc: Step size is\s+('+path_review.NUMBER+')', native)
    if (not reported or abs(path_review._number(reported[1])-int(steps[0])/100) > 1.e-7
            or 'RxnCyc: Invalid step size found!' not in native
            or not re.search(r'Integration scheme\s*=\s*DVV', native)
            or 'Error termination' not in native
            or re.search(r'Pt\s+\d+\s+Step|Point Number:|SCF Done:|NET REACTION COORDINATE', native)):
        raise ValueError('Native output does not prove a pre-propagation negative-time rejection.')
    bindings = checkpoint_inputs.inspect_dependencies(directory, task['input_name'],
                                                       (directory/task['input_name']).read_text())
    if len(bindings) != 1: raise ValueError('Expected one accepted dependency Hessian.')
    for binding in bindings:
        parent = Path(binding['source'])
        proof['files'][str(parent.relative_to(graph))] = binding['sha256']
        proof['files'][str((parent.parent/'execution.json').relative_to(graph))] = binding['source_execution_sha256']
    new_id = ident+'_positive_dt_v2'
    replacement = re.sub(r'\bStepSize=-([1-9]\d*)\b', r'StepSize=\1', text, count=1, flags=re.I)
    replacement = re.sub(r'(?im)^%chk=.*$', '%chk='+new_id+'.chk', replacement)
    new = deepcopy(task)
    new.update(id=new_id, input_name=new_id+'.com', stdin=new_id+'.com',
        stdout=new_id+'.log', stderr=new_id+'.err', input_template=replacement,
        required_outputs=[new_id+'.log', new_id+'.chk'],
        success_marker=dict(file=new_id+'.log', contains='Normal termination of Gaussian'),
        failure_markers=[dict(file=new_id+'.log', contains='Error termination')],
        native_probe_repair=dict(recipe='dvv_positive_time', source=ident,
            evidence=path_review.digest(proof), budget=1, scientific_approval=False,
            time_step_fs=abs(int(steps[0]))/100, spatial_path_length='not specified'))
    path_review.verify(proof)
    return dict(task=new, proof=proof, report=dict(failure='invalid_dvv_time_step'))


def dvv_departure_probe(graph, ident):
    """Bounded raw path after DVV declares a minimum at its first TS step.

    Gaussian Overlay 1 documents IOp(1/108=-1) as disabling only the DVV
    gradient-magnitude stopping test. Keep the force/velocity-angle test and
    integration controls; impose 64 points, with no endpoint optimization.
    This is a departure diagnostic, never a minimum or connection acceptance.
    """
    graph = Path(graph).resolve()
    _, task, outcome, proof = path_review.source(graph, ident, [ident+'.log'])
    text = task['input_template']
    route = re.search(r'IRC=\(([^()]*)\)', text, re.I)
    metadata = task.get('native_probe_repair', {})
    if (task.get('backend') != 'gaussian' or metadata.get('recipe') != 'dvv_positive_time'
            or outcome.get('status') != 'executed' or not route
            or '--Link1--' in text or re.search(r'\biop\s*\(', text, re.I)
            or any(k in task for k in ('result_parser','geometry_output','manual_fallback'))):
        raise ValueError('Departure probe requires an executed positive-time raw DVV repair.')
    options = route[1].lower().split(',')
    required = {'rcfc','gradientonly','dvv','massweighted','verytight','report=cartesians','maxcycle=100'}
    steps = re.findall(r'\bstepsize=([1-9]\d*)\b', route[1], re.I)
    bounds = re.findall(r'\bmaxpoints=([1-9]\d*)\b', route[1], re.I)
    if (len(options) != 10 or not required <= set(options) or len(steps) != 1
            or len(bounds) != 1 or int(bounds[0]) < 64
            or ('forward' in options) == ('reverse' in options)):
        raise ValueError('Unrecognized or exhausted DVV departure settings.')
    directory = graph/'tasks'/ident
    native = (directory/(ident+'.log')).read_text(errors='replace')
    report = path_review.gaussian_path(native.splitlines())
    dt = re.search(r'Step size\s*=\s*('+path_review.NUMBER+r')\s+fs', native)
    if (report['classification'] != 'early_gradient_stop' or report['points'] != 1
            or not 0 < report.get('coordinate', 0) <= .01
            or not 0 < report.get('gradient', 0) < .0015
            or not re.search(r'Integration scheme\s*=\s*DVV', native)
            or 'Current Structure is TS' not in native or not dt
            or abs(path_review._number(dt[1])-int(steps[0])/100) > 1.e-7):
        raise ValueError('Native evidence does not establish a first-step DVV gradient stop near the TS.')
    # Bind this successful input to the previously reviewed sign-only repair.
    ancestry = dvv_step_repair(graph, metadata['source'])
    saved = deepcopy(task); saved['resources'] = ancestry['task']['resources']
    if saved != ancestry['task']:
        raise ValueError('Positive-time source differs from its recorded sign repair.')
    bindings = checkpoint_inputs.inspect_dependencies(directory, task['input_name'],
                                                       (directory/task['input_name']).read_text())
    if len(bindings) != 1: raise ValueError('Departure probe requires one verified parent Hessian.')
    for binding in bindings:
        parent = Path(binding['source'])
        proof['files'][str(parent.relative_to(graph))] = binding['sha256']
        proof['files'][str((parent.parent/'execution.json').relative_to(graph))] = binding['source_execution_sha256']
    new_id = ident+'_departure_v3'
    replacement = re.sub(r'\bMaxPoints=[1-9]\d*\b', 'MaxPoints=64', text, count=1, flags=re.I)
    replacement = re.sub(r'(?im)^%chk=.*$', '%chk='+new_id+'.chk', replacement)
    replacement = replacement.rstrip()+' IOp(1/108=-1)\n\n'
    new = deepcopy(task)
    new.update(id=new_id, input_name=new_id+'.com', stdin=new_id+'.com',
        stdout=new_id+'.log', stderr=new_id+'.err', input_template=replacement,
        required_outputs=[new_id+'.log',new_id+'.chk'],
        success_marker=dict(file=new_id+'.log', contains='Normal termination of Gaussian'),
        failure_markers=[dict(file=new_id+'.log', contains='Error termination')],
        native_probe_repair=dict(recipe='dvv_departure', source=ident, budget=1,
            evidence=path_review.digest(proof), ancestor_evidence=path_review.digest(ancestry['proof']),
            scientific_approval=False, max_points=64, endpoint_acceptance=False,
            gradient_stop='disabled for bounded departure diagnostic',
            documentation='https://gaussian.com/overlay1/#iop_(1/108)'))
    path_review.verify(proof); path_review.verify(ancestry['proof'])
    return dict(task=new, proof=proof, report=report)


def electronic_failure(text):
    """Recognize two distinct failures before the first derivative starts.

    The CC test examines the final iteration table, not the preceding MP2
    table. Printed zero energy changes are not evidence of convergence.
    """
    if re.search(r'out.of.memory|no space left|time limit|insufficient memory', text, re.I):
        raise ValueError('Infrastructure failure needs separate recovery.')
    if (re.search(r'PROGRAM \* (?:FORCES|FREQUENCIES)|Numerical (?:gradient|hessian)', text, re.I)
            or len(re.findall(r'PROGRAM\s*\*\s*CCSD', text, re.I)) > 1):
        raise ValueError('A displaced-geometry failure needs separate review.')
    if 'No convergence in rhfpr' in text and 'RHF not converged' in text:
        rows = []
        for line in text.splitlines():
            f = line.split()
            if len(f) >= 10 and f[0].isdigit() and f[-1].startswith('diag'):
                try: rows.append((int(f[0]), path_review._number(f[3])))
                except ValueError: pass
        if (len(rows) < 20 or rows[-1][0] != 60 or min(x[1] for x in rows[-10:]) <= 1.e-6
                or median(x[1] for x in rows[-5:]) < .75*median(x[1] for x in rows[-10:-5])):
            raise ValueError('Insufficient native evidence for stagnant default ROHF.')
        return dict(recipe='rohf_level_shift', last_iteration=60,
                    density_tail=[x[1] for x in rows[-10:]])
    if ('CONVERGENCE NOT REACHED AFTER MAX. ITERATIONS' not in text
            or 'The problem occurs in ccsd' not in text):
        raise ValueError('No recognized native electronic failure.')
    rows = []
    for line in text.splitlines():
        if re.search(r'ITER\.\s+SQ.NORM\s+CORR.ENERGY', line): rows = []
        f = line.split()
        if len(f) == 11 and f[0].isdigit():
            try: rows.append((int(f[0]), path_review._number(f[6])+path_review._number(f[7])))
            except ValueError: pass
    if (len(rows) < 20 or rows[-1][0] != 39
            or [r[0] for r in rows[-10:]] != list(range(30,40))
            or not all(0 < r[1] < 1.e-9 for r in rows[-10:])
            or rows[-1][1] >= .5*rows[-10][1]
            or median(r[1] for r in rows[-5:]) >= .7*median(r[1] for r in rows[-10:-5])):
        raise ValueError('No sustained small-amplitude contraction supports a bounded CC extension.')
    return dict(recipe='contracting_cc100', last_iteration=39,
                variance_tail=[r[1] for r in rows[-10:]])


DERIVATIVES = {
    'gradient': 'forces,numerical,central,dstep=0.005,rstep=0.005,startcmd=rhf',
    'hessian_010': 'frequencies,numerical,step=0.01,new,print=0',
    'hessian_005': 'frequencies,numerical,step=0.005,new,print=0',
}


def curvature_repair(graph):
    """A single electronic gate then dependent fixed-geometry derivatives.

    No threshold relaxation, ignored errors, orbital import, optimization,
    geometry export or result parser. All three failed tasks must independently
    support the same repair; later failures remain for review.
    """
    from kinbot.anl import dispatch
    from ase.io import read
    import numpy as np
    graph, spec, _ = dispatch._load(graph)
    if (spec.get('intent', {}).get('integration') != 'none'
            or spec.get('intent', {}).get('scientific_approval') is not False
            or spec.get('intent', {}).get('native_probe_repair')
            or {t['id'] for t in spec['tasks']} != set(DERIVATIVES)):
        raise ValueError('Expected the original unaccepted three-task curvature audit.')
    proofs, reports, templates = [], [], []
    for task in spec['tasks']:
        ident = task['id']; directory = graph/'tasks'/ident
        names = [ident+'.out'] + ([ident+'.log'] if (directory/(ident+'.log')).exists() else [])
        _, saved, outcome, proof = path_review.source(graph, ident, names)
        record = json.loads((directory/'task.json').read_text())
        atoms = read(directory/'geometry.xyz')
        molecule = spec['molecule']
        if (record['molecule'] != dict(symbols=molecule['symbols'],
                charge=molecule.get('charge',0), multiplicity=molecule.get('multiplicity',1))
                or atoms.get_chemical_symbols() != molecule['symbols']
                or not np.allclose(atoms.positions, molecule['positions'], rtol=0, atol=5.1e-9)):
            raise ValueError('Diagnostic geometry or electronic state differs from the native attempt.')
        if outcome.get('status') != 'failed' or outcome.get('failure_phase') != 'execution':
            raise ValueError('Only failed electronic audit attempts can seed this repair.')
        text = task['input_template']
        basis = re.findall(r'(?m)^basis=([\w()-]+)$', text)
        if len(basis) != 1: raise ValueError('Ambiguous diagnostic basis.')
        body = f'basis={basis[0]}\nrhf\nuccsd(t)-f12b,scale_trip=1\n'
        expected = molpro_task(ident, body+DERIVATIVES[ident]+'\n', geometry_from='initial',
                               scf_convergence=1.e-12, cc_convergence=1.e-11)
        if (text != expected['input_template'] or task.get('geometry_from') != 'initial'
                or task.get('kind') != 'external' or task.get('backend') != 'molpro'
                or any(k in task for k in ('result_parser','geometry_output','manual_fallback','native_probe_repair'))
                or any(task.get(k) for k in ('depends_on','stage_files','files_from_env','geometry_record','restart'))):
            raise ValueError('Changed electronic controls or exhausted curvature settings.')
        native = (directory/(ident+'.out')).read_text(errors='replace')
        report = electronic_failure(native)
        new_text = text
        if report['recipe'] == 'contracting_cc100':
            # Reuse the shared CC MAXIT edit, with a stricter applicability
            # test here. The parser hint is only for planning, never installed
            # on a diagnostic task (which must not get automatic retries).
            hint = dict(saved, result_parser=dict(method='UCCSD(T)-F12B', file=ident+'.out'))
            planned = cc_iteration_recovery.options(hint, directory, outcome)
            if not planned or planned[1]['limit'] != 100:
                raise ValueError('Shared CC procedure has no unexhausted bounded edit.')
            for edit in planned[0]['input_edits']:
                new_text = new_text.replace(edit['old'], edit['new'])
        else:
            # Molpro's documented first legacy remedy for oscillating ROHF.
            # Keep 60 iterations and the original ROHF determinant definition.
            new_text = new_text.replace('\nrhf\n', '\n{rhf,maxit=60;shift,-1.0,-0.5}\n')
        proofs.append(proof); reports.append(report); templates.append((task,new_text))
    if len({r['recipe'] for r in reports}) != 1:
        raise ValueError('Mixed electronic failures require independent review.')
    recipe = reports[0]['recipe']
    prefix = templates[0][1].rsplit(DERIVATIVES[templates[0][0]['id']], 1)[0]
    gate = molpro_task('electronic', '', geometry_from='initial')
    gate['input_template'] = prefix.replace('\nbasis=', '\ngprint,orbitals,orben\nbasis=') + 'put,molden,electronic.molden\n'
    gate['resources'] = deepcopy(templates[0][0]['resources'])
    gate['required_outputs'].append('electronic.molden')
    tasks = [gate]
    for task,text in templates:
        new = deepcopy(task)
        new.update(input_template=text, depends_on=['electronic'])
        tasks.append(new)
    for task in tasks:
        task['failure_markers'] = [dict(file=task['id']+'.out', contains=s) for s in
            ('No convergence','NO CONVERGENCE','CONVERGENCE NOT REACHED','RHF not converged')]
        task['native_probe_repair'] = dict(recipe=recipe, budget=1, scientific_approval=False,
                                           evidence=path_review.digest(proofs))
    result = deepcopy(spec)
    result.update(name=spec['name']+'_electronic_v2', tasks=tasks)
    result['intent'] = dict(purpose='Bounded electronic recovery before fixed-geometry curvature audit.',
        source_proofs=proofs, native_probe_repair=recipe, integration='none', scientific_approval=False)
    dispatch.validate_spec(result)
    for proof in proofs: path_review.verify(proof)
    return dict(spec=result, proofs=proofs, reports=reports)


def quadratic_rohf_probe(graph):
    """One electronic-only SO trial after default and shifted ROHF stagnate.

    Reuse the reconnaissance determinant renderer. No derivatives are released:
    converged orbitals/energies still need comparison with the target reference.
    Native version, occupations and both exhausted inputs are evidence-bound.
    """
    from ase.data import atomic_numbers
    from ase.io import read
    import numpy as np
    from kinbot.anl import dispatch
    from kinbot.molpro_solver import rhf_determinant

    graph, spec, state = dispatch._load(graph)
    intent = spec.get('intent', {})
    parents = intent.get('source_proofs', [])
    if (intent.get('integration') != 'none' or intent.get('scientific_approval') is not False
            or intent.get('native_probe_repair') != 'rohf_level_shift'
            or {t['id'] for t in spec['tasks']} != {'electronic', *DERIVATIVES}
            or len(parents) != 3 or len({p['graph'] for p in parents}) != 1):
        raise ValueError('Quadratic trial requires the exhausted default/shifted ROHF audit.')
    for proof in parents: path_review.verify(proof)
    original = curvature_repair(parents[0]['graph'])
    expected = original['spec']
    normalized = deepcopy(spec)
    resources = {t['id']: t['resources'] for t in expected['tasks']}
    for task in normalized['tasks']: task['resources'] = resources[task['id']]
    if normalized != expected:
        raise ValueError('Shifted audit differs from its recorded recovery procedure.')
    for ident in DERIVATIVES:
        entry = state['tasks'].get(ident, {})
        if (entry.get('job_id') or entry.get('status', 'waiting') not in ('waiting', 'staged', 'blocked')
                or (graph/'tasks'/ident/'execution.json').exists()):
            raise ValueError('A derivative already started; review its evidence separately.')
    _, task, outcome, proof = path_review.source(graph, 'electronic', ['electronic.out'])
    if outcome.get('status') != 'failed' or outcome.get('failure_phase') != 'execution':
        raise ValueError('Expected a failed electronic gate.')
    directory = graph/'tasks/electronic'
    native = (directory/'electronic.out').read_text(errors='replace')
    report = electronic_failure(native)
    if report['recipe'] != 'rohf_level_shift':
        raise ValueError('Shifted ROHF did not show the required stagnation.')
    version = re.search(r'\bVersion\s+(\d{4})\.(\d+)', native, re.I)
    if not version or int(version[1]) < 2021:
        raise ValueError('Native evidence must identify Molpro 2021 or later for RHF,SO.')
    molecule = spec['molecule']; atoms = read(directory/'geometry.xyz')
    record = json.loads((directory/'task.json').read_text())
    definition = dict(symbols=molecule['symbols'], charge=molecule.get('charge', 0),
                      multiplicity=molecule.get('multiplicity', 1))
    if (record['molecule'] != definition or atoms.get_chemical_symbols() != molecule['symbols']
            or not np.allclose(atoms.positions, molecule['positions'], rtol=0, atol=5.1e-9)):
        raise ValueError('Electronic trial geometry/state differs from its native source.')
    electrons = sum(atomic_numbers[s] for s in molecule['symbols']) - definition['charge']
    spin = definition['multiplicity'] - 1
    alpha = re.findall(r'Final alpha occupancy:\s*(\d+)\s*$', native, re.M)
    beta = re.findall(r'Final beta\s+occupancy:\s*(\d+)\s*$', native, re.M)
    if (len(alpha) != 1 or len(beta) != 1
            or int(alpha[0]) + int(beta[0]) != electrons
            or int(alpha[0]) - int(beta[0]) != spin):
        raise ValueError('Native ROHF occupations do not establish the requested determinant counts.')
    old = '{rhf,maxit=60;shift,-1.0,-0.5}\n'
    if task['input_template'].count(old) != 1:
        raise ValueError('No unique exhausted level-shift input.')
    new = deepcopy(task)
    new['input_template'] = task['input_template'].replace(old,
        rhf_determinant('rhf,so,maxit=60,energy=1.d-12', electrons, spin, '2100.2'))
    sources = [*original['proofs'], proof]
    new['native_probe_repair'] = dict(recipe='rohf_quadratic', budget=1,
        scientific_approval=False, evidence=path_review.digest(sources))
    result = deepcopy(spec)
    result.update(name=spec['name']+'_quadratic_v3', tasks=[new])
    result['intent'] = dict(purpose='Quadratic ROHF electronic reference audit after DIIS and shift stagnation.',
        source_proofs=sources, native_probe_repair='rohf_quadratic', integration='none',
        scientific_approval=False, derivative_release='requires_reference_review')
    dispatch.validate_spec(result)
    for item in sources: path_review.verify(item)
    return dict(spec=result, proofs=sources, reports=[dict(report, recipe='rohf_quadratic',
        molpro_version=version[0], alpha=int(alpha[0]), beta=int(beta[0]), max_iterations=60)])
