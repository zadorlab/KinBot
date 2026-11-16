"""Finite repairs of failed native diagnostic probes, never PES acceptance.

Return new tasks/graphs with hashed native provenance. The owner applies them;
failed attempts remain intact. No repair may recursively repair its own trial.
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
