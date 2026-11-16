"""Evidence-bound, finite native path probes; never scientific acceptance.

Keep an IRC endpoint separate from a subsequent optimization. These helpers
propose additional tasks, preserving the source attempt and accepted Hessian.
The owning driver must serialize application with live_control.commit().
"""
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import re
from statistics import median

from kinbot.execution import file_sha256

NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[DEde][-+]?\d+)?'


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _number(text):
    return float(text.replace('D', 'E').replace('d', 'e'))


def gaussian_path(lines):
    """Read native point records, first IRC termination, and later geometry.

    Geometry connectivity below is only a review flag, never a bond assignment
    or a reason to accept/reject a reaction. Point counts alone prove no edge.
    """
    result = dict(points=0, stop='unknown', scientific_approval=False)
    section, rows, row, coordinates, mode, last_geometry = 0, [], None, [], None, None
    for line in lines:
        if 'CURRENT STRUCTURE' in line or 'Input orientation:' in line or 'Standard orientation:' in line:
            mode, coordinates = ('current' if 'CURRENT STRUCTURE' in line else 'orientation'), []
        fields = line.split()
        if mode:
            count = 5 if mode == 'current' else 6
            if len(fields) == count and fields[0].isdigit() and fields[1].isdigit():
                try: coordinates.append([int(fields[1]), *[_number(s) for s in fields[-3:]]])
                except ValueError: pass
            elif coordinates:
                last_geometry = coordinates
                if mode == 'current' and section == 0:
                    result['irc_geometry'] = coordinates
                mode = None
        if section == 0:
            match = re.search(r'Pt\s+(\d+)\s+Step number\s+(\d+)\s+out of a maximum of\s+(\d+)', line)
            if match:
                if row: rows.append(row)
                if rows and rows[-1]['point'] != int(match[1]): rows = []
                row = dict(point=int(match[1]), iteration=int(match[2]), limit=int(match[3]))
            for label, key in (('New End-Old End Dist.', 'dx'),
                               ('Recorrection delta-x convergence threshold', 'threshold')):
                if row and label in line:
                    found = re.search(NUMBER, line.split('=')[-1] if '=' in line else line.split(':')[-1])
                    if found: row[key] = _number(found[0])
            match = re.search(r'(?:Total number of points:|# OF POINTS ALONG THE PATH\s*=)\s*(\d+)', line)
            if match: result['points'] = int(match[1])
            match = re.search(r'NET REACTION COORDINATE UP TO THIS POINT\s*=\s*('+NUMBER+')', line)
            if match: result['coordinate'] = _number(match[1])
            if 'PES minimum detected' in line: result['stop'] = 'gradient_minimum'
            if 'Magnitude of the gradient =' in line:
                result['gradient'] = _number(line.split('=')[-1].strip())
            if 'Maximum number of corrector steps' in line: result['stop'] = 'corrector_limit'
            if 'maximum number of points' in line.lower(): result['stop'] = 'point_limit'
            if 'Error termination' in line: result['irc_termination'] = 'error'
            if 'Normal termination' in line: result['irc_termination'] = 'normal'
        if 'Normal termination' in line or 'Error termination' in line:
            if section == 0 and 'irc_geometry' not in result and last_geometry:
                result['irc_geometry'] = last_geometry
            if section > 0 and last_geometry: result['post_optimization_geometry'] = last_geometry
            section += 1
    if row: rows.append(row)
    result['last_corrector'] = rows
    result['classification'] = classify(result)
    first, last = result.get('irc_geometry'), result.get('post_optimization_geometry')
    if first and last:
        result['irc_components'] = components(first)
        result['optimized_components'] = components(last)
        result['possible_fragment_recapture'] = len(result['irc_components']) > len(result['optimized_components'])
    return result


def components(geometry):
    """Covalent-radius distance screen, with explicit atom-index membership."""
    from ase.data import covalent_radii
    pending, result = set(range(len(geometry))), []
    while pending:
        stack = [pending.pop()]; group = []
        while stack:
            i = stack.pop(); group.append(i)
            for j in list(pending):
                distance = math.dist(geometry[i][1:], geometry[j][1:])
                if distance < 1.25*(covalent_radii[geometry[i][0]]+covalent_radii[geometry[j][0]]):
                    pending.remove(j); stack.append(j)
        result.append(sorted(group))
    return sorted(result)


def classify(report):
    if report.get('stop') == 'gradient_minimum' and report.get('irc_termination') == 'normal':
        return 'early_gradient_stop' if report['points'] <= 1 else 'gradient_endpoint_review'
    if report.get('stop') != 'corrector_limit': return 'unclassified'
    rows = [r for r in report['last_corrector'] if 'dx' in r and r.get('threshold', 0) > 0]
    if len(rows) < 40: return 'insufficient_history'
    last = rows[-1]
    # An extension is justified only by recent contraction near the unchanged
    # threshold, not by the presence of a generic iteration-limit message.
    if (last['iteration'] == last['limit'] == 100 and 1 < last['dx']/last['threshold'] <= 2
            and median(r['dx'] for r in rows[-20:]) < .85*median(r['dx'] for r in rows[-40:-20])):
        return 'contracting_near_threshold'
    if median(r['dx'] for r in rows[-20:]) >= .85*median(r['dx'] for r in rows[-40:-20]):
        return 'stagnating_corrector'
    return 'unclassified'


def source(graph, ident, filenames):
    """Require matching task/job generations and an inactive allocation."""
    from kinbot.anl import dispatch
    graph, spec, state = dispatch._load(graph)
    task = next(t for t in spec['tasks'] if t['id'] == ident)
    entry = state['tasks'].get(ident, {})
    directory = graph/'tasks'/ident
    outcome = json.loads((directory/'execution.json').read_text())
    job = entry.get('job_id')
    if (not job or str(outcome.get('slurm_job_id')) != str(job)
            or entry.get('status') not in ('complete', 'failed', 'submitted')
            or outcome.get('status') not in ('executed', 'failed')):
        raise ValueError('Probe source has no matching finished allocation receipt.')
    if dispatch._job_active(str(job)): raise ValueError('Probe source allocation is still active.')
    dispatch._verify_stage_files(graph, task, entry)
    dispatch._verify_execution(graph, task, entry, outcome)
    names = set(filenames) | {'execution.json', 'task.json', 'job.slurm', task['input_name'], 'geometry.xyz'}
    paths = [directory/n for n in names]
    if any(p.parent != directory or p.resolve() != p or not p.is_file() for p in paths):
        raise ValueError('Probe evidence must be regular task files, not external paths.')
    return spec, task, outcome, dict(graph=str(graph), task=ident, job_id=str(job),
        files={str(p.relative_to(graph)):file_sha256(p) for p in paths})


def verify(proof):
    graph = Path(proof['graph'])
    for name, checksum in proof['files'].items():
        p = graph/name
        if p.resolve() != p or file_sha256(p) != checksum:
            raise ValueError('Native probe evidence changed; keep the previous decision for review.')
    from kinbot.anl import dispatch
    _, _, state = dispatch._load(graph)
    entry = state['tasks'].get(proof['task'], {})
    if str(entry.get('job_id')) != proof['job_id'] or dispatch._job_active(proof['job_id']):
        raise ValueError('Native probe source allocation changed or is active.')


def irc_probe(graph, ident, recipe):
    """One new raw-path task at the accepted parent, with a fixed trial budget."""
    from kinbot.anl import checkpoint_inputs, dispatch
    from kinbot.manual_fallback import logical_task
    graph = Path(graph).resolve()
    _, spec, _ = dispatch._load(graph)
    task = next(t for t in spec['tasks'] if t['id'] == ident)
    spec, task, outcome, proof = source(graph, ident, [task['stdout']])
    if task.get('backend') != 'gaussian' or task.get('native_path_probe'):
        raise ValueError('Only original Gaussian branch attempts can seed a probe.')
    # Reviewed numerical settings are based on the effective deck, not an
    # earlier logical deck hidden by a resource/manual fallback wrapper.
    logical_task(task)  # validates any wrapper without dropping its edits
    text = task['input_template']; first = text.split('--Link1--')[0]
    match = re.search(r'IRC=\(([^()]*)\)', first, re.I)
    if not match: raise ValueError('Expected a single explicit IRC route.')
    options = match[1].split(','); lower = [s.lower() for s in options]
    required = {'rcfc', 'gradientonly', 'eulerpc', 'massweighted', 'tight', 'report=cartesians', 'maxcycle=100'}
    if not required <= set(lower) or ('forward' in lower) == ('reverse' in lower):
        raise ValueError('Probe requires the documented RCFC EulerPC mass-weighted route.')
    if (len(lower) != 10 or not any(re.fullmatch(r'stepsize=-[1-9]\d*', s) for s in lower)
            or not any(re.fullmatch(r'maxpoints=[1-9]\d*', s) for s in lower)):
        raise ValueError('Unknown, repeated or exhausted IRC settings require separate review.')
    if len(re.findall(r'(?im)^%chk=', first)) != 1:
        raise ValueError('Probe needs one branch checkpoint.')
    directory = graph/'tasks'/ident
    _, imports = checkpoint_inputs.plan(directory, task['input_name'], (directory/task['input_name']).read_text())
    if len(imports) != 1: raise ValueError('Probe requires one verified dependency Hessian.')
    for binding in imports:
        p = Path(binding['source'])
        proof['files'][str(p.relative_to(graph))] = binding['sha256']
        proof['files'][str((p.parent/'execution.json').relative_to(graph))] = binding['source_execution_sha256']
    with (directory/task['stdout']).open(errors='replace') as stream: report = gaussian_path(stream)
    if recipe == 'pc150':
        if report['classification'] != 'contracting_near_threshold':
            raise ValueError('No near-threshold contraction supports an iteration extension.')
        options[lower.index('maxcycle=100')] = 'MaxCycle=150'
    elif recipe == 'dvv':
        if report['classification'] not in ('stagnating_corrector', 'early_gradient_stop', 'gradient_endpoint_review'):
            raise ValueError('No applicable stagnation or endpoint evidence for a DVV probe.')
        options[lower.index('eulerpc')] = 'DVV'
        options[lower.index('tight')] = 'VeryTight'
    else: raise ValueError('Unknown bounded path recipe.')
    new_id = ident+'_'+recipe+'_v1'
    first = first[:match.start(1)]+','.join(options)+first[match.end(1):]
    first = re.sub(r'(?im)^%chk=.*$', '%chk='+new_id+'.chk', first).rstrip()+'\n\n'
    new = deepcopy(task)
    for key in ('manual_fallback', 'geometry_output', 'result_parser', 'output_aliases'):
        new.pop(key, None)
    new.update(id=new_id, input_name=new_id+'.com', stdin=new_id+'.com',
        stdout=new_id+'.log', stderr=new_id+'.err', input_template=first,
        required_outputs=[new_id+'.log',new_id+'.chk'],
        success_marker=dict(file=new_id+'.log',contains='Normal termination of Gaussian'),
        failure_markers=[dict(file=new_id+'.log',contains='Error termination')],
        native_path_probe=dict(schema=1, source=ident, recipe=recipe, evidence=digest(proof),
                               budget=1, scientific_approval=False))
    verify(proof)
    return dict(proof=proof, report=report, task=new)


def molpro_modes(text):
    """Signed modes from Molpro's final native normal-mode table."""
    modes, kind = [], None
    for line in text.splitlines():
        if re.match(r'\s*Normal Modes\s*$', line): modes, kind = [], 'real'
        elif 'Normal Modes of imaginary frequencies' in line: kind = 'imaginary'
        elif 'Normal Modes of low/zero frequencies' in line: kind = 'zero'
        elif 'Wavenumbers [cm-1]' in line and kind in ('real','imaginary'):
            modes.extend((-1 if kind == 'imaginary' else 1)*_number(n)
                         for n in re.findall(NUMBER, line.split(']')[-1]))
    return modes


def curvature_probe(graph, ident):
    """Diagnose a rejected native Molpro saddle at its saved geometry.

    Two finite-difference steps and a central gradient use tighter electronic
    convergence. No optimization, no accepted geometry export, no mode gate
    relaxation. The source remains failed and its descendants remain blocked.
    """
    from ase.io import read
    from kinbot.anl.tasks import MOLPRO_HEADER, molpro_task
    spec, task, outcome, proof = source(graph, ident, ['optimized.xyz', ident+'.out'])
    if (task.get('backend') != 'molpro'
            or outcome.get('error') != 'Molpro transition state must contain exactly one imaginary mode.'):
        raise ValueError('Curvature probe requires a native harmonic saddle rejection.')
    text = task['input_template']
    # Accept only the generated, explicitly defined single-reference method.
    # This is a method applicability check, not a molecule/paper identifier.
    blocks = re.findall(r'(?im)^\{uccsd\(t\)-f12b,scale_trip=1,thrden=1\.d-8,thrvar=1\.d-10\}\s*$', text)
    basis = re.findall(r'(?im)^basis=([\w()-]+)\s*$', text)
    if len(blocks) != 2 or len(basis) != 1 or text.count('\nrhf\n') != 2:
        raise ValueError('Curvature probe requires an explicit generated ROHF-UCCSD(T)-F12b(T*) deck.')
    prefix = (MOLPRO_HEADER+'gthresh,energy=1.d-10,coeff=1.d-8,civec=1.d-8\n'
              +f'basis={basis[0]}\nrhf\n'+blocks[0].strip()+'\n')
    suffix = 'rhf\n'+blocks[1].strip()+'\nfrequencies,numerical\n'
    if (not text.startswith(prefix) or not text.endswith(suffix)
            or not re.fullmatch(r'\{optg,[^{}]*savexyz=\x27optimized.xyz\x27[^{}]*\}\s*',
                                text[len(prefix):-len(suffix)])):
        raise ValueError('Unknown electronic/geometry controls require separate curvature review.')
    directory = Path(graph)/'tasks'/ident
    modes = molpro_modes((directory/(ident+'.out')).read_text(errors='replace'))
    if sum(f < 0 for f in modes) <= 1:
        raise ValueError('Native table does not confirm multiple imaginary modes.')
    atoms = read(directory/'optimized.xyz')
    if atoms.get_chemical_symbols() != spec['molecule']['symbols'] or len(modes) != 3*len(atoms)-6:
        raise ValueError('Curvature seed atom order or mode count differs from the source.')
    molecule = deepcopy(spec['molecule']); molecule['positions'] = atoms.positions.tolist()
    tasks = []
    for name, end in [('gradient', 'forces,numerical,central,dstep=0.005,rstep=0.005,startcmd=rhf'),
                      ('hessian_010', 'frequencies,numerical,step=0.01,new,print=0'),
                      ('hessian_005', 'frequencies,numerical,step=0.005,new,print=0')]:
        native = molpro_task(name, f'basis={basis[0]}\nrhf\nuccsd(t)-f12b,scale_trip=1\n'+end+'\n',
            geometry_from='initial', scf_convergence=1.e-12, cc_convergence=1.e-11)
        native['resources'] = deepcopy(task['resources'])
        native['failure_markers'] = [dict(file=name+'.out', contains=s) for s in
            ('No convergence', 'NO CONVERGENCE', 'CONVERGENCE NOT REACHED')]
        # Each finite diagnostic is attempted once; no result parser promoting
        # it to an accepted TS, nor a CC-iteration fallback on these probes.
        tasks.append(native)
    verify(proof)
    return dict(proof=proof, modes=modes, spec=dict(schema=1,
        name=spec['name']+'_curvature_v1', molecule=molecule, limits=deepcopy(spec['limits']), tasks=tasks,
        intent=dict(purpose='Fixed-geometry numerical curvature/gradient audit of a rejected saddle.',
                    source=proof, scientific_approval=False, integration='none')))
