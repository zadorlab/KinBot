"""One isolated, tighter L2 minimum after native VPT2 framework exhaustion.

The original geometry and all its electronic-energy descendants remain intact.
Only VPT2 moves to the new minimum, after a separate harmonic validation. A new
geometry generation gets the whole symmetry ladder, even if its coordinates
round to the same values as the original minimum.
"""
from copy import deepcopy
import json
from pathlib import Path
import re

from kinbot.execution import atomic_json
from kinbot.anl import vpt2_recovery as symmetry


def _tasks(spec, task):
    """Only generated, explicitly understood routes can seed a new optimizer."""
    parser = task['result_parser']
    template = symmetry.replace_symmetry(task['input_template'], '')
    route = symmetry._ROUTE.search(template)
    if route is None or task.get('command') != ['g16']:
        raise ValueError('Tighter VPT2 geometry requires a supported Gaussian command/route.')
    remaining = route[0]
    expressions = [r'^\s*#[pnt]?\s*',
        re.escape(parser['method']+'/'+parser['basis']), r'\bFreq(?:uencies)?=Anharmonic\b',
        r'\bSCF\s*=\s*(?:XQC|\(XQC,Conver=\d+\))',
        r'\bEmpiricalDispersion='+re.escape(parser.get('dispersion', ''))+r'\b',
        r'\bIntegral=UltraFine\b']
    for expression in expressions:
        remaining = re.sub(expression, '', remaining, flags=re.I)
    if remaining.strip() or not re.search(r'\bIntegral=UltraFine\b', route[0], re.I):
        raise ValueError('Custom VPT2 route needs a reviewed geometry-refinement profile.')
    scf = re.search(r'\bConver=(\d+)', route[0], re.I)
    convergence = max(10, int(scf[1]) if scf else 0)
    tight = re.sub(r'\bSCF\s*=\s*(?:XQC|\(XQC,Conver=\d+\))',
                   f'SCF=(XQC,Conver={convergence})', template, flags=re.I)
    if not re.search(r'\bSCF=', tight, re.I):
        raise ValueError('Explicit SCF control is required for L2 refinement.')
    source = task.get('geometry_from', 'initial')
    original = next((t for t in spec['tasks'] if t['id'] == source), {})
    fmax = min(1.e-5, original.get('optimizer', {}).get('fmax', .0005)/10)
    ident = task['id']
    geom, harmonic = ident+'_refine_geometry', ident+'_refine_harmonic'
    if {geom, harmonic} & {t['id'] for t in spec['tasks']}:
        raise ValueError('VPT2 refinement task names already exist without a recovery receipt.')
    kwargs = dict(Symm='None', scf=f'(XQC,Conver={convergence})', integral='UltraFine')
    if parser.get('dispersion'):
        kwargs['EmpiricalDispersion'] = parser['dispersion']
    optimization = dict(id=geom, kind='ase_optimize', geometry_from=source,
        geometry_output='final.xyz', resources=deepcopy(task['resources']),
        profile=dict(calculator='gaussian', method=parser['method'], basis=parser['basis'],
                     command='g16', calculator_kwargs=kwargs, optimizer='sella'),
        optimizer=dict(fmax=fmax, steps=160, step_retries=0, order=0,
                       sella_kwargs={'internal': len(spec['molecule']['symbols']) > 2}))
    hessian = deepcopy(task)
    for key in ('recovery', 'manual_fallback', 'failure_policy', 'failure_skips', 'component_families'):
        hessian.pop(key, None)
    hessian.update(id=harmonic, geometry_from=geom, depends_on=[geom],
        input_template=re.sub(r'\bFreq(?:uencies)?=Anharmonic\b', 'Freq', tight, flags=re.I),
        result_parser=dict(parser, kind='gaussian_harmonic', natoms=len(spec['molecule']['symbols'])))
    replacement = deepcopy(task)
    replacement.pop('recovery', None)
    replacement.update(geometry_from=geom, depends_on=[*task.get('depends_on', []), harmonic],
                       input_template=tight)
    replacement['vpt2_geometry_refinement'] = dict(schema=1, geometry_task=geom,
        harmonic_task=harmonic, original_geometry_from=source, fmax=fmax,
        scf_convergence=convergence, maximum_geometry_refinements=1)
    if 'manual_fallback' in replacement:
        from kinbot.manual_fallback import logical_task, digest
        logical_task(task)
        replacement['manual_fallback']['patched_sha256'] = digest(replacement)
    return optimization, hessian, replacement


def recover(graph, spec, state):
    from kinbot.anl import dispatch, live_control
    from kinbot.anl.refinement import require_open
    from kinbot.anl.salvage import descendants
    from kinbot.anl.timeout_recovery import _job_id
    graph = Path(graph)
    if spec.get('intent', {}).get('stationary_point') != 'minimum':
        return spec, state  # No inferred minimum and no automatic saddle VPT2.
    for task in list(spec['tasks']):
        ident = task['id']
        entry = state['tasks'].get(ident, {})
        if (task.get('result_parser', {}).get('kind') != 'gaussian_vpt2'
                or entry.get('status') not in ('failed', 'unavailable')
                or entry.get('vpt2_refinement_exhausted')):
            continue
        directory = graph/'tasks'/ident
        native = directory/task['result_parser']['file']
        if not native.is_file():
            continue
        try:
            history, linear = symmetry.context(graph, task, entry)
            symmetry.plan(task, native.read_text(errors='replace'), history=history, linear=linear)
        except symmetry.VPT2RecoveryExhausted:
            pass
        except symmetry.VPT2RecoveryUnavailable:
            continue
        else:
            continue  # Ordinary fallbacks have precedence.
        job = _job_id(directory, ident, entry)
        if job and dispatch._job_active(job):
            continue
        dispatch._verify_stage_files(graph, task, entry)
        execution = json.loads((directory/'execution.json').read_text())
        dispatch._verify_execution(graph, task, entry, execution)
        if execution['status'] != 'failed' or (execution.get('slurm_job_id') and not job):
            raise ValueError('VPT2 refinement requires a verified failed allocation.')
        if task.get('vpt2_geometry_refinement'):
            entry.update(status='unavailable', vpt2_refinement_exhausted=True,
                unavailable_reason='Native framework failure persisted through both symmetry ladders and one tighter L2 minimum.',
                component_families=list(task.get('component_families', ['vpt2_correction'])))
            # Retain job ID, native evidence and all attempts for audit. This
            # explicit bounded terminal state also applies under hold_failures.
            continue
        require_open(graph)
        children = descendants(spec, ident)
        if any(state['tasks'].get(n, {}).get('status') not in
               (None, 'waiting', 'blocked', 'skipped') for n in children):
            raise RuntimeError('VPT2 results have downstream consumers; reopen their owner before refinement.')
        try:
            geometry, harmonic, replacement = _tasks(spec, task)
        except ValueError as error:
            entry['refinement_review'] = str(error)
            continue
        updated = deepcopy(spec)
        updated['tasks'][updated['tasks'].index(task)] = replacement
        updated['tasks'].extend((geometry, harmonic))
        dispatch.assign_partitions(updated)
        dispatch.validate_spec(updated)
        atomic_json(graph/'state.json', state)
        reset = {ident} | {n for n in children if state['tasks'].get(n, {}).get('status') in ('blocked', 'skipped')}
        try:
            live_control.commit(graph, updated, state, reset,
                {n: 'One tighter L2 minimum and harmonic Hessian; restart VPT2 symmetry ladder.' for n in reset})
        except RuntimeError as error:
            if 'legacy workers' not in str(error):
                raise
            continue
        _, spec, state = dispatch._load(graph)
    return spec, state
