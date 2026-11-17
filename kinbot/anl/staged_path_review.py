"""Evidence-bound release of an unattempted, pre-staged half-step IRC.

The owning controller supplies the admission allowlist. This module neither
submits nor changes tasks, parent checkpoints, scientific gates or attempts.
"""
from copy import deepcopy
from pathlib import Path
import re

from kinbot.anl import checkpoint_inputs, dispatch, path_review
from kinbot.anl.results import parse_result
from kinbot.execution import file_sha256


class Inapplicable(ValueError):
    """Native evidence calls for a different investigation, not this trial."""


def _route(task):
    if task.get('backend') != 'gaussian' or task.get('native_path_probe'):
        raise Inapplicable('Requires an original Gaussian IRC branch.')
    first = task.get('input_template', '').split('--Link1--')[0]
    found = re.search(r'IRC=\(([^()]*)\)', first, re.I)
    if not found:
        raise Inapplicable('No explicit IRC route.')
    options = found[1].lower().split(',')
    required = {'rcfc', 'gradientonly', 'eulerpc', 'massweighted', 'tight',
                'report=cartesians', 'maxcycle=100'}
    if (len(options) != 10 or not required <= set(options)
            or ('forward' in options) == ('reverse' in options)):
        raise Inapplicable('Requires the unchanged bounded EulerPC route.')
    step = next((int(s.split('=')[1]) for s in options
                 if re.fullmatch(r'stepsize=-[1-9]\d*', s)), None)
    points = next((int(s.split('=')[1]) for s in options
                   if re.fullmatch(r'maxpoints=[1-9]\d*', s)), None)
    if step is None or points is None:
        raise Inapplicable('Missing spatial step or point budget.')
    return step, points


def _rename(value, old, new):
    if isinstance(value, str): return value.replace(old, new)
    if isinstance(value, list): return [_rename(v, old, new) for v in value]
    if isinstance(value, dict): return {k:_rename(v, old, new) for k,v in value.items()}
    return value


def prepare(graph, coarse_id, fine_id):
    """Approve one untried spatial half-step using the finished coarse log.

    Corrector exhaustion supports changing the spatial increment. A completed
    multi-point path with subsequent fragment recapture supports a step-size
    sensitivity check. A one-point gradient stop needs a departure test.
    """
    graph, spec, state = dispatch._load(graph)
    tasks = {t['id']:t for t in spec['tasks']}
    coarse, fine = tasks[coarse_id], tasks[fine_id]
    step, points = _route(coarse)
    if _route(fine) != (step / 2, points * 2):
        raise Inapplicable('The candidate is not the untried half-step with equal path budget.')
    expected = _rename(deepcopy(coarse), coarse_id, fine_id)
    expected['input_template'] = re.sub(r'(?i)StepSize=-\d+', 'StepSize='+str(step//2), expected['input_template'])
    expected['input_template'] = re.sub(r'(?i)MaxPoints=\d+', 'MaxPoints='+str(points*2), expected['input_template'])
    if fine != expected:
        raise Inapplicable('Candidate changes more than branch filenames and spatial sampling.')
    entry = state['tasks'].get(fine_id, {})
    directory = graph/'tasks'/fine_id
    if (entry.get('status') != 'staged' or entry.get('job_id')
            or entry.get('attempt', 1) != 1 or entry.get('previous_attempt')
            or state.get('manual_retries', {}).get(fine_id)
            or (graph/'attempts'/fine_id).exists()
            or (directory/'execution.json').exists() or list(directory.glob('*.jobid'))
            or any((directory/f).exists() for f in fine.get('required_outputs', []))):
        raise Inapplicable('Candidate has already been attempted or is owned by another action.')
    first = fine['input_template'].split('--Link1--')[0]
    for other in spec['tasks']:
        if other['id'] == fine_id: continue
        prior = state['tasks'].get(other['id'], {})
        same_route = other.get('input_template', '').split('--Link1--')[0].replace(other['id'], fine_id) == first
        if (same_route and other.get('depends_on') == fine.get('depends_on')
                and other.get('geometry_from') == fine.get('geometry_from')
                and (prior.get('job_id') or prior.get('previous_attempt')
                     or prior.get('status') in ('complete','failed','submitted','submitting'))):
            raise Inapplicable('The same spatial sampling was already attempted under another task ID.')
    dispatch._verify_stage_files(graph, fine, entry)
    _, _, outcome, proof = path_review.source(graph, coarse_id, [coarse['stdout']])
    with (graph/'tasks'/coarse_id/coarse['stdout']).open(errors='replace') as stream:
        report = path_review.gaussian_path(stream)
    rows = report['last_corrector']
    if (outcome['status'] == 'failed' and report.get('irc_termination') == 'error'
            and report['stop'] == 'corrector_limit' and rows
            and rows[-1]['iteration'] == rows[-1]['limit'] == 100):
        reason = 'corrector_limit_untried_half_step'
    elif (outcome['status'] == 'executed' and report['classification'] == 'gradient_endpoint_review'
          and report.get('possible_fragment_recapture')):
        reason = 'raw_path_endpoint_sensitivity'
    else:
        raise Inapplicable('Native path requires separate review: '+report['classification'])
    bindings = checkpoint_inputs.inspect_dependencies(directory, fine['input_name'],
        (directory/fine['input_name']).read_text())
    old_directory = graph/'tasks'/coarse_id
    old_bindings = checkpoint_inputs.inspect_dependencies(old_directory, coarse['input_name'],
        (old_directory/coarse['input_name']).read_text())
    if len(bindings) != 1 or bindings != old_bindings:
        raise ValueError('The two spatial samplings do not share one accepted Hessian.')
    binding = bindings[0]
    parent = tasks[binding['source_task']]
    parser = parent.get('result_parser', {})
    if parser.get('kind') != 'gaussian_harmonic' or parser.get('stationary_point') != 'transition_state':
        raise ValueError('The parent must have an accepted transition-state harmonic analysis.')
    _, _, _, parent_proof = path_review.source(graph, parent['id'],
        [parent['stdout'], Path(binding['source']).name])
    spectrum = parse_result((graph/'tasks'/parent['id']/parent['stdout']).read_text(errors='replace'), parser)
    # Bind the target's scientific files, not job.slurm: the scheduler may
    # legitimately refresh only its launch wrapper immediately before sbatch.
    files = {name:file_sha256(directory/name) for name in ('task.json', 'geometry.xyz', fine['input_name'])}
    result = dict(schema=1, graph=str(graph), task=fine_id, definition=deepcopy(fine),
        files=files, proofs=[proof, parent_proof], reason=reason, budget=1,
        frequencies_cm_inverse=spectrum['frequencies_cm_inverse'],
        source=coarse_id, source_report=report, scientific_approval=False)
    verify(result)
    return result


def verify(release):
    """Permit normal task state evolution, but reject changed evidence/decks."""
    graph, spec, state = dispatch._load(release['graph'])
    task = next(t for t in spec['tasks'] if t['id'] == release['task'])
    if task != release['definition']:
        raise ValueError('Released spatial sampling definition changed.')
    directory = graph/'tasks'/task['id']
    for name, checksum in release['files'].items():
        path = directory/name
        if path.resolve() != path or file_sha256(path) != checksum:
            raise ValueError('Released spatial sampling files changed.')
    for proof in release['proofs']: path_review.verify(proof)
