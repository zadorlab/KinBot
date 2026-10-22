"""One-time recovery of verified native-interface defects, not QC failures.

The owner calls this under its drive lock. The existing live-control journal
archives only the selected attempt, fences old workers with the manifest
lock, and survives interruption. Scientific task/workflow identities do not
change; complete and allocated tasks are never reset.
"""
from copy import deepcopy
import json
import re

from ase.io import read

from kinbot.anl.cfour import basis_keyword, normalize_cfour_zmat


def _cfour_alias(task, directory, outcome):
    if task.get('backend') != 'cfour' or task.get('kind') != 'external':
        return None
    old = (directory/task['input_name']).read_text()
    if normalize_cfour_zmat(old) == old:
        return None
    # Restrict recovery to the documented AUG-PV/PCV spelling fix, not to
    # arbitrary edits or another normalization of the keyword block.
    aliases = [m for m in re.findall(r'(?im)^\s*BASIS\s*=\s*(\S+?)\s*[,)]?\s*$', old)
               if basis_keyword(m) != m]
    if len(aliases) != 1:
        return None
    if outcome is not None:
        native = directory/task.get('stdout', 'stdout.txt')
        if not native.is_file():
            return None
        text = native.read_text(errors='replace')
        if (outcome.get('status') != 'failed' or not re.search(
                rf'Option\s+{re.escape(aliases[0])}\s+for keyword BASIS\s+not known',
                text, re.I)):
            return None
    return 'cfour_augmented_basis_alias_v1'


def _molpro_geometry(task, directory, outcome):
    profile = task.get('profile', {})
    if (task.get('kind') != 'ase_optimize' or profile.get('calculator') != 'molpro'
            or profile.get('method') != 'CCSD(T)-F12b'
            or profile.get('calculator_kwargs', {}).get('effective_method') != 'CCSD-F12b'
            or outcome.get('status') != 'failed'
            or outcome.get('error') != 'Molpro output lacks the requested undisplaced F12b energy.'):
        return None
    from kinbot.ase_modules.calculators.molpro import parse_output, parse_numerical_gradient
    # This defect fails the very first force evaluation. Do not reset later
    # optimizer progress or treat a missing/invalid gradient as this defect.
    files = list(directory.glob(task['id']+'_step_*.out'))
    native = directory/(task['id']+'_step_0001.out')
    geometry = native.with_suffix('.xyz')
    if files != [native] or not geometry.is_file():
        return None
    text = native.read_text(errors='replace')
    if not re.search(r'^\s*!RHF-UCCSD-F12b\s+energy\s+', text, re.M | re.I):
        return None
    try:
        parse_output(native, method='CCSD(T)-F12b', effective_method='CCSD-F12b',
                     require_precise_energy=True)
        parse_numerical_gradient(native, geometry, read(directory/'geometry.xyz'))
    except ValueError:
        return None
    return 'molpro_ccsd_f12b_final_label_v1'


def _gaussian_linearity(task, directory, outcome, molecule):
    """Restage old linear VPT2 inputs with point-group detection preserved."""
    if (task.get('kind') != 'external' or task.get('backend') != 'gaussian'
            or task.get('result_parser', {}).get('kind') != 'gaussian_vpt2'):
        return None
    old = (directory/task['input_name']).read_text()
    none = r'(?i)\bSymm(?:etry)?\s*=\s*(?:\(\s*None\s*\)|None\b)'
    if not re.search(none, old):
        return None
    from kinbot.anl.dispatch import _render_input
    rendered = _render_input(task['input_template'], read(directory/'geometry.xyz'),
                             molecule, task['resources'], 'gaussian')
    # Rendering only makes this substitution for verified linear geometries.
    # All other native input bytes must match; no changed policy is smuggled
    # into a restart described as a representation-only interface repair.
    if rendered != re.sub(none, 'NoSymm', old) or rendered == old:
        return None
    if outcome is not None:
        native = directory/task['result_parser']['file']
        if not native.is_file() or outcome.get('status') != 'failed':
            return None
        text = native.read_text(errors='replace')
        if ('Molecule is non-linear but NTrRo=5' not in text
                or 'Error termination via Lnk1e' not in text
                or 'Normal termination of Gaussian' in text):
            return None
    return 'gaussian_vpt2_linear_framework_v1'


def recover(graph, spec, state):
    from kinbot.anl import dispatch, live_control
    from kinbot.anl.salvage import descendants
    from kinbot.anl.timeout_recovery import _job_id

    reset, reasons = set(), {}
    updated = deepcopy(state)
    for task in spec['tasks']:
        ident = task['id']
        entry = state['tasks'].get(ident, {})
        if entry.get('status') not in ('failed', 'unavailable', 'staged'):
            continue
        directory = graph/'tasks'/ident
        path = directory/'execution.json'
        outcome = json.loads(path.read_text()) if path.is_file() else None
        if entry['status'] == 'staged' and (entry.get('job_id') or outcome is not None):
            continue
        if entry['status'] != 'staged' and outcome is None:
            continue
        kind = (_cfour_alias(task, directory, outcome)
                or _gaussian_linearity(task, directory, outcome, spec['molecule'])
                or (_molpro_geometry(task, directory, outcome) if outcome else None))
        if not kind or kind in state.get('native_interface_repairs', {}).get(ident, {}):
            continue
        job = _job_id(directory, ident, entry)
        if entry['status'] != 'staged' and (not job or dispatch._job_active(job)):
            continue
        dispatch._verify_stage_files(graph, task, entry)
        if outcome is not None:
            dispatch._verify_execution(graph, task, entry, outcome)
        children = descendants(spec, ident)
        if any(state['tasks'].get(n, {}).get('status') in
               ('complete', 'submitted', 'submitting', 'staged') for n in children):
            # A reviewed dependency migration is required; never invalidate
            # a derived result automatically for an interface repair.
            continue
        reset.add(ident)
        reset.update(n for n in children if state['tasks'].get(n, {}).get('status')
                     in ('blocked', 'skipped'))
        reasons[ident] = kind
        updated.setdefault('native_interface_repairs', {}).setdefault(ident, {})[kind] = {
            'previous_job_id': job, 'previous_status': entry['status']}
    if reset:
        live_control.commit(graph, spec, updated, reset, reasons)
        return dispatch._load(graph)[2]
    return state
