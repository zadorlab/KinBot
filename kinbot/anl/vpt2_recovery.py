"""Bounded Gaussian framework recovery, selected from evidence and history.

These are symmetry-processing alternatives, not changes of geometry or theory.
The dispatcher owns archival, allocation fencing and result validation.
"""
import json
from pathlib import Path
import re

import numpy as np
from ase.io import read


class VPT2RecoveryUnavailable(RuntimeError):
    """No untried, diagnostic-specific framework recovery applies."""


_ROUTE = re.compile(r'(?m)^[ \t]*#[^\n]*(?:\n[^\S\n]*[^\s\n][^\n]*)*')
_SYMMETRY = re.compile(
    r'(?i)(?<!\S)(?:NoSymm(?:etry)?|Symm(?:etry)?\s*=\s*'
    r'(?:\([^)]*\)|PG\s*=\s*[A-Za-z0-9*]+|None))(?!\S)')
_PG = re.compile(r'PG=([A-Z0-9*]+)(,LOOSE)?\Z')


def symmetry(text):
    """Return a canonical setting and its span, only from the route section."""
    route = _ROUTE.search(text)
    if route is None:
        raise VPT2RecoveryUnavailable('Gaussian route section is missing.')
    matches = list(_SYMMETRY.finditer(text, route.start(), route.end()))
    if len(matches) > 1:
        raise VPT2RecoveryUnavailable('Input has multiple symmetry keywords.')
    if not matches:
        if re.search(r'(?i)(?<!\S)(?:No)?Symm(?:etry)?\b', route[0]):
            raise VPT2RecoveryUnavailable('Unrecognized symmetry options require review.')
        return 'DEFAULT', None
    match = matches[0]
    value = re.sub(r'\s+', '', match[0]).upper()
    if value in ('NOSYMM', 'NOSYMMETRY'):
        return 'NOSYMM', match.span()
    value = value.split('=', 1)[1].strip('()')
    if value != 'NONE' and not _PG.fullmatch(value):
        raise VPT2RecoveryUnavailable('Unrecognized symmetry options require review.')
    return value, match.span()


def replace_symmetry(template, keyword):
    _, span = symmetry(template)
    if span:
        return template[:span[0]] + keyword + template[span[1]:]
    # Insert after the route prefix, retaining compound Freq options verbatim.
    return re.sub(r'(?im)^([ \t]*#(?:[pnt](?=\s))?)(?=\s)',
                  lambda m: m[1] + ' ' + keyword, template, count=1)


def settings(task):
    """Read new receipts and conservatively migrate legacy retry receipts."""
    receipt = task.get('recovery', {})
    groups = receipt.get('attempted_framework_groups', [])
    if (not isinstance(groups, list) or any(not isinstance(g, str)
            or not re.fullmatch(r'[A-Za-z0-9*]+', g) for g in groups)):
        raise VPT2RecoveryUnavailable('Invalid VPT2 framework recovery history.')
    previous = receipt.get('attempted_symmetry_settings', [])
    if (not isinstance(previous, list) or any(not isinstance(s, str)
            or (s not in ('DEFAULT', 'NONE', 'NOSYMM') and not _PG.fullmatch(s))
            for s in previous)):
        raise VPT2RecoveryUnavailable('Invalid VPT2 symmetry recovery history.')
    tried = list(previous)
    tried.extend('PG=' + g.upper() for g in groups)
    # The old native-PG recovery was reachable only after the None workaround.
    # Its receipt did not explicitly record that earlier setting.
    if 'attempted_symmetry_settings' not in receipt and receipt.get('assigned_native_framework'):
        tried.append('NONE')
    tried.append(symmetry(task['input_template'])[0])
    return list(dict.fromkeys(tried))


def framework_history(task):
    groups = list(dict.fromkeys(_PG.fullmatch(s)[1] for s in settings(task)
                               if _PG.fullmatch(s)))
    current = _PG.fullmatch(symmetry(task['input_template'])[0])
    return groups, current[1] if current else None


def context(run_dir, task, entry):
    """Include actual decks in same-geometry/theory archives, also manual retries.

    Fixed archive layouts avoid walking potentially large native scratch trees.
    Receipts remain the fallback for archives predating journaled live controls.
    """
    root = Path(run_dir)
    directory = root / 'tasks' / task['id']
    tried = settings(task)
    tried.append(symmetry((directory / task['input_name']).read_text())[0])
    archives = root / 'attempts'
    paths = [*archives.glob(f'*/graphs/tasks/{task["id"]}/task.json'),
             *archives.glob(f'{task["id"]}/*/task.json')]
    for path in sorted(paths):
        record = json.loads(path.read_text())
        prior = record.get('task', {})
        if (record.get('geometry_sha256') != entry.get('geometry_sha256')
                or prior.get('id') != task['id']
                or prior.get('result_parser') != task.get('result_parser')):
            continue
        tried.extend(settings(prior))
        native = path.parent / prior['input_name']
        if native.is_file():
            tried.append(symmetry(native.read_text())[0])
    atoms = read(directory / 'geometry.xyz')
    singular = np.linalg.svd(atoms.positions - atoms.positions.mean(axis=0),
                             compute_uv=False)
    linear = len(atoms) > 1 and singular[1] < 1.e-7
    return list(dict.fromkeys(tried)), linear


def plan(task, output, *, history=(), linear=False):
    """Use established alternatives once each, only for the native mismatch.

    Prefer the symmetry bypass, then up to two native-reported groups and one
    corroborated Loose detection attempt, then the minimal C1 cap if untried.
    Linear geometries retain group detection via NoSymm instead of None/C1.
    """
    required = ('ERROR: Inconsistency found in framework group definition:',
                'Error termination via Lnk1e')
    changes = list(re.finditer(
        r'New:[ \t]*([A-Za-z0-9*]+)[ \t]*-[ \t]*Old:[ \t]*([A-Za-z0-9*]*)',
        output, re.I))
    if (not all(marker in output for marker in required) or not changes
            or 'Normal termination of Gaussian' in output):
        raise VPT2RecoveryUnavailable(
            f'{task["id"]}: native output is not the framework-group VPT2 failure.')
    change = changes[-1]
    new, old = (change[n].upper() for n in (1, 2))
    if new == old:
        raise VPT2RecoveryUnavailable('Native framework groups do not disagree.')
    current = symmetry(task['input_template'])[0]
    if linear and current == 'NONE':
        current = 'NOSYMM'  # The renderer preserves linear classification.
    tried = set(settings(task)) | set(history)
    bypass = 'NOSYMM' if linear else 'NONE'
    if bypass not in tried:
        return ('NoSymm' if linear else 'Symmetry=None'), False, change[0]
    pg = _PG.fullmatch(current)
    # Explicit native assignment is supported when the bypass left no group,
    # or the input's current group is the old framework reported by Gaussian.
    native_groups = {_PG.fullmatch(s)[1] for s in tried if _PG.fullmatch(s)
                     and _PG.fullmatch(s)[1] != 'C1'}
    if (((current == 'NONE' and not old) or current == 'NOSYMM'
            or (pg and old == pg[1]))
            and 'PG=' + new not in tried and len(native_groups) < 2
            and (not linear or new != 'C1')):
        return f'Symmetry=(PG={new})', True, change[0]
    if (pg and pg[1] == new and old and old != new
            and not any(s.endswith(',LOOSE') for s in tried)
            and re.search(r'(?im)^\s*Full point group\s+'+re.escape(old)+r'\s', output)
            and re.search(r'(?im)^\s*Framework definition based on atomic masses:\s*'
                          +re.escape(new)+r'\s*$', output)):
        return f'Symmetry=(PG={new},Loose)', True, change[0]
    if not linear and 'PG=C1' not in tried:
        return 'Symmetry=(PG=C1)', False, change[0]
    raise VPT2RecoveryUnavailable(
        f'{task["id"]}: applicable symmetry settings already attempted; '
        'bounded framework retries exhausted.')
