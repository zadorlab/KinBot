"""Preserve discovered wells across restarts; restore a saved controller inventory.

This edits only the discovery list. Existing PES owners admit restored workers
through their normal exploration limit, and still require current L2 evidence
before releasing ANL. No completion status or QC attempt is manufactured.
"""
import argparse
from contextlib import nullcontext
import json
import os
from pathlib import Path
import re
import tempfile

from kinbot import pes_phase
from kinbot.run_format import ensure_current_run
from kinbot.species_routing import is_species_name


def _names(path):
    values = path.read_text().splitlines() if path.exists() else []
    names = [v.strip() for v in values if v.strip()]
    if any(not is_species_name(n) for n in names):
        raise ValueError(f'{path}: invalid PES well name.')
    return list(dict.fromkeys(names))


def _validate(root, names):
    for name in names:
        directory = root/name
        ensure_current_run(directory)
        path = directory/(name+'.json')
        if path.is_symlink():
            raise ValueError(f'{path}: external worker inputs cannot restore discovery.')
        value = json.loads(path.read_text())
        if str(value.get('title')) != name or value.get('pes') != 1:
            raise ValueError(f'{path}: not a generated input for this PES well.')


def _write(root, names):
    fd, temporary = tempfile.mkstemp(prefix='.chemids.', dir=root)
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write('\n'.join(names)+'\n')
        os.replace(temporary, root/'chemids')
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def initialize(root, initial):
    """Seed a fresh PES, or keep its exact discovery order on restart."""
    root = Path(root).resolve()
    with pes_phase.locked(root):
        names = _names(root/'chemids')
        if initial not in names:
            names.append(initial)
        _validate(root, names)
        _write(root, names)
    return names


def summary_names(path):
    """Read only the explicit worker sections of an archived PES summary."""
    text = Path(path).read_text()
    total = re.search(r'^Total\s+(\d+)\s*$', text, re.M)
    names, section = [], False
    for line in text.splitlines():
        if line.strip() in ('Running:', 'Finished:', 'Waiting:'):
            section = True
        elif section and line.strip():
            name = line.strip()
            if not is_species_name(name):
                raise ValueError(f'{path}: invalid archived worker {name!r}.')
            names.append(name)
    if total is None or int(total[1]) != len(names) or len(set(names)) != len(names):
        raise ValueError(f'{path}: incomplete or duplicate archived PES inventory.')
    return names


def restore(root, summary, *, apply=False):
    """Append verified lost discoveries without reordering live admissions."""
    root = Path(root).resolve()
    expected = summary_names(summary)
    ensure_current_run(root)
    with pes_phase.locked(root) if apply else nullcontext():
        current = _names(root/'chemids')
        added = [n for n in expected if n not in current]
        _validate(root, [*current, *added])
        phase = root/pes_phase.NAME
        if added and phase.exists():
            record = json.loads(phase.read_text())
            if pes_phase.owner_alive(record) and record['phase'] != 'exploring':
                raise RuntimeError('PES already entered downstream planning; restore discovery '
                                   'during a controller handoff or a live exploration revisit.')
        result = dict(status='applied' if apply else 'reviewed', added=added,
                      retained=current, source=str(Path(summary).resolve()))
        if apply and added:
            _write(root, [*current, *added])
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args(argv)
    print(json.dumps(restore(args.root, args.summary, apply=args.apply), indent=2))


if __name__ == '__main__':
    main()
