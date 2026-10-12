"""Accept current calculation directories; do not migrate older results."""
import json
import os
from pathlib import Path
import re
import tempfile

from kinbot.rdkit_config import rdkit_runtime


MARKER = '.kinbot_run.json'
SCHEMA = 'kinbot.calculation.v4'
# Recognize both name lengths only when detecting existing calculation files.
# Version-one directories cannot be reused under the shorter naming format.
_SPECIES = re.compile(r'^\d+(?:-s(?:[0-9a-f]{16}|[0-9a-f]{64}))?$')
_JOB = re.compile(r'^\d+(?:-s(?:[0-9a-f]{16}|[0-9a-f]{64}))?_')


def _existing_results(directory, *, ignore_startup=False):
    for path in directory.iterdir():
        if ignore_startup and path.name in ('kinbot.log', 'pes.log'):
            continue
        if path.name in ('kinbot.db', 'kinbot.db-journal', 'kinbot.db-wal',
                         'chemids', 'kinbot.log', 'pes.log'):
            return path
        if path.is_file() and (path.suffix == '.pkl'
                or path.name.startswith('summary_') and path.suffix == '.out'
                or _JOB.match(path.name) and path.suffix in
                   {'.py', '.log', '.out', '.chk', '.fchk', '.com', '.inp'}):
            return path
        if path.is_dir() and path.name in ('conf', 'hir', 'hir_profiles', 'me',
                                           'perm', 'vrctst', 'aie', 'molpro', 'orca'):
            # Theory checks run after make_dirs, which creates empty trees
            # such as perm/conf. Those directories are not calculations.
            if (any(child.is_file() for child in path.rglob('*')) if ignore_startup
                    else any(path.iterdir())):
                return path
        if path.is_dir() and _SPECIES.fullmatch(path.name):
            # PES workers can already contain results before the parent has
            # written a database. Never authorize them as a fresh run.
            if (path.is_symlink() or (path / MARKER).exists()
                    or _existing_results(path, ignore_startup=ignore_startup) is not None):
                return path
    return None


def ensure_current_run(directory='.', *, create=False):
    """Check the format before any log rotation, result read or QC submission.

    Only a fresh directory can receive a new marker. Read-only PES processing
    must pass create=False. A restart requires the same RDKit version and modes.
    """
    directory = Path(directory)
    if directory.is_symlink():
        raise ValueError(f'{directory}: calculation-directory aliases are not supported; '
                         'use the actual current-format calculation directory.')
    marker = directory / MARKER
    expected = {'schema': SCHEMA, **rdkit_runtime()}
    if marker.exists():
        try:
            recorded = json.loads(marker.read_text())
        except (OSError, ValueError) as error:
            raise ValueError(f'{marker}: invalid calculation-format record; use a fresh directory.') from error
        runtime = dict(recorded) if isinstance(recorded, dict) else {}
        profiled = runtime.pop('profiled_theory', True)
        if runtime != expected or profiled is not True:
            raise ValueError(f'{marker}: incompatible calculation format or RDKit settings. '
                             f'Found {recorded}; required {expected}. Use the original environment '
                             'for a current-format restart, or start in a fresh directory.')
        return recorded
    if not directory.is_dir():
        if not create:
            raise ValueError(f'{directory}: current calculation directory does not exist.')
        directory.mkdir(parents=True, exist_ok=True)
    existing = _existing_results(directory)
    if existing is not None and marker.exists():
        # Another current worker can publish its marker and start writing
        # results after our first check. Validate that marker before refusal.
        return ensure_current_run(directory)
    if existing is not None or not create:
        detail = f' Existing calculation data: {existing}.' if existing is not None else ''
        raise ValueError(f'{directory}: missing {MARKER}; older calculations are incompatible.'
                         f'{detail} Start a new calculation in a fresh directory; no files were changed.')
    # Publish a complete marker without overwriting a concurrent creator.
    fd, temporary = tempfile.mkstemp(prefix='.kinbot_run_', dir=directory)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(expected, stream, indent=2, sort_keys=True)
            stream.write('\n')
        try:
            os.link(temporary, marker)
        except FileExistsError:
            return ensure_current_run(directory)
    finally:
        os.unlink(temporary)
    return expected


def record_profiled_run(directory='.'):
    """Remember profile use even if its custom record path is removed from input."""
    recorded = ensure_current_run(directory, create=True)
    if recorded.get('profiled_theory') is True:
        return
    recorded['profiled_theory'] = True
    directory = Path(directory)
    fd, temporary = tempfile.mkstemp(prefix='.kinbot_run_', dir=directory)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(recorded, stream, indent=2, sort_keys=True)
            stream.write('\n')
        os.replace(temporary, directory / MARKER)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
