"""One-time in-place repair of the reviewed downstream-failure runs.

The numerical fixes are general KinBot features. This migration preserves
accepted calculation caches and archives the obsolete attempts externally.
It never declares a failed native calculation complete.
"""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil

from kinbot.execution import atomic_json, file_sha256


class Archive:
    def __init__(self, root):
        self.root = root.resolve()
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        self.path = root.with_name(root.name + '.repair_archive.' + stamp)
        self.path.mkdir()
        self.entries = []

    def save(self, source, *, move=False, group='before'):
        source = Path(source)
        if not source.exists():
            return None
        if source.is_symlink() or not source.resolve().is_relative_to(self.root):
            raise ValueError(f'Refusing external/symlink migration: {source}')
        target = self.path / group / source.relative_to(self.root)
        if target.exists():
            if move:
                raise ValueError(f'Archive destination already exists: {target}')
            return target
        target.parent.mkdir(parents=True, exist_ok=True)
        if move:
            source.rename(target)
        elif source.is_dir():
            shutil.copytree(source, target)
        else:
            shutil.copy2(source, target)
        self.entries.append({'source': str(source), 'archive': str(target),
                             'action': 'move' if move else 'copy'})
        atomic_json(self.path / 'manifest.json', self.entries)
        return target


def update_electronic_inputs(root, record, archive):
    """Migrate reviewed electronic postprocessing, preserving numerical QC."""
    from kinbot.electronic_states import electronic_state_models
    package = Path(__file__).resolve().parents[1]
    updated = []
    archive.save(root / 'prepared.json')
    for filename in record['inputs']:
        case = Path(filename).parent
        reviewed = json.loads((package / f'{case.name}_anl1_f12'
                               / f'{case.name}.json').read_text())
        levels = reviewed['electronic_state_models']
        new_models = electronic_state_models(reviewed)
        existing = json.loads(Path(filename).read_text())
        for key, old in electronic_state_models(existing).items():
            new = new_models.get(key)
            if new is None or any(old.get(field) != new.get(field) for field in (
                    'smiles', 'ground_multiplicity', 'term_orbital_degeneracy')):
                raise ValueError('Electronic postprocessing migration requires '
                                 'unchanged spin and spatial term contracts.')
            old_so = old.get('spin_orbit_correction')
            new_so = new.get('spin_orbit_correction', {})
            if old_so is not None and any(old_so.get(field) != new_so.get(field)
                                         for field in ('value_hartree', 'reference')):
                raise ValueError('Existing SO correction differs from the reviewed model.')
        # PES regenerates its worker input on restart. Update existing child
        # and reference inputs as well so no restarted child sees stale levels.
        for path in case.rglob('*.json'):
            if path.is_symlink() or 'attempts' in path.relative_to(case).parts:
                continue
            try:
                data = json.loads(path.read_text())
            except (ValueError, OSError):
                continue
            if (not isinstance(data, dict) or 'theory_preset' not in data
                    or 'electronic_state_models' not in data
                    or data['electronic_state_models'] == levels):
                continue
            archive.save(path)
            data['electronic_state_models'] = levels
            atomic_json(path, data)
            updated.append(str(path.relative_to(root)))
    for name in record['sha256']:
        record['sha256'][name] = file_sha256(root / name)
    atomic_json(root / 'prepared.json', record)
    return updated


def archive_electronic_postprocessing(case, archive):
    """Retire derived energies/populations, keeping all accepted QC caches."""
    retired = []
    for base in case.glob('*/anl_automatic'):
        for name in ('profiles', 'formations', 'handoff.json', 'execution.json',
                     'final_report.json', 'cbh_thermochemistry_models.json',
                     'electronic_state_audit.json'):
            target = archive.save(base / name, move=True, group='old_postprocessing')
            if target is not None:
                retired.append(str(target))
        target = archive.save(base.parent / 'me', move=True, group='old_postprocessing')
        if target is not None:
            retired.append(str(target))
    return retired


def repair_graph(graph, archive):
    """Reparse native successes and migrate only uncompleted VPT2 tasks."""
    from kinbot.anl import dispatch
    for name in ('workflow.json', 'state.json'):
        archive.save(graph / name)
    for outcome in graph.glob('tasks/*/execution.json'):
        archive.save(outcome)
    # This reconciles calculations which finished after the old controller
    # failed, and recovers SH's native success using the corrected parser.
    dispatch.advance(graph, submit=False)
    _, spec, state = dispatch._load(graph)
    recovered = []
    for task in spec['tasks']:
        ident = task['id']
        status = state['tasks'].get(ident, {}).get('status')
        if status in ('submitted', 'submitting'):
            raise RuntimeError(f'{graph}/{ident}: unresolved submission; no restart launched.')
        # Optional native execution failures are recorded as unavailable by
        # advance(), so they must reach the same evidence-checked recovery.
        # Preflight unavailability has no native attempt to repair.
        if status == 'unavailable' and (
                task.get('result_parser', {}).get('kind') != 'gaussian_vpt2'
                or state['tasks'][ident].get('failure_phase') == 'preflight'
                or not (graph / 'tasks' / ident / 'execution.json').is_file()):
            continue
        if status not in ('failed', 'unavailable'):
            continue
        if task.get('result_parser', {}).get('kind') != 'gaussian_vpt2':
            raise RuntimeError(f'{graph}/{ident}: unreviewed failure; no restart launched.')
        # This function requires the specific native framework-group error;
        # it refuses unrelated Gaussian failures and archives before staging.
        old = dispatch.recover_gaussian_vpt2_symmetry(graph, ident)
        target = archive.save(old, move=True, group='failed_attempts')
        _, spec, state = dispatch._load(graph)
        state['tasks'][ident]['previous_attempt'] = str(target)
        atomic_json(graph / 'state.json', state)
        recovered.append(ident)

    _, spec, state = dispatch._load(graph)
    restage = []
    changed = False
    for task in spec['tasks']:
        ident = task['id']
        entry = state['tasks'].get(ident, {})
        if (task.get('result_parser', {}).get('kind') != 'gaussian_vpt2'
                or entry.get('status') in ('complete', 'unavailable')):
            continue
        template, count = re.subn(
            r'(?i)Symmetry\s*=\s*\(\s*PG\s*=\s*C1\s*\)',
            'Symmetry=None', task['input_template'])
        if not count:
            continue
        if entry.get('status') not in (None, 'staged', 'blocked', 'skipped'):
            raise RuntimeError(f'{graph}/{ident}: cannot migrate {entry.get("status")}.')
        directory = graph / 'tasks' / ident
        if directory.exists():
            if (directory / 'execution.json').exists():
                raise RuntimeError(f'{directory}: unreviewed execution record.')
            archive.save(directory, move=True, group='unstarted_templates')
            state['tasks'].pop(ident, None)
            restage.append(ident)
        task['input_template'] = template
        changed = True
    if changed:
        atomic_json(graph / 'workflow.json', spec)
        state['spec_sha256'] = file_sha256(graph / 'workflow.json')
        atomic_json(graph / 'state.json', state)
    state = dispatch.advance(graph, submit=False)
    if any(e['status'] in ('failed', 'blocked', 'submitted', 'submitting')
           for e in state['tasks'].values()):
        raise RuntimeError(f'{graph}: unresolved task state; no restart launched.')

    # These files remain hash-bound by successful execution records. Move
    # their content outside the live tree and preserve a read-only reference;
    # the completed task is never re-executed by this migration.
    for outcome in graph.glob('tasks/*/execution.failed.json'):
        if outcome.is_symlink():
            continue
        target = archive.save(outcome, move=True, group='parser_history')
        outcome.symlink_to(target)
    # Empty retry directories otherwise show up as old attempt folders.
    attempts = graph / 'attempts'
    if attempts.exists():
        for directory, _, _ in os.walk(attempts, topdown=False):
            path = Path(directory)
            if not any(path.iterdir()):
                path.rmdir()
    return {'graph': str(graph), 'recovered_vpt2': recovered,
            'updated_unstarted_vpt2': restage,
            'complete_tasks_preserved': sum(e['status'] == 'complete'
                                             for e in state['tasks'].values())}


def repair_before_resume(root, record):
    """Repair an inactive prepared handoff, leaving all accepted QC in place."""
    from kinbot import run_cleanup
    root = Path(root).resolve(strict=True)
    cases = [Path(path).parent for path in record['inputs']]
    if any(not case.is_relative_to(root) or case.is_symlink() for case in cases):
        raise ValueError('Prepared cases must be inside the run directory.')
    inventories = [run_cleanup.inspect(case) for case in cases]
    # Do not cancel useful in-flight QC for a saved-results restart. A live
    # controller could also rewrite the records during migration.
    if any(item['jobs'] or item['processes'] for item in inventories):
        raise RuntimeError('Saved-results restart requires both runs inactive. '
                           'Existing jobs/processes were left untouched.')
    archive = Archive(root)
    atomic_json(archive.path / 'inventory.json', inventories)
    archive.save(root / 'resume_provenance.json')
    result = {'archive': str(archive.path), 'graphs': [], 'vrc_archived': []}
    result['electronic_inputs_updated'] = update_electronic_inputs(root, record, archive)
    result['electronic_postprocessing_archived'] = []
    for case in cases:
        archive.save(case / 'resume_provenance.json')
        # A previous invocation may have updated both inputs, then stopped
        # partway through the first case. Finish invalidating derived data
        # even when this invocation finds the inputs already current.
        result['electronic_postprocessing_archived'].extend(
            archive_electronic_postprocessing(case, archive))
        # Copy database files only after all writers have stopped. Do not
        # remove unsuccessful reaction decisions or turn them into successes.
        for database in case.rglob('kinbot.db'):
            archive.save(database)
            for suffix in ('-wal', '-shm'):
                archive.save(Path(str(database) + suffix))
        for graph_file in sorted(case.glob('*/anl_automatic/**/workflow.json')):
            if 'attempts' in graph_file.relative_to(case).parts:
                continue
            result['graphs'].append(repair_graph(graph_file.parent, archive))
        for vrc in case.glob('*/vrctst'):
            # This migration is for the failed first scan, not a completed
            # VRC production run. Never silently retire accepted scan points.
            for state_file in vrc.glob('molpro_scan/*/state.json'):
                state = json.loads(state_file.read_text())
                if any(e.get('status') == 'complete' for e in state['tasks'].values()):
                    raise RuntimeError(f'{state_file}: accepted VRC point needs review before migration.')
            target = archive.save(vrc, move=True, group='invalid_vrc')
            result['vrc_archived'].append(str(target))
        historical_exits = {}
        from kinbot.qc import _worker_exit_status
        for marker in list(case.rglob('*.exitcode')):
            if marker.is_symlink() or 'anl_automatic' in marker.relative_to(case).parts:
                continue
            outcome = _worker_exit_status(marker)
            if outcome is None or outcome[0] == 0:
                continue
            # Keep the compact outcome marker and database rejection so a
            # failed TS cannot become accepted merely by removing its logs.
            # Native inputs/outputs of that failed attempt move out of view.
            archive.save(marker)
            historical_exits[str(marker.relative_to(case))] = file_sha256(marker)
            stem = marker.stem
            directories = [marker.parent]
            worker = next((p for p in marker.parents if (p / 'kinbot.db').is_file()), None)
            if worker is not None:
                directories.append(worker / 'perm' / marker.parent.relative_to(worker))
            for directory in directories:
                if not directory.is_dir():
                    continue
                for path in list(directory.iterdir()):
                    if path == marker or path.is_symlink():
                        continue
                    if (path.name.startswith(stem + '.')
                            or path.name.startswith(stem + '_sella.')
                            or path.name.startswith(stem + '_scf_recovery')
                            or path.name.startswith(stem + '_minimum_recovery')):
                        archive.save(path, move=True, group='rejected_qc')
        # Controller logs are not scientific inputs. Start new logs; retain
        # all old text unchanged outside the monitored directory.
        for directory, dirs, names in os.walk(case, followlinks=False):
            dirs[:] = [d for d in dirs if d not in ('attempts', 'orbitals', 'assets')
                       and not Path(directory, d).is_symlink()]
            for name in names:
                controller_log = (name in ('driver.log', 'driver.stdout', 'driver.stderr',
                    'kinbot.log', 'kinbot.out', 'kinbot.err', 'pes.log')
                    or re.fullmatch(r'(?:kinbot|pes)_\d{4}-.*\.log', name))
                controller_state = name in ('launch.json', 'pes_status.json',
                    'kinbot_status.json', 'pes_summary.txt', 'kinbot_monitor.out',
                    'kinbot_optimization_monitor.json', '.kinbot_post_reaction_branches.json')
                if controller_log or controller_state:
                    archive.save(Path(directory) / name, move=True, group='controller_history')
        atomic_json(case / 'resume_provenance.json', {
            'archive': str(archive.path), 'archived_failure_exits': historical_exits})
    atomic_json(archive.path / 'repair.json', result)
    atomic_json(root / 'resume_provenance.json', {'archive': str(archive.path),
        'status': 'prepared', 'time_utc': datetime.now(timezone.utc).isoformat()})
    return result
