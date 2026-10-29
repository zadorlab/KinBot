"""Read-only evidence capture for QC acceptance and rejected-search review.

Database snapshots and reports are written only outside the source run. No QC
object, recovery driver, scheduler command or saved Python worker is executed.
Findings are consistency questions, never proof of a bug or a useful retry.
"""
import argparse
from collections import Counter, defaultdict
from contextlib import closing
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3
import tarfile
import time

from ase.db.sqlite import SQLite3Database
from kinbot.frequencies import frequency_validation_error
from kinbot.run_pes_audit import rejected_searches

METADATA = {
    'input.json', 'resolved_parameters.json', 'pes_status.json', 'pes_summary.txt',
    'pes_exploration.json', 'kinbot_status.json', 'kinbot_monitor.out',
    'kinbot_optimization_monitor.json', 'accepted_l2_network.json',
    '.kinbot_search_state.json', '.kinbot_theory_jobs.json', '.kinbot_slurm_jobs.json',
    '.kinbot_revisit.json', '.kinbot_search_options.json', 'l1_screening.json',
    'l1_endpoint_checks.json', '.kinbot_resource_policy.json', '.kinbot_live_control.json',
    'launch.json', 'prepared.json', 'resume_provenance.json',
    '.kinbot_policy_admissions.json', '.kinbot_slurm_nodes.json',
}
RECEIPTS = ('.exitcode', '.result.json', '.scf_recovery.json',
            '.optimization_budget.json', '.optimization_resume.json')
NATIVE = re.compile(r'normal termination|error termination|SCF Done:|SCFChk:|'
    r'convergence criterion not met|convergence failure|no lower point found|'
    r'Frequencies\s+--|Zero-point correction|S\*\*2|maximum force|rms force|'
    r'Worker failure:|frequency.*fail|did not converge|^\s*Sella\s+\d|'
    r'Traceback|\b\w*(?:Error|Exception):|QPErr|No space left|Permission denied', re.I)


def pack(value):
    if hasattr(value, 'tolist'):
        return value.tolist()
    raise TypeError(type(value).__name__)


class ReadOnlyDatabase(SQLite3Database):
    def _connect(self):
        return sqlite3.connect(Path(self.filename).resolve().as_uri()+'?mode=ro', uri=True)


def database_snapshot(source, target):
    """SQLite backup API includes committed WAL content in a consistent copy."""
    started = time.monotonic()
    def progress(*_):
        if time.monotonic()-started > 30:
            raise TimeoutError('Database snapshot exceeded 30 seconds; no source changes made.')
    target.parent.mkdir(parents=True, exist_ok=True)
    with closing(sqlite3.connect(source.as_uri()+'?mode=ro', uri=True, timeout=3)) as src:
        src.execute('PRAGMA query_only=ON')
        with closing(sqlite3.connect(target)) as dst:
            src.backup(dst, pages=256, progress=progress, sleep=.05)
    history = defaultdict(list)
    for row in ReadOnlyDatabase(str(target), use_lock_file=False).select(sort='-id'):
        name = row.get('name')
        if not name:
            continue
        data = row.data
        history[name].append(dict(id=row.id, mtime_ase=row.mtime,
            status=data.get('status'), energy_eV=data.get('energy'),
            zpe_hartree=data.get('zpe'), frequencies=data.get('frequencies'),
            error=data.get('error'), natoms=row.natoms,
            positions=row.positions, numbers=row.numbers,
            calculator=row.get('calculator'), parameters=row.get('calculator_parameters')))
    return dict(history)


def accepted_checks(worker, parameters, rows, history):
    """Replay current data gates; cache-to-accepted-result identity is unverified."""
    targets = {(worker+'_well'+('_high' if parameters.get('high_level', 1) else ''), 0)} if rows else set()
    for row in rows:
        if row['code'] != -1:
            continue
        suffix = '_high' if parameters.get('high_level', 1) else ''
        if '_hom_sci_' not in row['name']:
            targets.add((row['name']+suffix, 1))
        targets.update((p+'_well'+suffix, 0) for p in row['products'])
    checks = []
    for job, order in sorted(targets):
        versions = history.get(job, [])
        problems = []
        if not versions:
            problems.append('expected cache label absent; selected result/reuse mapping needs inspection')
        else:
            latest = versions[0]
            if latest['status'] != 'normal':
                problems.append('latest expected cache row is not normal; selected-result provenance needs inspection')
            for key in ('energy_eV', 'zpe_hartree'):
                value = latest[key]
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                    problems.append('missing/nonfinite '+key)
            import numpy as np
            positions = np.asarray(latest['positions'])
            if latest['natoms'] < 1 or positions.shape != (latest['natoms'], 3) or not np.isfinite(positions).all():
                problems.append('malformed/nonfinite cached geometry')
            if latest['natoms'] > 1:
                try:
                    error = frequency_validation_error(latest['frequencies'], order,
                                                      parameters.get('imagfreq_threshold', 50.))
                except (TypeError, ValueError) as exc:
                    error = 'cannot validate frequency data: '+str(exc)
                if error:
                    problems.append(error)
        checks.append(dict(job=job, order=order, latest_row_id=versions[0]['id'] if versions else None,
            questions=problems, status='inspect_evidence' if problems else 'cache_data_gates_pass',
            selected_result_provenance_verified=False, native_acceptance_verified=False))
    return checks


def capture(root, output):
    root, output = Path(root).resolve(strict=True), Path(output).resolve()
    if output == root or output.is_relative_to(root):
        raise ValueError('Output must be outside the source run.')
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError('Output directory must be empty.')
    print('OUTPUT_DIRECTORY', output, flush=True)
    cases = [root] if (root/'input.json').is_file() else sorted(
        p for p in root.iterdir() if p.is_dir() and not p.is_symlink()
        and not p.name.startswith('.') and (p/'input.json').is_file())
    if not cases:
        raise ValueError('No PES input directories found.')
    issues, copied, results = [], [], []
    report = dict(schema=1, root=str(root), started_utc=datetime.now(timezone.utc).isoformat(),
        issues=issues, cases=results, source_modified=False, submits_jobs=False,
        global_snapshot=False, retry_success_inferred=False)
    report['rules_at_audit_time'] = {}
    for name in ('qc.py', 'irc.py', 'reaction_generator.py', 'frequencies.py',
                 'l1_screening.py', 'l1_endpoint.py', 'optimization_budget.py',
                 'minimum_recovery.py', 'ase_modules/calculators/gaussian.py',
                 'tpl/ase_sella_irc.tpl.py', 'tpl/ase_sella_opt_well.tpl.py'):
        code = Path(__file__).parent/name
        data = code.read_bytes()
        target = output/'audit_code'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        report['rules_at_audit_time'][name] = hashlib.sha256(data).hexdigest()
    report['historical_loaded_code_verified'] = False
    def read(path, limit=16*1024*1024, tail=False):
        try:
            if path.is_symlink() or not path.resolve().is_relative_to(root):
                raise ValueError('Symlink/outside-run source excluded')
            path = path.resolve(strict=True)
            before = path.stat()
            with path.open('rb') as stream:
                offset = max(0, before.st_size-limit) if tail else 0
                stream.seek(offset)
                data = stream.read(limit)
            after = path.stat()
            entry = dict(path=str(path.relative_to(root)), bytes=before.st_size,
                captured_bytes=len(data), byte_offset=offset, truncated=before.st_size>len(data),
                changed=(before.st_mtime_ns, before.st_size)!=(after.st_mtime_ns, after.st_size),
                sha256_captured=hashlib.sha256(data).hexdigest())
            if entry['changed']:
                issues.append(dict(path=entry['path'], issue='changed during read'))
            destination = output/'files'/entry['path']
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
            copied.append(entry)
            return data.decode(errors='replace')
        except (OSError, ValueError) as error:
            issues.append(dict(path=str(path), issue=str(error)))
            return ''
    def metadata(scope):
        found = {}
        paths = [p for p in scope.iterdir() if p.name in METADATA or
                 p.name == scope.name+'.json' or p.name.endswith(RECEIPTS)]
        for directory in ('.kinbot_timeout_allocations', '.kinbot_policy_submissions', '.kinbot_policy_history'):
            paths.extend((scope/directory).glob('*.json'))
        for path in sorted(paths):
            if path.is_file():
                value = read(path)
                try:
                    found[str(path.relative_to(scope))] = json.loads(value) if path.suffix == '.json' else value
                except ValueError:
                    issues.append(dict(path=str(path), issue='JSON could not be decoded'))
        return found
    report['root_metadata'] = metadata(root)
    for case in cases:
        item = dict(case=case.name, metadata=metadata(case) if case != root else report['root_metadata'], workers=[])
        results.append(item)
        item['mess_inputs'] = {p.name: read(p) for p in sorted((case/'me').glob('*.inp'))}
        workers = sorted(p for p in case.iterdir() if p.is_dir() and not p.is_symlink()
                         and ((p/'kinbot.db').is_file() or (p/(p.name+'.json')).is_file()))
        for worker in workers:
            print('READING', case.name, worker.name, flush=True)
            entry = dict(worker=worker.name, metadata=metadata(worker), rows=[], database={},
                         accepted_checks=[], rejected=[], native=[])
            item['workers'].append(entry)
            monitor = entry['metadata'].get('kinbot_monitor.out', '')
            for line in monitor.splitlines():
                fields = line.split()
                if len(fields) < 3:
                    continue
                try:
                    entry['rows'].append(dict(code=int(fields[0]), step=fields[1], name=fields[2], products=fields[3:]))
                except ValueError:
                    issues.append(dict(path=str(worker/'kinbot_monitor.out'), issue='malformed row: '+line))
            try:
                if (worker/'kinbot.db').is_file():
                    if (worker/'kinbot.db').is_symlink() or not (worker/'kinbot.db').resolve().is_relative_to(root):
                        raise ValueError('Symlink/outside-run database excluded')
                    db = output/'databases'/case.name/worker.name/'kinbot.db'
                    history = database_snapshot(worker/'kinbot.db', db)
                    entry['database'] = {name: dict(total_rows=len(rows), latest_rows=rows[:3])
                                         for name, rows in history.items()}
                    parameters = dict(item['metadata'].get('input.json', {}))
                    parameters.update(entry['metadata'].get(worker.name+'.json', {}))
                    entry['accepted_checks'] = accepted_checks(worker.name, parameters, entry['rows'], history)
            except Exception as error:
                issues.append(dict(path=str(worker/'kinbot.db'), issue=f'database audit: {error}'))
            entry['rejected'] = rejected_searches(worker, monitor, issues)
            # Fixed, shallow directories only; no traversal of scratch/ANL/sampling payloads.
            selected = set()
            for directory in (worker, worker/'conf', worker/'perm'):
                if not directory.is_dir() or directory.is_symlink():
                    continue
                if directory != worker:
                    entry['metadata'].update({str(Path(directory.name)/k): v for k,v in metadata(directory).items()})
                selected.update(p for p in directory.iterdir() if p.is_file() and
                    (p.suffix in {'.com', '.py', '.log', '.err', '.stdout', '.stderr'}))
            for pattern in ('*_minimum_recovery/attempt_*/*', 'conf/*_minimum_recovery/attempt_*/*',
                            '*_minimum_recovery/attempt_*/*/*', 'conf/*_minimum_recovery/attempt_*/*/*'):
                for path in sorted(worker.glob(pattern)):
                    if path.suffix in {'.com', '.log', '.py'}:
                        selected.add(path)
                    elif path.suffix == '.json':
                        value = read(path)
                        try:
                            entry['metadata'][str(path.relative_to(worker))] = json.loads(value)
                        except ValueError:
                            issues.append(dict(path=str(path), issue='recovery JSON could not be decoded'))
            # Preserve original inputs/outputs identified by SCF recovery receipts.
            for name, value in list(entry['metadata'].items()):
                if name.endswith('.scf_recovery.json') and isinstance(value, dict):
                    for file in value.get('preserved_files', {}).values():
                        path = Path(file.get('archive', ''))
                        path = path if path.is_absolute() else worker/path
                        if path.suffix in {'.com', '.log'}:
                            selected.add(path)
            for path in sorted(selected):
                native = path.suffix not in {'.com', '.py'}
                text = read(path, limit=65536 if native else 1024*1024, tail=native)
                lines = [s.strip() for s in text.splitlines() if NATIVE.search(s)]
                if lines:
                    entry['native'].append(dict(path=str(path.relative_to(root)), evidence=lines[-35:]))
            try:
                if (worker/'kinbot_monitor.out').is_file() and (worker/'kinbot_monitor.out').read_text() != monitor:
                    issues.append(dict(path=str(worker/'kinbot_monitor.out'), issue='search outcomes changed during worker capture'))
            except OSError as error:
                issues.append(dict(path=str(worker/'kinbot_monitor.out'), issue=str(error)))
    report['finished_utc'] = datetime.now(timezone.utc).isoformat()
    report['files'] = copied
    report['scope'] = ('Cache-gate replay plus exported evidence; no certification of scientific acceptance, '
        'pathway absence, or benefit from another retry. Native log captures are bounded tails; full DB histories are included.')
    (output/'report.json').write_text(json.dumps(report, indent=2, default=pack)+'\n')
    lines = [report['scope']]
    for case in results:
        outcomes = Counter(row['code'] for worker in case['workers'] for row in worker['rows'])
        checks = [c for w in case['workers'] for c in w['accepted_checks']]
        questions = [c for c in checks if c['questions']]
        lines.append(f"CASE {case['case']} search_codes={dict(outcomes)} cache_checks={len(checks)} questions={len(questions)}")
        for worker in case['workers']:
            counts = Counter(r['category'] for r in worker['rejected'])
            lines.append(f"WORKER {worker['worker']} rejected_search_evidence={dict(counts)}")
            for check in worker['accepted_checks']:
                if check['questions']:
                    lines.append('CONSISTENCY_QUESTION '+json.dumps(check))
    lines.append(f'CAPTURE_ISSUES {len(issues)}; truncated_file_captures={sum(p["truncated"] for p in copied)}')
    lines.extend(json.dumps(issue) for issue in issues)
    (output/'summary.txt').write_text('\n'.join(lines)+'\n')
    with (output/'evidence.txt').open('w') as stream:
        stream.write('\n'.join(lines)+'\n\nFULL STRUCTURED EVIDENCE\n')
        json.dump(report, stream, indent=2, default=pack)
    archive = output/'qc-decision-evidence.tar.gz'
    with tarfile.open(archive, 'w:gz') as tar:
        for path in sorted(output.iterdir()):
            if path != archive:
                tar.add(path, arcname=path.name)
    from kinbot.execution import file_sha256
    return dict(output=str(output), summary=str(output/'summary.txt'), text=str(output/'evidence.txt'),
        archive=str(archive), issues=len(issues), sha256=file_sha256(archive))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = capture(args.root, args.output)
    print(Path(result['summary']).read_text(), end='')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
