"""Read-only phase/workload evidence; never infer an ETA from file mtimes.

Uses only the standard library, one squeue query, current manifests and bounded
log excerpts. Does not import a dispatcher, open a QC database, hash checkpoints,
preflight, submit, repair, or execute native QC. Reports are written outside the
run directory. Counts are recorded observations, not scientific validation.
"""
import argparse
import ast
from collections import Counter, defaultdict
from datetime import datetime, timezone
import getpass
import json
from pathlib import Path
import re
import subprocess
import time


LIMIT = 256 * 1024
KEYS = ('vrc_tst_scan_points', 'vrc_tst_asymptote_distance',
        'vrc_tst_scan_gradient', 'vrc_tst_scan_fmax',
        'rotdpy_flux_parameters', 'pes', 'simultaneous_kinbot')


class Reader:
    def __init__(self):
        self.errors = []

    def text(self, path, *, head=False):
        try:
            with path.open('rb') as stream:
                if not head:
                    stream.seek(max(0, path.stat().st_size - LIMIT))
                return stream.read(LIMIT).decode(errors='replace')
        except FileNotFoundError:
            return ''
        except OSError as error:
            self.errors.append(f'{path}: {error}')
            return ''

    def json(self, path):
        try:
            if path.stat().st_size > 32 * 1024**2:
                raise ValueError('JSON larger than 32 MiB')
            return json.loads(path.read_text())
        except FileNotFoundError:
            return {}
        except (OSError, ValueError) as error:
            self.errors.append(f'{path}: {error}')
            return {}


def queue(root):
    try:
        result = subprocess.run(
            ['squeue', '-h', '-u', getpass.getuser(), '-o', '%i|%T|%M|%Z'],
            capture_output=True, text=True, check=False, timeout=20)
        if result.returncode:
            raise RuntimeError(result.stderr.strip())
        jobs = []
        for line in result.stdout.splitlines():
            if not line.strip():
                continue
            fields = line.split('|', 3)
            if len(fields) != 4:
                raise ValueError('Unrecognized squeue row: ' + line)
            ident, state, elapsed, work = (v.strip() for v in fields)
            if not work or work in ('(null)', 'N/A'):
                raise ValueError('Slurm omitted a job working directory')
            directory = Path(work).resolve()
            if directory.is_relative_to(root):
                jobs.append(dict(job_id=ident, state=state,
                                 elapsed=elapsed, directory=str(directory)))
        return jobs, None
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        return [], str(error)


def parameters(reader, case, worker):
    # The worker log prints resolved defaults, which a minimal input omits.
    result = {}
    for line in reader.text(worker/'kinbot.log', head=True).splitlines():
        for key in KEYS:
            match = re.match(r'\s*' + key + r':\s*(.*)$', line)
            if match:
                try:
                    result[key] = ast.literal_eval(match[1])
                except (ValueError, SyntaxError):
                    result[key] = match[1]
    result.update(reader.json(case/'input.json'))
    result.update(reader.json(worker/(worker.name + '.json')))
    return result


def surface_plan(params):
    plans = []
    model = params.get('vrc_tst_surface_model') or {}
    for channel in model.get('channels', []):
        for stage in channel.get('bottlenecks', []):
            controls = {**params, **stage.get('parameters', {})}
            grid = controls.get('vrc_tst_scan_points')
            # Never silently substitute today's defaults for a saved run.
            numeric = (isinstance(grid, list) and grid and all(
                isinstance(x, (int, float)) for x in grid))
            for state in stage.get('states', []):
                plans.append(dict(bond=channel.get('bond'), stage=stage.get('name'),
                    multiplicity=state.get('multiplicity'), root=state.get('root'),
                    relaxed_points=len(grid) if numeric else None,
                    scan_grid=grid, asymptote=controls.get('vrc_tst_asymptote_distance'),
                    gradient=controls.get('vrc_tst_scan_gradient'),
                    sampling_surfaces=len(stage.get('distances_angstrom', []))))
    return plans


def family(ident):
    if ident in ('l3_geometry', 'l3_geometry_tz', 'gaussian_vpt2'):
        return ident
    if ident.startswith('harmonic'):
        return 'harmonic'
    if ident in ('anl1_5z', 'anl1_6z'):
        return '5Z/6Z energies'
    if ident.startswith('ccsdt'):
        return 'higher-order CC'
    return 'other components'


def graph(reader, path, jobs):
    spec, state = reader.json(path/'workflow.json'), reader.json(path/'state.json')
    records = []
    for task in spec.get('tasks', []):
        ident = task['id']
        entry = state.get('tasks', {}).get(ident, {})
        work = path/'tasks'/ident
        row = dict(task=ident, status=entry.get('status', 'not_materialized'),
                   geometry_from=task.get('geometry_from'),
                   depends_on=task.get('depends_on', []),
                   job_id=entry.get('job_id'), error=entry.get('error'),
                   resources=task.get('resources'),
                   optimizer=task.get('optimizer'),
                   live=jobs.get(str(work.resolve())))
        if task.get('kind') == 'ase_optimize':
            lines = reader.text(work/'optimization.log').splitlines()
            steps = [s for s in lines if re.match(r'\s*(Sella|BFGS)\s+\d+\s', s)]
            row['optimizer_rows'] = steps[-6:]
            row['optimizer_first_observed_row'] = steps[0] if steps else None
            if steps:
                parts = steps[-1].split()
                try:
                    row.update(last_step=int(parts[1]), last_fmax=float(parts[4]))
                except (ValueError, IndexError):
                    pass
        # Read one current native output for wall-time text, never every force
        # evaluation or anything in retry/archive directories.
        candidates = [p for p in work.glob('*') if p.suffix in ('.out', '.log')
                      and p.name != 'optimization.log']
        if candidates:
            try:
                native = max(candidates, key=lambda p: p.stat().st_mtime)
                lines = reader.text(native).splitlines()
                row['native'] = dict(file=native.name,
                    age_seconds=round(time.time()-native.stat().st_mtime),
                    tail=lines[-8:], timing_lines=[s for s in lines if re.search(
                        r'Elapsed time:|REAL TIME|wall.*time|time.*wall', s, re.I)][-8:])
            except OSError as error:
                reader.errors.append(f'{work}: {error}')
        records.append(row)
    return dict(directory=str(path), molecule=spec.get('molecule'),
                intent=spec.get('intent'), tasks=records)


def collect(root, *, scheduler=None):
    root = Path(root).expanduser().resolve(strict=True)
    reader = Reader()
    live, scheduler_error = queue(root) if scheduler is None else scheduler
    jobs = {j['directory']: j for j in live}
    cases = []
    inputs = ([root/'input.json'] if (root/'input.json').exists()
              else sorted(root.glob('*/input.json')))
    for input_path in inputs:
        case = input_path.parent
        workers = []
        for directory in sorted(case.iterdir()):
            if not directory.is_dir() or not (directory/'kinbot.log').is_file():
                continue
            params = parameters(reader, case, directory)
            groups = {}
            for name, pattern in (
                    ('anl', 'anl_automatic/states/*/composite/workflow.json'),
                    ('vrc_scans', 'vrctst/molpro_scan/*/workflow.json')):
                groups[name] = [graph(reader, p.parent, jobs)
                                for p in sorted(directory.glob(pattern))]
            # Use the current correction dispatcher pointer, avoiding obsolete
            # dispatch directories from older input fingerprints.
            pointer = reader.json(directory/'vrctst/molpro/current_dispatch.json')
            correction_graph = None
            if pointer.get('run_dir'):
                current = Path(pointer['run_dir']).resolve()
                if current.is_relative_to(root):
                    correction_graph = graph(reader, current, jobs)
                else:
                    reader.errors.append(f'Correction graph outside run: {current}')
            workers.append(dict(worker=directory.name, **groups,
                surface_plan=surface_plan(params), controls={k: params.get(k) for k in KEYS},
                branches=reader.json(directory/'.kinbot_post_reaction_branches.json'),
                correction_graph=correction_graph,
                correction_files=[str(p.relative_to(directory))
                                  for p in sorted(directory.glob('vrctst/corr_*.json'))],
                network_request=reader.json(directory/'anl_automatic/network_request.json'),
                handoff_present=(directory/'anl_automatic/handoff.json').is_file()))
        cases.append(dict(case=case.name, workers=workers,
            pes_summary=reader.text(case/'pes_summary.txt'),
            rotd_inputs=[str(p.relative_to(case)) for p in sorted(case.glob('rotdPy/*.py'))],
            rotd_results=[str(p.relative_to(case)) for p in sorted(case.glob('rotdPy/*.rotdpy.json'))]))
    return dict(schema=1, root=str(root), captured_utc=datetime.now(timezone.utc).isoformat(),
                scheduler_error=scheduler_error, live_jobs=live, cases=cases,
                read_errors=reader.errors,
                note='Recorded status only. No ETA from mtimes, counts, or unconverged forces; future PES workers may add work.')


def render(report):
    out = ['Read-only phase evidence | ' + report['captured_utc'],
           f'Scheduler: {report["scheduler_error"] or "responding"}; '
           f'allocations={len(report["live_jobs"])}; read errors={len(report["read_errors"])}']
    for case in report['cases']:
        out.append('\n' + case['case'])
        out.append(f'  ANL graphs={sum(len(w["anl"]) for w in case["workers"])}; '
                   f'handoff files={sum(w["handoff_present"] for w in case["workers"])}')
        groups = defaultdict(Counter)
        for worker in case['workers']:
            for g in worker['anl']:
                for t in g['tasks']:
                    groups[family(t['task'])][t['status']] += 1
        for name, counts in sorted(groups.items()):
            out.append('  ' + name + ': ' + ', '.join(f'{k}={v}' for k, v in sorted(counts.items())))
        for worker in case['workers']:
            plan = worker['surface_plan']
            planned = (sum(s['relaxed_points'] for s in plan)
                       if plan and all(s['relaxed_points'] is not None for s in plan) else 'unknown')
            scans = [t for g in worker['vrc_scans'] for t in g['tasks']]
            counts = dict(Counter(t['status'] for t in scans))
            out.append(f'  {worker["worker"]} VRC: {len(plan)} configured leaves; '
                       f'planned relaxed points={planned}; current scan records={counts}')
            for g in worker['vrc_scans']:
                for t in g['tasks']:
                    if t['status'] in ('submitted', 'submitting', 'failed'):
                        opt = t.get('optimizer') or {}
                        out.append(f'    {Path(g["directory"]).name}: step={t.get("last_step")}, '
                                   f'fmax={t.get("last_fmax")}, target={opt.get("fmax")}; '
                                   f'job={(t.get("live") or {}).get("elapsed", "not in queue")}')
            correction = worker['correction_graph']
            out.append('  Correction task records: ' + str(dict(Counter(
                t['status'] for t in correction['tasks'])) if correction else 'not materialized')
                + f'; correction JSON files={len(worker["correction_files"])}')
        out.append('  PES: ' + ' | '.join(case['pes_summary'].splitlines()[:4]))
        out.append(f'  PES ROTD inputs={len(case["rotd_inputs"])}; results={len(case["rotd_results"])}')
    out.append('\nFull JSON includes individual graphs, dependencies and optimizer tails. '
               'Counts do not represent percent of wall time.')
    return '\n'.join(out)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    root, output = args.root.resolve(strict=True), args.output.resolve()
    if output.is_relative_to(root):
        parser.error('--output must be outside the run tree')
    if output.exists():
        parser.error('--output already exists; use a fresh report path')
    report = collect(root)
    with output.open('x') as stream:
        json.dump(report, stream, indent=2)
    print(render(report))
    print('Full report:', output)


if __name__ == '__main__':
    main()
