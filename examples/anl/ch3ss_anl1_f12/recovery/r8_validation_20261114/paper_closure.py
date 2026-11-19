"""Four bounded B2PLYP scans for unresolved inner pathways; existing dispatcher.

The frozen inputs retain source hashes and reusable separated-fragment energies.
Scan observations require path/geometry review before any saddle or rate-model
decision. This owner never retries the old CC panel or publishes MESS inputs.
"""
from contextlib import ExitStack
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import re
import time

import numpy as np

from kinbot.anl import dispatch
from kinbot.anl.tasks import resources
from kinbot.driver_status import run_recorded
from kinbot.execution import atomic_json

HERE = Path(__file__).resolve().parent
INPUTS = HERE/'paper_closure_inputs.json'
ROUTE = ('B2PLYP/cc-pVTZ EmpiricalDispersion=GD3BJ NoSymm '
         'Integral=UltraFine SCF=(XQC,Conver=8)')
TERMINAL = {'complete', 'failed', 'blocked', 'skipped', 'unavailable'}


def load_inputs():
    data = json.loads(INPUTS.read_text())
    if data['schema'] != 1 or len(data['scans']) != 4:
        raise ValueError('Unexpected closure plan.')
    for case in data['scans']:
        scan, mol = case['scan'], case['molecule']
        xyz = np.asarray(mol['positions'], float)
        first, second = scan['bond']
        if (xyz.shape != (len(mol['symbols']), 3) or not np.isfinite(xyz).all()
                or mol['charge'] != 0 or mol['multiplicity'] != 2
                or not 0 <= first < len(xyz) or not 0 <= second < len(xyz)
                or first == second or not 1 <= scan['steps'] <= 40
                or not 0 < abs(scan['step_angstrom']) <= .1):
            raise ValueError('Invalid frozen scan definition.')
        if abs(np.linalg.norm(xyz[first]-xyz[second])-scan['initial_distance_angstrom']) > 1.e-8:
            raise ValueError('Scan atom mapping changed.')
        if scan['initial_distance_angstrom']+scan['steps']*scan['step_angstrom'] <= 0:
            raise ValueError('Scan reaches an invalid bond distance.')
    return data


def graphs():
    result = {}
    for case in load_inputs()['scans']:
        scan = case['scan']
        first, second = [i+1 for i in scan['bond']]
        body = ('%mem={{WORK_MEMORY_MB}}MB\n%nprocshared={{CORES}}\n%chk=scan.chk\n'
            '#P '+ROUTE+' Opt=(ModRedundant,CalcFC,Tight,MaxCycle=150)\n\n'
            'Bounded inner-path scan\n\n{{CHARGE}} {{MULT}}\n{{CARTESIAN}}\n\n'
            f'B {first} {second} S {scan["steps"]} {scan["step_angstrom"]:.12f}\n\n')
        task = dict(id='scan', kind='external', backend='gaussian', geometry_from='initial',
            resources=resources('partition_max', max_cores=8, partition='month-long-cpu'),
            input_name='scan.com', input_template=body, command=['g16'], stdin='scan.com',
            stdout='scan.log', stderr='scan.err', required_outputs=['scan.log'],
            success_marker=dict(file='scan.log', contains='Normal termination of Gaussian'),
            failure_markers=[dict(file='scan.log', contains=x) for x in
                             ('Error termination', 'Optimization stopped.')])
        result[case['name']] = dict(schema=1, name=case['name'], molecule=case['molecule'],
            limits=dict(max_nodes=20, max_cores_per_node=8), tasks=[task],
            intent=dict(question=case['question'], scan=scan,
                input_template_sha256=hashlib.sha256(body.encode()).hexdigest(),
                seed_evidence=case['seed_evidence'], scientific_approval=False,
                integration='Review the constrained path before any saddle/kinetic use.'))
    return result


def _prepare(root):
    root = Path(root).resolve()
    source = load_inputs()
    frozen = root/'plan.json'
    if frozen.exists() and json.loads(frozen.read_text()) != source:
        raise ValueError('Existing closure inputs differ; use a distinct campaign directory.')
    atomic_json(frozen, source)
    specs = graphs()
    for name, spec in specs.items():
        folder = root/name
        if folder.exists():
            if not (folder/'state.json').is_file():
                # A partial stage is preserved for inspection, not deleted or
                # silently replaced while another actor may still be staging.
                raise RuntimeError(f'{folder}: incomplete staging; inspect before continuing.')
            _, saved, state = dispatch._load(folder)
            if saved['molecule'] != spec['molecule'] or saved.get('intent') != spec['intent']:
                raise ValueError(f'{folder}: frozen scientific inputs differ.')
            # Dispatcher protects resolved resource/input hashes, including any
            # documented bounded recovery, on repeated entry.
            for task in saved['tasks']:
                if task['id'] in state['tasks']:
                    dispatch._verify_stage_files(folder, task, state['tasks'][task['id']])
        else:
            dispatch.prepare(spec, folder)
    return specs


def execute(root, *, submit=False, interval=20):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        lock = stack.enter_context((root/'controller.lock').open('a'))
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print('Existing closure owner retains control.', flush=True)
            return 0
        specs = _prepare(root)
        if not submit:
            print(f'Staged {len(specs)} scans; no submission.', flush=True)
            return 0
        for name in specs:
            lock = stack.enter_context((root/name/'drive.lock').open('a'))
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

        def drive():
            while True:
                states, infrastructure = {}, {}
                for name in specs:
                    try:
                        state = dispatch.advance(root/name, submit=True)
                    except dispatch.SchedulerUnavailableError as error:
                        infrastructure[name] = str(error)
                        continue
                    states[name] = dispatch._summary(specs[name], state)
                done = not infrastructure and all(
                    value in TERMINAL for state in states.values() for value in state.values())
                atomic_json(root/'progress.json', dict(workflows=states,
                    infrastructure=infrastructure, execution_finished=done,
                    scientific_approval=False, mess_updated=False))
                if done:
                    collect(root)
                    return int(any(value != 'complete' for s in states.values() for value in s.values()))
                time.sleep(interval)
        return run_recorded(drive, root/'controller.json')


def scan_points(path, case):
    """Read converged native scan geometries with the correlated DH energy."""
    raw = Path(path).read_bytes()
    lines = raw.decode(errors='replace').splitlines()
    xyz = None
    numbers = None
    energy = None
    points = []
    for i, line in enumerate(lines):
        if ('Input orientation:' in line or 'Standard orientation:' in line):
            rows = []
            for row in lines[i+5:]:
                fields = row.split()
                if len(fields) != 6 or not fields[0].isdigit():
                    break
                rows.append(fields)
            if rows:
                numbers = [int(row[1]) for row in rows]
                xyz = np.array([row[3:6] for row in rows], float)
        match = re.search(r'(?<![A-Za-z0-9])E\(B2PLYP\)\s*=\s*([-+\d.DEde]+)', line)
        if match:
            energy = float(match[1].replace('D', 'E').replace('d', 'e'))
        if 'Optimization completed.' in line:
            from ase.data import atomic_numbers
            expected = [atomic_numbers[s] for s in case['molecule']['symbols']]
            if numbers != expected or xyz is None or energy is None:
                raise ValueError('Converged scan point lacks matching atoms or correlated energy.')
            first, second = case['scan']['bond']
            points.append(dict(line=i+1, energy_hartree=energy, positions=xyz.tolist(),
                distance_angstrom=float(np.linalg.norm(xyz[first]-xyz[second]))))
            energy = None
    scan = case['scan']
    targets = [scan['initial_distance_angstrom']+i*scan['step_angstrom'] for i in range(len(points))]
    coordinate_ok = all(abs(p['distance_angstrom']-target) < .003 for p,target in zip(points,targets))
    complete = (len(points) == scan['steps']+1 and coordinate_ok
                and any('Normal termination of Gaussian' in line for line in lines)
                and not any('Error termination' in line or 'Optimization stopped.' in line for line in lines))
    peaks = [i for i in range(1,len(points)-1)
             if points[i]['energy_hartree'] > max(points[i-1]['energy_hartree'],points[i+1]['energy_hartree'])]
    return dict(source=str(path), sha256=hashlib.sha256(raw).hexdigest(), points=points,
        complete_constrained_scan=complete, coordinate_targets_met=coordinate_ok,
        candidate_maximum_indices=peaks, scientific_approval=False,
        next_decision='Check H attachments, branch continuity, spin and correlated energies. '
          'A continuous inner maximum can seed an unconstrained saddle/frequency and two-sided path; '
          'a branch jump or failed scan is not a barrierless classification.')


def collect(root):
    root = Path(root).resolve()
    result = dict(scientific_approval=False, mess_updated=False, scans={})
    for case in load_inputs()['scans']:
        path = root/case['name']/'tasks/scan/scan.log'
        try:
            item = scan_points(path, case)
        except (OSError, ValueError) as error:
            item = dict(complete_constrained_scan=False, error=str(error))
        result['scans'][case['name']] = item
    atomic_json(root/'scan_observations.json', result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('plan','stage','run','collect'))
    parser.add_argument('root', type=Path)
    args = parser.parse_args(argv)
    if args.action == 'plan':
        print(json.dumps({c['name']:dict(question=c['question'],scan=c['scan'])
                          for c in load_inputs()['scans']},indent=2))
    elif args.action == 'collect':
        result = collect(args.root)
        print(json.dumps({k:dict(points=len(v.get('points',[])),
            complete_constrained_scan=v['complete_constrained_scan'])
            for k,v in result['scans'].items()},indent=2))
    else:
        return execute(args.root, submit=args.action == 'run')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
