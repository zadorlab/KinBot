"""Bounded R8 validation; observations only, never automatic MESS promotion.

Run the remaining L2 IRC panel through saddle_recheck. Independent native
Gaussian paths and targeted Molpro benchmarks use the existing QC dispatcher.
"""
from __future__ import annotations
import argparse
from collections import Counter
from contextlib import ExitStack
import fcntl
import hashlib
import json
from pathlib import Path
import time
import uuid

from kinbot.anl import dispatch
from kinbot.anl.tasks import molpro_task, resources
from kinbot.execution import atomic_json

HERE = Path(__file__).resolve().parent
INDEPENDENT = (
    'paper_ts1', 'paper_ts2', 'paper_ts3', 'paper_ts4', 'paper_ts5', 'paper_ts6',
    'native_ch3ss_rxn_12_1', 'native_ch3ss_rxn_12_2', 'native_ch3ss_rxn_54_1',
    'native_ch3ss_rxn_2_2', 'native_ch3s_rxn_1_2', 'native_reverse_ts1',
    'ch3s_h_loss_scale_1p6', 'reverse_26b8db7c', 'reverse_b18dcbc8')
# Only these three have a specific paper-versus-L2 connection disagreement that
# presently justifies the much more expensive coupled-cluster reaction paths.
CC_PATHS = ('si_TS2', 'si_TS3', 'si_TS5')
CC_SADDLES = (*CC_PATHS, 'si_TS6', 'native_ch3ss_rxn_12_1')
CC_MINIMA = ('si_W1', 'si_W2', 'si_W3', 'si_W4')
L2_ROUTE = 'B2PLYP/cc-pVTZ EmpiricalDispersion=GD3BJ NoSymm Integral=UltraFine SCF=(XQC,Conver=8)'
GAUSS_HEADER = '%mem={{WORK_MEMORY_MB}}MB\n%nprocshared={{CORES}}\n'


def molecule(point):
    return {key: point[key] for key in ('symbols', 'positions', 'charge', 'multiplicity')}


def workflow(name, point, tasks, purpose):
    return dict(schema=1, name=name, molecule=molecule(point),
                limits=dict(max_nodes=20, max_cores_per_node=8), tasks=tasks,
                intent=dict(purpose=purpose, scientific_approval=False,
                            integration='none; MESS remains paused'))


def gaussian(ident, body, *, dependencies=(), parser=None, checkpoint=False):
    task = dict(id=ident, kind='external', backend='gaussian', geometry_from='initial',
        resources=resources('partition_max', max_cores=8, partition='month-long-cpu'),
        input_name=ident+'.com', input_template=body, command=['g16'], stdin=ident+'.com',
        stdout=ident+'.log', stderr=ident+'.err', required_outputs=[ident+'.log'],
        success_marker=dict(file=ident+'.log', contains='Normal termination of Gaussian'),
        failure_markers=[dict(file=ident+'.log', contains='Error termination')],
        depends_on=list(dependencies))
    if checkpoint:
        task['required_outputs'].append(ident+'.chk')
    if parser:
        task['result_parser'] = dict(file=ident+'.log', **parser)
    return task


def native_irc(name, point, graph_dir, *, step_sizes=(10, 5)):
    """One Hessian, then two directions at each of two integration steps."""
    hessian = gaussian('hessian', GAUSS_HEADER + '%chk=hessian.chk\n#P '+L2_ROUTE+
        ' Freq=Numer\n\nIndependent L2 Hessian\n\n{{CHARGE}} {{MULT}}\n{{CARTESIAN}}\n\n',
        checkpoint=True, parser=dict(kind='gaussian_harmonic', method='B2PLYP', basis='cc-pVTZ',
            dispersion='GD3BJ', natoms=len(point['symbols']), stationary_point='transition_state'))
    tasks = [hessian]
    # Absolute OLDCHK is necessary: the dispatcher executes in isolated scratch.
    # It is read-only; each branch writes a separate checkpoint. Dependency
    # acceptance verifies the parent's published artifact hashes before staging.
    checkpoint = graph_dir/'tasks/hessian/hessian.chk'
    if any(c in str(checkpoint) for c in ('\n', '\r', '%')):
        raise ValueError('Unsupported checkpoint pathname.')
    for size in step_sizes:
        points = 6000 // size
        for direction in ('Forward', 'Reverse'):
            ident = direction.lower()+'_'+str(size)
            body = GAUSS_HEADER + f'%oldchk={checkpoint}\n%chk={ident}.chk\n#P '+L2_ROUTE+(
                f' Geom=AllCheck Guess=Read IRC=(RCFC,GradientOnly,EulerPC,MassWeighted,'
                f'{direction},StepSize=-{size},MaxPoints={points},MaxCycle=100,Tight,Report=Cartesians)\n\n'
                '--Link1--\n') + GAUSS_HEADER + f'%chk={ident}.chk\n#P '+L2_ROUTE+(
                ' Geom=AllCheck Guess=Read Opt=(Tight,CalcFC,MaxCycle=300) Freq=Numer\n\n')
            tasks.append(gaussian(ident, body, dependencies=('hessian',), checkpoint=True))
    return workflow(name, point, tasks,
        'Independent mass-weighted IRC and half-step check; inspect the pre-optimization path termination and final minima.')


def cc_body(basis):
    # The article/SI does not specify a/b or the triples scaling convention.
    # This is an explicit independent F12b(T*) benchmark, not an exact replica.
    return f'basis={basis}\nrhf\nuccsd(t)-f12b,scale_trip=1\n'


def cc_task(ident, body, **kwargs):
    return molpro_task(ident, body, geometry_from=kwargs.pop('geometry_from', 'initial'),
        max_cores=8, partition='month-long-cpu', **kwargs)


def cc_refinement(name, point, saddle, path_check):
    body = cc_body('cc-pVDZ-F12')
    body += ("{optg,method=qsd,root=2,maxit=100,gradient=1.d-5,energy=1.d-8,"
             "savexyz='optimized.xyz'\nprint,history,gradient\n}\n" if saddle else
             "{optg,maxit=100,gradient=1.d-5,energy=1.d-8,savexyz='optimized.xyz'\nprint,history,gradient\n}\n")
    # Re-evaluate at the optimized geometry before the one frequency section.
    body += 'rhf\nuccsd(t)-f12b,scale_trip=1\nfrequencies,numerical\n'
    opt = cc_task('refine', body, parser=dict(kind='molpro_harmonic', basis='cc-pVDZ-F12',
        reference='ROHF', method='CCSD(T)-F12b',
        stationary_point='transition_state' if saddle else 'minimum'))
    opt.update(geometry_output='optimized.xyz')
    opt['required_outputs'].append('optimized.xyz')
    # Molpro may terminate normally after an optimization step limit. Never
    # pass such a geometry to a path or CBS task as a converged structure.
    opt['failure_markers'] = [dict(file='refine.out', contains=x) for x in
        ('No convergence', 'NO CONVERGENCE', 'Maximum number of iterations reached')]
    tasks = [opt]
    for basis, label in (('cc-pVTZ-F12', 'tz'), ('cc-pVQZ-F12', 'qz')):
        tasks.append(cc_task(label, cc_body(basis), geometry_from='refine',
            parser=dict(kind='molpro_energy', method='CCSD(T)-F12b', basis=basis, reference='ROHF')))
    if path_check:
        for idir, label in ((1, 'forward'), (-1, 'reverse')):
            # Native QSD paths supply an independent higher-level descent. They
            # are not assumed equivalent to a converged mass-weighted IRC.
            path = cc_task(label, cc_body('cc-pVDZ-F12') +
                f"{{optg,method=qsdpath,idir={idir},maxit=100,stepmax=0.1,"
                "savexyz='path_end.xyz'\nprint,history,gradient\n}\n",
                geometry_from='refine')
            path.update(geometry_output='path_end.xyz')
            path['required_outputs'].append('path_end.xyz')
            tasks.append(path)
            end = cc_task(label+'_minimum', cc_body('cc-pVDZ-F12')+
                "{optg,maxit=100,gradient=1.d-5,energy=1.d-8,savexyz='minimum.xyz'\nprint,history,gradient\n}\n"
                'rhf\nuccsd(t)-f12b,scale_trip=1\nfrequencies,numerical\n',
                geometry_from=label, parser=dict(kind='molpro_harmonic', basis='cc-pVDZ-F12',
                    reference='ROHF', method='CCSD(T)-F12b', stationary_point='minimum'))
            end.update(geometry_output='minimum.xyz')
            end['required_outputs'].append('minimum.xyz')
            tasks.append(end)
    return workflow(name, point, tasks,
        'Explicit ROHF-UCCSD(T)-F12b(T*) comparison. Review convergence, T1, modes, path and basis dependence; not a certified single-reference surface.')


def single_points(name, point):
    tasks = [cc_task(label, cc_body(basis), parser=dict(kind='molpro_energy',
        method='CCSD(T)-F12b', basis=basis, reference='ROHF' if point['multiplicity'] > 1 else 'RHF'))
        for basis, label in (('cc-pVTZ-F12', 'tz'), ('cc-pVQZ-F12', 'qz'))]
    return workflow(name, point, tasks, 'Matched fixed-geometry T/Q benchmark; a single point does not validate a saddle or connection.')


def h_loss_guess(point, distance):
    """Independent product-side C-H addition seeds, not the exhausted parent stretch."""
    import numpy as np
    p = molecule(point)
    xyz = np.array(p['positions'], float)
    carbon = p['symbols'].index('C')
    hydrogens = [i for i, symbol in enumerate(p['symbols']) if symbol == 'H']
    leaving = max(hydrogens, key=lambda i: float(np.linalg.norm(xyz[i]-xyz[carbon])))
    vector = xyz[leaving]-xyz[carbon]
    xyz[leaving] = xyz[carbon]+vector*distance/np.linalg.norm(vector)
    p['positions'] = xyz.tolist()
    return p


def h_loss(name, point, graph_dir):
    graph = native_irc(name, point, graph_dir)
    task = graph['tasks'][0]
    task['input_template'] = task['input_template'].replace(' Freq=Numer',
        ' Opt=(TS,CalcFC,Tight,MaxCycle=200) Freq=Numer')
    graph['intent']['purpose'] = 'Product-side H-addition saddle search; valid first-order saddle required before either IRC. An unsuccessful search is not a channel failure.'
    return graph


def screening(name, point, graph_dir):
    if point['stationary_point'] == 'transition_state':
        graph = native_irc(name, point, graph_dir, step_sizes=(10,))
        graph['tasks'][0]['input_template'] = graph['tasks'][0]['input_template'].replace(
            ' Freq=Numer', ' Opt=(TS,CalcFC,Tight,MaxCycle=200) Freq=Numer')
        graph['intent']['purpose'] = 'Direct L2 check of a recorded L1 energy exclusion, with actual L2 endpoints.'
        return graph
    atomic = len(point['symbols']) == 1
    route = ' SP' if atomic else ' Opt=(Tight,CalcFC,MaxCycle=300) Freq=Numer'
    task = gaussian('minimum', GAUSS_HEADER+'#P '+L2_ROUTE+route+
        '\n\nL2 product threshold check\n\n{{CHARGE}} {{MULT}}\n{{CARTESIAN}}\n\n',
        parser=None if atomic else dict(kind='gaussian_harmonic', method='B2PLYP',
            basis='cc-pVTZ', dispersion='GD3BJ', natoms=len(point['symbols'])))
    return workflow(name, point, [task], 'L2 fragment energy and minimum check; atomic ZPE is exactly zero.')


def graphs(root):
    evidence = json.loads((HERE/'evidence.json').read_text())
    points = evidence['points']; result = {}
    for name in INDEPENDENT:
        ident = 'irc_'+name
        result[ident] = native_irc(ident, points[name], root/ident)
    for name in (*CC_SADDLES, *CC_MINIMA):
        ident = 'cc_'+name
        result[ident] = cc_refinement(ident, points[name], name in CC_SADDLES, name in CC_PATHS)
    # Include both observed low-energy complexes, the higher isomerization and
    # genuine H-loss saddles. Benchmark the endpoints in the same convention.
    names = ('paper_ts1', 'paper_ts2', 'paper_ts3', 'paper_ts4', 'paper_ts5', 'paper_ts6',
             'native_ch3ss_rxn_12_1', 'native_ch3ss_rxn_12_2', 'native_ch3ss_rxn_54_1',
             'ch3s_h_loss_scale_1p6', 'native_reverse_ts1',
             'paper_ts1_F', 'paper_ts1_R', 'paper_ts3_R', 'paper_ts5_F',
             'paper_ts5_R', 'ch3s_h_loss_scale_1p6_F', 'ch3s_h_loss_scale_1p6_R',
             'native_ch3s_rxn_1_1_R')
    for name in names:
        ident = 'sp_'+name
        result[ident] = single_points(ident, points[name])
    # The 1.6 parent-start attempt ended at a product-side minimum. Reuse that
    # heavy-atom framework, never rerun its exhausted parent-side seed.
    for distance in (1.5, 1.8):
        ident = 'h_loss_product_'+str(distance).replace('.', 'p')
        result[ident] = h_loss(ident, h_loss_guess(points['ch3ss_h_loss_scale_1p6'], distance), root/ident)
    for name, point in evidence['screening_points'].items():
        result[name] = screening(name, point, root/name)
    return result


def fingerprint():
    return {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in ('campaign.py', 'evidence.json', 'remaining_l2.json')}


def _stage_unlocked(root):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    receipt = root/'stage.json'
    signature = fingerprint()
    if receipt.exists() and json.loads(receipt.read_text())['sources'] != signature:
        raise ValueError('Staged scientific definitions changed; use a new campaign directory.')
    specs = graphs(root)
    # Freeze before preparation; a transient failure can resume graph by graph.
    atomic_json(receipt, dict(schema=1, sources=signature, workflows=list(specs), scientific_approval=False))
    for name, spec in specs.items():
        directory = root/name
        preparing = root/(name+'.preparing.json')
        if directory.exists() and not (directory/'state.json').exists():
            # prepare() never submits. Only recover an interrupted preparation
            # owned by this exact campaign, and preserve it for inspection.
            if (not preparing.is_file() or json.loads(preparing.read_text()) != spec
                    or (directory/'requested.json').exists()
                    or any(p.is_file() for pattern in ('execution.json', '*.exitcode', '*.jobid')
                           for p in directory.rglob(pattern))):
                raise ValueError(f'{name}: unowned partial preparation; preserve and inspect it.')
            archive = root/'interrupted_preparations'
            archive.mkdir(exist_ok=True)
            directory.rename(archive/(name+'.'+uuid.uuid4().hex))
        if directory.exists():
            plan = directory/'requested.json'
            _, saved, state = dispatch._load(directory)  # verifies workflow/state checksum
            if plan.exists():
                if json.loads(plan.read_text()) != spec:
                    raise ValueError(f'{name}: scientific definition changed.')
            else:
                # Recover a stop after prepare published a complete graph but
                # before our receipt. Only site-resolved resources may differ.
                from copy import deepcopy
                normalized = deepcopy(saved)
                expected = {t['id']:t for t in spec['tasks']}
                for t in normalized['tasks']:
                    t['resources'] = expected[t['id']]['resources']
                if normalized != spec:
                    raise ValueError(f'{name}: incomplete or changed preparation; preserve and inspect it.')
                atomic_json(plan, spec)
        else:
            atomic_json(preparing, spec)
            dispatch.prepare(spec, directory)
            atomic_json(directory/'requested.json', spec)
        preparing.unlink(missing_ok=True)
    return specs


def stage(root):
    root = Path(root).resolve(); root.mkdir(parents=True, exist_ok=True)
    with (root/'controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return _stage_unlocked(root)


def status(root):
    root = Path(root)
    result = {}
    for name in graphs(root.resolve()):
        directory = root/name
        if not (directory/'workflow.json').exists():
            result[name] = {'not_prepared': 1}; continue
        _, spec, state = dispatch._load(directory)
        result[name] = dict(Counter(dispatch._summary(spec, state).values()))
    return result


def run(root, interval=20):
    root = Path(root).resolve(); root.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        lock = stack.enter_context((root/'controller.lock').open('a'))
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print('Existing validation controller retains ownership.', flush=True); return
        specs = _stage_unlocked(root)
        for name in specs:
            handle = stack.enter_context((root/name/'drive.lock').open('a'))
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        while True:
            terminal = True; counts = {}; infrastructure = {}
            for name, spec in specs.items():
                try:
                    state = dispatch.advance(root/name, submit=True)
                except dispatch.SchedulerUnavailableError as error:
                    infrastructure[name] = str(error); terminal = False; continue
                current = dispatch._summary(spec, state)
                counts[name] = dict(Counter(current.values()))
                terminal &= all(s in ('complete', 'failed', 'blocked', 'skipped', 'unavailable') for s in current.values())
            atomic_json(root/'progress.json', dict(workflows=counts, infrastructure=infrastructure,
                execution_finished=terminal, scientific_approval=False))
            if terminal:
                atomic_json(root/'completion.json', dict(status='native_attempts_finished_review_required',
                    workflows=counts, scientific_approval=False, mess_updated=False))
                return
            time.sleep(interval)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('plan', 'stage', 'run', 'status'))
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    if args.action == 'plan':
        specs = graphs(args.root.resolve())
        print(json.dumps(dict(workflows=len(specs), jobs=sum(len(s['tasks']) for s in specs.values()),
            groups={name: len(spec['tasks']) for name,spec in specs.items()}), indent=2))
    elif args.action == 'stage':stage(args.root)
    elif args.action == 'run':run(args.root)
    else:print(json.dumps(status(args.root), indent=2))


if __name__ == '__main__':main()
