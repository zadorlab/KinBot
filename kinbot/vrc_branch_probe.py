"""Bounded fragment-seeded VRC diagnostic; never accepts production scan points.

The native profile, state average and final fragment orbitals come from a
verified production asymptote. Only separation and execution budget change.
One task uses the existing shared admission and bounded fragment-retry policy.
"""
from copy import deepcopy
import fcntl
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from ase.io import read as read_atoms

from kinbot.execution import file_sha256
from kinbot.vrc_scan_recovery import digest


def read(path):
    return json.loads(Path(path).read_text())


def specification(request, distance, *, steps=50, walltime='01:00:00'):
    from kinbot.anl.dispatch import _load, validate_spec
    from kinbot.run_policy import find
    from kinbot.run_walltime import seconds
    from kinbot.vrc_tst_scan import VTS

    request = Path(request).resolve(strict=True)
    saved = read(request)
    context = saved['context']
    if digest(context) != saved.get('context_sha256'):
        raise ValueError('Scan context hash is missing or changed.')
    policy = find(request.parent)
    if policy is None:
        raise ValueError('Diagnostic must share the production run resource policy.')
    root = policy.parent
    if (not np.isfinite(distance) or not 0. < distance < context['separation_angstrom']
            or type(steps) is not int or not 1 <= steps <= 100
            or not 0 < seconds(walltime) <= 7200):
        raise ValueError('Need an interior distance, 1–100 steps and at most two hours.')
    model = deepcopy(context['state_model'])
    selection = model['orbital_selection']
    checkpoint = Path(selection['path']).resolve(strict=True)
    graph = checkpoint.parents[2]
    directory = graph/'tasks/relaxed_scan'
    if (not graph.is_relative_to(root) or checkpoint.parent != directory
            or '_vts_fragments_asymptote_' not in graph.name
            or selection.get('mode') != 'molpro_record'
            or file_sha256(checkpoint) != selection.get('sha256')
            or (directory/'scan_rejection.json').exists()):
        raise ValueError('Need an intact, hash-bound production fragment checkpoint.')
    _, source, state = _load(graph)
    if state['tasks']['relaxed_scan']['status'] != 'complete':
        raise ValueError('Fragment source is not complete.')
    execution = read(directory/'execution.json')
    if execution.get('status') != 'executed':
        raise ValueError('Fragment execution did not succeed.')
    # Check the artifacts we use, without reading every historical force file.
    required = [directory/'task.json', directory/'final.xyz', checkpoint]
    for path in required:
        if execution.get('artifacts', {}).get(path.name) != file_sha256(path):
            raise ValueError('Fragment artifact changed or lacks provenance: '+path.name)
    task = deepcopy(read(directory/'task.json')['task'])
    if (len(source['tasks']) != 1 or task != source['tasks'][0]
            or task.get('manual_fallback') or task.get('kind') != 'ase_optimize'
            or task['profile']['calculator'] != 'molpro'):
        raise ValueError('Expected one unchanged production Molpro scan task.')
    profile = task['profile']
    parameters = context['parameters']
    if (profile['basis'] != parameters['vrc_tst_scan_basis']
            or profile['method'].lower() != parameters['vrc_tst_scan_method'].lower()):
        raise ValueError('Fragment method/basis differs from current scan context.')
    options = profile['calculator_kwargs']
    for key, parameter in (('caspt2_shift', 'vrc_tst_scan_shift'),
                           ('caspt2_ipea', 'vrc_tst_scan_ipea'),
                           ('gradient_mode', 'vrc_tst_scan_gradient'),
                           ('scf_convergence', 'vrc_tst_scf_convergence'),
                           ('correlation_convergence', 'vrc_tst_correlation_convergence')):
        if options.get(key) != parameters.get(parameter):
            raise ValueError('Fragment numerical contract differs from scan context: '+key)
    if task['optimizer']['fmax'] != parameters['vrc_tst_scan_fmax']:
        raise ValueError('Fragment force threshold differs from scan context.')
    without_selection = lambda value: {k:v for k,v in value.items() if k != 'orbital_selection'}
    if digest(without_selection(options['state_model'])) != digest(without_selection(model)):
        raise ValueError('Fragment electronic-state contract changed.')
    atoms = read_atoms(directory/'final.xyz')
    initial = np.asarray(context['fragment_geometry'], dtype=float)
    symbols = context['symbols']
    if (atoms.get_chemical_symbols() != symbols or initial.shape != atoms.positions.shape
            or not np.isfinite(initial).all()
            or not np.allclose(initial, atoms.positions, atol=1.e-7, rtol=0.)
            or any(source['molecule'][key] != context[key]
                   for key in ('symbols', 'charge', 'multiplicity'))):
        raise ValueError('Fragment coordinates or molecular identity changed.')
    maps, pair = context['maps'], context['coordinate']
    if (len(maps) != 2 or any(not m for m in maps)
            or sorted(i for m in maps for i in m) != list(range(len(atoms)))
            or len(pair) != 2 or pair[0] == pair[1]
            or not any(pair[0] in maps[i] and pair[1] in maps[1-i] for i in (0, 1))):
        raise ValueError('Fragment maps do not partition the selected coordinate.')
    if not np.isclose(np.linalg.norm(initial[pair[0]]-initial[pair[1]]),
                      context['separation_angstrom'], atol=.01, rtol=0.):
        raise ValueError('Fragment reference separation is inconsistent.')
    constraints = task['optimizer']['constraints']
    if constraints['fixed_bonds'] != [pair]:
        raise ValueError('Source scan constraint differs from review coordinate.')
    for minimum in constraints.get('minimum_bonds', []):
        if (len(minimum) != 3 or not np.isclose(minimum[2], context['separation_angstrom'])
                or not any(minimum[0] in maps[i] and minimum[1] in maps[1-i] for i in (0, 1))):
            raise ValueError('Unexpected source minimum-distance constraint.')
        minimum[2] = float(distance)
    positions = VTS._shift_scan_geometry(SimpleNamespace(scan_coo=pair, maps=maps), initial, distance)
    molecule = {key:deepcopy(context[key]) for key in ('symbols', 'charge', 'multiplicity')}
    molecule['positions'] = positions.tolist()
    options['state_model'] = model
    task['optimizer'].update(steps=steps, step_retries=0)
    task['resources'].update(walltime=walltime)
    task['resources'].pop('requested_walltime', None)
    proof = dict(request=str(request), context_sha256=saved['context_sha256'],
                 reaction=saved['reaction'], distance_angstrom=float(distance), maps=maps,
                 coordinate=pair, fragment_geometry=initial.tolist(),
                 source_hashes={str(p):file_sha256(p) for p in
                                [graph/'workflow.json', directory/'execution.json', *required]})
    spec = dict(schema=1, name='vrc-fragment-branch-probe', molecule=molecule,
                limits={'max_nodes':1}, tasks=[task], intent={
                    'vrc_branch_probe':proof, 'production_model_accepted':False,
                    # Reuse the dispatcher's established one-attempt/time-cap contract.
                    'vrc_fragment_retry':True, 'cap_walltime_to_explicit_partition':True})
    validate_spec(spec)
    return root, spec


def prepare(request, distance, **options):
    from kinbot.vrc_graph import select_graph
    root, spec = specification(request, distance, **options)
    parent = root/'diagnostics/vrc_branch_probes'
    parent.mkdir(parents=True, exist_ok=True)
    # select_graph's thread lock is sufficient inside production's single
    # owner; this standalone entry point also needs a cross-process lock.
    with (parent/'.prepare.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        return select_graph(parent, 'probe', spec)


def status(graph):
    from kinbot.anl.dispatch import _load
    graph, spec, state = _load(graph)
    proof = spec['intent']['vrc_branch_probe']
    directory = graph/'tasks/relaxed_scan'
    result = dict(graph=str(graph), reaction=proof['reaction'],
                  distance_angstrom=proof['distance_angstrom'],
                  task=state['tasks']['relaxed_scan'], production_model_accepted=False)
    outcome = directory/'execution.json'
    execution = read(outcome) if outcome.exists() else {}
    result['execution'] = {k:execution.get(k) for k in ('status', 'error')}
    if (directory/'optimization.log').exists():
        with (directory/'optimization.log').open('rb') as stream:
            stream.seek(max(0, (directory/'optimization.log').stat().st_size-1800))
            result['optimizer'] = stream.read().decode(errors='replace').splitlines()[-5:]
    if execution.get('status') == 'executed':
        final = directory/'final.xyz'
        if execution.get('artifacts', {}).get(final.name) != file_sha256(final):
            raise ValueError('Probe final geometry hash changed.')
        atoms = read_atoms(final)
        if atoms.get_chemical_symbols() != spec['molecule']['symbols']:
            raise ValueError('Probe atom order changed.')
        initial = np.asarray(proof['fragment_geometry'])
        result['energy_hartree'] = execution['details']['energy_hartree']
        result['pair_distances_angstrom'] = [dict(atoms=[i+1,j+1],
            symbols=[atoms[i].symbol,atoms[j].symbol],
            reference=float(np.linalg.norm(initial[i]-initial[j])),
            final=float(np.linalg.norm(atoms.positions[i]-atoms.positions[j])))
            for i in range(len(atoms)) for j in range(i+1,len(atoms))]
        result['review'] = 'Convergence is diagnostic evidence; inspect fragment/state identity before use.'
    return result


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run', 'status'])
    parser.add_argument('path', type=Path, help='Scan review JSON; for status, graph or production root.')
    parser.add_argument('--distance', type=float, default=6.)
    parser.add_argument('--steps', type=int, default=50)
    parser.add_argument('--walltime', default='01:00:00')
    args = parser.parse_args(argv)
    if args.action == 'status':
        graphs = ([args.path] if (args.path/'workflow.json').exists() else
                  sorted(p.parent for p in args.path.glob('diagnostics/vrc_branch_probes/probe_*/workflow.json')))
        print(json.dumps([status(g) for g in graphs], indent=2))
        return 0
    graph = prepare(args.path, args.distance, steps=args.steps, walltime=args.walltime)
    print('GRAPH', graph, flush=True)
    code = 0
    if args.action == 'run':
        from kinbot.anl.dispatch import drive_run, preflight_until_scheduler_ready
        preflight_until_scheduler_ready(graph)
        code, _ = drive_run(graph, hold_failures=False)
    print(json.dumps(status(graph), indent=2), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
