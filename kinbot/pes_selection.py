"""Outer-PES L2 pathway selection, before any new automatic ANL work.

Selection is a kinetics policy, not proof that two saddle geometries are the
same. Raw accepted chemistry and QC ownership are kept separately.
"""
from copy import copy, deepcopy
import json
import math
from pathlib import Path
import uuid

from kinbot.execution import atomic_json
from kinbot.anl.ts_reuse import digest

NAME = 'pes_l2_selection.json'


def choose(snapshots):
    """Use the same endpoint/pathway rule as final PES assembly, at L2 E+ZPE."""
    from kinbot.pes import select_summary_reaction
    paths, candidates, owners = {}, [], {}
    theories = set()
    for worker, evidence in sorted(snapshots.items()):
        if evidence.get('errors'):
            raise ValueError(f'{worker}: incomplete accepted L2 evidence: {evidence["errors"]}')
        network = evidence['network']
        stable = {s['routing_name']: s for s in network['stable_species']}
        for row in network['reactions']:
            name = row['name']
            if name in owners:
                raise ValueError(f'{name}: reaction name occurs in multiple workers.')
            owners[name] = worker
            if row['barrierless']:
                records = [stable[p] for p in row['products']]
            else:
                records = [network['transition_states'][row['transition_state']]]
                paths[name] = row.get('path_class')
            energy = 0.
            for record in records:
                theories.add(digest(record['l2_theory']))
                energy += record['energy_hartree'] + record['harmonic_zpe_hartree']
            if not math.isfinite(energy):
                raise ValueError(f'{name}: nonfinite L2 energy.')
            candidate = [row['reactant'], name, row['products'], energy]
            if row.get('product_complex'):
                # Final PES assembly prefers a complex-bearing route on an
                # exact energy tie. The last two fields carry that flag only.
                candidate.extend([0., 'vdW_selection'])
            candidates.append(candidate)
    if len(theories) > 1:
        raise ValueError('PES pathway selection requires a common L2 theory.')
    retained = []
    for candidate in candidates:
        select_summary_reaction(retained, candidate, paths)
    selected = {row[1] for row in retained}
    result = {}
    for worker, evidence in snapshots.items():
        names = {r['name'] for r in evidence['network']['reactions']}
        result[worker] = dict(network_sha256=digest(evidence['network']),
            selected=sorted(names & selected), omitted={},
            complexes=sorted({r['product_complex'] for r in evidence['network']['reactions']
                              if r['name'] in selected and r.get('product_complex')}))
        for candidate in candidates:
            if owners[candidate[1]] != worker or candidate[1] in selected:
                continue
            # Re-evaluate against the final retained set for an auditable
            # explanation, including reverse discoveries and placeholders.
            matches = []
            for winner in retained:
                trial = [deepcopy(winner)]
                select_summary_reaction(trial, deepcopy(candidate), paths)
                if all(r[1] != candidate[1] for r in trial):
                    matches.append(winner[1])
            if not matches:
                raise ValueError(f'{candidate[1]}: omitted route has no retained replacement.')
            result[worker]['omitted'][candidate[1]] = dict(
                retained=sorted(matches), reason='outer_pes_lowest_l2_pathway',
                energy_hartree=candidate[3])
    return dict(schema=1, status='selected', policy='lowest_l2', workers=result)


def publish(root, phase, workers):
    """PES owner only, under its exploration lock; no scheduler/QC changes."""
    root = Path(root)
    snapshots = {}
    managed, input_errors = [], []
    for name in sorted(workers):
        inp = root/name/(name+'.json')
        try:
            parameters = json.loads(inp.read_text())
        except (OSError, ValueError) as error:
            input_errors.append(f'{name}: {error}')
            continue
        if parameters.get('composite_method'):
            managed.append(name)
    if not managed and phase.get('anl_selection_policy') != 'lowest_l2':
        return
    try:
        if input_errors or set(managed) != set(workers):
            raise ValueError('Incomplete composite PES worker inventory: '+
                             '; '.join(input_errors or sorted(set(workers)-set(managed))))
        for name in managed:
            evidence = json.loads((root/name/'accepted_l2_network.json').read_text())
            if evidence.get('generation') != dict(token=phase['token'], epoch=phase['epoch']):
                raise ValueError(f'{name}: L2 snapshot belongs to another exploration generation.')
            snapshots[name] = evidence
        result = choose(snapshots)
    except (OSError, ValueError, KeyError, TypeError) as error:
        # A bad ANL input must not close the independent VRC review/execution
        # path. ANL fails at its planning boundary, before submitting QC.
        result = dict(schema=1, status='blocked', error=str(error))
    result.update(token=phase['token'], epoch=phase['epoch'])
    old = root/NAME
    if old.exists():
        previous = json.loads(old.read_text())
        if previous != result:
            atomic_json(root/'.pes_l2_selection_history'/(digest(previous)+'.json'), previous)
    atomic_json(old, result)


def selected_view(species, parameters):
    """Mask only downstream use; leave the exploration objects/statuses intact."""
    from kinbot import pes_phase
    phase = pes_phase._read()
    if phase is None or not parameters.get('composite_method'):
        return species
    root, worker = pes_phase._ROOT, pes_phase._WORKER
    selection = json.loads((root/NAME).read_text())
    if (selection.get('token'), selection.get('epoch')) != (phase['token'], phase['epoch']):
        raise ValueError('ANL selection belongs to another PES exploration generation.')
    if selection.get('status') != 'selected':
        raise ValueError('Outer PES L2 selection is blocked: '+selection.get('error', 'missing evidence'))
    entry = selection['workers'][worker]
    from kinbot.pes_snapshot import snapshot
    if digest(snapshot(species, parameters)['network']) != entry['network_sha256']:
        raise ValueError('Accepted L2 chemistry changed after outer PES pathway selection.')
    view = copy(species)
    view.reac_ts_done = list(species.reac_ts_done)
    view.anl_pathway_selection = deepcopy(entry)
    view.anl_accepted_source = species
    for i, reaction in enumerate(species.reac_obj):
        if view.reac_ts_done[i] == -1 and reaction.instance_name not in entry['selected']:
            view.reac_ts_done[i] = -2  # Selection omission, never a QC failure.
    return view


def task_names(plan, aliases, selection=None):
    """Return selected saddles plus the exact-identity owners they consume."""
    active = set(plan['network']['transition_states'])
    if selection is not None:
        active &= set(selection['selected'])
    return active | {aliases[n] for n in active if n in aliases}


def project(plan, selection):
    """Restrict stable/CBH work to the selected network, retaining its QC inputs.

    Keep all TS records until strict alias resolution; an omitted observation
    can still own the identical saddle needed by a selected observation.
    """
    if selection is None:
        return plan
    from kinbot.anl.cbh import generate_cbh_reaction, cbh_network_requirements, cbh_network_states
    result = deepcopy(plan)
    network = result['network']
    selected = [r for r in network['reactions'] if r['name'] in selection['selected']]
    routes = {network['initial'], *(p for r in selected for p in r['products'])}
    network['stable_species'] = [s for s in network['stable_species'] if s['routing_name'] in routes]
    network['route_to_smiles'] = {n:s for n,s in network['route_to_smiles'].items() if n in routes}
    if 'product_complexes' in network:
        network['product_complexes'] = {n:e for n,e in network['product_complexes'].items()
                                        if n in selection.get('complexes', ())}
    cbh = result['cbh']
    cbh['targets'] = {n:e for n,e in cbh['targets'].items() if n in routes}
    calculated, terminals = set(), set()
    for record in network['stable_species']:
        molecule = record['molecule']
        target = cbh['targets'][record['routing_name']]
        if target.get('terminal_atct'):
            calculated.add((record['smiles'], molecule['charge'], molecule['multiplicity']))
            terminals.add(record['smiles'])
        else:
            for rung in target['available_structural_rungs']:
                reaction = generate_cbh_reaction(record['smiles'], rung,
                    charge=molecule['charge'], multiplicity=molecule['multiplicity'])
                if reaction is None:
                    raise ValueError('Published CBH rung disappeared during pathway selection.')
                states = cbh_network_states(reaction)
                calculated.update((s,*state) for s,state in states.items())
                _, leaves = cbh_network_requirements(reaction)
                terminals.update(leaves)
    identity = lambda r: (r['smiles'],r['charge'],r['multiplicity'])
    cbh['calculated_states'] = [r for r in cbh['calculated_states'] if identity(r) in calculated]
    if {identity(r) for r in cbh['calculated_states']} != calculated:
        raise ValueError('Selected CBH states are absent from the frozen accepted plan.')
    cbh['terminal_atct_species'] = {s:r for s,r in cbh['terminal_atct_species'].items() if s in terminals}
    if set(cbh['terminal_atct_species']) != terminals:
        raise ValueError('Selected ATcT terminals are absent from the accepted plan.')
    thermo = calculated | {(s,*r['state']) for s,r in cbh['terminal_atct_species'].items()}
    cbh['thermochemistry_states'] = [r for r in cbh['thermochemistry_states'] if identity(r) in thermo]
    if {identity(r) for r in cbh['thermochemistry_states']} != thermo:
        raise ValueError('Selected thermochemistry states are absent from the accepted plan.')
    return result


def graph_keys(plan, needed_ts):
    from kinbot.anl.automatic import _state_key
    return ({'stable:'+_state_key(r['smiles'],r['charge'],r['multiplicity'])
             for r in plan['cbh']['calculated_states']} |
            {'complex:'+n for n in plan['network'].get('product_complexes',{})} | set(needed_ts))


def bind(base, plan, selection):
    """Invalidate derived assembly only when its selected chemistry changes."""
    from kinbot.anl.network_update import archive_assembly
    base = Path(base)
    value = dict(schema=1, selected=(sorted(plan['network']['transition_states'])
        if selection is None else sorted(selection['selected'])))
    path = base/'kinetic_selection.json'
    if selection is None and not path.exists():
        return  # Unmanaged/master-style execution has no new selection step.
    previous = json.loads(path.read_text()) if path.exists() else None
    journal = base/'selection_update.json'
    pending = json.loads(journal.read_text()) if journal.exists() else {}
    resuming = pending.get('status') == 'archiving'
    if resuming:
        if previous not in (pending['before'], pending['after']) or pending['after'] != value:
            raise ValueError('An interrupted pathway selection must finish before another revision.')
    if resuming or (previous != value and any((base/n).exists() for n in
            ('execution.json','handoff.json','profiles','formations',
             'final_report.json','cbh_thermochemistry_models.json'))):
        if not resuming:
            pending = dict(status='archiving', before=previous, after=value,
                           archive=str(base/'selection_history'/uuid.uuid4().hex))
            atomic_json(journal,pending)
        archive = Path(pending['archive'])
        archive_assembly(base, archive)
        atomic_json(archive/'selection.json',pending)
        atomic_json(path,value)
        atomic_json(journal,dict(pending,status='complete'))
        return
    atomic_json(path, value)
