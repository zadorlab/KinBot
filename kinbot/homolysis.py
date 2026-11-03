"""Check homolysis placeholders after ordinary reaction discovery.

The graph operation that proposes two fragments is not a saddle search.
Companions use the ordinary TS/frequency/IRC/L2 machinery, but never interpret
a failed or exhausted scan as proof that a channel is barrierless.
"""
import json
import math
from pathlib import Path

import numpy as np

from kinbot.execution import atomic_json
from kinbot.species_routing import routing_key


FAMILY = 'dissociation_saddle'


def export_notes(directory='.'):
    path=Path(directory)/'homolysis_checks.json'
    if not path.exists():return []
    return [f"{c['parent']}: provisional PST supported by a finite converged L1 scan; "
            'no global absence-of-saddle proof.' for c in json.loads(path.read_text())['checks']
            if c['outcome']=='sampled_path_supports_provisional_pst']


def ensure_export_ready(directory='.'):
    path = Path(directory)/'homolysis_checks.json'
    if not path.exists():
        return  # Legacy exports have no newly requested coverage revision.
    data = json.loads(path.read_text())
    pending = [c['parent'] for c in data['checks'] if c['outcome'] not in
               ('superseded_by_validated_saddle', 'parent_not_retained',
                'sampled_path_supports_provisional_pst')]
    if pending:
        raise RuntimeError(f'{path}: unresolved homolysis checks: {pending}; '
                           'old summaries cannot authorize a new MESS/PST export.')


def validate(selection):
    if isinstance(selection, bool):
        return selection
    if (not isinstance(selection, list) or not selection or
            any(not isinstance(n, str) or '_hom_sci_' not in n or
                '/' in n or '\\' in n for n in selection)):
        raise ValueError('homolytic_saddle_search must be true, false, or exact hom_sci names.')
    return sorted(set(selection))


def add_companions(finder):
    """Attach after ordinary selection/naming so no extra motif is invented."""
    from kinbot.reactions.reac_dissociation_saddle import DissociationSaddle
    selection = validate(finder.par.get('homolytic_saddle_search', True))
    species = finder.species
    for i, family in enumerate(list(species.reac_type)):
        if family != 'hom_sci':
            continue
        original = species.reac_name[i]
        if selection is False or isinstance(selection, list) and original not in selection:
            continue
        name = original.replace('_hom_sci_', '_dissociation_saddle_', 1)
        obj = DissociationSaddle(species, finder.qc, finder.par,
                                 list(species.reac_inst[i]), name)
        obj.homolysis_parent = original
        species.reac_obj.append(obj)
        species.reac_name.append(name)
        species.reac_type.append(FAMILY)
        species.reac_inst.append(list(species.reac_inst[i]))
        species.reac_step.append(0)
        species.reac_ts_done.append(0)
        species.reac_scan_energy.append([])
        species.reac_ts_geom.append(0)
        species.reac_ts_freq.append(0)


def _products(obj):
    return sorted(str(routing_key(p)) for p in obj.products)


def matching_saddles(species, parent):
    from kinbot.reaction_path import reaction_path_id
    target = _products(parent)
    if not target:
        return []
    matches = []
    for i, obj in enumerate(species.reac_obj):
        if species.reac_ts_done[i] != -1 or species.reac_type[i] == 'hom_sci':
            continue
        # The same separated products can also be reached by a rearrangement.
        # Only a validated IRC breaking this bond (without other connectivity
        # edits) covers this particular direct-cleavage proposal. Atom order
        # here is the common parent/IRC order, not a guessed SMILES mapping.
        endpoint = getattr(obj, 'irc_prod', None)
        if endpoint is None:
            continue
        before = np.asarray(species.bond) > 0
        after = np.asarray(endpoint.bond) > 0
        broken = np.argwhere(np.triu(before & ~after, 1)).tolist()
        formed = np.argwhere(np.triu(after & ~before, 1)).tolist()
        direct = broken == [sorted(map(int, parent.instance))] and not formed
        if direct and _products(obj) == target and reaction_path_id(obj):
            matches.append(obj.instance_name)
    return matches


def deferred(generator, index):
    """Finish ordinary discovery/reconstruction first; reuse any matching TS."""
    species = generator.species
    obj = species.reac_obj[index]
    if species.reac_type[index] != FAMILY:
        return False
    if any(s >= 0 and f != FAMILY for s, f in
           zip(species.reac_ts_done, species.reac_type)):
        return True
    parent_index = species.reac_name.index(obj.homolysis_parent)
    parent = species.reac_obj[parent_index]
    if species.reac_ts_done[index] != 0 or species.reac_step[index] != 0:
        return False  # Adopt an existing attempt; never retire live work here.
    if species.reac_ts_done[parent_index] != -1:
        species.reac_ts_done[index] = -999
        record(obj, outcome='parent_not_retained')
        return True
    matches = matching_saddles(species, parent)
    if matches:
        species.reac_ts_done[index] = -999
        record(obj, outcome='covered_by_validated_saddle', saddles=matches)
        return True
    return False


def record(obj, **changes):
    path = Path(obj.instance_name + '.dissociation_scan.json')
    result = json.loads(path.read_text()) if path.exists() else dict(
        schema=1, reaction=obj.instance_name, parent=obj.homolysis_parent,
        points=[], outcome='pending', proves_global_absence_of_saddle=False)
    result.update(changes)
    atomic_json(path, result)
    return result


def sampled_asymptote(generator, obj, data):
    """Finite-path scouting evidence, explicitly not a global no-saddle proof."""
    from kinbot import constants, l1_screening
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.electronic_states import apply_electronic_state_model
    parent = generator.species.reac_obj[generator.species.reac_name.index(obj.homolysis_parent)]
    end = StationaryPoint('scan_endpoint', obj.species.charge, obj.species.mult,
        atom=obj.species.atom, geom=np.array(data['points'][-1]['positions']))
    end.characterize()
    fragments, _ = end.start_multi_molecular()
    for fragment in fragments:
        apply_electronic_state_model(fragment, generator.par)
    if sorted(str(routing_key(p)) for p in fragments) != _products(parent):
        return dict(status='unresolved', reason='Scan endpoint does not match the proposed fragment states.')
    points = l1_screening.read_points(Path.cwd(),
        [str(routing_key(p))+'_well' for p in parent.products],
        l1_screening.theory(generator.par))
    tail = [p['energy_hartree'] for p in data['points'][-3:]]
    delta = (tail[-1]-sum(p['energy_hartree'] for p in points))*constants.AUtoKCAL
    span = (max(tail)-min(tail))*constants.AUtoKCAL
    supported = (len(tail)==3 and abs(delta)<=generator.par['homolytic_saddle_asymptote_tolerance']
                 and span<=generator.par['homolytic_saddle_flatness'])
    return dict(status='supports_provisional_pst' if supported else 'unresolved',
        fragment_points=points, endpoint_minus_fragments_kcal_mol=delta,
        final_three_span_kcal_mol=span,
        asymptote_tolerance_kcal_mol=generator.par['homolytic_saddle_asymptote_tolerance'],
        flatness_tolerance_kcal_mol=generator.par['homolytic_saddle_flatness'],
        scope='Converged finite L1 path only; no exhaustive saddle or electronic-state search.')


def advance(generator, index):
    """One bounded scan transition, using the normal durable submission owner.

    Return False after the final saddle submission, so the usual frequency
    validation, IRC, product checks and high-level refinement take over.
    """
    from kinbot.reac_family import carry_out_reaction
    species, qc = generator.species, generator.qc
    obj = species.reac_obj[index]
    step = species.reac_step[index]
    if step == obj.max_step + 1:
        return False
    status = qc.check_qc(obj.instance_name)
    if status == 'running':
        return True
    if status == 'error':
        record(obj, outcome='qc_failed', step=step)
        species.reac_ts_done[index] = -999
        return True
    data = record(obj)
    if step > 0 and data['outcome'] != 'saddle_submitting':
        if status != 'normal':
            return True  # A staged/admitted result may not yet be visible.
        err, energy = qc.get_qc_energy(obj.instance_name)
        geom_err, geom = qc.get_qc_geom(obj.instance_name, species.natom)
        if err or geom_err:
            return True
        geom = np.asarray(geom, dtype=float)
        if not math.isfinite(energy) or geom.shape != (species.natom, 3) or not np.isfinite(geom).all():
            record(obj, outcome='invalid_scan_point', step=step)
            species.reac_ts_done[index] = -999
            return True
        completed_step = step - 1
        distance = float(np.linalg.norm(geom[obj.instance[0]]-geom[obj.instance[1]]))
        target = obj.start_distance + completed_step * obj.par['homolytic_saddle_scan_step']
        if abs(distance-target) > 0.001:
            record(obj, outcome='scan_constraint_not_met', step=step,
                   distance_angstrom=distance, target_angstrom=target)
            species.reac_ts_done[index] = -999
            return True
        if not any(p['step'] == completed_step for p in data['points']):
            from kinbot.search_state import _job
            data['points'].append(dict(step=completed_step, energy_hartree=float(energy),
                positions=geom.tolist(), job_id=_job(qc, obj.instance_name),
                distance_angstrom=float(np.linalg.norm(geom[obj.instance[0]]-geom[obj.instance[1]]))))
            data = record(obj, points=data['points'])
            generator.search_checkpoint.scan_energy(index, energy)
        energies = [p['energy_hartree'] for p in data['points']]
        if len(energies) >= 3 and energies[-2] > max(energies[-3], energies[-1]):
            # A maximum supplies a guess, never an accepted saddle. Preserve
            # this intent before admission: pool-full/restart must not turn
            # max_step into an exhausted scan or submit it twice.
            data = record(obj, outcome='saddle_submitting', peak_step=step-2)
            species.reac_step[index] = step = obj.max_step
            generator.search_checkpoint.save(index)
        elif step >= obj.max_step:
            try:
                asymptote = sampled_asymptote(generator, obj, data)
            except (OSError, ValueError, KeyError) as error:
                asymptote = dict(status='unresolved', reason=str(error))
            record(obj, outcome='no_maximum_on_sampled_path', asymptote=asymptote)
            species.reac_ts_done[index] = -999
            return True
    if data['outcome'] == 'saddle_submitting':
        species.reac_step[index] = step = obj.max_step
        obj.scan_saddle_geometry = next(p['positions'] for p in data['points']
                                        if p['step'] == data['peak_step'])
    species.reac_step[index] = carry_out_reaction(obj, step,
        generator.par['qc_command'], bimol=generator.par['bimol'])
    if species.reac_step[index] == -1:
        record(obj, outcome='geometry_or_submission_failed')
        species.reac_ts_done[index] = -999
    return True


def finish(species):
    """Write outcomes and prevent unresolved assumptions reaching rate export."""
    checks = []
    for obj in species.reac_obj:
        if not hasattr(obj, 'homolysis_parent'):
            continue
        parent_i = species.reac_name.index(obj.homolysis_parent)
        parent = species.reac_obj[parent_i]
        matches = matching_saddles(species, parent) if species.reac_ts_done[parent_i] == -1 else []
        data = record(obj)
        if species.reac_ts_done[parent_i] != -1:
            outcome = 'parent_not_retained'
        elif matches:
            outcome = 'superseded_by_validated_saddle'
        elif (data['outcome']=='no_maximum_on_sampled_path'
              and data.get('asymptote',{}).get('status')=='supports_provisional_pst'):
            outcome = 'sampled_path_supports_provisional_pst'
        else:
            outcome = 'unresolved'
        checks.append(dict(parent=obj.homolysis_parent, search=obj.instance_name,
                           outcome=outcome, scan_outcome=data['outcome'], saddles=matches))
    if not checks:
        return
    atomic_json('homolysis_checks.json', dict(schema=1, checks=checks))
    unresolved = [c['parent'] for c in checks if c['outcome'] == 'unresolved']
    if unresolved:
        raise RuntimeError('Homolysis inner-barrier checks remain unresolved: ' +
            ', '.join(unresolved) + '. See homolysis_checks.json and native scan/TS/IRC evidence. '
            'No new MESS/PST export was authorized; a failed or finite scan is not proof of barrierlessness.')
