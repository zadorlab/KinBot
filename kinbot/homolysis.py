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


def repair_constraint_target(qc, obj, entry):
    """One owner-driven repair of an old scan frozen at the approximate guess.

    Preserve all completed points. Only a normal, strictly converged legacy
    point whose actual frozen coordinate explains the rejection is eligible.
    No SCF/optimizer failure or exhausted saddle search enters this procedure.
    """
    import ast
    import shutil
    from kinbot.execution import file_sha256
    path = Path(obj.instance_name+'.dissociation_scan.json')
    if not path.is_file():
        return None
    data = json.loads(path.read_text())
    if data.get('outcome') == 'constraint_target_repair_submitting':
        return data['constraint_target_repair']['step']  # journal-before-checkpoint interruption
    if data.get('outcome') != 'scan_constraint_not_met' or data.get('constraint_target_repair'):
        return None
    step = int(entry['step'])-1
    if (not 0 <= step < obj.max_step or data.get('step') != step+1
            or any(p['step'] >= step for p in data.get('points', []))
            or qc.check_qc(obj.instance_name) != 'normal'):
        return None
    source = Path(obj.instance_name+'.py')
    if not source.is_file():
        return None
    values = {}
    try:
        for node in ast.parse(source.read_text()).body:
            if not isinstance(node, ast.Assign):
                continue
            for target in node.targets:
                if not isinstance(target, ast.Name):
                    continue
                if target.id in ('strict_scan', 'constraint_targets'):
                    values[target.id] = ast.literal_eval(node.value)
                elif target.id == 'input_geom':
                    values[target.id] = ast.literal_eval(node.value.args[0])
                elif target.id == 'base_0_fix':
                    values[target.id] = ast.literal_eval(node.value.generators[0].iter)
        guess = np.asarray(values['input_geom'], dtype=float)
        if (values.get('strict_scan') is not True or 'constraint_targets' in values
                or values['base_0_fix'] != [[int(i)+1 for i in obj.instance]]
                or guess.shape != (obj.species.natom, 3) or not np.isfinite(guess).all()):
            return None
    except (SyntaxError, ValueError, TypeError, KeyError, AttributeError):
        return None
    frozen = float(np.linalg.norm(guess[obj.instance[0]]-guess[obj.instance[1]]))
    target = obj.start_distance + step*obj.par['homolytic_saddle_scan_step']
    error, geom = qc.get_qc_geom(obj.instance_name, obj.species.natom)
    geom = np.asarray(geom, dtype=float)
    if error or geom.shape != guess.shape or not np.isfinite(geom).all():
        return None
    actual = float(np.linalg.norm(geom[obj.instance[0]]-geom[obj.instance[1]]))
    if (abs(target-data['target_angstrom']) > 1.e-8
            or abs(actual-data['distance_angstrom']) > 1.e-5
            or abs(actual-frozen) > 1.e-4 or abs(frozen-target) <= .001):
        return None
    archive = Path('.kinbot_scan_repairs')/obj.instance_name/'explicit_target_v1'
    archive.mkdir(parents=True, exist_ok=True)
    evidence = {}
    for suffix in ('.py', '.com', '.log', '.traj', '_sella.log',
                   '.optimization_budget.json', '.dissociation_scan.json'):
        old = Path(obj.instance_name+suffix)
        if not old.is_file():
            continue
        saved = archive/old.name
        checksum = file_sha256(old)
        if saved.exists() and file_sha256(saved) != checksum:
            raise RuntimeError(f'Constraint-repair evidence changed: {saved}')
        if not saved.exists():
            shutil.copy2(old, saved)
        evidence[str(old)] = dict(archive=str(saved), sha256=checksum)
    repair = dict(step=step, reason='legacy worker froze approximate input instead of requested coordinate',
                  frozen_angstrom=frozen, requested_angstrom=target, evidence=evidence,
                  previous_outcome=data['outcome'])
    record(obj, outcome='constraint_target_repair_submitting', constraint_target_repair=repair)
    return step


def export_notes(directory='.'):
    path=Path(directory)/'homolysis_checks.json'
    if not path.exists():return []
    notes = []
    for c in json.loads(path.read_text())['checks']:
        if c['outcome'] == 'sampled_path_supports_provisional_pst':
            notes.append(f"{c['parent']}: provisional PST supported by a finite converged L1 scan; "
                         'no global absence-of-saddle proof.')
        elif c['outcome'] in ('unresolved', 'no_validated_saddle'):
            notes.append(f"{c['parent']}: no validated replacement saddle from the additional search "
                         f"({c.get('scan_outcome', 'not recorded')}); retain the capture-model "
                         'proposal unless superseded during pathway selection. This is not proof '
                         'of barrierlessness or a request to retry the search.')
    return notes


def unresolved_checks(directory='.'):
    """Return inconclusive evidence, not automatic export holds or retries."""
    path = Path(directory)/'homolysis_checks.json'
    if not path.exists():
        return []  # Legacy exports have no newly requested coverage revision.
    data = json.loads(path.read_text())
    return [c for c in data['checks'] if c['outcome'] not in
               ('superseded_by_validated_saddle', 'parent_not_retained',
                'sampled_path_supports_provisional_pst')]


def pending_checks(directory='.'):
    """Only unfinished discovery blocks export; a terminal nonfinding does not.

    Legacy ``unresolved`` ledgers were published after reaction generation
    finished. Preserve their evidence without perpetuating the old blanket
    hold. A newer active search checkpoint takes precedence over that ledger.
    The scan's ``saddle_submitting`` is an intent, not its final search status.
    """
    checks = unresolved_checks(directory)
    state_path = Path(directory)/'.kinbot_search_state.json'
    state = json.loads(state_path.read_text()) if checks and state_path.exists() else {}
    pending = []
    for c in checks:
        status = state.get('reactions', {}).get(c.get('search'), {}).get(
            'status', c.get('search_status'))
        if c['outcome'] == 'pending' or (status is not None and status >= 0):
            pending.append(c)
    return pending


def ensure_export_ready(directory='.'):
    path = Path(directory)/'homolysis_checks.json'
    pending = [c['parent'] for c in pending_checks(directory)]
    if pending:
        raise RuntimeError(f'{path}: unfinished homolysis searches: {pending}; '
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
    """Use the exporter rule: an accepted tight path with the same endpoints.

    This replaces a kinetic placeholder; it does not prove that every direct
    cleavage coordinate has a saddle. Requiring one particular bond edit here
    disagreed with final pathway selection and launched redundant companions.
    """
    from kinbot.reaction_path import reaction_path_id
    target = _products(parent)
    if not target:
        return []
    matches = []
    for i, obj in enumerate(species.reac_obj):
        if species.reac_ts_done[i] != -1 or species.reac_type[i] == 'hom_sci':
            continue
        if _products(obj) == target and reaction_path_id(obj):
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
    repair = data.get('constraint_target_repair', {})
    repairing = data['outcome'] == 'constraint_target_repair_submitting'
    if repairing and step > repair['step']:
        data = record(obj, outcome='scanning')
        repairing = False
    if step > 0 and data['outcome'] != 'saddle_submitting' and not repairing:
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
    """Publish discovery outcomes; terminal nonfindings do not block kinetics."""
    checks = []
    for index, obj in enumerate(species.reac_obj):
        if not hasattr(obj, 'homolysis_parent'):
            continue
        parent_i = species.reac_name.index(obj.homolysis_parent)
        parent = species.reac_obj[parent_i]
        matches = matching_saddles(species, parent) if species.reac_ts_done[parent_i] == -1 else []
        data = record(obj)
        if species.reac_ts_done[index] >= 0:
            outcome = 'pending'
        elif species.reac_ts_done[parent_i] != -1:
            outcome = 'parent_not_retained'
        elif matches:
            outcome = 'superseded_by_validated_saddle'
        elif (data['outcome']=='no_maximum_on_sampled_path'
              and data.get('asymptote',{}).get('status')=='supports_provisional_pst'):
            outcome = 'sampled_path_supports_provisional_pst'
        else:
            outcome = 'no_validated_saddle'
        checks.append(dict(parent=obj.homolysis_parent, search=obj.instance_name,
                           outcome=outcome, scan_outcome=data['outcome'], saddles=matches,
                           search_status=species.reac_ts_done[index]))
    if not checks:
        return
    atomic_json('homolysis_checks.json', dict(schema=1, checks=checks))
    return checks
