"""Read-only L1 energy/vdW selection audit of an existing PES directory.

This is an evidence report, never a retry, network edit, or seed manifest.
Unknown evidence stays in the possible retained set. Connectivity refers only
to the saved accepted graph, not to reactions an unfinished search may discover.
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sqlite3

from kinbot import constants, l1_screening as l1


def read(path, default=None):
    return json.loads(path.read_text()) if path.is_file() else (default if default is not None else {})


def unknown(reason):
    return dict(decision='unknown', reason=str(reason))


def combined(*items):
    values = [x['decision'] for x in items if x is not None]
    return 'exclude' if 'exclude' in values else 'unknown' if 'unknown' in values else 'retain'


def search_entries(worker):
    record = read(worker/'.kinbot_search_state.json')
    if record:
        return record.get('reactions', {})
    # Old monitors are evidence only. Do not import/publish a checkpoint here.
    path = worker/'kinbot_monitor.out'
    entries = {}
    if path.is_file():
        for line in path.read_text().splitlines():
            fields = line.split()
            if len(fields) >= 3:
                try:
                    entries[fields[2]] = dict(status=int(fields[0]), step=int(fields[1]), legacy=True)
                except ValueError:
                    pass
    return entries


def species(point, par):
    import numpy as np
    from kinbot.stationary_pt import StationaryPoint
    positions = np.asarray(point['positions'], dtype=float)
    if (positions.shape != (len(point['symbols']), 3) or not np.isfinite(positions).all()
            or any(np.linalg.norm(a-b) < .3 for i, a in enumerate(positions) for b in positions[:i])):
        raise ValueError('Invalid L1 geometry.')
    obj = StationaryPoint(point['name'], par['charge'], par['mult'],
                          atom=point['symbols'], geom=positions)
    obj.characterize()
    return obj


def product_side(worker, name, par, profile, parent):
    """Use the same two-endpoint identity test as IRC.irc2stationary_pt."""
    from kinbot.stereo_identity import require_supported_identity, identity_matches
    from kinbot.species_routing import routing_name
    from kinbot.electronic_states import apply_electronic_state_model
    if par.get('bimol') or par['charge'] != 0:
        raise ValueError('Bimolecular/charged endpoint assignment needs a separate identity audit.')
    points = l1.read_points(worker, [name+'_IRC_F_prod', name+'_IRC_R_prod'], profile,
                            geometry=True, require_zpe=False)
    reference = require_supported_identity(parent)
    objects = [species(p, par) for p in points]
    matches = [identity_matches(reference, require_supported_identity(o),
                               par.get('optical_population', 'specified')) for o in objects]
    if sum(matches) != 1:
        raise ValueError('The two L1 IRC endpoints do not identify exactly one parent well.')
    index = matches.index(False)
    fragments, _ = objects[index].start_multi_molecular()
    for fragment in fragments:
        apply_electronic_state_model(fragment, par)
    return points[index], sorted(routing_name(p) for p in fragments)


def depth_decision(endpoint, fragments, threshold):
    if Counter(endpoint['symbols']) != Counter(s for p in fragments for s in p['symbols']):
        return unknown('L1 complex and isolated fragments have different stoichiometry.')
    evidence = dict(threshold_kcal_mol=threshold, endpoint=endpoint, fragments=fragments)
    if any(not l1.finite_number(p.get('energy_hartree')) for p in [endpoint, *fragments]):
        return dict(unknown('Nonfinite electronic energy.'), **evidence)
    electronic = (sum(p['energy_hartree'] for p in fragments)
                  - endpoint['energy_hartree'])*constants.AUtoKCAL
    def comparison(depth, basis):
        return dict(decision='retain' if depth > threshold else 'exclude',
                    depth_kcal_mol=depth, threshold_kcal_mol=threshold, basis=basis)
    evidence['electronic_only'] = comparison(electronic, 'electronic_energy_without_ZPE')
    if any(not l1.finite_number(p.get('zpe_hartree')) for p in [endpoint, *fragments]):
        # QC.get_qc_zpe historically returns zero for an absent value. Show
        # that counterfactual explicitly, never as a matched E+ZPE decision.
        if (endpoint.get('zpe_status') == 'missing'
                and all(l1.finite_number(p.get('zpe_hartree')) for p in fragments)):
            old_depth = electronic + sum(p['zpe_hartree'] for p in fragments)*constants.AUtoKCAL
            evidence['legacy_zero_zpe_diagnostic'] = comparison(
                old_depth, 'fragment_E_plus_ZPE_minus_endpoint_E_with_missing_ZPE_set_to_zero')
        return dict(unknown('L1 endpoint/fragment ZPE unavailable; electronic and legacy depths are diagnostic only.'),
                    **evidence)
    depth = electronic + (sum(p['zpe_hartree'] for p in fragments)
                          - endpoint['zpe_hartree'])*constants.AUtoKCAL
    return dict(comparison(depth, 'matched_E_plus_ZPE'), **evidence)


def reachability(workers, initial):
    """Alternate paths and reversible single-product edges; never split fragments.

    The possible graph keeps unknown energies. A disconnected component that
    was already disconnected in the source is not attributed to this screen.
    """
    names = {w['worker'] for w in workers}
    graphs = {key: defaultdict(set) for key in ('before', 'possible', 'confirmed')}
    heights = {w['worker']: w['well']['decision'] for w in workers}
    for worker in workers:
        a = worker['worker']
        for channel in worker['channels']:
            targets = channel['product_ids']
            if not channel['accepted_l2'] or channel['operator_rejected'] or len(targets) != 1:
                continue
            b = targets[0]
            if b not in names:
                continue
            for key in graphs:
                statuses = [heights[a], heights[b], channel['direct_decision']]
                if key == 'possible' and 'exclude' in statuses:
                    continue
                if key == 'confirmed' and any(s != 'retain' for s in statuses):
                    continue
                graphs[key][a].add(b)
                graphs[key][b].add(a)
    def visit(graph):
        found, pending = {initial}, [initial]
        while pending:
            for name in graph[pending.pop()] - found:
                found.add(name)
                pending.append(name)
        return found
    before, possible, confirmed = (visit(graphs[k]) for k in ('before', 'possible', 'confirmed'))
    for worker in workers:
        name = worker['worker']
        if name == initial:
            status = 'initial'
        elif heights[name] == 'exclude':
            status = 'above_limit'
        elif name in confirmed:
            status = 'reachable_with_retained_edges'
        elif name in possible:
            status = 'reachable_with_unknown_edges'
        elif name in before:
            status = 'disconnected_by_screen_in_saved_graph'
        else:
            status = 'not_connected_in_saved_graph'
        worker['reachability'] = status
        # Disconnection is a candidate only: an unfinished search may add paths.
        for channel in worker['channels']:
            channel['decision'] = combined(worker['well'], dict(decision=channel['direct_decision']))
            channel['cascade_candidate'] = status == 'disconnected_by_screen_in_saved_graph'
            vrc = channel['vrc']
            if vrc['kind'] == 'none':
                vrc['decision'] = 'not_applicable'
            else:
                vrc['decision'] = combined(dict(decision=channel['decision']),
                    channel['vdw'] if vrc['kind'] == 'vdw_exit' else None)
            vrc['cascade_candidate'] = channel['cascade_candidate']


def work_inventory(worker, channels):
    """Current graph ledgers only, excluding retry archives and scratch trees."""
    from kinbot.reaction_decisions import owns
    graphs, issues = [], []
    paths = list(worker.glob('anl_automatic/states/*/composite/state.json'))
    for folder in ('molpro_scan', 'molpro', 'molpro_pivots'):
        paths.extend(worker.glob('vrctst/'+folder+'/*/state.json'))
    for path in sorted(set(paths)):
        try:
            state = read(path)
            counts = Counter(t.get('status', 'unknown') for t in state.get('tasks', {}).values())
            owners = [c['reaction'] for c in channels if owns(path.parent.name, c['reaction'])]
            graphs.append(dict(path=str(path.relative_to(worker)), task_statuses=dict(counts),
                               dedicated_vrc_reactions=owners))
        except (OSError, ValueError, TypeError, AttributeError) as error:
            issues.append(dict(path=str(path), error=str(error)))
    return dict(graphs=graphs, issues=issues,
                note='Ledger counts only; complete artifacts still need seed validation. Shared tasks are not assigned to individual reactions.')


def review_inventory(worker, name, configs):
    path = worker/'vrctst/orbital_review'/name/'request.json'
    request = read(path)
    matches = []
    for source, config in configs:
        model = config.get('vrc_tst_surface_models', {}).get(worker.name, {})
        for channel in model.get('channels', []):
            if channel.get('reaction') == name or (channel.get('reaction') is None
                    and request.get('bond') and sorted(channel.get('bond', [])) == sorted(request['bond'])):
                matches.append(dict(source=source, bond=channel['bond'],
                                    reaction=channel.get('reaction'), reviewed=model.get('reviewed')))
    return dict(request=str(path) if request else None, request_status=request.get('status', 'not_recorded'),
                configured_models=matches,
                gate='Existing models/receipts are inventory, not a new orbital approval; retained channels need seed-contract validation.')


def review_plans(plans, case, workers):
    lookup = {(w['worker'], c['reaction']): c for w in workers for c in w['channels']}
    result = []
    for path in plans:
        path = Path(path).resolve(strict=True)
        for item in read(path).get('screens', []):
            if item.get('case') != case.name:
                continue
            recipe = (path.parent/item['recipe']).resolve(strict=True)
            match = read(recipe).get('request_match', {})
            channel = lookup.get((match.get('parent'), match.get('reaction')))
            result.append(dict(plan=str(path), screen=item['name'], recipe=str(recipe),
                parent=match.get('parent'), reaction=match.get('reaction'),
                vrc_decision=channel['vrc']['decision'] if channel else 'unknown',
                cascade_candidate=channel['cascade_candidate'] if channel else False,
                reason=None if channel else 'No matching reaction in the saved audit inventory; do not drop this review.'))
    return result


def report(case, limit=150., vdw_depth=.5, plans=()):
    from kinbot.parameters import Parameters
    from kinbot.species_routing import matches_name, connectivity_name
    l1.validate(limit)
    if isinstance(vdw_depth, bool) or not math.isfinite(vdw_depth) or vdw_depth < 0:
        raise ValueError('vdW depth must be finite and nonnegative.')
    case = Path(case).resolve(strict=True)
    output = dict(schema=2, case=str(case), captured_utc=datetime.now(timezone.utc).isoformat(),
        limit_kcal_mol=limit, vdw_depth_kcal_mol=vdw_depth, workers=[], issues=[],
        limitations=[
            'Selection uses L1 electronic energy plus ZPE, relative to the initial L1 well; no L1/L2 mixing.',
            'Electronic-only and missing-ZPE-zero depth comparisons are diagnostic only and never exclude channels.',
            'Unknown evidence is retained for review. L1 filtering is approximate, not a kinetic flux criterion.',
            'A shallow complex loses its separate exit VRC channel, not automatically its preceding tight-TS reaction.',
            'Reachability uses reversible single-product edges in the saved accepted graph. Disconnection is a candidate, not proof for an unfinished exploration.',
            'Work counts are current graph ledgers, not CPU savings or validated restartable tasks. Shared CBH/fragment work can still be needed.',
            'No run files, database rows, decisions, controllers, or jobs are changed.'])
    try:
        ref = l1.reference(case)
    except (OSError, ValueError, sqlite3.Error) as error:
        ref = None
        output['issues'].append(dict(reference=str(error)))
    output['reference'] = ref
    case_config = read(case/'input.json')
    for worker in sorted(case.iterdir()):
        inp = worker/(worker.name+'.json')
        if not worker.is_dir() or not inp.is_file():
            continue
        try:
            par = Parameters(str(inp)).par
            overrides = read(worker/'.kinbot_search_options.json')
            par.update(overrides)
            profile = l1.theory(par)
        except (OSError, ValueError) as error:
            output['issues'].append(dict(worker=worker.name, error=str(error)))
            continue
        cache = {}
        def points(names):
            for name in names:
                if name not in cache:
                    cache[name] = l1.read_points(worker, [name], profile, geometry=True)[0]
            return [cache[n] for n in names]
        def measure(names, comparable=True):
            try:
                if ref is None:
                    return unknown('Initial L1 reference unavailable.')
                return l1.evaluate(ref, profile, points(names), limit,
                                   comparable=comparable, charge=par['charge'])
            except (OSError, ValueError, sqlite3.Error) as error:
                return unknown(error)
        well = measure([worker.name+'_well'])
        try:
            parent = species(points([worker.name+'_well'])[0], par)
        except (OSError, ValueError, sqlite3.Error) as error:
            parent = None
        entries = search_entries(worker)
        snapshot = read(worker/'accepted_l2_network.json')
        reactions = {r['name']: r for r in snapshot.get('network', {}).get('reactions', [])}
        requests = {p.parent.name for p in (worker/'vrctst/orbital_review').glob('*/request.json')}
        decisions = read(worker/'.kinbot_reaction_decisions.json').get('reactions', {})
        channels = []
        for name in sorted(set(entries) | set(reactions) | requests):
            saved, entry = reactions.get(name, {}), entries.get(name, {})
            family = saved.get('family', entry.get('signature', {}).get('family', 'unknown'))
            hom = family in ('hom_sci', 'hs') or saved.get('barrierless', False)
            comparable = l1.comparable_family(par, profile, family)
            ts = None if hom else measure([name], comparable)
            product_ids = saved.get('products', [])
            endpoint, endpoint_points = unknown('No verified product-side L1 endpoint.'), []
            vdw = dict(decision='not_applicable')
            identity_issue = None
            if not hom and (saved or entry.get('status') in (-1, 2, 3, 4, 5, 6)):
                try:
                    if parent is None:
                        raise ValueError('No valid parent L1 geometry.')
                    point, inferred = product_side(worker, name, par, profile, parent)
                    if product_ids and Counter(inferred) != Counter(product_ids):
                        raise ValueError('L1 endpoint fragment identities differ from accepted L2 products.')
                    product_ids = inferred
                    endpoint_points = [point]
                    endpoint = (l1.evaluate(ref, profile, endpoint_points, limit,
                                           comparable=comparable, charge=par['charge'])
                                if ref else unknown('Initial L1 reference unavailable.'))
                except (OSError, ValueError, sqlite3.Error) as error:
                    identity_issue = str(error)
                    endpoint = unknown(error)
            isolated = measure([p+'_well' for p in product_ids]) if product_ids else unknown('No verified product mapping.')
            if len(product_ids) == 2 and not hom:
                try:
                    if not endpoint_points or not comparable:
                        raise ValueError(identity_issue or 'Unverified L1 complex or theory.')
                    vdw = depth_decision(endpoint_points[0], points([p+'_well' for p in product_ids]), vdw_depth)
                    historical = vdw.get('legacy_zero_zpe_diagnostic')
                    if historical is None and vdw['decision'] != 'unknown':
                        historical = vdw
                    if historical is not None:
                        old_limit = par.get('vdW_detection')
                        if l1.finite_number(old_limit):
                            vdw['threshold_only_replay'] = dict(
                                old_threshold_kcal_mol=old_limit,
                                old_decision='retain' if historical['depth_kcal_mol'] > old_limit else 'exclude',
                                requested_decision=historical['decision'],
                                basis=historical['basis'], diagnostic_only=True)
                except (OSError, ValueError, sqlite3.Error) as error:
                    vdw = unknown(error)
            if identity_issue:
                products = unknown(identity_issue)
            elif hom or len(product_ids) != 2 or vdw['decision'] == 'exclude':
                products = isolated
            elif vdw['decision'] == 'retain':
                products = endpoint
            else:
                # Until depth is known, a lower complex may still be retained.
                products = unknown('Cannot choose complex versus separated-product energy without L1 depth.')
            configured = []
            for key in ('vrc_tst_scan', 'vrc_tst_noscan'):
                mapping = par.get(key, {})
                configured.extend(mapping.get(worker.name, mapping.get(connectivity_name(worker.name), [])))
            existing_vrc = bool(name in requests or matches_name(name, configured) or
                                (len(product_ids) == 2 and (hom or saved.get('product_complex'))))
            potential = bool(len(product_ids) == 2 and not hom and vdw['decision'] == 'retain')
            kind = 'homolytic' if hom and (existing_vrc or len(product_ids) == 2) else 'vdw_exit' if existing_vrc or potential else 'none'
            operator = decisions.get(name, {}).get('decision') == 'reject'
            direct = 'exclude' if operator else combined(ts, products)
            channels.append(dict(reaction=name, family=family, accepted_l2=name in reactions,
                search_status=entry.get('status'), operator_rejected=operator, barrier=ts,
                product_ids=product_ids, products=products, isolated_products=isolated,
                product_complex=endpoint, accepted_complex=saved.get('product_complex'), vdw=vdw,
                direct_decision=direct, vrc=dict(kind=kind, existing=existing_vrc,
                    review=review_inventory(worker, name, [(str(inp), par), (str(case/'input.json'), case_config)]))))
        output['workers'].append(dict(worker=worker.name, well=well, channels=channels,
            configured_vdw_detection=par.get('vdW_detection'), snapshot_errors=snapshot.get('errors', []),
            work=work_inventory(worker, channels)))
    if ref:
        reachability(output['workers'], ref['worker'])
    else:
        reachability(output['workers'], '')
    output['review_plans'] = review_plans(plans, case, output['workers'])
    output['counts'] = counts(output['workers'])
    return output


def counts(workers):
    channels = [c for w in workers for c in w['channels']]
    accepted = [c for c in channels if c['accepted_l2'] and not c['operator_rejected']]
    vrc = [c for c in channels if c['vrc']['kind'] != 'none']
    completed = Counter()
    by_selection = defaultdict(Counter)
    for w in workers:
        by_name = {c['reaction']: c for c in w['channels']}
        for g in w['work']['graphs']:
            completed.update(g['task_statuses'])
            owners = [by_name[n] for n in g['dedicated_vrc_reactions']]
            category = 'retained_or_shared_or_unknown'
            if w['well']['decision'] == 'exclude':
                category = 'under_above_limit_worker'
            elif owners and all(c['vrc']['decision'] == 'exclude' for c in owners):
                category = 'dedicated_to_excluded_vrc'
            elif w['reachability'] == 'disconnected_by_screen_in_saved_graph':
                category = 'disconnected_worker_candidate'
            by_selection[category].update(g['task_statuses'])
    return dict(workers=len(workers), well_energies=dict(Counter(w['well']['decision'] for w in workers)),
        reachability=dict(Counter(w['reachability'] for w in workers)),
        all_searches=len(channels), accepted_reactions=len(accepted),
        accepted_reaction_decisions=dict(Counter(c['decision'] for c in accepted)),
        accepted_reaction_cascade_candidates=sum(c['cascade_candidate'] and c['decision'] != 'exclude' for c in accepted),
        vrc_decisions=dict(Counter(c['vrc']['decision'] for c in vrc)),
        existing_vrc_decisions=dict(Counter(c['vrc']['decision'] for c in vrc if c['vrc']['existing'])),
        vdw_diagnostic_decisions={key: dict(Counter(c['vdw'][key]['decision'] for c in channels
            if key in c['vdw'])) for key in ('electronic_only', 'legacy_zero_zpe_diagnostic')},
        vdw_threshold_only_replay=dict(Counter(
            c['vdw']['threshold_only_replay']['old_decision']+'->'+c['vdw']['threshold_only_replay']['requested_decision']
            for c in channels if 'threshold_only_replay' in c['vdw'])),
        vrc_cascade_candidates=sum(c['cascade_candidate'] and c['vrc']['decision'] != 'exclude' for c in vrc),
        vrc_review_status_kept_or_unknown=dict(Counter(c['vrc']['review']['request_status'] for c in vrc if c['vrc']['decision'] != 'exclude')),
        vdw_depth_decisions=dict(Counter(c['vdw']['decision'] for c in channels if c['vdw']['decision'] != 'not_applicable')),
        recorded_task_statuses=dict(completed),
        recorded_task_statuses_by_selection={k: dict(v) for k,v in sorted(by_selection.items())})


def summarize(report):
    print('CASE', report['case'])
    print('COUNTS', json.dumps(report['counts'], sort_keys=True))
    for issue in report['issues']:
        print('ISSUE', json.dumps(issue, sort_keys=True))
    for w in report['workers']:
        print('WELL', w['worker'], w['well']['decision'], w['well'].get('relative_kcal_mol'), w['reachability'])
        for c in w['channels']:
            if c['accepted_l2'] or c['vrc']['kind'] != 'none':
                print('CHANNEL', json.dumps(dict(reaction=c['reaction'], accepted_l2=c['accepted_l2'],
                    decision=c['decision'], cascade_candidate=c['cascade_candidate'],
                    ts_kcal=(c['barrier'] or {}).get('relative_kcal_mol'),
                    products_kcal=c['products'].get('relative_kcal_mol'),
                    vdw_depth=c['vdw'].get('depth_kcal_mol'), vdw=c['vdw']['decision'],
                    electronic_depth=c['vdw'].get('electronic_only'),
                    legacy_zero_zpe_depth=c['vdw'].get('legacy_zero_zpe_diagnostic'),
                    threshold_only_replay=c['vdw'].get('threshold_only_replay'),
                    vrc=c['vrc']['decision'], review=c['vrc']['review']['request_status'],
                    reasons=[v['reason'] for v in [c['barrier'], c['products'], c['vdw']] if v and v.get('reason')]), sort_keys=True))
    for item in report.get('review_plans', []):
        print('PLANNED_REVIEW', json.dumps(item, sort_keys=True))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', help='PES case directory, or saved report with --summarize')
    parser.add_argument('--limit', type=float, default=150.)
    parser.add_argument('--vdw-depth', type=float, default=.5)
    parser.add_argument('--review-plan', action='append', default=[])
    parser.add_argument('--summarize', action='store_true')
    args = parser.parse_args(argv)
    if args.summarize:
        summarize(read(Path(args.path)))
    else:
        print(json.dumps(report(args.path, args.limit, args.vdw_depth, args.review_plan), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
