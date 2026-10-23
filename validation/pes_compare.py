"""Read-only PES comparison against an explicitly supplied reference network.

No scheduler, dispatcher, database writer, or recovery is called.
"""
import argparse
import ast
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sqlite3

from kinbot.run_pes_audit import input_label, read, database_record


def load_reference(reference):
    if isinstance(reference, (str, Path)):
        reference = json.loads(Path(reference).read_text())
    if not isinstance(reference, dict) or reference.get('schema') != 1:
        raise ValueError('Expected a schema-1 reference network.')
    nodes, edges = reference['nodes'], reference['connections']
    for name, row in nodes.items():
        endpoint(row['smiles'])
        if row.get('kind') == 'complex' and row.get('net_endpoint') not in nodes:
            raise ValueError(f'{name}: complex requires a declared net endpoint.')
    for name, row in edges.items():
        if (len(row['endpoints']) != 2 or any(n not in nodes for n in row['endpoints'])
                or row['kind'] not in ('association', 'barrierless', 'tight')):
            raise ValueError(f'{name}: invalid reference connection.')
    return reference


def connectivity(smiles):
    """Canonical explicit-H graph, retaining charge but ignoring bond order."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Cannot read molecular identity: {smiles}')
    mol = Chem.AddHs(mol)
    Chem.RemoveStereochemistry(mol)
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in mol.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetNumRadicalElectrons(0)
        atom.SetNoImplicit(True)
    return tuple(sorted(Chem.MolToSmiles(mol, canonical=True,
        allHsExplicit=True, allBondsExplicit=True).split('.')))


def endpoint(smiles):
    return tuple(sorted(s for label in smiles for s in connectivity(label)))


def _tail(path, errors):
    try:
        if not path.exists():
            return ''
        with path.open('rb') as f:
            size = path.stat().st_size
            f.seek(max(0, size-8*1024*1024))
            return f.read().decode(errors='replace')
    except OSError as error:
        errors.append(dict(file=str(path), error=str(error)))
        return ''


def _cached_labels(worker, names, errors):
    """Read existing product geometries without creating SQLite/ASE locks."""
    database = worker/'kinbot.db'
    if not database.is_file() or not names:
        return {}
    from ase.db.sqlite import SQLite3Database
    class ReadOnly(SQLite3Database):
        def _connect(self):
            return sqlite3.connect(database.resolve().as_uri()+'?mode=ro',
                                  uri=True, timeout=3)
    db = ReadOnly(str(database), use_lock_file=False)
    found = {}
    for name in sorted(names):
        try:
            for suffix in ('_well_high', '_well'):
                rows = list(db.select(name=name+suffix, sort='-id', limit=1))
                if not rows:
                    continue
                row = rows[0]
                # Never retrieve an older normal row after a newer failure.
                if row.data.get('status') == 'normal':
                    atoms = row.toatoms()
                    parameters = dict(row.get('calculator_parameters') or {})
                    structure = [v for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.positions)
                                 for v in [symbol, *map(float, pos)]]
                    label, _ = input_label(dict(structure=structure,
                        charge=parameters.get('charge', 0), mult=parameters.get('mult', 2)))
                    found[name] = dict(smiles=label, source=str(database), source_job=name+suffix,
                                       source_row_id=row.id, identity='geometry_connectivity_only')
                    found[name]['record'] = database_record(database, name+suffix,
                        dict(charge=parameters.get('charge', 0), mult=parameters.get('mult', 1)), label)
                    if found[name]['record'] is not None:
                        found[name]['record']['cache_level'] = 'L2' if suffix == '_well_high' else 'L1'
                break
        except Exception as error:
            errors.append(dict(file=str(database), product=name, error=str(error)))
    return found


def snapshot(root, *, reference):
    reference = load_reference(reference)
    case = Path(root).resolve(strict=True)
    if not (case/'input.json').is_file():
        raise ValueError('Expected a PES case directory containing input.json.')
    errors, labels, workers, rows, complexes = [], {}, {}, [], []
    cache_evidence, native_records = {}, {}
    phase = read(case/'pes_exploration.json', errors)
    selection = read(case/'pes_l2_selection.json', errors)
    selection_current = (selection.get('status') == 'selected' and
        all(selection.get(k) == phase.get(k) for k in ('token', 'epoch')))
    for w in sorted(case.glob('[0-9]*')):
        inp = w/(w.name+'.json')
        if not inp.is_file():
            continue
        params = read(inp, errors)
        try:
            label, _ = input_label(params)
            labels[w.name] = label
        except Exception as error:
            errors.append(dict(file=str(inp), error=str(error)))
        workers[w.name] = phase.get('workers', {}).get(w.name, {}).get('status', 'not_reported')
        if (w/'kinbot.db').is_file():
            job = w.name+('_well_high' if params.get('high_level', 1) else '_well')
            try:
                native_records[f'{w.name}/well'] = dict(
                    source=str(w/'kinbot.db'),
                    status='latest_native_cache; acceptance_not_revalidated',
                    record=database_record(w/'kinbot.db', job, params, labels.get(w.name)))
            except Exception as error:
                errors.append(dict(file=str(w/'kinbot.db'), job=job, error=str(error)))
        network, records = {}, {}
        for p in (w/params.get('anl_automatic_directory', 'anl_automatic')/'network_request.json',
                  w/'accepted_l2_network.json'):
            saved = read(p, errors).get('network', {})
            labels.update(saved.get('route_to_smiles', {}))
            records.update(saved.get('transition_states', {}))
            for category in ('species', 'product_complexes'):
                for key, value in saved.get(category, {}).items():
                    native_records[f'{w.name}/{category}/{key}'] = dict(source=str(p), record=value,
                        status='saved_record; current_acceptance_separate')
            network.update({r['name']: r for r in saved.get('reactions', [])})
            for name, value in saved.get('product_complexes', {}).items():
                if value.get('record', {}).get('smiles'):
                    labels[name] = value['record']['smiles']
        # Current monitor acceptance wins over older snapshots/plans.
        monitor = w/'kinbot_monitor.out'
        text = _tail(monitor, errors)
        hints = {}
        for line in _tail(w/'kinbot.log', errors).splitlines():
            found = re.search(r'reaction (\S+) leads to products (\[[^\]]*\])', line)
            if found:
                try:
                    hints[found[1]] = list(map(str, ast.literal_eval(found[2])))
                except (ValueError, SyntaxError):
                    continue
        if not monitor.is_file():
            errors.append(dict(file=str(monitor), error='Current search acceptance unavailable'))
        selected = selection.get('workers', {}).get(w.name, {})
        for line in text.splitlines():
            a = line.split()
            if len(a) < 3:
                continue
            try:
                code = int(a[0])
            except ValueError:
                continue
            name = a[2]
            old = network.get(name, {})
            products = a[3:] if code == -1 else hints.get(name, old.get('products', []))
            row = dict(name=name, worker=w.name, reactant=w.name, products=products,
                status='accepted' if code == -1 else 'rejected' if code == -999 else 'in_progress',
                step=a[1], endpoint_evidence='current_monitor' if code == -1 else 'candidate_only',
                barrierless='_hom_sci_' in name,
                selection='not_yet_selected' if not selection_current else
                    'retained' if name in selected.get('selected', []) else
                    'omitted_by_selection' if name in selected.get('omitted', {}) else 'not_yet_selected')
            if (code == -1 and old.get('product_complex')
                    and sorted(old.get('products', [])) == sorted(products)):
                row['product_complex'] = old['product_complex']
                complexes.append(dict(route=old['product_complex'], products=products, reaction=name))
            if code == -1 and not row['barrierless']:
                rec = records.get(old.get('transition_state')) if sorted(old.get('products', [])) == sorted(products) else None
                origin = 'saved_network; current_endpoints_match'
                if rec is None and (w/'kinbot.db').is_file():
                    try:
                        job = name+('_high' if params.get('high_level', 1) else '')
                        rec = database_record(w/'kinbot.db', job, params, labels.get(w.name))
                        origin = 'latest_native_cache; acceptance_not_revalidated'
                    except Exception as error:
                        errors.append(dict(file=str(w/'kinbot.db'), job=name, error=str(error)))
                row['native_record_key'] = f'{w.name}/transition_states/{name}'
                native_records[row['native_record_key']] = dict(source=origin, record=rec)
            rows.append(row)
        if text != _tail(monitor, errors):
            errors.append(dict(file=str(monitor), error='Search status changed during read; repeat check'))
    # A worker publishes its accepted network at the end of exploration.
    # Earlier monitor rows can still be identified from existing native caches.
    for worker in workers:
        missing = {p for row in rows if row['worker'] == worker for p in row['products'] if p not in labels}
        products = {p for row in rows if row['worker'] == worker and row['status'] == 'accepted'
                    for p in row['products']}
        cached = _cached_labels(case/worker, missing | products, errors)
        labels.update({name: row['smiles'] for name, row in cached.items() if name not in labels})
        cache_evidence.update(cached)
        for name, row in cached.items():
            native_records[f'{worker}/product_cache/{name}'] = dict(
                source=row['source'], record=row['record'], status='latest_native_cache; acceptance_not_revalidated')
    keys = {}
    for name, label in labels.items():
        try:
            keys[name] = connectivity(label)
        except ValueError as error:
            errors.append(dict(file=name, error=str(error)))
    def route_key(name):
        if name in keys:
            return keys[name]
        candidates = {value for route, value in keys.items() if route.split('-s')[0] == name}
        return next(iter(candidates)) if len(candidates) == 1 else None
    def end_key(names):
        parts = [route_key(n) for n in names]
        return tuple(sorted(x for part in parts for x in part)) if parts and all(parts) else None
    ref_nodes, ref_edges = reference['nodes'], reference['connections']
    ref_keys = {name: endpoint(v['smiles']) for name, v in ref_nodes.items()}
    complex_nodes = {n for n, v in ref_nodes.items() if v.get('kind') == 'complex'}
    def net_key(name):
        return ref_keys[ref_nodes[name]['net_endpoint']] if name in complex_nodes else ref_keys[name]
    complex_matches = {name: [c for c in complexes if end_key(c['products']) == net_key(name)]
                       for name in complex_nodes}
    accepted = [r for r in rows if r['status'] == 'accepted']
    observed = {end_key([r['reactant']]) for r in accepted} | {end_key(r['products']) for r in accepted}
    pending = {end_key(r['products']) for r in rows if r['status'] == 'in_progress'}
    if phase.get('initial_worker'):
        observed.add(route_key(str(phase['initial_worker'])))
    nodes = []
    for name, definition in ref_nodes.items():
        key = ref_keys[name]
        matches = [route for route, value in keys.items() if value == key]
        matched_endpoints = sorted({tuple(sorted(ns)) for r in accepted
            for ns in ([r['reactant']], r['products']) if end_key(ns) == key})
        if name in complex_nodes:
            matches = [c['route'] for c in complex_matches[name]]
            status = 'complex_candidate; reference_orientation_unverified' if matches else 'not_found_in_accepted_complex_records'
            matched_endpoints = []
        else:
            status = 'accepted_connectivity' if key in observed else (
                'well_discovered; acceptance_pending' if any(m in workers for m in matches)
                else 'candidate_endpoint_in_progress' if key in pending
                else 'not_found_in_current_accepted_network')
        nodes.append(dict(reference=name, definition=definition, status=status,
                          routing_keys=matches, endpoint_routing_keys=matched_endpoints))
    evidence = {name: [] for name in ref_edges}
    extra, extra_pending, unknown = [], [], []
    for row in rows:
        left, right = end_key([row['reactant']]), end_key(row['products'])
        if left is None or right is None:
            unknown.append(row)
            continue
        matched = []
        for name, definition in ref_edges.items():
            if definition['kind'] == 'association':
                continue
            if {left, right} == {net_key(n) for n in definition['endpoints']}:
                matched.append(name)
                evidence[name].append(dict(row, match='endpoints_only; TS_identity_unverified',
                    path_type_matches=row['barrierless'] == (definition['kind'] == 'barrierless'),
                    complex_recorded=bool(row.get('product_complex'))))
        if not matched and row['status'] in ('accepted', 'in_progress'):
            target = extra if row['status'] == 'accepted' else extra_pending
            target.append(dict(row, reactant_smiles=labels.get(row['reactant']),
                               product_smiles=[labels.get(p, '?') for p in row['products']]))
    edges = []
    for name, definition in ref_edges.items():
        required = [n for n in definition['endpoints'] if n in complex_nodes]
        candidates = [c for n in required for c in complex_matches[n]]
        ev = evidence[name]
        if definition['kind'] == 'association':
            status = 'complex_connection_recorded; geometry_unverified' if candidates else 'not_found'
            ev = candidates
        else:
            status = ('accepted_endpoint_and_path_type_match' if any(
                r['status'] == 'accepted' and r['path_type_matches'] for r in ev)
                else 'accepted_endpoints; different_path_type' if any(r['status'] == 'accepted' for r in ev)
                else 'candidate_in_progress' if any(r['status'] == 'in_progress' for r in ev)
                else 'rejected_candidate_only' if ev else 'not_found_in_current_records')
        edges.append(dict(reference=name, definition=definition, status=status,
            required_complexes={n: bool(complex_matches[n]) for n in required}, evidence=ev))
    fragments = {k for name, value in ref_keys.items() if name not in complex_nodes for k in value}
    accepted_routes = {p for r in accepted for p in [r['reactant'], *r['products']]}
    additions = [dict(route=n, smiles=labels.get(n), status='accepted_connectivity')
                 for n in sorted(accepted_routes) if route_key(n) is not None
                 and route_key(n) not in ref_keys.values()
                 and route_key(n) not in {(k,) for k in fragments}]
    pending_wells = [dict(route=n, smiles=labels.get(n), status=workers[n])
                     for n in sorted(workers) if n not in accepted_routes
                     and route_key(n) is not None and route_key(n) not in ref_keys.values()]
    return dict(reference=reference, case=str(case), exploration=phase,
        captured_utc=datetime.now(timezone.utc).isoformat(), selection_current=selection_current,
        comparison='Element and H connectivity only; no TS, conformer, spin or energy equivalence asserted. '
                   'Complexes use explicitly declared net endpoints and require separate geometry validation. '
                   'Not found is a snapshot observation, not a claim that a channel cannot exist; '
                   'unmapped endpoints may conceal further matches.',
        reference_nodes=nodes, reference_connections=edges, additional_accepted_species=additions,
        additional_accepted_connections=extra, unresolved_endpoints=unknown,
        additional_pending_wells=pending_wells, additional_in_progress_connections=extra_pending,
        cache_identity_evidence=cache_evidence, native_records=native_records,
        current_searches=rows, route_to_smiles=labels, selection=selection,
        search_counts=dict(Counter(r['status'] for r in rows)), errors=errors)


def render(report):
    out = [report['reference'].get('title', 'PES reference comparison')+' | read-only',
           report['comparison'], 'Searches: '+str(report['search_counts'])]
    for row in report['reference_nodes']:
        routes = row['routing_keys'] or row['endpoint_routing_keys']
        out.append(f"{row['reference']} {row['definition']['label']}: {row['status']}; routes={routes}")
    for row in report['reference_connections']:
        out.append(f"{row['reference']} {' <-> '.join(row['definition']['endpoints'])}: {row['status']}")
        for node, found in row['required_complexes'].items():
            if not found:
                out.append(f'  {node} complex not recorded; a net endpoint match does not establish this full path.')
        for ev in row['evidence']:
            out.append('  '+ev.get('name', ev.get('reaction', '?'))+' '+ev.get('selection', ''))
    for key in ('additional_accepted_species', 'additional_accepted_connections',
                'additional_pending_wells', 'additional_in_progress_connections'):
        out.append(key+':')
        out.extend('  '+json.dumps(row) for row in report[key])
    out.append(f"Unmapped endpoints: {len(report['unresolved_endpoints'])}; read errors: {len(report['errors'])}; full evidence in --json")
    out.extend('  READ CHECK '+json.dumps(e) for e in report['errors'])
    return '\n'.join(out)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('case', type=Path)
    p.add_argument('--reference', required=True, type=Path)
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--json', action='store_true')
    mode.add_argument('--evidence', action='store_true', help='Readable comparison plus compact native energies, frequencies and coordinates.')
    args = p.parse_args()
    report = snapshot(args.case, reference=args.reference)
    print(json.dumps(report, indent=2) if args.json else render(report))
    if args.evidence:
        fields = ('smiles', 'source_job', 'source_row_id', 'cache_level', 'energy_hartree',
                  'harmonic_zpe_hartree', 'frequencies_cm_inverse', 'l2_theory')
        for key, item in sorted(report['native_records'].items()):
            record = item.get('record') or {}
            record = record.get('record', record)
            native = {k: record[k] for k in fields if k in record}
            mol = record.get('molecule', {})
            native['molecule'] = {k: mol[k] for k in ('symbols', 'positions', 'charge', 'multiplicity') if k in mol}
            print('NATIVE', key, json.dumps(dict(source=item['source'], **native), separators=(',', ':')))
        for row in report['unresolved_endpoints']:
            print('UNMAPPED', json.dumps(row, separators=(',', ':')))


if __name__ == '__main__':
    main()
