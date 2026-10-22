"""Read-only comparison of the CH3SS network with D4CP04873C, Fig. 6.

This compares element/H connectivity and unordered endpoints, not energies,
electronic terms, stereochemistry, or proof of identical transition states.
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

from kinbot.run_pes_audit import input_label, read


SOURCE = 'https://pubs.rsc.org/en/content/articlehtml/2026/cp/d4cp04873c'
# Connectivity from section 3.3.2; labels/relative energies from the supplied
# Fig. 6. Radical/bond-order conventions are intentionally not identity tests.
NODES = {
    'W1': ('CH2S...SH complex', ['C=S', '[SH]'], -4.77),
    'W2': ('CH2S2H', ['[CH2]SS'], -18.46),
    'W3': ('HSCH2S', ['SC[S]'], -22.82),
    'W4': ('CH3S2', ['CS[S]'], -43.02),
    'P1': ('CH2S + SH', ['C=S', '[SH]'], 0.),
    'P2': ('CHS + H2S', ['[CH]=S', 'S'], 1.49),
    'P3': ('SCHSH + H', ['S=CS', '[H]'], 16.60),
    'P4': ('CH3 + S2', ['[CH3]', '[S][S]'], 3.57),
}
EDGES = {
    'entrance': ('P1', 'W1', None), 'TS1': ('W1', 'P2', 7.48),
    'TS2': ('W2', 'W1', -4.22), 'TS3': ('W3', 'W1', -.17),
    'TS4': ('W4', 'W2', 7.38), 'TS5': ('W3', 'P3', 21.62),
    'TS6': ('W2', 'W3', 18.20), 'scission': ('W4', 'P4', None),
}


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
                break
        except Exception as error:
            errors.append(dict(file=str(database), product=name, error=str(error)))
    return found


def snapshot(root):
    root = Path(root).resolve(strict=True)
    case = root/'ch3ss' if (root/'ch3ss/input.json').is_file() else root
    if not (case/'input.json').is_file():
        raise ValueError('Expected the production root or its CH3SS case directory.')
    errors, labels, workers, rows, complexes = [], {}, {}, [], []
    cache_evidence = {}
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
        network = {}
        for p in (w/params.get('anl_automatic_directory', 'anl_automatic')/'network_request.json',
                  w/'accepted_l2_network.json'):
            saved = read(p, errors).get('network', {})
            labels.update(saved.get('route_to_smiles', {}))
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
            rows.append(row)
        if text != _tail(monitor, errors):
            errors.append(dict(file=str(monitor), error='Search status changed during read; repeat check'))
    # A worker publishes its accepted network at the end of exploration.
    # Earlier monitor rows can still be identified from existing native caches.
    for worker in workers:
        missing = {p for row in rows if row['worker'] == worker for p in row['products'] if p not in labels}
        cached = _cached_labels(case/worker, missing, errors)
        labels.update({name: row['smiles'] for name, row in cached.items()})
        cache_evidence.update(cached)
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
    reference = {name: endpoint(v[1]) for name, v in NODES.items() if name != 'W1'}
    reverse = {value: name for name, value in reference.items()}
    accepted = [r for r in rows if r['status'] == 'accepted']
    observed = {end_key([r['reactant']]) for r in accepted} | {end_key(r['products']) for r in accepted}
    pending = {end_key(r['products']) for r in rows if r['status'] == 'in_progress'}
    if phase.get('initial_worker'):
        observed.add(route_key(str(phase['initial_worker'])))
    complex_matches = [c for c in complexes if end_key(c['products']) == reference['P1']]
    nodes = []
    for name, (label, smi, energy) in NODES.items():
        matches = [route for route, key in keys.items() if key == endpoint(smi)]
        if name == 'W1':
            status = 'complex_candidate; paper_orientation_unverified' if complex_matches else 'not_found_in_accepted_complex_records'
            matches = [c['route'] for c in complex_matches]
        else:
            status = 'accepted_connectivity' if reference[name] in observed else (
                'well_discovered; acceptance_pending' if any(m in workers for m in matches)
                else 'candidate_endpoint_in_progress' if reference[name] in pending
                else 'not_found_in_current_accepted_network')
        matched_endpoints = sorted({tuple(sorted(ns)) for r in accepted
            for ns in ([r['reactant']], r['products']) if end_key(ns) == endpoint(smi)}) if name != 'W1' else []
        nodes.append(dict(paper=name, label=label, paper_E0_kcal_mol=energy,
                          status=status, routing_keys=matches, endpoint_routing_keys=matched_endpoints))
    evidence = {name: [] for name in EDGES}
    extra, extra_pending, unknown = [], [], []
    for row in rows:
        left_key, right_key = end_key([row['reactant']]), end_key(row['products'])
        left, right = reverse.get(left_key), reverse.get(right_key)
        if left_key is None or right_key is None:
            unknown.append(row)
            continue
        matched = []
        for name, (a, b, energy) in EDGES.items():
            if name == 'entrance':
                continue
            pair = {a.replace('W1', 'P1'), b.replace('W1', 'P1')}
            if {left, right} == pair:
                matched.append(name)
                evidence[name].append(dict(row, match='endpoints_only; TS_identity_unverified',
                    path_type_matches=row['barrierless'] == (energy is None),
                    complex_recorded=bool(row.get('product_complex'))))
        if not matched and row['status'] in ('accepted', 'in_progress'):
            target = extra if row['status'] == 'accepted' else extra_pending
            target.append(dict(row, reactant_smiles=labels.get(row['reactant']),
                               product_smiles=[labels.get(p, '?') for p in row['products']]))
    edges = [dict(paper=name, endpoints=[a,b], paper_E0_kcal_mol=energy,
        status=('complex_connection_recorded; geometry_unverified' if complex_matches else 'not_found')
            if name == 'entrance' else ('accepted_endpoint_and_path_type_match' if any(
                r['status']=='accepted' and r['path_type_matches'] for r in evidence[name])
            else 'accepted_endpoints; different_path_type' if any(r['status']=='accepted' for r in evidence[name])
            else 'candidate_in_progress' if any(r['status']=='in_progress' for r in evidence[name])
            else 'rejected_candidate_only' if evidence[name] else 'not_found_in_current_records'),
        W1_required=(a == 'W1' or b == 'W1'), W1_complex_recorded=bool(complex_matches),
        evidence=complex_matches if name == 'entrance' else evidence[name])
        for name,(a,b,energy) in EDGES.items()]
    paper_fragments = {k for value in reference.values() for k in value}
    accepted_routes = {p for r in accepted for p in [r['reactant'], *r['products']]}
    additions = [dict(route=n, smiles=labels.get(n), status='accepted_connectivity')
                 for n in sorted(accepted_routes) if route_key(n) is not None
                 and route_key(n) not in reference.values()
                 and route_key(n) not in {(k,) for k in paper_fragments}]
    pending_wells = [dict(route=n, smiles=labels.get(n), status=workers[n])
                     for n in sorted(workers) if n not in accepted_routes
                     and route_key(n) is not None and route_key(n) not in reference.values()]
    return dict(source=SOURCE, figure=6, case=str(case), exploration=phase.get('phase'),
        captured_utc=datetime.now(timezone.utc).isoformat(), selection_current=selection_current,
        comparison='Element and H connectivity only; no TS, conformer, spin or energy equivalence asserted. '
                   'Edges involving W1 use net fragment endpoints; a separate complex record is required '
                   'and its paper orientation remains unverified. '
                   'Not found is a snapshot observation, not a claim that a channel cannot exist; '
                   'unmapped endpoints may conceal further matches.',
        paper_nodes=nodes, paper_connections=edges, additional_accepted_species=additions,
        additional_accepted_connections=extra, unresolved_endpoints=unknown,
        additional_pending_wells=pending_wells, additional_in_progress_connections=extra_pending,
        cache_identity_evidence=cache_evidence,
        search_counts=dict(Counter(r['status'] for r in rows)), errors=errors)


def render(report):
    out = ['CH3SS versus D4CP04873C Fig. 6 | read-only', report['comparison'],
           'Searches: '+str(report['search_counts'])]
    for row in report['paper_nodes']:
        routes = row['routing_keys'] or row['endpoint_routing_keys']
        out.append(f"{row['paper']} {row['label']}: {row['status']}; routes={routes}")
    for row in report['paper_connections']:
        out.append(f"{row['paper']} {' <-> '.join(row['endpoints'])}: {row['status']}")
        if row['W1_required'] and not row['W1_complex_recorded']:
            out.append('  W1 complex not recorded; a net endpoint match does not establish this full paper path.')
        for ev in row['evidence']:
            detail = '' if 'path_type_matches' not in ev else (
                ' tight-TS candidate' if not ev['barrierless'] else ' barrierless')
            out.append('  '+ev.get('name', ev.get('reaction', '?'))+' '+ev.get('selection', '')+detail)
    out.append('Additional accepted species:')
    out.extend(f"  {r['route']} {r['smiles']}" for r in report['additional_accepted_species'])
    out.append('Additional accepted connections (searches, not unique saddles):')
    out.extend(f"  {r['name']}: {r['reactant_smiles']} -> {' + '.join(r['product_smiles'])}; {r['selection']}"
               for r in report['additional_accepted_connections'])
    out.append('Additional wells still awaiting search acceptance:')
    out.extend(f"  {r['route']} {r['smiles']}; {r['status']}" for r in report['additional_pending_wells'])
    out.append('Additional in-progress endpoint candidates (not accepted discoveries):')
    out.extend(f"  {r['name']}: {r['reactant_smiles']} -> {' + '.join(r['product_smiles'])}"
               for r in report['additional_in_progress_connections'])
    counts = dict(Counter(r['status'] for r in report['unresolved_endpoints']))
    out.append(f"Unmapped endpoints by search status: {counts}; read errors: {len(report['errors'])}; full evidence in --json")
    for error in report['errors'][:10]:
        out.append(f"  READ CHECK {error['file']}: {error['error']}")
    return '\n'.join(out)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('--json', action='store_true')
    args = p.parse_args()
    report = snapshot(args.root)
    print(json.dumps(report, indent=2) if args.json else render(report))


if __name__ == '__main__':
    main()
