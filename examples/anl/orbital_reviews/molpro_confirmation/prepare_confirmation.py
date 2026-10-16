"""Prepare a fixed-geometry transfer check using ordinary KinBot recon inputs.

All chemistry, points, state weights and site controls are in the catalogue.
This adapter submits nothing. Run on the host containing the declared WFU files.
"""
import argparse
from copy import deepcopy
import json
from pathlib import Path

from kinbot.execution import atomic_json, file_sha256
from kinbot.molpro_reconnaissance import prepare, _load


def configurations(catalog):
    settings = catalog['settings']
    configs = []
    for case in catalog['cases']:
        electrons, orbitals = case['candidate_active_space']
        source = case['seed']
        closed = (source['active_electrons'] - electrons) // 2
        discarded = source['active_orbitals'] - closed - orbitals
        source_inactive = (case['electron_count'] - source['active_electrons']) // 2
        swaps = [[a - source_inactive, b - source_inactive]
                 for a, b in case['candidate_partitions'][f'{electrons},{orbitals}']]
        for basis in settings['bases']:
            for stage in ('caspt2', 'mrci'):
                if stage == 'mrci' and basis['basis'] != settings['mrci_basis']:
                    continue
                points = case['confirmation_points']
                if stage == 'mrci':
                    points = [p for p in points if p['name'] in settings['mrci_points']]
                title = case['id'] + '_' + basis['name'] + '_' + stage
                model = deepcopy(case['state_model'])
                model['correlation_options'] = deepcopy(settings['solver_options']['correlation'])
                model['correlation_options']['reference_policy'] = settings['reference_policies'][stage]
                cfg = {
                    'schema': 1, 'title': title,
                    'symbols': case['geometry']['symbols'],
                    'geometry_angstrom': points[0]['geometry_angstrom'],
                    'charge': 0, 'multiplicity': 2, 'bases': [basis],
                    'ppn': settings['processes'], 'memory_mw': settings['memory_mw_per_rank'],
                    'memory_mb': settings['memory_mb'], 'partition': settings['partition'],
                    'walltime': settings['walltime'],
                    'scf_convergence': settings['solver_options']['scf']['energy_convergence'],
                    'scf_options': {'solver': 'so_sci',
                        'max_iterations': settings['solver_options']['scf']['max_iterations']},
                    'casscf_options': {**settings['solver_options']['casscf'],
                        'solver': 'so_sci', 'exact_diagonalization': orbitals <= 8,
                        'quasi_newton': True},
                    'probe_active_space': True, 'probe_method': 'casscf',
                    'orbital_initialization': 'propagated',
                    'anchors': [{'name': p['name'], 'orbital_chain': title,
                                 'geometry_angstrom': p['geometry_angstrom']} for p in points],
                    'initial_records': {basis['name'] + '__' + points[0]['name']: {
                        'path': source['path'], 'sha256': source['sha256'],
                        'basis': source['basis'], 'project_orbitals': True}},
                    'probe_definition': {'active_electrons': electrons,
                        'active_orbitals': orbitals, 'state_model': model,
                        'chemical_note': case['chemical_note'],
                        'orbital_partition': {'source_active_electrons': source['active_electrons'],
                            'source_active_orbitals': source['active_orbitals'],
                            'close_lowest': closed, 'discard_highest': discarded,
                            'orbital_swaps': swaps}},
                    'correlation_checks': ([{'method': f'caspt2-f12({electrons},{orbitals})',
                        'shift': settings['caspt2_shift'], 'ipea': settings['caspt2_ipea']}]
                        if stage == 'caspt2' else [{'method': f'mrci-f12+q({electrons},{orbitals})'}]),
                    'confirmation_provenance': {'catalog_case': case['id'], 'stage': stage,
                        'source_orbital_review_sha256': source['reviewed_xml_sha256'],
                        'point_evidence': [p['source'] for p in points],
                        'basis_projection_is_only_a_starting_guess': True,
                        'production_approved': False}}
                configs.append(cfg)
    return configs


def build(catalog_path, root):
    catalog_path = Path(catalog_path).resolve()
    root = Path(root).resolve()
    if root.exists():
        raise FileExistsError('Choose a fresh confirmation directory; resume existing runs in place.')
    catalog = json.loads(catalog_path.read_text())
    # Check every externally reviewed file before creating a partly usable run.
    for case in catalog['cases']:
        source = case['seed']
        if file_sha256(source['path']) != source['sha256']:
            raise ValueError('Reviewed seed changed: ' + source['path'])
    root.mkdir(parents=True)
    from kinbot.run_format import ensure_current_run
    ensure_current_run(root, create=True)
    atomic_json(root / 'catalog.json', catalog)
    records = []
    for cfg in configurations(catalog):
        path = root / 'inputs' / (cfg['title'] + '.json')
        path.parent.mkdir(exist_ok=True)
        atomic_json(path, cfg)
        _load(path)
        manifest = prepare(path, root / cfg['title'])
        m = json.loads(manifest.read_text())
        m['slurm_excluded_nodes'] = catalog['execution']['excluded_nodes']
        m['startup_retries'] = catalog['execution']['startup_retries']
        m['solver_options'] = catalog['settings']['solver_options']
        atomic_json(manifest, m)
        records.append({'title': cfg['title'], 'configuration': str(path), 'configuration_sha256': file_sha256(path),
                        'manifest': str(manifest), 'jobs': len(m['jobs'])})
    atomic_json(root / 'confirmation.json', {'schema': 1, 'status': 'prepared_review_required',
        'catalog_sha256': file_sha256(root / 'catalog.json'), 'systems': records,
        'execution': {**catalog['execution'], 'memory_mb': catalog['settings']['memory_mb']},
        'production_approved': False})
    print(f'Prepared {sum(r["jobs"] for r in records)} points in {len(records)} chains: {root}')
    return root / 'confirmation.json'


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('--catalog', type=Path,
                   default=Path(__file__).with_name('geometry_and_seed_catalog.json'))
    args = p.parse_args()
    build(args.catalog, args.root)
