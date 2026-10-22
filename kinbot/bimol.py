"""Low-cost candidate generation and validation for bimolecular reactions.

Reaction families provide encounter geometries and formed/broken bond lists.
This module applies the target-free force bias, refines saddle candidates,
clusters equivalent candidates, validates unique representatives by L0 IRC,
and promotes validated structures to KinBot's normal L1 workflow.
"""

import copy
import hashlib
import json
import logging
import math
import shutil
import time
import traceback
from importlib import import_module
from pathlib import Path

import numpy as np
import rmsd
from ase import Atoms
from ase.constraints import ExternalForce
from ase.data import atomic_numbers, covalent_radii, vdw_radii
from ase.io import write
from ase.optimize import BFGS
from sella import IRC as SellaIRC
from sella import Sella

from kinbot.calculation import load_calculation_record, selected_calculation_job
from kinbot import constants, frequencies, geometry
from kinbot.stationary_pt import StationaryPoint


GOLDEN_ANGLE = math.pi * (3. - math.sqrt(5.))
logger = logging.getLogger('KinBot')
DEFAULTS = {
    'sequential': True,
    'pose_angles': [0., 20., 40., 60.],
    'pose_azimuths': 4,
    'encounter_extra_distance': 1.,
    'collision_scale': .7,
    'force_max': 3.,
    'force_levels': 31,
    'force_bisections': 5,
    'break_scale': 1.,
    'form_scale': 1.,
    'path_fmax': .1,
    'path_steps': 80,
    'max_candidates': 4,
    'candidate_rmsd': .08,
    'irc_fmax': .01,
    'irc_steps': 100,
    'endpoint_fmax': .03,
    'endpoint_steps': 300,
}


def same_stereopath(first, second):
    """Compare bimolecular reaction-site classes during candidate clustering.

    These labels describe the enumerated abstraction sites. The shared TS
    pathway metadata additionally protects configurations during refinement.
    """
    return (getattr(first, 'stereopath_id', None)
            == getattr(second, 'stereopath_id', None))


def search_job(name, attempt=0, sequential=False):
    """Return the L0 candidate-search job label for a reaction channel."""
    suffix = '_sequential' if sequential else ''
    suffix += f'_try{attempt:04d}' if attempt else ''
    return str(Path('bimol', name, f'{name}_bimol_guess{suffix}'))


def validation_job(directory):
    """Return the L0 IRC-validation job label for a unique candidate."""
    directory = Path(directory)
    return str(directory / f'{directory.name}_bimol_irc')


def validate_search_kwargs(kwargs):
    """Validate user overrides for the bimolecular search."""
    unknown = set(kwargs) - set(DEFAULTS)
    if unknown:
        raise ValueError(f'Unknown bimol_search_kwargs: {sorted(unknown)}')
    settings = DEFAULTS.copy()
    settings.update(kwargs)
    if not isinstance(settings['sequential'], bool):
        raise ValueError('bimol_search_kwargs.sequential must be true or false.')
    positive = [key for key in DEFAULTS if key not in {'pose_angles', 'sequential'}]
    if any(settings[key] <= 0 for key in positive):
        raise ValueError('Bimolecular search values must be positive.')
    if settings['force_levels'] < 2:
        raise ValueError('force_levels must be at least two.')
    if any(angle < 0 or angle > 180 for angle in settings['pose_angles']):
        raise ValueError('pose_angles must be between 0 and 180 degrees.')
    return settings


def _search_settings(overrides):
    overrides = overrides or {}
    settings = validate_search_kwargs({key: value
                                        for key, value in overrides.items()
                                        if key in DEFAULTS})
    settings.update({key: value for key, value in overrides.items()
                     if key not in DEFAULTS})
    return settings


def create_calculator(qc, module, calculator_class, kwargs, command, ppn,
                      model, task, device, compile_model=False):
    """Construct the configured ASE calculator inside a queued worker."""
    if qc == 'fc':
        from fairchem.core import FAIRChemCalculator
        from kinbot.fairchem_utils import load_predictor
        return FAIRChemCalculator(load_predictor(model, device, compile_model=compile_model),
                                  task_name=task)
    calculator_type = getattr(import_module(module), calculator_class)
    kwargs = kwargs.copy()
    if calculator_class == 'ORCA':
        from kinbot.ase_modules.calculators.orca import OrcaProfile
        kwargs['profile'] = OrcaProfile(command=kwargs['profile'])
    calculator = calculator_type(**kwargs)
    if qc == 'gauss':
        calculator.command = f'{command} < PREFIX.com > PREFIX.log'
    elif qc == 'qchem':
        calculator.command = (f'{command} -nt {ppn} -save '
                              'PREFIX.inp PREFIX.out PREFIX.sv')
    return calculator


def _json(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def _merge_status(record, updates):
    """Return whether any saved value actually changed."""
    lifecycle = record.get('lifecycle', {})
    if all(key in record and record[key] == value
           and key in lifecycle and lifecycle[key] == value
           for key, value in updates.items()):
        return False
    record.update(updates)
    record.setdefault('lifecycle', {}).update(updates)
    return True


def _update_index_batch(path, key, changes):
    """Read and write a shared index at most once for a group of updates."""
    path = Path(path)
    if not path.is_file():
        return
    data = json.loads(path.read_text())
    changed = False
    for item in data.get(key, []):
        updates = changes.get(item.get('id'))
        if updates is not None:
            changed = _merge_status(item, updates) or changed
    if changed:
        _json(path, data)


def _update_index(path, key, item_id, updates):
    _update_index_batch(path, key, {item_id: updates})


def update_bimol_status(reaction, *, index_updates=None, **updates):
    """Persist one promoted candidate's lifecycle in every L0 manifest."""
    if not hasattr(reaction, 'search_parent'):
        return
    if 'final_status' in updates:
        reaction.bimol_final_status = updates['final_status']
    candidate = reaction.search_candidate
    files = [Path(candidate['directory'], 'result.json')]
    candidate_id = getattr(reaction, 'l0_candidate_id', None)
    l0_ts_id = getattr(reaction, 'l0_ts_id', None)
    if candidate_id:
        files.append(Path('bimol', 'unique_candidates', candidate_id,
                          'metadata.json'))
    if l0_ts_id:
        files.append(Path('bimol', 'validated_l0', l0_ts_id,
                          'metadata.json'))
    for result_file in files:
        if not result_file.is_file():
            continue
        result = json.loads(result_file.read_text())
        if _merge_status(result, updates):
            _json(result_file, result)
    for key, item_id in (('unique_candidates', candidate_id),
                         ('validated_l0', l0_ts_id)):
        if not item_id:
            continue
        path = f'bimol/{key}/index.json'
        if index_updates is None:
            _update_index(path, key, item_id, updates)
        else:
            index_updates.setdefault((path, key), {}).setdefault(item_id, {}).update(updates)


def load_selected_species_properties(optimization, qc):
    """Synchronize an optimized species to its selected calculation."""
    job = selected_calculation_job(optimization)
    species = optimization.species
    rows = list(qc.db.select(name=job))
    row = rows[-1] if rows else None
    # Optimize has already prepared the selected Hessian and thermal modes.
    # Loading the same row again discards a native-backend Hessian and restores
    # the unprojected frequencies. Keep that completed model when its complete
    # calculation record is unchanged; a newer row must still be loaded.
    unchanged = (row is not None
                 and getattr(species, 'source_job', None) == row.name
                 and getattr(species, 'source_row_id', None) == row.id
                 and np.array_equal(species.geom, row.positions)
                 and row.data.get('energy') is not None
                 and species.energy == row.data['energy'] * constants.EVtoHARTREE
                 and species.zpe == row.data.get('zpe')
                 and np.array_equal(species.freq, row.data.get('frequencies')))
    if unchanged:
        optimization.selected_job = job
        return job

    load_calculation_record(species, qc, job)
    optimization._projection_source = None
    if (optimization.par['multi_conf_tst'] == 0
            and optimization.par['rotor_scan']
            and not optimization.just_high
            and not optimization.defer_hir):
        from kinbot.conformer_records import hessian_record
        from kinbot.hindered_rotors import recover_hir_model, use_harmonic_model
        if not len(species.hess):
            # Read existing native Gaussian/Q-Chem output, without submitting
            # an unpolled frequency job from this completed-result loader.
            record = hessian_record(qc, job, species.geom, species.atom)
            if record.get('hessian') is not None:
                species.hess = np.asarray(record['hessian']).tolist()
                species.hessian_source_job = job
                species.optical_hessian_reference = record['hessian_reference']
        if getattr(species, 'hir', None) is not None:
            recover_hir_model(species, qc, optimization.par)
        else:
            use_harmonic_model(species, optimization.par,
                               'The selected calculation has no HIR scans')
    for attr in ('kinbot_freqs', 'reduced_freqs'):
        setattr(species, attr, frequencies.thermochemical_frequencies(
            getattr(species, attr), species.wellorts,
            optimization.par.get('imagfreq_threshold', 50.)))
    optimization.selected_job = job
    return job


def build_bimolecular_ts(reaction, qc):
    """Build the L1 TS topology before final conformer selection and IRC."""
    species = reaction.species
    err, geom = qc.get_qc_geom(reaction.instance_name, species.natom)
    if err:
        raise ValueError(f'Cannot read L1 TS geometry for {reaction.instance_name}.')
    ts = StationaryPoint(
        reaction.instance_name, species.charge, species.mult,
        atom=species.atom, geom=geom, wellorts=1)
    from kinbot.stereo_identity import canonical_identity
    ts.optical_reference = canonical_identity(species)
    load_calculation_record(ts, qc, reaction.instance_name)
    data = reaction.get_bimol_search_data()
    union = np.array(species.bond, dtype=int)
    product = np.array(species.bond01, dtype=int)
    validation = getattr(reaction, 'l0_validation', {})
    observed = next((endpoint for endpoint in validation.get('endpoints', [])
                     if endpoint.get('state') != 'reactant'), None)
    if observed is not None and 'bonds' in observed:
        product.fill(0)
        for first, second in observed['bonds']:
            product[first, second] = product[second, first] = 1
            union[first, second] = union[second, first] = max(
                1, union[first, second])
    else:
        for first, second in data['formed_bonds']:
            union[first, second] = union[second, first] = max(
                1, union[first, second])
            product[first, second] = product[second, first] = 1
        for first, second in data['broken_bonds']:
            product[first, second] = product[second, first] = 0
    ts.bond = union
    ts.bonds = []
    for reactant_matrix in species.bonds or [species.bond]:
        matrix = np.array(reactant_matrix, dtype=int)
        for first in range(species.natom):
            for second in range(first + 1, species.natom):
                if product[first, second]:
                    matrix[first, second] = matrix[second, first] = max(
                        1, matrix[first, second])
        if not any(np.array_equal(matrix, existing)
                   for existing in ts.bonds):
            ts.bonds.append(matrix)
    ts.bond01 = (ts.bond > 0).astype(int)
    ts.reac_bond = product - np.asarray(species.bond01)
    ts.rads = [np.asarray(radical).copy()
               for radical in getattr(species, 'rads', [])]
    ts.find_cycle()
    ts.characterize(bond_mx=ts.bond)
    endpoints = _l0_endpoint_species(reaction)
    if not endpoints:
        raise ValueError(f'Missing configured L0 endpoints for {reaction.instance_name}.')
    from kinbot.reaction_path import prepare_stereopath
    prepare_stereopath(ts, *endpoints)
    # The final IRC, including after a rotor refinement, determines products.
    ts.bimol_reactant_only = True
    reaction.ts = ts
    return ts


def check_l1_ts_refinement(reaction, par):
    """Check the L0-to-L1 change before starting TS conformer calculations."""
    from kinbot.optimize import same_ts_geometry, valid_ts_frequencies
    from kinbot.reaction_path import prepare_stereopath
    from kinbot.stereo_identity import configured_geometry_allowed

    ts = reaction.ts
    reference = copy.copy(ts)
    reference.geom = np.asarray(reaction.search_candidate['geometry'])
    prepare_stereopath(reference, *ts.ts_endpoint_graphs)
    refined = StationaryPoint('refined', ts.charge, ts.mult,
                              atom=ts.atom, geom=ts.geom)
    refined.bond_mx()
    # Keep the TS union topology, including forming and breaking bonds, on
    # the L0 reference. Saved L0 candidates have no comparable normal mode,
    # so use the ordinary close-geometry arm without another Hessian job.
    if not same_ts_geometry(reference, refined):
        raise ValueError(f'L1 TS refinement changed L0 bond lengths or formed '
                         f'extra bonds for {reaction.instance_name}.')
    if not valid_ts_frequencies(ts.freq, par['imagfreq_threshold']):
        raise ValueError(f'L1 TS refinement has invalid imaginary frequencies '
                         f'for {reaction.instance_name}.')
    if not configured_geometry_allowed(reference, ts.geom, par.get('optical_population', 'specified')):
        raise ValueError(f'L1 TS refinement changed the configured stereochemistry '
                         f'for {reaction.instance_name}.')


def update_ts_topology_from_irc(reaction):
    """Set TS topology to the union of its observed final IRC endpoints."""
    reactant = reaction.irc_reactant
    product = reaction.irc_prod
    ts = reaction.ts
    reactant_matrices = reactant.bonds or [reactant.bond]
    product_matrices = product.bonds or [product.bond]
    ts.bonds = []
    for reactant_matrix in reactant_matrices:
        for product_matrix in product_matrices:
            matrix = np.maximum(reactant_matrix, product_matrix)
            if not any(np.array_equal(matrix, existing)
                       for existing in ts.bonds):
                ts.bonds.append(matrix)
    ts.bond = np.maximum(reactant.bond, product.bond)
    ts.bond01 = (ts.bond > 0).astype(int)
    ts.reac_bond = (np.asarray(product.bond01)
                    - np.asarray(reactant.bond01))
    ts.find_cycle()
    ts.characterize(bond_mx=ts.bond)
    from kinbot.reaction_path import prepare_stereopath
    prepare_stereopath(ts, reactant, product)
    reaction.actual_bond_changes = {
        'formed': [[first, second] for first in range(ts.natom)
                   for second in range(first + 1, ts.natom)
                   if ts.reac_bond[first][second] > 0],
        'broken': [[first, second] for first in range(ts.natom)
                   for second in range(first + 1, ts.natom)
                   if ts.reac_bond[first][second] < 0],
    }
    assign_observed_site(reaction)


def assign_observed_site(reaction):
    """Use observed reacting atoms for grouping, retaining the search target."""
    if not hasattr(reaction, 'get_bimol_search_data'):
        return
    def edges(pairs):
        return {tuple(sorted(pair)) for pair in pairs}
    actual = reaction.actual_bond_changes
    templates = reaction.species.reac_obj or [reaction]
    for template in templates:
        if not getattr(template, 'bimol_search', False):
            continue
        data = template.get_bimol_search_data()
        choices = [dict(instance=list(template.instance), **{
            key: data[key] for key in ('formed_bonds', 'broken_bonds')})]
        choices.extend(data.get('equivalent_bond_changes', []))
        for choice in choices:
            if (edges(actual['formed']) != edges(choice['formed_bonds'])
                    or edges(actual['broken']) != edges(choice['broken_bonds'])):
                continue
            if not matches_expected_product(reaction.species, reaction.irc_prod.geom,
                    choice['formed_bonds'], choice['broken_bonds']):
                continue
            observed = dict(instance=list(choice['instance']), family_name=template.family_name)
            for name in ('stereopath_id', 'site_degeneracy', 'site_members', 'site_relations'):
                if hasattr(template, name):
                    observed[name] = copy.deepcopy(getattr(template, name))
            reaction.observed_site = observed
            return
    # A different reaction has no demonstrated multiplicity of equivalent
    # sites. Count its observed route once, rather than borrowing the target's.
    atoms = sorted({atom for kind in actual.values() for pair in kind for atom in pair})
    if not atoms:
        atoms = list(range(reaction.ts.natom))  # A stereo-only reaction changes no bonds.
    members = {f'atom_{i}': [atom] for i, atom in enumerate(atoms)}
    reaction.observed_site = dict(instance=atoms, family_name='bimolecular',
        stereopath_id='observed:' + reaction.ts.stereopath_id,
        site_degeneracy=1, site_members=members,
        site_relations={role: 'singleton' for role in members})


def _l0_endpoint_identities(reaction, validation=None):
    """Configured evidence gates early selection, not final reaction identity."""
    if validation is None:
        validation = getattr(reaction, 'l0_validation', {})
    endpoints = validation.get('endpoints', [])
    if len(endpoints) != 2 or sum(p.get('state') == 'reactant' for p in endpoints) != 1:
        return None
    identities = [p.get('configured_identity', {}) for p in endpoints]
    if any(p.get('status') != 'assigned' for p in identities):
        return None
    return tuple(sorted((endpoint.get('state'), identity['id'])
                        for endpoint, identity in zip(endpoints, identities)))


def _l0_endpoint_species(reaction):
    """Rebuild the saved L0 endpoints in the common TS atom order."""
    endpoints = getattr(reaction, 'l0_validation', {}).get('endpoints', [])
    result = []
    for endpoint in sorted(endpoints, key=lambda p: p.get('state') != 'reactant'):
        identity = endpoint.get('configured_identity', {})
        if identity.get('status') != 'assigned' or 'geometry_angstrom' not in endpoint:
            return []
        point = StationaryPoint('l0_endpoint', identity['charge'], identity['multiplicity'],
                                atom=reaction.species.atom,
                                geom=np.asarray(endpoint['geometry_angstrom'], dtype=float))
        point.characterize()
        result.append(point)
    return result if len(result) == 2 else []


def final_irc_job_name(reaction):
    """Keep final IRC results specific to the selected TS geometry and level."""
    level = 'l2' if getattr(reaction, 'par', {}).get('high_level') else 'l1'
    source = getattr(reaction.ts, 'source_job', None)
    qc = getattr(reaction, 'qc', None)
    if source and qc is not None:
        level = qc.for_job(source).theory_level
    origin_hash = hashlib.sha1(np.round(
        reaction.ts.geom, 6).tobytes()).hexdigest()[:10]
    # An L2 refinement can keep the starting coordinates. Its IRC must still
    # be calculated at L2, rather than reusing an L1 path with the same hash.
    return f'{reaction.instance_name}_final_{level}_{origin_hash}'


def has_complete_irc_endpoint(qc, job, level):
    """Reuse endpoint properties only when the saved result has the final theory."""
    backend = qc.for_job(job)
    if backend.theory_level != level:
        return False
    rows = list(qc.db.select(name=job))
    if not rows:
        return False
    data = rows[-1].data
    return (data.get('status') == 'normal'
            and (backend.qc != 'fc' or data.get('endpoint_converged') is True)
            and all(data.get(key) is not None
                    for key in ('energy', 'zpe', 'frequencies'))
            and len(data['frequencies']) > 0)


def release_final_ts_clusters(species, par):
    """Cluster selected L1 TSs, then release only representatives to IRC."""
    changes = {}
    try:
        return _release_final_ts_clusters(species, par, changes)
    finally:
        for (path, key), updates in changes.items():
            _update_index_batch(path, key, updates)


def _release_final_ts_clusters(species, par, index_updates):
    from kinbot.bimol_network import geometry_signature, symmetry_rmsd, unique_observation_id

    def record(reaction, **updates):
        update_bimol_status(reaction, index_updates=index_updates, **updates)

    def numerical_duplicate(member, representative):
        return (abs(member.ts.energy + member.ts.zpe - representative.ts.energy - representative.ts.zpe)
                * constants.AUtoKCAL <= par['bimol_energy_tolerance']
                and symmetry_rmsd(member.ts, representative.ts,
                    first_instance=member.instance, second_instance=representative.instance,
                    allow_reflection=(getattr(member, 'allow_reflection', True)
                                      and getattr(representative, 'allow_reflection', True))) <= par['bimol_ts_rmsd'])

    def update_archive(members, representative):
        state = (representative.instance_name, tuple(
            (member.instance_name, int(species.reac_ts_done[index]))
            for index, member in members))
        if getattr(representative, '_final_archive_state', None) == state:
            return
        for _, member in members:
            member.final_conformer_role = ('selected' if member is representative else
                'duplicate' if numerical_duplicate(member, representative)
                else 'archived conformer')
        statuses = {member.instance_name: (index, member) for index, member in members}
        for _, member in members:
            for entry in getattr(member, 'final_conformer_archive', []):
                index, observation = statuses[entry['observation']]
                entry['role'] = observation.final_conformer_role
                entry['final_irc_workflow_status'] = int(species.reac_ts_done[index])
            record(member, final_conformer_archive=getattr(
                member, 'final_conformer_archive', []),
                final_conformer_role=member.final_conformer_role)
        representative._final_archive_state = state

    indices = [index for index, reaction in enumerate(species.reac_obj)
               if getattr(reaction, 'bimol_search', False)
               and not getattr(reaction, 'pose_search_template', False)]
    if not indices:
        return False
    clustered = [index for index in indices if hasattr(
        species.reac_obj[index], 'final_ts_cluster_id')]
    if clustered:
        clusters = {}
        for index in clustered:
            reaction = species.reac_obj[index]
            clusters.setdefault(reaction.final_ts_cluster_id, []).append(
                (index, reaction))
        for cluster_id, members in clusters.items():
            successful = [item for item in members
                          if species.reac_ts_done[item[0]] in {2, 3, 4, 5, -1}]
            waiting = [item for item in members
                       if species.reac_ts_done[item[0]] == 8]
            active = [item for item in members
                      if species.reac_ts_done[item[0]] == 7]
            if successful:
                representative = min(
                    successful, key=lambda item: item[1].final_ts_rank)[1]
                for index, reaction in waiting:
                    # A duplicate can already have validated endpoints from
                    # an earlier run. Read them through the ordinary workflow
                    # so its lower RC/PC remains a candidate for selection.
                    # The geometry-specific names prevent reuse after a TS
                    # change. Do not start work for uncalculated duplicates.
                    job = final_irc_job_name(reaction)
                    jobs = [f'{job}_IRC_{direction}{suffix}'
                            for direction in ('F', 'R')
                            for suffix in ('', '_prod')]
                    qc = getattr(reaction, 'qc', None)
                    completed = (qc is not None and all(
                        qc.check_qc(name) == 'normal' for name in jobs))
                    if completed and par.get('high_level'):
                        # Poll first: FairChem's completed worker files are
                        # copied into the result database by check_qc.
                        completed = all(
                            has_complete_irc_endpoint(qc, f'{job}_IRC_{direction}_prod', 'l2')
                            or qc.check_qc(f'{job}_IRC_{direction}_prod_high') == 'normal'
                            for direction in ('F', 'R'))
                    if completed:
                        reaction.final_irc_name = job
                        reaction.final_irc_attempted = True
                        species.reac_ts_done[index] = 7
                        record(reaction, final_status='restoring_final_irc')
                        logger.info('Reading completed IRC and endpoint calculations for %s.',
                                    reaction.instance_name)
                        continue
                    species.reac_ts_done[index] = -998
                    reaction.final_representative = (
                        representative.instance_name)
                    record(
                        reaction, final_cluster_id=cluster_id,
                        final_representative=representative.instance_name,
                        final_status=('clustered_duplicate' if getattr(
                            reaction, 'final_conformer_role', '') == 'duplicate'
                            else 'archived_conformer'))
                update_archive(members, representative)
                continue
            if active or not waiting:
                continue
            attempted = [reaction for index, reaction in members
                         if species.reac_ts_done[index] == -999
                         and getattr(reaction, 'final_irc_attempted', False)]
            replacement = None
            for item in sorted(waiting,
                               key=lambda value: value[1].final_ts_rank):
                reaction = item[1]
                if all(same_stereopath(reaction, previous)
                       and symmetry_rmsd(
                        reaction.ts, previous.ts,
                        first_instance=reaction.instance,
                        second_instance=previous.instance,
                        allow_reflection=(
                            getattr(reaction, 'allow_reflection', True)
                            and getattr(previous,
                                        'allow_reflection', True))) > 1.e-5
                       for previous in attempted):
                    replacement = item
                    break
            if replacement is None:
                for index, reaction in waiting:
                    species.reac_ts_done[index] = -998
                    record(
                        reaction, final_cluster_id=cluster_id,
                        final_status='untried_exact_duplicate')
                continue
            index, reaction = replacement
            reaction.final_irc_name = final_irc_job_name(reaction)
            reaction.final_irc_attempted = True
            species.reac_ts_done[index] = 7
            for _, member in members:
                member.final_representative = reaction.instance_name
            update_archive(members, reaction)
            record(
                reaction, final_cluster_id=cluster_id,
                final_representative=reaction.instance_name,
                final_status='retrying_final_irc')
        # A failed representative can trigger another pose search. Its new
        # selected TS must still reach IRC when other clusters already exist.
        indices = [index for index in indices if index not in clustered]
        if not indices:
            return True

    if any(species.reac_ts_done[index] not in {6, -999}
           for index in indices):
        return False
    ready = [species.reac_obj[index] for index in indices
             if species.reac_ts_done[index] == 6]
    tolerance = par['bimol_energy_tolerance'] * constants.KCALtoHARTREE
    groups = []
    single_conformer = not par.get('multi_conf_tst', 0)
    for reaction in sorted(ready, key=lambda item: (
            item.ts.energy + item.ts.zpe, item.instance_name)):
        match = None
        for group in groups:
            representative = group[0]
            if reaction.family_name != representative.family_name:
                continue
            if not same_stereopath(reaction, representative):
                continue
            same_site = (tuple(reaction.instance)
                         == tuple(representative.instance))
            identity = _l0_endpoint_identities(reaction)
            if identity is None or identity != _l0_endpoint_identities(representative):
                # Unknown or incompatible configured endpoints need their own
                # final validation. They are not yet independent rate channels.
                continue
            if single_conformer and same_site:
                match = group
                break
            if abs((reaction.ts.energy + reaction.ts.zpe)
                   - (representative.ts.energy + representative.ts.zpe)) \
                    > tolerance:
                continue
            if symmetry_rmsd(
                    reaction.ts, representative.ts,
                    first_instance=reaction.instance,
                    second_instance=representative.instance,
                    allow_reflection=(
                        getattr(reaction, 'allow_reflection', True)
                        and getattr(representative,
                                    'allow_reflection', True))) \
                    <= par['bimol_ts_rmsd']:
                match = group
                break
        if match is None:
            groups.append([reaction])
        else:
            match.append(reaction)

    used_ids = set()
    for group in groups:
        ranked = sorted(group, key=lambda item: (
            item.ts.energy + item.ts.zpe, item.instance_name))
        representative = ranked[0]
        archive = []
        for member in ranked:
            duplicate = member is not representative and numerical_duplicate(member, representative)
            member.final_conformer_role = ('selected' if member is representative else
                                           'duplicate' if duplicate else 'archived conformer')
            archive.append({
                'observation': member.instance_name,
                'role': member.final_conformer_role,
                'geometry_angstrom': np.asarray(member.ts.geom).tolist(),
                'atoms': list(map(str, member.ts.atom)),
                'raw_frequencies_cm-1': list(map(float, getattr(member.ts, 'freq', []))),
                'electronic_energy_hartree': float(member.ts.energy),
                'zpe_hartree': float(member.ts.zpe),
                'source_job': getattr(member.ts, 'source_job', None),
                'conformer_inventory': [dict(record.as_dict(),
                    observation_member_id=f'{member.instance_name}:{record.member_id}')
                    for record in getattr(member.ts, 'conformer_inventory', ())],
                'inventory_atom_mapping': 'same order as this observation atoms',
                'optical_population_scope': getattr(member.ts, 'optical_counting_scope', None),
                'contains_graph_cycle': bool(np.any(getattr(member.ts, 'cycle', []))),
                'equilibration_established': None,
                'independent_rate_contribution': False,
            })
        names = sorted(item.instance_name for item in group)
        identity = '|'.join([
            representative.family_name,
            getattr(representative, 'stereopath_id', ''),
            min(geometry_signature(item.ts) for item in group),
            repr(_l0_endpoint_identities(representative)),
            repr(sorted(tuple(map(int, item.instance)) for item in group)),
            # This is a selection-observation ID, not a chemical channel ID.
            # Coarse geometric hashes cannot distinguish every retained group.
            repr(names),
        ])
        cluster_id = unique_observation_id(f'final_{hashlib.sha1(identity.encode()).hexdigest()[:10]}', used_ids)
        instances = sorted({tuple(instance) for item in group
                            for instance in getattr(
                                item, 'search_instances', [item.instance])})
        for rank, member in enumerate(ranked):
            member.final_ts_cluster_id = cluster_id
            member.final_ts_members = names
            member.final_ts_rank = rank
            member.final_representative = representative.instance_name
            member.final_conformer_archive = archive
            member.final_conformer_policy = ('MC-RRHO members; numerical duplicates removed'
                if not single_conformer else 'lowest selected conformer; other geometries archived')
            member.search_instances = [list(item) for item in instances]
            member_index = species.reac_obj.index(member)
            if member is representative:
                member.final_irc_name = final_irc_job_name(member)
                member.final_irc_attempted = True
                species.reac_ts_done[member_index] = 7
                final_status = 'awaiting_final_irc'
            else:
                species.reac_ts_done[member_index] = 8
                final_status = 'final_irc_fallback'
            record(
                member, final_cluster_id=cluster_id,
                final_representative=representative.instance_name,
                final_ts_rank=rank, final_status=final_status,
                final_conformer_role=member.final_conformer_role,
                final_conformer_archive=archive,
                final_conformer_policy=member.final_conformer_policy)
    return True


def unit_vector(vector):
    vector = np.asarray(vector, dtype=float)
    norm = np.linalg.norm(vector)
    if norm < 1.e-12:
        raise ValueError('Cannot normalize a zero vector.')
    return geometry.unit_vector(vector)


def _perpendiculars(axis):
    axis = unit_vector(axis)
    trial = np.array([1., 0., 0.])
    if abs(np.dot(axis, trial)) > .8:
        trial = np.array([0., 1., 0.])
    first = unit_vector(np.cross(axis, trial))
    return first, unit_vector(np.cross(axis, first))


def rotate_between(source, target):
    source, target = unit_vector(source), unit_vector(target)
    cross = np.cross(source, target)
    dot = float(np.clip(np.dot(source, target), -1., 1.))
    if np.linalg.norm(cross) > 1.e-12:
        return geometry.rotation_matrix(cross, math.acos(dot))
    if dot > 0.:
        return np.eye(3)
    return geometry.rotation_matrix(_perpendiculars(source)[0], math.pi)


def open_directions(fragment, geom, center):
    """Return geometric open-valence directions at a KinBot radical center."""
    neighbors = np.flatnonzero(fragment.bond01[center])
    if not len(neighbors):
        return [np.array([0., 0., 1.])]
    vectors = [unit_vector(geom[index] - geom[center])
               for index in neighbors]
    if len(vectors) == 1:
        angle = math.radians(104.5 if fragment.atom[center] == 'O' else 109.47)
        return [unit_vector(math.cos(angle) * vectors[0]
                            + math.sin(angle)
                            * _perpendiculars(vectors[0])[0])]
    if len(vectors) == 2:
        direction = vectors[0] + vectors[1]
        if np.linalg.norm(direction) < 1.e-8:
            direction = _perpendiculars(vectors[0])[0]
        else:
            direction = -unit_vector(direction)
        return [direction, -direction]
    covariance = sum(np.outer(vector, vector) for vector in vectors)
    direction = np.linalg.eigh(covariance)[1][:, 0]
    return [direction, -direction]


def approach_directions(axis, settings):
    first, second = _perpendiculars(axis)
    for ring, angle in enumerate(settings['pose_angles']):
        count = 1 if abs(angle) < 1.e-8 else settings['pose_azimuths']
        theta = math.radians(angle)
        for azimuth in range(count):
            phi = 2. * math.pi * azimuth / count + ring * GOLDEN_ANGLE
            yield unit_vector(math.cos(theta) * axis + math.sin(theta)
                              * (math.cos(phi) * first
                                 + math.sin(phi) * second))


def _bond_limit(first, second):
    for key in (first + second, second + first):
        if key in constants.st_bond:
            return float(constants.st_bond[key])
    return float(1.2 * (covalent_radii[atomic_numbers[first]]
                        + covalent_radii[atomic_numbers[second]]))


def _bond_graph(atoms):
    return {(i, j) for i in range(len(atoms)) for j in range(i + 1, len(atoms))
            if atoms.get_distance(i, j) < _bond_limit(atoms[i].symbol,
                                                       atoms[j].symbol)}


def _expected_graphs(species, formed_bonds, broken_bonds):
    graph = set()
    shift = 0
    for fragment in (species.fragA, species.fragB):
        graph.update((i + shift, j + shift)
                     for i in range(fragment.natom)
                     for j in range(i + 1, fragment.natom)
                     if fragment.bond01[i][j])
        shift += fragment.natom
    product = graph.copy()
    for bond in broken_bonds:
        product.remove(tuple(sorted(bond)))
    product.update(tuple(sorted(bond)) for bond in formed_bonds)
    return graph, product


def _reaction_state(atoms, formed_bonds, broken_bonds):
    graph = _bond_graph(atoms)
    formed = {tuple(sorted(bond)) for bond in formed_bonds}
    broken = {tuple(sorted(bond)) for bond in broken_bonds}
    if not (graph & formed) and broken <= graph:
        return 'reactant'
    if formed <= graph and not (graph & broken):
        return 'product'
    return 'ambiguous'


def _collision_free(atoms, size_a, formed_bonds, scale):
    intended = {tuple(sorted(bond)) for bond in formed_bonds}
    for i in range(size_a):
        for j in range(size_a, len(atoms)):
            if (i, j) in intended:
                continue
            radii = []
            for number in (atoms.numbers[i], atoms.numbers[j]):
                radius = vdw_radii[number]
                radii.append(radius if np.isfinite(radius)
                             else 1.5 * covalent_radii[number])
            if atoms.get_distance(i, j) < scale * sum(radii):
                return False
    return True


def _bias_constraints(strength, formed_bonds, broken_bonds, settings):
    constraints = [ExternalForce(*bond,
                                 strength * settings['break_scale'])
                   for bond in broken_bonds]
    constraints.extend(ExternalForce(*bond,
                                     -strength * settings['form_scale'])
                       for bond in formed_bonds)
    return constraints


def _distances(atoms, formed_bonds, broken_bonds):
    return {
        'formed': [float(atoms.get_distance(*bond)) for bond in formed_bonds],
        'broken': [float(atoms.get_distance(*bond)) for bond in broken_bonds],
    }


def _trial(start, strength, formed_bonds, broken_bonds, calculator, settings,
           directory):
    class ProductFormed(RuntimeError):
        pass

    atoms = start.copy()
    atoms.calc = calculator
    atoms.set_constraint(_bias_constraints(
        strength, formed_bonds, broken_bonds, settings))
    frames = []
    first_product = None

    def monitor():
        nonlocal first_product
        frame = atoms.copy()
        frame.set_constraint()
        frame.info.update({
            'force': strength,
            'energy': float(atoms.get_potential_energy(
                apply_constraint=False)),
        })
        frames.append(frame)
        if _reaction_state(atoms, formed_bonds, broken_bonds) == 'product':
            first_product = frame.copy()
            raise ProductFormed

    optimizer = BFGS(atoms, logfile=str(directory / 'opt.log'), maxstep=.05)
    optimizer.attach(monitor)
    try:
        converged = bool(optimizer.run(fmax=settings['path_fmax'],
                                       steps=settings['path_steps']))
    except ProductFormed:
        converged = False
    record = {
        'force': strength,
        'energy': float(atoms.get_potential_energy(apply_constraint=False)),
        'distances': _distances(atoms, formed_bonds, broken_bonds),
        'state': _reaction_state(atoms, formed_bonds, broken_bonds),
        'converged': converged,
    }
    atoms.set_constraint()
    write(directory / 'path.xyz', frames)
    _json(directory / 'result.json', record)
    return {'force': strength, 'atoms': atoms.copy(), 'frames': frames,
            'product': record['state'] == 'product',
            'first_product': first_product, 'record': record}


def _tangent(frames, index):
    before = frames[max(0, index - 1)].positions
    after = frames[min(len(frames) - 1, index + 1)].positions
    return unit_vector((after - before).ravel())


def _force_path(encounter, formed_bonds, broken_bonds, calculator, settings,
                directory):
    directory.mkdir(parents=True, exist_ok=True)
    trials = []

    def run(start, force, phase):
        trial_dir = directory / f'{len(trials) + 1:03d}_{phase}'
        trial_dir.mkdir(exist_ok=True)
        result = _trial(start, force, formed_bonds, broken_bonds, calculator,
                        settings, trial_dir)
        trials.append(result)
        return result

    low = {'force': 0., 'atoms': encounter.copy(), 'frames': [],
           'product': False}
    high = None
    for force in np.linspace(0., settings['force_max'],
                             settings['force_levels'])[1:]:
        trial = run(low['atoms'], float(force), 'coarse')
        if trial['product']:
            high = trial
            break
        low = trial
    if high is not None and not low['product']:
        for _ in range(settings['force_bisections']):
            force = .5 * (low['force'] + high['force'])
            trial = run(low['atoms'], force, 'bisect')
            if trial['product']:
                high = trial
            else:
                low = trial
    if high is None or len(high['frames']) < 2:
        return [], {'status': 'no connectivity change',
                    'trials': [trial['record'] for trial in trials]}

    frames = high['frames']
    raw = []

    def add(index, source, source_frames=frames):
        if len(source_frames) > 1:
            raw.append({'atoms': source_frames[index].copy(),
                        'mode': _tangent(source_frames, index),
                        'source': source})

    add(len(frames) - 2, 'pre_connectivity')
    if not low['product']:
        add(len(low['frames']) - 1, 'lower_force', low['frames'])
    energies = [frame.info['energy'] for frame in frames]
    for index in range(1, len(frames) - 1):
        if energies[index] >= max(energies[index - 1], energies[index + 1]):
            add(index, 'energy_maximum')
    coordinate = [
        sum(frame.get_distance(*bond) for bond in broken_bonds)
        - sum(frame.get_distance(*bond) for bond in formed_bonds)
        for frame in frames]
    add(int(np.argmax(np.abs(np.diff(coordinate)))), 'pre_discontinuity')
    if high['first_product'] is not None:
        left = encounter if low['product'] else low['atoms']
        difference = high['first_product'].positions - left.positions
        if np.linalg.norm(difference) > 1.e-10:
            midpoint = left.copy()
            midpoint.positions += .5 * difference
            raw.append({'atoms': midpoint,
                        'mode': unit_vector(difference.ravel()),
                        'source': 'midpoint'})

    candidates = []
    for candidate in raw:
        if all(rmsd.kabsch_rmsd(
                   candidate['atoms'].positions, kept['atoms'].positions,
                   translate=True)
               >= settings['candidate_rmsd'] for kept in candidates):
            candidates.append(candidate)
        if len(candidates) == settings['max_candidates']:
            break
    metadata = {
        'status': 'candidates generated',
        'critical_force': [low['force'], high['force']],
        'trials': [trial['record'] for trial in trials],
        'candidates': [candidate['source'] for candidate in candidates],
    }
    _json(directory / 'path.json', metadata)
    return candidates, metadata


def _endpoint_configuration(atoms):
    """Record the supplied total charge/spin without guessing fragment states."""
    from kinbot.stereo_identity import canonical_identity
    if 'charge' not in atoms.info or 'spin' not in atoms.info:
        return {'status': 'unavailable', 'reason': 'endpoint total charge/spin not supplied'}
    point = StationaryPoint('l0_endpoint', atoms.info['charge'], atoms.info['spin'],
                            atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    try:
        point.characterize()
        return canonical_identity(point)
    except (ValueError, RuntimeError, KeyError) as error:
        return {'status': 'unsupported', 'reason': str(error)}


def _endpoint(atoms, calculator, settings, expected, directory):
    atoms = atoms.copy()
    atoms.calc = calculator
    optimizer = Sella(atoms, order=0, internal=False,
                      trajectory=str(directory / 'opt.traj'),
                      logfile=str(directory / 'opt.log'))
    error = None
    try:
        converged = bool(optimizer.run(
            fmax=settings['endpoint_fmax'],
            steps=settings['endpoint_steps']))
    except (RuntimeError, ValueError) as exc:
        converged = False
        error = str(exc)
    graph = _bond_graph(atoms)
    state = 'reactant' if graph == expected[0] else (
        'product' if graph == expected[1] else 'unexpected')
    write(directory / 'endpoint.xyz', atoms, format='xyz')
    result = {
        'converged': converged,
        'state': state,
        'bonds': [list(bond) for bond in sorted(graph)],
        'distances': _distances(atoms, expected[2], expected[3]),
        'geometry_angstrom': atoms.positions.tolist(),
        'configured_identity': _endpoint_configuration(atoms),
    }
    if error:
        result['error'] = error
    return result, atoms


def _irc_endpoints(ts, calculator, settings, expected, directory):
    """Follow both L0 IRC directions and optimize their endpoint basins."""
    atoms = ts.copy()
    atoms.calc = calculator
    irc = SellaIRC(atoms, dx=.1, eta=1.e-4, gamma=0., keep_going=True,
                   logfile=str(directory / 'irc.log'))
    # Required by Sella 2.5.0 with ASE 3.28; newer Sella overrides this method.
    if irc.gradient_converged.__module__.startswith('ase.'):
        irc.gradient_converged = lambda gradient=None: irc.converged()
    frames = []
    irc.attach(lambda: frames.append(atoms.copy()))
    endpoints = []
    for direction in ('forward', 'reverse'):
        endpoint_dir = directory / f'irc_{direction}'
        endpoint_dir.mkdir(exist_ok=True)
        frames.clear()
        irc.nsteps = 0
        error = None
        try:
            irc_converged = bool(irc.run(
                fmax=settings['irc_fmax'], steps=settings['irc_steps'],
                direction=direction))
        except (RuntimeError, ValueError) as exc:
            irc_converged = False
            error = str(exc)
        if not frames or not np.array_equal(
                frames[-1].positions, atoms.positions):
            frames.append(atoms.copy())
        write(endpoint_dir / 'path.xyz', frames)
        endpoint, _ = _endpoint(
            atoms, calculator, settings, expected, endpoint_dir)
        endpoint.update({'direction': direction,
                         'irc_converged': irc_converged})
        if error:
            endpoint['irc_error'] = error
        endpoints.append(endpoint)
    return endpoints


def _refine_candidate(candidate, charge, mult, calculator, settings,
                      directory):
    directory.mkdir(parents=True, exist_ok=True)
    atoms = candidate['atoms'].copy()
    atoms.info.update({'charge': charge, 'spin': mult})
    atoms.calc = calculator
    optimizer = Sella(atoms, order=1, internal=False,
                      trajectory=str(directory / 'sella.traj'),
                      logfile=str(directory / 'sella.log'))
    # Sella currently expects v0 in its projected free-coordinate basis.
    optimizer.pes.v0 = (unit_vector(candidate['mode'])
                        @ optimizer.pes.get_Ufree())
    error = None
    try:
        converged = bool(optimizer.run(
            fmax=settings.get('sella_fmax', .005),
            steps=settings.get('sella_steps', 300)))
    except (RuntimeError, ValueError) as exc:
        converged = False
        error = str(exc)
    result = {'converged': converged, 'saddle': False}
    if error:
        result['error'] = error
    if converged:
        optimizer.pes.diag()
        free = optimizer.pes.get_Ufree()
        projected = optimizer.pes.get_HL_projected(free)
        result['sella_order'] = int(np.count_nonzero(projected.evals < 0.))
        if result['sella_order'] == 1:
            result['energy'] = float(atoms.get_potential_energy())
            result['saddle'] = True
    write(directory / 'candidate.xyz', atoms, format='xyz')
    _json(directory / 'result.json', result)
    return result, atoms


def matches_expected_product(species, geom, formed_bonds, broken_bonds, alternatives=()):
    """Recognize H transfer at the enumerated site or its equivalent atoms."""
    graph = _bond_graph(Atoms(species.atom, positions=geom))
    choices = [dict(formed_bonds=formed_bonds, broken_bonds=broken_bonds), *alternatives]
    return any(graph == _expected_graphs(species, choice['formed_bonds'], choice['broken_bonds'])[1]
               for choice in choices)


def validate_candidate(species, geom, formed_bonds, broken_bonds, calculator,
                       workdir, search_kwargs=None):
    """Validate one unique L0 saddle by following both IRC directions."""
    settings = _search_settings(search_kwargs)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    atoms = Atoms(species.atom, positions=geom)
    atoms.info.update({'charge': species.charge, 'spin': species.mult})
    expected = _expected_graphs(species, formed_bonds, broken_bonds)
    expected += (formed_bonds, broken_bonds)
    endpoints = _irc_endpoints(
        atoms, calculator, settings, expected, workdir)
    for endpoint in endpoints:
        if endpoint['state'] == 'reactant':
            if _requested_reactant_configuration(species, {'endpoints': [endpoint]},
                    settings.get('optical_population', 'specified')):
                continue
            endpoint['state'] = 'unexpected'
        matched = matches_expected_product(species, endpoint['geometry_angstrom'],
            formed_bonds, broken_bonds, settings.get('equivalent_bond_changes', ()))
        endpoint['state'] = 'product' if matched else 'unexpected'
    states = [endpoint['state'] for endpoint in endpoints]
    observed_product = next((endpoint for endpoint in endpoints
                             if endpoint['state'] != 'reactant'), None)
    result = {
        'endpoints': endpoints,
        'valid': (all(endpoint['converged'] for endpoint in endpoints)
                  and states.count('reactant') == 1),
        'product_matches_expected': (
            observed_product is not None
            and observed_product['state'] == 'product'),
    }
    _json(workdir / 'validation.json', result)
    return result


def _requested_reactant_configuration(species, validation, population):
    """Compare endpoint fragments with the requested configured reactants."""
    from kinbot.stereo_identity import identity_matches, require_supported_identity
    endpoint = next((item for item in validation['endpoints']
                     if item['state'] == 'reactant'), None)
    if endpoint is None or 'geometry_angstrom' not in endpoint:
        return False
    geometry = np.asarray(endpoint['geometry_angstrom'])
    offset = 0
    original_matches, mirror_matches = [], []
    for fragment in (species.fragA, species.fragB):
        observed = StationaryPoint('l0_reactant', fragment.charge, fragment.mult,
            atom=fragment.atom, geom=geometry[offset:offset + fragment.natom])
        observed.characterize()
        reference = require_supported_identity(fragment)
        observation = require_supported_identity(observed)
        original_matches.append(identity_matches(reference, observation, 'specified'))
        mirror_matches.append(identity_matches(dict(reference, id=reference['mirror_id']),
                                               observation, 'specified'))
        offset += fragment.natom
    # Match the same global-mirror population used by final IRC and counting.
    # Reflecting just one of two chiral reactants is a different population.
    return all(original_matches) or (population == 'racemic' and all(mirror_matches))


def find_ts_candidates(species, conformers_a, conformers_b, instance,
                       formed_bonds, broken_bonds, encounter_generator,
                       calculator, workdir, search_kwargs=None, start_after=None):
    """Search one KinBot bimolecular channel for L0 saddle candidates.

    Args:
        species: Combined bimolecular :class:`StationaryPoint`.
        conformers_a: Optimized geometries of ``species.fragA``.
        conformers_b: Optimized geometries of ``species.fragB``.
        instance: Combined indices defining the reaction-family instance.
        formed_bonds: Atom-index pairs attracted by the force bias.
        broken_bonds: Atom-index pairs repelled by the force bias.
        encounter_generator: Reaction-family encounter-pose generator.
        calculator: Unbiased ASE calculator used for paths and refinement.
        workdir: Directory for auditable search output.
        search_kwargs: Overrides for :data:`DEFAULTS`.

    Returns:
        ``(status, geometries, metadata)``. Sequential mode returns the first
        saddle connected to the requested reactant and a cursor for the remaining trials.
        Exhaustive mode returns every locally refined first-order saddle.
    """
    settings = _search_settings(search_kwargs)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    sequential = settings['sequential']
    summary = {'instance': list(instance), 'attempts': [], 'observations': [],
               'sequential': sequential, 'start_after': start_after,
               'conformer_pair_shape': [len(conformers_a), len(conformers_b)],
               'sampling_complete': False,
               'sampling_policy': ('first reactant-connected IRC-validated pose; TS conformers sampled subsequently'
                                   if sequential else 'all encounter poses')}
    if sequential and start_after is not None and (workdir / 'search.json').is_file():
        saved = json.loads((workdir / 'search.json').read_text())
        # Keep completed IRC checks when a failed L1 calculation reopens the
        # remaining poses. Later saved poses do not belong to this continuation.
        summary['observations'] = [item for item in saved.get('observations', [])
            if 'cursor' in item and tuple(item['cursor']) <= tuple(start_after)]
    saddles = []
    for index_a, geom_a in enumerate(conformers_a):
        for index_b, geom_b in enumerate(conformers_b):
            pair_dir = workdir / f'c{index_a:04d}_{index_b:04d}'
            try:
                encounters = encounter_generator(
                    species, geom_a, geom_b, instance, settings)
                for encounter, pose_data in encounters:
                    pose_key = (index_a, index_b, pose_data['pose'])
                    if start_after is not None and pose_key < tuple(start_after[:3]):
                        continue
                    pose = dict(pose_data)
                    pose['conformers'] = [index_a, index_b]
                    pose_dir = pair_dir / f'p{pose["pose"]:03d}'
                    pose_dir.mkdir(parents=True, exist_ok=True)
                    summary['attempts'].append(pose)
                    try:
                        write(pose_dir / 'encounter.xyz', encounter)
                        pose['collision_free'] = _collision_free(
                            encounter, species.fragA.natom, formed_bonds,
                            settings['collision_scale'])
                        if not pose['collision_free']:
                            _json(workdir / 'search.json', summary)
                            continue
                        encounter.info.update({'charge': species.charge,
                                               'spin': species.mult})
                        saved_candidates = pose_dir / 'force_candidates.json'
                        if (start_after is not None and pose_key == tuple(start_after[:3])
                                and saved_candidates.is_file()):
                            saved = json.loads(saved_candidates.read_text())
                            candidates = [dict(atoms=Atoms(species.atom, positions=item['geometry']),
                                               mode=np.asarray(item['mode']), source=item['source'])
                                          for item in saved['candidates']]
                            path = saved['path']
                        else:
                            candidates, path = _force_path(
                                encounter, formed_bonds, broken_bonds, calculator,
                                settings, pose_dir / 'force')
                            if sequential:
                                _json(saved_candidates, {'path': path, 'candidates': [
                                    dict(geometry=item['atoms'].positions.tolist(),
                                         mode=np.asarray(item['mode']).tolist(), source=item['source'])
                                    for item in candidates]})
                        pose['path'] = path['status']
                    except Exception:
                        pose.update({'status': 'error',
                                     'traceback': traceback.format_exc()})
                        _json(workdir / 'search.json', summary)
                        continue
                    for candidate_index, candidate in enumerate(candidates):
                        cursor = [*pose_key, candidate_index]
                        if start_after is not None and tuple(cursor) <= tuple(start_after):
                            continue
                        candidate_dir = pose_dir / (
                            f'candidate_{candidate_index:02d}')
                        try:
                            result, ts = _refine_candidate(
                                candidate, species.charge, species.mult,
                                calculator, settings, candidate_dir)
                        except Exception:
                            candidate_dir.mkdir(parents=True, exist_ok=True)
                            result = {'converged': False, 'saddle': False,
                                      'traceback': traceback.format_exc()}
                            _json(candidate_dir / 'result.json', result)
                            continue
                        if result['saddle']:
                            observation = {
                                'conformers': [index_a, index_b],
                                'pose': pose['pose'],
                                'candidate': candidate_index,
                                'directory': str(candidate_dir),
                                'geometry': ts.positions.tolist(),
                                'result': result,
                            }
                            if sequential:
                                observation['cursor'] = cursor
                                previous = next((item for item in summary['observations']
                                    if _completed_l0_irc(item.get('validation', {}))
                                    and _same_l0_saddle(species, observation, item,
                                        instance, instance,
                                        settings.get('bimol_ts_rmsd', .1),
                                        settings.get('bimol_energy_tolerance', .25),
                                        settings.get('allow_reflection',
                                                     not np.any(species.chiral)))), None)
                                if previous is not None:
                                    # Keep the original endpoint atom order in its
                                    # own files; a symmetry match may permute atoms.
                                    observation['irc_reused_from'] = previous['directory']
                                    summary['observations'].append(observation)
                                    continue
                                validation_dir = candidate_dir / 'validation'
                                try:
                                    validation = validate_candidate(species, ts.positions,
                                        formed_bonds, broken_bonds, calculator,
                                        validation_dir, settings)
                                    validation['reactant_configuration_matches'] = (
                                        _requested_reactant_configuration(species, validation,
                                            settings.get('optical_population', 'specified')))
                                    _json(validation_dir / 'validation.json', validation)
                                except Exception:
                                    validation = {'valid': False, 'traceback': traceback.format_exc()}
                                observation.update(validation=validation,
                                    validation_source=str(validation_dir), cursor=cursor)
                                summary['observations'].append(observation)
                                if not (validation.get('valid')
                                        and validation.get('reactant_configuration_matches')):
                                    continue
                            saddles.append(observation)
                            if sequential:
                                summary.update(status='normal', candidates=saddles,
                                               next_cursor=cursor, exhausted=False)
                                _json(workdir / 'search.json', summary)
                                return 0, [ts.positions.copy()], summary
                    _json(workdir / 'search.json', summary)
            except Exception:
                summary['attempts'].append({
                    'conformers': [index_a, index_b],
                    'status': 'encounter_error',
                    'traceback': traceback.format_exc(),
                })
                _json(workdir / 'search.json', summary)
    summary.update({'status': 'normal' if saddles else 'error',
                    'candidates': saddles, 'exhausted': True,
                    'sampling_complete': True})
    _json(workdir / 'search.json', summary)
    geometries = [np.asarray(candidate['geometry']) for candidate in saddles]
    return (0 if geometries else -1), geometries, summary


def _completed_l0_irc(validation):
    """Reuse completed endpoint checks, including ordinary IRC point limits."""
    endpoints = validation.get('endpoints', [])
    return (not validation.get('traceback') and len(endpoints) == 2
            and 'reactant_configuration_matches' in validation
            and all(endpoint.get('converged') and not endpoint.get('error')
                    and not endpoint.get('irc_error') for endpoint in endpoints))


def _same_l0_saddle(species, first, second, first_instance, second_instance,
                    rmsd_threshold, energy_tolerance, allow_reflection):
    """Compare L0 energy and geometry while preserving reacting atom roles."""
    from kinbot.bimol_network import symmetry_rmsd

    tolerance = energy_tolerance * constants.KCALtoHARTREE / constants.EVtoHARTREE
    energy = first.get('result', {}).get('energy')
    other = second.get('result', {}).get('energy')
    if energy is not None and other is not None and abs(energy - other) > tolerance:
        return False
    return symmetry_rmsd(species, species, first['geometry'], second['geometry'],
        first_instance, second_instance, allow_reflection=allow_reflection) <= rmsd_threshold


def cluster_candidates(species, entries, rmsd_threshold, energy_tolerance):
    """Cluster symmetry-equivalent L0 saddles across reaction channels."""
    ordered = sorted(entries, key=lambda item: (
        item['candidate'].get('result', {}).get('energy', np.inf),
        item['reaction'].instance_name,
        item['candidate'].get('directory', '')))
    groups = []
    for entry in ordered:
        match = None
        for group in groups:
            for member in group:
                if (entry['reaction'].family_name
                        != member['reaction'].family_name):
                    continue
                if not same_stereopath(
                        entry['reaction'], member['reaction']):
                    continue
                if _same_l0_saddle(
                    species, entry['candidate'], member['candidate'],
                    entry['reaction'].instance, member['reaction'].instance,
                    rmsd_threshold, energy_tolerance,
                    allow_reflection=(getattr(
                        entry['reaction'], 'allow_reflection', True)
                        and getattr(member['reaction'],
                                    'allow_reflection', True))):
                    match = group
                    break
            if match is not None:
                break
        if match is None:
            groups.append([entry])
        else:
            match.append(entry)
    return groups


def _group_members(group):
    return [{
        'reaction': member['reaction'].instance_name,
        'instance': list(member['reaction'].instance),
        'source': member['candidate'].get('directory'),
        'energy': member['candidate'].get('result', {}).get('energy'),
    } for member in group]


def _candidate_digest(species, group, representative, append):
    from kinbot.bimol_network import geometry_signature

    candidate = representative['candidate']
    energy = candidate.get('result', {}).get('energy')
    identity = '|'.join(
        [geometry_signature(species, candidate['geometry']), str(energy)]
        + sorted(member['reaction'].instance_name for member in group))
    if append:
        # Distinct retries can have the same rounded geometry signature.
        identity += '|' + candidate['directory'] + '|' + repr(candidate.get('cursor'))
    return hashlib.sha1(identity.encode()).hexdigest()[:10]


def export_unique_candidates(species, groups,
                             directory='bimol/unique_candidates', append=False):
    """Write one unvalidated representative for every L0 saddle cluster."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    index = []
    for group in groups:
        representative = min(group, key=lambda item: item['candidate'].get(
            'result', {}).get('energy', np.inf))
        candidate = representative['candidate']
        energy = candidate.get('result', {}).get('energy')
        digest = _candidate_digest(species, group, representative, append)
        candidate_id = f'candidate_{digest}'
        target = directory / candidate_id
        target.mkdir(exist_ok=True)
        atoms = Atoms(species.atom, positions=candidate['geometry'])
        write(target / 'candidate.xyz', atoms, format='xyz')
        members = _group_members(group)
        metadata = {
            'id': candidate_id,
            'representative': representative['reaction'].instance_name,
            'instance': list(representative['reaction'].instance),
            'energy': energy,
            'members': members,
        }
        _json(target / 'metadata.json', metadata)
        index.append(metadata)
        representative['cluster_id'] = candidate_id.removeprefix('candidate_')
        representative['validation_directory'] = str(target)
        representative['members'] = members
    _write_candidate_index(directory / 'index.json', 'unique_candidates', index, append)
    return groups


def _write_candidate_index(path, key, records, append):
    if append and path.is_file():
        old = json.loads(path.read_text())[key]
        replacements = {item['id']: item for item in records}
        records = [replacements.pop(item['id'], item) for item in old] + list(replacements.values())
    _json(path, {key: records})


def export_validated_l0(species, representatives,
                        directory='bimol/validated_l0', append=False):
    """Export unique L0 saddles with validated L0 IRC connectivity."""
    directory = Path(directory)
    if directory.exists() and not append:
        shutil.rmtree(directory)
    directory.mkdir(parents=True, exist_ok=True)
    index = []
    for representative in representatives:
        ts_id = f'l0_ts_{representative["cluster_id"]}'
        target = directory / ts_id
        target.mkdir(exist_ok=True)
        candidate = representative['candidate']
        write(target / 'ts.xyz', Atoms(
            species.atom, positions=candidate['geometry']), format='xyz')
        source = Path(representative.get('validation_source', representative['validation_directory']))
        validation = representative['validation']
        for endpoint in validation['endpoints']:
            direction = endpoint['direction']
            state = endpoint['state']
            path = source / f'irc_{direction}' / 'path.xyz'
            endpoint_file = source / f'irc_{direction}' / 'endpoint.xyz'
            if path.is_file():
                shutil.copy2(path, target / f'irc_{direction}.xyz')
            if endpoint_file.is_file():
                shutil.copy2(endpoint_file, target / f'{state}.xyz')
        metadata = {
            'id': ts_id,
            'level': 'L0',
            'final_status': 'l0_validated',
            'representative': representative['reaction'].instance_name,
            'instance': list(representative['reaction'].instance),
            'energy': candidate.get('result', {}).get('energy'),
            'members': representative['members'],
            'validation': validation,
        }
        _json(target / 'metadata.json', metadata)
        index.append(metadata)
        representative['unique_id'] = ts_id
    _write_candidate_index(directory / 'index.json', 'validated_l0', index, append)


def _job_data(qc, job):
    rows = list(qc.db.select(name=job))
    return rows[-1].data if rows else {}


_REACTION_FIELDS = ('reac_type', 'reac_inst', 'reac_step', 'reac_scan_energy',
                    'reac_ts_done', 'reac_ts_geom', 'reac_ts_freq', 'reac_name', 'reac_obj')


def _initialize_pose_searches(species, par):
    """Keep one independent continuation for each enumerated reaction site."""
    species.bimol_pose_searches = []
    for index, reaction in enumerate(species.reac_obj):
        if not getattr(reaction, 'bimol_search', False):
            continue
        values = {field: getattr(species, field)[index] for field in _REACTION_FIELDS}
        species.bimol_pose_searches.append(dict(reaction=reaction, values=values,
            attempt=0, cursor=None, active=None, done=False))
        reaction.pose_search_template = True
        species.reac_ts_done[index] = -998
    if species.bimol_pose_searches:
        logger.info('Bimolecular poses will be tried sequentially for each reaction site. '
                    'One reactant-connected L0 IRC-validated TS starts each TS conformer search; '
                    'remaining poses are kept for failed refinement or final IRC.')
        if not par.get('conformer_search'):
            logger.warning('Sequential poses with conformer_search=0 select one TS structure '
                           'per pathway; other conformers and ring basins are not exhaustively sampled.')


def _promote_pose(species, par, qc, state, candidate):
    parent = state['reaction']
    representative = dict(reaction=parent, candidate=candidate,
        validation=candidate['validation'], validation_source=candidate['validation_source'])
    export_unique_candidates(species, [[representative]], append=True)
    export_validated_l0(species, [representative], append=True)
    name = f'{parent.instance_name}_{representative["unique_id"]}'
    variant = type(parent)(species, qc, par, parent.instance, name)
    variant.guess_geom = np.asarray(candidate['geometry'])
    variant.search_candidate = candidate
    variant.l0_validation = candidate['validation']
    variant.search_parent = parent.instance_name
    variant.l0_ts_id = representative['unique_id']
    variant.l0_candidate_id = f'candidate_{representative["cluster_id"]}'
    variant.search_members = representative['members']
    variant.search_instances = [list(parent.instance)]
    variant.sequential_pose_search = True
    variant.pose_search_attempt = state['attempt']
    variant.pose_search_cursor = candidate['cursor']
    variant.pose_search_policy = {
        'mode': 'sequential reactant-connected IRC gate',
        'attempt': state['attempt'], 'accepted_cursor': candidate['cursor'],
        'remaining_poses': 'available for failed refinement or final IRC',
        'ts_conformer_search': bool(par.get('conformer_search')),
        'conformer_and_ring_basin_completeness': 'not established',
    }
    variant.prod_done, variant.valid_prod = 0, []
    for attribute in ('stereopath_id', 'site_degeneracy', 'site_relations',
                      'site_members', 'allow_reflection'):
        if hasattr(parent, attribute):
            setattr(variant, attribute, getattr(parent, attribute))
    values = dict(state['values'], reac_step=0, reac_ts_done=0, reac_scan_energy=[],
                  reac_ts_geom=0, reac_ts_freq=0, reac_name=name, reac_obj=variant)
    index = len(species.reac_obj)
    for field in _REACTION_FIELDS:
        getattr(species, field).append(values[field])
    update_bimol_status(variant, search_policy='sequential validated poses',
                       pose_search_attempt=state['attempt'], pose_search_cursor=candidate['cursor'])
    return index


def advance_pose_searches(species, par, qc):
    """Poll L0 workers; reopen only the pathway whose TS or IRC failed."""
    pending = False
    ready = []
    for state in getattr(species, 'bimol_pose_searches', []):
        if state['done']:
            continue
        if state['active'] is not None:
            index = state['active']
            status = species.reac_ts_done[index]
            if status >= 0:
                continue
            previous = species.reac_obj[index]
            ts_rejected = getattr(previous, 'bimol_final_status', None) in {
                'final_barrier_rejected', 'ts_frequency_failed'}
            unexpected = (getattr(previous, 'final_irc_processed', False)
                          and not getattr(previous, 'product_matches_expected', True))
            if not unexpected and (status != -999 or (
                    getattr(previous, 'final_irc_processed', False) and not ts_rejected)):
                # A later fragment/HIR failure does not justify another TS search.
                state['done'] = True
                continue
            state['attempt'] += 1
            state['active'] = None
            logger.info('Trying remaining poses for %s after %s %s.',
                        state['reaction'].instance_name, previous.instance_name,
                        'found another product' if unexpected else 'failed')
        reaction = state['reaction']
        job = search_job(reaction.instance_name, state['attempt'], sequential=True)
        status = qc.check_qc(job)
        if status == 0:
            qc.qc_bimol_search(reaction, attempt=state['attempt'], start_after=state['cursor'])
            pending = True
            continue
        if status == 'running':
            pending = True
            continue
        search = _job_data(qc, job).get('search', {})
        candidates = search.get('candidates', [])
        if status != 'normal' or not candidates:
            state['done'] = True
            logger.warning('No further reactant-connected TS found for %s (%s).', reaction.instance_name,
                           'poses exhausted' if search.get('exhausted') else 'search worker failed')
            continue
        candidate = candidates[0]
        validation = candidate.get('validation', {})
        if not all(validation.get(key) for key in ('valid',
                                                   'reactant_configuration_matches')):
            logger.warning('Saved pose search %s lacks requested-reactant IRC validation; '
                           'this pathway is not promoted. Use a fresh directory for the new policy.', job)
            state['done'] = True
            continue
        state['cursor'] = candidate['cursor']
        entry = dict(reaction=reaction, candidate=candidate)
        digest = _candidate_digest(species, [entry], entry, append=True)
        job = f'{reaction.instance_name}_l0_ts_{digest}'
        ready.append((state, entry, qc.check_qc(job)))

    # A later-listed site may have finished first in the previous run. Restore
    # its L1 job before an equivalent earlier-listed site can start another one.
    for state, entry, existing in sorted(ready, key=lambda item: item[2] == 0):
        reaction, candidate = entry['reaction'], entry['candidate']
        validation = candidate['validation']
        # Compare with previous promotions, including equivalent reacting sites.
        # Different configured products still need their own L1 calculation.
        identity = _l0_endpoint_identities(reaction, validation)
        tried = [dict(reaction=previous, candidate=previous.search_candidate)
                 for previous in species.reac_obj
                 if hasattr(previous, 'search_candidate')
                 and (getattr(previous, 'search_parent', None) == reaction.instance_name
                      or (identity is not None
                          and identity == _l0_endpoint_identities(previous)))]
        groups = cluster_candidates(species, tried + [entry],
            par.get('bimol_ts_rmsd', .1), par.get('bimol_energy_tolerance', .25))
        group = next(group for group in groups if any(member is entry for member in group))
        previous = next((member['reaction'] for member in group if member is not entry), None)
        if previous is not None and previous.search_parent != reaction.instance_name:
            if existing != 0:
                # Restore existing L1 work on restart. It may have a lower TS
                # or complex; ordinary final selection must still see it.
                previous = None
            else:
                members = previous.search_members
                for member in _group_members([entry]):
                    if member not in members:
                        members.append(member)
                previous.search_instances = [list(instance) for instance in sorted({
                    tuple(member['instance']) for member in members})]
                update_bimol_status(previous, members=members,
                                    search_instances=previous.search_instances)
                state['active'] = species.reac_obj.index(previous)
                pending = True
                logger.info('Using %s for equivalent L0 saddle at %s. Both reaction '
                            'sites retain their remaining poses if refinement or IRC fails.',
                            previous.instance_name, candidate['directory'])
                continue
        if previous is not None:
            logger.info('Skipping L0 saddle at %s: equivalent to %s, already promoted '
                        'to L1. Continuing with the remaining poses.',
                        candidate['directory'], previous.instance_name)
            state['attempt'] += 1
            pending = True
            continue
        state['active'] = _promote_pose(species, par, qc, state, candidate)
    return pending


def prepare_bimolecular_reactions(species, par, qc):
    """Run generic L0 search, clustering, validation, and L1 promotion."""
    if par.get('bimol_search_kwargs', {}).get('sequential', True):
        _initialize_pose_searches(species, par)
        advance_pose_searches(species, par, qc)
        return
    indices = [index for index, reaction in enumerate(species.reac_obj)
               if getattr(reaction, 'bimol_search', False)]
    if not indices:
        return

    while True:
        pending = False
        for index in indices:
            reaction = species.reac_obj[index]
            status = qc.check_qc(search_job(reaction.instance_name))
            if status == 0:
                qc.qc_bimol_search(reaction)
                pending = True
            elif status == 'running':
                pending = True
        if not pending:
            break
        time.sleep(1)

    entries = []
    channels_with_candidates = set()
    for index in indices:
        reaction = species.reac_obj[index]
        job = search_job(reaction.instance_name)
        if qc.check_qc(job) != 'normal':
            continue
        candidates = _job_data(qc, job).get('search', {}).get(
            'candidates', [])
        if candidates:
            channels_with_candidates.add(index)
        entries.extend({'index': index, 'reaction': reaction,
                        'candidate': candidate}
                       for candidate in candidates)

    groups = cluster_candidates(
        species, entries, par.get('bimol_ts_rmsd', .1),
        par.get('bimol_energy_tolerance', .25))
    export_unique_candidates(species, groups)
    representatives = [next(entry for entry in group
                            if 'cluster_id' in entry)
                       for group in groups]

    while True:
        pending = False
        for representative in representatives:
            job = validation_job(representative['validation_directory'])
            status = qc.check_qc(job)
            if status == 0:
                qc.qc_bimol_irc(representative['reaction'], representative)
                pending = True
            elif status == 'running':
                pending = True
        if not pending:
            break
        time.sleep(1)

    validated = []
    candidate_index_file = Path('bimol/unique_candidates/index.json')
    candidate_index = json.loads(candidate_index_file.read_text())
    for representative in representatives:
        job = validation_job(representative['validation_directory'])
        validation = _job_data(qc, job).get(
            'validation', {'valid': False})
        representative['validation'] = validation
        metadata_file = Path(
            representative['validation_directory'], 'metadata.json')
        metadata = json.loads(metadata_file.read_text())
        metadata['validation'] = validation
        _json(metadata_file, metadata)
        candidate_id = f'candidate_{representative["cluster_id"]}'
        for item in candidate_index['unique_candidates']:
            if item['id'] == candidate_id:
                item['validation'] = validation
                break
        if validation.get('valid'):
            validated.append(representative)
    _json(candidate_index_file, candidate_index)
    export_validated_l0(species, validated)

    promoted = {index: [] for index in indices}
    representative_indices = set()
    for representative in representatives:
        representative_indices.add(representative['index'])
    for representative in validated:
        promoted[representative['index']].append(representative)

    fields = ('reac_type', 'reac_inst', 'reac_step', 'reac_scan_energy',
              'reac_ts_done', 'reac_ts_geom', 'reac_ts_freq', 'reac_name',
              'reac_obj')
    expanded = {field: [] for field in fields}
    for index, reaction in enumerate(species.reac_obj):
        values = {field: getattr(species, field)[index] for field in fields}
        if index not in indices:
            for field in fields:
                expanded[field].append(values[field])
            continue
        if promoted[index]:
            for representative in promoted[index]:
                candidate = representative['candidate']
                name = f'{reaction.instance_name}_{representative["unique_id"]}'
                variant = type(reaction)(species, qc, par,
                                         reaction.instance, name)
                variant.guess_geom = np.asarray(candidate['geometry'])
                variant.search_candidate = candidate
                variant.l0_validation = representative['validation']
                variant.search_parent = reaction.instance_name
                variant.l0_ts_id = representative['unique_id']
                variant.l0_candidate_id = (
                    f'candidate_{representative["cluster_id"]}')
                variant.search_members = representative['members']
                variant.search_instances = sorted({
                    tuple(member['instance'])
                    for member in representative['members']})
                for attribute in ('stereopath_id', 'site_degeneracy',
                                  'site_relations', 'site_members',
                                  'allow_reflection'):
                    if hasattr(reaction, attribute):
                        setattr(variant, attribute,
                                getattr(reaction, attribute))
                candidate_values = values.copy()
                candidate_values.update({
                    'reac_step': 0,
                    'reac_ts_done': 0,
                    'reac_scan_energy': [],
                    'reac_ts_geom': 0,
                    'reac_ts_freq': 0,
                    'reac_name': name,
                    'reac_obj': variant,
                })
                for field in fields:
                    expanded[field].append(candidate_values[field])
        elif (index not in channels_with_candidates
              or index in representative_indices):
            values['reac_ts_done'] = -999
            for field in fields:
                expanded[field].append(values[field])

    for field in fields:
        setattr(species, field, expanded[field])
