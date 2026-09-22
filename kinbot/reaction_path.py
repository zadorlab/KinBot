"""Stereochemical pathway classes for locally distinguished reaction motifs.

These are reaction metadata, not species identifiers or rate multipliers.
The two endpoint graphs annotate which bond is broken/formed. Virtual atom
substitution on those graphs at the TS geometry retains fleeting stereo.
It also makes the class independent of atom numbering and search direction.
"""
import copy
import hashlib
import json
import logging

import numpy as np

from kinbot.stereo_identity import canonical_identity, optical_scope
from kinbot.stereochemistry import (refine_equivalence_group, virtually_labelled,
                                    configuration_erased_graph)


class StereoAssignmentUnavailable(ValueError):
    """A demonstrated pathway split lacks a supported stereo assignment."""


def endpoint_snapshot(species):
    """Retain the IRC graph/order before product reuse or optimization changes it."""
    from kinbot.stationary_pt import StationaryPoint
    snapshot = StationaryPoint('irc_endpoint', species.charge, species.mult,
        atom=copy.deepcopy(species.atom), geom=copy.deepcopy(species.geom))
    # Keep the legacy chirality methods available when optional RDKit is
    # absent, without copying QC handles or the live reaction-object graph.
    names = ('bond', 'bond01', 'bonds', 'rads', 'atom_eqv', 'atomid',
             'isotopes', 'formal_charges')
    for name in names:
        if hasattr(species, name):
            setattr(snapshot, name, copy.deepcopy(getattr(species, name)))
    return snapshot


def set_endpoint_populations(ts, reactants, products):
    """A global TS mirror must belong to both connected populations."""
    ts.ts_endpoint_identities = tuple(tuple(canonical_identity(p) for p in side)
                                      for side in (reactants, products))
    if not hasattr(ts, 'optical_reference') and len(reactants) == 1:
        ts.optical_reference = canonical_identity(reactants[0])


def _configured_pair(endpoints, geom):
    identities = [canonical_identity(endpoint, geom) for endpoint in endpoints]
    if any(item['status'] != 'assigned' for item in identities):
        return None, None
    return tuple(item['id'] for item in identities), tuple(item['mirror_id'] for item in identities)


def prepare_ts_context(ts, reactant, product):
    """Save endpoint population and spectator configuration before TS searches."""
    set_endpoint_populations(ts, [reactant], [product])
    changed = np.flatnonzero(np.any(np.asarray(reactant.bond) != product.bond, axis=1))
    ts.ts_endpoint_graphs = tuple(endpoint_snapshot(p) for p in (reactant, product))
    ts.ts_reference_geometry = np.array(ts.geom, copy=True)
    endpoints = tuple(copy.copy(p) for p in ts.ts_endpoint_graphs)
    for endpoint in endpoints:
        endpoint.stereo_ignored_atoms = tuple(map(int, changed))
    ts.configuration_endpoints = endpoints
    ts.configuration_reference, ts.configuration_mirror = _configured_pair(endpoints, ts.geom)
    ts.configuration_ignored_atoms = tuple(map(int, changed))


def _tagged_pair(species, geom):
    endpoints = [copy.copy(endpoint) for endpoint in species.stereopath_endpoints]
    for endpoint in endpoints:
        # A forming/breaking double bond can be nonplanar in a TS. Its
        # fictitious E/Z tag must not depend on which substituent RDKit visits.
        endpoint.stereo_ignored_bonds = getattr(species, 'stereopath_ignored_bonds', ())
    roles = getattr(species, 'stereopath_joint_roles', None)
    if roles is not None:
        identities = [canonical_identity(virtually_labelled(endpoint, roles), geom)
                      for endpoint in endpoints]
        if any(item['status'] != 'assigned' for item in identities):
            raise StereoAssignmentUnavailable('Cannot assign a demonstrated joint stereochemical pathway.')
        return (tuple(sorted(item['id'] for item in identities)),
                tuple(sorted(item['mirror_id'] for item in identities)))
    atoms = getattr(species, 'stereopath_atoms', None)
    if atoms is None:
        atoms = [species.stereopath_hydrogen]  # existing saved objects
    pairs = []
    mirrors = []
    for atom in atoms:
        identities = [canonical_identity(endpoint, geom, tagged_atom=atom)
                      for endpoint in endpoints]
        if any(identity['status'] != 'assigned' for identity in identities):
            raise StereoAssignmentUnavailable('Cannot assign a demonstrated stereochemical pathway. '
                                              'A supported canonical stereo assignment is required.')
        pairs.append(tuple(sorted(identity['id'] for identity in identities)))
        mirrors.append(tuple(sorted(identity['mirror_id'] for identity in identities)))
    # Preserve the original single-transfer class convention.
    return ((pairs[0], mirrors[0]) if len(pairs) == 1
            else (tuple(sorted(pairs)), tuple(sorted(mirrors))))



def _reaction_roles(reactant, product, changed):
    """Number chemical roles without relying on atom order or search direction."""
    descriptions = {}
    for atom in changed:
        edits = []
        for neighbor in changed:
            before, after = int(reactant.bond[atom, neighbor]), int(product.bond[atom, neighbor])
            if before != after:
                edits.append((str(reactant.atom[neighbor]), min(before, after), max(before, after)))
        descriptions[atom] = (str(reactant.atom[atom]), tuple(sorted(edits)))
    numbers = {value: i for i, value in enumerate(sorted(set(descriptions.values())))}
    return {atom: numbers[value] for atom, value in descriptions.items()}


def _has_joint_stereo_split(endpoints, roles):
    """Detect diastereotopy created by the other reacting sites' fixed roles.

    Conditioning also separates ordinary graph positions. Only compare
    alternatives with the same isotope-labelled graph after removing stereo;
    a difference between those mirror families demonstrates a stereo split.
    """
    for endpoint in endpoints:
        for atom in roles:
            group = next((g for g in endpoint.atom_eqv if atom in g), [atom])
            if len(group) < 2 or len(refine_equivalence_group(endpoint, group)) > 1:
                continue  # already covered by the original single-site classes
            previous = {other: role for other, role in roles.items() if other != atom}
            available = [candidate for candidate in group if candidate not in previous]
            if len(available) < 2:
                continue
            view = virtually_labelled(endpoint, previous)
            identities = [canonical_identity(view, tagged_atom=candidate) for candidate in available]
            if any(item['status'] != 'assigned' for item in identities):
                continue
            families = {}
            for identity in identities:
                graph = configuration_erased_graph(identity)
                families.setdefault(graph, set()).add(identity['mirror_family_id'])
            if any(len(values) > 1 for values in families.values()):
                return True
    return False


def prepare_stereopath(ts, reactant, product):
    """Protect all TS spectators; classify demonstrated reacting-site splits.

    Site classification is element-neutral and allows additional bond edits,
    as in HO2 elimination. Ordinary unsplit paths keep lowest-barrier policy.
    """
    prepare_ts_context(ts, reactant, product)
    for key in ('stereopath_atoms', 'stereopath_hydrogen', 'stereopath_endpoints',
                'stereopath_product_identity', 'stereopath_reference',
                'stereopath_mirror', 'stereopath_id', 'stereopath_metadata',
                'stereopath_joint_roles', 'stereopath_ignored_bonds'):
        ts.__dict__.pop(key, None)
    delta = np.asarray(product.bond) - np.asarray(reactant.bond)
    broken = np.argwhere(np.triu(delta, 1) < 0)
    formed = np.argwhere(np.triu(delta, 1) > 0)
    changed = sorted(map(int, set(broken.ravel()) | set(formed.ravel())))
    split_atoms = [int(atom) for atom in changed if any(
        len(refine_equivalence_group(endpoint, group)) > 1
        for endpoint in (reactant, product) for group in endpoint.atom_eqv
        if atom in group and len(group) > 1)]
    roles = _reaction_roles(reactant, product, changed)
    joint = len(split_atoms) > 1 or _has_joint_stereo_split((reactant, product), roles)
    if not split_atoms and not joint:
        from kinbot.stereo_identity import legacy_stereo_warning
        for side in ts.ts_endpoint_identities:
            for identity in side:
                if identity['status'] != 'assigned':
                    legacy_stereo_warning(ts, identity.get('reason'))
        # An examined ordinary path is not missing legacy classification.
        # All ordinary routes retain one lowest-barrier class per endpoint pair.
        ts.stereopath_id = 'ordinary'
        ts.stereopath_metadata = {
            'schema': 'kinbot.stereopath.v1', 'id': ts.stereopath_id,
            'site_relation': 'no demonstrated stereochemical split',
            'atom_index_base': 0, 'broken_bonds': broken.tolist(),
            'formed_bonds': formed.tolist(),
            'counting': 'lowest barrier among ordinary routes to these endpoints',
        }
        return
    if joint:
        ts.stereopath_joint_roles = roles
        split_atoms = changed
    ts.stereopath_atoms = split_atoms
    if not joint and len(split_atoms) == 1 and reactant.atom[split_atoms[0]] == 'H':
        ts.stereopath_hydrogen = split_atoms[0]
    ts.stereopath_endpoints = tuple(endpoint_snapshot(p) for p in (reactant, product))
    ts.stereopath_ignored_bonds = [list(map(int, pair)) for pair in np.argwhere(np.triu(delta != 0, 1))
                                   if max(reactant.bond[tuple(pair)], product.bond[tuple(pair)]) > 1]
    ts.stereopath_product_identity = canonical_identity(product)
    own, mirror = _tagged_pair(ts, ts.geom)
    ts.stereopath_reference, ts.stereopath_mirror = own, mirror
    prefix = 'htransfer:' if hasattr(ts, 'stereopath_hydrogen') else 'stereopath:'
    ts.stereopath_id = prefix + hashlib.sha256(repr(min(own, mirror)).encode()).hexdigest()
    ts.stereopath_metadata = {
        'schema': 'kinbot.stereopath.v1', 'id': ts.stereopath_id,
        'reacting_atoms': split_atoms, 'atom_index_base': 0,
        'broken_bonds': broken.tolist(), 'formed_bonds': formed.tolist(),
        'site_relation': 'joint diastereotopic' if joint else 'diastereotopic',
        'counting': 'one explicit route; no additional route-number multiplier',
    }
    if hasattr(ts, 'stereopath_hydrogen'):
        ts.stereopath_metadata['transferred_hydrogen'] = ts.stereopath_hydrogen


def path_geometry_allowed(species, geom, population=None):
    """Protect spectator configuration and any separately classified pathway."""
    if not getattr(species, 'wellorts', 0) and not hasattr(species, 'stereopath_reference'):
        from kinbot.stereo_identity import configured_geometry_allowed
        return configured_geometry_allowed(species, geom, population)
    population = population or getattr(species, 'optical_population', 'specified')
    mirror_allowed = optical_scope(species, population)['mirror_allowed']
    reference = getattr(species, 'configuration_reference', None)
    if reference is not None:
        own, _ = _configured_pair(species.configuration_endpoints, geom)
        allowed = [reference]
        if mirror_allowed:
            allowed.append(species.configuration_mirror)
        if own not in allowed:
            return False
    if getattr(species, 'ts_endpoint_graphs', ()):
        from kinbot.stereo_identity import endpoint_configuration_allowed
        references = [species.ts_reference_geometry]
        if mirror_allowed:
            references.append(species.ts_reference_geometry * [-1., 1., 1.])
        try:
            if not any(all(endpoint_configuration_allowed(endpoint, reference, geom)
                           for endpoint in species.ts_endpoint_graphs) for reference in references):
                return False
        except (ImportError, ValueError, RuntimeError):
            # Ordinary unsupported configurations retain their existing guard.
            if getattr(species, 'configuration_reference', None) is not None:
                raise
    if hasattr(species, 'stereopath_reference'):
        own, _ = _tagged_pair(species, geom)
        allowed = [species.stereopath_reference]
        if mirror_allowed:
            allowed.append(species.stereopath_mirror)
        return own in allowed
    return True


def reaction_path_id(reaction):
    ts = reaction.ts
    if not hasattr(ts, 'stereopath_metadata'):
        return None
    if not path_geometry_allowed(ts, ts.geom):
        raise ValueError(f'Selected TS {reaction.instance_name} changed stereochemical pathway.')
    return getattr(ts, 'stereopath_id', None)


def same_path_class(left, right):
    """Lowest-barrier filtering stays within a demonstrated pathway class."""
    if (left is None) != (right is None):
        return False
    return left == right


def reject_invalid_pathway(species, index, par=None):
    """Exclude an unusable channel, preserving unrelated successful chemistry."""
    if species.reac_ts_done[index] != -1:
        return False
    reaction = species.reac_obj[index]
    try:
        population = (par or {}).get('optical_population', getattr(species, 'optical_population', 'specified'))
        if population == 'racemic':
            products = ([opt.species for opt in getattr(reaction, 'prod_opt', ())]
                        or reaction.products)
            identities = [canonical_identity(point) for point in products]
            if sum(item.get('is_chiral_configuration', False) for item in identities) > 1:
                raise ValueError('Independent racemic fragments are not the two correlated '
                                 'product configurations of a global mirror pair; this product '
                                 'population is not represented by the current writer')
        if species.reac_type[index] != 'hom_sci':
            reaction_path_id(reaction)
    except ValueError as error:
        logging.getLogger('KinBot').warning(
            '%s: omitted reaction from rates: %s. Calculation files retained; '
            'the exported network is incomplete.', reaction.instance_name, error)
        reaction.stereochemical_rejection = str(error)
        species.reac_ts_done[index] = -999
        return True
    return False


def compare_pathways(existing_name, existing_path, existing_energy,
                     candidate_name, candidate_path, candidate_energy, *,
                     candidate_complex=False):
    """Select barriers with matching endpoints: replace, keep, or distinct.

    A homolytic-scission placeholder is superseded by an optimized saddle,
    before comparing stereochemical saddle classes. Within an actual class,
    retain the lowest barrier; PES can retain its equal-energy complex tie rule.
    """
    existing_hom = 'hom_sci' in existing_name
    candidate_hom = 'hom_sci' in candidate_name
    if existing_hom or candidate_hom:
        return 'replace' if existing_hom and not candidate_hom else 'keep'
    if (existing_path is None) != (candidate_path is None):
        logging.getLogger('KinBot').warning(
            'Pathway classification missing for %s; retaining the classified route %s '
            'for these endpoints. The unclassified observation is retained in the '
            'calculation files, not counted as an additional pathway.',
            existing_name if existing_path is None else candidate_name,
            candidate_name if existing_path is None else existing_name)
        # Never let an unclassified route merge two known distinct classes.
        return 'replace' if existing_path is None else 'keep'
    if not same_path_class(existing_path, candidate_path):
        return 'distinct'
    if candidate_energy < existing_energy:
        return 'replace'
    if candidate_energy == existing_energy and candidate_complex:
        return 'replace'
    return 'keep'


def summary_path_line(reaction):
    if reaction_path_id(reaction) is None:
        return None
    return '# kinbot_stereopath ' + json.dumps(
        dict(reaction.ts.stereopath_metadata, reaction=reaction.instance_name), sort_keys=True)


def read_summary_paths(lines):
    paths = {}
    for line in lines:
        if line.startswith('# kinbot_stereopath '):
            record = json.loads(line[len('# kinbot_stereopath '):])
            if (record.get('schema') != 'kinbot.stereopath.v1'
                    or not isinstance(record.get('id'), str) or not record['id']
                    or not isinstance(record.get('reaction'), str) or not record['reaction']):
                raise ValueError('Unsupported or incomplete reaction-path metadata in summary.')
            name, key = record['reaction'], record['id']
            if name in paths and paths[name] != key:
                raise ValueError(f'Contradictory stereochemical pathway records for {name}.')
            paths[name] = key
    return paths
