"""Select bimolecular rate models without discarding calculation observations.

Conformer counting is completed here, while the original graphs and Hessians
are still available. The MESS reader needs only the resulting saved properties.
"""
import copy
from dataclasses import replace
import logging

import networkx as nx
import numpy as np

from kinbot.calculation import array_fingerprint
from kinbot.conformer_counting import evaluate_members
from kinbot.conformer_records import ConformerRecord, retain
from kinbot.molecular_symmetry import chemical_graph
from kinbot.reaction_path import compare_pathways, reaction_path_id
from kinbot.stereo_identity import canonical_identity
from kinbot import frequencies

logger = logging.getLogger('KinBot')


def population_route(left, right):
    """Identify a route and its global mirror, not independent endpoint flips."""
    return min((left[0], right[0]), (left[1], right[1]))


def _resonance_states(point, order):
    matrices = list(getattr(point, 'bonds', ())) or [point.bond]
    radicals = list(getattr(point, 'rads', ()))
    return {tuple(np.asarray(matrix)[np.ix_(order, order)].ravel()) +
            tuple(np.asarray(radicals[i])[order] if i < len(radicals) else ())
            for i, matrix in enumerate(matrices)}


def atom_order(target, source):
    """Map an observation into the reference order, preserving reacting roles."""
    if (target.natom, target.charge, target.mult) != (source.natom, source.charge, source.mult):
        raise ValueError('Conformer atoms, charge or multiplicity differ.')
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        chemical_graph(target), chemical_graph(source),
        node_match=lambda a, b: a['label'] == b['label'],
        edge_match=lambda a, b: a['label'] == b['label'])
    reference = _resonance_states(target, list(range(target.natom)))
    for mapping in matcher.isomorphisms_iter():
        order = [mapping[i] for i in range(target.natom)]
        if _resonance_states(source, order) == reference:
            return order
    raise ValueError('Conformer chemical graphs or reacting-bond roles differ.')


def ensemble(points, population):
    """Count all observed conformers together; retain each proper duplicate once."""
    reference = copy.copy(min(points, key=lambda p: p.energy + p.zpe))
    reference.bimol_reactant_only = False
    records = []
    for point in points:
        order = atom_order(reference, point)
        existing = [r for r in getattr(point, 'conformer_inventory', ())
                    if r.retained and r.status == 'valid']
        if not existing:
            existing = [ConformerRecord(str(point.name), 0,
                getattr(point, 'source_job', None), 'valid',
                geometry=tuple(map(tuple, point.geom)),
                electronic_energy_hartree=point.energy, zpe_hartree=point.zpe,
                zero_energy_hartree=point.energy + point.zpe,
                frequencies_cm1=tuple(point.freq),
                hessian=tuple(map(tuple, getattr(point, 'hess', ()))) or None,
                hessian_reference=getattr(point, 'optical_hessian_reference', None))]
        for old in existing:
            geom = np.asarray(old.geometry)[order]
            hessian, hessref = old.hessian, copy.deepcopy(old.hessian_reference)
            if hessian is not None and (not hessref or
                    hessref.get('geometry_sha256') != array_fingerprint(old.geometry) or
                    hessref.get('hessian_sha256') != array_fingerprint(hessian) or
                    hessref.get('atoms') != list(map(str, point.atom)) or
                    hessref.get('source_job') != old.source_job):
                hessian, hessref = None, None
            if hessian is not None:
                indices = [3 * i + j for i in order for j in range(3)]
                hessian = np.asarray(hessian)[np.ix_(indices, indices)]
                if hessref is not None:
                    hessref.update(geometry_sha256=array_fingerprint(geom),
                        hessian_sha256=array_fingerprint(hessian),
                        atoms=list(map(str, reference.atom)), original_atom_indices=order)
                hessian = tuple(map(tuple, hessian))
            index = len(records)
            records.append(replace(old, index=index,
                member_id=f'{reference.name}:observed:{index}', geometry=tuple(map(tuple, geom)),
                hessian=hessian, hessian_reference=hessref, retained=False,
                duplicate_of=None, exclusion_reason=None, mirror_partner_ids=(),
                remaining_optical_weight=None, mirror_coverage=None,
                optical_evidence={'warnings': (old.optical_evidence or {}).get('warnings', [])}))
    records, groups = evaluate_members(reference, records, population)
    accepted = [records[i] for group in groups for i in group]
    if not accepted:
        raise ValueError('No conformer remains in the requested configured population.')
    retain(reference, records, [r.index for r in accepted])
    reference.conformer_index = [r.index for r in accepted]
    reference.conformer_geom = [np.asarray(r.geometry) for r in accepted]
    reference.conformer_freq = [list(r.frequencies_cm1) for r in accepted]
    reference.conformer_zeroenergy = [r.zero_energy_hartree for r in accepted]
    reference.conformer_energy = [r.electronic_energy_hartree for r in accepted]
    reference.conformer_representation = 'MC-RRHO'
    reference.hir = None
    reference.rotor_projection = None
    reference.reduced_freqs = list(reference.freq)
    return reference


def prepare_models(network):
    """One complex population per fragment pair; one model per chemical path."""
    from kinbot.bimol_network import _species_data
    mc = bool(network.par.get('multi_conf_tst'))
    population = network.par.get('optical_population', 'specified')
    warnings, complexes, pathways = [], [], []
    if not hasattr(network, 'species'):
        return dict(complexes=[], pathways=[], warnings=['Only raw calculation evidence is available.'])
    aliases = {key: key for key in network.nodes}
    unsupported = set()
    configurations = {}
    if population == 'racemic':
        families = {}
        for key, node in network.nodes.items():
            if node['kind'] != 'bimolecular':
                continue
            identities = [canonical_identity(p) for p in node['fragment_species']]
            configurations[key] = tuple(tuple(sorted((p.chemid, p.charge, p.mult, i[field])
                for p, i in zip(node['fragment_species'], identities))) for field in ('id', 'mirror_id'))
            if sum(i['is_chiral_configuration'] for i in identities) > 1:
                unsupported.add(key)
                warnings.append(f'{key}: two independent chiral reactant populations require separate configured-pair calculations.')
                continue
            signature = tuple(sorted((p.chemid, p.charge, p.mult, i['mirror_family_id'])
                                     for p, i in zip(node['fragment_species'], identities)))
            families.setdefault(signature, []).append(key)
        for keys in families.values():
            chosen = (network.reference_node if network.reference_node in keys else
                      min(keys, key=lambda k: (network.nodes[k]['energy'], k)))
            for key in keys:
                aliases[key] = chosen
    groups = {}
    for original_transition in network.transitions:
        transition = dict(original_transition)
        reaction = transition['reaction']
        try:
            path = reaction_path_id(reaction)
            original = (reaction.reactant_asymptote_id, reaction.product_asymptote_id)
            if any(network.nodes[key]['kind'] != 'bimolecular' for key in original):
                raise ValueError('This H-abstraction rate model requires two separated fragments on each side.')
            if any(key in unsupported for key in original):
                continue
            endpoints = tuple(aliases[key] for key in original)
            transition['_population_route'] = (population_route(*(configurations[key] for key in original))
                                                if population == 'racemic' else None)
            if path is None:
                raise ValueError('Missing stereochemical pathway information.')
            for group in groups.setdefault(endpoints, []):
                first = group[0]
                if first['_population_route'] != transition['_population_route']:
                    continue
                decision = compare_pathways(first['reaction'].instance_name,
                    reaction_path_id(first['reaction']), first['energy'],
                    reaction.instance_name, path, transition['energy'])
                if decision == 'distinct':
                    continue
                group.append(transition)
                if decision == 'replace':
                    group.insert(0, group.pop())
                break
            else:
                groups[endpoints].append([transition])
        except ValueError as error:
            warnings.append(f"{transition['id']}: {error}")
    complex_groups = {}
    # Compare every optimized complex, including observations merged earlier
    # as geometrical duplicates, but use only the lowest E + ZPE for rates.
    for reaction in [r for i, r in enumerate(network.species.reac_obj)
                     if network.species.reac_ts_done[i] == -1]:
        for role, flag, attr, endpoint in (
                ('rc', 'do_reac_vdW', 'irc_reactant_opt', reaction.reactant_asymptote_id),
                ('pc', 'do_vdW', 'irc_prod_opt', reaction.product_asymptote_id)):
            if getattr(reaction, flag, False):
                if endpoint not in unsupported:
                    complex_groups.setdefault((role, aliases[endpoint]), []).append(getattr(reaction, attr).species)
    complex_ids = {}
    for (role, endpoint), points in sorted(complex_groups.items()):
        try:
            point = min(points, key=lambda p: p.energy + p.zpe)
            if any(f <= 0 or not np.isfinite(f) for f in
                   frequencies.thermochemical_frequencies(point.freq, False,
                       network.par.get('imagfreq_threshold', 50.))):
                raise ValueError('Complex has an invalid minimum frequency; use separated fragments instead.')
            identifier = f'{role}_{endpoint}'
            complexes.append(dict(id=identifier, role=role, asymptote=endpoint,
                stationary_point=_species_data(point),
                sampling='lowest E + ZPE optimized IRC complex; single harmonic conformer'))
            complex_ids[role, endpoint] = identifier
        except ValueError as error:
            warnings.append(f'{role} {endpoint}: {error}')
    for endpoints, path_groups in sorted(groups.items()):
        for group in path_groups:
            best = group[0]
            try:
                observations = [reaction for entry in group
                                for reaction in entry.get('_observations', [entry['reaction']])]
                point = ensemble([reaction.ts for reaction in observations], population) if mc else best['species']
                # Different complex records can keep equivalent TS observations
                # apart. Check the reacting sites of every observation used here.
                counting = network._site_counting(dict(best, _observations=observations))
                if counting['unlabelled_model_conversion']['additional_site_multiplier'] != 1.:
                    raise ValueError('Reaction-site counting has conflicting or missing observations.')
                pathways.append(dict(id=best['id'], reactant=endpoints[0], product=endpoints[1],
                    reactant_complex=complex_ids.get(('rc', endpoints[0])),
                    product_complex=complex_ids.get(('pc', endpoints[1])),
                    path_id=reaction_path_id(best['reaction']), site_counting=counting,
                    population_route=best['_population_route'],
                    observations=[entry['id'] for entry in group],
                    stationary_point=_species_data(point)))
            except ValueError as error:
                warnings.append(f"{best['id']}: {error}")
    for message in warnings:
        logger.warning('Bimolecular rate model omitted: %s Calculation observations retained.', message)
    return dict(complexes=complexes, pathways=pathways, warnings=warnings,
        asymptote_aliases=aliases,
        complex_population_assumption='one lowest E + ZPE complex per configured fragment pair and role',
        conformer_treatment='MC-RRHO' if mc else 'lowest E + ZPE within each chemical pathway')
