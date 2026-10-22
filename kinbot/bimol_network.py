"""Canonical reaction network for bimolecular KinBot searches."""

import hashlib
import copy
import json
import logging
import shutil
import tempfile
from pathlib import Path

import networkx as nx
import numpy as np
from kinbot.stereo_identity import canonical_identity
import rmsd
from ase import Atoms
from ase.io import write

from kinbot import constants, kb_path
from kinbot.bimol import same_stereopath
from kinbot.thermochemistry import thermochemistry_evidence, json_record
from kinbot.counting_contract import site_counting
from kinbot.reaction_path import compare_pathways, reaction_path_id


def _digest(value):
    return hashlib.sha1(value.encode()).hexdigest()[:10]


def unique_observation_id(base, used):
    """Allocate a technical ID; coarse/hash identity is never proof of a basin."""
    value, suffix = base, 2
    while value in used:
        value = f'{base}_{suffix}'
        suffix += 1
    used.add(value)
    return value


def _atomic_text(path, value):
    """Replace one text artifact only after its complete contents are ready."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
                mode='w', dir=path.parent,
                prefix=f'.{path.name}.', delete=False) as handle:
            handle.write(value)
            temporary = Path(handle.name)
        temporary.replace(path)
    except Exception:
        if temporary is not None and temporary.exists():
            temporary.unlink()
        raise


def _graph(species):
    """Build a graph using KinBot's characterized chemical information."""
    graph = nx.Graph()
    for index, element in enumerate(species.atom):
        radicals = tuple(sorted(
            int(radical_pattern[index])
            for radical_pattern in getattr(species, 'rads', [])
        ))
        graph.add_node(
            index, element=str(element),
            atomid=int(species.atomid[index]), radicals=radicals)
    for first in range(species.natom):
        for second in range(first):
            if not species.bond01[first][second]:
                continue
            orders = {int(matrix[first][second])
                      for matrix in getattr(species, 'bonds', [])
                      if matrix[first][second]}
            if species.bond[first][second]:
                orders.add(int(species.bond[first][second]))
            graph.add_edge(first, second, orders=tuple(sorted(orders)))
    return graph


def symmetry_rmsd(first, second, first_geom=None, second_geom=None,
                  first_instance=None, second_instance=None,
                  allow_reflection=True):
    """Return the lowest RMSD over graph automorphisms and reflections."""
    if (first.natom != second.natom or first.charge != second.charge
            or first.mult != second.mult):
        return np.inf
    if first_instance is not None and (second_instance is None
            or len(first_instance) != len(second_instance)):
        return np.inf
    left = first.geom if first_geom is None else np.asarray(first_geom)
    right = second.geom if second_geom is None else np.asarray(second_geom)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _graph(first), _graph(second),
        node_match=lambda a, b: (
            a['element'] == b['element']
            and a['atomid'] == b['atomid']
            and a['radicals'] == b['radicals']),
        edge_match=lambda a, b: a['orders'] == b['orders'])
    best = np.inf
    for mapping in matcher.isomorphisms_iter():
        if first_instance is not None and any(
                mapping[a] != b
                for a, b in zip(first_instance, second_instance)):
            continue
        reordered = np.asarray([right[mapping[i]]
                                for i in range(first.natom)])
        value = rmsd.kabsch_rmsd(left, reordered, translate=True)
        if allow_reflection:
            value = min(value, rmsd.kabsch_rmsd(
                left, reordered * [-1., 1., 1.], translate=True))
        best = min(best, value)
        if best < 1.e-5:
            break
    return best


def geometry_signature(species, geom=None):
    geom = species.geom if geom is None else np.asarray(geom)
    distances = []
    for first in range(species.natom):
        for second in range(first):
            pair = ''.join(sorted((str(species.atom[first]),
                                   str(species.atom[second]))))
            distances.append((pair, round(float(np.linalg.norm(
                geom[first] - geom[second])), 2)))
    return _digest(json.dumps(sorted(distances)))


def _species_data(species):
    if getattr(species, 'bimol_reactant_only', False):
        species = copy.copy(species)
        species.bimol_reactant_only = False
    data = {
        'atom': [str(atom) for atom in species.atom],
        'geometry': np.asarray(species.geom).tolist(),
        'energy': float(species.energy),
        'zpe': float(species.zpe),
        'frequencies': [float(freq) for freq in
                        getattr(species, 'reduced_freqs', [])],
        'charge': int(species.charge),
        'mult': int(species.mult),
        'name': str(species.name),
        'chemid': str(species.chemid),
        'bond': np.asarray(species.bond).tolist(),
        # Legacy fields above retain their values; units are local to this
        # stationary-point record, not the network's relative kcal/mol scale.
        'energy_unit': 'hartree', 'geometry_unit': 'angstrom',
        'frequency_unit': 'cm^-1',
        'frequencies_semantics': 'thermochemical frequencies after any HIR projection',
        'thermochemistry': thermochemistry_evidence(species),
    }
    for attribute in ('sigma_ext', 'nopt', 'optical_isomers',
                      'proper_mirror_rmsd'):
        value = getattr(species, attribute, None)
        if value is not None and value != -1:
            data[attribute] = (int(value) if attribute in {
                'nopt', 'optical_isomers'} else float(value))
    data['stereochemical_identity'] = canonical_identity(species)
    if hasattr(species, 'conformer_symmetry'):
        data['conformer_symmetry'] = species.conformer_symmetry
    return json_record(data)


class BimolecularNetwork:
    """Nodes and pathways produced by one bimolecular KinBot run."""

    def __init__(self, species, par):
        self.species = species
        self.par = par
        self.reference = species.energy + species.zpe
        self.nodes = {}
        self.transitions = []
        self.connections = []
        self._build()

    @staticmethod
    def _fragment_key(fragments):
        from kinbot.stereo_routing import refuse_routing
        fragments = list(fragments)
        keys = []
        for fragment in fragments:
            identity = canonical_identity(fragment)
            if identity['status'] != 'assigned':
                refuse_routing(
                    'Canonical bimolecular endpoint identity is unavailable or unsupported: '
                    + identity.get('reason', identity['status']), fragments)
            keys.append(f'{fragment.chemid}:q{fragment.charge}:m{fragment.mult}'
                        + ':stereo=' + identity['id'])
        return tuple(sorted(keys))

    def _relative_energy(self, species):
        return ((species.energy + species.zpe - self.reference)
                * constants.AUtoKCAL)

    def _equivalence_class(self, atom):
        fragments = ((self.species.fragA, 0),
                     (self.species.fragB, self.species.fragA.natom))
        for fragment_index, (fragment, offset) in enumerate(fragments):
            local = atom - offset
            if 0 <= local < fragment.natom:
                group = next(group for group in fragment.atom_eqv
                             if local in group)
                return ((fragment_index,
                         tuple(index + offset for index in group)), len(group))
        raise IndexError(f'Atom {atom} does not belong to either fragment.')

    def _site_degeneracy(self, reactions):
        sites = {}
        for reaction in reactions:
            if hasattr(reaction, 'site_degeneracy'):
                key = (reaction.family_name,
                       getattr(reaction, 'stereopath_id',
                               tuple(reaction.instance)))
                sites.setdefault(key, int(reaction.site_degeneracy))
                continue
            instances = getattr(reaction, 'search_instances',
                                [reaction.instance])
            for instance in instances:
                get_site = getattr(reaction, 'get_bimol_site',
                                   lambda value: value)
                classes = [self._equivalence_class(atom)
                           for atom in get_site(instance)]
                key = tuple(item[0] for item in classes)
                count = int(np.prod([item[1] for item in classes]))
                sites.setdefault((reaction.family_name, key), count)
        return sum(sites.values())

    def _add_bimolecular(self, fragments):
        key = self._fragment_key(fragments)
        node_id = f'b_{_digest("_".join(key))}'
        if node_id not in self.nodes:
            energy = (sum(fragment.energy + fragment.zpe
                          for fragment in fragments) - self.reference)
            self.nodes[node_id] = {
                'id': node_id,
                'kind': 'bimolecular',
                'fragments': list(key),
                'fragment_species': list(fragments),
                'energy': energy * constants.AUtoKCAL,
            }
        return node_id

    def _add_products(self, fragments):
        """Add an observed product as a well or a separated-fragment node."""
        if len(fragments) != 1:
            return self._add_bimolecular(fragments)
        species = fragments[0]
        key = self._fragment_key(fragments)[0]
        node_id = f'w_{_digest(key)}'
        if node_id not in self.nodes:
            self.nodes[node_id] = {
                'id': node_id,
                'kind': 'well',
                'energy': self._relative_energy(species),
                'species': species,
            }
        return node_id

    def _cluster_complexes(self, candidates, role):
        groups = []
        candidates.sort(key=lambda item: (
            item['fragments'], geometry_signature(item['species']),
            item['reaction'].instance_name))
        for candidate in candidates:
            match = None
            for group in groups:
                representative = group[0]
                if candidate['fragments'] != representative['fragments']:
                    continue
                energy_difference = abs(
                    self._relative_energy(candidate['species'])
                    - self._relative_energy(representative['species']))
                if energy_difference > self.par['bimol_energy_tolerance']:
                    continue
                if symmetry_rmsd(candidate['species'],
                                 representative['species']) \
                        <= self.par['bimol_basin_rmsd']:
                    match = group
                    break
            if match is None:
                groups.append([candidate])
            else:
                match.append(candidate)

        used_ids = set(self.nodes)
        created = []
        for group in groups:
            representative = min(
                group, key=lambda item: item['species'].energy
                + item['species'].zpe)
            fragment_hash = _digest('_'.join(representative['fragments']))
            geometry_hash = min(geometry_signature(item['species'])
                                for item in group)
            node_id = unique_observation_id(f'{role}_{fragment_hash}_{geometry_hash}', used_ids)
            self.nodes[node_id] = {
                'id': node_id,
                'kind': 'complex',
                'role': role,
                'fragments': list(representative['fragments']),
                'energy': self._relative_energy(representative['species']),
                'species': representative['species'],
                'members': sorted(item['reaction'].instance_name
                                  for item in group),
                'observations': [{'observation': item['reaction'].instance_name,
                                  'stationary_point': _species_data(item['species'])} for item in group],
                'energy_conflicts': [], 'conflict_evidence': [],
            }
            created.append((self.nodes[node_id], group))
            for item in group:
                setattr(item['reaction'], f'{role}_basin_id', node_id)

        for index, (first, left_group) in enumerate(created):
            for second, right_group in created[:index]:
                if first['fragments'] != second['fragments']:
                    continue
                witnesses = []
                for left in left_group:
                    for right in right_group:
                        delta = abs(self._relative_energy(left['species'])-self._relative_energy(right['species']))
                        if delta <= self.par['bimol_energy_tolerance']:
                            continue
                        value = symmetry_rmsd(left['species'], right['species'])
                        if value <= self.par['bimol_basin_rmsd']:
                            witnesses.append({'observations': [left['reaction'].instance_name,
                                                              right['reaction'].instance_name],
                                              'rmsd_angstrom': value, 'energy_difference_kcal_mol': delta})
                if witnesses:
                    for target, other in ((first, second), (second, first)):
                        target['energy_conflicts'].append(other['id'])
                        target['conflict_evidence'].append({'kind': 'energy', 'other': other['id'],
                                                           'witnesses': witnesses})
                        target['basin_identity_status'] = 'unresolved energy observations'

    def _add_connection(self, first, second, role):
        key = tuple(sorted((first, second)))
        if any(connection['key'] == key for connection in self.connections):
            return
        self.connections.append({
            'id': f'{role}_{_digest("_".join(key))}',
            'kind': 'capture',
            'role': role,
            'reactant': first,
            'product': second,
            'model': None,
            'key': key,
        })

    def _cluster_transitions(self, reactions):
        # Keep the original search target available for sequential retries.
        # Network grouping and site counting describe the bonds actually found.
        observed = []
        for reaction in reactions:
            if getattr(reaction, 'observed_site', None):
                reaction = copy.copy(reaction)
                for key, value in reaction.observed_site.items():
                    setattr(reaction, key, copy.deepcopy(value))
                reaction.search_instances = [list(reaction.instance)]
            observed.append(reaction)
        reactions = observed
        groups = []
        ordered = sorted(reactions, key=lambda reaction: (
            reaction.family_name, geometry_signature(reaction.ts),
            reaction.instance_name))
        for reaction in ordered:
            match = None
            for group in groups:
                representative = group[0]
                if reaction.family_name != representative.family_name:
                    continue
                if not same_stereopath(reaction, representative):
                    continue
                if (reaction.rc_basin_id, reaction.pc_basin_id) != (
                        representative.rc_basin_id, representative.pc_basin_id):
                    continue
                energy_difference = abs(
                    self._relative_energy(reaction.ts)
                    - self._relative_energy(representative.ts))
                if energy_difference > self.par['bimol_energy_tolerance']:
                    continue
                value = symmetry_rmsd(
                    reaction.ts, representative.ts,
                    reaction.ts.geom, representative.ts.geom,
                    reaction.instance, representative.instance,
                    allow_reflection=(
                        getattr(reaction, 'allow_reflection', True)
                        and getattr(representative,
                                    'allow_reflection', True)))
                if value <= self.par['bimol_ts_rmsd']:
                    match = group
                    break
            if match is None:
                groups.append([reaction])
            else:
                match.append(reaction)

        used_ids = {transition['id'] for transition in self.transitions}
        for group in groups:
            representative = min(
                group, key=lambda reaction: reaction.ts.energy
                + reaction.ts.zpe)
            archive, seen = [], set()
            for reaction in sorted(group, key=lambda r: r.instance_name):
                for entry in getattr(reaction, 'final_conformer_archive', []):
                    serialized = json.dumps(json_record(entry), sort_keys=True, allow_nan=False)
                    if serialized not in seen:
                        seen.add(serialized)
                        archive.append(copy.deepcopy(entry))
            endpoints = (representative.rc_basin_id,
                         representative.pc_basin_id)
            family = representative.family_name
            geometry_hash = min(geometry_signature(reaction.ts)
                                for reaction in group)
            site_instances = sorted({tuple(instance) for reaction in group
                                     for instance in getattr(
                                         reaction, 'search_instances',
                                         [reaction.instance])})
            members = {reaction.instance_name for reaction in group}
            for reaction in group:
                members.update(getattr(reaction, 'final_ts_members', []))
                members.update(member['reaction'] for member in getattr(
                    reaction, 'search_members', []))
            stereopath = getattr(representative, 'stereopath_id', '')
            identity = '_'.join(
                endpoints + (family, stereopath, geometry_hash))
            edge_id = unique_observation_id(f'ts_{_digest(identity)}', used_ids)
            self.transitions.append({
                'id': edge_id,
                'kind': 'transition_state',
                'family': family,
                'reactant': endpoints[0],
                'product': endpoints[1],
                'energy': self._relative_energy(representative.ts),
                'degeneracy': self._site_degeneracy(group),
                'site_instances': [list(instance)
                                   for instance in site_instances],
                'stereopath_ids': sorted({
                    reaction.stereopath_id for reaction in group
                    if hasattr(reaction, 'stereopath_id')}),
                'species': representative.ts,
                'reaction': representative,
                'members': sorted(members),
                'product_matches_expected': all(getattr(
                    reaction, 'product_matches_expected', True)
                    for reaction in group),
                'actual_bond_changes': getattr(
                    representative, 'actual_bond_changes', None),
                'parallel_group': _digest('_'.join(endpoints + (family,))),
                'conformer_selection_policy': getattr(representative, 'final_conformer_policy', None),
                'pose_search_policy': getattr(representative, 'pose_search_policy', None),
                'product_spin_selection': getattr(representative, 'product_spin_selection', None),
                'conformer_archive': archive,
                'observations': [{'observation': reaction.instance_name,
                    'reactant': reaction.rc_basin_id, 'product': reaction.pc_basin_id,
                    'stationary_point': _species_data(reaction.ts)} for reaction in group],
                'endpoint_conflicts': [],
                'energy_conflicts': [],
                'conflict_evidence': [],
                '_observations': list(group),
            })

        # Different endpoint observations are preserved, but nearly identical
        # saddles must not automatically become two independent flux terms.
        for index, first in enumerate(self.transitions):
            for second in self.transitions[:index]:
                if first['family'] != second['family']:
                    continue
                evidence = {'endpoint': [], 'energy': []}
                for left in first['_observations']:
                    for right in second['_observations']:
                        if not same_stereopath(left, right):
                            continue
                        value = symmetry_rmsd(left.ts, right.ts, first_instance=left.instance,
                            second_instance=right.instance,
                            allow_reflection=(getattr(left, 'allow_reflection', True)
                                              and getattr(right, 'allow_reflection', True)))
                        if value > self.par['bimol_ts_rmsd']:
                            continue
                        delta = abs(self._relative_energy(left.ts)-self._relative_energy(right.ts))
                        # Different optimized complexes can lead to the same
                        # fragments. Only a different configured reactant or
                        # product identity makes the TS endpoints disagree.
                        kind = ('energy' if delta > self.par['bimol_energy_tolerance'] else
                                'endpoint' if (left.reactant_asymptote_id, left.product_asymptote_id)
                                != (right.reactant_asymptote_id, right.product_asymptote_id) else None)
                        if kind:
                            evidence[kind].append({'observations': [left.instance_name, right.instance_name],
                                'rmsd_angstrom': value, 'energy_difference_kcal_mol': delta})
                for kind, witnesses in evidence.items():
                    if witnesses:
                        first[f'{kind}_conflicts'].append(second['id'])
                        second[f'{kind}_conflicts'].append(first['id'])
                        for target, other in ((first, second), (second, first)):
                            target['conflict_evidence'].append({'kind': kind, 'other': other['id'],
                                                                'witnesses': witnesses})
        for transition in self.transitions:
            transition['flux_assignment'] = ('unresolved observation conflicts; no independent flux assigned'
                if transition['endpoint_conflicts'] or transition['energy_conflicts']
                else 'not yet assigned by a rate model')

    def _build(self):
        reactions = [reaction for index, reaction in
                     enumerate(self.species.reac_obj)
                     if self.species.reac_ts_done[index] == -1]
        reactant_fragments = [self.species.fragA, self.species.fragB]
        reactant_id = self._add_bimolecular(reactant_fragments)
        self.reference_node = reactant_id
        rc_candidates = []
        pc_candidates = []
        for reaction in reactions:
            product_fragments = [opt.species for opt in reaction.prod_opt]
            reaction.reactant_asymptote_id = reactant_id
            reaction.product_asymptote_id = self._add_products(
                product_fragments)
            if reaction.do_reac_vdW:
                species = reaction.irc_reactant_opt.species
                rc_candidates.append({
                    'fragments': self._fragment_key(reactant_fragments),
                    'species': species,
                    'reaction': reaction,
                })
            else:
                reaction.rc_basin_id = reactant_id
            if reaction.do_vdW:
                species = reaction.irc_prod_opt.species
                pc_candidates.append({
                    'fragments': self._fragment_key(product_fragments),
                    'species': species,
                    'reaction': reaction,
                })
            else:
                reaction.pc_basin_id = reaction.product_asymptote_id

        self._cluster_complexes(rc_candidates, 'rc')
        self._cluster_complexes(pc_candidates, 'pc')
        for reaction in reactions:
            if reaction.rc_basin_id != reaction.reactant_asymptote_id:
                self._add_connection(reaction.reactant_asymptote_id,
                                     reaction.rc_basin_id, 'entrance')
            if reaction.pc_basin_id != reaction.product_asymptote_id:
                self._add_connection(reaction.pc_basin_id,
                                     reaction.product_asymptote_id, 'exit')
        self._cluster_transitions(reactions)
        self._select_hir_models()

    def _select_hir_models(self):
        """Select one HIR model per chemical pathway, retaining every observation.

        Use the same barrier comparison as direct MESS and PES. Encounter
        geometries do not define separate chemical reactions. This is a
        single-conformer approximation, not a sum over all ring conformers.
        """
        self.hir_selection = None
        for transition in self.transitions:
            transition.pop('hir_selection', None)
        if not self.par.get('rotor_scan') or self.par.get('multi_conf_tst'):
            return
        groups = {}
        for transition in sorted(self.transitions, key=lambda t: t['reaction'].instance_name):
            reaction = transition['reaction']
            selection = {'selected': False, 'representative': None}
            transition['hir_selection'] = selection
            try:
                endpoints = (reaction.reactant_asymptote_id, reaction.product_asymptote_id)
                if any(node not in self.nodes for node in endpoints):
                    raise ValueError('Configured reactant or product identity is missing.')
                path = reaction_path_id(reaction)
                if path is None or not np.isfinite(transition['energy']):
                    raise ValueError('The TS lacks a pathway class or finite energy.')
            except (AttributeError, ValueError) as error:
                selection['reason'] = str(error)
                logging.getLogger('KinBot').warning(
                    '%s: no HIR rate model selected: %s. Calculation records retained.',
                    reaction.instance_name, error)
                continue
            selection.update(reactant=endpoints[0], product=endpoints[1], stereopath_id=path)
            for group in groups.setdefault(endpoints, []):
                previous = group[0]
                decision = compare_pathways(
                    previous['reaction'].instance_name, previous['hir_selection']['stereopath_id'],
                    previous['energy'], reaction.instance_name, path, transition['energy'])
                if decision == 'distinct':
                    continue
                group.append(transition)
                if decision == 'replace':
                    group.insert(0, group.pop())
                break
            else:
                groups[endpoints].append([transition])
        for pathways in groups.values():
            for group in pathways:
                representative = group[0]['id']
                for transition in group:
                    transition['hir_selection'].update(
                        selected=transition['id'] == representative,
                        representative=representative)
        self.hir_selection = {
            'policy': 'lowest E + ZPE per configured reactant/product pair and TS stereochemical pathway',
            'approximation': 'one selected conformer with HIR; other ring conformations are not summed',
            'selected_transition_ids': sorted(t['id'] for t in self.transitions
                                              if t['hir_selection']['selected']),
        }

    def _site_counting(self, transition):
        reaction = transition['reaction']
        # Rate models select one lowest E + ZPE complex per fragment pair.
        # Other complex energies do not change the equivalent reacting sites.
        conflict = bool(transition.get('endpoint_conflicts') or transition.get('energy_conflicts')
                        or len(transition.get('stereopath_ids', [])) > 1)
        attributes = ('site_members', 'site_relations', 'site_degeneracy', 'stereopath_id')
        expected = json_record({key: getattr(reaction, key, None) for key in attributes})
        expected_class = None
        for observation in [reaction, *transition.get('_observations', [])]:
            supplied = json_record({key: getattr(observation, key, None) for key in attributes})
            checked = site_counting(supplied['site_degeneracy'], supplied['site_members'],
                                    supplied['site_relations'], supplied['stereopath_id'])
            if checked['unlabelled_model_conversion']['status'] == 'unresolved':
                conflict = True
                continue
            # Clustering has already checked the pathway and geometrical match.
            # Equivalent reacting atoms can have different indices in those
            # observations. Compare their site classes, not the atom numbers.
            site_class = (supplied['stereopath_id'], supplied['site_degeneracy'],
                          supplied['site_relations'],
                          {role: len(atoms) for role, atoms in supplied['site_members'].items()})
            if expected_class is None:
                expected_class = site_class
            conflict |= site_class != expected_class
        return site_counting(transition.get('degeneracy'), expected['site_members'],
                             expected['site_relations'], expected['stereopath_id'],
                             conflict=conflict)

    def as_dict(self):
        nodes = []
        for node in sorted(self.nodes.values(), key=lambda item: item['id']):
            data = {key: value for key, value in node.items()
                    if key not in {'species', 'fragment_species'}}
            if node['kind'] in {'complex', 'well'}:
                data['stationary_point'] = _species_data(node['species'])
            else:
                data['fragment_data'] = [
                    _species_data(fragment)
                    for fragment in node['fragment_species']]
            nodes.append(data)
        transitions = []
        for transition in sorted(self.transitions,
                                 key=lambda item: item['id']):
            data = {key: value for key, value in transition.items()
                    if key not in {'species', 'reaction', '_observations'}}
            data['stationary_point'] = _species_data(transition['species'])
            data['site_counting'] = self._site_counting(transition)
            transitions.append(data)
        connections = [
            {key: value for key, value in connection.items() if key != 'key'}
            for connection in sorted(self.connections,
                                     key=lambda item: item['id'])]
        from kinbot.bimol_rate_models import prepare_models
        return json_record({'schema_version': 3, 'energy_unit': 'kcal/mol',
                'rate_model_ready': False,
                'rate_models': prepare_models(self),
                'hir_selection': getattr(self, 'hir_selection', None),
                'reference_node': self.reference_node, 'nodes': nodes,
                'transitions': transitions, 'connections': connections})

    def write_json(self, filename='bimol_network.json'):
        _atomic_text(filename, json.dumps(self.as_dict(), indent=2, allow_nan=False) + '\n')

    def write_unique_ts(self, directory='bimol/unique_ts'):
        """Write one final L1 structure for each clustered transition."""
        directory = Path(directory)
        directory.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(
            prefix=f'.{directory.name}.stage-', dir=directory.parent))
        index = []
        try:
            for transition in sorted(self.transitions,
                                     key=lambda item: item['id']):
                reaction = transition['reaction']
                target = staging / transition['id']
                target.mkdir()
                write(target / 'ts.xyz', Atoms(
                    reaction.ts.atom, positions=reaction.ts.geom),
                    format='xyz')
                irc_name = getattr(
                    reaction, 'final_irc_name', reaction.instance_name)
                for label, direction in (('forward', 'F'),
                                         ('reverse', 'R')):
                    source = Path(f'{irc_name}_IRC_{direction}.xyz')
                    if source.is_file():
                        shutil.copy2(source, target / f'irc_{label}.xyz')
                for label, endpoint in (
                        ('reactant', reaction.irc_reactant),
                        ('product', reaction.irc_prod)):
                    write(target / f'irc_{label}_endpoint.xyz', Atoms(
                        endpoint.atom, positions=endpoint.geom), format='xyz')

                endpoint_metadata = {}
                for label in ('reactant', 'product'):
                    node_id = transition[label]
                    node = self.nodes[node_id]
                    endpoint_metadata[label] = {
                        'node_id': node_id,
                        'kind': node['kind'],
                        'irc_direction': getattr(
                            reaction,
                            ('irc_reactant_direction' if label == 'reactant'
                             else 'irc_prod_direction'), None),
                    }
                    if node['kind'] in {'complex', 'well'}:
                        point = node['species']
                        endpoint_metadata[label]['stationary_point'] = _species_data(point)
                        write(target / f'{label}.xyz', Atoms(
                            point.atom, positions=point.geom), format='xyz')
                    else:
                        endpoint_metadata[label]['fragment_data'] = [
                            _species_data(point) for point in node['fragment_species']]
                        for fragment_index, point in enumerate(
                                node['fragment_species'], 1):
                            write(target / (
                                f'{label}_fragment_{fragment_index}.xyz'),
                                Atoms(point.atom, positions=point.geom),
                                format='xyz')
                metadata = {
                    'schema_version': 3,
                    'id': transition['id'],
                    'level': ('L2' if self.par.get('high_level', 0)
                              else 'L1'),
                    'irc_level': ('L2' if self.par.get('high_level', 0)
                                  else 'L1'),
                    'representative': reaction.instance_name,
                    'instance': list(reaction.instance),
                    'stereopath_id': getattr(
                        reaction, 'stereopath_id', None),
                    'site_degeneracy': getattr(
                        reaction, 'site_degeneracy', None),
                    'site_relations': getattr(
                        reaction, 'site_relations', None),
                    'site_members': getattr(
                        reaction, 'site_members', None),
                    'site_counting': self._site_counting(transition),
                    'allow_reflection': getattr(
                        reaction, 'allow_reflection', True),
                    'selected_ts_job': getattr(
                        reaction, 'ts_source_job',
                        getattr(reaction.ts, 'source_job', None)),
                    'final_irc_job': irc_name,
                    'product_spin_selection': getattr(reaction, 'product_spin_selection', None),
                    'electronic_energy_hartree': float(reaction.ts.energy),
                    'zpe_hartree': float(reaction.ts.zpe),
                    'frequencies_cm-1': [float(value) for value in
                                         reaction.ts.freq],
                    'frequencies_semantics': 'raw harmonic frequencies of the selected calculation',
                    'thermochemistry': thermochemistry_evidence(reaction.ts),
                    'hessian_available': bool(len(getattr(
                        reaction.ts, 'hess', []))),
                    'symmetry': {
                        key: value for key, value in _species_data(
                            reaction.ts).items()
                        if key in {'sigma_ext', 'nopt', 'optical_isomers',
                                   'proper_mirror_rmsd',
                                   'conformer_symmetry'}},
                    'product_matches_expected': getattr(
                        reaction, 'product_matches_expected', True),
                    'actual_bond_changes': getattr(
                        reaction, 'actual_bond_changes', None),
                    'relative_energy_kcal_mol': transition['energy'],
                    'endpoints': endpoint_metadata,
                    'members': transition['members'],
                    'conformer_selection_policy': transition.get('conformer_selection_policy'),
                    'pose_search_policy': transition.get('pose_search_policy'),
                    'conformer_archive': transition.get('conformer_archive', []),
                    'observations': transition.get('observations', []),
                    'endpoint_conflicts': transition.get('endpoint_conflicts', []),
                    'energy_conflicts': transition.get('energy_conflicts', []),
                    'conflict_evidence': transition.get('conflict_evidence', []),
                    'flux_assignment': transition.get('flux_assignment', 'not yet assigned by a rate model'),
                    'hir_selection': transition.get('hir_selection'),
                }
                metadata = json_record(metadata)
                (target / 'metadata.json').write_text(
                    json.dumps(metadata, indent=2, allow_nan=False) + '\n')
                index.append(metadata)
            (staging / 'index.json').write_text(json.dumps(
                {'schema_version': 3, 'unique_ts': index,
                 'hir_selection': getattr(self, 'hir_selection', None)},
                indent=2, allow_nan=False) + '\n')

            backup_root = Path(tempfile.mkdtemp(
                prefix=f'.{directory.name}.backup-',
                dir=directory.parent))
            backup = backup_root / directory.name
            try:
                if directory.exists():
                    directory.replace(backup)
                staging.replace(directory)
            except Exception:
                if not directory.exists() and backup.exists():
                    backup.replace(directory)
                raise
            finally:
                if backup_root.exists():
                    shutil.rmtree(backup_root)
        except Exception:
            if staging.exists():
                shutil.rmtree(staging)
            raise

    def write_pesviewer(self, filename='pesviewer.inp'):
        from kinbot.postprocess import make_xyz

        directory = Path('xyz')
        directory.mkdir(exist_ok=True)
        wells = []
        bimolecular = []
        for node in sorted(self.nodes.values(), key=lambda item: item['id']):
            if node['kind'] in {'complex', 'well'}:
                wells.append(f'{node["id"]} {node["energy"]:.2f}')
                species = node['species']
                make_xyz(species.atom, species.geom, node['id'],
                         str(directory))
            else:
                bimolecular.append(f'{node["id"]} {node["energy"]:.2f}')
                for index, fragment in enumerate(node['fragment_species'], 1):
                    make_xyz(fragment.atom, fragment.geom,
                             f'{node["id"]}{index}', str(directory))
        transition_lines = []
        for transition in sorted(self.transitions,
                                 key=lambda item: item['id']):
            transition_lines.append(
                f'{transition["id"]} {transition["energy"]:.2f} '
                f'{transition["reactant"]} {transition["product"]}')
            species = transition['species']
            make_xyz(species.atom, species.geom, transition['id'],
                     str(directory))
        connection_lines = [
            f'{connection["id"]} {connection["reactant"]} '
            f'{connection["product"]}' for connection in sorted(
                self.connections, key=lambda item: item['id'])]
        with open(f'{kb_path}/tpl/pesviewer.inp.tpl') as handle:
            template = handle.read()
        _atomic_text(filename, template.format(
            wells='\n'.join(wells),
            bimolecs='\n'.join(bimolecular),
            ts='\n'.join(transition_lines),
            barrierless='\n'.join(connection_lines),
            id=self.species.chemid,
        ))
