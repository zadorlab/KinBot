"""Optical counting in the represented rigid-body/torsional model.

One global chemical atom mapping must match every anchored rigid part. Only
accepted bridge rotors supply independent motions. This is a coordinate-model
test, not a calculation of inversion dynamics or vibrational overlap.
"""
from itertools import combinations
import logging
import numpy as np
import networkx as nx

from kinbot import constants
from kinbot.molecular_symmetry import (chemical_graph, _physical_role_graphs,
                                      _preserves_resonance, OPTICAL_RMSD_TOLERANCE)
from kinbot.stereo_identity import optical_scope, canonical_identity

METHOD = 'anchored-parts-v1'
# Measured scan witnesses must still identify the same calculated structure.
# Keep this limit separate from the part-model tolerance: loosening that model
# must not turn nearby barrier geometries into measured reference/mirror minima.
SCAN_MIRROR_RMSD_TOLERANCE = .1
logger = logging.getLogger('KinBot')


def unresolved_optical_default(result):
    """Keep geometric uncertainty explicit while permitting weight-one output."""
    if (result.get('status') == 'unresolved' and result.get('reason') in (
            'Rigid mirror comparison is numerically undetermined.',
            'The represented motions do not establish mirror coverage or retained handedness.',
            'Uncertain explicit mirror coverage.',
            'Multiple relaxed rotor coordinates independently contain the same mirror; overlapping coverage is unresolved.',
            'The mirror witness coincides with the reference under the assigned rotor periods; coverage is not established.')):
        return dict(result, remaining_multiplier=1., fallback='unresolved_symmetry')
    return result


def _fit_errors(left, right):
    left = left - left.mean(axis=0)
    right = right - right.mean(axis=0)
    u, _, vt = np.linalg.svd(left.T @ right)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    return np.linalg.norm(left @ u @ correction @ vt - right, axis=1)


def _fit(left, right):
    errors = _fit_errors(left, right)
    return float(np.sqrt(np.mean(errors**2))), float(max(errors, default=0.))


class _Parts:
    def __init__(self, species, rotors):
        self.species = species
        self.graph = chemical_graph(species)
        cut = self.graph.copy()
        self.axes = {}
        bridges = {frozenset(edge) for edge in nx.bridges(cut)}
        for rotor in rotors:
            axis = frozenset(rotor['axis'])
            domain = rotor.get('represented_domain_degrees')
            sigma = rotor.get('sigma_int')
            if (axis not in bridges or axis in self.axes or not rotor.get('usable')
                    or not domain or domain[0] != 0 or not 0 < domain[1] <= 360
                    or sigma is None or sigma <= 0):
                raise ValueError('Optical coverage requires distinct usable bridge-rotor domains.')
            self.axes[axis] = (float(sigma), float(domain[1]))
        cut.remove_edges_from(tuple(axis) for axis in self.axes)
        # Separate fragments still belong to ONE rigid-body model. Their
        # relative pose cannot be changed by independently aligning them.
        if not nx.is_connected(self.graph) and self.axes:
            raise ValueError('Disconnected assemblies require an explicit relative-motion model.')
        components = list(nx.connected_components(cut)) if self.axes else [set(cut)]
        self.blocks = []
        for component in components:
            anchored = set(component)
            for axis in self.axes:
                if axis & component:
                    anchored.update(axis)
            self.blocks.append(tuple(sorted(anchored)))
        # Small bonded neighbourhoods protect local handedness from dilution
        # by distant atoms. No enumeration of arbitrary molecular quadruples.
        self.cores = []
        self.neighborhoods = {}
        for block in self.blocks:
            self.cores.append(block)
            self.neighborhoods[block] = []
            for atom in block:
                local = tuple(sorted({atom} | (set(self.graph[atom]) & set(block))))
                self.neighborhoods[block].append(local)
                if 4 <= len(local) < len(block) and local not in self.cores:
                    self.cores.append(local)
            for atom in block:
                for path in self._four_atom_paths(atom, set(block)):
                    local = tuple(sorted(path))
                    if local not in self.neighborhoods[block]:
                        self.neighborhoods[block].append(local)

    def compare(self, left, right, tolerance, *, reflected=True, observations=None):
        """Return match/distinct/undetermined, with evidence under one map.

        Use one proper alignment per anchored part, then require the RMSD
        on every atom's bonded neighbourhood to pass. A distant tail cannot
        dilute a local error. Independent minimum RMSDs are lower bounds;
        a failed positive fit without a lower-bound proof stays undetermined.
        """
        n = len(self.species.atom)
        left, right = np.asarray(left, float), np.asarray(right, float)
        if any(g.shape != (n, 3) or not np.all(np.isfinite(g)) for g in (left, right)):
            raise ValueError('Optical counting requires finite coordinates for every atom.')
        if reflected:
            right = right * [-1., 1., 1.]
        obs = None if observations is None else [np.asarray(g, float) for g in observations]
        if obs is not None and any(g.shape != (n, 3) or not np.all(np.isfinite(g)) for g in obs):
            raise ValueError('A successful HIR point lacks finite geometry.')
        pairs = {frozenset(pair) for block in self.blocks for pair in combinations(block, 2)}
        graphs = []
        for geom in (left, right):
            graph = nx.Graph()
            graph.add_nodes_from(self.graph.nodes(data=True))
            for i in range(n):
                for j in range(i):
                    pair = frozenset((i, j))
                    graph.add_edge(i, j, label=self.graph.get_edge_data(i, j, {}).get('label'),
                        axis=self.axes.get(pair),
                        distance=float(np.linalg.norm(geom[i]-geom[j])) if pair in pairs else None)
            graphs.append(graph)

        bound = 2*np.sqrt(max(len(group) for groups in self.neighborhoods.values()
                             for group in groups))*tolerance + 1.e-10
        def edge_match(a, b):
            return (a['label'] == b['label'] and a['axis'] == b['axis']
                    and ((a['distance'] is None and b['distance'] is None)
                         or (a['distance'] is not None and b['distance'] is not None
                             and (obs is not None or abs(a['distance']-b['distance']) <= bound))))

        witnesses = []
        # A local orientation obstruction can prune the rest of a mapping
        # before permutations of remote equivalent methyl hydrogens expand it.
        def obstructed(core, mapping):
            a, b = left[list(core)], right[[mapping[i] for i in core]]
            if _fit(a, b)[0] <= tolerance:
                return False
            if obs is None:
                return True
            # Relaxation is allowed. Require a local signed-volume witness,
            # not identical shapes, in all retained scan observations.
            # Neighbourhoods have <=5 atoms in the ordinary valence cases;
            # large ring blocks are tested through consecutive bonded paths.
            quartets = combinations(core, 4) if len(core) <= 5 else (
                path for start in core for path in self._four_atom_paths(start, set(core)))
            for quartet in quartets:
                quartet = tuple(quartet)
                x = left[list(quartet)]
                y = right[[mapping[i] for i in quartet]]
                volume = np.linalg.det(x[1:]-x[0])
                if volume * np.linalg.det(y[1:]-y[0]) >= 0 or _fit(x, y)[0] <= tolerance:
                    continue
                values = [g[list(quartet)] for g in obs]
                if all(np.linalg.det(v[1:]-v[0])*volume >
                       1.e-8 * max(np.linalg.norm(v[:, None]-v[None, :], axis=-1).max(), 1.)**3
                       * abs(volume) for v in values):
                    example = dict(atom_indices=list(quartet),
                                   mapped_atom_indices=[mapping[i] for i in quartet])
                    if example not in witnesses:
                        witnesses.append(example)
                    return True
            return False

        owner = self
        class Matcher(nx.algorithms.isomorphism.GraphMatcher):
            def semantic_feasibility(self, a, b):
                if not super().semantic_feasibility(a, b):
                    return False
                mapping = dict(self.core_1)
                mapping[a] = b
                return not any(a in core and all(i in mapping for i in core)
                               and obstructed(core, mapping) for core in owner.cores)

        ambiguous = False
        best = float('inf')
        for target in _physical_role_graphs(graphs[1], self.species):
            matcher = Matcher(graphs[0], target,
                node_match=lambda a, b: a['label'] == b['label'], edge_match=edge_match)
            for mapping in matcher.isomorphisms_iter():
                order = [mapping[i] for i in range(n)]
                if not _preserves_resonance(self.species, order):
                    continue
                if any(tuple(sorted(mapping[i] for i in block)) not in self.blocks for block in self.blocks):
                    continue
                maximum, local_error = 0., 0.
                for block in self.blocks:
                    errors = _fit_errors(left[list(block)], right[[mapping[i] for i in block]])
                    maximum = max(maximum, float(max(errors)))
                    local_error = max(local_error, float(np.sqrt(np.mean(errors**2))))
                    for local in self.neighborhoods[block]:
                        local_error = max(local_error, float(np.sqrt(np.mean(
                            errors[[block.index(i) for i in local]]**2))))
                best = min(best, local_error)
                if local_error <= tolerance:
                    return dict(status='match', atom_mapping=order,
                                maximum_displacement_angstrom=maximum,
                                largest_local_rmsd_angstrom=local_error,
                                rigid_parts=[list(b) for b in self.blocks])
                ambiguous = True
        return dict(status='undetermined' if ambiguous else 'distinct',
                    largest_local_rmsd_angstrom=best if np.isfinite(best) else None,
                    rejected_mapping_examples=witnesses,
                    rigid_parts=[list(b) for b in self.blocks])

    def _four_atom_paths(self, start, allowed):
        def extend(path):
            if len(path) == 4:
                yield tuple(path)
            else:
                for node in self.graph[path[-1]]:
                    if node in allowed and node not in path:
                        yield from extend(path + [node])
        yield from extend([start])


def rigid_mirror(species, geometry=None, tolerance=OPTICAL_RMSD_TOLERANCE):
    geom = np.asarray(species.geom if geometry is None else geometry, float)
    result = compare_rigid(species, geom, geom, tolerance)
    result.update(method=METHOD, mirror_states={'match': 1, 'distinct': 2}.get(result['status']))
    return result


def compare_rigid(species, left, right, tolerance=OPTICAL_RMSD_TOLERANCE, *, reflected=True):
    """Protect assigned configurations before approximate geometric matching."""
    own, other = canonical_identity(species, left), canonical_identity(species, right)
    if (own.get('status') == other.get('status') == 'assigned'
            and own['id'] != other['mirror_id' if reflected else 'id']):
        return dict(status='distinct', reason='Different assigned configurations.')
    return _Parts(species, []).compare(left, right, tolerance, reflected=reflected)


def evaluate_optical(species, *, geometry=None, rotors=(), population=None,
                     tolerance=OPTICAL_RMSD_TOLERANCE, energy_tolerance=1.):
    """One remaining optical decision for a rigid conformer or HIR model."""
    scope = optical_scope(species, population or getattr(species, 'optical_population', 'specified'))
    geom = np.asarray(species.geom if geometry is None else geometry, float)
    if scope['identity'].get('status') != 'assigned':
        from kinbot.stereo_identity import legacy_stereo_warning
        from kinbot.symmetry import conformer_symmetry_numbers
        legacy_stereo_warning(species, scope['identity'].get('reason', 'unknown reference'))
        legacy = conformer_symmetry_numbers(species, geom)['nopt']
        return dict(method=METHOD, status='legacy_unverified', population_scope=scope,
                    total_optical_states=None, allowed_global_mirror_states=None,
                    remaining_multiplier=float(legacy),
                    reason='Stereochemical population not established; retained the legacy optical convention.')
    rigid = rigid_mirror(species, geom, tolerance)
    size = rigid['mirror_states']
    result = dict(method=METHOD, status='unresolved', population_scope=scope,
                  total_optical_states=size, allowed_global_mirror_states=None,
                  remaining_multiplier=None, rigid_mirror=rigid,
                  reason='Rigid mirror comparison is numerically undetermined.')
    allowed = size if scope['mirror_allowed'] else 1
    result['allowed_global_mirror_states'] = allowed
    if not rotors and (size == 1 or not scope['mirror_allowed']):
        result.update(status='resolved', remaining_multiplier=1., reason='No additional allowed mirror.')
        return result
    if not rotors:
        if size == 2:
            result.update(status='resolved', remaining_multiplier=2., reason='One rigid geometry omits its allowed mirror.')
        return result
    observations, omitted = [], []
    measured = None
    for rotor in rotors:
        points = [p for p in rotor['points'] if p['status'] == 'successful']
        for p in rotor['points']:
            if p['status'] != 'successful':
                omitted.append(dict(rotor_index=rotor['index'], point_index=p['index']))
        if not points:
            result['reason'] = 'A represented rotor has no successful observations.'
            return result
        for p in points:
            g, e = p.get('geometry_angstrom'), p.get('electronic_energy_hartree')
            if (g is None or np.shape(g) != geom.shape or not np.all(np.isfinite(g))
                    or e is None or not np.isfinite(e)):
                result['reason'] = 'A successful HIR point lacks finite geometry or energy.'
                return result
            observations.append(g)
        measured_tolerance = min(tolerance, SCAN_MIRROR_RMSD_TOLERANCE)
        anchors = [p for p in points if compare_rigid(species, geom, p['geometry_angstrom'], measured_tolerance,
                                                     reflected=False)['status'] == 'match']
        partners = [p for p in points if compare_rigid(species, geom, p['geometry_angstrom'], measured_tolerance)['status'] == 'match']
        coverage = rotor['mirror_coverage'] = dict(status='unresolved', covered=None, witnesses=[],
            mirror_observed_in_full_scan=bool(partners) and size == 2,
            counting_domain_degrees=[0., 360.], potential_period_degrees=rotor['represented_domain_degrees'][1],
            symmetry_quotient=rotor['sigma_int'], symmetry_periodicity_verified=False)
        if size == 2 and partners and not scope['mirror_allowed']:
            result['reason'] = 'HIR contains a global mirror outside the specified population.'
            return result
        if size != 2:
            continue
        for a in anchors:
            for b in partners:
                delta = abs(a['electronic_energy_hartree']-b['electronic_energy_hartree'])*constants.AUtoKCAL
                if delta >= energy_tolerance:
                    message = (f'Measured HIR mirror geometries for {species.name}, rotor '
                               f'{rotor["index"]}, points {a["index"]}/{b["index"]} differ by '
                               f'{delta:g} kcal/mol. Retaining the measured potential; '
                               'this energy difference does not create another optical state.')
                    logger.warning(message)
                    result.setdefault('warnings', []).append(message)
                period = rotor['represented_domain_degrees'][1]
                separation = (b['angle_offset_degrees']-a['angle_offset_degrees']) % period
                if min(separation, period-separation) < 1.e-6:
                    coverage.setdefault('unresolved_coordinate_witnesses', []).append(
                        {'reference_point': a['index'], 'mirror_point': b['index'],
                         'reason': 'These observations share the same represented rotor coordinate.'})
                    continue  # another measured pair or the part comparison may establish coverage
                if all(0 <= p['angle_offset_degrees'] < 360 for p in (a, b)):
                    witness = dict(rotor_index=rotor['index'], reference_point=a['index'],
                                    mirror_point=b['index'], energy_difference_kcal_mol=delta)
                    measured = measured or witness
                    coverage['witnesses'].append(dict(witness, within_serialized_potential_period=
                        all(0 <= p['angle_offset_degrees'] < period for p in (a, b))))
                    coverage.update(status='observed_pair', covered=2)
    if size == 1 or not scope['mirror_allowed']:
        result.update(status='resolved', remaining_multiplier=1., states_covered_by_hir=1,
                      reason='No additional allowed mirror.')
        return result
    if measured:
        if sum(r['mirror_coverage'].get('status') == 'observed_pair' for r in rotors) > 1:
            result['reason'] = 'Multiple relaxed rotor coordinates independently contain the same mirror; overlapping coverage is unresolved.'
            return result
        result.update(status='assumed' if result.get('warnings') else 'resolved',
                      remaining_multiplier=1., states_covered_by_hir=2,
                      measured_coverage=measured, reason='The accepted scan contains the selected structure and its mirror.')
        return result
    if rigid.get('reason') == 'Different assigned configurations.':
        reference = canonical_identity(species, geom)
        identities = [canonical_identity(species, g) for g in observations]
        if all(i.get('status') == 'assigned' and i['id'] == reference['id'] for i in identities):
            result.update(status='resolved', remaining_multiplier=2., states_covered_by_hir=1,
                          reason='The represented torsions retain the assigned fixed configuration.')
        else:
            result['reason'] = 'Scan configuration changes do not establish coverage of the selected mirror.'
        return result
    try:
        parts = _Parts(species, rotors)
    except ValueError as error:
        if str(error) != 'Disconnected assemblies require an explicit relative-motion model.':
            raise
        result.update(reason='The represented motions do not establish mirror coverage or retained handedness.',
                      coordinate_coverage={'status': 'unsupported', 'reason': str(error)})
        return result
    comparison = parts.compare(geom, geom, tolerance, observations=observations)
    result['coordinate_coverage'] = dict(comparison, additional_qc_calculations=0,
        coupled_potential_accuracy_established=False, successful_observations=len(observations),
        omitted_points=omitted)
    if comparison['status'] == 'match':
        from ase import Atoms
        source = Atoms(species.atom, positions=geom)
        target = Atoms(species.atom, positions=(geom * [-1., 1., 1.])[comparison['atom_mapping']])
        changes = []
        for rotor in rotors:
            dihed = rotor['dihedral']
            angle = (target.get_dihedral(*dihed)-source.get_dihedral(*dihed)) % 360.
            period = rotor['represented_domain_degrees'][1]
            changes.append(dict(rotor_index=rotor['index'], angle_offset_degrees=angle,
                                potential_period_degrees=period))
        result['coordinate_coverage']['coordinate_changes'] = changes
        if all(min(step['angle_offset_degrees'] % step['potential_period_degrees'],
                   step['potential_period_degrees'] - step['angle_offset_degrees'] % step['potential_period_degrees'])
               < 1.e-6 for step in changes):
            result['reason'] = 'The mirror witness coincides with the reference under the assigned rotor periods; coverage is not established.'
            return result
        result.update(status='resolved', remaining_multiplier=1., states_covered_by_hir=size,
                      reason='The represented independent torsions include the mirror.')
    elif comparison['status'] == 'distinct' and comparison['rejected_mapping_examples']:
        result.update(status='resolved', remaining_multiplier=2., states_covered_by_hir=1,
                      reason='Retained rigid-part handedness is absent from the represented torsional motions.')
    else:
        result['reason'] = 'The represented motions do not establish mirror coverage or retained handedness.'
    return result


def bind_assumption(species, parameters):
    """Use an exact public state name; never a QC filename or a wildcard."""
    from kinbot.species_routing import routing_name
    name = str(species.name) if getattr(species, 'wellorts', 0) else routing_name(species)
    assumption = parameters.get('optical_factor_assumptions', {}).get(name)
    species.optical_factor_assumption = (dict(assumption, state=name) if assumption else None)
