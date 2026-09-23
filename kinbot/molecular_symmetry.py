"""Spatial atom mappings for duplicate removal and optical coverage.

Graph permutations are candidate atom mappings, never symmetry numbers by
themselves. Distance bounds prune the search before proper Kabsch alignment.
This module does not assign rotational symmetry numbers or pathway degeneracy.
Rotational symmetry is assigned by the legacy rules in kinbot.symmetry.
"""
import networkx as nx
import numpy as np
from ase.data import atomic_numbers, covalent_radii


OPTICAL_RMSD_TOLERANCE = .1


def reaction_bond_changes(species):
    """Signed endpoint bond changes for a TS, or None for an ordinary species."""
    changes = getattr(species, 'reac_bond', None)
    if not getattr(species, 'wellorts', 0) or changes is None:
        return None
    changes = np.asarray(changes)
    if changes.shape != (len(species.atom), len(species.atom)):
        raise ValueError('TS reacting-bond matrix must contain every atom.')
    return changes


def chemical_graph(species):
    atoms = list(map(str, species.atom))
    n = len(atoms)
    bond = getattr(species, 'bond', None)
    if bond is None:
        geom = np.asarray(species.geom)
        radii = np.array([covalent_radii[atomic_numbers[a]] for a in atoms])
        distance = np.linalg.norm(geom[:, None] - geom[None, :], axis=-1)
        bond = ((distance < 1.25 * (radii[:, None] + radii[None, :]))
                & (distance > .1)).astype(int)
    matrices = list(getattr(species, 'bonds', [])) or [bond]
    radicals = list(getattr(species, 'rads', []))
    isotopes = getattr(species, 'isotopes', [0] * n)
    charges = getattr(species, 'formal_charges', [0] * n)
    changes = reaction_bond_changes(species)
    graph = nx.Graph()
    for index, atom in enumerate(atoms):
        graph.add_node(index, label=(atom, int(isotopes[index]), int(charges[index]),
                       tuple(sorted(int(rad[index]) for rad in radicals))))
    for i in range(n):
        for j in range(i):
            orders = tuple(sorted(int(matrix[i][j]) for matrix in matrices))
            if any(orders) or bond[i][j] or (changes is not None and changes[i, j]):
                # A chosen Kekule structure must not distinguish symmetry-
                # equivalent atoms. Keep TS union-only edges distinct, and
                # preserve literal orders separately for molecule builders.
                label = (bool(bond[i][j]) and not any(orders), orders)
                if changes is not None:
                    label += (int(changes[i, j]),)
                graph.add_edge(i, j, label=label, order=int(bond[i][j]))
    return graph


def _physical_role_graphs(graph, species):
    """Physical TS symmetry may reverse the entire reaction direction.

    A proper rotation of H--H--H exchanges its forming and breaking bonds.
    Permit that global exchange, never independent swaps at individual bonds.
    The original directed graph and reaction-path metadata remain unchanged.
    """
    yield graph
    changes = reaction_bond_changes(species)
    if changes is not None and np.any(changes):
        reversed_graph = graph.copy()
        for _, _, edge in reversed_graph.edges(data=True):
            label = edge.get('label')
            if label is not None:
                edge['label'] = (*label[:-1], -label[-1])
        yield reversed_graph


def graph_equivalence_classes(species):
    """Atom orbits preserving graph labels, including forming/breaking bonds.

    Rooted graph matching stops at the first compatible mapping for each pair;
    it does not enumerate every permutation of equivalent methyl hydrogens.
    """
    graph = chemical_graph(species)
    labels = [-1] * len(species.atom)
    for atom in graph:
        if labels[atom] >= 0:
            continue
        labels[atom] = atom
        left = graph.copy()
        nx.set_node_attributes(left, False, 'root')
        left.nodes[atom]['root'] = True
        for other in range(atom + 1, len(labels)):
            if labels[other] >= 0 or graph.nodes[atom]['label'] != graph.nodes[other]['label']:
                continue
            right = graph.copy()
            nx.set_node_attributes(right, False, 'root')
            right.nodes[other]['root'] = True
            for target in _physical_role_graphs(right, species):
                matcher = nx.algorithms.isomorphism.GraphMatcher(
                    left, target,
                    node_match=lambda a, b: (a['label'], a['root']) == (b['label'], b['root']),
                    edge_match=lambda a, b: a['label'] == b['label'])
                if any(_preserves_resonance(species, [mapping[i] for i in range(len(labels))])
                       for mapping in matcher.isomorphisms_iter()):
                    labels[other] = atom
                    break
    return labels


def reaction_atom_labels(species):
    """Refine legacy atom classes only for TSs with endpoint bond changes."""
    if reaction_bond_changes(species) is None:
        return species.atomid
    orbits = graph_equivalence_classes(species)
    classes = {}
    return [classes.setdefault((atomid, orbit), len(classes))
            for atomid, orbit in zip(species.atomid, orbits)]


def _preserves_resonance(species, order):
    matrices = list(getattr(species, 'bonds', []))
    radicals = list(getattr(species, 'rads', []))
    if not matrices:
        return True
    # Per-edge order sets are only a pruning aid: preserve their correlations.
    states = {tuple(np.asarray(matrix).ravel()) +
              tuple(radicals[i] if i < len(radicals) else [])
              for i, matrix in enumerate(matrices)}
    mapped = {tuple(np.asarray(matrix)[np.ix_(order, order)].ravel()) +
              tuple(np.asarray(radicals[i])[order] if i < len(radicals) else [])
              for i, matrix in enumerate(matrices)}
    return mapped == states


def proper_rmsd(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    left = left - left.mean(axis=0)
    right = right - right.mean(axis=0)
    u, _, vt = np.linalg.svd(left.T @ right)
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(u @ vt)
    return float(np.sqrt(np.mean(np.sum((left @ u @ correction @ vt-right)**2, axis=1))))


def spatial_mappings(species, left, right, tolerance=.05):
    """Yield chemically compatible maps that may align within RMSD tolerance."""
    base = chemical_graph(species)
    graphs = []
    for coordinates in (left, right):
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape != (len(species.atom), 3) or not np.all(np.isfinite(coordinates)):
            raise ValueError('Symmetry requires finite coordinates for every atom.')
        graph = nx.Graph()
        graph.add_nodes_from(base.nodes(data=True))
        for i in range(len(coordinates)):
            for j in range(i):
                graph.add_edge(i, j, label=base.get_edge_data(i, j, {}).get('label'),
                               distance=float(np.linalg.norm(coordinates[i]-coordinates[j])))
        graphs.append(graph)
    # Any alignment with RMSD <= tolerance satisfies this pair-distance bound.
    bound = 2 * np.sqrt(len(species.atom)) * tolerance + 1.e-10
    for target in _physical_role_graphs(graphs[1], species):
        matcher = nx.algorithms.isomorphism.GraphMatcher(graphs[0], target,
            node_match=lambda a, b: a['label'] == b['label'],
            edge_match=lambda a, b: a['label'] == b['label']
            and abs(a['distance']-b['distance']) <= bound)
        for mapping in matcher.isomorphisms_iter():
            order = [mapping[index] for index in range(len(species.atom))]
            if not _preserves_resonance(species, order):
                continue
            yield order


def equivalent_geometry(species, left, right, tolerance=.05, reflected=False):
    """Duplicate/coverage check under proper rotations and chemical mappings."""
    right = np.asarray(right, dtype=float)
    if reflected:
        right = right * [-1., 1., 1.]
    return any(proper_rmsd(left, right[mapping]) <= tolerance
               for mapping in spatial_mappings(species, left, right, tolerance))


def geometric_mirror_states(species, geom=None, tolerance=OPTICAL_RMSD_TOLERANCE):
    """Identify a geometry's mirror pair, independently of rotational sigma."""
    from kinbot.optical import rigid_mirror
    return rigid_mirror(species, geom, tolerance)['mirror_states']
