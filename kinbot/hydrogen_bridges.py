"""Three-centre C--H--C connections in positive-ion well geometries.

The connectivity graph retains both contacts. Ordinary valence descriptions
attach H to one carbon at a time, so reaction motifs never give H two single
bonds. Geometry recognition does not establish that a structure is a minimum.
"""
import itertools

import networkx as nx
import numpy as np


def find_bridges(species, bond):
    """Recognize compact, shared H contacts without extending ordinary bonds.

    Both C--H contacts must be elongated: a normal C--H bond followed by an
    internal hydrogen bond is not a shared proton. Carbon atoms can be across
    a ring or in different fragments; an existing C--C bond is not required.
    """
    known = {tuple(bridge) for bridge in getattr(species, 'hydrogen_bridges', ())}
    transition_state = getattr(species, 'wellorts', 0)
    if species.charge <= 0 or transition_state and not known:
        return []
    bridges = []
    for hydrogen, element in enumerate(species.atom):
        if element != 'H':
            continue
        if any(bond[hydrogen, other] and symbol != 'C'
               for other, symbol in enumerate(species.atom)):
            continue
        contacts = [carbon for carbon, symbol in enumerate(species.atom)
                    if symbol == 'C' and species.dist[hydrogen, carbon] < 1.5]
        if len(contacts) != 2:
            continue
        first, second = contacts
        if transition_state and (first, hydrogen, second) not in known:
            continue
        lengths = species.dist[hydrogen, contacts]
        if min(lengths) < 1.16 or abs(lengths[0]-lengths[1]) > .20:
            continue
        # H must lie between the carbons along their connecting direction;
        # proximity to two carbons on the same side is not a bridge.
        direction = species.geom[second]-species.geom[first]
        projection = np.dot(species.geom[hydrogen]-species.geom[first], direction)
        if not 0. < projection < np.dot(direction, direction):
            continue
        bridges.append((first, hydrogen, second))
    return bridges


def spectator_bridges(reactant, product):
    """Keep C-H-C contacts confirmed at both endpoints of a reaction."""
    confirmed = {tuple(bridge) for bridge in getattr(product, 'hydrogen_bridges', ())}
    return [tuple(bridge) for bridge in getattr(reactant, 'hydrogen_bridges', ())
            if tuple(bridge) in confirmed
            and all(reactant.bond[carbon, bridge[1]] == product.bond[carbon, bridge[1]] > 0
                    for carbon in (bridge[0], bridge[2]))]


def copy_ts_bridges(source, target):
    """Do not rediscover H-transfer TS contacts as minimum bridges on a copy."""
    if source.charge > 0 and getattr(source, 'wellorts', 0):
        target.wellorts = source.wellorts
        target.hydrogen_bridges = [tuple(bridge) for bridge in getattr(source, 'hydrogen_bridges', ())]


def localized_bonds(bond, bridges):
    """Yield the two localized attachments of each shared hydrogen."""
    for choices in itertools.product((0, 1), repeat=len(bridges)):
        matrix = bond.copy()
        for (first, hydrogen, second), choice in zip(bridges, choices):
            matrix[hydrogen, :] = matrix[:, hydrogen] = 0
            carbon = (first, second)[choice]
            matrix[hydrogen, carbon] = matrix[carbon, hydrogen] = 1
        # Invalid contacts must not manufacture hypervalent carbon states.
        yield matrix


def bridge_edges(species):
    return {frozenset((carbon, hydrogen))
            for first, hydrogen, second in getattr(species, 'hydrogen_bridges', ())
            for carbon in (first, second)}


def rotation_axis(species, first, second):
    """A rotor must separate the physical connection when its axis is cut."""
    graph = nx.from_numpy_array(np.asarray(species.bond) > 0)
    if not graph.has_edge(first, second):
        return False
    graph.remove_edge(first, second)
    return not nx.has_path(graph, first, second)


def rotational_symmetry(species):
    """Count proper rotations of a bridged geometry using existing mappings."""
    from kinbot.molecular_symmetry import spatial_mappings, proper_rmsd

    return sum(proper_rmsd(species.geom, species.geom[order]) <= .05
               for order in spatial_mappings(species, species.geom, species.geom))
