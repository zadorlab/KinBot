"""Read-only stereochemical helpers for reaction-path bookkeeping.

The temporary atom labels used here never modify the atoms or geometries sent
to an electronic-structure backend.  They only refine reaction-site classes
that KinBot's graph equivalence would otherwise merge.
"""

import copy
import hashlib

import networkx as nx
import numpy as np
import rmsd
from kinbot.stereo_identity import canonical_identity
from kinbot.molecular_symmetry import OPTICAL_RMSD_TOLERANCE


def _component(species, atom):
    """Keep disconnected bimolecular fragments distinct during matching."""
    if not hasattr(species, 'fragA'):
        return 'single'
    fragment = (species.fragA if atom < species.fragA.natom
                else species.fragB)
    return (f'{fragment.chemid}:q{fragment.charge}:m{fragment.mult}')


def _graph(species, tagged_atom=None):
    graph = nx.Graph()
    for index, element in enumerate(species.atom):
        radicals = tuple(sorted(
            int(pattern[index]) for pattern in getattr(species, 'rads', [])))
        graph.add_node(
            index, element=str(element), atomid=int(species.atomid[index]),
            radicals=radicals, component=_component(species, index),
            tagged=(index == tagged_atom))
    bonds = getattr(species, 'bonds', []) or [species.bond]
    for first in range(species.natom):
        for second in range(first):
            if not species.bond01[first][second]:
                continue
            orders = {int(matrix[first][second]) for matrix in bonds
                      if matrix[first][second]}
            if species.bond[first][second]:
                orders.add(int(species.bond[first][second]))
            graph.add_edge(first, second, orders=tuple(sorted(orders)))
    return graph


def _node_match(first, second):
    return (first['element'] == second['element']
            and first['atomid'] == second['atomid']
            and first['radicals'] == second['radicals']
            and first['component'] == second['component']
            and first['tagged'] == second['tagged'])


def _tagged_species(species, atom):
    """Return a copy carrying one virtual isotope-like graph label."""
    # A shallow copy is intentional: OpenBabel objects attached to a
    # StationaryPoint are not pickleable, and calc_chiral only writes the
    # copied ``chiral`` attribute while reading the structural arrays.
    tagged = copy.copy(species)
    tagged.atomid = list(tagged.atomid)
    tagged.atomid[atom] = max(int(value) for value in tagged.atomid) + 1
    tagged.atom_eqv = [list(group) for group in tagged.atom_eqv
                       if atom not in group]
    original_group = next(
        (list(group) for group in species.atom_eqv if atom in group), [atom])
    remainder = [index for index in original_group if index != atom]
    if remainder:
        tagged.atom_eqv.append(remainder)
    tagged.atom_eqv.append([atom])
    tagged.atom_uniq = [group[0] for group in tagged.atom_eqv]
    tagged.calc_chiral()
    return tagged


def stereotopic_relation(species, first, second):
    """Classify two graph-equivalent atoms as homo-, enantio-, or diastereotopic.

    The test uses KinBot's atom identifiers and chirality labels after
    virtually distinguishing each site in turn.  Existing stereocentres are
    therefore retained, while the actual chemical structure is untouched.
    """
    if first == second:
        return 'homotopic'
    left_id = canonical_identity(species, tagged_atom=first)
    right_id = canonical_identity(species, tagged_atom=second)
    if left_id['status'] == right_id['status'] == 'assigned':
        if left_id['id'] == right_id['id']:
            return 'homotopic'
        if left_id['mirror_id'] == right_id['id']:
            return 'enantiotopic'
        return 'diastereotopic'
    # Preserve the pre-existing local chirality refinement when RDKit is not
    # installed; this fallback is not advertised as a canonical identity.
    left = _tagged_species(species, first)
    right = _tagged_species(species, second)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _graph(left, first), _graph(right, second),
        node_match=_node_match,
        edge_match=lambda a, b: a['orders'] == b['orders'])
    enantiotopic = False
    for mapping in matcher.isomorphisms_iter():
        left_chiral = np.sign(np.asarray(left.chiral, dtype=float))
        right_chiral = np.sign(np.asarray(
            [right.chiral[mapping[index]] for index in range(species.natom)],
            dtype=float))
        active = (left_chiral != 0) | (right_chiral != 0)
        if not np.any(active) or np.array_equal(
                left_chiral[active], right_chiral[active]):
            return 'homotopic'
        if np.array_equal(left_chiral[active], -right_chiral[active]):
            enantiotopic = True
    return 'enantiotopic' if enantiotopic else 'diastereotopic'


def refine_equivalence_group(species, group):
    """Split one graph-equivalence group only at diastereotopic boundaries."""
    group = [int(atom) for atom in group]
    if len(group) < 2:
        return [{'members': [atom], 'representative': atom, 'relation': 'homotopic'}
                for atom in group]
    identities = [canonical_identity(species, tagged_atom=atom) for atom in group]
    if all(identity['status'] == 'assigned' for identity in identities):
        # Virtual substitution is needed once per atom, not twice per pair.
        # A mirror family is exactly one homo/enantiotopic class; unrelated
        # families stay diastereotopic. Insertion order preserves representatives.
        families = {}
        for atom, identity in zip(group, identities):
            item = families.setdefault(identity['mirror_family_id'],
                {'members': [], 'identities': set()})
            item['members'].append(atom)
            item['identities'].add(identity['id'])
        return [{'members': item['members'], 'representative': item['members'][0],
                 'relation': ('enantiotopic' if len(item['identities']) > 1 else
                              'diastereotopic' if len(families) > 1 else 'homotopic')}
                for item in families.values()]
    # Retain the existing noncanonical fallback when RDKit cannot assign scope.
    classes = []
    for atom in group:
        for item in classes:
            relations = [stereotopic_relation(species, atom, member)
                         for member in item['members']]
            if all(relation != 'diastereotopic' for relation in relations):
                item['members'].append(atom)
                item['relations'].extend(relations)
                break
        else:
            classes.append({'members': [atom], 'relations': []})
    split = len(classes) > 1
    for item in classes:
        item['representative'] = item['members'][0]
        item['relation'] = ('enantiotopic'
                            if 'enantiotopic' in item['relations']
                            else ('diastereotopic' if split
                                  else 'homotopic'))
        del item['relations']
    return classes


def stereotopic_classes(species, atoms):
    """Return refined classes for selected atoms, preserving original groups."""
    selected = set(int(atom) for atom in atoms)
    classes = []
    for group_index, group in enumerate(species.atom_eqv):
        subset = [int(atom) for atom in group if atom in selected]
        if not subset:
            continue
        for class_index, item in enumerate(
                refine_equivalence_group(species, subset)):
            item.update({'group_index': group_index,
                         'class_index': class_index})
            classes.append(item)
    return classes


def stereotopic_site_signature(species, atom):
    """Return a permutation- and mirror-invariant label for one site."""
    identity = canonical_identity(species, tagged_atom=atom)
    if identity['status'] == 'assigned':
        return identity['mirror_family_id']
    tagged = _tagged_species(species, atom)
    graph = _graph(tagged, atom)
    distances = nx.single_source_shortest_path_length(graph, atom)

    def descriptors(mirror):
        values = []
        for index, data in graph.nodes(data=True):
            handedness = int(np.sign(tagged.chiral[index])) * mirror
            values.append((
                int(distances.get(index, species.natom + 1)),
                data['component'], data['element'], str(data['atomid']),
                data['radicals'], data['tagged'], handedness))
        return tuple(sorted(values))

    # A global reflection changes every handedness, but does not create a new
    # reaction-site label. Relative handedness still distinguishes
    # diastereotopic sites.
    canonical = min(descriptors(1), descriptors(-1))
    return hashlib.sha1(repr(canonical).encode()).hexdigest()[:12]



def virtually_labelled(species, labels):
    """Return a temporary isotope-labelled view; QC atoms remain untouched."""
    view = copy.copy(species)
    view.isotopes = list(getattr(species, 'isotopes', [0] * species.natom))
    offset = max(view.isotopes + [1000]) + 1
    for atom, role in labels.items():
        view.isotopes[int(atom)] = offset + int(role)
    return view


def motif_identity(species, motif):
    """Identify an ordered joint selection, retaining its relative stereo."""
    view = virtually_labelled(species, {atom: role for role, atom in enumerate(motif)})
    return canonical_identity(view)


def configuration_erased_graph(identity):
    """Discard stereo while retaining virtual labels and graph connectivity."""
    from rdkit import Chem
    graphs = []
    for text in identity['canonical_graphs']:
        mol = Chem.MolFromSmiles(text, sanitize=False)
        Chem.RemoveStereochemistry(mol)
        graphs.append(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
    return tuple(sorted(graphs))


def reaction_atom_equivalence(species):
    """Refine an entire reaction motif without changing molecular atom_eqv.

    A heavy-atom branch can be pruned before its terminal transferred atom is
    visited. Use the same substitution test for every atom, including groups
    transferred by intra_R_migration, additions and bond scission.
    """
    return [item['members'] for group in species.atom_eqv
            for item in refine_equivalence_group(species, group)]


def stereotopic_hydrogen_equivalence(species):
    """Refine only H equivalence groups for selected H-transfer finders."""
    refined = []
    for group in species.atom_eqv:
        if (len(group) < 2
                or any(species.atom[atom] != 'H' for atom in group)):
            refined.append(list(group))
            continue
        refined.extend(item['members']
                       for item in refine_equivalence_group(species, group))
    return refined


def proper_mirror_rmsd(species, geom=None):
    """Return the best proper-rotation RMSD to the structure's mirror image."""
    left = np.asarray(species.geom if geom is None else geom, dtype=float)
    mirrored = left * np.array([-1., 1., 1.])
    from kinbot.molecular_symmetry import proper_rmsd, spatial_mappings
    best = np.inf
    # Use the same chemistry/reaction-role mappings as rate counting. An
    # infinite search bound requests the minimum, rather than a yes/no test.
    for mapping in spatial_mappings(species, left, mirrored, tolerance=np.inf):
        best = min(best, proper_rmsd(left, mirrored[mapping]))
        if best < 1.e-5:
            break
    return float(best)


def optical_isomers(species, tolerance=OPTICAL_RMSD_TOLERANCE):
    """Count the mirror pair represented by one stationary geometry."""
    from kinbot.optical import rigid_mirror
    mirror_rmsd = proper_mirror_rmsd(species)
    return rigid_mirror(species, tolerance=tolerance)['mirror_states'], mirror_rmsd


def assign_optical_metadata(species, tolerance=OPTICAL_RMSD_TOLERANCE):
    """Record optical information without changing legacy ``nopt``."""
    count, mirror_rmsd = optical_isomers(species, tolerance=tolerance)
    species.optical_isomers = count
    species.proper_mirror_rmsd = mirror_rmsd
    return {
        'sigma_ext': float(species.sigma_ext),
        'nopt': int(species.nopt),
        'optical_isomers': species.optical_isomers,
        'proper_mirror_rmsd': float(mirror_rmsd),
    }


def conformer_symmetry(species, tolerance=OPTICAL_RMSD_TOLERANCE):
    """Export legacy rotational numbers and separate optical observations."""
    from kinbot.symmetry import conformer_symmetry_numbers
    from kinbot.molecular_symmetry import geometric_mirror_states

    records = []
    geometries = getattr(species, 'conformer_geom', [])
    indices = getattr(species, 'conformer_index', [])
    for offset, geom in enumerate(geometries):
        index = int(indices[offset] if offset < len(indices) else offset)
        if index < 0:
            continue
        numbers = conformer_symmetry_numbers(species, geom)
        mirrors = geometric_mirror_states(species, geom, tolerance)
        data = dict(index=index, sigma_ext=numbers['sigma_ext'],
                    nopt=numbers['nopt'], mirror_states=mirrors,
                    optical_isomers=mirrors)
        records.append(data)
    return records


def same_stereopath(first, second):
    """Allow legacy unlabeled reactions to cluster, but protect new labels."""
    left = getattr(first, 'stereopath_id', None)
    right = getattr(second, 'stereopath_id', None)
    if left is None and right is None:
        return True
    return left == right
