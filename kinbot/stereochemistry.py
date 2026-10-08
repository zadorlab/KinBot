"""Read-only stereochemical helpers for reaction-path bookkeeping.

The temporary atom labels used here never modify the atoms or geometries sent
to an electronic-structure backend.  They only refine reaction-site classes
that KinBot's graph equivalence would otherwise merge.
"""

import copy
import numpy as np
from kinbot.stereo_identity import require_supported_identity
from kinbot.molecular_symmetry import OPTICAL_RMSD_TOLERANCE


def stereotopic_relation(species, first, second):
    """Compare virtual substitutions at two graph-equivalent atoms."""
    if first == second:
        return 'homotopic'
    left = require_supported_identity(species, tagged_atom=first)
    right = require_supported_identity(species, tagged_atom=second)
    if left['id'] == right['id']:
        return 'homotopic'
    if left['mirror_id'] == right['id']:
        return 'enantiotopic'
    return 'diastereotopic'


def refine_equivalence_group(species, group):
    """Split one graph-equivalence group only at diastereotopic boundaries."""
    group = [int(atom) for atom in group]
    if len(group) < 2:
        return [{'members': [atom], 'representative': atom, 'relation': 'homotopic'}
                for atom in group]
    identities = [require_supported_identity(species, tagged_atom=atom) for atom in group]
    # Each mirror family contains the homotopic/enantiotopic alternatives.
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
    return require_supported_identity(species, tagged_atom=atom)['mirror_family_id']



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
    return require_supported_identity(view)


def configuration_erased_graph(identity):
    """Discard stereo while retaining virtual labels and graph connectivity."""
    from rdkit import Chem
    graphs = []
    for text in identity['canonical_graphs']:
        mol = Chem.MolFromSmiles(text, sanitize=False)
        Chem.RemoveStereochemistry(mol)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)  # identity-only axial orientation
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
