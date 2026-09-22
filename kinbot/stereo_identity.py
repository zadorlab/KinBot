"""Canonical identities for supplied molecular graphs with explicit 3-D stereo.

The identity covers RDKit's tetrahedral and double-bond stereochemistry. It
does not identify transition states, torsional basins or atropisomers. The
charge and multiplicity belong to the identity. Formal charges are included
when supplied, but electronic localization is not inferred. The canonical
strings encode configured graphs, not necessarily sanitized Lewis structures.
No coordinates or graph arrays are modified.
"""
import hashlib
import copy
import logging
import numpy as np


def legacy_stereo_warning(species, reason=None):
    """Report unsupported stereo once while retaining the legacy treatment.

    The identity remains unsupported: a warning must not manufacture an R/S
    assignment or claim that a requested racemate was completely represented.
    """
    reason = reason or canonical_identity(species).get('reason', 'unknown configuration scope')
    if getattr(species, 'stereo_fallback_reason', None) != reason:
        logging.getLogger('KinBot').warning(
            '%s: advanced stereochemistry is unavailable (%s); using the legacy '
            'symmetry treatment. Complete racemic counting is not established.',
            getattr(species, 'name', 'species'), reason)
    species.stereo_fallback_reason = reason
    species.stereo_routing_status = 'legacy; unsupported stereochemistry; unverified'


def _three_ring_benzenoid(mol):
    """Recognize the unsubstituted anthracene/phenanthrene graph class.

    Ordinary bending does not require a new configured identity for these
    skeletons. This narrow exemption is not a general helicene classifier.
    """
    rings = [set(ring) for ring in mol.GetRingInfo().AtomRings()]
    if len(rings) != 3 or any(len(ring) != 6 for ring in rings):
        return False
    shared = [rings[i] & rings[j] for i in range(3) for j in range(i)]
    if sorted(map(len, shared)) != [0, 2, 2]:
        return False
    if any(atoms and mol.GetBondBetweenAtoms(*sorted(atoms)) is None for atoms in shared):
        return False
    core = set.union(*rings)
    if len(core) != 14 or mol.GetNumAtoms() != 24:
        return False
    for atom in mol.GetAtoms():
        orders = sorted(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
        if atom.GetIdx() in core:
            if atom.GetAtomicNum() != 6 or orders != [1., 1., 2.]:
                return False
        elif (atom.GetAtomicNum() != 1 or orders != [1.]
              or atom.GetNeighbors()[0].GetIdx() not in core):
            return False
    return True


def _check_supported_configuration(mol):
    """Conservative boundary for fixed stereo outside the supported tags.

    This detects potential axes/frameworks; it does not infer their barriers
    or assert that every flagged structure is a stable atropisomer.
    """
    import networkx as nx
    from rdkit import Chem
    ranks = Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=False)
    doubles = nx.Graph()
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 2:
            doubles.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    for component in nx.connected_components(doubles):
        ends = [i for i in component if doubles.degree(i) == 1]
        substituents = [[ranks[a.GetIdx()] for a in mol.GetAtomWithIdx(i).GetNeighbors()
                         if a.GetIdx() not in component] for i in ends]
        if (len(component) > 2 and len(ends) == 2
                and all(len(set(group)) >= 2 for group in substituents)):
            raise ValueError('cumulene/axial configuration is outside the supported stereo scope')
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 4 or (atom.GetDegree() and atom.GetAtomicNum()
                                   not in {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}):
            raise ValueError('non-tetrahedral/coordination configuration is outside the supported stereo scope')
    for bond in mol.GetBonds():
        ends = (bond.GetBeginAtom(), bond.GetEndAtom())
        if (not bond.IsInRing() and all(atom.IsInRing() for atom in ends)
                and all(any(b.GetBondTypeAsDouble() == 2 for b in atom.GetBonds())
                        for atom in ends)
                and all(len({ranks[a.GetIdx()] for a in atom.GetNeighbors()
                             if a.GetIdx() not in (ends[0].GetIdx(), ends[1].GetIdx())}) >= 2
                        for atom in ends)):
            raise ValueError('potential biaryl atropisomer requires an explicit configuration model')
    pi_rings = sum(any(mol.GetBondWithIdx(i).GetBondTypeAsDouble() == 2 for i in ring)
                   for ring in mol.GetRingInfo().BondRings())
    if pi_rings >= 3 and not _three_ring_benzenoid(mol):
        raise ValueError('polycyclic pi-framework configuration requires a helical/planar stereo model')


def _molecules(species, geom, tagged_atom=None):
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    from kinbot.molecular_symmetry import chemical_graph

    n = len(species.atom)
    bonds = list(getattr(species, 'bonds', []))
    if not bonds:
        graph = chemical_graph(species)
        bond = np.zeros((n, n), dtype=int)
        for i, j, data in graph.edges(data=True):
            bond[i, j] = bond[j, i] = data['order']
        bonds = [bond]
    isotopes = list(getattr(species, 'isotopes', [0] * n))
    charges = list(getattr(species, 'formal_charges', [0] * n))
    if tagged_atom is not None:
        isotopes[tagged_atom] = max(isotopes + [1000]) + 1
    result = []
    for matrix in bonds:
        mol = Chem.RWMol()
        for i, element in enumerate(species.atom):
            atom = Chem.Atom(str(element))
            atom.SetNoImplicit(True)
            atom.SetIsotope(int(isotopes[i]))
            atom.SetFormalCharge(int(charges[i]))
            # KinBot's rads are valence remainders, not assigned unpaired
            # electrons: a three-bonded O has -1; an alkoxide O has +1.
            # Neither is an RDKit radical count. Canonicalize the supplied
            # graph without guessing local charges or electron populations.
            mol.AddAtom(atom)
        for i in range(n):
            for j in range(i):
                order = int(matrix[i][j])
                if order:
                    mol.AddBond(i, j, {1: Chem.BondType.SINGLE,
                                     2: Chem.BondType.DOUBLE,
                                     3: Chem.BondType.TRIPLE}[order])
        mol = mol.GetMol()
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
        _check_supported_configuration(mol)
        conf = Chem.Conformer(n)
        conf.Set3D(True)
        for i, xyz in enumerate(geom):
            conf.SetAtomPosition(i, Point3D(*map(float, xyz)))
        mol.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(mol, confId=0, replaceExistingTags=True)
        # TS spectator guards use endpoint graphs while allowing reacting
        # centres to lose or change configuration. This is an in-memory view.
        ignored = set(getattr(species, 'stereo_ignored_atoms', ()))
        for index in ignored:
            mol.GetAtomWithIdx(int(index)).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        ignored_bonds = {frozenset(map(int, pair))
                         for pair in getattr(species, 'stereo_ignored_bonds', ())}
        for bond in mol.GetBonds():
            if (bond.GetBeginAtomIdx() in ignored or bond.GetEndAtomIdx() in ignored
                    or frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())) in ignored_bonds):
                bond.SetStereo(Chem.BondStereo.STEREONONE)
                bond.SetBondDir(Chem.BondDir.NONE)
        for pair in ignored_bonds:
            for index in pair:
                for adjacent in mol.GetAtomWithIdx(index).GetBonds():
                    adjacent.SetBondDir(Chem.BondDir.NONE)
        # Canonical SMILES normalizes atom numbering, not the supplied set of
        # resonance structures. Retain the entire provided ensemble explicitly.
        result.append(mol)
    return result


def _strings(species, geom, tagged_atom=None):
    from rdkit import Chem
    return tuple(sorted({Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True, allHsExplicit=True)
                         for mol in _molecules(species, geom, tagged_atom)}))


def endpoint_configuration_allowed(endpoint, reference_geom, observed_geom):
    """Compare configurations present in an endpoint and the discovery TS.

    Atom order and graph are fixed here. A centre becoming planar may lose its
    tag; it is not required to retain a fictitious R/S assignment. A definite
    opposite assignment still belongs to another configured endpoint channel.
    """
    actual = _molecules(endpoint, endpoint.geom)
    reference = _molecules(endpoint, reference_geom)
    observed = _molecules(endpoint, observed_geom)
    for stable, before, after in zip(actual, reference, observed):
        for a, b, c in zip(stable.GetAtoms(), before.GetAtoms(), after.GetAtoms()):
            tags = tuple(int(atom.GetChiralTag()) for atom in (a, b, c))
            if all(tags) and tags[1] != tags[2]:
                return False
        for a, b, c in zip(stable.GetBonds(), before.GetBonds(), after.GetBonds()):
            tags = tuple(int(bond.GetStereo()) for bond in (a, b, c))
            if all(tags) and tags[1] != tags[2]:
                return False
    return True


def canonical_identity(species, geom=None, *, tagged_atom=None):
    """Return a scoped stereo identity, or an explicit unsupported status."""
    if getattr(species, 'wellorts', 0):
        return {'status': 'unsupported', 'reason': 'transition-state graph'}
    coordinates = np.asarray(species.geom if geom is None else geom, dtype=float)
    if coordinates.shape != (len(species.atom), 3) or not np.all(np.isfinite(coordinates)):
        raise ValueError('Canonical stereo identity requires finite atom coordinates.')
    try:
        own = _strings(species, coordinates, tagged_atom)
        mirror = _strings(species, coordinates * [-1., 1., 1.], tagged_atom)
    except ImportError:
        return {'status': 'unavailable', 'reason': 'install kinbot[stereo] (RDKit)'}
    except (ValueError, KeyError, RuntimeError) as error:
        return {'status': 'unsupported', 'reason': str(error)}
    charge = int(getattr(species, 'charge', 0))
    multiplicity = int(getattr(species, 'mult', 1))
    def key(value):
        return hashlib.sha256(repr((value, charge, multiplicity)).encode()).hexdigest()
    return {'status': 'assigned', 'schema': 'kinbot.stereo.v1',
            'scope': 'tetrahedral and double-bond; provided resonance ensemble',
            'id': key(own), 'mirror_id': key(mirror),
            'mirror_family_id': key(min(own, mirror)),
            'is_chiral_configuration': own != mirror,
            'canonical_graphs': own, 'charge': charge, 'multiplicity': multiplicity,
            'electronic_localization': 'not inferred from valence remainders',
            'formal_charge_localization': ('specified' if hasattr(species, 'formal_charges')
                                           else 'not supplied')}


def optical_scope(species, population='specified'):
    """Whether a global mirror is in this species' declared population.

    Saddles require compatible populations on both endpoints. Unknown
    identity does not authorize adding an unobserved enantiomer.
    """
    if population not in ('specified', 'racemic'):
        raise ValueError('optical_population must be specified or racemic')
    reference = getattr(species, 'optical_reference', None)
    if reference is None:
        # Freeze the declared input configuration before a selected geometry
        # changes. TS callers supply the reactant reference explicitly.
        reference = canonical_identity(species)
        species.optical_reference = copy.deepcopy(reference)
    identity = reference
    known = identity.get('status') == 'assigned'
    mirror_allowed = bool(known and (population == 'racemic'
                                    or not identity['is_chiral_configuration']))
    endpoints = getattr(species, 'ts_endpoint_identities', ())
    if endpoints and population == 'specified':
        from collections import Counter
        known_endpoints = all(item.get('status') == 'assigned' for side in endpoints for item in side)
        if known_endpoints:
            own = [tuple(sorted(Counter(item['id'] for item in side).items())) for side in endpoints]
            mirrored = [tuple(sorted(Counter(item['mirror_id'] for item in side).items())) for side in endpoints]
            # An elementary channel is undirected. Reflection may preserve
            # each side or exchange R -> S with S -> R within the same channel.
            mirror_allowed = known and sorted(own) == sorted(mirrored)
        else:
            mirror_allowed = False
    product = getattr(species, 'stereopath_product_identity', None)
    if product is not None and not endpoints and population == 'specified':
        # For an explicitly classified route, a global mirror reaching a
        # different configured product belongs to that product's channel.
        mirror_allowed &= (product.get('status') == 'assigned'
                           and not product['is_chiral_configuration'])
    return {'population': population, 'identity': identity,
            'mirror_allowed': mirror_allowed}


def configured_geometry_allowed(species, geom, population=None):
    """Check stable configuration, or an explicitly annotated TS pathway."""
    if getattr(species, 'wellorts', 0):
        from kinbot.reaction_path import path_geometry_allowed
        return path_geometry_allowed(species, geom, population)
    population = population or getattr(species, 'optical_population', 'specified')
    reference = optical_scope(species, population)['identity']
    observed = canonical_identity(species, geom)
    matches = identity_matches(reference, observed, population)
    if matches is None:
        # Unknown scope is rejected by strict MC counting, not by ordinary
        # connectivity validation in legacy non-MC runs.
        return True
    return matches


def identity_matches(reference, observed, population='specified'):
    """Compare configured identities; None means the assignment is unknown.

    A racemic population contains the global mirror, not other diastereomers.
    Callers retain their existing fallback when either assignment is unknown.
    """
    if reference.get('status') != 'assigned' or observed.get('status') != 'assigned':
        return None
    allowed = {reference['id']}
    if population == 'racemic':
        allowed.add(reference['mirror_id'])
    return observed['id'] in allowed
