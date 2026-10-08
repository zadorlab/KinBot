"""Axial assignments on disposable identity graphs, never on QC atoms.

RDKit supplies the atom equivalences and SMILES ordering. Its ordinary 3-D
assignment does not encode allene handedness. Identity-only atom-map labels
record that missing orientation while retaining the complete molecular graph.
The signs are not CIP R/S descriptors. They are invariant to atom numbering
and proper rotation; reflection reverses the axial signs.
"""
import numpy as np
import networkx as nx


# Coplanar terminal groups have no resolved axial orientation. Five degrees
# also prevents small out-of-plane residuals in planar radical resonance
# drawings from becoming a new axial configuration. This is not a barrier.
PLANAR_SINE = np.sin(np.deg2rad(5.))


def assign_axial_stereo(mol, geom, aromatic_forms, ignored_atoms=(), ignored_bonds=()):
    """Assign cumulene and potential biaryl orientations from explicit atoms.

    Even cumulenes carry handedness; odd cumulenes carry terminal E/Z order.
    Potential biaryls use the same substituted-ring criterion as the former
    scope refusal. Their rotation barrier is not inferred: specified mode
    retains one orientation, whereas racemic mode permits its global mirror.
    A symmetry-equivalent pair of terminal arms does not define an axis.
    """
    from rdkit import Chem

    probe = Chem.Mol(mol)
    Chem.SetAromaticity(probe)
    xyz = np.asarray(geom, float)
    ignored_atoms = set(map(int, ignored_atoms))
    ignored_bonds = {frozenset(map(int, pair)) for pair in ignored_bonds}
    doubles = nx.Graph()
    for bond in probe.GetBonds():
        if bond.GetBondTypeAsDouble() == 2:
            doubles.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    candidates = []
    for component in nx.connected_components(doubles):
        ends = [i for i in component if doubles.degree(i) == 1]
        if (len(component) < 3 or len(ends) != 2
                or any(probe.GetAtomWithIdx(i).GetDegree() != 2
                       for i in component if i not in ends)):
            continue
        path = nx.shortest_path(doubles, *ends)
        edges = {frozenset(pair) for pair in zip(path, path[1:])}
        if any(edges <= form for form in aromatic_forms()):
            continue
        # Individual C=C tags are not terminal cumulene stereochemistry.
        for a, b in zip(path, path[1:]):
            bond = mol.GetBondBetweenAtoms(a, b)
            bond.SetStereo(Chem.BondStereo.STEREONONE)
            bond.SetBondDir(Chem.BondDir.NONE)
        axial = len(path) % 2 == 1
        middle = [path[len(path)//2]] if axial else path[len(path)//2-1:len(path)//2+1]
        candidates.append((path, middle, 'cumulene' if axial else 'cumulene_ez'))

    for bond in probe.GetBonds():
        ends = (bond.GetBeginAtom(), bond.GetEndAtom())
        if (bond.GetBondType() == Chem.BondType.SINGLE and not bond.IsInRing()
                and all(atom.IsInRing() for atom in ends)
                and all(any(b.GetIsAromatic() or b.GetBondTypeAsDouble() == 2
                            for b in atom.GetBonds()) for atom in ends)):
            path = [a.GetIdx() for a in ends]
            candidates.append((path, path, 'biaryl'))

    # Assign independent axes first. Their labels can distinguish the two
    # arms of another axis. Add each round together and keep earlier signs:
    # an axis must not change meaning when another axis is encountered.
    assigned = set()
    while True:
        ranks = list(Chem.CanonicalRankAtoms(probe, breakTies=False, includeChirality=True))
        assignments = []
        for number, (path, centers, kind) in enumerate(candidates):
            if number in assigned:
                continue
            edges = {frozenset(pair) for pair in zip(path, path[1:])}
            if ignored_atoms.intersection(path) or edges.intersection(ignored_bonds):
                continue
            arms = [[a.GetIdx() for a in probe.GetAtomWithIdx(i).GetNeighbors()
                     if a.GetIdx() not in path] for i in (path[0], path[-1])]
            if any(len(side) != 2 or ranks[side[0]] == ranks[side[1]] for side in arms):
                continue
            # Both arms fix the terminal plane. Their difference also handles a
            # bent cumulene without interpreting its central bond angle as chirality.
            ordered = [sorted(side, key=lambda i: ranks[i]) for side in arms]
            axis = xyz[path[-1]] - xyz[path[0]]
            length = np.linalg.norm(axis)
            if length < 1.e-8:
                continue
            axis /= length
            vectors = [xyz[high] - xyz[low] for low, high in ordered]
            vectors = [v - np.dot(v, axis)*axis for v in vectors]
            norms = [np.linalg.norm(v) for v in vectors]
            if min(norms) < 1.e-8:
                continue
            left, right = [v/length for v, length in zip(vectors, norms)]
            if kind == 'cumulene_ez':
                value = np.dot(left, right)
                # Perpendicular terminal planes do not define terminal E/Z.
                if abs(value) <= PLANAR_SINE:
                    continue
                label = 3 if value > 0 else 4
            else:
                value = np.dot(axis, np.cross(left, right))
                if abs(value) <= PLANAR_SINE:
                    continue
                label = (1 if value > 0 else 2) + (4 if kind == 'biaryl' else 0)
            assignments.append((number, centers, label))
        if not assignments:
            break
        for number, centers, label in assignments:
            assigned.add(number)
            for index in centers:
                mol.GetAtomWithIdx(index).SetAtomMapNum(label)
                probe.GetAtomWithIdx(index).SetAtomMapNum(label)
    return bool(assigned)
