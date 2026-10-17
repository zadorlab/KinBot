"""Keep a supplied reference identity across bond-order perception.

This does not canonicalize unrelated database entries into one species. It
validates a particular requested reference against its geometry-derived
SMILES, whose localized radical/bond orders may differ on hypervalent atoms.
"""
from rdkit import Chem
import networkx as nx


def reference_smiles_match(requested, perceived):
    """Match explicit-H connectivity, isotope, charge and radical count.

    Exact isomeric SMILES matches remain preferred. Nonidentical stereo
    declarations are not reconciled by this connectivity-only comparison.
    Proton transfers, different hydrogen counts and connectivity changes
    are rejected. Spin multiplicity is checked separately by the caller.
    """
    molecules = [Chem.MolFromSmiles(value) for value in (requested, perceived)]
    if any(mol is None for mol in molecules):
        return False
    if Chem.MolToSmiles(molecules[0]) == Chem.MolToSmiles(molecules[1]):
        return True
    if any(c in value for value in (requested, perceived) for c in '@/\\'):
        return False
    molecules = [Chem.AddHs(mol) for mol in molecules]
    if (Chem.GetFormalCharge(molecules[0]) != Chem.GetFormalCharge(molecules[1])
            or sum(a.GetNumRadicalElectrons() for a in molecules[0].GetAtoms())
            != sum(a.GetNumRadicalElectrons() for a in molecules[1].GetAtoms())):
        return False
    graphs = []
    for mol in molecules:
        graph = nx.Graph()
        for atom in mol.GetAtoms():
            graph.add_node(atom.GetIdx(), element=atom.GetAtomicNum(),
                           isotope=atom.GetIsotope())
        graph.add_edges_from((b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                             for b in mol.GetBonds())
        graphs.append(graph)
    return nx.is_isomorphic(*graphs, node_match=lambda a, b: a == b)
