"""Retain a connected incidental product without inventing an abstraction rate."""
import copy
from types import SimpleNamespace

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from kinbot import bimol
from kinbot.bimol_network import BimolecularNetwork
from kinbot.reactions.reac_abstraction import Abstraction
from kinbot.stationary_pt import StationaryPoint


def _properties(point, energy):
    point.energy, point.zpe = energy, .01
    point.freq = [500.] * (3 * point.natom - 6)
    point.reduced_freqs = list(point.freq)
    point.sigma_ext = 1


def test_ethyl_recombination_remains_a_well_with_no_h_abstraction_rate():
    # These constructed structures exercise bookkeeping only. No saddle point
    # or barrier for radical recombination is established by this test.
    molecule = Chem.AddHs(Chem.MolFromSmiles('CCCC'))
    assert AllChem.EmbedMolecule(molecule, randomSeed=124) == 0
    AllChem.UFFOptimizeMolecule(molecule)
    disconnected = Chem.RWMol(molecule)
    disconnected.RemoveBond(1, 2)
    groups = Chem.GetMolFrags(disconnected.GetMol())
    assert len(groups) == 2
    order = [atom for group in groups for atom in group]
    atoms = [molecule.GetAtomWithIdx(atom).GetSymbol() for atom in order]
    geometry = molecule.GetConformer().GetPositions()[order]
    product = StationaryPoint('butane', 0, 1, atom=atoms, geom=geometry)
    product.characterize()
    size = len(groups[0])
    first = StationaryPoint('ethyl_a', 0, 2, atom=atoms[:size], geom=geometry[:size])
    second = StationaryPoint('ethyl_b', 0, 2,
        atom=atoms[size:], geom=geometry[size:] + [10., 0., 0.])
    for fragment in (first, second):
        fragment.characterize()
        _properties(fragment, -79.)
    species = StationaryPoint('two_ethyl_singlet', 0, 1, atom=atoms,
        geom=np.vstack((first.geom, second.geom)), fragA=first, fragB=second)
    species.characterize()
    species.energy, species.zpe = -158., .02
    donor = order.index(1)
    abstractor = order.index(2)
    hydrogens = [atom for atom in range(size)
                 if atoms[atom] == 'H' and species.bond01[donor, atom]]
    hydrogen = hydrogens[0]
    reaction = Abstraction(species, None, {}, [abstractor, hydrogen, donor],
        'ethyl_h_abstraction_found_recombination')
    reaction.site_degeneracy = 2
    reaction.site_members = {'abstractor': [abstractor], 'transferred_atom': hydrogens}
    reaction.site_relations = {'abstractor': 'singleton', 'transferred_atom': 'enantiotopic'}
    reaction.stereopath_id = 'requested_h_transfer'
    reaction.ts = copy.deepcopy(product)
    reaction.ts.name, reaction.ts.wellorts = 'constructed_ts', 1
    reaction.irc_reactant, reaction.irc_prod = species, product
    reaction.product_matches_expected = False
    reaction.prod_opt = [SimpleNamespace(species=product)]
    reaction.do_reac_vdW = reaction.do_vdW = False
    species.reac_obj, species.reac_ts_done = [reaction], [-1]
    _properties(product, -159.)
    _properties(reaction.ts, -157.9)
    reaction.ts.freq[0] = -400.
    reaction.ts.reduced_freqs = list(reaction.ts.freq)
    bimol.update_ts_topology_from_irc(reaction)

    network = BimolecularNetwork(species, dict(rotor_scan=1, multi_conf_tst=0,
        optical_population='specified', bimol_energy_tolerance=.25,
        bimol_ts_rmsd=.1, bimol_basin_rmsd=.15))
    data = network.as_dict()
    assert len(data['transitions']) == 1
    transition = data['transitions'][0]
    assert not transition['product_matches_expected']
    assert transition['family'] == 'bimolecular' and transition['degeneracy'] == 1
    product_node = next(node for node in data['nodes']
                        if node['id'] == transition['product'])
    assert product_node['kind'] == 'well'
    assert product_node['stationary_point']['name'] == 'butane'
    assert data['rate_models']['pathways'] == []
    assert any('two separated fragments on each side' in warning
               for warning in data['rate_models']['warnings'])
    assert reaction.instance == [abstractor, hydrogen, donor]
    assert reaction.site_degeneracy == 2
