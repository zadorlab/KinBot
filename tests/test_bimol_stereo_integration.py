"""Exercise the shared stereochemistry checks through bimolecular handoffs.

These are deterministic structural fixtures, not optimized saddle points.
No QC calculations are submitted.
"""
import json
from types import SimpleNamespace

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from kinbot import bimol, postprocess
from kinbot.counting_contract import counting_view
from kinbot.parameters import Parameters
from kinbot.reaction_finder_bimol import ReactionFinderBimol
from kinbot.reaction_path import (path_geometry_allowed, reaction_path_id,
                                  read_summary_paths)
from kinbot.reactions.reac_abstraction import Abstraction
from kinbot.species_routing import routing_key
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity


def donor_and_hydrogen(smiles='CO'):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=312) == 0
    donor = StationaryPoint('donor', 0, 1,
        atom=[atom.GetSymbol() for atom in mol.GetAtoms()],
        geom=mol.GetConformer().GetPositions())
    donor.characterize()
    hydrogen = StationaryPoint('H', 0, 2, atom=['H'], geom=np.array([[10., 0., 0.]]))
    hydrogen.characterize()
    species = StationaryPoint('reactants', 0, 2,
        atom=list(donor.atom) + ['H'],
        geom=np.vstack([donor.geom, hydrogen.geom]), fragA=donor, fragB=hydrogen)
    species.characterize()
    for point in (donor, hydrogen, species):
        point.energy, point.zpe = -1., .01
    species.energy, species.zpe = -2., .02
    return species


def abstraction(species, donor=0, hydrogen=None):
    if hydrogen is None:
        hydrogen = next(i for i, symbol in enumerate(species.atom)
                        if symbol == 'H' and species.bond[donor, i])
    acceptor = species.natom - 1
    product_geom = species.geom.copy()
    product_geom[hydrogen] = product_geom[acceptor] + [0., .74, 0.]
    product = StationaryPoint('product', 0, 2, atom=species.atom, geom=product_geom)
    product.characterize()
    reaction = Abstraction(species, None, {'scan_step': 20},
                           [acceptor, hydrogen, donor], f'abstraction_{hydrogen}')
    reaction.stereopath_id = 'enumerated-site-class'
    reaction.l0_validation = {'endpoints': [
        {'state': state, 'configured_identity': canonical_identity(point),
         'geometry_angstrom': point.geom.tolist(),
         'bonds': np.argwhere(np.triu(point.bond01, 1)).tolist()}
        for state, point in (('reactant', species), ('product', product))]}
    row = SimpleNamespace(name=reaction.instance_name, id=1, positions=species.geom,
        data={'energy': -50., 'zpe': .02,
              'frequencies': [-800.] + [500.] * (3 * species.natom - 7),
              'hess': np.eye(3 * species.natom)})
    reaction.qc = SimpleNamespace(qc='fc',
        db=SimpleNamespace(select=lambda **kwargs: iter([row])),
        get_qc_geom=lambda *args, **kwargs: (0, species.geom.copy()))
    reaction.irc_reactant, reaction.irc_prod = species, product
    return reaction


def test_saved_endpoints_protect_ts_and_survive_summary_and_network(tmp_path, monkeypatch):
    species = donor_and_hydrogen()
    reaction = abstraction(species)
    ts = bimol.build_bimolecular_ts(reaction, reaction.qc)
    assert ts.stereopath_metadata['broken_bonds'] == [[0, reaction.instance[1]]]
    assert path_geometry_allowed(ts, ts.geom)
    assert reaction.stereopath_id == 'enumerated-site-class'
    assert reaction_path_id(reaction) == 'ordinary'

    # Final IRC observations replace the discovery evidence before export.
    bimol.update_ts_topology_from_irc(reaction)
    assert reaction_path_id(reaction) == 'ordinary'
    reaction.products = reaction.irc_prod.start_multi_molecular()[0]
    for point in reaction.products:
        point.energy, point.zpe = -1.1, .01
    reaction.prod_opt = [SimpleNamespace(species=point) for point in reaction.products]
    reaction.do_reac_vdW = reaction.do_vdW = False
    species.reac_obj, species.reac_ts_done = [reaction], [-1]
    species.reac_inst, species.reac_name = [reaction.instance], [reaction.instance_name]
    species.reac_type = ['abstraction']
    par = {'bimol': 1, 'high_level': 0, 'optical_population': 'specified',
           'bimol_energy_tolerance': .25, 'bimol_basin_rmsd': .15, 'bimol_ts_rmsd': .1}
    monkeypatch.chdir(tmp_path)
    postprocess.create_summary_file(species, reaction.qc, par)
    assert species.reac_ts_done == [-1]
    summary = (tmp_path / f'summary_{routing_key(species)}.out').read_text()
    assert 'SUCCESS' in summary
    assert read_summary_paths(summary.splitlines()) == {reaction.instance_name: 'ordinary'}
    postprocess.createPESViewerInput(species, reaction.qc, par)
    network = json.loads((tmp_path / 'bimol_network.json').read_text())
    assert len(network['transitions']) == 1
    assert network['transitions'][0]['stereopath_ids'] == ['enumerated-site-class']


def test_bimolecular_diastereotopic_sites_keep_distinct_ts_guards(tmp_path):
    species = donor_and_hydrogen('CC[C@H](O)C')
    hydrogens = [i for i, symbol in enumerate(species.atom)
                 if symbol == 'H' and species.bond[1, i]]
    assert len(hydrogens) == 2
    input_path = tmp_path / 'input.json'
    input_path.write_text(json.dumps({'barrier_threshold': 50.}))
    finder = ReactionFinderBimol(species, Parameters(input_path).par, None)
    site_keys, path_keys = [], []
    for hydrogen in hydrogens:
        reaction = abstraction(species, donor=1, hydrogen=hydrogen)
        reaction.stereopath_id = finder._site_signature(hydrogen)
        ts = bimol.build_bimolecular_ts(reaction, reaction.qc)
        site_keys.append(reaction.stereopath_id)
        path_keys.append(reaction_path_id(reaction))
        swapped = ts.geom.copy()
        swapped[hydrogens] = swapped[hydrogens[::-1]]
        # A working TS may lead to another product configuration. The final
        # IRC assigns that channel; rate counting keeps its products separate.
        assert path_geometry_allowed(ts, swapped)
        assert not path_geometry_allowed(counting_view(ts), swapped)
        assert path_geometry_allowed(ts, ts.geom)
    assert len(set(site_keys)) == len(set(path_keys)) == 2


def test_missing_l0_evidence_does_not_create_an_unprotected_ts():
    reaction = abstraction(donor_and_hydrogen())
    del reaction.l0_validation
    with pytest.raises(ValueError, match='Missing configured L0 endpoints'):
        bimol.build_bimolecular_ts(reaction, reaction.qc)


def test_bimol_site_comparison_keeps_distinct_and_unknown_classes_separate():
    first, second = SimpleNamespace(), SimpleNamespace()
    assert bimol.same_stereopath(first, second)
    first.stereopath_id = 'site-A'
    assert not bimol.same_stereopath(first, second)
    second.stereopath_id = 'site-B'
    assert not bimol.same_stereopath(first, second)
