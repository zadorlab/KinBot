"""Preserve TS conformer guards and the discovery-only unimolecular IRC policy."""
import copy
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from ase.build import molecule
import numpy as np
import pytest

from kinbot import geometry
from kinbot.conformers import Conformers
from kinbot.parameters import Parameters
from kinbot.reaction_generator import ReactionGenerator
from kinbot.reaction_path import endpoint_snapshot
from kinbot.stationary_pt import StationaryPoint


@pytest.mark.parametrize('failure', [ImportError('no RDKit'), ValueError('unsupported assignment')])
def test_unassigned_diastereotopic_channel_does_not_stop_other_reactions(tmp_path, monkeypatch, caplog, failure):
    from test_reaction_paths import peroxy, transfer
    from kinbot.species_routing import routing_name
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps({'barrier_threshold': 100.}))
    par = Parameters('input.json', show_warnings=False).par
    par.update(pes=0, high_level=0, L3_calc=0, delete_intermediate_files=1)
    parent = peroxy()
    reactions = []
    # H10 on the adjacent CH2 is diastereotopic; H7 on terminal CH3 is not.
    for hydrogen, donor in ((9, 1), (6, 0)):
        ts, product = transfer(parent, hydrogen, donor=donor)
        product.reset_order = Mock()  # keep this structural fixture's supplied endpoint graph
        reactions.append(SimpleNamespace(instance_name=ts.name, species=parent,
            instance=[5, hydrogen], products=[product], prod_opt=[], do_vdW=False,
            irc_prod=product, irc_product_reference=endpoint_snapshot(product)))
    parent.reac_obj, parent.reac_inst = reactions, [r.instance for r in reactions]
    parent.reac_name = [r.instance_name for r in reactions]
    parent.reac_type = ['intra_H_migration'] * 2
    parent.reac_ts_done, parent.reac_step = [3, 3], [0, 0]
    evidence = Path(reactions[0].instance_name + '.log')
    evidence.write_text('saved TS calculation\n')
    qc = Mock()
    qc.get_qc_energy.return_value = (0, -100.)
    qc.get_qc_zpe.return_value = (0, .1)
    qc.get_qc_freq.return_value = (0, [-1000.] + [500.] * (3 * parent.natom - 7))

    def optimize(point, *args):
        point.reduced_freqs = list(point.freq)
        return SimpleNamespace(species=point, shir=1, shigh=1, do_optimization=Mock())

    with patch('kinbot.stereo_identity._strings', side_effect=failure), \
            patch('kinbot.reaction_generator.Optimize', side_effect=optimize) as opt, \
            patch('kinbot.reaction_generator.time.sleep'), \
            patch('kinbot.reaction_generator.postprocess.createPESViewerInput'), \
            patch.object(ReactionGenerator, 'delete_files', side_effect=AssertionError('evidence deleted')):
        products = {routing_name(r.products[0]) + '_well': r.products[0] for r in reactions}
        qc.get_qc_geom.side_effect = lambda job, *args, **kwargs: (
            (0, products[job].geom.copy(), products[job].atom.copy())
            if kwargs.get('reorder') else (0, parent.geom.copy()))
        ReactionGenerator(parent, par, qc, 'input.json').generate()
    assert parent.reac_ts_done == [-999, -1]
    assert 'Cannot assign a demonstrated' in reactions[0].stereochemical_rejection
    assert not hasattr(reactions[0], 'ts_opt')
    assert reactions[1].ts.stereopath_id == 'ordinary'
    assert opt.call_count == 2  # only the ordinary TS and its product
    assert 'exported network is incomplete' in caplog.text
    assert 'using the legacy symmetry treatment' in caplog.text
    assert evidence.read_text() == 'saved TS calculation\n'
    qc.submit_qc.assert_not_called()


class TestTSConformerGuards(unittest.TestCase):
    def test_existing_and_reacting_bond_lengths_are_checked(self):
        point = SimpleNamespace(natom=3, wellorts=1,
            geom=np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]),
            bond=np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]))
        for distance, accepted in ((1.05, True), (1.2, False)):
            with self.subTest(distance=distance):
                candidate = copy.deepcopy(point)
                candidate.geom[2, 0] = 1. + distance
                self.assertEqual(bool(geometry.equal_geom(point, candidate, .1)), accepted)

    def test_new_nonreacting_bond_is_rejected_even_with_identical_lengths(self):
        point = SimpleNamespace(natom=3, wellorts=1,
            geom=np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.]]),
            bond=np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]))
        candidate = copy.deepcopy(point)
        candidate.bond[0, 2] = candidate.bond[2, 0] = 1
        self.assertFalse(geometry.equal_geom(point, candidate, .1))

    def test_conformer_polling_uses_the_ten_percent_geometry_guard(self):
        atoms = molecule('CH3OH')
        point = StationaryPoint('ts', 0, 1, atom=atoms.get_chemical_symbols(),
                                geom=atoms.positions, wellorts=1)
        point.characterize()
        search = Conformers.__new__(Conformers)
        search.species, search.semi_emp = point, 0
        search.get_job_name = lambda *args, **kwargs: 'conf/ts_0000'
        search.qc = Mock()
        search.qc.get_qc_geom.return_value = (0, point.geom)
        with patch('kinbot.conformers.geometry.equal_geom', return_value=False) as guard:
            self.assertEqual(search.test_conformer(0)[1], 1)
            self.assertEqual(guard.call_args.args[2], .1)


class TestDiscoveryOnlyIRC(unittest.TestCase):
    def test_completed_ts_conformers_keep_the_discovery_channel_without_new_ircs(self):
        with TemporaryDirectory() as directory:
            previous = Path.cwd()
            try:
                os.chdir(directory)
                Path('input.json').write_text(json.dumps({'barrier_threshold': 50.}))
                par = Parameters('input.json', show_warnings=False).par
                par.update(bimol=0, high_level=0, pes=0, L3_calc=0,
                           multi_conf_tst=1, rotor_scan=0)
                atoms = molecule('H2O')
                point = StationaryPoint('water', 0, 1,
                    atom=atoms.get_chemical_symbols(), geom=atoms.positions)
                point.characterize()
                discovery_product = copy.deepcopy(point)
                discovery_product.freq = discovery_product.reduced_freqs = [100., 200., 300.]
                # A one-fragment IRC product is the same object later optimized.
                product = discovery_product
                reaction = SimpleNamespace(instance_name='test_ts', species=point,
                    instance=[0, 1], products=[product], prod_opt=[], do_vdW=False,
                    irc_prod=discovery_product, irc=Mock(),
                    # This fixture enters stage 3, after stage 2 captured it.
                    irc_product_reference=endpoint_snapshot(discovery_product))
                point.reac_obj, point.reac_inst = [reaction], [[0, 1]]
                point.reac_name, point.reac_type = ['test_ts'], ['test_family']
                point.reac_ts_done, point.reac_step = [3], [0]
                qc = Mock()
                qc.get_qc_geom.side_effect = lambda name, *args, **kwargs: (
                    (0, atoms.positions[[1, 2, 0]], ['H', 'H', 'O'])
                    if kwargs.get('reorder') else (0, atoms.positions.copy()))
                qc.get_qc_energy.return_value = (0, -76.)
                qc.get_qc_zpe.return_value = (0, .01)
                qc.get_qc_freq.return_value = (0, [-1000., 200., 300.])

                def optimize(species, parameters, calculator):
                    # Retain two distinct members; neither launches an IRC.
                    if species.wellorts:
                        species.conformer_index = [0, 1]
                        species.conformer_geom = [species.geom.copy(), species.geom.copy() + .01]
                    species.reduced_freqs = list(species.freq)
                    return SimpleNamespace(species=species, shir=1, shigh=1,
                                           do_optimization=Mock())

                with patch('kinbot.reaction_generator.Optimize', side_effect=optimize) as opt, \
                        patch('kinbot.reaction_generator.IRC', side_effect=AssertionError('Unexpected IRC')), \
                        patch('kinbot.reaction_generator.postprocess.createPESViewerInput'), \
                        patch('kinbot.reaction_generator.time.sleep'):
                    ReactionGenerator(point, par, qc, 'input.json').generate()
                self.assertEqual(point.reac_ts_done, [-1])
                self.assertEqual(point.reac_obj, [reaction])
                self.assertIs(reaction.irc_prod, discovery_product)
                self.assertEqual(list(discovery_product.atom), ['H', 'H', 'O'])
                # Original TS coordinates are still O,H,H. The reordered
                # product must not invent an H-H bond in that coordinate order.
                np.testing.assert_array_equal(reaction.ts.bond, point.bond)
                np.testing.assert_array_equal(reaction.ts.reac_bond, np.zeros((3, 3)))
                self.assertEqual(reaction.ts.conformer_index, [0, 1])
                self.assertEqual(opt.call_count, 2)
                self.assertEqual(reaction.irc.mock_calls, [])
                self.assertFalse(Path('final_ts').exists())
            finally:
                os.chdir(previous)


if __name__ == '__main__':
    unittest.main()
