"""Final representative IRCs do not become an IRC per retained MC member."""
import copy
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from ase.build import molecule

from kinbot.parameters import Parameters
from kinbot.reaction_generator import ReactionGenerator
from kinbot.stationary_pt import StationaryPoint


class TestBimolecularIRCPolicy(unittest.TestCase):
    def test_ordinary_irc_keeps_existing_endpoint_identification(self):
        endpoint = SimpleNamespace(name='formaldehyde_product')
        irc = Mock()
        irc.check_irc.return_value = ['normal', 'error']
        irc.irc2stationary_pt.return_value = endpoint
        reaction = SimpleNamespace(instance_name='formaldehyde_h2_elimination')
        species = SimpleNamespace(reac_ts_done=[1])
        generator = ReactionGenerator(species, {'bimol': 0}, Mock(), 'input.json')
        with patch('kinbot.reaction_generator.IRC', return_value=irc):
            generator._advance_irc(0, reaction, 80.)
        irc.irc2stationary_pt.assert_called_once_with()
        self.assertIs(reaction.irc_prod, endpoint)
        self.assertEqual(species.reac_ts_done, [2])

    def test_product_stages_keep_selected_ts_and_members_without_more_ircs(self):
        with TemporaryDirectory() as directory:
            previous = Path.cwd()
            try:
                os.chdir(directory)
                Path('input.json').write_text(json.dumps({'barrier_threshold': 50.}))
                par = Parameters('input.json', show_warnings=False).par
                par.update(bimol=1, high_level=0, pes=0, L3_calc=0,
                           multi_conf_tst=1, rotor_scan=0)
                atoms = molecule('H2O')
                point = StationaryPoint('reactants', 0, 1,
                    atom=atoms.get_chemical_symbols(), geom=atoms.positions)
                point.characterize()
                product = copy.deepcopy(point)
                product.freq = product.reduced_freqs = [100., 200., 300.]
                selected = copy.deepcopy(point)
                selected.name, selected.wellorts = 'selected_ts', 1
                selected.freq = selected.reduced_freqs = [-1000., 200., 300.]
                selected.energy, selected.zpe = -76., .01
                selected.conformer_index = [0, 1]
                selected.conformer_geom = [selected.geom.copy(), selected.geom.copy() + .01]
                point.fragA = SimpleNamespace(energy=-70., zpe=.005)
                point.fragB = SimpleNamespace(energy=-6., zpe=.005)
                ts_opt = SimpleNamespace(species=selected, shir=1, shigh=1,
                                          do_optimization=Mock(), defer_hir=True,
                                          defer_cleanup=True,
                                          delete_files=Mock())
                reaction = SimpleNamespace(instance_name='test_ts', species=point,
                    instance=[0, 1], products=[product], prod_opt=[], do_vdW=False,
                    do_reac_vdW=False, irc_prod=copy.deepcopy(product),
                    ts=selected, ts_opt=ts_opt, irc=Mock(), final_irc_processed=True)
                endpoint = reaction.irc_prod
                point.reac_obj, point.reac_inst = [reaction], [[0, 1]]
                point.reac_name, point.reac_type = ['test_ts'], ['abstraction']
                point.reac_ts_done, point.reac_step = [3], [0]
                qc = Mock()
                qc.get_qc_geom.return_value = (0, atoms.positions.copy(), atoms.get_chemical_symbols())

                def optimize(species, parameters, calculator):
                    self.assertFalse(species.wellorts, 'Selected TS must not be rebuilt')
                    return SimpleNamespace(species=species, shir=1, shigh=1,
                        selected_job='product_well', do_optimization=Mock())

                with patch('kinbot.reaction_generator.Optimize', side_effect=optimize), \
                        patch('kinbot.reaction_generator.IRC', side_effect=AssertionError('Unexpected IRC')), \
                        patch('kinbot.reaction_generator.postprocess.createPESViewerInput'), \
                        patch('kinbot.reaction_generator.time.sleep'), \
                        patch.object(ReactionGenerator, '_classify_bimol_complexes'), \
                        patch.object(ReactionGenerator, '_finalize_bimol_complexes'), \
                        patch('kinbot.bimol.prepare_bimolecular_reactions'), \
                        patch('kinbot.bimol.release_final_ts_clusters'), \
                        patch('kinbot.bimol_network.BimolecularNetwork'), \
                        patch('kinbot.bimol.update_bimol_status'):
                    ReactionGenerator(point, par, qc, 'input.json').generate()
                self.assertEqual(point.reac_ts_done, [-1])
                self.assertEqual(point.reac_obj, [reaction])
                self.assertIs(reaction.ts, selected)
                self.assertIs(reaction.ts_opt.species, reaction.ts)
                self.assertIs(reaction.irc_prod, endpoint)
                self.assertEqual(reaction.ts.conformer_index, [0, 1])
                self.assertEqual(reaction.irc.mock_calls, [])
                self.assertFalse(Path('final_ts').exists())
                ts_opt.delete_files.assert_not_called()
                par.update(L3_calc=1, delete_intermediate_files=1)
                # Only terminal export runs on this second call. Sources must
                # still exist when final artifacts are written.
                def export():
                    ts_opt.delete_files.assert_not_called()
                    self.assertTrue(ts_opt.defer_hir)
                    self.assertTrue(ts_opt.defer_cleanup)

                with patch('kinbot.bimol_network.BimolecularNetwork') as network, \
                        patch('kinbot.bimol.update_bimol_status') as update_status, \
                        patch('kinbot.reaction_generator.postprocess.createPESViewerInput'):
                    network.return_value.write_unique_ts.side_effect = OSError('export failed')
                    with self.assertRaises(OSError):
                        ReactionGenerator(point, par, qc, 'input.json').generate()
                    ts_opt.delete_files.assert_not_called()
                    self.assertTrue(ts_opt.defer_hir)
                    self.assertIn('ready_for_export', [c.kwargs.get('final_status') for c in update_status.call_args_list])
                    self.assertNotIn('final_exported', [c.kwargs.get('final_status') for c in update_status.call_args_list])
                    update_status.reset_mock()
                    network.return_value.write_unique_ts.side_effect = export
                    ReactionGenerator(point, par, qc, 'input.json').generate()
                    self.assertEqual(update_status.call_args_list[-1].kwargs['final_status'], 'final_exported')
                ts_opt.delete_files.assert_called_once_with()
                self.assertFalse(ts_opt.defer_hir)
                self.assertFalse(ts_opt.defer_cleanup)
                ts_opt.do_optimization.assert_not_called()
            finally:
                os.chdir(previous)


if __name__ == '__main__':
    unittest.main()
