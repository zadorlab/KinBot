"""Configured names agree across jobs, restarts, object reuse and MESS/PES."""
import copy
import json
import logging
import io
import sys
from contextlib import redirect_stdout
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
import numpy as np
from ase import Atoms
from kinbot import constants, symmetry
from kinbot.conformers import Conformers
from kinbot.calculation import load_calculation_record
from kinbot.mess import MESS, finalize_mc_mess
from kinbot.parameters import Parameters
from kinbot.pes import get_energy, write_input, create_mess_input
from kinbot import pes, postprocess
from kinbot.qc import QuantumChemistry
from kinbot.reaction_generator import ReactionGenerator
from kinbot.species_routing import routing_key, routing_name, connect, resolve_job, prepare_qc_routing, input_species, is_species_name, expand_pes_names, configured_result_matches, mess_filename
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity, optical_scope
from kinbot.stereo_routing import StereoRoutingError, guard_well_job
from kinbot.vrc_tst_scan import VTS, fragment_routing_state
from kinbot.fragments import Fragment


def point(smiles):
    p = StationaryPoint('fixture', 0, 1, smiles=smiles)
    p.characterize()
    p.name = routing_name(p)
    p.energy, p.zpe = -100., .01
    p.freq = p.reduced_freqs = [100.] * (3*p.natom-6)
    symmetry.calculate_symmetry(p)
    return p


def input_data(p):
    return dict(charge=p.charge, mult=p.mult,
        structure=[v for a, xyz in zip(p.atom, p.geom) for v in [str(a), *map(float, xyz)]])


class TestSpeciesRouting(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        for directory in ('conf', 'hir', 'me'):
            Path(directory).mkdir()
        Path('input.json').write_text(json.dumps(dict(barrier_threshold=100.,
            high_level=0, rotor_scan=0, multi_conf_tst=0, conformer_search=0,
            me=0, uq=0, epsilon=100., sigma=3., queuing='local')))
        self.par = Parameters('input.json', show_warnings=False).par
        self.a = point('C[C@H](O)[C@H](F)C')
        self.b = point('C[C@H](O)[C@@H](F)C')
        self.qc = QuantumChemistry(self.par)
        self.qc.submit_qc = Mock()
        self.addCleanup(patch.stopall)
        patch('kinbot.pes.logger', logging.getLogger('KinBot'), create=True).start()

    def test_reaction_exclusions_match_legacy_and_exact_configured_names_before_qc(self):
        from kinbot.reaction_finder import ReactionFinder
        from test_reaction_paths import peroxy
        seed = peroxy()
        suffix = '_intra_H_migration_6_10'
        legacy = str(seed.chemid) + suffix
        exact = routing_name(seed) + suffix
        for reflected in (False, True):
            for selector in (legacy, exact):
                with self.subTest(reflected=reflected, selector=selector):
                    parent = copy.deepcopy(seed)
                    if reflected:
                        parent.geom *= [-1., 1., 1.]
                        parent.__dict__.pop('optical_reference', None)
                    par = dict(self.par, families=['intra_H_migration'], ringrange=[5, 6],
                               skip_reactions=[selector], delete_intermediate_files=0)
                    ReactionFinder(parent, par, None).find_reactions()
                    reaction = next(r for r in parent.reac_obj if r.instance_name.endswith(suffix))
                    parent.reac_obj, parent.reac_inst = [reaction], [reaction.instance]
                    parent.reac_name, parent.reac_type = [reaction.instance_name], ['intra_H_migration']
                    parent.reac_ts_done, parent.reac_step = [0], [0]
                    qc = Mock()
                    qc.check_qc.side_effect = RuntimeError('unexcluded reaction reached QC')
                    with patch('kinbot.reaction_generator.time.sleep'):
                        generator = ReactionGenerator(parent, par, qc, 'input.json')
                        if reflected and selector == exact:
                            with self.assertRaisesRegex(RuntimeError, 'unexcluded reaction'):
                                generator.generate()
                            qc.check_qc.assert_called_once_with(reaction.instance_name)
                        else:
                            generator.generate()
                            self.assertEqual(parent.reac_ts_done, [-999])
                            self.assertEqual(qc.mock_calls, [])

    def record(self, name, p, energy=-100., reference=False):
        if reference:
            guard_well_job(self.qc, p, p.geom, name)
        return self.qc.db.write(Atoms(p.atom, positions=p.geom), name=name,
            data={'status': 'normal', 'energy': energy / constants.EVtoHARTREE,
                  'charge': p.charge, 'multiplicity': p.mult,
                  'zpe': .01, 'frequencies': p.freq, 'hess': np.eye(3*p.natom)})

    def test_ordinary_names_and_connectivity_ids_are_unchanged(self):
        for smi in ('CO', 'CC', '[CH3]O', 'O'):
            p = point(smi)
            self.assertEqual(routing_key(p), p.chemid)
        self.assertEqual(self.a.chemid, self.b.chemid)
        self.assertNotEqual(routing_key(self.a), routing_key(self.b))
        self.assertNotIn('_', routing_name(self.a))

    def test_atom_order_and_rigid_motion_keep_the_key_but_mirrors_do_not(self):
        p = self.a
        order = np.arange(p.natom)[::-1]
        q = StationaryPoint('permuted', p.charge, p.mult, atom=p.atom[order],
            geom=(p.geom @ np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]) + 3)[order])
        q.characterize()
        self.assertEqual(routing_key(p), routing_key(q))
        q = copy.copy(p)
        q.geom = p.geom * [-1, 1, 1]
        self.assertNotEqual(routing_key(p), routing_key(q))
        optical_scope(p, 'racemic')
        optical_scope(q, 'racemic')
        self.assertNotEqual(routing_key(p), routing_key(q))
        self.assertTrue(optical_scope(p, 'racemic')['mirror_allowed'])

    def test_achiral_configurations_with_stereo_tags_also_get_separate_keys(self):
        trans, cis = point('F/C=C/F'), point('F/C=C\\F')
        self.assertFalse(canonical_identity(trans)['is_chiral_configuration'])
        self.assertEqual(trans.chemid, cis.chemid)
        self.assertNotEqual(routing_key(trans), routing_key(cis))
        self.assertTrue(is_species_name(routing_name(cis)))

    def test_separate_names_cannot_double_count_the_same_full_racemate(self):
        mirror = point('C[C@@H](O)[C@@H](F)C')
        par = dict(self.par, multi_conf_tst=1, optical_population='racemic')
        writer = MESS(par, self.a)
        writer.well_names = {routing_key(self.a): 'w1', routing_key(mirror): 'w2'}
        blocks = [writer.write_well(p, 0., 1., 0) for p in (self.a, mirror)]
        self.assertIn('! kinbot_racemic_population', blocks[0])
        self.assertEqual(finalize_mc_mess(blocks[0]), blocks[0])
        with self.assertRaisesRegex(ValueError, 'Overlapping racemic populations'):
            finalize_mc_mess(''.join(blocks))
        # Repeated appearances of the same fragment population are legitimate.
        finalize_mc_mess(blocks[0] + blocks[0])
        writer.par['optical_population'] = 'specified'
        blocks = [writer.write_well(p, 0., 1., 0) for p in (self.a, mirror)]
        finalize_mc_mess(''.join(blocks))

    def test_object_reuse_merges_only_the_same_configuration(self):
        repeat = copy.copy(self.a)
        fragments = [self.a, self.b, repeat]
        ReactionGenerator.equate_identical(None, fragments)
        self.assertIs(fragments[2], self.a)
        self.assertIs(fragments[1], self.b)
        unique = [self.a]
        candidates = [self.b, repeat]
        ReactionGenerator.equate_unique(None, candidates, unique)
        self.assertEqual(len(unique), 2)
        self.assertIs(candidates[1], self.a)

    def test_backend_writers_and_conformer_readers_use_the_same_namespace(self):
        for backend in ('gauss', 'qchem', 'nwchem', 'fc', 'orca'):
            par = dict(self.par, qc=backend, use_sella=int(backend != 'gauss'))
            qc = QuantumChemistry(par)
            qc.submit_qc = Mock()
            for p in (self.a, self.b):
                base = routing_name(p)
                qc.qc_opt(p, p.geom)
                self.assertEqual(qc.submit_qc.call_args.args[0], base+'_well')
                qc.qc_conf(p, p.geom, 0)
                confs = Conformers(p, par, qc)
                self.assertEqual(qc.submit_qc.call_args.args[0], confs.get_job_name(0))
                qc.qc_ring_conf(p, p.geom, [], [], 0, 0)
                self.assertEqual(qc.submit_qc.call_args.args[0], f'conf/{base}_r0000_0000')
                qc.qc_hir(p, p.geom, 0, 0, [[1, 2, 4, 5]], 0)
                self.assertEqual(qc.submit_qc.call_args.args[0], f'hir/{base}_hir_0_00')

    def test_completed_legacy_well_can_be_read_without_rewriting_its_files(self):
        p = self.a
        legacy, key = str(p.chemid), routing_name(p)
        self.record(legacy+'_well', p)
        Path(legacy+'_well.py').write_text('# original script\n')
        Path(legacy+'_well.log').write_text('done\n')
        self.qc.qc_opt(p, p.geom)
        self.qc.submit_qc.assert_not_called()
        self.assertEqual(Path(legacy+'_well.py').read_text(), '# original script\n')
        self.assertEqual(resolve_job(self.qc.db, key+'_well'), legacy+'_well')
        self.assertAlmostEqual(self.qc.get_qc_energy(key+'_well')[1], -100.)
        np.testing.assert_array_equal(self.qc.get_qc_geom(key+'_well', p.natom)[1], p.geom)
        self.assertEqual(self.qc.get_qc_freq(key+'_well', p.natom)[1], p.freq)
        self.assertEqual(self.qc.get_qc_zpe(key+'_well')[1], .01)
        load_calculation_record(p, self.qc, key+'_well')
        self.assertEqual(p.source_job, legacy+'_well')
        self.assertEqual(p.requested_source_job, key+'_well')
        self.record(f'conf/{legacy}_0000', p)
        self.assertEqual(resolve_job(self.qc.db, f'conf/{key}_0000'), f'conf/{key}_0000')
        # A new process reads the same recorded alias.
        self.assertEqual(resolve_job(connect('kinbot.db'), key+'_well'), legacy+'_well')

    def test_original_input_licenses_children_and_new_results_supersede_aliases(self):
        p = self.a
        legacy, key = str(p.chemid), routing_name(p)
        self.record(legacy+'_well', p, reference=True)
        for suffix in ('_0000', '_low'):
            self.record('conf/'+legacy+suffix, p)
        prepare_qc_routing(self.qc, p)
        self.assertEqual(resolve_job(self.qc.db, f'conf/{key}_0000'), f'conf/{legacy}_0000')
        self.assertEqual(list(self.qc.db.select(name=f'conf/{key}_low'))[-1].name,
                         f'conf/{legacy}_low')
        self.record(key+'_well', p, -99.)
        self.assertEqual(resolve_job(self.qc.db, key+'_well'), key+'_well')
        self.assertEqual(len(list(self.qc.db.raw.select(name=legacy+'_well'))), 1)

    def test_ambiguous_pending_restart_requires_explicit_saved_input(self):
        p = self.a
        legacy = str(p.chemid)
        Path(legacy+'_well.py').write_text('# unidentified pending input\n')
        prepare_qc_routing(self.qc, p)
        configured = routing_name(p) + '_well'
        self.assertEqual(resolve_job(self.qc.db, configured), configured)
        self.assertEqual(Path(legacy+'_well.py').read_text(), '# unidentified pending input\n')
        Path('original.json').write_text(json.dumps(input_data(p)))
        self.qc.par['stereo_legacy_inputs'] = {legacy: 'original.json'}
        prepare_qc_routing(self.qc, p)
        self.assertEqual(resolve_job(self.qc.db, routing_name(p)+'_well'), legacy+'_well')

    def test_other_legacy_isomer_does_not_block_a_new_namespace(self):
        self.record(str(self.a.chemid)+'_well', self.a, reference=True)
        prepare_qc_routing(self.qc, self.b)
        key = routing_name(self.b)
        self.assertEqual(resolve_job(self.qc.db, key+'_well'), key+'_well')
        self.qc.qc_opt(self.b, self.b.geom)
        self.assertEqual(self.qc.submit_qc.call_args.args[0], key+'_well')

    def test_later_original_input_can_verify_child_jobs_and_auxiliary_writers(self):
        p = self.a
        legacy, key = str(p.chemid), routing_name(p)
        self.record(legacy+'_well', p)
        self.record(f'conf/{legacy}_0000', p)
        prepare_qc_routing(self.qc, p)
        self.assertEqual(resolve_job(self.qc.db, f'conf/{key}_0000'), f'conf/{key}_0000')
        Path(f'{legacy}.json').write_text(json.dumps(input_data(p)))
        for folder in ('aie', 'vrctst'):
            Path(folder).mkdir()
        for job in (f'aie/{legacy}_AIE0_0', f'aie/{legacy}_AIE1_0', f'vrctst/{legacy}_vts'):
            Path(job+'.py').write_text('# pending legacy job\n')
            self.record(job, p)
            Path(job+'.log').write_text('done\n')
        prepare_qc_routing(self.qc, p)
        self.assertEqual(resolve_job(self.qc.db, f'conf/{key}_0000'), f'conf/{legacy}_0000')
        self.qc.qc_aie(p, p.geom, '0')
        self.qc.qc_vts_frag(p)
        self.qc.submit_qc.assert_not_called()
        self.assertFalse(Path(f'aie/{key}_AIE0_0.py').exists())

    def test_pes_connectivity_selectors_expand_and_keep_declared_scope(self):
        p = self.a
        optical_scope(p, 'racemic')
        p.geom = p.geom * [-1, 1, 1]
        self.par['optical_population'] = 'racemic'
        Path('input.json').write_text(json.dumps(self.par))
        for species in (p, self.b):
            write_input('input.json', species, 100., None, '.', 2)
        names = expand_pes_names('.', [str(p.chemid)])
        self.assertCountEqual(names, [routing_name(p), routing_name(self.b)])
        self.assertEqual(expand_pes_names('.', [routing_name(p)]), [routing_name(p)])
        pes.write_input_keep('input.json', routing_name(p), '.')
        restored = input_species(f'{routing_name(p)}/{routing_name(p)}.json')
        self.assertEqual(routing_name(restored), routing_name(p))
        self.assertEqual(restored.optical_population, 'racemic')

    def test_changed_pes_directory_alias_preserves_input_and_refuses_reuse(self):
        p = self.a
        legacy = Path(str(p.chemid))
        legacy.mkdir()
        filename = legacy / f'{p.chemid}.json'
        filename.write_text(json.dumps(input_data(p)))
        write_input('input.json', p, 100., None, '.', 2)
        filename.write_text(json.dumps(input_data(self.b)))
        with self.assertRaises(StereoRoutingError):
            write_input('input.json', p, 100., None, '.', 2)
        self.assertEqual(json.loads(filename.read_text()), input_data(self.b))

    def test_late_guard_reference_does_not_license_old_children(self):
        p = self.a
        legacy, key = str(p.chemid), routing_name(p)
        self.record(legacy+'_well', p)
        guard_well_job(self.qc, p, p.geom, legacy+'_well')
        self.record(f'conf/{legacy}_0000', p)
        prepare_qc_routing(self.qc, p)
        self.assertEqual(resolve_job(self.qc.db, key+'_well'), legacy+'_well')
        self.assertEqual(resolve_job(self.qc.db, f'conf/{key}_0000'), f'conf/{key}_0000')

    def test_changed_legacy_result_cannot_hide_behind_a_recorded_alias(self):
        p = self.a
        legacy = str(p.chemid)
        self.record(legacy+'_well', p, reference=True)
        prepare_qc_routing(self.qc, p)
        self.record(legacy+'_well', self.b)
        with self.assertRaises(StereoRoutingError):
            prepare_qc_routing(self.qc, p)

    def test_invalidated_alias_starts_a_new_canonical_job(self):
        p = self.a
        legacy, key = str(p.chemid), routing_name(p)
        self.record(legacy+'_well', p, reference=True)
        self.record(legacy+'_well_high', p)
        Path(legacy+'_well_high.log').write_text('done\n')
        prepare_qc_routing(self.qc, p)
        self.qc.invalidate_qc(key+'_well_high')
        self.qc.qc_opt(p, p.geom, high_level=1)
        self.assertEqual(self.qc.submit_qc.call_args.args[0], key+'_well_high')
        self.assertTrue(Path(key+'_well_high.py').exists())
        self.assertEqual(list(self.qc.db.select(name=legacy+'_well_high'))[-1].data.status, 0)
        self.assertTrue(list(Path('.').glob(legacy+'_well_high.log.restart_*')))

    def test_legacy_conformer_summary_copies_the_actual_completed_source(self):
        p = self.a
        legacy, key = str(p.chemid), routing_name(p)
        self.record(legacy+'_well', p, reference=True)
        self.record(f'conf/{legacy}_0000', p)
        Path(legacy+'_well.log').write_text('parent done\n')
        Path(f'conf/{legacy}_0000.log').write_text('conformer done\n')
        prepare_qc_routing(self.qc, p)
        for status, source in ((0, 'conformer'), (1, 'parent')):
            conf = Conformers(p, self.par, self.qc)
            conf.conf, conf.conf_status = 1, [status]
            self.assertEqual(conf.check_conformers()[0], 1)
            self.assertEqual(Path(f'conf/{key}_low.log').read_text(), source+' done\n')

    def test_result_without_electronic_state_cannot_license_a_requested_charge(self):
        p = self.a
        old = str(p.chemid)+'_well'
        self.qc.db.write(Atoms(p.atom, positions=p.geom), name=old,
            data={'status': 'normal', 'energy': -100., 'zpe': .01})
        requested = copy.copy(p)
        requested.charge = 2
        prepare_qc_routing(self.qc, requested)
        configured = routing_name(requested) + '_well'
        self.assertEqual(resolve_job(self.qc.db, configured), configured)
        self.assertEqual(len(list(self.qc.db.select(name=old))), 1)

    def test_late_reference_cannot_supply_unreported_original_electronic_state(self):
        p = self.a
        old = str(p.chemid)+'_well'
        self.qc.db.write(Atoms(p.atom, positions=p.geom), name=old,
            data={'status': 'normal', 'energy': -100., 'zpe': .01})
        requested = copy.copy(p)
        requested.charge = 2
        guard_well_job(self.qc, requested, requested.geom, old)
        prepare_qc_routing(self.qc, requested)
        configured = routing_name(requested) + '_well'
        self.assertEqual(resolve_job(self.qc.db, configured), configured)

    def test_ordinary_legacy_smiles_keeps_working_without_optional_rdkit(self):
        p = point('CO')
        old = str(p.chemid)
        Path(old).mkdir()
        Path(old, old+'.json').write_text(json.dumps({'smiles': 'CO', 'charge': 0, 'mult': 1}))
        with patch.dict(sys.modules, {'rdkit': None}):
            self.assertEqual(expand_pes_names('.', [old]), [old])

    def test_charged_exact_alias_keeps_its_trusted_reference_for_pes(self):
        p = self.a
        p.charge, p.mult = 1, 2
        p.calc_chemid()
        old = str(p.chemid)+'_well'
        self.record(old, p)
        prepare_qc_routing(self.qc, p)
        row = list(self.qc.db.select(name=old))[-1]
        self.assertTrue(configured_result_matches(self.qc.db, routing_name(p), row))

    def test_unspecified_saved_smiles_cannot_license_a_legacy_prefix(self):
        p = self.a
        old = str(p.chemid)
        Path(old+'_well.py').write_text('# old pending job\n')
        Path(old+'.json').write_text(json.dumps({'smiles': 'CC(O)C(F)C', 'mult': 1, 'charge': 0}))
        prepare_qc_routing(self.qc, p)
        configured = routing_name(p) + '_well'
        self.assertEqual(resolve_job(self.qc.db, configured), configured)
        self.assertEqual(Path(old+'_well.py').read_text(), '# old pending job\n')
        self.assertFalse(list(self.qc.db.select(name='stereo_route/'+routing_name(p))))

    def test_mirror_selected_racemic_legacy_input_licenses_its_declared_namespace(self):
        p = self.a
        optical_scope(p, 'racemic')
        p.geom = p.geom * [-1, 1, 1]
        data = dict(input_data(p), optical_population='racemic', stereo_reference=p.optical_reference)
        Path(str(p.chemid)+'.json').write_text(json.dumps(data))
        self.qc.par['optical_population'] = 'racemic'
        self.record(str(p.chemid)+'_well', p)
        prepare_qc_routing(self.qc, p)
        self.assertEqual(resolve_job(self.qc.db, routing_name(p)+'_well'), str(p.chemid)+'_well')

    def test_actual_pes_entrypoint_restores_reference_before_writing(self):
        p = self.a
        optical_scope(p, 'racemic')
        p.geom = p.geom * [-1, 1, 1]
        data = dict(self.par, **input_data(p), optical_population='racemic', stereo_reference=p.optical_reference)
        Path('input.json').write_text(json.dumps(data))
        def capture(_input, species, *args):
            self.assertEqual(routing_name(species), routing_name(p))
            self.assertEqual(species.optical_population, 'racemic')
            raise RuntimeError('stop before job submission')
        with patch.object(sys, 'argv', ['pes', 'input.json']), \
                patch.object(pes, 'write_input', side_effect=capture), \
                patch.object(pes, 'config_log', return_value=logging.getLogger('KinBot')), \
                redirect_stdout(io.StringIO()), \
                self.assertRaisesRegex(RuntimeError, 'stop before job submission'):
            pes.main()

    def test_configured_cache_requires_compatible_atom_indexing(self):
        p = self.a
        self.record(str(p.chemid)+'_well', p, reference=True)
        order = np.arange(p.natom)[::-1]
        q = StationaryPoint('reordered', p.charge, p.mult, atom=p.atom[order], geom=p.geom[order])
        q.characterize()
        self.assertEqual(routing_name(p), routing_name(q))
        self.qc.qc_opt(q, q.geom)
        configured = routing_name(q) + '_well'
        self.assertEqual(resolve_job(self.qc.db, configured), configured)
        self.assertEqual(self.qc.submit_qc.call_args.args[0], configured)

    def test_legacy_no_kinbot_requires_summary_regeneration_without_deletion(self):
        p = self.a
        old, key = str(p.chemid), routing_name(p)
        Path(old).mkdir()
        Path(old, old+'.json').write_text(json.dumps(input_data(p)))
        original = Path(old, f'summary_{old}.out')
        original.write_text('legacy results\n')
        write_input('input.json', p, 100., None, '.', 2)
        with self.assertRaisesRegex(ValueError, 'legacy PES summary needs'):
            pes.get_wells(key)
        with self.assertRaisesRegex(ValueError, 'legacy PES summary needs'):
            pes.postprocess(self.par, [key], 'all', [], p.mass)
        with patch('kinbot.pes.subprocess.Popen', return_value=SimpleNamespace(pid=1)), \
                patch('kinbot.pes.time.sleep'):
            pes.submit_job(key, self.par)
        self.assertEqual(original.read_text(), 'legacy results\n')

    def test_vrc_verified_sources_and_selectors_follow_the_same_namespace(self):
        p = self.a
        old, key = str(p.chemid), routing_name(p)
        Path('vrctst').mkdir()
        self.record(old+'_well', p, reference=True)
        legacy = f'vrctst/{old}_vts'
        self.record(legacy, p)
        Path(legacy+'.log').write_text('done\n')
        Path(legacy+'.chk').write_text('checkpoint\n')
        prepare_qc_routing(self.qc, p)
        reaction = SimpleNamespace(instance_name=key+'_hom_sci_1_2', products=[p], species=p)
        p.reac_obj = [reaction]
        vts = VTS(p, self.par, self.qc)
        self.assertEqual(vts.configured_reactions([old+'_hom_sci_1_2']), [reaction.instance_name])
        vts.scan_reac[reaction.instance_name] = reaction
        with patch('kinbot.vrc_tst_scan.which', return_value='/bin/example'), \
                patch('kinbot.vrc_tst_scan.Popen') as process:
            process.return_value.communicate.return_value = (b'', b'')
            vts.save_products([reaction.instance_name])
            self.assertEqual(process.call_args.kwargs['args'], ['formchk', legacy+'.chk', legacy+'.fchk'])
        Path(legacy+'.cube').write_text('cube evidence\n')
        p.parent = '.'
        with patch('builtins.open', side_effect=RuntimeError('read cube')) as handle, \
                self.assertRaisesRegex(RuntimeError, 'read cube'):
            Fragment.pp_from_homo(p, 0)
        self.assertEqual(handle.call_args.args[0], './'+legacy+'.cube')
        scan = f'vrctst/{old}_hom_sci_1_2_vts_pt00'
        self.record(scan, p)
        Path(scan+'.log').write_text('done\n')
        job = self.qc.qc_vts(reaction, p.geom, 0, [], False, p.geom)
        self.assertEqual(job, f'vrctst/{key}_hom_sci_1_2_vts_pt00')
        self.assertEqual(self.qc.get_qc_geom(job, p.natom)[0], 0)
        asymptote = f'vrctst/{old}_hom_sci_1_2_vts_pt_asymptote'
        self.record(asymptote, p)
        Path(asymptote+'.log').write_text('done\n')
        reaction.do_vdW, reaction.scan_coo, reaction.maps, reaction.equiv = False, [0, 1], [[], [1]], []
        with patch('kinbot.vrc_tst_scan.time.sleep'):
            vts.do_scan([reaction.instance_name], noscan=True)
        self.qc.submit_qc.assert_not_called()

    def test_vrc_correction_json_restores_fragment_scope_and_legacy_selectors(self):
        parent = self.b
        root = routing_name(parent)
        Path(root, 'vrctst').mkdir(parents=True)
        p = self.a
        p.charge, p.mult = 1, 2
        p.calc_chemid()
        optical_scope(p, 'racemic')
        p.optical_population = 'racemic'
        p.geom = p.geom * [-1, 1, 1]
        other = point('[H]')
        products = [routing_name(p), routing_name(other)]
        name = root+'_hom_sci_1_2'
        data = {'dist': [1.], 'e_inf_samp': -100., 'ra': [[0], [0]],
                'unique': [[[0]], [[0]]], 'frags_atom': [list(p.atom), list(other.atom)],
                'frags_geom': [p.geom.tolist(), other.geom.tolist()],
                'frags_mult': [p.mult, other.mult],
                'frags_routing': [fragment_routing_state(p), fragment_routing_state(other)]}
        filename = Path(root, 'vrctst', 'corr_'+name+'.json')
        filename.write_text(json.dumps(data))
        par = dict(self.par, rotdpy_dist=[10.],
                   vrc_tst_scan={str(parent.chemid): [str(parent.chemid)+'_hom_sci_1_2']})
        def verify(**kwargs):
            restored = kwargs['fragments'][0]
            self.assertEqual(restored.charge, 1)
            self.assertEqual(restored.mult, 2)
            self.assertEqual(restored.optical_population, 'racemic')
            self.assertEqual(routing_name(restored), routing_name(p))
            raise RuntimeError('verified configured surface input')
        with patch.object(Fragment, '_instances', []), \
                patch('kinbot.pes.pp_settings.create_all_surf_for_dist', side_effect=verify), \
                self.assertRaisesRegex(RuntimeError, 'verified configured surface input'):
            pes.create_rotdpy_inputs(par, [[root, name, products, 0.]], [])
        data['frags_atom'][0] = list(parent.atom)
        data['frags_geom'][0] = parent.geom.tolist()
        data['frags_mult'][0] = parent.mult
        data['frags_routing'][0] = fragment_routing_state(parent)
        filename.write_text(json.dumps(data))
        with patch.object(Fragment, '_instances', []), \
                self.assertRaisesRegex(ValueError, 'disagree with the requested configured endpoints'):
            pes.create_rotdpy_inputs(par, [[root, name, products, 0.]], [])
        data.pop('frags_routing')
        filename.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'lacks original state and population'):
            pes.create_rotdpy_inputs(par, [[root, name, products, 0.]], [])

    def test_termolecular_artifacts_do_not_exceed_filesystem_name_limits(self):
        products = [self.a, self.a, self.b]
        key = '_'.join(sorted(routing_name(p) for p in products))
        self.assertGreater(len(key+'_0000.mess'), 255)
        writer = MESS(self.par, self.a)
        writer.termolec_names[key] = 't1'
        block = writer.write_termol(products, None, 0)
        filename = mess_filename(key, 0)
        self.assertLess(len(filename), 255)
        self.assertEqual(Path(filename).read_text(), block)
        self.assertIn(key, block)

    def test_pes_inputs_and_energies_round_trip_two_diastereomers(self):
        for i, p in enumerate((self.a, self.b)):
            key = routing_name(p)
            write_input('input.json', p, 100., None, '.', 2)
            restored = input_species(f'{key}/{key}.json')
            self.assertEqual(routing_key(restored), routing_key(p))
            db = connect(f'{key}/kinbot.db')
            db.write(Atoms(p.atom, positions=p.geom), name=key+'_well',
                data={'energy': (-100.+i) / constants.EVtoHARTREE, 'zpe': .01, 'status': 'normal'})
            self.assertAlmostEqual(get_energy([key], key, 0, 0)[0], -100.+i)
        key = routing_name(self.a)
        db = connect(f'{key}/kinbot.db')
        db.write(Atoms(self.b.atom, positions=self.b.geom), name=key+'_well',
            data={'energy': -200. / constants.EVtoHARTREE, 'zpe': .01, 'status': 'normal'})
        with self.assertRaises(ValueError):
            get_energy([key], key, 0, 0)

    def test_racemic_selected_mirror_preserves_declared_name_across_pes_restart(self):
        p = self.a
        optical_scope(p, 'racemic')
        p.optical_population = 'racemic'
        key = routing_name(p)
        p.geom = p.geom * [-1, 1, 1]
        self.par['optical_population'] = 'racemic'
        Path('input.json').write_text(json.dumps(self.par))
        write_input('input.json', p, 100., None, '.', 2)
        restored = input_species(f'{key}/{key}.json')
        self.assertEqual(routing_name(restored), key)
        self.assertEqual(canonical_identity(restored)['id'], p.optical_reference['mirror_id'])
        db = connect(f'{key}/kinbot.db')
        qc = SimpleNamespace(db=db, par=self.par)
        guard_well_job(qc, restored, restored.geom, key+'_well')
        db.write(Atoms(restored.atom, positions=restored.geom), name=key+'_well',
            data={'status': 'normal', 'energy': -100./constants.EVtoHARTREE, 'zpe': .01})
        self.assertAlmostEqual(get_energy([key], key, 0, 0)[0], -100.)
        data = json.loads(Path(f'{key}/{key}.json').read_text())
        data['optical_population'] = 'specified'
        Path('invalid.json').write_text(json.dumps(data))
        with self.assertRaises(StereoRoutingError):
            input_species('invalid.json')

    def test_verified_legacy_pes_directory_keeps_original_input_and_results(self):
        p = self.a
        old, key = str(p.chemid), routing_name(p)
        Path(old).mkdir()
        original = json.dumps(input_data(p))
        Path(f'{old}/{old}.json').write_text(original)
        Path(f'{old}/result.log').write_text('old result\n')
        write_input('input.json', p, 100., None, '.', 2)
        self.assertTrue(Path(key).is_symlink())
        self.assertEqual(Path(f'{old}/{old}.json').read_text(), original)
        self.assertEqual(Path(f'{key}/result.log').read_text(), 'old result\n')
        self.assertEqual(routing_key(input_species(f'{key}/{key}.json')), routing_key(p))

    def test_final_direct_and_pes_keep_diastereomers_and_select_lowest_per_endpoint(self):
        # Synthetic channel observations isolate namespace and selection behavior;
        # these are not a physical reaction-rate benchmark.
        well = point('CC')
        base = routing_name(well)
        Path(base).mkdir()
        Path(base, 'me').mkdir()
        reactions = []
        for label, product, barrier in (('low', self.a, 30.), ('other', self.b, 40.),
                                         ('high', self.a, 50.)):
            ts = copy.copy(well)
            ts.name, ts.wellorts = base + '_test_' + label, 1
            ts.energy = well.energy + barrier / constants.AUtoKCAL
            ts.freq = ts.reduced_freqs = [-1000.] + well.freq[1:]
            reactions.append(SimpleNamespace(instance_name=ts.name, ts=ts,
                products=[product], prod_opt=[SimpleNamespace(species=product)],
                do_vdW=False, mp2=0))
        well.reac_obj, well.reac_ts_done = reactions, [-1]*3
        well.reac_type, well.reac_inst = ['test']*3, [None]*3
        well.reac_name = [r.instance_name for r in reactions]
        par = dict(self.par, **input_data(well), smiles='', me=2)
        db = connect(f'{base}/kinbot.db')
        for p, job in [(well, base+'_well'),
                       (self.a, routing_name(self.a)+'_well'),
                       (self.b, routing_name(self.b)+'_well')]+[(r.ts, r.ts.name) for r in reactions]:
            db.write(Atoms(p.atom, positions=p.geom), name=job,
                data={'status': 'normal', 'energy': p.energy / constants.EVtoHARTREE,
                      'zpe': p.zpe, 'frequencies': p.freq})
        previous = Path.cwd()
        try:
            os.chdir(base)
            writer = MESS(dict(par, pes=0), well)
            writer.write_input(None)
            direct = Path('me/mess_0000.inp').read_text()
            self.assertEqual(len(writer.well_names), 3)
            writer = MESS(dict(par, pes=1), well)
            writer.write_input(None)
            postprocess.create_summary_file(well, SimpleNamespace(qc='fc'), par)
        finally:
            os.chdir(previous)
        with patch('kinbot.pes.copy_from_kinbot'), \
             patch('kinbot.pes.create_pesviewer_input'), patch('kinbot.pes.create_rotdpy_inputs'), \
             patch('kinbot.pes.t1_analysis'), patch('kinbot.pes.get_l3energy', return_value=(0, -1)):
            pes.postprocess(par, [base], 'all', [], well.mass)
        deferred = Path('me/mess_0000.inp').read_text()
        for output in (direct, deferred):
            self.assertIn(routing_name(self.a), output)
            self.assertIn(routing_name(self.b), output)
            # Achiral entrance also reaches each product's global mirror.
            self.assertEqual(sum(line.strip().startswith('Barrier ') for line in output.splitlines()), 4)
            self.assertEqual(output.count('derived by global reflection from'), 4)
            self.assertNotIn(base+'_test_high', output)
            self.assertIn(base+'_test_low', output)


if __name__ == '__main__':
    unittest.main()
