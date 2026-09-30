"""Configured names agree across jobs, restarts, object reuse and MESS/PES."""
import copy
import hashlib
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
from ase.db import connect
from kinbot.run_format import ensure_current_run
from kinbot import constants, symmetry
from kinbot.conformers import Conformers
from kinbot.calculation import load_calculation_record
from kinbot.mess import MESS, finalize_mc_mess
from kinbot.parameters import Parameters
from kinbot.pes import get_energy, write_input, create_mess_input
from kinbot import pes, postprocess
from kinbot.qc import QuantumChemistry
from kinbot.reaction_generator import ReactionGenerator
from kinbot.reaction_path import set_endpoint_populations
from kinbot.species_routing import routing_key, routing_name, input_species, is_species_name, expand_pes_names, configured_result_matches, mess_filename, same_species
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
        ensure_current_run(create=True)
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

    def test_reaction_exclusions_match_connectivity_and_exact_configured_names_before_qc(self):
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
        identity = canonical_identity(self.a)
        self.assertEqual(len(identity['id']), 64)
        self.assertEqual(routing_name(self.a), f"{self.a.chemid}-s{identity['id'][:16]}")

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
        self.record(routing_name(p)+'_well', p, reference=True)
        order = np.arange(p.natom)[::-1]
        q = StationaryPoint('reordered', p.charge, p.mult, atom=p.atom[order], geom=p.geom[order])
        q.characterize()
        self.assertEqual(routing_name(p), routing_name(q))
        with self.assertRaisesRegex(StereoRoutingError, 'atom indexing'):
            self.qc.qc_opt(q, q.geom)
        self.qc.submit_qc.assert_not_called()

    def test_vrc_current_sources_and_name_options_follow_the_same_namespace(self):
        p = self.a
        old, key = str(p.chemid), routing_name(p)
        Path('vrctst').mkdir()
        self.record(key+'_well', p, reference=True)
        source = f'vrctst/{key}_vts'
        self.record(source, p)
        Path(source+'.log').write_text('done\n')
        Path(source+'.chk').write_text('checkpoint\n')
        reaction = SimpleNamespace(instance_name=key+'_hom_sci_1_2', products=[p], species=p)
        p.reac_obj = [reaction]
        vts = VTS(p, self.par, self.qc)
        self.assertEqual(vts.configured_reactions([old+'_hom_sci_1_2']), [reaction.instance_name])
        vts.scan_reac[reaction.instance_name] = reaction
        with patch('kinbot.vrc_tst_scan.which', return_value='/bin/example'), \
                patch('kinbot.vrc_tst_scan.Popen') as process:
            process.return_value.communicate.return_value = (b'', b'')
            vts.save_products([reaction.instance_name])
            self.assertEqual(process.call_args.kwargs['args'], ['formchk', source+'.chk', source+'.fchk'])
        Path(source+'.cube').write_text('cube evidence\n')
        p.parent = '.'
        with patch('builtins.open', side_effect=RuntimeError('read cube')) as handle, \
                self.assertRaisesRegex(RuntimeError, 'read cube'):
            Fragment.pp_from_homo(p, 0)
        self.assertEqual(handle.call_args.args[0], './'+source+'.cube')
        self.qc.submit_qc.assert_not_called()

    def test_vrc_correction_json_restores_fragment_scope_and_connectivity_names(self):
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
        with self.assertRaisesRegex(ValueError, 'require the current fragment identities'):
            pes.create_rotdpy_inputs(par, [[root, name, products, 0.]], [])

    def test_termolecular_artifacts_do_not_exceed_filesystem_name_limits(self):
        products = [self.a, self.a, self.b]
        key = '_'.join(sorted(routing_name(p) for p in products))
        self.assertLess(len(key+'_0000.mess'), 255)
        writer = MESS(self.par, self.a)
        writer.termolec_names[key] = 't1'
        block = writer.write_termol(products, None, 0)
        filename = mess_filename(key, 0)
        self.assertEqual(filename, key+'_0000.mess')
        self.assertLess(len(filename), 255)
        self.assertEqual(Path(filename).read_text(), block)
        self.assertIn(key, block)
        self.assertTrue(mess_filename('very-long-network-' * 20, 0).startswith('network-'))

    def test_pes_inputs_and_energies_round_trip_two_diastereomers(self):
        for i, p in enumerate((self.a, self.b)):
            key = routing_name(p)
            write_input('input.json', p, 100., None, '.', 2)
            restored = input_species(f'{key}/{key}.json')
            self.assertEqual(routing_key(restored), routing_key(p))
            db = connect(f'{key}/kinbot.db')
            guard_well_job(SimpleNamespace(db=db), p, p.geom, key+'_well')
            db.write(Atoms(p.atom, positions=p.geom), name=key+'_well',
                data={'energy': (-100.+i) / constants.EVtoHARTREE, 'zpe': .01, 'status': 'normal'})
            self.assertAlmostEqual(get_energy([key], key, 0, 0)[0], -100.+i)
        key = routing_name(self.a)
        db = connect(f'{key}/kinbot.db')
        db.write(Atoms(self.b.atom, positions=self.b.geom), name=key+'_well',
            data={'energy': -200. / constants.EVtoHARTREE, 'zpe': .01, 'status': 'normal'})
        with self.assertRaises(ValueError):
            get_energy([key], key, 0, 0)

    def test_short_name_collision_cannot_reuse_a_different_full_identity(self):
        # Force the prefix collision; the remaining hash and actual chemical
        # graphs still distinguish the two real diastereomers.
        sha256 = hashlib.sha256
        selected_ids = {canonical_identity(p)['id'] for p in (self.a, self.b)}
        def colliding_hash(value):
            digest = sha256(value).hexdigest()
            return SimpleNamespace(hexdigest=lambda: 'a'*16 + digest[16:]
                                   if digest in selected_ids else digest)
        with patch('kinbot.stereo_identity.hashlib.sha256', side_effect=colliding_hash):
            first, second = copy.copy(self.a), copy.copy(self.b)
            for p in (first, second):
                p.__dict__.pop('optical_reference', None)
                optical_scope(p, 'specified')
            key = routing_name(first)
            self.assertEqual(key, routing_name(second))
            self.assertNotEqual(first.optical_reference['id'], second.optical_reference['id'])
            self.assertFalse(same_species(first, second))
            row_id = self.record(key+'_well', first, reference=True)
            row = self.qc.db.get(id=row_id)
            self.assertTrue(configured_result_matches(self.qc.db, key, row))
            self.assertEqual(len(self.qc.db.get(name='stereochemistry/'+key+'_well').data.identity['id']), 64)
            wrong_id = self.record(key+'_well', second)
            self.assertFalse(configured_result_matches(self.qc.db, key, self.qc.db.get(id=wrong_id)))
            with self.assertRaisesRegex(StereoRoutingError, 'different full stereoisomer identity'):
                guard_well_job(self.qc, second, second.geom, key+'_well')
            write_input('input.json', first, 100., None, '.', 2)
            filename = Path(key, key+'.json')
            original = filename.read_bytes()
            with self.assertRaisesRegex(StereoRoutingError, 'different full stereoisomer identity'):
                write_input('input.json', second, 100., None, '.', 2)
            self.assertEqual(filename.read_bytes(), original)

            for worker, p, energy in [('111', first, -100.), ('222', second, -90.)]:
                ensure_current_run(worker, create=True)
                db = connect(f'{worker}/kinbot.db')
                guard_well_job(SimpleNamespace(db=db), p, p.geom, key+'_well')
                db.write(Atoms(p.atom, positions=p.geom), name=key+'_well', data={
                    'status': 'normal', 'energy': energy / constants.EVtoHARTREE,
                    'zpe': .01, 'charge': p.charge, 'multiplicity': p.mult})
            for workers in (['111', '222'], ['222', '111']):
                with self.assertRaisesRegex(ValueError, 'different complete stereoisomer identities'):
                    get_energy(workers, key, 0, 0)

            # With no full input reference, a matching filename is insufficient.
            self.qc.db.delete([self.qc.db.get(name='stereochemistry/'+key+'_well').id])
            self.assertFalse(configured_result_matches(self.qc.db, key, row))

        # Reflected models need separate names even without a second QC job.
        identity = canonical_identity(self.a)
        selected_ids = {identity['id'], identity['mirror_id']}
        mirrored = copy.copy(self.a)
        mirrored.__dict__.pop('optical_reference', None)
        with patch('kinbot.stereo_identity.hashlib.sha256', side_effect=colliding_hash):
            with self.assertRaisesRegex(ValueError, 'same 16-character name suffix'):
                routing_name(mirrored)

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

    def test_final_direct_and_pes_keep_diastereomers_and_select_lowest_per_endpoint(self):
        # Synthetic channel observations isolate namespace and selection behavior;
        # these are not a physical reaction-rate benchmark.
        well = point('CC')
        base = routing_name(well)
        ensure_current_run(base, create=True)
        Path(base, 'me').mkdir()
        reactions = []
        for label, product, barrier in (('low', self.a, 30.), ('other', self.b, 40.),
                                         ('high', self.a, 50.)):
            ts = copy.copy(well)
            ts.name, ts.wellorts = base + '_test_' + label, 1
            ts.stereopath_id = 'ordinary'
            ts.stereopath_metadata = {'schema': 'kinbot.stereopath.v1', 'id': 'ordinary'}
            set_endpoint_populations(ts, [well], [product])
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
            if not getattr(p, 'wellorts', 0):
                guard_well_job(SimpleNamespace(db=db), p, p.geom, job)
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
