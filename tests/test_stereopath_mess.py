"""Final direct/PES files sum stereochemical routes, not duplicate observations."""
import copy
import json
import logging
import math
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot import constants, pes, postprocess, symmetry
from kinbot.mess import MESS, union_stereochemical_barriers
from kinbot.parameters import Parameters
from kinbot.species_routing import routing_name
from kinbot.stereo_routing import StereoRoutingError
from tests.counting_fixtures import methanol_data, saved_point
from test_reaction_paths import peroxy, transfer
from test_mess_conformers import values


class TestStereopathMESS(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        Path('input.json').write_text(json.dumps(dict(barrier_threshold=100.,
            high_level=0, rotor_scan=0, multi_conf_tst=0, conformer_search=0,
            me=0, uq=0, epsilon=100., sigma=3., queuing='local')))
        self.par = Parameters('input.json', show_warnings=False).par
        self.addCleanup(patch.stopall)
        patch('kinbot.pes.logger', logging.getLogger('KinBot'), create=True).start()

    def observations(self):
        p = peroxy()
        p.name = routing_name(p)
        routes = []
        for hydrogen, barrier in ((9, 30.), (10, 32.), (9, 35.)):
            ts, q = transfer(p, hydrogen)
            # Give the endpoint a geometry with its stated O-H connectivity.
            # These remain structural fixtures, not calculated reaction energies.
            direction = q.geom[5]-q.geom[4]
            q.geom[hydrogen] = q.geom[5] + .97*direction/np.linalg.norm(direction)
            q.name = routing_name(q)
            q.energy = p.energy - 5./constants.AUtoKCAL
            symmetry.calculate_symmetry(q)
            ts.name = p.name + f'_intra_H_migration_6_{hydrogen+1}_{int(barrier)}'
            ts.energy = p.energy + barrier/constants.AUtoKCAL
            routes.append(SimpleNamespace(instance_name=ts.name, ts=ts,
                products=[q], prod_opt=[SimpleNamespace(species=q)], do_vdW=False, mp2=0))
        # Both paths reach the same configured endpoint.
        self.assertEqual(routing_name(routes[0].products[0]), routing_name(routes[1].products[0]))
        p.reac_obj, p.reac_ts_done = routes, [-1]*3
        p.reac_type, p.reac_inst = ['intra_H_migration']*3, [None]*3
        p.reac_name = [r.instance_name for r in routes]
        return p, routes

    def test_actual_direct_summary_and_pes_selection_keep_two_routes(self):
        self.direct_summary_and_pes_selection()

    def test_ordinary_and_stereo_routes_survive_direct_summary_and_pes(self):
        self.direct_summary_and_pes_selection(ordinary=True)

    def direct_summary_and_pes_selection(self, ordinary=False):
        p, routes = self.observations()
        if ordinary:
            for route in (routes[0], routes[2]):
                route.ts.stereopath_id = 'ordinary'
                route.ts.stereopath_metadata = dict(route.ts.stereopath_metadata,
                                                    id='ordinary', site_relation='ordinary')
        base = routing_name(p)
        Path(base, 'me').mkdir(parents=True)
        par = dict(self.par, smiles='', charge=0, mult=2, me=2,
                   structure=[v for atom, xyz in zip(p.atom, p.geom)
                              for v in (str(atom), *map(float, xyz))])
        db = connect(f'{base}/kinbot.db')
        # Last-row lookup needs one electronic-energy record per named calculation.
        points = [(p, base+'_well'), (routes[0].products[0], routing_name(routes[0].products[0])+'_well')]
        points += [(r.ts, r.instance_name) for r in routes]
        for point, name in points:
            db.write(Atoms(point.atom, positions=point.geom), name=name,
                     data={'status': 'normal', 'energy': point.energy/constants.EVtoHARTREE,
                           'zpe': point.zpe, 'frequencies': point.freq})
        root = Path.cwd()
        try:
            os.chdir(base)
            MESS(dict(par, pes=0), p).write_input(None)
            direct = Path('me/mess_0000.inp').read_text()
            MESS(dict(par, pes=1), p).write_input(None)
            postprocess.create_summary_file(p, SimpleNamespace(qc='fc'), par)
        finally:
            os.chdir(root)
        def assemble(jobs, task='all', names=()):
            with patch('kinbot.pes.copy_from_kinbot'), \
                 patch('kinbot.pes.create_pesviewer_input'), patch('kinbot.pes.create_rotdpy_inputs'), \
                 patch('kinbot.pes.t1_analysis'), patch('kinbot.pes.get_l3energy', return_value=(0, -1)):
                pes.postprocess(par, jobs, task, list(names), p.mass)
            return Path('me/mess_0000.inp').read_text()
        deferred = assemble([base])
        lowest = assemble([base], 'lowestpath', [base, routes[0].products[0].name])
        for output in (direct, deferred, lowest):
            barrier = output[output.index('  Barrier'):]
            self.assertEqual(sum(line.strip().startswith('Barrier ') for line in output.splitlines()), 1)
            self.assertIn('Union ! 2 stereochemical pathways', barrier)
            self.assertEqual(values(barrier, 'ZeroEnergy'), [30., 32.])
            self.assertEqual(values(barrier, 'WellDepth'), [30., 35., 32., 37.])
            self.assertEqual(values(barrier, 'SymmetryFactor'), [1., 1.])
            self.assertEqual(barrier.count('! kinbot_stereopath'), 2)
            self.assertNotIn(routes[2].instance_name, barrier)
        # Adding two equal TS terms doubles their sum; unequal barriers must
        # instead keep their own Boltzmann factors. No arbitrary factor of 2.
        rt = .00198720425864083*800.  # R*T in kcal/mol at 800 K
        explicit = sum(math.exp(-energy/rt) for energy in values(direct[direct.index('  Barrier'):], 'ZeroEnergy'))
        self.assertAlmostEqual(explicit/math.exp(-30./rt), 1.+math.exp(-2./rt))

        # Discover the first route again from QOOH. Its slightly higher
        # observation must not become a third independent contribution.
        q = copy.deepcopy(routes[0].products[0])
        reverse = copy.copy(routes[0])
        reverse.ts = copy.deepcopy(routes[0].ts)
        reverse.ts.energy = p.energy + 33./constants.AUtoKCAL
        reverse.instance_name = reverse.ts.name = q.name + '_intra_H_migration_reverse'
        reverse.products, reverse.prod_opt = [p], [SimpleNamespace(species=p)]
        q.reac_obj, q.reac_ts_done, q.reac_type = [reverse], [-1], ['intra_H_migration']
        q.reac_inst, q.reac_name = [None], [reverse.instance_name]
        Path(q.name, 'me').mkdir(parents=True)
        db = connect(f'{q.name}/kinbot.db')
        for point, name in [(q, q.name+'_well'), (reverse.ts, reverse.instance_name)]:
            db.write(Atoms(point.atom, positions=point.geom), name=name,
                     data={'status': 'normal', 'energy': point.energy/constants.EVtoHARTREE,
                           'zpe': point.zpe, 'frequencies': point.freq})
        try:
            os.chdir(q.name)
            MESS(dict(par, pes=1), q).write_input(None)
            postprocess.create_summary_file(q, SimpleNamespace(qc='fc'), par)
        finally:
            os.chdir(root)
        output = assemble([base, q.name])
        barrier = output[output.index('  Barrier'):]
        self.assertEqual(values(barrier, 'ZeroEnergy'), [30., 32.])
        self.assertNotIn(reverse.instance_name, barrier)
        # An old summary without the classification cannot license extra flux.
        filename = Path(q.name, f'summary_{q.name}.out')
        filename.write_text('\n'.join(line for line in filename.read_text().splitlines()
                                     if not line.startswith('# kinbot_stereopath')))
        output = assemble([base, q.name])
        barrier = output[output.index('  Barrier'):]
        self.assertEqual(values(barrier, 'ZeroEnergy'), [30., 32.])
        self.assertNotIn(reverse.instance_name, barrier)
        self.assertIn('SUCCESS', filename.read_text())

    def test_mc_ensembles_are_nested_once_and_preserve_each_routes_reference(self):
        self.mc_nested_routes()

    def test_mc_ordinary_and_stereo_ensembles_are_nested_once(self):
        self.mc_nested_routes(ordinary=True)

    def mc_nested_routes(self, ordinary=False):
        p, routes = self.observations()
        if ordinary:
            for route in (routes[0], routes[2]):
                route.ts.stereopath_id = 'ordinary'
                route.ts.stereopath_metadata = dict(route.ts.stereopath_metadata, id='ordinary')
        Path('me').mkdir()
        self.par.update(multi_conf_tst=1, conformer_search=1, pes=0)
        for route in routes:
            ts = route.ts
            ts.conformer_index = [0]
            ts.conformer_geom = [ts.geom.copy()]
            ts.conformer_zeroenergy = [ts.energy+ts.zpe]
            ts.conformer_freq = [ts.freq]
        MESS(self.par, p).write_input(None)
        output = Path('me/mess_0000.inp').read_text()
        barrier = output[output.index('  Barrier'):]
        self.assertEqual(values(barrier, 'ZeroEnergy'), [30., 32.])
        self.assertEqual(barrier.count('number of species in union is 1'), 2)
        self.assertEqual(values(barrier, 'SymmetryFactor'), [1., 1.])

    def test_racemic_weights_do_not_double_the_normalized_unimolecular_rate(self):
        p, routes = self.observations()
        Path('me').mkdir()
        self.par.update(optical_population='racemic', pes=0)
        MESS(self.par, p).write_input(None)
        output = Path('me/mess_0000.inp').read_text()
        barrier = output[output.index('  Barrier'):]
        writer = MESS(self.par, p)
        parent_divisor = writer._parent_symmetry(p)
        self.assertEqual(parent_divisor, .5)
        self.assertEqual(values(barrier, 'SymmetryFactor'), [.5, .5])
        self.assertEqual([parent_divisor/d for d in values(barrier, 'SymmetryFactor')], [1., 1.])

        # With explicit mirrors inside each MC sum, every geometry has weight
        # one. Two geometries / a two-enantiomer parent gives the same weight
        # as one representative with its allowed optical factor of two.
        self.par.update(multi_conf_tst=1, conformer_search=1)
        for route in routes:
            ts = route.ts
            ts.conformer_index = [0, 1]
            ts.conformer_geom = [ts.geom.copy(), ts.geom * [-1., 1., 1.]]
            ts.conformer_zeroenergy = [ts.energy+ts.zpe]*2
            ts.conformer_freq = [ts.freq, ts.freq]
        MESS(self.par, p).write_input(None)
        output = Path('me/mess_0000.inp').read_text()
        barrier = output[output.index('  Barrier'):]
        self.assertEqual(values(barrier, 'ZeroEnergy'), [30., 30., 32., 32.])
        self.assertEqual(values(barrier, 'SymmetryFactor'), [1.]*4)
        self.assertEqual(parent_divisor * sum(1./d for d in values(barrier, 'SymmetryFactor')), 2.)

    def test_union_requires_distinct_classes_and_preserves_opposite_orientation(self):
        def block(name, key, left='w1', right='w2', energy=30.):
            label = '' if key is None else f'! kinbot_stereopath {key}\n'
            return label+f'  Barrier {name} {left} {right}\n  RRHO\n ZeroEnergy[kcal/mol] {energy}\n End\n'
        first = block('a', 'path_a')
        second = block('b', 'path_b', 'w2', 'w1', 32.)
        combined = union_stereochemical_barriers([first, second])
        self.assertEqual(len(combined), 1)
        self.assertEqual(values(combined[0], 'ZeroEnergy'), [30., 32.])
        self.assertEqual(union_stereochemical_barriers([first]), [first])
        for second in (block('b', 'path_a'), block('b', None)):
            with self.assertRaisesRegex(ValueError, 'distinct classified'):
                union_stereochemical_barriers([first, second])

    def test_existing_mirror_coverage_contract_is_used_for_classified_routes(self):
        writer = MESS(self.par, peroxy())
        # Reuse measured methanol TS/scan evidence as a coverage fixture. This
        # does not classify methanol's homotopic methyl sites as diastereotopic.
        for run, divisor in [('MeOH_H', .5), ('MeOH_H_rotor', 1.)]:
            ts = saved_point(methanol_data(run)['saddles'][0])
            ts.stereopath_metadata = {'id': 'coverage-fixture'}
            self.assertEqual(writer._parent_symmetry(ts), divisor)
            self.assertEqual(ts.mess_optical_counting['remaining_multiplier'], 1./divisor)
        ts.hir.scan_reference = None
        self.assertEqual(writer._parent_symmetry(ts), .5)
        self.assertEqual(ts.reduced_freqs, ts.freq)
        self.assertEqual(ts.rotor_projection['internal_rank'], 0)
        writer.par['rotor_scan'] = 1
        self.assertIn('kept as a harmonic oscillator', writer.make_rotors(ts, 1.))

    def test_pes_path_queries_preserve_parallel_routes(self):
        reactions = [['a', 'x', ['b'], 30.], ['a', 'y', ['b'], 32.],
                     ['b', 'z', ['c'], 40.]]
        conn, bars = pes.get_connectivity(['a', 'b', 'c'], [], reactions)
        self.assertEqual(bars[0, 1], 30.)
        paths = pes.get_all_pathways(['a', 'b', 'c'], [], reactions, ['a', 'c'], conn)
        self.assertEqual([[r[1] for r in route] for route in paths], [['x', 'z'], ['y', 'z']])


if __name__ == '__main__':
    unittest.main()
