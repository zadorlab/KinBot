"""Repeated channels must reuse a product's complete optimized state."""
import copy
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np

from kinbot.calculation import array_fingerprint
from kinbot.parameters import Parameters
from kinbot.reaction_generator import ReactionGenerator
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint


def point(smiles, mult=1):
    species = StationaryPoint('fixture', 0, mult, smiles=smiles)
    species.characterize()
    species.energy, species.zpe = -100., .01
    species.start_energy, species.start_zpe = -100., .01
    species.freq = species.reduced_freqs = [100.] * (3 * species.natom - 6)
    return species


class TestProductOptimizationReuse(unittest.TestCase):
    def test_channels_share_final_results_without_overwriting_them_with_initial_wells(self):
        # Both channels break the C-O bond of ethanol. The optimizer mock
        # selects a different geometry and records the corresponding HIR input.
        for delayed in (False, True):
            with self.subTest(second_channel_after_completion=delayed), TemporaryDirectory() as tmp:
                previous = Path.cwd()
                try:
                    os.chdir(tmp)
                    Path('input.json').write_text(json.dumps({'barrier_threshold': 100.}))
                    par = Parameters('input.json', show_warnings=False).par
                    par.update(pes=0, high_level=0, L3_calc=0)
                    parent = point('CCO')
                    fragments = [point('C[CH2]', 2), point('[OH]', 2)]
                    initial = {routing_name(p) + '_well': p for p in fragments}
                    reactions = []
                    for i in range(2):
                        endpoint = copy.copy(parent)
                        endpoint.start_multi_molecular = Mock(return_value=(
                            [copy.copy(p) for p in fragments], None))
                        reactions.append(SimpleNamespace(instance_name=f'hom_sci_{i}',
                            species=parent, instance=[1, 2], products=[], prod_opt=[],
                            do_vdW=False, irc_prod=endpoint))
                    parent.reac_obj = reactions
                    parent.reac_inst = [[1, 2], [1, 2]]
                    parent.reac_type = ['hom_sci', 'hom_sci']
                    parent.reac_name = [r.instance_name for r in reactions]
                    parent.reac_step = [0, 0]
                    parent.reac_ts_done = [2, 2]
                    qc = Mock()
                    qc.db.select.return_value = []
                    active_channel = [None]

                    def submit(species, geom):
                        # Only the first channel may advance while its
                        # optimizer is finishing in the delayed variant.
                        active_channel[0] = next(i for i, r in enumerate(reactions)
                            if any(species is p for p in r.products))

                    def read_geom(name, *args, **kwargs):
                        p = initial[name]
                        pending = (delayed and active_channel[0] == 1
                                   and parent.reac_ts_done[0] != -1)
                        result = (1 if pending else 0, p.geom.copy())
                        return (*result, p.atom.copy()) if kwargs.get('reorder') else result

                    qc.qc_opt.side_effect = submit
                    qc.get_qc_geom.side_effect = read_geom
                    qc.get_qc_energy.side_effect = lambda name: (0, initial[name].energy)
                    qc.get_qc_zpe.side_effect = lambda name: (0, initial[name].zpe)
                    optimizers = []

                    def optimize(species, *args):
                        original = species.geom.copy()
                        species.geom = original + np.array([.1, .2, .3])
                        species.energy -= .02
                        species.zpe = .02
                        species.freq = [200.] * len(species.freq)
                        species.reduced_freqs = list(species.freq)
                        species.source_job = 'conf/' + routing_name(species) + '_0003'
                        species.source_row_id = 123
                        species.hir = SimpleNamespace(species=species,
                            scan_reference={'geometry_sha256': array_fingerprint(species.geom)})
                        opt = SimpleNamespace(species=species, shir=1, shigh=1,
                                              do_optimization=Mock())
                        optimizers.append(opt)
                        return opt

                    with patch('kinbot.reaction_generator.Optimize', side_effect=optimize), \
                            patch('kinbot.reaction_generator.time.sleep'), \
                            patch('kinbot.reaction_generator.postprocess.createPESViewerInput'):
                        ReactionGenerator(parent, par, qc, 'input.json').generate()

                    self.assertEqual(parent.reac_ts_done, [-1, -1])
                    self.assertEqual(len(optimizers), 2)  # one per distinct product
                    for j, opt in enumerate(optimizers):
                        for reaction in reactions:
                            self.assertIs(reaction.products[j], opt.species)
                            self.assertIs(reaction.prod_opt[j], opt)
                        p = opt.species
                        np.testing.assert_allclose(p.geom, fragments[j].geom + [.1, .2, .3])
                        self.assertEqual(p.energy, -100.02)
                        self.assertEqual(p.zpe, .02)
                        self.assertEqual(p.freq, [200.] * len(p.freq))
                        self.assertEqual(p.hir.scan_reference['geometry_sha256'],
                                         array_fingerprint(p.geom))
                        self.assertEqual(p.source_row_id, 123)
                finally:
                    os.chdir(previous)

    def test_initial_copies_preserve_identical_fragment_multiplicity(self):
        fragment = point('[CH3]', 2)
        copies = ReactionGenerator.initial_product_copies([fragment, fragment])
        self.assertIs(copies[0], copies[1])
        self.assertIsNot(copies[0], fragment)
        copies[0].geom[0, 0] += .2
        copies[0].bond[:] = 0
        self.assertNotEqual(copies[0].geom[0, 0], fragment.geom[0, 0])
        self.assertTrue(np.any(fragment.bond))


if __name__ == '__main__':
    unittest.main()
