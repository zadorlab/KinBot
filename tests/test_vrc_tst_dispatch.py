"""VRC correction jobs remain distinct and complete before rotdPy handoff."""

import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

from ase import Atoms
from ase.db import connect
import numpy as np

from kinbot.vrc_tst_scan import VTS


class _FakeMolpro:
    def __init__(self, species, parameters):
        self.species = species

    def create_molpro_input(self, name='', VTS=False, sample=False):
        assert VTS
        Path(f'vrctst/molpro/{name}.inp').write_text(
            f'***,{name}\n! sample={sample}\n')

    def get_molpro_energy(self, key, name='', VTS=False):
        assert key == 'MYENERGY'
        assert VTS
        if not Path(f'vrctst/molpro/{name}.out').is_file():
            return 0, -1.
        return 1, (-10. if name.endswith('_fr') else -11.)


def test_noscan_uses_distinct_sampling_and_high_level_asymptotes():
    with TemporaryDirectory() as temporary:
        previous = Path.cwd()
        os.chdir(temporary)
        try:
            Path('vrctst/molpro').mkdir(parents=True)
            atoms = Atoms('H2', positions=[[0., 0., 0.], [0., 0., 30.]])
            connect('kinbot.db').write(
                atoms, name='vrctst/r_vts_pt_asymptote_fr',
                data={'status': 'normal'})
            well = SimpleNamespace(
                chemid='parent', geom=np.asarray([[0., 0., 0.],
                                                   [0., 0., .74]]))
            parameters = {
                'vrc_tst_scan_points': [2.5],
                'vrc_tst_scan_molpro_key': 'MYENERGY',
            }
            vts = VTS(well, parameters, None)
            product_a = SimpleNamespace(
                atom=['H'], geom=[[0., 0., 0.]], mult=2)
            product_b = SimpleNamespace(
                atom=['H'], geom=[[0., 0., 30.]], mult=2)
            vts.scan_reac['r'] = SimpleNamespace(
                usym=[[[0]], [[0]]], equiv=[[0], [1]],
                maps=[np.asarray([0]), np.asarray([1])],
                products=[product_a, product_b])
            observed = set()

            def dispatch(pending):
                observed.update(pending)
                for name in pending:
                    Path(f'vrctst/molpro/{name}.out').write_text(
                        'Molpro calculation terminated\n')

            with patch('kinbot.vrc_tst_scan.Molpro', _FakeMolpro), \
                    patch.object(vts, '_dispatch_molpro_corrections', dispatch):
                vts.energies(['r'], noscan=True)

            assert observed == {
                'r_vts_pt_asymptote_fr', 'r_vts_pt_asymptote'}
            correction = json.loads(
                Path('vrctst/corr_r.json').read_text())
            assert correction['dist'] == [30]
            assert correction['e_samp'] == [0.]
            assert correction['e_high'] == [0.]
            assert correction['e_inf_samp'] == -10.
            assert correction['e_inf_high'] == -11.
        finally:
            os.chdir(previous)
