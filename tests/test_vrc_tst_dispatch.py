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
import pytest

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
                'vrc_tst_sample_method': 'caspt2(2,2)',
                'vrc_tst_sample_basis': 'vdz',
                'vrc_tst_high_method': 'mrci+q(2,2)',
                'vrc_tst_high_basis': 'vtz',
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
                    source = Path(f'vrctst/molpro/{name}.inp')
                    Path(f'vrctst/molpro/{name}.input.sha256').write_text(
                        __import__('hashlib').sha256(
                            source.read_bytes()).hexdigest() + '\n')

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
            assert correction['levels']['trusted_correction'] == {
                'method': 'mrci+q(2,2)', 'basis': 'vtz'}
        finally:
            os.chdir(previous)


def test_molpro_dispatch_normalizes_numpy_charge_and_multiplicity():
    with TemporaryDirectory() as temporary:
        previous = Path.cwd()
        os.chdir(temporary)
        try:
            Path('vrctst/molpro').mkdir(parents=True)
            Path('vrctst/molpro/r.inp').write_text(
                '***,VRC correction\nMolpro calculation terminated\n')
            species = SimpleNamespace(
                atom=['H', 'H'],
                geom=np.asarray([[0., 0., 0.], [0., 0., .74]]),
                charge=np.float64(0.), mult=np.int64(1))
            parameters = {
                'vrc_tst_walltime': '00:10:00',
                'single_point_ppn': 4,
                'vrc_tst_min_stack_mw': 50,
                'queue_name': 'short',
                'vrc_tst_max_nodes': 1,
            }
            vts = VTS(SimpleNamespace(), parameters, None)
            captured = {}

            def prepare_dispatch(spec, run_dir):
                from kinbot.anl.dispatch import _atoms
                _atoms(spec['molecule'])
                json.dumps(spec)
                captured['spec'] = spec
                output = run_dir / 'tasks' / 'r' / 'r.out'
                output.parent.mkdir(parents=True)
                output.write_text('Molpro calculation terminated\n')

            with patch('kinbot.anl.dispatch.prepare', prepare_dispatch), \
                    patch('kinbot.anl.dispatch.preflight'), \
                    patch('kinbot.vrc_tst_scan.subprocess.run',
                          return_value=SimpleNamespace(returncode=0)):
                vts._dispatch_molpro_corrections({'r': species})

            molecule = captured['spec']['molecule']
            assert type(molecule['charge']) is int
            assert type(molecule['multiplicity']) is int
            assert Path('vrctst/molpro/r.out').read_text() == \
                'Molpro calculation terminated\n'
        finally:
            os.chdir(previous)


def test_molpro_dispatch_surfaces_native_task_failure():
    with TemporaryDirectory() as temporary:
        previous = Path.cwd()
        os.chdir(temporary)
        try:
            Path('vrctst/molpro').mkdir(parents=True)
            Path('vrctst/molpro/r.inp').write_text('***,failed VRC point\n')
            species = SimpleNamespace(
                atom=['H', 'H'], geom=[[0., 0., 0.], [0., 0., .74]],
                charge=0., mult=1)
            parameters = {
                'vrc_tst_walltime': '00:10:00', 'single_point_ppn': 2,
                'vrc_tst_min_stack_mw': 50, 'queue_name': 'short',
                'vrc_tst_max_nodes': 1,
            }
            vts = VTS(SimpleNamespace(), parameters, None)
            state = {'tasks': {'r': {'status': 'failed'}}}

            def prepare_dispatch(_, run_dir):
                task_dir = run_dir / 'tasks' / 'r'
                task_dir.mkdir(parents=True)
                (task_dir / 'execution.json').write_text(json.dumps({
                    'status': 'failed',
                    'error': 'PMPI_Init: OFI endpoint open failed',
                }))

            with patch('kinbot.anl.dispatch.prepare', prepare_dispatch), \
                    patch('kinbot.anl.dispatch.preflight'), \
                    patch('kinbot.anl.dispatch._load',
                          side_effect=lambda run_dir: (Path(run_dir), {}, state)), \
                    patch('kinbot.vrc_tst_scan.subprocess.run',
                          return_value=SimpleNamespace(returncode=1)):
                with pytest.raises(RuntimeError, match='PMPI_Init.*OFI endpoint'):
                    vts._dispatch_molpro_corrections({'r': species})
        finally:
            os.chdir(previous)
