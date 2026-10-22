"""Two-surface selection, native MESS agreement and MC/capture counting."""
import copy
from contextlib import contextmanager
import os
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from kinbot import constants
from kinbot.bimol_mess import BimolecularMESS
from kinbot.bimol_two_ts import header, read_rates, variational
from tests.test_bimol_mess import add_complex, methanol_records, parameters


@contextmanager
def workdir(path):
    previous = Path.cwd()
    path.mkdir(exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def fixture(par, *, rc=5., pc=2., mc=False, hir=False):
    records = methanol_records(mc=mc, hir=hir)
    records['rate_models']['pathways'] = records['rate_models']['pathways'][:1]
    writer = BimolecularMESS(par, records)
    path = records['rate_models']['pathways'][0]
    for side, depth in (('reactant', rc), ('product', pc)):
        if depth is None:
            continue
        identifier = side + '_complex'
        zero = sum(writer.ground(d) for d in writer.nodes[path[side]]['fragment_data'])
        add_complex(records, identifier, 'rc' if side == 'reactant' else 'pc', path[side],
                    zero - depth / constants.AUtoKCAL)
        path[identifier] = identifier
    return records


def set_product_energy(records, par, relative_kcal):
    """Shift one fragment's complete energy description for an analytic case."""
    writer = BimolecularMESS(par, records)
    endpoint = records['rate_models']['pathways'][0]['product']
    node = next(n for n in records['nodes'] if n['id'] == endpoint)
    current = sum(writer.ground(d) for d in node['fragment_data'])
    shift = writer.reference + relative_kcal / constants.AUtoKCAL - current
    data = node['fragment_data'][0]
    data['energy'] += shift
    thermo = data['thermochemistry']
    thermo['electronic_energy_hartree'] += shift
    members = list(thermo.get('conformer_inventory', []))
    if thermo.get('rrho_representative_counting'):
        members.append(thermo['rrho_representative_counting'])
    for member in members:
        for key in ('electronic_energy_hartree', 'zero_energy_hartree'):
            if member.get(key) is not None:
                member[key] += shift


class TestTwoTSSelection(unittest.TestCase):
    def test_asymptote_selection_ignores_binding_depth_and_reaction_direction(self):
        with TemporaryDirectory() as directory:
            par = parameters(directory)
            for mc in (False, True):
                for rc, pc in ((5., 2.), (2., 5.), (None, 5.), (5., None)):
                    for energy, expected in ((-1., 'reactant'), (1., 'product')):
                        with self.subTest(mc=mc, rc=rc, pc=pc, energy=energy):
                            records = fixture(par, rc=rc, pc=pc, mc=mc)
                            set_product_energy(records, par, energy)
                            writer = BimolecularMESS(par, records)
                            path = writer.models['pathways'][0]
                            original = writer._outer_choice(path)
                            self.assertEqual(original['side'], expected)
                            self.assertAlmostEqual(original['asymptote_energy_kcal_mol'], max(0., energy))
                            reverse = dict(path, reactant=path['product'], product=path['reactant'],
                                           reactant_complex=path['product_complex'], product_complex=path['reactant_complex'])
                            self.assertEqual(writer._outer_choice(reverse)['endpoint'], original['endpoint'])
                            self.assertEqual(writer._depths(reverse), writer._depths(path)[::-1])
            records = fixture(par, rc=2., pc=5.)
            set_product_energy(records, par, 0.)
            writer = BimolecularMESS(par, records)
            path = writer.models['pathways'][0]
            reverse = dict(path, reactant=path['product'], product=path['reactant'],
                           reactant_complex=path['product_complex'], product_complex=path['reactant_complex'])
            self.assertEqual(writer._outer_choice(path)['endpoint'], writer._outer_choice(reverse)['endpoint'])
            path.update(reactant_complex=None, product_complex=None)
            self.assertIsNone(writer._outer_choice(path))

    def test_selection_uses_zpe_and_not_only_electronic_energy(self):
        with TemporaryDirectory() as directory:
            par = parameters(directory)
            records = fixture(par)
            set_product_energy(records, par, -1.)
            writer = BimolecularMESS(par, records)
            path = writer.models['pathways'][0]
            self.assertEqual(writer._outer_choice(path)['side'], 'reactant')
            point = writer.nodes[path['product']]['fragment_data'][0]
            point['zpe'] += 2 / constants.AUtoKCAL
            point['thermochemistry']['zpe_hartree'] = point['zpe']
            self.assertEqual(writer._outer_choice(path)['side'], 'product')

    def test_no_well_models_require_no_complex_frequencies_or_symmetry(self):
        with TemporaryDirectory() as directory, workdir(Path(directory)):
            par = parameters(directory)
            par['me_skip_vdW'] = 1
            records = fixture(par, rc=2., pc=5.)
            for complex_model in records['rate_models']['complexes']:
                complex_model['stationary_point'].pop('thermochemistry')
            writer = BimolecularMESS(par, records)
            writer.write_input()
            self.assertFalse(writer.warnings)
            self.assertEqual(writer.mess_jobs[0]['kind'], 'two_ts')
            selected = writer.mess_jobs[0]['channels'][0]['outer_pst_selection'][0]
            self.assertEqual(selected['side'], 'reactant')
            text = Path('me/mess_two_ts_0000_states.inp').read_text()
            self.assertNotIn('Well reactant_complex', text)
            self.assertNotIn('Well product_complex', text)
            self.assertEqual(text.count('WellDepth[kcal/mol]'), 2)

    def test_higher_asymptote_need_not_have_its_own_complex(self):
        with TemporaryDirectory() as directory, workdir(Path(directory)):
            par = parameters(directory)
            par['me_skip_vdW'] = 1
            writer = BimolecularMESS(par, fixture(par, rc=None, pc=5.))
            writer.write_input()
            self.assertFalse(writer.warnings)
            self.assertEqual(writer.mess_jobs[0]['kind'], 'two_ts')
            selection = writer.mess_jobs[0]['channels'][0]['outer_pst_selection'][0]
            self.assertEqual(selection['side'], 'reactant')
            self.assertIsNone(selection['complex'])
            self.assertIn('higher separated-fragment E + ZPE', Path('me/mess_two_ts_0000_states.inp').read_text())

    def test_vdw_on_keeps_explicit_wells_and_three_surfaces(self):
        with TemporaryDirectory() as directory, workdir(Path(directory)):
            par = parameters(directory)
            writer = BimolecularMESS(par, fixture(par))
            writer.write_input()
            self.assertFalse(writer.warnings)
            self.assertEqual(writer.mess_jobs[0]['kind'], 'master_equation')
            text = Path('me/mess_0000.inp').read_text()
            self.assertEqual(text.count('Core PhaseSpaceTheory'), 2)
            self.assertIn('Well reactant_complex', text)
            self.assertIn('Well product_complex', text)

    def test_input_only_queue_script_runs_the_state_and_rate_stages(self):
        with TemporaryDirectory() as directory, workdir(Path(directory)):
            par = parameters(directory)
            par.update(me_skip_vdW=1, run_me=0, queuing='slurm', queue_name='local')
            writer = BimolecularMESS(par, fixture(par))
            writer.write_input()
            writer.run()
            script = next(Path('me').glob('run_mess_two_ts_*.sbatch')).read_text()
            self.assertIn('kinbot.bimol_two_ts', script)
            self.assertIn('--state-command', script)
            self.assertIn('--mess-command', script)


@unittest.skipUnless(shutil.which('mess') and shutil.which('mess-abstraction'), 'MESS executables required')
class TestTwoTSMESS(unittest.TestCase):
    def test_reversing_reaction_preserves_rates_and_outer_selection(self):
        with TemporaryDirectory() as directory, workdir(Path(directory)):
            par = parameters(directory)
            par.update(me_skip_vdW=1, ModelEnergyLimit=100.)
            forward = fixture(par, rc=2., pc=5., hir=True)
            reverse = copy.deepcopy(forward)
            path = reverse['rate_models']['pathways'][0]
            path['reactant'], path['product'] = path['product'], path['reactant']
            path['reactant_complex'], path['product_complex'] = path['product_complex'], path['reactant_complex']
            reverse['reference_node'] = path['reactant']
            outputs, choices = [], []
            for i, records in enumerate((forward, reverse)):
                with workdir(Path(directory) / str(i)):
                    writer = BimolecularMESS(par, records)
                    writer.write_input()
                    writer.run()
                    job = writer.mess_jobs[0]
                    choices.append(job['outer_pst_endpoint'])
                    outputs.append(read_rates(Path('me') / job['channels'][0]['rate_output'], par['TemperatureList']))
            self.assertEqual(choices[0], choices[1])
            np.testing.assert_allclose(outputs[0], outputs[1][:, [0, 2, 1]], rtol=0.005)

    def test_single_ts_matches_native_variational_including_tunneling(self):
        with TemporaryDirectory() as directory, workdir(Path(directory)):
            par = parameters(directory)
            par.update(me_skip_vdW=1, ModelEnergyLimit=100.)
            for hir, rc, pc, product_energy in ((False, 5., 2., -1.), (True, 2., 5., 1.00173)):
                with self.subTest(hir=hir, rc=rc), workdir(Path(directory) / str(hir)):
                    records = fixture(par, rc=rc, pc=pc, hir=hir)
                    set_product_energy(records, par, product_energy)
                    writer = BimolecularMESS(par, records)
                    writer.write_input()
                    writer.run()
                    job = writer.mess_jobs[0]
                    path = writer.models['pathways'][0]
                    outer = writer.nodes[path[writer._outer_choice(path)['side']]]['fragment_data']
                    import json
                    config = json.loads(Path('me', job['stem'] + '_states.json').read_text())
                    self.assertEqual(config['outer_threshold_cm1'], 0.)
                    self.assertAlmostEqual(config['energy_reference_hartree'], sum(writer.ground(d) for d in outer))
                    text = header('native', par['TemperatureList'], par['ModelEnergyLimit'])
                    text += ''.join(writer.bimolecular(dict(writer.nodes[path[s]], id=n), writer.reference)
                                    for s, n in (('reactant', 'R'), ('product', 'P')))
                    text += 'Barrier native\n' + variational(writer.rrho(path['stationary_point'], writer.reference,
                        depths=writer._depths(path)), outer=writer.pst(outer, writer.reference, path['stationary_point']['mult']))
                    Path('native.inp').write_text(text)
                    subprocess.run([par['mess_abstraction_command'], 'native.inp'], check=True,
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    actual = read_rates(Path('me') / job['channels'][0]['rate_output'], par['TemperatureList'])
                    native = read_rates('native.out', par['TemperatureList'])
                    np.testing.assert_allclose(actual, native, rtol=0.005)

    def test_mc_sum_equals_identical_weighted_ts_without_duplicate_capture(self):
        with TemporaryDirectory() as directory, workdir(Path(directory)):
            par = parameters(directory)
            par.update(me_skip_vdW=1, ModelEnergyLimit=100.)
            original = fixture(par)
            writer = BimolecularMESS(par, original)
            ts = original['rate_models']['pathways'][0]['stationary_point']
            member = copy.deepcopy(writer.members(ts)[0])
            mc = copy.deepcopy(original)
            thermo = mc['rate_models']['pathways'][0]['stationary_point']['thermochemistry']
            thermo.update(conformer_representation='MC-RRHO', hir={'rotors': []},
                conformer_inventory=[dict(member, retained=True, status='valid', member_id=f'fixture_{i}') for i in range(2)])
            ts['thermochemistry']['optical_counting']['remaining_multiplier'] *= 2
            observed = []
            for name, records in (('weighted', original), ('mc', mc)):
                with workdir(Path(directory) / name):
                    current = BimolecularMESS(par, records)
                    current.write_input()
                    current.run()
                    job = current.mess_jobs[0]
                    observed.append(read_rates(Path('me') / job['channels'][0]['rate_output'], par['TemperatureList']))
            np.testing.assert_allclose(observed[0], observed[1], rtol=2e-4)

    def test_shared_entrance_keeps_channels_and_detailed_balance(self):
        with TemporaryDirectory() as directory, workdir(Path(directory)):
            par = parameters(directory)
            par.update(me_skip_vdW=1, ModelEnergyLimit=100.)
            records = methanol_records()
            base = BimolecularMESS(par, records)
            add_complex(records, 'rc', 'rc', records['reference_node'], base.reference - 5 / constants.AUtoKCAL)
            for path in records['rate_models']['pathways']:
                path['reactant_complex'] = 'rc'
            writer = BimolecularMESS(par, records)
            writer.write_input()
            writer.run()
            self.assertEqual(len(writer.mess_jobs), 1)
            job = writer.mess_jobs[0]
            self.assertEqual(len(job['channels']), 2)
            for i, channel in enumerate(job['channels']):
                path = next(p for p in writer.models['pathways'] if p['id'] == channel['pathways'][0])
                text = header(f'direct{i}', par['TemperatureList'], par['ModelEnergyLimit'])
                text += ''.join(writer.bimolecular(dict(writer.nodes[key], id=name), writer.reference)
                                for name, key in channel['endpoints'].items())
                text += 'Barrier direct\n' + writer.rrho(path['stationary_point'], writer.reference, depths=writer._depths(path))
                Path(f'direct{i}.inp').write_text(text)
                subprocess.run([par['mess_abstraction_command'], f'direct{i}.inp'], check=True,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                direct = read_rates(f'direct{i}.out', par['TemperatureList'])
                actual = read_rates(Path('me') / channel['rate_output'], par['TemperatureList'])
                # Barrier changes alter the two rates equally; their ratio is
                # fixed by the same independent separated-fragment models.
                np.testing.assert_allclose(actual[:, 1] / actual[:, 2], direct[:, 1] / direct[:, 2], rtol=3e-5)
            outer = np.loadtxt(Path('me') / (job['stem'] + '_outer.dat'))
            from kinbot.bimol_two_ts import read_states
            import json
            config = json.loads(Path('me', job['stem'] + '_states.json').read_text())
            energy, combined = read_states(Path('me', job['stem'] + '_combined.micro'), config['grid'], 1)
            bound = np.interp(energy, outer[:, 0], outer[:, 1], left=0.)
            self.assertTrue(np.all(combined[:, 0] <= bound * 1.00002))


if __name__ == '__main__':
    unittest.main()
