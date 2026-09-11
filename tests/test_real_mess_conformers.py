"""Replay saved THF calculations through the complete direct and PES writers.

The fixture contains two actual isolated-conformer records, not an invented
reaction or a validated kinetic model. No QC job or MESS executable is run.
"""
import json
import logging
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from ase import Atoms

from kinbot import constants, pes, symmetry
from kinbot.mess import MESS
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint

FIXTURE = Path(__file__).parent / 'reference' / 'thf_mc_rrho.json'


def ground_energies(text):
    return [float(line.split()[1]) for line in text.splitlines()
            if line.strip().startswith('ZeroEnergy[kcal/mol]')]


def write_example(directory):
    """Save both final MESS inputs, preserving the fixture's calculated values."""
    fixture = json.loads(FIXTURE.read_text())
    records = fixture['conformers']
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    previous = Path.cwd()
    previous_logger = getattr(pes, 'logger', None)
    try:
        os.chdir(directory)
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'smiles': fixture['smiles'],
            'charge': 0, 'mult': 1, 'high_level': 0, 'qc': 'fc',
            'fc_model_path': fixture['model'],
            'method': fixture['model'], 'high_level_method': fixture['model'],
            'basis': '', 'high_level_basis': '',
            'multi_conf_tst': 1, 'conformer_search': 1, 'rotor_scan': 0,
            'me': 0, 'uq': 0, 'epsilon': 100., 'sigma': 3.,
        }))
        par = Parameters('input.json', show_warnings=False).par
        parent = records[0]
        species = StationaryPoint('THF', 0, 1, atom=parent['atoms'],
                                  geom=np.array(parent['geometry_angstrom']))
        species.characterize()
        symmetry.calculate_symmetry(species)
        species.energy = parent['electronic_energy_ev'] * constants.EVtoHARTREE
        species.zpe = parent['zpe_hartree']
        species.freq = list(parent['frequencies_cm-1'])
        species.reduced_freqs = list(species.freq)
        species.conformer_index = [0, 1]
        species.conformer_geom = [np.array(row['geometry_angstrom']) for row in records]
        species.conformer_freq = [row['frequencies_cm-1'] for row in records]
        species.conformer_energy = [row['electronic_energy_ev'] * constants.EVtoHARTREE
                                    for row in records]
        species.conformer_zeroenergy = [energy + row['zpe_hartree']
                                        for energy, row in zip(species.conformer_energy, records)]
        Path('me').mkdir(exist_ok=True)
        par['pes'] = 0
        writer = MESS(par, species)
        writer.write_input(None)
        direct = Path('me/mess_0000.inp').read_text()
        Path('direct_mess.inp').write_text(direct)
        par['pes'] = 1
        block = writer.write_well(species, 0., 1., 0)
        name = str(species.chemid)
        Path(name).mkdir(exist_ok=True)
        Path(f'{name}/{name}_0000.mess').write_text(block)
        pes.logger = logging.getLogger('KinBot')
        mass = Atoms(parent['atoms']).get_masses().sum()
        pes.create_mess_input(par, [name], [], [], [], [],
                              {name: 0.}, {}, {name: name}, mass, False)
        combined = Path('me/mess_0000.inp').read_text()
        Path('pes_mess.inp').write_text(combined)
        return direct, combined, species.conformer_zeroenergy
    finally:
        os.chdir(previous)
        if previous_logger is None and hasattr(pes, 'logger'):
            delattr(pes, 'logger')
        elif previous_logger is not None:
            pes.logger = previous_logger


class TestRealConformerOutput(unittest.TestCase):
    def test_saved_thf_member_grounds_reach_both_final_inputs(self):
        with TemporaryDirectory() as directory:
            direct, combined, zero = write_example(directory)
        expected = [round((energy - zero[0]) * constants.AUtoKCAL, 2)
                    for energy in zero]
        self.assertGreater(expected[1], 0.)
        self.assertEqual(ground_energies(direct), expected)
        self.assertEqual(ground_energies(combined), expected)
        for output in (direct, combined):
            self.assertEqual(output.count('End ! RRHO'), 2)
            self.assertEqual(output.count('Frequencies[1/cm]'), 2)
            self.assertNotIn('{zeroenergy}', output)


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--write-example':
        write_example(sys.argv[2])
    else:
        unittest.main()
