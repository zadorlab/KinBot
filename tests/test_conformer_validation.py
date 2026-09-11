"""Exercise completed conformer records without submitting QC calculations."""
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock

from ase.build import molecule

from kinbot import constants
from kinbot.conformers import Conformers
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


class TestConformerValidation(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(directory.name)
        Path('conf').mkdir()
        Path('input.json').write_text('{"barrier_threshold": 50}')
        self.par = Parameters('input.json', show_warnings=False).par
        atoms = molecule('CH3CH2OH')
        self.species = StationaryPoint('ethanol', 0, 1,
            atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        self.species.characterize()
        self.species.name = str(self.species.chemid)
        self.qc = QuantumChemistry(self.par)
        self.qc.qc = 'fc'  # Read database results without native log copying.
        self.qc.check_qc = Mock(return_value='normal')
        self.atoms = atoms
        self.nfreq = 3 * len(atoms) - 6

    def record(self, name, energy, frequencies, zpe=0.):
        self.qc.db.write(self.atoms, name=name, data={
            'energy': energy / constants.EVtoHARTREE, 'zpe': zpe,
            'frequencies': frequencies, 'status': 'normal'})

    def check_pair(self, wellorts, higher_frequencies):
        self.species.wellorts = wellorts
        if wellorts:
            self.species.name = 'ethanol_ts_fixture'
        parent = self.species.name if wellorts else f'{self.species.name}_well'
        valid = ([-1000.] + [100.] * (self.nfreq - 1)
                 if wellorts else [100.] * self.nfreq)
        self.species.freq = valid
        search = Conformers(self.species, self.par, self.qc)
        self.record(parent, -100., valid)
        self.record(search.get_job_name(0), -100., valid)
        self.record(search.get_job_name(1), -99.999, higher_frequencies)
        search.conf = 2
        result = search.check_conformers()
        return search, result

    def test_higher_well_with_imaginary_frequency_is_excluded_before_population(self):
        search, result = self.check_pair(0, [-100.] + [100.] * (self.nfreq - 1))
        self.assertEqual(result[7], [0, 1])
        self.assertEqual(search.find_unique(*result[4:8], temp=300., boltz=.001)[-1], [0])

    def test_ts_requires_one_reaction_coordinate_with_existing_soft_mode_allowance(self):
        for initial, accepted in [([], False), ([-1000.], True),
                                  ([-1000., -100.], False),
                                  ([-1000., -20.], True),
                                  ([-1000., -20., -10.], False)]:
            with self.subTest(initial=initial):
                frequencies = initial + [100.] * (self.nfreq - len(initial))
                _, result = self.check_pair(1, frequencies)
                self.assertEqual(result[7], [0, 0 if accepted else 1])

    def test_empty_and_nonfinite_spectra_are_excluded(self):
        for wellorts in (0, 1):
            for frequencies in ([], [float('nan')] * self.nfreq):
                with self.subTest(wellorts=wellorts, frequencies=frequencies):
                    _, result = self.check_pair(wellorts, frequencies)
                    self.assertEqual(result[7], [0, 1])

    def test_existing_small_imaginary_well_mode_is_still_allowed(self):
        _, result = self.check_pair(0, [-20.] + [100.] * (self.nfreq - 1))
        self.assertEqual(result[7], [0, 0])


if __name__ == '__main__':
    unittest.main()
