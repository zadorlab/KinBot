import os
import re

import numpy as np

from ase.units import Hartree, Bohr
from ase.io.orca import write_orca
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError


class ORCA(FileIOCalculator):
    implemented_properties = ['energy', 'forces']

    if 'ORCA_COMMAND' in os.environ:
        command = os.environ['ORCA_COMMAND'] + ' PREFIX.inp > PREFIX.out'
    else:
        command = 'orca PREFIX.inp > PREFIX.out'

    default_parameters = dict(
        charge=0, mult=1,
        task='gradient',
        orcasimpleinput='tightscf PBE def2-SVP',
        orcablocks='%scf maxiter 200 end')

    def __init__(self, restart=None,
                 ignore_bad_restart_file=FileIOCalculator._deprecated,
                 label='orca', atoms=None, **kwargs):
        """ ASE interface to ORCA 4
        by Ragnar Bjornsson, Based on NWchem interface but simplified.
        Only supports energies and gradients (no dipole moments,
        orbital energies etc.) for now.

        For more ORCA-keyword flexibility, method/xc/basis etc.
        keywords are not used. Instead, two keywords:

            orcasimpleinput: str
                What you'd put after the "!" in an orca input file.

            orcablock: str
                What you'd put in the "% ... end"-blocks.

        are used to define the ORCA simple-inputline and the ORCA-block input.
        This allows for more flexible use of any ORCA method or keyword
        available in ORCA instead of hardcoding stuff.

        Point Charge IO functionality added by A. Dohn.
        """
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        self.pcpot = None

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        p.write(self.label + '.ase')
        p['label'] = self.label
        if self.pcpot:  # also write point charge file and add things to input
            p['pcpot'] = self.pcpot

        write_orca(atoms, **p)

    def read(self, label):
        FileIOCalculator.read(self, label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        with open(self.label + '.inp') as fd:
            for line in fd:
                if line.startswith('geometry'):
                    break
            symbols = []
            positions = []
            for line in fd:
                if line.startswith('end'):
                    break
                words = line.split()
                symbols.append(words[0])
                positions.append([float(word) for word in words[1:]])

        self.parameters = Parameters.read(self.label + '.ase')
        self.read_results()

    def read_results(self):
        self.read_energy()
        if self.parameters.task.find('gradient') > -1:
            self.read_forces()

    def read_energy(self):
        """Read Energy from ORCA output file."""
        with open(self.label + '.out', mode='r', encoding='utf-8') as fd:
            text = fd.read()
        # Energy:
        re_energy = re.compile(r"FINAL SINGLE POINT ENERGY.*\n")
        re_not_converged = re.compile(r"Wavefunction not fully converged")
        found_line = re_energy.search(text)
        if found_line and not re_not_converged.search(found_line.group()):
            self.results['energy'] = float(found_line.group().split()[-1]) * Hartree

    def read_forces(self):
        """Read Forces from ORCA output file."""
        with open(f'{self.label}.engrad', 'r') as fd:
            lines = fd.readlines()
        getgrad = False
        gradients = []
        tempgrad = []
        for i, line in enumerate(lines):
            if line.find('# The current gradient') >= 0:
                getgrad = True
                gradients = []
                tempgrad = []
                continue
            if getgrad and "#" not in line:
                grad = line.split()[-1]
                tempgrad.append(float(grad))
                if len(tempgrad) == 3:
                    gradients.append(tempgrad)
                    tempgrad = []
            if '# The at' in line:
                getgrad = False
        self.results['forces'] = -np.array(gradients) * Hartree / Bohr

    def embed(self, mmcharges=None, **parameters):
        """Embed atoms in point-charges (mmcharges)
        """
        self.pcpot = PointChargePotential(mmcharges, label=self.label)
        return self.pcpot


class PointChargePotential:
    def __init__(self, mmcharges, label=None, positions=None, directory=None):
        """ Point Charge Potential Interface to ORCA """
        if positions is not None:
            self.set_positions(positions)
        if directory is None:
            directory = os.getcwd()

        self.directory = directory + os.sep
        self.mmcharges = mmcharges
        self.label = label

    def set_positions(self, positions):
        self.positions = positions

    def set_charges(self, mmcharges):
        self.q_p = mmcharges

    def write_mmcharges(self, filename):
        pc_file = open(os.path.join(self.directory,
                                    filename + '.pc'), 'w')

        pc_file.write('{0:d}\n'.format(len(self.mmcharges)))
        for [pos, pc] in zip(self.positions, self.mmcharges):
            [x, y, z] = pos
            pc_file.write('{0:12.6f} {1:12.6f} {2:12.6f} {3:12.6f}\n'
                          .format(pc, x, y, z))

        pc_file.close()

    def get_forces(self, calc):
        ''' reads forces on point charges from .pcgrad file '''
        with open(os.path.join(self.directory, self.label + '.pcgrad'),
                  'r', encoding='utf-8') as fd:
            lines = fd.readlines()
        numpc = int(lines[0])
        forces = np.zeros((numpc, 3))
        for i in range(numpc):
            [fx, fy, fz] = [float(f) for f in lines[i + 1].split()]
            forces[i, :] = fx, fy, fz

        return -forces * Hartree / Bohr
