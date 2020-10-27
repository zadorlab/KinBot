import os
import pkg_resources
import numpy as np
import logging

from kinbot import constants


class Molpro:
    """
    Class to write and read molpro file and to run molpro
    """
    def __init__(self, species, par):
        self.species = species
        self.par = par

    def create_molpro_input(self, bls=0, name='', shift_vec=None, natom1=None):
        """
        Create the input for molpro based on the template,
        which is either the one in the system, or provided
        by the user.
        : bls is whether it is a bls rection
        : name is an alternative name for the file
        : shift_vec is for bls to define the direction of shift in prod scan
        : natom1 is the number of atoms in the second fragment
        """
        if bls == 1 and shift_vec is None:
            tpl_file = self.par['barrierless_saddle_single_point_template']
        elif bls == 1 and shift_vec is not None:
            tpl_file = self.par['barrierless_saddle_prod_single_point_template']
        elif self.par['single_point_template'] == '':
            tpl_file = pkg_resources.resource_filename('tpl', 'molpro.tpl')
        else:
            tpl_file = self.par['single_point_template']

        with open(tpl_file) as f:
            file = f.read()

        fname = self.get_name(name)

        geom = ''
        nelectron = 0
        for i, at in enumerate(self.species.atom):
            x, y, z = self.species.geom[i]
            geom += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
            nelectron += constants.znumber[at]

        nelectron -= self.species.charge
        symm = self.molpro_symm()
        spin = self.species.mult - 1

        if bls == 0:
            with open('molpro/' + fname + '.inp', 'w') as outf:
                outf.write(file.format(name=fname,
                                       natom=self.species.natom,
                                       geom=geom,
                                       nelectron=nelectron,
                                       symm=symm,
                                       spin=spin,
                                       charge=self.species.charge
                                       ))

        else:
            closed = (nelectron - self.par['barrierless_saddle_nelectron']) / 2
            if closed.is_integer() is not True:
                logging.warning("The number of closed orbitals is not an integer,\n\
                             the CASPT2-like calculation will crash, but\n\
                             KinBot carries on for now. Revise your input,\n\
                             barrierless_saddle_nelectron is incorrect.")
            else:
                closed = int(closed)

            occ = closed + self.par['barrierless_saddle_norbital']

            if shift_vec is None:
                with open('molpro/' + fname + '.inp', 'w') as outf:
                    outf.write(file.format(name=fname,
                                           natom=self.species.natom,
                                           geom=geom,
                                           nelectron=nelectron,
                                           symm=symm,
                                           spin=spin,
                                           charge=self.species.charge,
                                           state=self.par['barrierless_saddle_nstate'],
                                           closed=closed,
                                           occ=occ
                                           ))
            else:
                shift_vec = shift_vec / np.linalg.norm(shift_vec) * 0.5  # step of 0.5 A
                geom0 = ''
                scancoo = ''
                scanstart = ''
                shift = ''
                for i, at in enumerate(self.species.atom):
                    x, y, z = self.species.geom[i]
                    if i < self.species.natom - natom1:
                        geom0 += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
                    else:
                        scancoo += '{} s{} s{} s{}\n'.format(at, 3 * i, 3 * i + 1, 3 * i + 2)
                        scanstart += 's{} = {:.8f}\ns{}= {:.8f}\ns{}= {:.8f}\n'.\
                                     format(3 * i, x, 3 * i + 1, y, 3 * i + 2, z)
                        shift += 's{0} = s{0} + {1:.8f}\n'.format(3 * i, shift_vec[0])
                        shift += 's{0} = s{0} + {1:.8f}\n'.format(3 * i + 1, shift_vec[1])
                        shift += 's{0} = s{0} + {1:.8f}\n'.format(3 * i + 2, shift_vec[2])
                with open('molpro/' + fname + '.inp', 'w') as outf:
                    outf.write(file.format(name=fname,
                                           natom=self.species.natom,
                                           geom=geom0,
                                           scanstart=scanstart,
                                           scancoo=scancoo,
                                           shift=shift,
                                           nelectron=nelectron,
                                           symm=symm,
                                           spin=spin,
                                           charge=self.species.charge,
                                           state=self.par['barrierless_saddle_nstate'],
                                           closed=closed,
                                           occ=occ
                                           ))
        return 0

    def get_molpro_energy(self, key, name=''):
        """
        Verify if there is a molpro output file and if yes, read the energy
        key is the keyword for the energy we want to read
        returns 1, energy if successful
        returns 0, -1 if the energy was not there
        A non-object-oriented version is used in pes.py
        """
        fname = self.get_name(name)
        status = os.path.exists('molpro/' + fname + '.out')
        if status:
            with open('molpro/' + fname + '.out') as f:
                lines = f.readlines()

            for index, line in enumerate(reversed(lines)):
                if ('SETTING ' + key) in line:
                    return 1, float(line.split()[3])
        else:
            return 0, -1

    def create_molpro_submit(self, name=''):
        """
        write a pbs file for the molpro input file
        """
        fname = self.get_name(name)

        # open the template head and template
        molpro_head = pkg_resources.resource_filename('tpl', self.par['queuing'] + '.tpl')
        with open(molpro_head) as f:
            tpl_head = f.read()
        molpro_tpl = pkg_resources.resource_filename('tpl', self.par['queuing'] + '_molpro.tpl')
        with open(molpro_tpl) as f:
            tpl = f.read()
        # substitution
        with open('molpro/' + fname + '.' + self.par['queuing'], 'w') as f:
            if self.par['queuing'] == 'pbs':
                f.write((tpl_head + tpl).format(
                        name=fname,
                        ppn=self.par['single_point_ppn'],
                        queue_name=self.par['queue_name'],
                        dir='molpro',
                        command=self.par['single_point_command']))
            elif self.par['queuing'] == 'slurm':
                f.write((tpl_head + tpl).format(
                        name=fname,
                        ppn=self.par['single_point_ppn'],
                        queue_name=self.par['queue_name'],
                        dir='molpro',
                        command=self.par['single_point_command'],
                        slurm_feature=self.par['slurm_feature']))

        return 0

    def molpro_symm(self):
        if np.array_equal(self.species.atom, ['O']) and self.species.mult == 3:
            return 4
        if np.array_equal(self.species.atom, ['O', 'O']) \
                and self.species.mult == 3:
            return 4
        if np.array_equal(sorted(self.species.atom), sorted(['O', 'H'])) \
                and self.species.mult == 3:
            return 2
        return 1

    def get_name(self, name):
        if name != '':
            fname = name
        elif self.species.wellorts:
            fname = self.species.name
        else:
            fname = str(self.species.chemid)
        return fname
