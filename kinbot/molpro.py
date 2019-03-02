###################################################
##                                               ##
## This file is part of the KinBot code v2.0     ##
##                                               ##
## The contents are covered by the terms of the  ##
## BSD 3-clause license included in the LICENSE  ##
## file, found at the root.                      ##
##                                               ##
## Copyright 2018 National Technology &          ##
## Engineering Solutions of Sandia, LLC (NTESS). ##
## Under the terms of Contract DE-NA0003525 with ##
## NTESS, the U.S. Government retains certain    ##
## rights to this software.                      ##
##                                               ##
## Authors:                                      ##
##   Judit Zador                                 ##
##   Ruben Van de Vijver                         ##
##                                               ##
###################################################
import os
import pkg_resources

from kinbot import constants


class Molpro:
    """
    Class to write and read molpro file and to run molpro
    """
    def __init__(self, species, par):
        self.species = species
        self.par = par
        # self.qc = qc

    def create_molpro_input(self):
        """
        Create the input for molpro based on the template,
        which is either the one in the system, or provided
        by the user.
        """
        if self.par.par['single_point_template'] == '':
            tpl_file = pkg_resources.resource_filename('tpl', 'molpro.tpl')
        else:
            tpl_file = self.par.par['single_point_template']  
        with open(tpl_file) as f:
            file = f.read()

        fname = str(self.species.chemid)
        if self.species.wellorts:
            fname = self.species.name

        geom = ''
        nelectron = 0
        for i, at in enumerate(self.species.atom):
            x, y, z = self.species.geom[i]
            geom += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
            nelectron += constants.znumber[at]

        nelectron -= self.species.charge

        with open('molpro/' + fname + '.inp', 'w') as outf:
            outf.write(file.format(name=fname,
                                   natom=self.species.natom,
                                   geom=geom,
                                   nelectron=nelectron,
                                   spin=self.species.mult - 1,
                                   charge=self.species.charge
                                   ))


    def get_molpro_energy(self, key='MYENERGY'):
        """
        Verify if there is a molpro output file and if yes, read the energy
        key is the keyword for the energy we want to read
        returns 1, energy if successful
        returns 0, -1 if the energy was not there
        A non-object-oriented version is used in pes.py
        """
        fname = str(self.species.chemid)
        if self.species.wellorts:
            fname = self.species.name

        status = os.path.exists('molpro/' + fname + '.out')
        if status:
            with open('molpro/' + fname + '.out') as f:
                lines = f.readlines()

            for index, line in enumerate(reversed(lines)):
                if ('SETTING ' + key) in line:
                    return 1, float(line.split()[2])
        else:
            return 0, -1

    def create_molpro_submit(self):
        """
        write a pbs file for the molpro input file
        """

        fname = str(self.species.chemid)
        if self.species.wellorts:
            fname = self.species.name
        
        # open the template head and template
        molpro_head = pkg_resources.resource_filename('tpl', self.par.par['queuing'] + '.tpl')
        with open(molpro_head) as f:
            tpl_head = f.read()
        molpro_tpl = pkg_resources.resource_filename('tpl', self.par.par['queuing'] + '_molpro.tpl')
        with open(molpro_tpl) as f:
            tpl = f.read()
        # substitution
        with open('molpro/' + fname + '.' + self.par.par['queuing'], 'w' ) as f:
            f.write((tpl_head + tpl).format(name=fname, ppn=self.par.par['single_point_ppn'], queue_name=self.par.par['queue_name'], dir='molpro'))

        #command = ['qsub', 'run_molpro.pbs']
        #process = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        #out, err = process.communicate()
        #out = out.decode()
        #pid = out.split('\n')[0].split('.')[0]

        return 0



