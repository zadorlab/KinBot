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

import constants


class Molpro:
    """
    Class to write and read molpro file and to run molpro
    """
    def __init__(self, species, par, qc):
        self.species = species
        self.par = par
        self.qc = qc

    def create_molpro_input(self):
        """
        Create the input for molden
        """
        tpl_file = pkg_resources.resource_filename('tpl', 'molpro.tpl')
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
                                   electron=nelectron,
                                   spin=self.species.mult - 1,
                                   charge=self.species.charge
                                   ))

    def get_molpro_energy(self):
        """
        Verify if there is a molpro output file and if yes, read the energy
        """
        fname = str(self.species.chemid)
        if self.species.wellorts:
            fname = self.species.name

        status = os.path.exists('molpro/' + fname + '.out')
        if status:
            with open('molpro/' + fname + '.out') as f:
                lines = f.readlines()

            for index, line in enumerate(reversed(lines)):
                if 'SETTING MYENA' in line:
                    return 1, float(line.split()[2])
        else:
            return 0, -1

    def run(self):
        """
        TODO
        """
        pass
