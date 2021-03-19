from __future__ import division
from __future__ import print_function
import random
import numpy as np

"""
Functions for the generation and postprocessing of uncertainty analysis
"""


class UQ:

    def __init__(self, par):
        self.par = par
        self.wellUQ = par['well_uq']
        self.barUQ = par['barrier_uq']
        self.freqUQ = par['freq_uq']
        self.imagfreqUQ = par['imagfreq_uq']
        self.hirUQ = par['hir_uq']
        self.relaxPowerUQ = par['energy_relaxation_power_uq']
        self.relaxFactorUQ = par['energy_relaxation_factor_uq']
        self.eWellUQ = par['epsilon_well_uq']
        self.sWellUQ = par['sigma_well_uq']
        self.uq_iter = 0

    def calc_factor(self, propertyType, species, uq_iter, runUQ):

        # runUQ = 0, treat all iterations as uq iter = 1
        #     runUQ = 0 during kinbot procedure for freq, and energies during PES runs
        #     PES will alter energies, frequencies at a later point.

        # runUQ = 1, treat all interations as uqiter = 0 to n
	#     runUQ = 1 for all energies, freq, HIR during kinbot only runs.
        #     runUQ = 1 for HIR during PES runs, because they are consistent between kinbot & pes mess files
        #     TODO: COULD runUQ = 1 for frequencies? during both pes and kinbot runs?

        if self.uq_iter != uq_iter:  # new iteration
            with open('uqtk.data', 'a') as f:
                f.write('')  # new line
        # DO WE NEED THIS IF/ELSE STATEMENT???
        if self.par['pes'] == 1 and runUQ == 0:
            uq_iter = 0
        elif self.par['pes'] == 1 and runUQ == 1:
            uq_iter = uq_iter


        if uq_iter == 0:
            if propertyType == 'freq' or propertyType == 'imagfreq' or propertyType == 'rotor' or propertyType == 'relax_factor' or propertyType == 'e_well' or propertyType == 's_well':
                factor = 1
                normfactor = 0
            else:
                factor = 0
                normfactor = 0
                # normfactor = factor / self.wellUQ

            self.write_uqtk_header(species, propertyType, uq_iter)
            self.write_uqtk_data(propertyType, normfactor, species, uq_iter)

            return factor
        else:
            if propertyType == 'energy':
                factor = random.uniform(-self.wellUQ, self.wellUQ)
                normfactor = factor / self.wellUQ

            elif propertyType == 'barrier':
                factor = random.uniform(-self.barUQ, self.barUQ)
                normfactor = factor / self.barUQ

            elif propertyType == 'relax_power':
                factor = random.uniform(-self.relaxPowerUQ, self.relaxPowerUQ)
                normfactor = factor / self.relaxPowerUQ
            
            elif propertyType == 'freq':
                factor = np.exp(random.uniform(np.log(1./self.freqUQ), np.log(self.freqUQ)))
                normfactor = np.log(factor) / np.log(self.freqUQ)

            elif propertyType == 'imagfreq':
                factor = np.exp(random.uniform(np.log(1./self.imagfreqUQ), np.log(self.imagfreqUQ)))
                normfactor = np.log(factor) / np.log(self.imagfreqUQ)

            elif propertyType == 'rotor':
                factor = np.exp(random.uniform(np.log(1./self.hirUQ), np.log(self.hirUQ)))
                normfactor = np.log(factor) / np.log(self.hirUQ)
                if normfactor < -1 or normfactor > 1:
                    e = "ERROR"
                else:
                    e = ''
                fi=open("rotors.txt", 'a')
                fi.write("{}\t{}\t{}\t{}\n".format(e, uq_iter, factor, normfactor))
                fi.close()

            elif propertyType == 'relax_factor':
                factor = np.exp(random.uniform(np.log(1./self.relaxFactorUQ), np.log(self.relaxFactorUQ)))
                normfactor = np.log(factor) / np.log(self.relaxFactorUQ)

            elif propertyType == 's_well':
                factor = np.exp(random.uniform(np.log(1./self.sWellUQ), np.log(self.sWellUQ)))
                normfactor = np.log(factor) / np.log(self.sWellUQ)

            elif propertyType == 'e_well':
                factor = np.exp(random.uniform(np.log(1./self.eWellUQ), np.log(self.eWellUQ)))
                normfactor = np.log(factor) / np.log(self.eWellUQ)

            self.write_uqtk_header(species, propertyType, uq_iter)
            self.write_uqtk_data(propertyType, normfactor, species, uq_iter)

        return factor

    def write_uqtk_header(self, species, propertyType, uq_iter):
        with open('uqtk.data', 'a') as fi:
            fi.write("{} {} {}\n".format(species, propertyType, uq_iter))
        return 0

    def write_uqtk_data(self, propertyType, normfactor, species, uq_iter):
        with open('uqtk.data', 'a') as fi:
            fi.write("{}\n".format(normfactor))

        return 0

    def format_uqtk_data(self):
        with open('uqtk.data', 'r') as fi:
            data = {}
            for line in fi:
                line = line.split(' ')
                line.pop()
                data_point = line[0] + ' ' + line[1]
                if data_point not in data:
                    data[data_point] = [line[2]]
                else:
                    data[data_point].append(line[2])

        uq_n = self.par['uq_n']
        names = []
        vals = []
        uq_counter = 0
        while uq_counter < uq_n:
            for item in data.items():
                names.append(item[0])
                vals.append(item[1])
                uq_counter += 1

        vals = np.array(vals)
        vals = vals.transpose()
        with open('fuqtk.data', 'w') as f:
            for item in vals:
                for i in item:
                    f.write(i + '\t')
                f.write('\n')

        return 0

    def pes_freq_uqtk_data(self, parent, reaction_items):
        for key in parent:
            with open(key + '/uqtk.data') as f:
                for line in f:
                    if line[0] in reaction_items:
                        if line[1] == 'freq' or line[1] == 'imagfreq':
                            uqfi = open('uqtk.data', 'a')
                            uqfi.write(line)
                            uqfi.close()

        f.close()
