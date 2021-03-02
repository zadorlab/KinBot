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
        self.uq_iter = 0


    def calc_factor(self, propertyType, species, uq_iter, runUQ):
        # runUQ = 0, treat all iterations as uq iter = 1
        # runUQ = 1, treat all interations as uqiter = 0 to n
        if self.uq_iter != uq_iter:  # new iteration
            with open('uqtk.data', 'a') as f:
                f.write('')  # new line
        if self.par['pes'] == 1 and runUQ == 0:
            uq_iter = 0
        elif self.par['pes'] == 1 and runUQ == 1:
            uq_iter = uq_iter

        if uq_iter == 0:
            if propertyType == 'freq' or propertyType == 'imagfreq':
                factor = 1
                normfactor = 0
                
            else:
                factor = 0
                normfactor = factor /self.wellUQ
        
            self.write_uqtk_header(species, propertyType)
            self.write_uqtk_data(propertyType, normfactor, species, uq_iter)

            return factor
        else:
            if propertyType == 'energy':
                factor = random.uniform(-self.wellUQ, self.wellUQ)
                normfactor = factor / self.wellUQ

            elif propertyType == 'freq':
                factor = np.exp(random.uniform(np.log(1./self.freqUQ), np.log(self.freqUQ)))
                normfactor = np.log(factor) / np.log(self.freqUQ)
        
            elif propertyType == 'imagfreq':
                factor = np.exp(random.uniform(np.log(1./self.imagfreqUQ), np.log(self.imagfreqUQ)))
                normfactor = np.log(factor) / np.log(self.imagfreqUQ)

            elif propertyType == 'barrier':
                factor = random.uniform(-self.barUQ, self.barUQ)
                normfactor = factor / self.barUQ

            self.write_uqtk_header(species, propertyType)
            self.write_uqtk_data(propertyType, normfactor, species, uq_iter)

        return factor


    def write_uqtk_header(self, species, propertyType):
        with open('uqtk.data', 'a') as fi:
            fi.write("{} {} ".format(species, propertyType))
        return 0


    def write_uqtk_data(self, propertyType, normfactor, species, uq_iter):
        with open('uqtk.data', 'a') as fi:
            fi.write("{} \n".format(normfactor))

        return 0


    def format_uqtk_data(self):
        with open('uqtk.data', 'r') as fi:
            data = {}
            for line in fi:
                line = line.split(' ')
                line.pop()
                data_point = line[0] + ' '  + line[1] 
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
                    print(reaction_items + '\n' + line[0])
                    if line[0] in reaction_items:
                        print('{} in rxnitems'.format(line[0]))
                        print(line[1])
                        if line[1] == 'freq' or line[1] == 'imagfreq':
                            uqfi = open('uqtk.data', 'a')
                            uqfi.write(line)
                            uqfi.close()

        f.close()

