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
        self.imagFreqUQ = par['imagfreq_uq']
        print("UQ init: {} {} {} {}".format(self.wellUQ, self.barUQ, self.freqUQ, self.imagFreqUQ))

    def calc_factor(self, propertyType, species, uq_iter):
        if uq_iter == 0:
            if propertyType == 'freq' or propertyType == 'imagfreq':
                factor = 1
            else:
                factor = 0
        
            self.write_file(propertyType, factor, species, uq_iter)

            return factor
        
        if propertyType == 'energy':
            factor = random.uniform(-self.wellUQ, self.wellUQ)

        elif propertyType == 'freq':
            factorMin = 1.0 / self.freqUQ
            factor = random.uniform(factorMin, self.freqUQ)

        elif propertyType == 'imagFreq':
            factorMin = 1.0 / self.imagFreqUQ
            factor = random.uniform(factorMin, self.imagFreqUQ)

        elif propertyType == 'barrier':
            factor = random.uniform(-self.barUQ, self.barUQ)

        self.write_file(propertyType, factor, species, uq_iter)

        return factor

    def write_file(self, propertyType, propertyValue, species, uq_iter):
        uqfile = 'uq_' + propertyType + '.txt'
        with open(uqfile, 'a') as fi:
            fi.write("{}, {}, {}\n".format(uq_iter, species, propertyValue))
        fi.close()
        return 0
