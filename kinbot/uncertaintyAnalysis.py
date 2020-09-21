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
        self.uq_iter = 0


    def calc_factor(self, propertyType, species, uq_iter):
        if self.uq_iter != uq_iter:  # new iteration
            with open('uqtk.dat', 'a') as f:
                f.write('')  # new line

        if uq_iter == 0:
            if propertyType == 'freq' or propertyType == 'imagfreq':
                factor = 1
            else:
                factor = 0
            self.write_uqtk_header(species, propertyType)
            return factor
        
        if propertyType == 'energy':
            factor = random.uniform(-self.wellUQ, self.wellUQ)
            normfactor = factor / self.wellUQ

        elif propertyType == 'freq':
            # Amanda, this is the right way to sample from a log normal distribution
            factor = np.exp(random.uniform(np.log(1./self.freqUQ), np.log(self.freqUQ)))
            normfactor = np.log(factor) / np.log(self.freqUQ)

        elif propertyType == 'imagFreq':
            factor = np.exp(random.uniform(np.log(1./self.imagfreqUQ), np.log(self.imagfreqUQ)))
            normfactor = np.log(factor) / np.log(self.imagfreqUQ)

        elif propertyType == 'barrier':
            factor = random.uniform(-self.barUQ, self.barUQ)
            normfactor = factor / self.barUQ

        self.write_uqtk_data(normfactor)

        return factor


    def write_uqtk_header(self, species, propertyType):
        with open('uqtk.dat', 'a') as fi:
            fi.write("{}/{}".format(species, propertyType.rstrip('\n')))
        return 0


    def write_uqtk_data(self, normfactor)
        with open('uqtk.dat', 'a') as fi:
            fi.write("{}".format(normfactor.rstrip('\n')))
        return 0
