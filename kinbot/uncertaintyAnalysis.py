from __future__ import division
from __future__ import print_function
import random
import numpy as np

"""
Functions for the generation and postprocessing of uncertainty analysis
"""


class UQ:

    # random.seed in pes.py, if pes off then inside kb.py, make parameter
    def calc_energyUQ(self, factor):
        factor = random.uniform(-factor, factor)

        return factor

    def calc_freqUQ(self, factor):
        factorMax = factor
        factorMin = 1.0 / factor
        factor = random.uniform(factorMin, factorMax)

        return factor

    def norm_energy(self, energyList, speciesType, names, n):
        fi = 'normalizedData.log'
        fio = open(fi, 'a')
        fio.write("\n{} UQ Normalization".format(speciesType))
        names = names
        print(names)
        data = energyList
        if len(names) != 0 and len(energyList) == 0:
            elements = n * len(names)
            energyList = np.zeros(elements)
        if len(energyList) == 0:
            fio.write("\nNo {} data".format(speciesType))
            return 0
        
        # number of entries > n means more than one structure
        if len(energyList) > n:
            nSpecies = int(len(energyList) / n)
            nCount = 0
            totalSpecies = len(energyList)
            eArray = []
            while nCount < nSpecies:
                eArray.append(energyList[nCount:totalSpecies:nSpecies])
                nCount += 1
            for i, species in enumerate(eArray):
                norm = []
                npData = species
                originalEnergy = npData[0]
                low = np.min(npData)
                high = np.max(npData)
                for e in npData:
                    norm_e = 2 * ((e - low) / (high - low)) - 1
                    norm.append(norm_e)
                npNormVals = np.array(norm)
                if len(npNormVals) == n:
                    i = 0
                    print(names)
                    for v, val in enumerate(npNormVals):
                        norm_energy = ('{:.4}'.format(float(npNormVals[v])))
                        en = ('{:.4}'.format(float(npData[v])))
                        x = v % 5
                        rxn = (names[x])
                        fio.write("\n{}\t{}\t{}".format(rxn, en, norm_energy))
                        
        else:
            # stats = []  # originalEnergy, min, max, avg, stDev
            normVals = []
            npData = np.array(data)
            if len(npData) == 0:
                norm = []
                npNormVals = []
            for i, val_i in enumerate(npData):
                norm = []
                # stats = []
                originalEnergy = npData[0]
                # stats.append(originalEnergy)
                low = np.min(npData)
                # stats.append(low)
                high = np.max(npData)
                # stats.append(high)
                # avg = np.mean(npData)
                # stats.append(avg)
                # stDev = np.std(npData)
                # stats.append(stDev)
                for e in data:
                    norm_e = 2 * ((e - low) / (high - low)) - 1
                    norm.append(norm_e)
                npNormVals = np.array(norm)
            print("data: {}".format(npData))
            print("names: {}".format(names))
            if len(npNormVals) == n:
                for k, val_k in enumerate(npNormVals):
                    norm_energy = ('{:.4}'.format(str(float(npNormVals[k]))))
                    energy = ('{:.4}'.format(float(str(npData[k])[1: -1])))
                    fio.write("\n{}\t{}\t{}".format(names[k], energy, norm_energy))
        return 0

    def norm_imagfreq(self, imagfreqList, n):

        return norm_imagfreq_vals

    def norm_posFreq(self, posFreqList, n):

        return norm_posFreq_vals
