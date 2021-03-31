from __future__ import division
from __future__ import print_function
import os
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

    def calc_rotor_factor(self, propertyType, species, uq_iter, runUQ, name):
        string = []
        for i, char in enumerate(name):
            string.append(char)
        string.reverse()
        stop_parsing = 0
        name_array = []
        for char in string:
            if stop_parsing == 0:
                if char != " ":
                    name_array.append(char)
                elif char == " ":
                    stop_parsing = 1
            else:
                pass
        name_array.reverse()
        name_string = ''.join(name_array)

        if uq_iter == 0:
            factor = 1
            normfactor = 0
            self.write_uqtk_data(propertyType, normfactor, name_string, uq_iter)

            return factor
        else:
            factor = np.exp(random.uniform(np.log(1./self.hirUQ), np.log(self.hirUQ)))
            normfactor = np.log(factor) / np.log(self.hirUQ)
            if normfactor < -1 or normfactor > 1:
                e = "ERROR"
            else:
                e = ''
            self.write_uqtk_data(propertyType, normfactor, name_string, uq_iter)

        return factor

    def calc_factor(self, propertyType, species, uq_iter, runUQ):

        # runUQ = 0, treat all iterations as uq iter = 1
        #     runUQ = 0 during kinbot procedure for freq, and energies during PES runs
        #     PES will alter energies, frequencies at a later point.

        # runUQ = 1, treat all interations as uqiter = 0 to n
	    # runUQ = 1 for all energies, freq, HIR during kinbot only runs.
        #     runUQ = 1 for HIR during PES runs, because they are consistent between kinbot & pes mess files
        #     TODO: COULD runUQ = 1 for frequencies? during both pes and kinbot runs?

        if uq_iter == 0:
            if propertyType == 'freq' or propertyType == 'imagfreq' or propertyType == 'rotor' or propertyType == 'relax_factor' or propertyType == 'e_well' or propertyType == 's_well':
                factor = 1
                normfactor = 0
            else:
                factor = 0
                normfactor = 0

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

            elif propertyType == 'relax_factor':
                factor = np.exp(random.uniform(np.log(1./self.relaxFactorUQ), np.log(self.relaxFactorUQ)))
                normfactor = np.log(factor) / np.log(self.relaxFactorUQ)

            elif propertyType == 's_well':
                factor = np.exp(random.uniform(np.log(1./self.sWellUQ), np.log(self.sWellUQ)))
                normfactor = np.log(factor) / np.log(self.sWellUQ)

            elif propertyType == 'e_well':
                factor = np.exp(random.uniform(np.log(1./self.eWellUQ), np.log(self.eWellUQ)))
                normfactor = np.log(factor) / np.log(self.eWellUQ)

            self.write_uqtk_data(propertyType, normfactor, species, uq_iter)

        return factor

    def write_uqtk_data(self, propertyType, normfactor, species, uq_iter):
        file = "uq_" + str(species) + "_" + str(propertyType) + ".txt"
        with open(file, 'a') as fi:
            fi.write("{} | {}\n".format(uq_iter, normfactor))
            if uq_iter == self.par["uq_n"] + 1:
                fi.write("\n")

        return 0

    def format_uqtk_data(self):
        normalized_data = []
        for file in os.listdir("./"):
            if file.startswith("uq_"):
                parameters = []
                name=file
                parameters.append(name[3:-4:])
                read_file = open(file, 'r')
                for i, line in enumerate(read_file):
                    line_chars = []
                    for char in line:
                        line_chars.append(char)
                    line_chars = line_chars[2::]
                    truncated_line_chars = []
                    for i, char in enumerate(line_chars):
                        if line_chars[i] == ' ':
                            pass
                        elif line_chars[i] == '|':
                            pass
                        elif line_chars[i] == '\n':
                            pass
                        else:
                            truncated_line_chars.append(char)
                    final_string = ''.join(truncated_line_chars)
                    parameters.append(final_string)
                normalized_data.append(parameters)
        for x in normalized_data:
            print(x[0], len(x))
        normalized_data_cols = len(normalized_data[0])
        normalized_data_rows = len(normalized_data)
        print(normalized_data_cols, normalized_data_rows)
        fi = open("normalization.txt", 'w')
        row = 0 #i
        col = 0 #j
        while col < normalized_data_cols:
            norm_data = []
            for i, row in enumerate(normalized_data):
                norm_data.append(row[col])
            norm_data_string = " ".join(norm_data)
            fi.write(norm_data_string)
            fi.write("\n")
            col = col + 1
        fi.close()
        return 0

    def pes_freq_uqtk_data(self, parent, pes_reactions):
    # TO DO
    # ADD ALL normalization files together and delete columns that don't match pesviewer files?
        all_reaction = []
        for key in parent:
            with open(key + '/normalization.txt') as f:
                for line in f:
                    if line[0] in reaction_items:
                        if line[1] == 'freq' or line[1] == 'imagfreq':
                            uqfi = open('normalization.txt', 'a')
                            uqfi.write(line)
                            uqfi.close()

        f.close()
