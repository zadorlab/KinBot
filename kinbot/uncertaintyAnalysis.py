import random
import numpy as np

"""
Functions for the generation and postprocessing of uncertainty analysis
"""


class UQ:

    def __init__(self, par):
        self.uq_n = par['uq_n']
        self.limit = {'energy': par['well_uq'],
                      'barrier': par['barrier_uq'],
                      'freq': par['freq_uq'],
                      'imagfreq': par['imagfreq_uq'],
                      'epsilon': par['epsilon_uq'],
                      'sigma': par['sigma_uq'],
                      'enrelfact': par['enrelfact_uq'],
                      'enrelpow': par['enrelpow_uq'],
                      'pstsymm': par['pstsymm_uq'],
                      }

        self.additive = ['energy', 'barrier']

    def calc_factor(self, parameter, uq_iter):
        #if self.uq_iter != uq_iter:  # new iteration
        #    with open('uqtk.data', 'a') as f:
        #        f.write('')  # new line


        if uq_iter == 0:
            if parameter in self.additive:
                factor = 0
            else:
                factor = 1
            normfactor = 0
        
            return factor


        if parameter in self.additive:
            factor = random.uniform(-self.limit[parameter], self.limit[parameter])
            normfactor = factor / self.limit[parameter]
        else:
            factor = np.exp(random.uniform(np.log(1. / self.limit[parameter]), np.log(self.limit[parameter])))
            normfactor = np.log(factor) / np.log(self.limit[parameter])

        #self.write_uqtk_header(species, parameter)
        #self.write_uqtk_data(parameter, normfactor, species, uq_iter)

        return factor


    def write_uqtk_header(self, species, parameter):
        with open('uqtk.data', 'a') as fi:
            fi.write("{}{} ".format(species, parameter))
        return 0


    def write_uqtk_data(self, parameter, normfactor, uq_iter):
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
        
        names = []
        vals = []
        uq_counter = 0
        while uq_counter < self.uq_n:
            for item in data.items():
                names.append(item[0])
                vals.append(item[1])
                uq_counter += 1

        vals = np.array(vals, dtype=object)
        vals = vals.transpose()
        with open('fuqtk.data', 'w') as f:
            for item in vals:
                for i in item:
                    f.write(i + '\t')
                f.write('\n')

        return 0
