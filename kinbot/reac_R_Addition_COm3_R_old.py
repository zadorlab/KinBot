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
import numpy as np
import copy
import time

from vector import *
from qc import *
from constants import *
from kinbot import *
from stationary_pt import *
from geom import *
from qc import *
import par



def do_R_Addition_COm3_R(species, instance, step, instance_name):
    """
    Carry out the reaction in reverse, which is an alpha scission.
    The number of steps of 0.1 Angstrom is controlled by a parameter.
    """

    if step > 0:
        if check_qc(instance_name) != 'normal': return step
    
    if step == 0:
        if os.path.isfile(instance_name + '.log'):
            err, freq = read_qc_freq(instance_name, par.natom)
            if err == 0 and freq[0] != 0.:
                step = par.scan_step + 1
                return step
            else:
                os.remove(instance_name + '.log')
        with open(par.tpldir + 'gauss_ts_scan_start.tpl') as f:
            lines = f.readlines()
    elif step < par.scan_step:
        with open(par.tpldir + 'gauss_ts_scan_mid.tpl') as f:
            lines = f.readlines()
    else:
        with open(par.tpldir + 'gauss_ts_search_end.tpl') as f:
            lines = f.readlines()    
            
    f = open(instance_name + '.com', 'w')
    for line in lines:
        line = line.replace('_chk_', instance_name)
        line = line.replace('_ppn_', str(par.ppn))
        line = line.replace('_charge_', str(par.charge))
        line = line.replace('_multiplicity_', str(par.mult))
        line = line.replace('_level_', str(par.level))
        if line.split():
            if line.split()[0] == '_geom_':
                for atom in range(par.natom):
                    f.write(par.atom[atom] + '  ')
                    [f.write(str(species.geom[atom][i]) + '  ') for i in range(3)]
                    f.write('\n')
            else:
                f.write(line)
        else:
            f.write(line)
    
    if step < par.scan_step:
        f.write(str(instance[0] + 1) + ' ' + str(instance[1] + 1) + ' += 0.1 F\n')
              
    
    f.write('\n\n')                                    
    f.close()
   
    step += submit_qc(instance_name, 0)

    return step


