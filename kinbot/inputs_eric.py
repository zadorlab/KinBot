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
import os,sys,shutil,time
#import logging
#import datetime
sys.dont_write_bytecode = True

import imp

from thread_kinbot import *
from cheminfo import *
#try to import pybel to generate the structure from the smiles using OpenBabel
try:
    import pybel
except ImportError:
    pass


def main(workdir):
    """
    Main method for generating many kinbot runs to get a large dataset of initial ts structures
    """
    dir = os.path.expanduser(workdir)
    
    #read the .dat file
    f = open('{}smi.dat'.format(dir))
    par = imp.load_source('par', '', f)
    
    #make a sdf file for visualization
    output = pybel.Outputfile("sdf", dir + "species.sdf",overwrite=True)
    for name in par.smiles:
        smi = par.smiles[name]
        obmol = pybel.readstring("smi",smi)
        output.write(obmol)
    output.close()
    
    #list with the jobs that need to be done
    jobs = []
    
    #iterate the input files
    for name in par.smiles:
        #name = input_file.replace('.inp','') #name of the calculation
        test_dir = dir + name #location where the calculations will be done
        if not os.path.exists(test_dir):
            os.mkdir(test_dir)
                    
        #copy the input file to the working directory
        write_input_file(par,name,par.smiles[name],test_dir + '/input.inp')
        job = workdir + name + '/'
        jobs.append(job)
    
    run_threads(jobs, 'eric', max_running = 3)

def write_input_file(par,name,smi,file_name):
    f = open(file_name,'w+')
    f.write('title = \'%s\'\n\n'%name)
    f.write('method = \'%s\'\n'%par.method)
    f.write('basis = \'%s\'\n'%par.basis)
    f.write('qc = \'%s\'\n'%par.qc)
    f.write('conformer_search = %i\n'%par.conformer_search)
    f.write('reaction_search = %i\n'%par.reaction_search)
    f.write('barrier_threshold = %.1f\n'%par.barrier_threshold)
    f.write('ga = %i\n'%par.ga)
    f.write('ngen = %i\n'%par.ngen)
    f.write('ppn = %i\n\n'%par.ppn)
    f.write('eric = 1\n\n')
    
    obmol = pybel.readstring('smi',smi)
    obmol.OBMol.AddHydrogens()
    
    charge = 0
    f.write('charge = %i\n'%charge)
    mult = obmol.spin
    f.write('mult = %i\n'%mult)
    natom = len(obmol.atoms)
    f.write('natom = %i\n'%natom)
    f.write('smiles = \'%s\'\n'%smi)
    

if __name__ == "__main__":
    workdir = '~/eric/'
    main(workdir)