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
sys.path.insert(0,os.path.expanduser('~/ml-kinbot/code/kinbot'))

import imp

from thread_kinbot import *
#try to import pybel to generate the structure from the smiles using OpenBabel
try:
    import pybel
except ImportError:
    pass


def main(workdir):
    """
    Main method for the testing of kinbot
    """
    wipe_jobs = 1 #this boolean tells if previous jobs should be deleted (default should be 1)
    
    #read the .dat file
    f = open('test_kinbot.dat')
    par = imp.load_source('par', '', f)

    #list with the jobs that need to be done
    jobs = []
    
    #iterate the input files
    for name in par.jobs:
        #name = input_file.replace('.inp','') #name of the calculation
        test_dir = os.path.expanduser(workdir) + name #location where the calculations will be done
        if not os.path.exists(test_dir):
            os.mkdir(test_dir)
        else:
            if wipe_jobs: #check if the previous calculations need to be deleted
                for f in os.listdir(test_dir):
                    if os.path.isfile(test_dir + '/' + f):
                        os.remove(test_dir + '/' + f) 
                    
        #copy the input file to the working directory
        write_input_file(par,name,test_dir + '/input.inp')
        #shutil.copy(os.path.expanduser(workdir + 'input_files/') + input_file,test_dir + '/input.inp')
        job = workdir + name + '/'
        jobs.append(job)
    
    #run_threads(jobs, 'test_kinbot')

def write_input_file(par,name,file_name):
    f = open(file_name,'w+')
    f.write('title = \'%s\'\n\n'%name)
    f.write('method = \'%s\'\n'%par.method)
    f.write('basis = \'%s\'\n'%par.basis)
    f.write('qc = \'%s\'\n'%par.qc)
    f.write('conformer_search = %i\n'%par.conformer_search)
    f.write('reaction_search = %i\n'%par.reaction_search)
    f.write('barrier_threshold = %.1f\n'%par.barrier_threshold)
    f.write('families = [%s]\n'%','.join(["'%s'"%fi for fi in par.jobs[name][1]]))
    f.write('ga = %i\n'%par.ga)
    f.write('ngen = %i\n'%par.ngen)
    f.write('ppn = %i\n\n'%par.ppn)
    
    smi = par.jobs[name][0]
    obmol = pybel.readstring('smi',smi)
    obmol.OBMol.AddHydrogens()
    
    charge = 0
    f.write('charge = %i\n'%charge)
    mult = obmol.spin
    f.write('mult = %i\n'%mult)
    natom = len(obmol.atoms)
    f.write('natom = %i\n'%natom)
    f.write('smiles = \'%s\'\n'%smi)
    
    if name in par.structures: 
        f.write('structure = %s\n\n'%par.structures[name])
    
    
    

if __name__ == "__main__":
    workdir = '~/ml-kinbot/tests/kinbot/'
    main(workdir)