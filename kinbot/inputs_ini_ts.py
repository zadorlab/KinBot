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
import os,sys,shutil,time,json
import imp

import thread_kinbot 
import cheminfo
from parameters import Parameters

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
    
    #read the .json file
    f = '{}smi.json'.format(dir)
    par = Parameters(f)
    
    #make a sdf file for visualization
    output = pybel.Outputfile("sdf", dir + "species.sdf",overwrite=True)
    for name in sorted(par.par['smiles'].keys()):
        smi = par.par['smiles'][name]
        obmol = pybel.readstring("smi",smi)
        output.write(obmol)
    output.close()
    
    #list with the jobs that need to be done
    # keys: directories of the jobs
    # values: names of the json input files
    jobs = {}
    
    #iterate the input files
    for name in sorted(par.par['smiles'].keys()):
        #name = input_file.replace('.inp','') #name of the calculation
        test_dir = dir + name #location where the calculations will be done
        if not os.path.exists(test_dir):
            os.mkdir(test_dir)
        
        #input file
        inpfile = name + '.json'
        
        #copy the input file to the working directory
        write_input_file(par,name,par.par['smiles'][name],test_dir + '/' + inpfile)
        job = workdir + name + '/'
        jobs[job] = inpfile
    
    thread_kinbot.run_threads(jobs, 'ini_ts', max_running = 10)

def write_input_file(par,name,smi,file_name):
    #make a new parameters instance and overwrite some keys
    par2 = Parameters(par.input_file)
    
    obmol = pybel.readstring('smi',smi)
    obmol.OBMol.AddHydrogens()
    
    par2.par['title'] = name
    par2.par['charge'] = 0
    par2.par['mult'] = obmol.spin
    par2.par['smiles'] = smi
    
    with open(file_name,'w') as outfile:
        json.dump(par2.par,outfile,indent = 4, sort_keys = True)

if __name__ == "__main__":
    workdir = '~/ini_ts/'
    main(workdir)
