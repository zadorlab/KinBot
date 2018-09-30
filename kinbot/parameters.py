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
import os
import json
import logging

import imp
import numpy as np

from cheminfo import *

from ase.db import connect

class Parameters:
    """
    This class initiates all parameters to their defaults and reads in the 
    user-defined variables, which overwrite the defaults
    """
    def __init__(self,file):
        """
        Initialize all the variable and read the file which is the user input
        file
        """
        #user input file
        self.input_file = file
        
        self.par = {
            # GENERAL INFO
            #title of the current calculations
            'title' : '',
            
            # WHICH STEPS TO TAKE
            #Do a reaction search
            'reaction_search' : 1,
            #Which reaction families to include in the search
            'families' : ['all'],
            #Threshold above which barriers are deemed unimportant
            'barrier_threshold' : 0.,
            #Number of 0.1 Angstrom steps in bond scans
            'scan_step' : 30,
            #Do a full PES scan instead of one well
            'pes' : 0,
            #Perform high level optimization and freq calculation
            'high_level' : 0,
            #Do a conformational search
            'conformer_search' : 1,
            #Do a hindered rotor scan
            'rotor_scan' : 0,
            #Number of points along the rotor scan
            'nrotation' : 12,
            #Do master equation calculations
            'me' : 0,
            
            # QUANTUM CHEMISTRY INFO
            #Which quantum chemistry code to use
            'qc' : 'gauss', # or nwchem
            #nwchem-specific parameter
            'methodclass' : 'dft', # or scf or mp2  
            #Command for gaussian
            'gaussian_command' : 'g09',
            #Command for NWChem
            'nwchem_command' : 'nwchem',
            #Quantum chemistry method to use
            'method' : 'b3lyp',
            #Basis set to use
            'basis' : '6-31G',
            #Quantum chemistry method to use for high-level
            'high_level_method' : 'M062X',
            #Basis set to use for high-level
            'high_level_basis' : '6-311++G(d,p)',
            
            # INPUT SPECIES INFOR
            #SMILES of the species
            'smiles' : '',
            #geometry of the species
            'structure' : [],
            #Charge of the species
            'charge' : 0,
            #Multiplicity of the species
            'mult' : 0,
            #List of the elements
            'atom' : [],
            #Number of atoms
            'natom' : 0,
            
            # COMPUTATIONAL ENVIRONEMNT
            #Directory with the templates
            'tpldir' : os.path.expanduser('~/KinBot/tpl/'),
            #Which queuing system to use
            'queuing' : 'pbs', # or slurm
            #Scratch directory
            'scratch' : '/scratch/jzador',
            #User name
            'username' : 'jzador',
            #Name of the queue
            'queue_name' : 'medium',
            #E.g. the type of node or anything that comes with -C in SLURM
            'slurm_feature' : 'knl',
            #Number of cores to run the qc jobs on
            'ppn' : 1,
            #This many spaces can be used for numbering files, e.g., in ga
            'zf' : 4,

            # MASTER EQUATION
            #Which ME code to use:
            'me_code' : 'mess', # or mesmer
            #MESS specific keywords
            'mess_command' : 'mess',
            'TemperatureList' : [500+100*i for i in range(16)],
            'PressureList' : [760],
            'EnergyStepOverTemperature' : .2,
            'ExcessEnergyOverTemperature' : 30,
            'ModelEnergyLimit' : 400,
            'CalculationMethod' : 'direct', # or low-eigenvalue
            'ChemicalEigenvalueMax' : 0.2,
            'EnergyRelaxationFactor' : 200,
            'EnergyRelaxationPower' : .85,
            'EnergyRelaxationExponentCutoff' : 15,
            'Epsilons' : [7.08,310.387],
            'Sigmas' : [2.576,6.000],
            'Masses' : [4.0,87.0],
            #MESMER specific keywords
            'mesmer_command' : 'mesmer',
        }
        #Read the user input and overwrite the user-defined parameters
        self.read_user_input()

    def read_user_input(self):
        """
        Read the user input file and overwrite the default values
        """
        print self.input_file
        with open(self.input_file) as json_file:
            print json_file
            user_data = json.load(json_file)
            for key in user_data:
                if key in self.par:
                    self.par[key] = user_data[key]
                else:
                    logging.error('KinBot does not recognize option {} with value {}'.format(key,user_data[key]))

    def print_parameters(self):
        """
        Make a string out of the parameters
        """
        s = ''
        for key in self.par:
            s += '{}\t{}\n'.format(key,self.par[key])
        return s

def read_input(inputfile):
    """
    Read the main input file.
    Whenever possible, parameters have defaults.
    """

    global par
        
    f = open(os.path.expanduser('~/KinBot/kinbot/default_par.dat'))
    par = imp.load_source('par', '', f)
    print par.title
    f.close()

    f = open(inputfile)    
    par = imp.load_source('par', '', f)
    print par
    print par.title
    f.close()

    global level

    if par.basis == 'none':
        level = par.method # for composite methods in Gaussian
    else:
        level = par.method + '/' + par.basis # for Gaussian

    if not hasattr(par,'structure'):
        #generate the structure from the smiles
        obmol, par.structure = generate_3d_structure(par.smiles)
        if not hasattr(par,'natom'):
            par.natom = len(obmol.atoms)

    structure = np.reshape(par.structure, (par.natom,4))
    
    global atom

    atom = structure[:,0]
    geom = structure[:,1:4].astype(float)
    
    #ase database object
    global db
    db = connect('kinbot.db')

    return geom



def main():
    """
    This module reads the input file.
    """



if __name__ == "__main__":
    main()
    
