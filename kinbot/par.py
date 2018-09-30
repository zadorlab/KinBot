import os
#import random

import imp
import numpy as np

from cheminfo import *

from ase.db import connect

def read_input(inputfile):
    """
    Read the main input file.
    Whenever possible, parameters have defaults.
    """
    

    global par
        
    f = open(os.path.expanduser('~/ml-kinbot/code/kinbot/default_par.dat'))
    par = imp.load_source('par', '', f)
    f.close()

    f = open(inputfile)    
    par = imp.load_source('par', '', f)
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
    
