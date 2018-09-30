"""
Template to run ase to calculate the frequencies of a well using Gaussian
KinBot needs to pass to the template: 
1. A label for the calculation
2. The number of cores
3. The kwargs for Gaussian
4. The atom vector
5. The geometry

"""

import os, sys, re

import numpy as np

import ase
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.vibrations import Vibrations
from ase.db import connect
from ase.io import read


label = '{label}'
kwargs = {kwargs}

Gaussian.command = 'g09 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)

atom = {atom}
geom = {geom}


mol = Atoms(symbols = atom, positions = geom)
mol.set_calculator(calc)
try:
    e = mol.get_potential_energy() # use the Gaussian optimizer (task optimize)
    
    #read the frequencies
    with open('{label}.log') as f:
        lines = f.readlines()

    natom = len([at for at in atom if at !='X']) #filter out the dummy atoms
    if natom == 1:
        freq = []
    else:
        freq = []
        for line in lines:
            if re.search('Frequencies', line) != None:
                if natom == 2:
                    freq.append(np.array(line.split()[2]).astype(float))
                    break
                else:
                    f = np.array(line.split()[2:5]).astype(float)
                    freq.extend(f)

    #read the zpe
    for line in reversed(lines):  
        if re.search('Zero-point correction=', line) != None:
            zpe = float(line.split()[2])
            break 

    dummy = {dummy}
    for d in dummy:
        #remove the dummy atom to write to the database
        mol.pop()
        
    db = connect('kinbot.db')
    db.write(mol,name = label,data={{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status' : 'normal'}})
except RuntimeError, e: 
    db = connect('kinbot.db')
    db.write(mol, name = label, data = {{'status' : 'error'}})

# add the finished stamp
f = open(label + '.log','a')
f.write('done\n')
f.close()

