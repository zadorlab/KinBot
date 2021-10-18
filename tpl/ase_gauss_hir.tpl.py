"""
Template to run ase to do one step of a hindered rotor calculation
KinBot needs to pass to the template: 
1. A label for the calculation
2. The number of cores
3. The kwargs for Gaussian
4. The atom vector
5. The geometry
6. The Gaussian command

"""

import os, sys, re

import numpy as np

import ase
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.db import connect


label = '{label}'
kwargs = {kwargs}

calc = Gaussian(**kwargs)
calc.command = '{qc_command} < PREFIX.com > PREFIX.log'

atom = {atom}
geom = {geom}


mol = Atoms(symbols = atom, positions = geom)
mol.set_calculator(calc)

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    #Positions (geom) updated in ase/ases/io/gaussian.py code
    dummy = {dummy}
    for d in dummy:
        #remove the dummy atom to write to the database
        mol.pop()
    db = connect('{working_dir}/kinbot.db')
    db.write(mol, name = label, data = {{'energy': e,'status' : 'normal'}})
except RuntimeError: 
    # in case of fail, try again with final geometry
    try:
        #read the geometry from the output file
        outfile = '{label}.log'
        with open(outfile) as f:
            lines = f.readlines()
        for index, line in enumerate(reversed(lines)):
            if re.search('Input orientation:', line) != None:
                for n in range(len(mol)):
                    geom[n][0:3] = np.array(lines[-index+4+n].split()[3:6]).astype(float)
                break
        dummy = {dummy}
        for i,d in enumerate(dummy):
            geom[-(i+1)][0:3] = d[0:3]
        mol.positions = geom
        e = mol.get_potential_energy() # use the Gaussian optimizer
        #read the geometry from the output file
        outfile = '{label}.log'
        with open(outfile) as f:
            lines = f.readlines()
        for index, line in enumerate(reversed(lines)):
            if re.search('Input orientation:', line) != None:
                for n in range(len(mol)):
                    geom[n][0:3] = np.array(lines[-index+4+n].split()[3:6]).astype(float)
                break
        for i,d in enumerate(dummy):
            geom[-(i+1)][0:3] = d[0:3]
        mol.positions = geom
        
        for d in dummy:
            #remove the dummy atom to write to the database
            mol.pop()
        db = connect('{working_dir}/kinbot.db')
        db.write(mol, name = label, data = {{'energy': e,'status' : 'normal'}})
    except:
        db = connect('{working_dir}/kinbot.db')
        db.write(mol, name = label, data = {{'status' : 'error'}})

f = open(label + '.log','a')
f.write('done\n')
f.close()
