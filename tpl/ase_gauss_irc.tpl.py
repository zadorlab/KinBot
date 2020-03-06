"""
Template to run ircs with ase using Gaussian
KinBot needs to pass to the template: 
1. A label for the calculation
2. The number of cores
3. The kwargs for Gaussian
4. The atom vector
5. The geometry
6. Tha Gaussian command
"""

import os, sys, re

from math import pi
import numpy as np

import ase
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.optimize import BFGS
from ase.db import connect
from ase.constraints import FixInternals


label = '{label}'
kwargs = {kwargs}

Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)

atom = {atom}
geom = {geom}

mol = Atoms(symbols = atom, positions = geom)
mol.set_calculator(calc)

success = 1

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    #Positions (geom) updated in ase/ases/io/gaussian.py code    
    geom = mol.positions
    db = connect('{working_dir}/kinbot.db')
    db.write(mol, name = label, data = {{'energy': e,'status' : 'normal'}})
except:
    #read the geometry from the output file
    outfile = '{label}.log'
    with open(outfile) as f:
        lines = f.readlines()
    geom = np.zeros((len(mol),3))
    new_geom = 0
    for index, line in enumerate(reversed(lines)):
        if re.search('Input orientation:', line) != None:
            for n in range(len(mol)):
                geom[n][0:3] = np.array(lines[-index+4+n].split()[3:6]).astype(float)
            new_geom = 1
            break
    if new_geom:
        mol.positions = geom
        db = connect('{working_dir}/kinbot.db')
        db.write(mol, name = label, data = {{'status' : 'normal'}}) #although there is an error, continue from the final geometry
    else:
        success = 0
        db = connect('{working_dir}/kinbot.db')
        db.write(mol, name = label, data = {{'status' : 'error'}})

if success:
    # start the product optimization
    pr_kwargs = {prod_kwargs}
    label = '{label}_prod'
    
    calc = Gaussian(**pr_kwargs)
    mol = Atoms(symbols = atom, positions = geom)
    mol.set_calculator(calc)
    try:
        e = mol.get_potential_energy() # use the Gaussian optimizer
        #read the geometry from the output file
        outfile = '{label}_prod.log'
        with open(outfile) as f:
            lines = f.readlines()
        geom = np.zeros((len(mol),3))
        for index, line in enumerate(reversed(lines)):
            if re.search('Input orientation:', line) != None:
                for n in range(len(mol)):
                    geom[n][0:3] = np.array(lines[-index+4+n].split()[3:6]).astype(float)
                break
        mol.positions = geom
        db = connect('{working_dir}/kinbot.db')
        db.write(mol, name = label, data = {{'energy': e,'status' : 'normal'}})
    except RuntimeError: 
        #read the geometry from the output file
        outfile = '{label}_prod.log'
        with open(outfile) as f:
            lines = f.readlines()
        geom = np.zeros((len(mol),3))
        new_geom = 0
        for index, line in enumerate(reversed(lines)):
            if re.search('Input orientation:', line) != None:
                for n in range(len(mol)):
                    geom[n][0:3] = np.array(lines[-index+4+n].split()[3:6]).astype(float)
                new_geom = 1
                break
        if new_geom:
            mol.positions = geom
            db = connect('{working_dir}/kinbot.db')
            db.write(mol, name = label, data = {{'status' : 'normal'}}) #although there is an error, continue from the final geometry
        else:
            db = connect('{working_dir}/kinbot.db')
            db.write(mol, name = label, data = {{'status' : 'error'}})


f = open(label + '.log','a')
f.write('done\n')
f.close()
