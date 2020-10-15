"""
Template to run ase to optimize a ts using Gaussian
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

Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)

atom = {atom}
geom = {geom}


mol = Atoms(symbols = atom, positions = geom)
mol.set_calculator(calc)
try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    #Positions (geom) updated in ase/ases/io/gaussian.py code    
    #read the frequencies
    natom = len(mol)
    if natom == 1:
        freq = np.array([])
    elif natom == 2:
        freq = np.array([0.])
    else:
        freq = np.zeros(3 * natom - 6) 
    freq_all = np.zeros(3 * natom) 
    with open('{label}.log') as f:
        lines = f.readlines()
    nfreq = 0
    for line in lines:
        if re.search('Frequencies', line) != None:
            if natom == 2:
                freq[0] = np.array(line.split()[2]).astype(float)
                break
            else:
                freq[nfreq:nfreq+3] = np.array(line.split()[2:5]).astype(float)
                nfreq = nfreq + 3

    #read the zpe
    for line in reversed(lines):  
        if re.search('Zero-point correction=', line) != None:
            zpe = float(line.split()[2])
            break 

    db = connect('{working_dir}/kinbot.db')
    db.write(mol,name = label,data={{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status' : 'normal'}})
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
        mol.positions = geom
        #read the frequencies
        natom = len(mol)
        if natom == 1:
            freq = np.array([])
        elif natom == 2:
            freq = np.array([0.])
        else:
            freq = np.zeros(3 * natom - 6) 
        freq_all = np.zeros(3 * natom) 
        with open('{label}.log') as f:
            lines = f.readlines()
        nfreq = 0
        for line in lines:
            if re.search('Frequencies', line) != None:
                if natom == 2:
                    freq[0] = np.array(line.split()[2]).astype(float)
                    break
                else:
                    freq[nfreq:nfreq+3] = np.array(line.split()[2:5]).astype(float)
                    nfreq = nfreq + 3

        #read the zpe
        for line in reversed(lines):  
            if re.search('Zero-point correction=', line) != None:
                zpe = float(line.split()[2])
                break 

        db = connect('{working_dir}/kinbot.db')
        db.write(mol,name = label,data={{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status' : 'normal'}})
    except:
        db = connect('{working_dir}/kinbot.db')
        db.write(mol, name = label, data = {{'status' : 'error'}})

f = open(label + '.log','a')
f.write('done\n')
f.close()
