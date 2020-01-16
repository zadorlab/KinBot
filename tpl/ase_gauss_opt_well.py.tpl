"""
Template to run ase to optimize a well using Gaussian
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
from ase.optimize import BFGS
from ase.db import connect

label = '{label}'
kwargs = {kwargs}

Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)

atom = {atom}
geom = {geom}


mol = Atoms(symbols = atom, positions = geom)
mol.set_calculator(calc)

method = '{{}}/{{}}'.format(kwargs['method'], kwargs['basis'])
natom = len([at for at in atom if at !='X']) #filter out the dummy atoms

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    """
    #read the geometry from the output file
    outfile = '{label}.log'
    #fc = open("doneCheck.log", 'a')
    #fc.write("Entering into gDone check\n")
    #gDone=outfile.closed()
    #while gDone == False:
    #    fc.write(gDone)
    #    time.sleep(1)
    #    gDone=outfile.closed()
    #fc.closed()        
    with open(outfile) as f:
        lines = f.readlines()
    for index, line in enumerate(reversed(lines)):
        if re.search('Input orientation:', line) != None:
            for n in range(len(mol)):
                geom[n][0:3] = np.array(lines[-index+4+n].split()[3:6]).astype(float)
            break
    mol.positions = geom
    
    """
    #read the frequencies
    with open('{label}.log') as f:
        lines = f.readlines()

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
    zpe = -1
    for line in reversed(lines):  
        if re.search('Zero-point correction=', line) != None:
            zpe = float(line.split()[2])
            break 
    dummy = {dummy}
    for d in dummy:
        #remove the dummy atom to write to the database
        mol.pop()
    db = connect('{working_dir}/kinbot.db')
    db.write(mol, name = label, data = {{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status' : 'normal'}})

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
        
        #read the frequencies
        with open('{label}.log') as f:
            lines = f.readlines()

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
        zpe = -1
        for line in reversed(lines):  
            if re.search('Zero-point correction=', line) != None:
                zpe = float(line.split()[2])
                break 
        for d in dummy:
            #remove the dummy atom to write to the database
            mol.pop()
        db = connect('{working_dir}/kinbot.db')
        db.write(mol, name = label, data = {{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status' : 'normal'}})
    except:
        db = connect('{working_dir}/kinbot.db')
        db.write(mol, name = label, data = {{'status' : 'error'}})

f = open(label + '.log','a')
f.write('done\n')
f.close()
