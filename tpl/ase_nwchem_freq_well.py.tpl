"""
Template to run ase to calculate the frequencies of a well using NWChem
KinBot needs to pass to the template: 
1. A label for the calculation
2. The number of cores
3. The kwargs for NWChem
4. The atom vector
5. The geometry

"""

import os, sys, re

import numpy as np

import ase
from ase import Atoms
from ase.calculators.nwchem import NWChem
from ase.vibrations import Vibrations
from ase.db import connect
from ase.io import read


label = '{label}'
kwargs = {kwargs}

NWChem.command = 'mpirun -np {ppn} -path /usr/local/bin nwchem PREFIX.nw > PREFIX.out'
calc = NWChem(**kwargs)

atom = {atom}
geom = {geom}


mol = Atoms(symbols = atom, positions = geom)
mol.set_calculator(calc)
try:
    e = mol.get_potential_energy() # use the NWChem optimizer (task optimize)
    
    #read the frequencies
    natom = len(mol)
    if natom == 1:
        freq = np.array([])
    elif natom == 2:
        freq = np.array([0.])
    else:
        freq = np.zeros(3 * natom - 6) 
    freq_all = np.zeros(3 * natom) 
    with open('{label}.out') as f:
        lines = f.readlines()
    nfreq = 0
    for line in lines:
        if re.search('P.Frequency', line) != None:
            nval = len(line.split()) - 1
            freq_all[nfreq:nval+nfreq] = np.array(line.split()[1:nval+1]).astype(float)
            nfreq += nval
    i = 0
    for fr in range(nfreq):
        if abs(freq_all[fr]) > 5:
            freq[i] = freq_all[fr]
            i += 1 

    #read the zpe
    for line in reversed(lines):  
        if re.search('Zero-Point correction to Energy', line) != None:
            zpe = float(line.split()[8])
            break  

    db = connect('kinbot.db')
    db.write(mol,name = label,data={{'energy': e, 'frequencies': np.asarray(freq), 'zpe':zpe, 'status' : 'normal'}})
except RuntimeError, e: 
    db = connect('kinbot.db')
    db.write(mol, name = label, data = {{'status' : 'error'}})



# add the finished stamp
f = open(label + '.out','a')
f.write('done\n')
f.close()

