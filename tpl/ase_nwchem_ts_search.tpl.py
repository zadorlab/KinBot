"""
Template to run ase search for a ts using NWChem
KinBot needs to pass to the template: 
1. A label for the calculation
2. The number of cores
3. The kwargs for NWChem
4. The atom vector
5. The geometry
6. The constraints for the optimization
    a. Fix: The coordinates to fix at their current value
    b. Change: The coordinates to change and fix at the new value
    c. Release: The coordinates to release (only for gaussian)
"""

import os, sys, re

from math import pi
import numpy as np

import ase
from ase import Atoms
from ase.calculators.nwchem import NWChem
from ase.db import connect


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
    #read the geometry from the output file
    outfile = '{label}.out'
    with open(outfile) as f:
        lines = f.readlines()
    for index, line in enumerate(reversed(lines)):
        if re.search('Output coordinates in angstroms', line) != None:
            for n in range(len(atom)):
                geom[n][0:3] = np.array(lines[-index+3+n].split()[3:6]).astype(float)
            break
    mol.positions = geom
    
    # do another single point calculation because the geometry has changed
    # and ase clears the energies and forces in this case.
    #del kwargs['task']
    #calc = NWChem(**kwargs)
    #mol.set_calculator(calc)
    #mol.get_potential_energy() # use the NWChem optimizer (task optimize)
    
    db = connect('kinbot.db')
    db.write(mol, name = label, data = {{'energy':e, 'status' : 'normal'}})
except RuntimeError, e: 
    db = connect('kinbot.db')
    db.write(mol, name = label, data = {{'status' : 'error'}})


f = open(label + '.out','a')
f.write('done\n')
f.close()
