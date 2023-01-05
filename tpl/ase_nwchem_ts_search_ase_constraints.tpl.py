"""
Template to run ase to optimize a well using NWChem
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

from math import pi

from ase import Atoms
from ase.calculators.nwchem import NWChem
from ase.optimize import BFGS
from ase.db import connect
from kinbot.ase_modules.constraints import FixInternals

label = '{label}'
kwargs = {kwargs}

NWChem.command = 'mpirun -np {ppn} -path /usr/local/bin nwchem PREFIX.nw > PREFIX.out'
calc = NWChem(**kwargs)

atom = {atom}
geom = {geom}

mol = Atoms(symbols=atom, positions=geom)
mol.set_calculator(calc)

# apply the constraints:
fix = {fix}
change = {change}
release = {release}

bonds = []
angles = []
dihedrals = []

for f in fix:
    if len(f) == 3:
        dist = mol.get_distance(f[0] - 1, f[1] - 1)  # ase enumerates atoms starting from zero
        bonds.append([dist, [f[0] - 1, f[1] - 1]])
    if len(f) == 3:
        angle = mol.get_angle(f[0] - 1, f[1] - 1, f[2] - 1)  # ase enumerates atoms starting from zero
        angle *= pi / 180.
        angles.append([angle, [f[0] - 1, f[1] - 1, f[2] - 1]])
    if len(f) == 4:
        dih = mol.get_dihedral(f[0] - 1, f[1] - 1, f[2] - 1, f[3] - 1)  # ase enumerates atoms starting from zero
        dih *= pi / 180.
        dihedrals.append([dih, [f[0] - 1, f[1] - 1, f[2] - 1, f[3] - 1]])
for c in change:
    if len(c) == 3:
        bonds.append([c[2], [c[0] - 1, c[1] - 1]])
    if len(c) == 4:
        angle = c[3] * pi / 180.
        angles.append([angle, [c[0] - 1, c[1] - 1, c[2] - 1]])
    if len(c) == 5:
        dih = c[4] * pi / 180.
        dihedrals.append([dih, [c[0] - 1, c[1] - 1, c[2] - 1, c[3] - 1]])

cons = FixInternals(bonds=bonds, angles=angles, dihedrals=dihedrals)
mol.set_constraint(cons)

cons.adjust_positions(mol, mol.positions)

try:
    dyn = BFGS(mol, trajectory='%s.traj' % label)
    dyn.run(fmax=0.05)
    db = connect('kinbot.db')
    db.write(mol, name=label, data={{'status': 'normal'}})
except RuntimeError as e:
    print('error')
    db = connect('kinbot.db')
    db.write(mol, name=label, data={{'status': 'error'}})

"""
try:
    mol.get_potential_energy() # use the NWChem optimizer (task optimize)
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
    db = connect('kinbot.db')
    db.write(mol, name = label, data = {{'status' : 'normal'}})
except RuntimeError, e: 
    db = connect('kinbot.db')
    db.write(mol, name = label, data = {{'status' : 'error'}})

"""
f = open(label + '.out', 'a')
f.write('done\n')
f.close()
