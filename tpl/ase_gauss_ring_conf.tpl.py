import numpy as np
import ase
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.optimize.pcobfgs import PCOBFGS
from ase.db import connect

db = connect('{working_dir}/kinbot.db')

mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}

Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.set_calculator(calc)

fix = {fix}
change = {change}

bonds = []
angles = []
dihedrals = []

for fi in fix:
    if len(fi) == 2:
        #careful: atom indices in the fix lists start at 1
        bondlength = mol.get_distance(fi[0]-1, fi[1]-1)
        bonds.append([bondlength,[fi[0]-1, fi[1]-1]])
    if len(fi) == 3:
        #careful: atom indices in the fix lists start at 1
        angle = mol.get_angle(fi[0]-1, fi[1]-1, fi[2]-1) * np.pi / 180
        angles.append([angle,[fi[0]-1, fi[1]-1, fi[2]-1]])
    if len(fi) == 4:
        #careful: atom indices in the fix lists start at 1
        dihed = mol.get_dihedral(fi[0]-1, fi[1]-1, fi[2]-1, fi[3]-1) * np.pi / 180
        dihedrals.append([dihed,[fi[0]-1, fi[1]-1, fi[2]-1, fi[3]-1]])
for ci in change:
    if len(ci) == 3:
        #careful: atom indices in the fix lists start at 1
        bondlength = ci[2]
        bonds.append([bondlength,[ci[0]-1, ci[1]-1]])
    if len(ci) == 4:
        #careful: atom indices in the fix lists start at 1
        angle = ci[3] * np.pi / 180
        angles.append([angle,[ci[0]-1, ci[1]-1, ci[2]-1]])
    if len(ci) == 5:
        #careful: atom indices in the fix lists start at 1
        dihed = ci[4] * np.pi / 180
        dihedrals.append([dihed,[ci[0]-1, ci[1]-1, ci[2]-1, ci[3]-1]])

dyn = PCOBFGS(mol,
              trajectory='{label}.traj',
              bonds=bonds,
              angles=angles,
              dihedrals=dihedrals,
              force_consistent=False)

try:
    dyn.run(fmax=0.01, steps=400)
    e = mol.get_potential_energy()
    data = {{'energy': e, 'status': 'normal'}}
except RuntimeError:
    data = {{'status': 'error'}}

db.write(mol, name='{label}', data=data)

with open(f'{label}.log','a') as f:
    f.write('done\n')
