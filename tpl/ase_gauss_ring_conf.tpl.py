from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.optimize import LBFGS
from ase.constraints import FixInternals
from ase.db import connect
from kinbot import reader_gauss

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

mol = Atoms(symbols={atom}, positions={geom})
kwargs = {kwargs}

Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

bonds, angles, dihedrals = reader_gauss.constraint(mol, {fix}, {change})
constraints = FixInternals(bonds=bonds, angles_deg=angles, dihedrals_deg=dihedrals)
mol.set_constraint(constraints)

dyn = LBFGS(atoms=mol, trajectory='ringopt.traj')
try:
    dyn.run(fmax=0.01, steps=400)
    e = mol.get_potential_energy()
    iowait(logfile, 'gauss')
    data = {{'energy': e, 'status': 'normal'}}
except (RuntimeError, ValueError):
    data = {{'status': 'error'}}

db.write(mol, name=label, data=data)

with open(logfile,'a') as f:
    f.write('done\n')
