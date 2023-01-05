from ase import Atoms
from ase.optimize import LBFGS
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot.ase_modules.constraints import FixInternals  # New
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

mol = Atoms(symbols={atom}, positions={geom})
kwargs = {kwargs}

Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

bonds, angles, dihedrals = reader_gauss.constraint(mol, {fix}, {change})
# if version.parse(ase.__version__) >= version.parse("3.21"):
constraints = FixInternals(bonds=bonds, angles_deg=angles,
                           dihedrals_deg=dihedrals)
# else:
#     constraints = FixInternals(bonds, angles, dihedrals)
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
