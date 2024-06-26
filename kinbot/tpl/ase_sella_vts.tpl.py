import numpy as np
from ase import Atoms
from ase.db import connect
import time
import os
import copy
from sella import Sella, Constraints, Internals
from kinbot.ase_modules.calculators.{code} import {Code}

import rmsd

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

mol = Atoms(symbols={atom}, positions={init_geom})
kwargs = {kwargs}
kwargs['chk'] = '{label}'.replace('vrctst/', '')
mol.calc = {Code}(**kwargs)
e = mol.get_potential_energy()
kwargs['guess'] = 'Read'
mol.calc = {Code}(**kwargs)
mol_prev = copy.deepcopy(mol)

# RELAXED
scan_coo = {scan_coo}
scan_coo_equiv = {equiv}
scan_dist = np.linalg.norm(mol.positions[scan_coo[0]] - mol.positions[scan_coo[1]])

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')
sella_kwargs = {sella_kwargs}
cons = Constraints(mol)
cons.fix_bond((scan_coo[0], scan_coo[1]))
# do not allow other equivalent point pairs to be closer than the one scanned
for sci in scan_coo_equiv[0]:
    for scj in scan_coo_equiv[1]:
        if [sci, scj] == scan_coo or [scj, sci] == scan_coo:
            continue
        cons.fix_bond((sci, scj), target=scan_dist, comparator='gt')
internals = Internals(mol)
internals.find_all_bonds()
# add scanned bond if missing from internals - may not have any effect in the end
hit = False
for bond in internals.internals['bonds']:
    if (bond.indices[0] == scan_coo[0] and bond.indices[1] == scan_coo[1]) or \
       (bond.indices[1] == scan_coo[0] and bond.indices[0] == scan_coo[1]):
        hit = True 
        break
if not hit:
    internals.add_bond((scan_coo[0], scan_coo[1]))

internals.find_all_angles()
internals.find_all_dihedrals()
opts = []
opt = Sella(mol, 
            order=0, 
            constraints=cons,
            internal=False,
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            delta0=1.3e-3,
            **sella_kwargs)
opts.append(opt)
mol.calc.label = '{label}'
bonds = {bonds}
distances = np.array([np.linalg.norm(mol.positions[bond[0]] - mol.positions[bond[1]]) for bond in bonds])

energies = []
model = 0
fmax = 1.e-4
while 1:
    ok = True
    last = True  # take the last geometry, otherwise the one before that
    if {asymptote}:
        break
    try:
        for i in opts[model].irun(fmax=fmax, steps=100):
            # due to dummy atom, constraint is lost
            if abs(np.linalg.norm(mol.positions[scan_coo[0]] - mol.positions[scan_coo[1]]) - scan_dist) > 0.01:
                ok = False
                print('constraint lost')
                break
            if rmsd.kabsch_rmsd(np.array({init_geom}), mol.positions, translate=True) > {scan_deviation}:
                last = False
                print('rmsd is too large, optimization is stopped')
                break
            curr_distances = np.array([np.linalg.norm(mol.positions[bond[0]] - mol.positions[bond[1]]) for bond in bonds])
            ratio = curr_distances / distances
            if any([True if ri > 1.1 else False for ri in ratio]):
                last = False
                print('a bond is more than 10% stretched, optimization is stopped')
                #break
            mol_prev = copy.deepcopy(mol)
            e = mol.get_potential_energy() 
            energies.append(e)
            # when forces don't fully converge, but energy doesn't change anymore
            if len(energies) > 11 and np.std(energies[-10:]) < 1.e-7:
                break
    except (RuntimeError, AssertionError):
        ok = False
        pass

    if not ok:
        mol = copy.deepcopy(mol_prev)
        cons = Constraints(mol)
        cons.fix_bond((scan_coo[0], scan_coo[1]))
        for sci in scan_coo_equiv[0]:
            for scj in scan_coo_equiv[1]:
                if [sci, scj] == scan_coo or [scj, sci] == scan_coo:
                    continue
                cons.fix_bond((sci, scj), target=scan_dist, comparator='lt')
        # give up internal coo optimization and losen convergence
        fmax *= 3.
        opt = Sella(mol,
                    order=0,
                    constraints=cons,
                    internal=False,
                    trajectory='{label}.traj', 
                    logfile='{label}_sella.log',
                    **sella_kwargs)
        opts.append(opt)
        model += 1
    elif model == 3:
        break
    else:
        break

if not last:
    mol = mol_prev

db.write(mol, name='{label}',
         data={{'energy': e, 'status': 'normal'}})

with open('{label}.log', 'a') as f:
    f.write('done\n')

# FROZEN

# create fragment geometries from the relaxed scan and save scan atoms' position and index
coo_A = np.array([mol.positions[i] for i in {frag_maps}[0]])
coo_B = np.array([mol.positions[i] for i in {frag_maps}[1]])
try:
    scan_atom_A = {frag_maps}[0].index(scan_coo[0])
    scan_atom_B = {frag_maps}[1].index(scan_coo[1])
except ValueError:
    scan_atom_A = {frag_maps}[0].index(scan_coo_equiv[1])
    scan_atom_B = {frag_maps}[1].index(scan_coo_equiv[0])
scan_pos_A = coo_A[scan_atom_A]
scan_pos_B = coo_B[scan_atom_B]

# Move frozen fragments to own centroids (translation only)
coo_A_fr = np.array({froz_A_geom})
coo_B_fr = np.array({froz_B_geom})

# set pivot atom to 1
A_weight = copy.copy({frag_bonds_0}[scan_atom_A])
B_weight = copy.copy({frag_bonds_1}[scan_atom_B])
A_weight[scan_atom_A] = 1
B_weight[scan_atom_B] = 1 

# rotate fragments with weight
coo_A_fr, _ = rmsd.kabsch_weighted_fit(coo_A_fr, coo_A, W=A_weight)
coo_B_fr, _ = rmsd.kabsch_weighted_fit(coo_B_fr, coo_B, W=B_weight)

# move fragments back to scan position
coo_A_fr -= coo_A_fr[scan_atom_A] - scan_pos_A
coo_B_fr -= coo_B_fr[scan_atom_B] - scan_pos_B

# reassemble full geometry from frozen fragments
frozen_geom = np.zeros((len({atom}), 3))
for fi, fa in enumerate({frag_maps}[0]):
    frozen_geom[fa] = coo_A_fr[fi]
for fi, fb in enumerate({frag_maps}[1]):
    frozen_geom[fb] = coo_B_fr[fi]
    
# here creating the run for the frozen version of this geometry based on the fully optimized one
label = '{label}_fr'
logfile = '{label}_fr.log'
kwargs.pop('opt', None)
kwargs.pop('addsec', None)
kwargs['label'] = label
calc = Gaussian(**kwargs)
mol = Atoms(symbols={atom}, positions=frozen_geom)
mol.calc = calc
try:
    e = mol.get_potential_energy()  # single point energy
    db.write(mol, name=label, data={{'energy': e, 'status': 'normal'}})
except:
    db.write(mol, name=label, data={{'energy': 10., 'status': 'normal'}})

time.sleep(1)  
with open(logfile, 'a') as f:
    f.write('done\n')
