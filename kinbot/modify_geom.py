###################################################
##                                               ##
## This file is part of the KinBot code v2.0     ##
##                                               ##
## The contents are covered by the terms of the  ##
## BSD 3-clause license included in the LICENSE  ##
## file, found at the root.                      ##
##                                               ##
## Copyright 2018 National Technology &          ##
## Engineering Solutions of Sandia, LLC (NTESS). ##
## Under the terms of Contract DE-NA0003525 with ##
## NTESS, the U.S. Government retains certain    ##
## rights to this software.                      ##
##                                               ##
## Authors:                                      ##
##   Judit Zador                                 ##
##   Ruben Van de Vijver                         ##
##                                               ##
###################################################
"""
This class modifies a given geometry according to a set of coordinates that 
need to have a new value assigned. 

The optimization is done based on interatomic distances only
The deviations of the the distances are weighted by the inverse of the distance itself

"""

import sys, os, copy

import numpy as np

sys.dont_write_bytecode = True

import ase
from ase import Atoms
from ase.io import read,write,Trajectory
from ase.calculators.singlepoint import SinglePointCalculator

from BFGS import *
from geom import *
from stationary_pt import *
from cheminfo import *
from constants import *
from motif import *
from zmat import *


class cost_function():
    def __init__(self,coords):
        self.coords = coords
    def eval(self,x):
        """
        x is a vector of length 3N with N the number of atoms
        containing the cartesian coordinates [x1, y1, z1, x2, ..., xN, yN, zN]
        """
        e = 0
        for coord in self.coords:
            i = coord[0]
            j = coord[1]
            d = coord[2]
            weight = coord[3]
            dc = ( x[3*i] - x[3*j] )**2 + ( x[3*i+1] - x[3*j+1] )**2 + ( x[3*i+2] - x[3*j+2] )**2
            e += ( (dc - d) * weight )**2
        return e
        
    def gradient(self,x):
        grad = np.zeros(len(x))
        for coord in self.coords:
            i = coord[0]
            j = coord[1]
            d = coord[2]
            weight = coord[3]
            dc = ( x[3*i] - x[3*j] )**2 + ( x[3*i+1] - x[3*j+1] )**2 + ( x[3*i+2] - x[3*j+2] )**2
            
            grad[3*i] += 2 * ( (dc - d) * weight ) * 2 * (x[3*i] - x[3*j])
            grad[3*i+1] += 2 * ( (dc - d) * weight ) * 2 * (x[3*i+1] - x[3*j+1])
            grad[3*i+2] += 2 * ( (dc - d) * weight ) * 2 * (x[3*i+2] - x[3*j+2])
            
            grad[3*j] += 2 * ( (dc - d) * weight ) * 2 * (x[3*i] - x[3*j]) * -1
            grad[3*j+1] += 2 * ( (dc - d) * weight ) * 2 * (x[3*i+1] - x[3*j+1]) * -1
            grad[3*j+2] += 2 * ( (dc - d) * weight ) * 2 * (x[3*i+2] - x[3*j+2]) * -1
            
        return grad

def append_geom(natom,step,new_e,atom,x_new,grad,atoms_list,f_out = None):
    
    if f_out != None:
        f_out.write('%s\nPoint  %i  Energy=  %f\n'%(natom,step,new_e))
        for at in range(natom):
            f_out.write(atom[at] + ' ')
            [f_out.write(str(np.reshape(x_new,(natom,3))[at][i]) + '  ') for i in range(3)]
            f_out.write('\n')
        step += 1

    atoms = Atoms(symbols = atom, positions = np.reshape(x_new, (natom,3)))
    calc = SinglePointCalculator(atoms, energy= new_e, forces= 10. * np.reshape(grad,(natom,3)))
    atoms.set_calculator(calc)
    atoms_list.append(atoms)
    
    return step

def modify_coordinates(species,name,geom,changes,bond,natom,atom):
    """
    Geom is the geometry (n x 3 matrix with n the number of atoms)
    in cartesian coordinates
    
    Changes is a list of lists, each list containing the coordinates 
    and their new value (atom indices start at 0):
    To change a bond length: [atom1, atom2, bond_length]
    To change a bond angle: [neighbor1, central_atom, neighbor2, 
                             angle_in_degrees]
    To change a dihedral angle: [atom1, atom2, atom3, atom4, 
                                 dihedarl_angle_in_degrees]
    
    Bond is the bond matrix of the molecule
    """
    
    status = 1
    step = 1
    atoms_list = []
    
    count = 0
    fname = '{}_{}.xyz'.format(name,count)
    while os.path.exists(fname):
        count += 1
        fname = '{}_{}.xyz'.format(name,count)
    #f_out = open(fname,'w')
    f_out = None
    
    new_geom = copy.deepcopy(geom)
    step = append_geom(natom,step,0.,atom,new_geom,np.zeros((natom*3)),atoms_list,f_out = f_out)
    #change dihedrals, if necessary
    for ci in changes:
        if len(ci) == 5:
            zmat_atom, zmat_ref, zmat, zmatorder = make_zmat_from_cart(species, ci[:-1], natom, atom, new_geom, 2)
            orig_dih = zmat[3][2]
            new_dih = ci[-1]
            dih_diff = new_dih - orig_dih
            zmat[3][2] += dih_diff
            for i in range(4, natom):
                if zmat_ref[i][2] == 4:
                    zmat[i][2] += dih_diff
                if zmat_ref[i][2] == 1:
                    zmat[i][2] += dih_diff
            new_geom = make_cart_from_zmat(zmat, zmat_atom, zmat_ref, natom, atom, zmatorder)
            step = append_geom(natom,step,0.,atom,new_geom,np.zeros((natom*3)),atoms_list,f_out = f_out)
        #change angles, if necessary
        if len(ci) == 4:
            # original angle in radians
            orig_angle = calc_angle(new_geom[ci[0]],new_geom[ci[1]],new_geom[ci[2]])
            new_angle = np.radians(ci[-1]) #new angle in radians
            
            v1 = new_geom[ci[0]] - new_geom[ci[1]]
            v2 = new_geom[ci[2]] - new_geom[ci[1]]
            rot_ax = [0.,0.,0.]
            
            #create a vector perpendicular to v1 and v2
            #verify if points are collinear
            if np.linalg.norm(np.cross(v1,v2)) == 0:
                #rotate around any axis perpendicular to the axis along the three points:
                if v1[0] != 0 or v1[1] != 0:
                    rot_ax = [v1[1],-v1[0],0.]
                elif v1[0] != 0 or v1[2] != 0:
                    rot_ax = [v1[2],0.,-v1[0]]
                else: 
                    rot_ax = [1.,0.,0.]
            else:
                rot_ax = np.cross(v1,v2)
            
            rot_ax = rot_ax/np.linalg.norm(rot_ax)
            #rotate all the atoms on the side of the last atom
            st, ats, ats2 = divide_atoms(ci[2],ci[1],bond,natom,atom)
            if not st:
                status = 0
                break
            for atj in ats:
                new_geom[atj] = perform_rotation(new_geom[atj],new_geom[ci[1]],rot_ax,new_angle-orig_angle)
                step = append_geom(natom,step,1.,atom,new_geom,np.zeros((natom*3)),atoms_list,f_out = f_out)
    
    
    
    #change bond lengths and angles, if necessary
    if any([len(ci) == 3 or len(ci) == 4 for ci in changes]):
        """
        #first only change the dedicated bond lengths
        coords = get_coords(natom,bond,new_geom,changes,1)
        #optimize the geometry to meet the coords list
        x0 = np.reshape(new_geom,3*natom)
        cost_fct = cost_function(coords)
        opt = bfgs()
        x_opt, x_i = opt.optimize(cost_fct,x0)
        
        new_geom = np.reshape(x_opt,(natom, 3))
        for xi in x_i:
            gi = np.reshape(xi,(natom, 3))
            step = append_geom(natom,step,0.,atom,gi,np.zeros((natom*3)),atoms_list,f_out = f_out)
        """
        coords = get_coords(natom,bond,new_geom,changes,0)
        #optimize the geometry to meet the coords list
        x0 = np.reshape(new_geom,3*natom)
        cost_fct = cost_function(coords)
        opt = bfgs()
        x_opt, x_i = opt.optimize(cost_fct,x0)
        
        new_geom = np.reshape(x_opt,(natom, 3))
        for xi in x_i:
            gi = np.reshape(xi,(natom, 3))
            step = append_geom(natom,step,2.,atom,gi,np.zeros((natom*3)),atoms_list,f_out = f_out)
    
    #write(fname.replace('.xyz','.traj'),atoms_list)
    #f_out.close()
    
    return new_geom

def get_coords(natom,bond,geom,changes,mode):
    """
    list the (N*(N-1)) / 2 possible bond lengths and the value we optimize to
    mode = 0: include all interatomic distances
    mode = 1: only include the interatomic distances that need to be changed
    """
    coords = []
    for i in range(natom-1):
        for j in range(i+1, natom):
            is_change = 0 # the distance from i to j needs to be changed
            is_in_change = 0 # i or j needs to be changed, don't include this bond length
            change = []
            for ci in changes:
                if [i,j] == sorted([ci[0],ci[-2]]):
                    is_change = 1
                    change = ci
                if i in ci or j in ci:
                    is_in_change = 1
            if is_change:
                if len(change) == 3:
                    coords.append([i,j,change[-1]**2,1.,0])
                elif len(change) == 4:
                    #calculate the bond length that corresponds to the new angle
                    b1 = np.linalg.norm(geom[i]-geom[change[1]])
                    b2 = np.linalg.norm(geom[j]-geom[change[1]])
                    a = np.radians(change[-1])
                    d = b1**2 + b2**2 - 2*b1*b2*np.cos(a)
                    coords.append([i,j,d,10])
            else:
                if mode == 0:
                    d = np.linalg.norm(geom[i]-geom[j])**2
                    if np.sqrt(d) < 4.: #use a cutoff of 4 angstroms
                        if bond[i][j] > 0:
                            coords.append([i,j,d,1./d,1]) #use a larger weight for bonds
                        else:
                            if not is_in_change:
                                coords.append([i,j,d,.5/d,1])
    return coords

def divide_atoms(ati,atj,bond,natom,atom,forbidden = []):
    """
    This method divides the atoms in a molecule in two sets, 
    which are separated by a bond
    In the case of rings, the atoms are equally divided in the two sets, 
    which will change the bond length of the bond furthest away from 
    the given bond.
    Be careful when using this method for cyclic structures!
    """
    status = 1
    if bond[ati,atj] == 0:
        return 0, [ati], []
    
    #Get all the atoms on the side of ati
    visited = [ati]
    forbidden.append(atj)
    division = [ati]
    
    # check for cycles and cut them in half
    for ring_size in range(3,natom+1):
        motif = ['X' for at in range(ring_size)]
        inst = start_motif(motif,natom,bond,atom,-1,[])
        for ins in inst:
            if bond[ins[0]][ins[-1]] > 0:
                #cycle found
                if ins[0] == ati and ins[-1] == atj:
                    forbidden.append(ins[ring_size / 2])
                if ins[0] == atj and ins[-1] == ati:
                    forbidden.append(ins[- ring_size/2 - 1])
        if len(inst) == 0:
            break

    get_neighbors(ati,visited,forbidden,division,bond,natom)
    division2 = [x for x in range(natom) if x not in division]
    
    return status, division,division2

def get_neighbors(ati,visited,forbidden,division,bond,natom):
    for atj in range(natom):
        if not atj in visited and not atj in forbidden:
            if bond[atj,ati] > 0:
                division.append(atj)
                visited.append(atj)
                get_neighbors(atj,visited,forbidden,division,bond,natom)

def perform_rotation(at,center,axis,angle):
    #move center to origin: 
    at -= center
    
    #create rotation matrix: 
    a = math.cos(angle/2)
    b,c,d = -axis*math.sin(angle/2)
    aa,bb,cc,dd = a*a,b*b,c*c,d*d
    bc,ad,ac,ab,bd,cd = b*c,a*d,a*c,a*b,b*d,c*d
    rot_matrix = ([[aa+bb-cc-dd, 2*(bc+ad),   2*(bd-ac)],
                   [2*(bc-ad),   aa+cc-bb-dd, 2*(cd+ab)],
                   [2*(bd+ac),   2*(cd-ab),   aa+dd-bb-cc]])
    
    #perform the rotation: 
    at = np.dot(rot_matrix,at)

    #put the center back to its original coordinates:
    at += center

    return at
    
def main():
    smi = 'S[CH2]'
    obmol, structure = generate_3d_structure(smi)
    
    mult = 2
    charge = 0
    natom = len(obmol.atoms)
    
    structure = np.reshape(structure, (natom,4))
    atom = structure[:,0]
    geom = structure[:,1:4].astype(float)

    well0 = stationary_pt('well0')
    well0.geom = geom
    well0.characterize(natom,atom,2,0)
    well0.bond_mx(natom,atom)
    bond = well0.bond
    
    geom = [np.array(gi) for gi in geom]
    changes = [
    [1, 2, 1.0972959660511175],
    [0, 2, 2.4604284873750513],
    ]
    geom = modify_coordinates(well0,'test',geom,changes,bond,natom,atom)
    """
    changes = [
    [0,1,2,90.],
    [1,2,3,90.],
    ]
    geom = modify_coordinates(well0,'test',geom,changes,bond,natom,atom)
    """

if __name__ == "__main__":
    main()