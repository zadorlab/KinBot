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
import sys
import os
import copy
import numpy as np
import re
import logging 

sys.dont_write_bytecode = True

from stationary_pt import *
from conformers import *
from zmat import *
from par import *
from modify_geom import *
from geom import *

def create_zmat(mol):
    """
    Create the zmat of the molecule mol
    """
    zmat_atom, zmat_ref, zmat, zmatorder = make_zmat_from_cart_all_dihedrals(mol.bond,mol.cycle,mol.dihed,mol.conf_dihed,par.natom,par.atom,mol.geom,0)

    zfile = open('zmat_test.out','w')
    
    write_zmat_molden(zfile,zmat_atom,zmat_ref,zmat)
    
    #make 1D rotations
    for i,rotor in enumerate(mol.dihed):
        #get the appropriate dihedral position
        for j in range(3,par.natom):
            indices = [zmatorder[j],zmatorder[zmat_ref[j][0]-1],zmatorder[zmat_ref[j][1]-1],zmatorder[zmat_ref[j][2]-1]]
            if indices == rotor or indices == rotor[::-1]:
                break

        file = '/home/rvandev/rotors/rotor_' + str(i) + '_' + '_'.join([str(ri) for ri in rotor]) + '.xyz'
        f = open(file,'w')
        #do 36 iterations of 10 degrees
        for k in range(36):
            zmat_copy = copy.deepcopy(zmat)
            zmat_copy[j][2] += float(k * 10)
            cart = make_cart_from_zmat(zmat_copy, zmat_atom, zmat_ref, par.natom, par.atom, zmatorder)
            cart = translate_and_rotate(cart,par.atom,rotor[1],rotor[2])
            
            f.write(write_cart(cart,par.atom))
        f.close()

def translate_and_rotate(cart,atom,i,j):
    """
    translate the molecule as such that the first rotor atom i is the center of mass
    and the ij vector is along the x axis 
    """
    
    #translate the molecule:
    trans = copy.deepcopy(cart[i])
    for ci in cart:
        ci -= trans
        
    #rotate the molecule
    end_vec = [0.,0.,1.]
    angle = calc_angle(cart[j], cart[i], end_vec)
    if angle != 0:
        axis = np.cross(cart[j],end_vec)
        axis = axis/np.linalg.norm(axis)
        a = math.cos(angle/2)
        b,c,d = -axis*math.sin(angle/2)
        aa,bb,cc,dd = a*a,b*b,c*c,d*d
        bc,ad,ac,ab,bd,cd = b*c,a*d,a*c,a*b,b*d,c*d
        rot_matrix = ([[aa+bb-cc-dd, 2*(bc+ad),   2*(bd-ac)],
                       [2*(bc-ad),   aa+cc-bb-dd, 2*(cd+ab)],
                       [2*(bd+ac),   2*(cd-ab),   aa+dd-bb-cc]])
        for i in range(len(cart)):
            cart[i] = np.dot(rot_matrix,cart[i])

    return cart

def write_cart(geom,atom):
    s = '%i\n'%len(geom)
    s += 'Energy = 1.0\n'
    
    for index in range(len(geom)):
        s += '%s %.6f %.6f %.6f\n'%(atom[index],geom[index][0],geom[index][1],geom[index][2])
    return s

def main():
    """
    Generate a molecule and have a look at its zmat
    
    Next, generate conformers and start 1D rotor scans
    """
    
    logging.basicConfig(filename='zmat_test.log',level=logging.INFO)
    logging.info('Starting KinBot in %s at %s'%(workdir,datetime.datetime.now()))
    
    mol = stationary_pt('well0')
    mol.geom = par.read_input('zmat_test_input.inp')
    
    file = open('zmat_test.xyz','w')
    file.write(write_cart(mol.geom,par.atom))
    file.write('\n\n\n')
    file.close()
    
    mol.characterize(par.natom, par.atom, par.mult, par.charge)

    logging.debug(mol.cycle_chain)
    logging.debug(mol.geom)

    create_zmat(mol)

if __name__ == "__main__":
    main()
    