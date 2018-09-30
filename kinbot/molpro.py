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
import os, sys

from constants import *
import par

def create_molpro_input(species, natom, atom, mult, charge, wellorts):
    """
    Create the input for molden
    
    species: stationary point object
    natom: number of atoms in the stationary point
    atom: elements of all the atoms in the stationary point
    mult: multiplicity of the stationary point
    charge: charge of the stationary point
    wellorts: 1 for transition states, 0 for molecules
    """
    
    with open(par.tpldir + 'molpro.tpl') as f:
        file = f.read()
    
    fname = str(species.chemid)
    if wellorts: fname = species.name
    
    geom = ''
    nelectron = 0
    for i,at in enumerate(atom):
        x,y,z = species.geom[i]
        geom += '{} {:.8f} {:.8f} {:.8f}\n'.format(at,x,y,z)
        nelectron += znumber[at]
    
    nelectron -= charge
    
    outf = open('molpro/' + fname + '.inp','w')
    outf.write(file.format( name = fname,
                            natom = natom,
                            atom = atom,
                            geom = geom,
                            nelectron = nelectron,
                            spin = mult - 1,
                            charge = charge
                            ))
    outf.close()

def get_molpro_energy(species,wellorts):
    """
    Verify if there is a molpro output file and if yes, read the energy
    """
    fname = str(species.chemid)
    if wellorts: fname = species.name
    
    status = os.path.exists('molpro/' + fname + '.out')
    
    if status:
        with open('molpro/' + fname + '.out') as f:
            lines = f.readlines()
        
        for index, line in enumerate(reversed(lines)):
            if 'SETTING MYENA' in line:
                return 1, float(line.split()[2])
    else:
        return 0, -1
