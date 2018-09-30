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
import os,sys
import time
import logging
import numpy as np
import matplotlib.pyplot as plt

from constants import *
from stationary_pt import *
from zmat import *
from qc import *
import par

def generate_hir_geoms(species, natom, atom, mult, charge, cart, wellorts):
    species.hir_status = []
    species.hir_energies = []
    species.hir_geoms = []
    while len(species.hir_status) < len(species.dihed):
        species.hir_status.append([-1 for i in range(par.nrotation)])
        species.hir_energies.append([-1 for i in range(par.nrotation)])
        species.hir_geoms.append([[] for i in range(par.nrotation)])

    for rotor in range(len(species.dihed)):
        cart = np.asarray(cart)       
        zmat_atom, zmat_ref, zmat, zmatorder = make_zmat_from_cart(species, rotor, natom, atom, cart, 0)

        #first element has same geometry ( TODO: this shouldn't be recalculated)
        cart_new = make_cart_from_zmat(zmat, zmat_atom, zmat_ref, natom, atom, zmatorder)

        fi = [(zi+1) for zi in zmatorder[:4]]
        qc_hir(species,cart_new,wellorts,natom,atom,mult,charge,rotor,0,[fi])
        for ai in range(1,par.nrotation):
            ang = 360. / float(par.nrotation)
            zmat[3][2] += ang
            for i in range(4, natom):
                if zmat_ref[i][2] == 4:
                    zmat[i][2] += ang 
                if zmat_ref[i][2] == 1:
                    zmat[i][2] += ang 
            cart_new = make_cart_from_zmat(zmat, zmat_atom, zmat_ref, natom, atom, zmatorder)
            qc_hir(species,cart_new,wellorts,natom,atom,mult,charge,rotor,ai,[fi])
    return 0

def test_hir(species,natom,atom,mult,charge,wellorts):
    for rotor in range(len(species.dihed)):
        for ai in range(par.nrotation):
            if species.hir_status[rotor][ai] == -1:
                if wellorts:
                    job = 'hir/' + species.name + '_hir_' + str(rotor) + '_' + str(ai).zfill(2)
                else:
                    job = 'hir/' + str(species.chemid) + '_hir_' + str(rotor) + '_' + str(ai).zfill(2)
                err, geom = get_qc_geom(job, natom)
                if err == 1: #still running
                    continue
                elif err == -1: #failed
                    species.hir_status[rotor][ai] = 1
                    species.hir_energies[rotor][ai] = -1
                    species.hir_geoms[rotor][ai] = geom
                else:
                    
                    #check if all the bond lenghts are within 15% or the original bond lengths
                    if equal_geom(species.bond,species.geom,geom,0.15):
                        err, energy = get_qc_energy(job)
                        species.hir_status[rotor][ai] = 0
                        species.hir_energies[rotor][ai] = energy
                        species.hir_geoms[rotor][ai] = geom
                    else:
                        species.hir_status[rotor][ai] = 1
                        species.hir_energies[rotor][ai] = -1
                        species.hir_geoms[rotor][ai] = geom
    return 0

    
def check_hir(species, natom, atom, mult, charge, wellorts, wait = 0):
    """
    Check for hir calculations and optionally wait for them to finish
    """
    while 1:
        #check if all the calculations are finished
        test_hir(species,natom,atom,mult,charge,wellorts)
        if all([all([test >= 0 for test in status]) for status in species.hir_status]):
            for rotor in range(len(species.dihed)):
                if wellorts:
                    job = species.name + '_hir_' + str(rotor)
                else:
                    job = str(species.chemid) + '_hir_' + str(rotor)
                
                angles = [i * 2 * np.pi / float(par.nrotation) for i in range(par.nrotation)]
                #write profile to file
                write_profile(species,rotor,job,atom,natom)
                species.hir_fourier.append(fourier_fit(job,angles,species.hir_energies[rotor],species.hir_status[rotor],plot_fit = 0))
            return 1
        else:
            if wait:
                time.sleep(1)
            else:
                return 0

def write_profile(species,rotor,job,atom,natom):
    """
    Write a molden-readable file with the HIR scan (geometries and energies)
    """
    file = open('hir/' + job + '.xyz','w')
    for i in range(par.nrotation):
        s = str(natom) + '\n'
        s += 'energy = ' + str(species.hir_energies[rotor][i]) + '\n'
        for j,at in enumerate(atom):
            x,y,z = species.hir_geoms[rotor][i][j]
            s += '{} {:.8f} {:.8f} {:.8f}\n'.format(at,x,y,z)
        file.write(s)
    file.close()
    
def fourier_fit(job,angles,energies,status,plot_fit = 0):
    """
    Create a alternative fourier formulation of a hindered rotor
    (Vanspeybroeck et al.)
    
    profile, the angles are in radians and the eneries in 
    kcal per mol
    plot_fit: plot the profile and the fit to a png
    """
    n_terms = 6 #the number of sine and cosine terms
    
    ang = [angles[i] for i in range(len(status)) if status[i] == 0]
    ens = [(energies[i] - energies[0])*AUtoKCAL for i in range(len(status)) if status[i] == 0]
    
    if len(ens) < par.par.nrotation - 2:
        #more than two points are off
        logging.warning("Hindered rotor potential has more than 2 failures for " + job)
    
    X = np.zeros((len(ang), 2 * n_terms))
    for i,ai in enumerate(ang):
        for j in range(n_terms):
            X[i][j] = (1 - np.cos((j+1) * ai))
            X[i][j+n_terms] = np.sin((j+1) * ai)

    A = np.linalg.lstsq(X,np.array(ens))[0]
    
    for i,si in enumerate(status):
        if si == 1:
            energies[i] = energies[0] + get_fit_value(A,angles[i])/AUtoKCAL
    
    if plot_fit:
        #fit the plot to a png file
        plt.plot(ang,ens,'ro')
        fit_angles = [i * 2. * np.pi / 360 for i in range(360)]
        fit_energies = [get_fit_value(A,ai) for ai in fit_angles]
        plt.plot(fit_angles,fit_energies)
        plt.xlabel('Dihedral angle [radians]')
        plt.ylabel('Energy [kcal/mol]')
        plt.savefig('hir_profiles/{}.png'.format(job))
        plt.clf()
    return A

def get_fit_value(A,ai):
    """
    Get the fitted energy
    """
    e = 0.
    n_terms = (len(A)) / 2
    for j in range(n_terms):
        e += A[j] * ( 1 - np.cos((j+1) * ai))
        e += A[j+n_terms] * np.sin((j+1) * ai)
    
    return e
    

def main():
    """
    Calculate the 1D hindered rotor profiles
    Create a fourier fit representation of the profile
    """



if __name__ == "__main__":
    main()

