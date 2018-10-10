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

import constants
import geometry
import zmatrix

class HIR:
    """
    Class that does all the steps for the HIR calculations of one species
    """
    def __init__(self,species,qc,par):
        """
        species: instance of StationaryPoint
        qc: instance of QuantumChemistry
        par: instance of Parameters
        """
        self.species = species
        self.qc = qc
        
        #number of points along one scan
        self.nrotation = par.par['nrotation']
        
        # -1 (not finished), 0 (successful) or 1 (failed) for each HIR scan point
        self.hir_status = []
        # energies of all the HIR scan points
        self.hir_energies = []
        # Fourier fit of each scan
        self.hir_fourier = []
        # all the geometries of the HIR scan points
        self.hir_geoms = []

    def generate_hir_geoms(self,cart):
        """
        Generate the initial geometries of the points along the scans
        """
        #re-initialize the lists incase of a restart of the HIR scans
        self.hir_status = []
        self.hir_energies = []
        self.hir_geoms = []
        
        while len(self.hir_status) < len(self.species.dihed):
            self.hir_status.append([-1 for i in range(self.nrotation)])
            self.hir_energies.append([-1 for i in range(self.nrotation)])
            self.hir_geoms.append([[] for i in range(self.nrotation)])

        for rotor in range(len(self.species.dihed)):
            cart = np.asarray(cart)       
            zmat_atom, zmat_ref, zmat, zmatorder = zmatrix.make_zmat_from_cart(self.species, rotor, cart, 0)

            #first element has same geometry ( TODO: this shouldn't be recalculated)
            cart_new = zmatrix.make_cart_from_zmat(zmat, zmat_atom, zmat_ref, self.species.natom, self.species.atom, zmatorder)

            fi = [(zi+1) for zi in zmatorder[:4]]
            self.qc.qc_hir(self.species,cart_new,rotor,0,[fi])
            for ai in range(1,self.nrotation):
                ang = 360. / float(self.nrotation)
                zmat[3][2] += ang
                for i in range(4, self.species.natom):
                    if zmat_ref[i][2] == 4:
                        zmat[i][2] += ang 
                    if zmat_ref[i][2] == 1:
                        zmat[i][2] += ang 
                cart_new = zmatrix.make_cart_from_zmat(zmat, zmat_atom, zmat_ref, self.species.natom, self.species.atom, zmatorder)
                self.qc.qc_hir(self.species,cart_new,rotor,ai,[fi])
        return 0

    def test_hir(self):
        for rotor in range(len(self.species.dihed)):
            for ai in range(self.nrotation):
                if self.hir_status[rotor][ai] == -1:
                    if self.species.wellorts:
                        job = 'hir/' + self.species.name + '_hir_' + str(rotor) + '_' + str(ai).zfill(2)
                    else:
                        job = 'hir/' + str(self.species.chemid) + '_hir_' + str(rotor) + '_' + str(ai).zfill(2)
                    err, geom = self.qc.get_qc_geom(job, self.species.natom)
                    if err == 1: #still running
                        continue
                    elif err == -1: #failed
                        self.hir_status[rotor][ai] = 1
                        self.hir_energies[rotor][ai] = -1
                        self.hir_geoms[rotor][ai] = geom
                    else:
                        
                        #check if all the bond lenghts are within 15% or the original bond lengths
                        if geometry.equal_geom(self.species.bond,self.species.geom,geom,0.15):
                            err, energy = self.qc.get_qc_energy(job)
                            self.hir_status[rotor][ai] = 0
                            self.hir_energies[rotor][ai] = energy
                            self.hir_geoms[rotor][ai] = geom
                        else:
                            self.hir_status[rotor][ai] = 1
                            self.hir_energies[rotor][ai] = -1
                            self.hir_geoms[rotor][ai] = geom
        return 0

        
    def check_hir(self, wait = 0):
        """
        Check for hir calculations and optionally wait for them to finish
        """
        while 1:
            #check if all the calculations are finished
            self.test_hir()
            if all([all([test >= 0 for test in status]) for status in self.hir_status]):
                for rotor in range(len(self.species.dihed)):
                    if self.species.wellorts:
                        job = self.species.name + '_hir_' + str(rotor)
                    else:
                        job = str(self.species.chemid) + '_hir_' + str(rotor)
                    
                    angles = [i * 2 * np.pi / float(self.nrotation) for i in range(self.nrotation)]
                    #write profile to file
                    self.write_profile(rotor,job)
                    self.hir_fourier.append(self.fourier_fit(job,angles,rotor,plot_fit = 0))
                return 1
            else:
                if wait:
                    time.sleep(1)
                else:
                    return 0

    def write_profile(self,rotor,job):
        """
        Write a molden-readable file with the HIR scan (geometries and energies)
        """
        file = open('hir/' + job + '.xyz','w')
        for i in range(self.nrotation):
            s = str(self.species.natom) + '\n'
            s += 'energy = ' + str(self.hir_energies[rotor][i]) + '\n'
            for j,at in enumerate(self.species.atom):
                x,y,z = self.hir_geoms[rotor][i][j]
                s += '{} {:.8f} {:.8f} {:.8f}\n'.format(at,x,y,z)
            file.write(s)
        file.close()
        
    def fourier_fit(self,job,angles,rotor,plot_fit = 0):
        """
        Create a alternative fourier formulation of a hindered rotor
        (Vanspeybroeck et al.)
        
        profile, the angles are in radians and the eneries in 
        kcal per mol
        plot_fit: plot the profile and the fit to a png
        """
        energies = self.hir_energies[rotor]
        status = self.hir_status[rotor]
        
        n_terms = 6 #the number of sine and cosine terms
        
        ang = [angles[i] for i in range(len(status)) if status[i] == 0]
        ens = [(energies[i] - energies[0])*constants.AUtoKCAL for i in range(len(status)) if status[i] == 0]
        
        if len(ens) < self.nrotation - 2:
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
                energies[i] = energies[0] + self.get_fit_value(A,angles[i])/constants.AUtoKCAL
        
        if plot_fit:
            #fit the plot to a png file
            plt.plot(ang,ens,'ro')
            fit_angles = [i * 2. * np.pi / 360 for i in range(360)]
            fit_energies = [self.get_fit_value(A,ai) for ai in fit_angles]
            plt.plot(fit_angles,fit_energies)
            plt.xlabel('Dihedral angle [radians]')
            plt.ylabel('Energy [kcal/mol]')
            plt.savefig('hir_profiles/{}.png'.format(job))
            plt.clf()
        return A

    def get_fit_value(self,A,ai):
        """
        Get the fitted energy
        """
        e = 0.
        n_terms = (len(A)) / 2
        for j in range(n_terms):
            e += A[j] * ( 1 - np.cos((j+1) * ai))
            e += A[j+n_terms] * np.sin((j+1) * ai)
        return e
        
