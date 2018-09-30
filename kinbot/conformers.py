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
import random
import time
import numpy as np

from modify_geom import *
from cheminfo import *
from stationary_pt import *
from zmat import *
from qc import *
from geom import *
import par


def generate_ring_conformers(species, natom, atom, mult, charge, cart, wellorts):
    """
    Generate the conformers of a cyclic structure
    by randomly sampling the dihedrals of the ring
    """
    
    confs = []
    
    all_dihs = []
    all_cycdih = []
    
    #iterate the different rings in the species
    for cyc in species.cycle_chain:
        if len(cyc) > 3: #three membered rings don't have conformers
            dihs = [] #list of the ring dihedrals
            for i,at in enumerate(cyc):
                dihs.append([cyc[i-3],cyc[i-2],cyc[i-1],cyc[i]])
            
            #define the flatness of the ring by the sum of the absolute values of the
            #dihedrals along the ring divided by the number of atoms in the ring
            cycdih = 0.
            for dih in dihs:
                cycdih += np.abs(calc_dihedral(cart[dih[0]], cart[dih[1]], cart[dih[2]], cart[dih[3]])[0])
            cycdih /= float(len(cyc))
            #randomly select N-3 dihedrals, with N the number of dihedrals in the ring
            all_dihs.extend(random.sample(dihs,len(dihs) - 3))
            all_cycdih.append(cycdih)
        else:
            all_cycdih.append(0.0)
    
    return start_ring_conformer_search( species, natom, atom, mult, charge, [cart], 
                                        wellorts, all_dihs, all_cycdih, [], 0, 0)

def start_ring_conformer_search(species, natom, atom, mult, charge, carts, wellorts, all_dihs, 
                                all_cycdih, dih_values, dih_index, conf_index):
    

    if dih_index == len(all_dihs):
        return carts

    for cart in carts:
        #define the dihedrals that have been sampled already and fix those
        fix = []
        for i in range(dih_index):
            dih = all_dihs[i]
            fix.append([dih[0] + 1,dih[1] + 1,dih[2] + 1,dih[3] + 1])
        
        #select the proper dihedral change for the current ring
        cyc_index = 0
        for i, cyc in enumerate(species.cycle_chain):
            if all([di in cyc for di in all_dihs[dih_index]]):
                cyc_index = i
                break
        cycdih = all_cycdih[cyc_index]
    
        #optimize three searches, with - cycdih, 0 and + cycdih
        dih = all_dihs[dih_index]
        d = calc_dihedral(cart[dih[0]], cart[dih[1]], cart[dih[2]], cart[dih[3]])[0]
        
        for new_dih in [-cycdih, 0., + cycdih]:
            relaxed_scan = []
            diff = (new_dih - d) 
            nsteps = int(np.abs(diff)) / 10 + 1 #take steps of 10 degrees maximum
            if nsteps < 1:
                nsteps = 1
            relaxed_scan.append([dih[0] + 1,dih[1] + 1,dih[2] + 1,dih[3] + 1, nsteps, diff / nsteps])
            qc_ring_conf(species, cart, wellorts, natom, atom, mult, charge, relaxed_scan, fix, conf_index)
            conf_index += 1
        
    
    #wait for all the current runs to finish
    status, new_carts = check_conformers(species, conf_index, natom, atom, mult, charge,wait = 1, ts = wellorts, ring = 1)
    dih_index += 1
    new_carts = start_ring_conformer_search(species, natom, atom, mult, charge, new_carts, 
                                wellorts, all_dihs, all_cycdih, dih_values, dih_index, conf_index)
    
    return new_carts

    
    
    
def start_ring_conformer    (species, natom, atom, mult, charge, cart, 
                            wellorts, all_dihs, all_cycdih, dih_values, 
                            dih_index = 0, conf_index = 0):
    """
    Combine all the ring conformer deviations
    """
    if dih_index == len(all_dihs): 
        relaxed_scan = []
        for i, dih in enumerate(all_dihs):
            d = calc_dihedral(cart[dih[0]], cart[dih[1]], cart[dih[2]], cart[dih[3]])[0]
            diff = (dih_values[i] - d) 
            nsteps = int(diff) / 10 +1
            relaxed_scan.append([dih[0] + 1,dih[1] + 1,dih[2] + 1,dih[3] + 1, nsteps, diff / nsteps])
        

        
        #qc_ring_conf(species, cart, wellorts, natom, atom, mult, charge, relaxed_scan, conf_index)
        conf_index += 1
        return conf_index
    
    
    dih_values.append(0.)
    new_dih_values = [di for di in dih_values]
    conf_index = start_ring_conformer(  species, natom, atom, mult, charge, cart, 
                                        wellorts, all_dihs, all_cycdih, new_dih_values, 
                                        dih_index = dih_index + 1, conf_index = conf_index)
    
    cyc_index = 0
    for i, cyc in enumerate(species.cycle_chain):
        if all([di in cyc for di in all_dihs[dih_index]]):
            cyc_index = i
            break
    cycdih = all_cycdih[cyc_index]
    
    dih_values[-1] = cycdih
    new_dih_values = [di for di in dih_values]
    conf_index = start_ring_conformer(  species, natom, atom, mult, charge, cart, 
                                        wellorts, all_dihs, all_cycdih, new_dih_values, 
                                        dih_index = dih_index + 1, conf_index = conf_index)
    

    dih_values[-1] = -cycdih
    new_dih_values = [di for di in dih_values]
    conf_index = start_ring_conformer(  species, natom, atom, mult, charge, cart, 
                                        wellorts, all_dihs, all_cycdih, new_dih_values, 
                                        dih_index = dih_index + 1, conf_index = conf_index)
    
    return conf_index

def generate_conformers(species, natom, atom, mult, charge, rotor, cart, conf, wellorts):
    """
    Generate guesses for all of the canonical conformers.
    This is a recursive routine to generate them.
    rotor: the rotor number in the order it was discovered
    conf: counter of the conformers
    wellorts = 0 for wells, 1 for ts's
    """
    
    if len(species.conf_dihed) > 4:
        return generate_conformers_random_sampling(species, natom, atom, mult, charge, cart, conf, wellorts)
    
    if rotor == len(species.conf_dihed): 
        qc_conf(species, cart, wellorts, natom, atom, mult, charge, conf)
        conf += 1
        return conf

    cart = np.asarray(cart)                                    
    zmat_atom, zmat_ref, zmat, zmatorder = make_zmat_from_cart(species, rotor, natom, atom, cart, 1)

    rotor += 1
    cart0 = make_cart_from_zmat(zmat, zmat_atom, zmat_ref, natom, atom, zmatorder)   
    conf = generate_conformers(species, natom, atom, mult, charge, rotor, cart0, conf, wellorts)
    
    zmat[3][2] += 120.
    for i in range(4, natom):
        if zmat_ref[i][2] == 4:
            zmat[i][2] += 120. 
        if zmat_ref[i][2] == 1:
            zmat[i][2] += 120. 
    cart1 = make_cart_from_zmat(zmat, zmat_atom, zmat_ref, natom, atom, zmatorder)
    conf = generate_conformers(species, natom, atom, mult, charge, rotor, cart1, conf, wellorts)
    
    zmat[3][2] += 120.  
    for i in range(4, natom):
        if zmat_ref[i][2] == 4:
            zmat[i][2] += 120. 
        if zmat_ref[i][2] == 1:
            zmat[i][2] += 120.       
    cart2 = make_cart_from_zmat(zmat, zmat_atom, zmat_ref, natom, atom, zmatorder)
    conf = generate_conformers(species, natom, atom, mult, charge, rotor, cart2, conf, wellorts)
    
    return conf

def generate_conformers_random_sampling(species, natom, atom, mult, charge, ini_cart, conf, wellorts, nconfs = 100):
    """
    Generate a random sampling of each dihedral for a number nconfs of conformers
    """
    for ni in range(nconfs):
        cart = copy.deepcopy(ini_cart)
        if ni == 0:
            sample = [0. for di in species.conf_dihed]
        else:
            sample = [random.choice([0.,120., 240.]) for di in species.conf_dihed]
        for rotor in range(len(species.conf_dihed)):
            zmat_atom, zmat_ref, zmat, zmatorder = make_zmat_from_cart(species, rotor, natom, atom, cart, 1)
            zmat[3][2] += sample[rotor]
            for i in range(4, natom):
                if zmat_ref[i][2] == 4:
                    zmat[i][2] += sample[rotor]
                if zmat_ref[i][2] == 1:
                    zmat[i][2] += sample[rotor]
            cart = make_cart_from_zmat(zmat, zmat_atom, zmat_ref, natom, atom, zmatorder)
        qc_conf(species, cart, wellorts, natom, atom, mult, charge, conf)
        conf += 1
    return conf


def test_conformer(species, conf, natom, atom, mult, charge, ts = 0, ring = 0):
    """
    Test whether a conformer has the same bond matrix as the original structure.
    Returns the conformer object and -1 if not yet finished, 0 if same, and 1 if not.
    """
    r = ''
    if ring: r = 'r'
    
    if ts:
        job = 'conf/' + species.name + '_' + r + str(conf).zfill(par.zf)
    else:
        job = 'conf/' + str(species.chemid) + '_' + r + str(conf).zfill(par.zf)
    
    status, geom = get_qc_geom(job, natom)
    if status == 1: #still running
        return np.zeros((natom,3)), -1
    elif status == -1: #conformer search failed
        return np.zeros((natom,3)), 1
    else:
        #check if all the bond lenghts are withing 10% or the original bond lengths
        if equal_geom(species.bond,species.geom,geom,0.10):
            return geom, 0
        else:
            return np.zeros((natom,3)), 1
        """
        if not ts:
            temp = stationary_pt('temporary_species')
            temp.geom = geom
            temp.characterize(natom, atom, mult, charge)
            if all([all([temp.bond[i][j] == species.bond[i][j] for j in range(natom)]) for i in range(natom)]):
                return geom, 0
            else:
                return np.zeros((natom,3)), 1
        else:
            # TODO: find a way to compare ts geometries
            return geom, 0
        """

def check_conformers(species, conf, natom, atom, mult, charge,wait = 0, ts = 0, ring = 0):
    """
    Check if the conformer optimizations finished.
    Test them, and submit frequency calculations.
    Then select the lowest energy one.
    
    returns:
    *status: 0 if still running, 1 if finished
    * if ring  = 0: geometry of lowest energy conformer
      if ring  = 1: geometries of all the conformers
    
    wait: wait for all the conformer calculations to finish before returning anything
    ts: 1 for transition states and 0 for wells
    ring: 1 for cyclic searches, 0 for open searches
    """ 
    
    if ring:
        if len(species.cyc_conf_status) < conf:
            for i in range(len(species.cyc_conf_status),conf):
                species.cyc_conf_status.append(-1)
        status = species.cyc_conf_status
    else:
        if len(species.conf_status) < conf:
            for i in range(len(species.conf_status),conf):
                species.conf_status.append(-1)
        status = species.conf_status
    
    r = ''
    if ring: r = 'r'
    
    while 1:
        #check if conformational search is finished
        for i,si in enumerate(status):
            if si == -1:
                status[i] = test_conformer(species,i,natom,atom,mult,charge, ts = ts, ring = ring)[1]
        if all([si >= 0 for si in status]):
            lowest_energy = species.energy
            lowest_e_geom = species.geom
            final_geoms = [] # list of all final conformer geometries
            geoms = [species.geom] # list used for intermediate ring conformer generation
            energies = []
            for ci in range(conf): 
                si = status[ci]
                if si == 0: #this is a valid confomer
                    if ts:
                        job = 'conf/' + species.name + '_' + r + str(ci).zfill(par.zf)
                    else:
                        job = 'conf/' + str(species.chemid) + '_' + r + str(ci).zfill(par.zf)
                    err, energy = get_qc_energy(job)
                    err, geom = get_qc_geom(job,natom)
                    geoms.append(geom)
                    final_geoms.append(geom)
                    energies.append(energy)
                    if energy < lowest_energy:
                        lowest_energy = energy
                        lowest_e_geom = geom
                else:
                    energies.append(0.)
                    final_geoms.append(np.zeros((natom,3)))
            if ring:
                return 1, geoms
            else:
                write_profile(species,ts,status,final_geoms,energies,atom,natom)
                return 1, lowest_e_geom
        else:
            if wait:
                time.sleep(1)
            else:
                return 0, np.zeros((natom,3))

def write_profile(species,ts,status,final_geoms,energies,atom,natom):
    """
    Write a molden-readable file with the CONF analysis (geometries and energies)
    """
    if ts:
        file = open('conf/' + species.name + '.xyz','w')
    else:
        file = open('conf/' + str(species.chemid) + '.xyz','w')


    for i,st in enumerate(status):
        s = str(natom) + '\n'
        s += 'energy = ' + str(energies[i]) + '\n'
        for j,at in enumerate(atom):
            x,y,z = final_geoms[i][j]
            s += '{} {:.8f} {:.8f} {:.8f}\n'.format(at,x,y,z)
        if st == 0: #valid conformer:
            file.write(s)

    file.close()

def main():
    """
    This code generates conformers.
    """

if __name__ == "__main__":
    main()
