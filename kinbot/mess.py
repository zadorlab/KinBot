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
import os
import numpy as np
import re
import subprocess
import time
import copy

import constants
import frequencies
"""
from constants import *
from frequencies import *
from qc import *
from rate_tables import *
import par
"""

class MESS:
    """
    Class that reads and write all MESS files and runs MESS
    """
    def __init__(self,par):
        self.par = par

    def write_header(self,reactant):
        """
        Create the header block for MESS
        """
        with open(self.par.par['tpldir'] + 'mess_header.tpl') as f:
            tpl = f.read()
        
        header = tpl.format(TemperatureList = ' '.join([str(ti) for ti in self.par.par['TemperatureList']]),
                            PressureList = ' '.join([str(pi) for pi in self.par.par['PressureList']]),
                            EnergyStepOverTemperature = self.par.par['EnergyStepOverTemperature'],
                            ExcessEnergyOverTemperature = self.par.par['ExcessEnergyOverTemperature'],
                            ModelEnergyLimit = self.par.par['ModelEnergyLimit'],
                            CalculationMethod = self.par.par['CalculationMethod'],
                            ChemicalEigenvalueMax = self.par.par['ChemicalEigenvalueMax'],
                            Reactant = reactant,
                            EnergyRelaxationFactor = self.par.par['EnergyRelaxationFactor'],
                            EnergyRelaxationPower = self.par.par['EnergyRelaxationPower'],
                            EnergyRelaxationExponentCutoff = self.par.par['EnergyRelaxationExponentCutoff'],
                            Epsilons = ' '.join([str(ei) for ei in self.par.par['Epsilons']]),
                            Sigmas = ' '.join([str(si) for si in self.par.par['Sigmas']]),
                            Masses = ' '.join([str(mi) for mi in self.par.par['Masses']]))
        return header
    def write(self,species):
        """
        #generate names for the mess calculations
        if len(obj.products) == 1:
            if not str(obj.products[0].chemid) in well_short_names:
                short_name = 'w'+str(len(well_short_names)+1)
                well_short_names[str(obj.products[0].chemid)] = short_name
            else:
                short_name = well_short_names[str(obj.products[0].chemid)]
            obj.products[0].short_name = short_name
        else:
            name = '_'.join(sorted([str(prod.chemid) for prod in obj.products]))
            if not name in bimol_short_names:
                bimol_short_names[name] = 'b'+str(len(bimol_short_names)+1)
            for prod in obj.products:
                if not str(prod.chemid) in fragment_short_names:
                    name = 'f' + str(len(fragment_short_names)+1)
                    fragment_short_names[str(prod.chemid)] = name
                else:
                    name = fragment_short_names[str(prod.chemid)]
                prod.short_name = name
        ts_short_names[instance_name] = 'ts' + str(len(ts_short_names)+1)
        
        #write the mess block for the ts
        write_mess_barrier(self.species, index, ts, obj.products, self.species.natom, self.species.atom, self.species.mult, self.species.charge, well_short_names, bimol_short_names, ts_short_names)
        
        #write the mess block for the products
        if len(obj.products) > 1:
            write_mess_bimol(obj.products,self.species, bimol_short_names, fragment_short_names)
        else:
            st_pt = obj.products[0]
            write_mess_well(st_pt, st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge, self.species, well_short_names)

        if self.par.par['me']:
            #write a complete mess file
            write_mess_all(self.species, barriers, products, well_short_names, bimol_short_names, fragment_short_names, ts_short_names)
            if not os.path.exists('mess/all.out'):
                run_mess()
            all_rates = read_mess()
            write_high_p_rates(all_rates, well_short_names,bimol_short_names)
        """
        return 0

    
    def write_all(self,species, barriers, products, well_short_names, bimol_short_names, fragment_short_names, ts_short_names):
        """
        Take the different blocks of the reactant, barrier, product and header
        and merge them into one file
        """
        
        #write the short names to a file:
        names_out = open('MESS_names.out','w')
        names_out.write('WELLS:\n')
        for w in well_short_names:
            names_out.write(w + '\t' + well_short_names[w] + '\n')
        names_out.write('BIMOLECULAR PRODUCTS:\n')
        for b in bimol_short_names:
            names_out.write(b + '\t' + bimol_short_names[b] + '\n')
        names_out.write('FRAGMENTS:\n')
        for f in fragment_short_names:
            names_out.write(f + '\t' + fragment_short_names[f] + '\n')
        names_out.write('BARRIERS:\n')
        for ts in ts_short_names:
            names_out.write(ts + '\t' + ts_short_names[ts] + '\n')
            
        names_out.close()
        
        #filter ts's with the same reactants and products:
        ts_unique = {} #key: ts name, value: [prod_name,energy]
        for index, instance in enumerate(species.reac_inst):
            if species.reac_ts_done[index] == -1:
                prod_name = '_'.join([str(pi.chemid) for pi in products[index]])
                energy = barriers[index].energy
                new = 1
                remove = []
                for ts in ts_unique:
                    if ts_unique[ts][0] == prod_name:
                        #check for the barrier with the lowest energy
                        if ts_unique[ts][1] > energy:
                            #remove the current barrier
                            remove.append(ts)
                        else:
                            new = 0
                for ts in remove:
                    ts_unique.pop(ts,None)
                if new:
                    ts_unique[species.reac_name[index]] = [prod_name,energy]
        
        #write the mess input file
        reactant = well_short_names[str(species.chemid)]
        header = write_mess_header(reactant)
        wells = ''
        for well in well_short_names:
            wells += open(well + '.mess','r').read()
            wells += '\n!****************************************\n'
        bimols = ''
        for bimol in bimol_short_names:
            bimols += open(bimol + '.mess','r').read()
            bimols += '\n!****************************************\n'
        tss = ''
        for ts in ts_unique:
            tss += open(ts + '.mess', 'r').read()
            tss += '\n!****************************************\n'

        dummy = open(self.par.par['tpldir'] + 'mess_dummy.tpl','r').read()
        dum = dummy.format(barrier = 'tsd', reactant = reactant, dummy = 'd1')
        
        if not os.path.exists('mess/'):
            os.mkdir('mess/')
        
        f_out = open('mess/all.inp','w')
        f_out.write(header + '\n!****************************************\n')
        f_out.write(wells)
        f_out.write(bimols)
        f_out.write(tss)
        f_out.write(dum)
        f_out.write('\n!****************************************\nEnd ! end kinetics\n')
        f_out.close()
        
    def write_bimol(self,species_list, well0, bimol_short_names, fragment_short_names):
        """ 
        Create the block for MESS for a bimolecular product.
        well0: reactant on this PES (zero-energy reference)
        """ 

        with open(self.par.par['tpldir'] + 'mess_bimol.tpl') as f:
            tpl = f.read()   
        with open(self.par.par['tpldir'] + 'mess_fragment.tpl') as f:
            fragment_tpl = f.read()
        with open(self.par.par['tpldir'] + 'mess_hinderedrotor.tpl') as f:
            rotor_tpl = f.read()
        with open(self.par.par['tpldir'] + 'mess_atom.tpl') as f:
            atom_tpl = f.read()
        
        fragments = ''
        for species in species_list:
            if species.natom > 1:
                rotors = []
                if self.par.par['rotor_scan']:
                    for i,rot in enumerate(species.dihed):
                        group = ' '.join([str(pi+1) for pi in frequencies.partition(species,rot,species.natom)[0][1:]])
                        axis = '{} {}'.format(str(rot[1]+1),str(rot[2]+1))
                        rotorsymm = species.sigma_int[rot[1]][rot[2]]
                        nrotorpot = species.hir.nrotation / rotorsymm
                        ens = species.hir_energies[i]
                        rotorpot = [(ei - ens[0])*constants.AUtoKCAL for ei in ens]
                        rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot[:species.hir.nrotation / rotorsymm]])
                        rotors.append(rotor_tpl.format( group = group,
                                                        axis = axis,
                                                        rotorsymm = rotorsymm,
                                                        nrotorpot = nrotorpot,
                                                        rotorpot = rotorpot
                                                        ))
                rotors = '\n'.join(rotors)
                
                freq = ''
                for i,fr in enumerate(species.reduced_freqs):
                    if i == 0:
                        freq += '{:.4f}'.format(fr)
                    elif i > 0 and i%3 == 0:
                        freq += '\n            {:.4f}'.format(fr)
                    else:
                        freq += '    {:.4f}'.format(fr)
                    
                geom = ''
                for i,at in enumerate(species.atom):
                    if i > 0:
                        geom += '            '
                    x,y,z = species.geom[i]
                    geom += '{} {:.6f} {:.6f} {:.6f}\n'.format(at,x,y,z)
                
                energy = (  ( species.energy + species.zpe )- ( well0.energy + well0.zpe) ) * constants.AUtoKCAL
                
                fragments += fragment_tpl.format(   chemid = fragment_short_names[str(species.chemid)] + ' ! ' + str(species.chemid),
                                                    natom = species.natom,
                                                    geom = geom,
                                                    symm = float(species.sigma_ext) / float(species.nopt),
                                                    nfreq = len(species.reduced_freqs),
                                                    freq = freq,
                                                    hinderedrotor = rotors,
                                                    nelec = 1,
                                                    charge = species.charge,
                                                    mult = species.mult,
                                                    )
                fragments += '\n'
            else:
                # atom template
                fragments += atom_tpl.format(   chemid = fragment_short_names[str(species.chemid)] + ' ! ' + str(species.chemid),
                                                element = species.atom[0],
                                                nelec = 1,
                                                charge = species.charge,
                                                mult = species.mult
                                                )
                fragments += '\n'
        
        if par.pes:
            energy = '{ground_energy}'
        else:
            energy = ( sum([sp.energy for sp in species_list]) + sum([sp.zpe for sp in species_list]) -
                     ( well0.energy + well0.zpe ) ) * constants.AUtoKCAL

        name = '_'.join(sorted([str(species.chemid) for species in species_list]))
        bimol = tpl.format( chemids = bimol_short_names[name] + ' ! ' + name,
                            fragments = fragments,
                            ground_energy = energy
                            )

        f = open(name + '.mess', 'w')
        f.write(bimol)
        f.close()


    def write_well(self,species, well0, well_short_names):
        """ 
        Create the block for MESS for a well.
        
        well0: reactant on this PES (zero-energy reference)
        
        """ 
        with open(self.par.par['tpldir'] + 'mess_well.tpl') as f:
            tpl = f.read()   
        with open(self.par.par['tpldir'] + 'mess_hinderedrotor.tpl') as f:
            rotor_tpl = f.read()   
        
        rotors = []
        if self.par.par['rotor_scan']:
            for i,rot in enumerate(species.dihed):
                group = ' '.join([str(pi+1) for pi in frequencies.partition(species,rot,species.natom)[0][1:]])
                axis = '{} {}'.format(str(rot[1]+1),str(rot[2]+1))
                rotorsymm = species.sigma_int[rot[1]][rot[2]]
                nrotorpot = species.hir.nrotation / rotorsymm
                ens = species.hir.hir_energies[i]
                rotorpot = [(ei - ens[0])*constants.AUtoKCAL for ei in ens]
                rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot[:species.hir.nrotation / rotorsymm]])
                rotors.append(rotor_tpl.format( group = group,
                                                axis = axis,
                                                rotorsymm = rotorsymm,
                                                nrotorpot = nrotorpot,
                                                rotorpot = rotorpot
                                                ))
        rotors = '\n'.join(rotors)
        
        freq = ''
        for i,fr in enumerate(species.reduced_freqs):
            if i == 0:
                freq += '{:.4f}'.format(fr)
            elif i > 0 and i%3 == 0:
                freq += '\n            {:.4f}'.format(fr)
            else:
                freq += '    {:.4f}'.format(fr)
            
        geom = ''
        for i,at in enumerate(species.atom):
            if i > 0:
                geom += '            '
            x,y,z = species.geom[i]
            geom += '{} {:.6f} {:.6f} {:.6f}\n'.format(at,x,y,z)
        
        if self.par.par['pes']:
            energy = '{zeroenergy}'
        else:
            energy = (  ( species.energy + species.zpe )- ( well0.energy + well0.zpe) ) * constants.AUtoKCAL
        
        mess = tpl.format(  chemid = well_short_names[str(species.chemid)] + ' ! ' + str(species.chemid),
                            natom = species.natom,
                            geom = geom,
                            symm = float(species.sigma_ext) / float(species.nopt),
                            nfreq = len(species.reduced_freqs),
                            freq = freq,
                            hinderedrotor = rotors,
                            nelec = 1,
                            charge = species.charge,
                            mult = species.mult,
                            zeroenergy = energy
                            )
        
        f = open(str(species.chemid) + '.mess', 'w')
        f.write(mess)
        f.close()


    def write_barrier(self, species, index, ts, prods, natom, atom, mult, charge, well_short_names, bimol_short_names, ts_short_names):
        """ 
        Create the block for a MESS barrier.
        
        species: well0 of this PES search
        index: reaction index of the current barrier
        ts: stationary_pt object of the current barrier
        prods: products of the current reaction (list of stationary_pt objects)
        natom: number of atoms in the ts
        atom: list of element symbols for the atoms in the ts
        mult: multiplicity of the ts
        charge: charge of the ts
        """ 
        
        

        with open(self.par.par['tpldir'] + 'mess_ts.tpl') as f:
            tpl = f.read()   
        with open(self.par.par['tpldir'] + 'mess_hinderedrotor.tpl') as f:
            rotor_tpl = f.read()   
        with open(self.par.par['tpldir'] + 'mess_tunneling.tpl') as f:
            tun_tpl = f.read()   
        
        rotors = []
        if self.par.par['rotor_scan']:
            for i,rot in enumerate(ts.dihed):
                group = ' '.join([str(pi+1) for pi in frequencies.partition(ts,rot,natom)[0][1:]])
                axis = '{} {}'.format(str(rot[1]+1),str(rot[2]+1))
                rotorsymm = ts.sigma_int[rot[1]][rot[2]]
                nrotorpot = par.nrotation / rotorsymm
                ens = ts.hir_energies[i]
                rotorpot = [(ei - ens[0])*constants.AUtoKCAL for ei in ens]
                rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot[:par.nrotation / rotorsymm]])
                rotors.append(rotor_tpl.format( group = group,
                                                axis = axis,
                                                rotorsymm = rotorsymm,
                                                nrotorpot = nrotorpot,
                                                rotorpot = rotorpot
                                                ))
        rotors = '\n'.join(rotors)
        
        freq = ''
        for i,fr in enumerate(ts.reduced_freqs[1:]):
            if i == 0:
                freq += '{:.4f}'.format(fr)
            elif i > 0 and i%3 == 0:
                freq += '\n            {:.4f}'.format(fr)
            else:
                freq += '    {:.4f}'.format(fr)
            
        geom = ''
        for i,at in enumerate(atom):
            if i > 0:
                geom += '            '
            x,y,z = ts.geom[i]
            geom += '{} {:.6f} {:.6f} {:.6f}\n'.format(at,x,y,z)
            
        
        barriers = [
            ((ts.energy + ts.zpe) - (species.energy + species.zpe)) * constants.AUtoKCAL,
            ((ts.energy + ts.zpe) - sum([(pr.energy + pr.zpe) for pr in prods]))*constants.AUtoKCAL,
        ]
        if any([bi < 0 for bi in barriers]):
            tun = ''
        else:
            tun = tun_tpl.format(   cutoff = min(barriers),
                                    imfreq = -ts.reduced_freqs[0],
                                    welldepth1 = barriers[0],
                                    welldepth2 = barriers[1])
        
        if par.pes:
            energy = '{zeroenergy}'
        else:
            energy = (  ( ts.energy + ts.zpe )- ( species.energy + species.zpe) ) * constants.AUtoKCAL
        
        if len(prods) == 1:
            prod_name = well_short_names[str(prods[0].chemid)]
        else:
            long_name = '_'.join(sorted([str(pi.chemid) for pi in prods]))
            prod_name = bimol_short_names[long_name]
        
        
        mess = tpl.format(  rxn_name = ts_short_names[species.reac_name[index]],
                            chemid_reac = well_short_names[str(species.chemid)],
                            chemid_prod = prod_name,
                            long_rxn_name = species.reac_name[index],
                            natom = natom,
                            geom = geom,
                            symm = float(ts.sigma_ext) / float(ts.nopt),
                            nfreq = len(ts.reduced_freqs) - 1,
                            freq = freq,
                            hinderedrotor = rotors,
                            tunneling = tun,
                            
                            nelec = 1,
                            charge = charge,
                            mult = mult,
                            zeroenergy = energy
                            )
        
        f = open(species.reac_name[index] + '.mess', 'w')
        f.write(mess)
        f.close()

    def write_high_p_rates(all_rates,well_short_names, bimol_short_names):
        f_out = open('rates.out','w')
        if len(all_rates) > 0:
            f_out.write('temperatures\t\t')
            f_out.write('\t'.join([str(ti) for ti in all_rates[0].temperatures]))
            f_out.write('\n')

            for rates in all_rates:
                true_rxn = 1
                
                if rates.reactant in well_short_names.values():
                    for name in well_short_names:
                        if well_short_names[name] == rates.reactant:
                            reac_name = name
                elif rates.reactant in bimol_short_names.values():
                    for name in bimol_short_names:
                        if bimol_short_names[name] == rates.reactant:
                            reac_name = name
                else:
                    true_rxn = 0
                
                if rates.product in well_short_names.values():
                    for name in well_short_names:
                        if well_short_names[name] == rates.product:
                            prod_name = name
                elif rates.product in bimol_short_names.values():
                    for name in bimol_short_names:
                        if bimol_short_names[name] == rates.product:
                            prod_name = name
                else:
                    true_rxn = 0

                if true_rxn:
                    s = reac_name + '\t'
                    s += prod_name + '\t'
                    s += '\t'.join([str(ri) for ri in rates.rates[-1]])
                    s += '\n'
                    f_out.write(s)

        
        f_out.close()

    def read(self):
        lines = open('mess/all.out','r').read().split('\n')
        index = 0
        all_rates = []
        while index < len(lines):
            line = lines[index]
            if line.startswith('Temperature-Pressure Rate Tables'):
                #start reading this section
                while index < len(lines):
                    line = lines[index]
                    if '->' in line:
                        if len(line.split('->')) > 1:
                            reactant = line.split('->')[0]
                            product = line.split('->')[1]
                            rates = Rates(reactant,product)
                            rates.temperatures = [float(ti) for ti in lines[index+2].split()[1:]]
                            #read_temperatures()
                            rates.pressures = []
                            #rates.temperatures = par.TemperatureList
                            #rates.pressures = copy.deepcopy(par.PressureList)
                            #rates.temperatures = [300+100*j for j in range(18)]
                            #rates.pressures = [76, 760, 7600]
                            #rates.pressures.append(float('inf'))
                            index += 3 #go to the start of the rates
                            while len(lines[index]) > 0:
                                line = lines[index]
                                if line.split()[0] == 'O-O':
                                    rates.pressures.append(float('inf'))
                                else:
                                    rates.pressures.append(float(line.split()[0]))
                                
                                r = []
                                for j in range(len(rates.temperatures)):
                                    rj = line[7+13*j:7+13*(j+1)]
                                    if '***' in rj:
                                        r.append(0.)
                                    else:
                                        r.append(float(rj))
                                rates.rates.append(r)
                                
                                #r = line.split()[1:]
                                #for j in range(len(r)):
                                #    if r[j] == '***': r[j] = 0.
                                #rates.rates.append([float(ri) for ri in r])
                                index += 1
                            #i -= 1
                            all_rates.append(rates)
                    index += 1
            index += 1
        return all_rates

        
    def run(self):
        """
        write a pbs file for the mess/all.inp mess input file
        submit the pbs file to the queue
        wait for the mess run to finish
        """
        with open(self.par.par['tpldir'] + 'pbs_mess.tpl') as f:
            tpl = f.read()
        pbs = open('run_mess.pbs','w')
        pbs.write(tpl.format(name = 'mess', ppn = par.ppn, queue_name = par.queue_name, dir = 'mess'))
        pbs.close()
        
        command = ['qsub','run_mess.pbs']
        process = subprocess.Popen(command,shell=False,stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = process.communicate()
        pid = out.split('\n')[0].split('.')[0]
        
        while 1:
            devnull = open(os.devnull, 'w')
            command = 'qstat -f | grep ' + '"Job Id: ' + pid + '"' + ' > /dev/null'
            if int(subprocess.call(command, shell = True, stdout=devnull, stderr=devnull)) == 0: 
                time.sleep(1)
            else:
                break
        
        return 0
