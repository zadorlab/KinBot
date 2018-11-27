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
from __future__ import division 
import sys
import os
import numpy as np
import re
import subprocess
import time
import copy
import xml.etree.cElementTree as ET
import xml.dom.minidom as minidom

import pkg_resources

import constants
import frequencies
import cheminfo

class MESMER:
    """
    Class that reads and write all MESMER files and runs MESMER
    """
    def __init__(self,par,species):
        self.par = par
        self.species = species

    def write_header(self):
        """
        Create the header block for MESMER
        """
        #Read the header template
        header_file = pkg_resources.resource_filename('tpl', 'mess_header.tpl')
        with open(header_file) as f:
            tpl = f.read()
        
        header = tpl.format(TemperatureList = ' '.join([str(ti) for ti in self.par.par['TemperatureList']]),
                            PressureList = ' '.join([str(pi) for pi in self.par.par['PressureList']]),
                            EnergyStepOverTemperature = self.par.par['EnergyStepOverTemperature'],
                            ExcessEnergyOverTemperature = self.par.par['ExcessEnergyOverTemperature'],
                            ModelEnergyLimit = self.par.par['ModelEnergyLimit'],
                            CalculationMethod = self.par.par['CalculationMethod'],
                            ChemicalEigenvalueMax = self.par.par['ChemicalEigenvalueMax'],
                            Reactant = self.well_names[self.species.chemid],
                            EnergyRelaxationFactor = self.par.par['EnergyRelaxationFactor'],
                            EnergyRelaxationPower = self.par.par['EnergyRelaxationPower'],
                            EnergyRelaxationExponentCutoff = self.par.par['EnergyRelaxationExponentCutoff'],
                            Epsilons = ' '.join([str(ei) for ei in self.par.par['Epsilons']]),
                            Sigmas = ' '.join([str(si) for si in self.par.par['Sigmas']]),
                            Masses = ' '.join([str(mi) for mi in self.par.par['Masses']]))
        return header
    
    def write_input(self):
        """
        write the input for all the wells, bimolecular products and barriers
        both in a separate file, as well as in one large ME file
        """
        #create short names for all the species, bimolecular products and barriers
        #header = self.write_header()

        #filter ts's with the same reactants and products:
        ts_unique = {} #key: ts name, value: [prod_name,energy]
        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] == -1:
                prod_name = '_'.join([str(pi.chemid) for pi in reaction.products])
                energy = reaction.ts.energy
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
                    ts_unique[reaction.instance_name] = [prod_name,energy]

        root = ET.Element(  'me:mesmer',{'xmlns':'http://www.xml-cml.org/schema',
                            'xmlns:me':'http://www.chem.leeds.ac.uk/mesmer',
                            'xmlns:xsi':'http://www.w3.org/2001/XMLSchema-instance'})
        
        ET.SubElement(root,'me:title').text = str(self.species.chemid)
        mollist = ET.SubElement(root,'moleculeList')
        reaclist = ET.SubElement(root,'reactionList')
        
        #write the mess input for the different blocks
        wells = [self.species.chemid]
        self.write_well(self.species,mollist)
        
        for index,reaction in enumerate(self.species.reac_obj):
            if reaction.instance_name in ts_unique:
                for st_pt_opt in reaction.prod_opt:
                    st_pt = st_pt_opt.species
                    if not st_pt.chemid in wells:
                        self.write_well(st_pt,mollist)
                self.write_barrier(reaction,mollist,reaclist)
        
        #add a bath gas molecule
        molecule = ET.SubElement(mollist, 'molecule', {'id':'He'})
        propertylist = ET.SubElement(molecule, 'propertyList')
        #add epsilon
        epsilon = ET.SubElement(propertylist, 'property', {'dictRef':'me:epsilon'})
        ET.SubElement(epsilon, 'scalar').text = '48.0'
        #add sigma
        sigma = ET.SubElement(propertylist, 'property', {'dictRef':'me:sigma'})
        ET.SubElement(sigma, 'scalar').text = '3.90'
        #add the molecular weight
        mw = ET.SubElement(propertylist, 'property', {'dictRef':'me:mw'})
        ET.SubElement(mw, 'scalar').text = '28.0'
        
        #add the conditions
        conditions = ET.SubElement(root,'me:conditions')
        ET.SubElement(conditions, 'me:bathGas').text = 'He'
        pts = ET.SubElement(conditions,'me:PTs')
        for p in self.par.par['PressureList']:
            for t in self.par.par['TemperatureList']:
                ET.SubElement(pts,'me:PTpair',{'units':'Torr','p':str(p),'t':str(t),'precision':'d'})
        
        #add the model parameters
        modelparams = ET.SubElement(root,'me:modelParameters')
        ET.SubElement(modelparams,'me:grainSize',{'units':'cm-1'}).text = '50'
        ET.SubElement(modelparams,'me:energyAboveTheTopHill').text = '30.0'
        
        #add the control parameters
        control = ET.SubElement(root,'me:control')
        ET.SubElement(control,'me:printSpeciesProfile')
        ET.SubElement(control,'me:testRateConstants')
        ET.SubElement(control,'me:eigenvalues').text = '3'
        
        st = ET.tostring(root,'utf-8')
        st = minidom.parseString(st)
        fout = open('test.xml','w')
        fout.write(st.toprettyxml(indent = ' '))
        fout.close()

        return 0

    def write_well(self,species,mollist):
        """ 
        Create the block for MESS for a well.
        
        well0: reactant on this PES (zero-energy reference)
        
        """ 
        molecule = ET.SubElement(mollist, 'molecule', {'id':str(species.chemid)})
        atomarray = ET.SubElement(molecule, 'atomArray')
        for i,at in enumerate(species.atom):
            args = {'id':'a{}'.format(i+1)}
            args['elementType'] = at
            args['x3'] = '{:.8f}'.format(species.geom[i][0])
            args['y3'] = '{:.8f}'.format(species.geom[i][1])
            args['z3'] = '{:.8f}'.format(species.geom[i][2])
            at = ET.SubElement(atomarray, 'atom', args)
        
        bondarray = ET.SubElement(molecule, 'bondArray')
        bond_id = 1
        bond_ref = {}
        for i in range(species.natom-1):
            for j in range(i+1,species.natom):
                if species.bond[i][j] > 0:
                    bond_name = 'b{}'.format(bond_id)
                    bond_ref['{}_{}'.format(i,j)] = bond_name
                    args = {'id':'b{}'.format(bond_id)}
                    args['atomRefs2']="a{} a{}".format(i+1,j+1)
                    args['order']="{}".format(species.bond[i][j])
                    b = ET.SubElement(bondarray,'bond',args)
                    bond_id += 1
        propertylist = ET.SubElement(molecule, 'propertyList')
        
        #add the zpe
        if self.par.par['pes']:
            energy = '{zeroenergy}'
        else:
            energy = (  ( species.energy + species.zpe )- ( self.species.energy + self.species.zpe) ) * constants.AUtoKCAL
        zpe = ET.SubElement(propertylist, 'property', {'dictRef':'me:ZPE'})
        ET.SubElement(zpe, 'scalar', {'units':'kcal/mol'}).text = str(energy)
        
        #add the multiplicity
        mult = ET.SubElement(propertylist, 'property', {'dictRef':'me:spinMultiplicity'})
        ET.SubElement(mult, 'scalar').text = str(species.mult)
        
        #add the external symmetry number
        sigma = ET.SubElement(propertylist, 'property', {'dictRef':'me:symmetryNumber'})
        ET.SubElement(sigma, 'scalar').text = str(species.sigma_ext / species.nopt)
        
        #add the vibrational frequencies
        vibs = ET.SubElement(propertylist, 'property', {'dictRef':'me:vibFreqs'})
        ET.SubElement(vibs, 'array', {'units':'cm-1'}).text = ' '.join([str(fi) for fi in species.kinbot_freqs])
        
        #add the rotor potentials
        if self.par.par['rotor_scan'] and len(species.dihed) > 0:
            qmrotors = ET.SubElement(molecule, 'me:DOSCMethod',{'xsi:type':'me:QMRotors'})
            for i,rot in enumerate(species.dihed):
                hinderedrotor = ET.SubElement(qmrotors, 'me:ExtraDOSCMethod',{'xsi:type':'me:HinderedRotorQM1D'})
                ET.SubElement(hinderedrotor, 'me:bondRef').text = bond_ref[ '_'.join(sorted([str(rot[1]),str(rot[2])])) ]
                pot = ET.SubElement(hinderedrotor, 'me:HinderedRotorPotential', {'format':'numerical','units':'kcal/mol','expansionSize':'7','useSineTerms':'yes'})
                ens = species.hir.hir_energies[i]
                rotorpot = [(ei - ens[0])*constants.AUtoKCAL for ei in ens]
                for j,ri in enumerate(rotorpot):
                    ET.SubElement(pot, 'me:PotentialPoint', {'angle':str(360/species.hir.nrotation*j),'potential':str(ri)})
                ET.SubElement(pot, 'me:PotentialPoint', {'angle':str(360.0),'potential':str(rotorpot[0])})
                ET.SubElement(hinderedrotor, 'me:periodicity').text = str(species.sigma_int[rot[1]][rot[2]])
        else:
            qmrotors = ET.SubElement(molecule, 'me:DOSCMethod',{'name':'QMRotors'})
        
        #add the energy transfer model
        etransfer = ET.SubElement(molecule, 'me:energyTransferModel',{'xsi:type':'me:ExponentialDown'})
        ET.SubElement(etransfer, 'me:deltaEDown', {'units':'cm-1'}).text = '-1'


    def write_barrier(self, reaction, mollist, reaclist):
        """ 
        Create the block for a MESS barrier.
        """ 
        
        rxn = ET.SubElement(reaclist, 'reaction', {'id':reaction.instance_name})
        
        #add the reactant refs
        reactant = ET.SubElement(rxn, 'reactant')
        ET.SubElement(reactant, 'molecule', {'ref':str(self.species.chemid), 'role':'modelled'})
        
        # add the product refs
        if len(reaction.products) == 1:
            product = ET.SubElement(rxn, 'product')
            ET.SubElement(product, 'molecule', {'ref':str(reaction.products[0].chemid), 'role':'modelled'})
        else:
            for st_pt in reaction.products:
                product = ET.SubElement(rxn, 'product')
                ET.SubElement(product, 'molecule', {'ref':str(st_pt.chemid), 'role':'sink'})
        
        #add the transition state ref
        ts = ET.SubElement(rxn, 'me:transitionState')
        ET.SubElement(ts, 'molecule', {'ref':reaction.instance_name, 'role':'transitionState'})
        
        #add the MCRCMethod
        ET.SubElement(rxn, 'me:MCRCMethod',{'name':'RRKM'})
        
        #add the tunneling
        barriers = [
            ((reaction.ts.energy + reaction.ts.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL,
            ((reaction.ts.energy +reaction.ts.zpe) - sum([(opt.species.energy + opt.species.zpe) for opt in reaction.prod_opt]))*constants.AUtoKCAL,
        ]
        if all([bi > 0 for bi in barriers]):
            ET.SubElement(rxn, 'me:tunneling',{'name':'Eckart'})
        
        #add the transition state to the molecule list
        molecule = ET.SubElement(mollist, 'molecule', {'id':reaction.instance_name}) 
        atomarray = ET.SubElement(molecule, 'atomArray')
        for i,at in enumerate(reaction.ts.atom):
            args = {'id':'a{}'.format(i+1)}
            args['elementType'] = at
            args['x3'] = '{:.8f}'.format(reaction.ts.geom[i][0])
            args['y3'] = '{:.8f}'.format(reaction.ts.geom[i][1])
            args['z3'] = '{:.8f}'.format(reaction.ts.geom[i][2])
            at = ET.SubElement(atomarray, 'atom', args)
        
        bondarray = ET.SubElement(molecule, 'bondArray')
        bond_id = 1
        bond_ref = {}
        for i in range(reaction.ts.natom-1):
            for j in range(i+1,reaction.ts.natom):
                if reaction.ts.bond[i][j] > 0:
                    bond_name = 'b{}'.format(bond_id)
                    bond_ref['{}_{}'.format(i,j)] = bond_name
                    args = {'id':'b{}'.format(bond_id)}
                    args['atomRefs2']="a{} a{}".format(i+1,j+1)
                    args['order']="{}".format(reaction.ts.bond[i][j])
                    b = ET.SubElement(bondarray,'bond',args)
                    bond_id += 1
        propertylist = ET.SubElement(molecule, 'propertyList')
        
        #add the zpe
        if self.par.par['pes']:
            energy = '{zeroenergy}'
        else:
            energy = (  ( reaction.ts.energy + reaction.ts.zpe )- ( self.species.energy + self.species.zpe) ) * constants.AUtoKCAL
        zpe = ET.SubElement(propertylist, 'property', {'dictRef':'me:ZPE'})
        ET.SubElement(zpe, 'scalar', {'units':'kcal/mol'}).text = str(energy)
        
        #add the multiplicity
        mult = ET.SubElement(propertylist, 'property', {'dictRef':'me:spinMultiplicity'})
        ET.SubElement(mult, 'scalar').text = str(reaction.ts.mult)
        
        #add the external symmetry number
        sigma = ET.SubElement(propertylist, 'property', {'dictRef':'me:symmetryNumber'})
        ET.SubElement(sigma, 'scalar').text = str(reaction.ts.sigma_ext / reaction.ts.nopt)
        
        #add the vibrational frequencies
        vibs = ET.SubElement(propertylist, 'property', {'dictRef':'me:vibFreqs'})
        ET.SubElement(vibs, 'array', {'units':'cm-1'}).text = ' '.join([str(fi) for fi in reaction.ts.kinbot_freqs[1:]])
        
        #add the imaginary frequencies
        imfreq = ET.SubElement(propertylist, 'property', {'dictRef':'me:imFreqs'})
        ET.SubElement(imfreq, 'scalar', {'units':'cm-1'}).text = str(-reaction.ts.kinbot_freqs[0])
        
        #add the rotor potentials
        if self.par.par['rotor_scan'] and len(reaction.ts.dihed) > 0:
            qmrotors = ET.SubElement(molecule, 'me:DOSCMethod',{'xsi:type':'me:QMRotors'})
            for i,rot in enumerate(reaction.ts.dihed):
                hinderedrotor = ET.SubElement(qmrotors, 'me:ExtraDOSCMethod',{'xsi:type':'me:HinderedRotorQM1D'})
                ET.SubElement(hinderedrotor, 'me:bondRef').text = bond_ref[ '_'.join(sorted([str(rot[1]),str(rot[2])])) ]
                pot = ET.SubElement(hinderedrotor, 'me:HinderedRotorPotential', {'format':'numerical','units':'kcal/mol','expansionSize':'7','useSineTerms':'yes'})
                ens = reaction.ts.hir.hir_energies[i]
                rotorpot = [(ei - ens[0])*constants.AUtoKCAL for ei in ens]
                for j,ri in enumerate(rotorpot):
                    ET.SubElement(pot, 'me:PotentialPoint', {'angle':str(360/reaction.ts.hir.nrotation*j),'potential':str(ri)})
                ET.SubElement(pot, 'me:PotentialPoint', {'angle':str(360.0),'potential':str(rotorpot[0])})
                ET.SubElement(hinderedrotor, 'me:periodicity').text = str(reaction.ts.sigma_int[rot[1]][rot[2]])
        else:
            qmrotors = ET.SubElement(molecule, 'me:DOSCMethod',{'name':'QMRotors'})
        
        #add the energy transfer model
        etransfer = ET.SubElement(molecule, 'me:energyTransferModel',{'xsi:type':'me:ExponentialDown'})
        ET.SubElement(etransfer, 'me:deltaEDown', {'units':'cm-1'}).text = '-1'

    def run(self):
        """
        write a pbs file for the me/all.xml input file
        submit the pbs file to the queue
        wait for the mess run to finish
        """
        #open the template
        pbs_file = pkg_resources.resource_filename('tpl', 'pbs_mesmer.tpl')
        with open(pbs_file) as f:
            tpl = f.read()
        pbs = open('run_mesmer.pbs','w')
        pbs.write(tpl.format(name = 'mesmer', ppn = self.par.par['ppn'], queue_name = self.par.par['queue_name'], dir = 'me'))
        pbs.close()
        
        command = ['qsub','run_mesmer.pbs']
        process = subprocess.Popen(command,shell=False,stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = process.communicate()
        out = out.decode()
        pid = out.split('\n')[0].split('.')[0]
        
        while 1:
            devnull = open(os.devnull, 'w')
            command = 'qstat -f | grep ' + '"Job Id: ' + pid + '"' + ' > /dev/null'
            if int(subprocess.call(command, shell = True, stdout=devnull, stderr=devnull)) == 0: 
                time.sleep(1)
            else:
                break
        
        return 0
