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
This file is part of the KinBot software.

This script reads the reaction mechanism in ~/hept/nc7_16mech.dat.txt,
selects all the hydrogen migration reactions and for each of the reactants
of these reactions, a kinbot run is started. Optionally, the reaction
search is limited to the reactions in the ~/hept/1percent.txt selectivity
file. 
The resulting rate coefficients calculated by KinBot through MESS are finally
compared to the rates in the original mechanism.
"""

import sys,os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

sys.dont_write_bytecode = True

from read_mech import *
from thread_kinbot import *
from cheminfo import *
from stationary_pt import *
import par

def postprocess(dir,jobs,hmigrs,job_names):
    """
    Read the output of all the jobs and put them in a dictionary 
    to compare them to the mechanism reactions
    """
    
    # dictionary with the name of the reactant as key and
    # the inchis of the reactants and products and the Arrhenius parameters are value
    reactions = {}
    
    #iterate all the jobs
    for name in job_names:
        job = job_names[name][0] #name of the calculation
        reactions[name] = []
        
        #iterate the files and look for the summary file to get the chemid of the reactant
        for f in os.listdir(os.path.expanduser(job)):
            if 'summary' in f:
                reac_chemid = f.split('_')[1].split('.')[0]
        #iterate the files and look for the rates.out file to read the reactions and their kinetics
        for f in os.listdir(os.path.expanduser(job)):
            if 'rates.out' in f:
                rate_lines = open(os.path.expanduser(job) + f).read().split('\n')
                #default temperatures
                temps = [300+100*i for i in range(18)]
                if 'temperatures' in rate_lines[0]: #check for custom temperatures
                    temps = [float(ti) for ti in rate_lines[0].split()[2:]] #read the custom temperatures
                for line in rate_lines:
                    if len(line) > 0:
                        if not 'temperatures' in line:
                            pieces = line.split()
                            react_name = pieces[0]
                            #the reactants of the reactions are written by their chemid divided by '_'
                            reactant_chemids = react_name.split('_') 
                            #only keep this reaction is the reactant corresponds to the current job
                            #(rates are written in both direction, so the reverse ones are omitted here)
                            if len(reactant_chemids ) == 1 and reac_chemid in reactant_chemids:
                                #create the inchis of the reactant
                                reactant = [create_inchi(job,chemid) for chemid in react_name.split('_')]
                                prod_name = pieces[1]
                                #create the inchis of the products
                                products = [create_inchi(job,chemid) for chemid in prod_name.split('_')]
                                #give the reaction a name based on the chemids
                                rxn_name = react_name + '_' + prod_name
                                #read the reaction rate coefficient at different temperatures
                                rates = {}
                                for i,ti in enumerate(temps):
                                    if ti >= 500:
                                        rates[ti] = float(pieces[2+i])
                                #add the reaction to the reactions dictionary
                                reactions[name].append([reactant,products,rates,rxn_name])
    
    #compare the kinbot reactions and the mechanism reactions
    compare_reactions(dir,reactions,hmigrs,job_names)
    
def compare_reactions(dir,kb_reactions,hmigrs,job_names):
    """
    Compare the mechanism and kb reactions
    """
    
    #output file to write the summary to
    f_iden_out = open('%s/identical_rxns.out'%dir,'w+')
    
    #output file with the arrhenius parameters and chemkin names of the species
    #to directly be implemented in a chemkin input file
    f_out_chemkin = open('%s/chemkin.inp'%dir,'w+')
    #chemkin names of the species
    chemkin_names = read_chemkin_names(dir)
    
    #convert the hmigration reaction of the mechanism to an inchi representation for
    #easy comparison 
    hmigrs_inchi = []
    for hmigr in hmigrs:
        re_inchi = [create_inchi_from_smi(re_smi) for re_smi in hmigr[0]]
        pr_inchi = [create_inchi_from_smi(pr_smi) for pr_smi in hmigr[1]]
        hmigrs_inchi.append([re_inchi,pr_inchi,hmigr[2]])
    
    #iterate the reactions in the mechanism
    for i, hmigr in enumerate(hmigrs_inchi):
        new = 1
        for name in job_names:
            if job_names[name][2] == hmigr[0][0]:
                for index,kb_rxn in enumerate(kb_reactions[name]):
                    #check if a kinbot reaction and a mechanism reaction are identical
                    if equal_rxns(kb_rxn[:2],hmigr[:2]):
                        new = 0
                        #if identical, regress the rate coefficients into a modified Arrhenius form
                        kinbot_arrh = regress_arrh(kb_rxn)
                        #make a figure consisting of the reaction depiction and the kinetics from kinbot 
                        #and from the original mechanism
                        make_fig(dir,name,index,kb_rxn,kinbot_arrh,hmigr, '.'.join(hmigrs[i][0]), '.'.join(hmigrs[i][1]))
                        break
                break
        #write the reaction to a summary file
        if new:
            line = '{}\t{}\t{}\t'.format(name,''.join(hmigrs[i][0]),''.join(hmigrs[i][1]))
            line += '{}\t{}\t'.format(0,1)
            line += '{:.2E}\t{:.2f}\t{:.2f}\t{:.2E}\t{:.2f}\t{:.2f}\n'.format(hmigr[2][0],hmigr[2][1],hmigr[2][2],0.,0.,0.)
        else:
            line = '{}\t{}\t{}\t'.format(name,''.join(hmigrs[i][0]),''.join(hmigrs[i][1]))
            line += '{}\t{}\t'.format(1,0)
            line += '{:.2E}\t{:.2f}\t{:.2f}\t{:.2E}\t{:.2f}\t{:.2f}\n'.format(hmigr[2][0],hmigr[2][1],hmigr[2][2],kinbot_arrh[0],kinbot_arrh[1],kinbot_arrh[2])
            
            re_names = get_names(chemkin_names,hmigr[0])
            pr_names = get_names(chemkin_names,hmigr[1])
            chemkin_line = '{}<=>{}\t'.format('+'.join(re_names),'+'.join(pr_names))
            chemkin_line += '{:.2E}\t{:.2f}\t{:.2f}\n'.format(kinbot_arrh[0],kinbot_arrh[1],kinbot_arrh[2])
            f_out_chemkin.write(chemkin_line)
        f_iden_out.write(line)
        
    f_out_chemkin.close()
    f_iden_out.close()

def get_names(chemkin_names,inchis):
    """
    Get the names of the species from the original mechanism
    """
    names = ['' for inchi in inchis]
    
    for chemkin_name in chemkin_names:
        chemkin_inchi = chemkin_names[chemkin_name]
        for i,inchi in enumerate(inchis):
            if inchi == chemkin_inchi:
                names[i] = chemkin_name
    return names
    
def read_chemkin_names(dir):
    """
    Read the chemkin names from the dir/smi.txt file.
    """
    chemkin_names = {}
    lines = open('{}/smi.txt'.format(dir),'r').readlines()
    for line in lines:
        if len(line) > 0:
            pieces = line.split()
            name = pieces[0]
            smi = pieces[1]
            smi = smi.replace('singlet','')
            smi = smi.replace('triplet','')
            chemkin_inchi = create_inchi_from_smi(smi)
            chemkin_names[name] = chemkin_inchi
    return chemkin_names

def regress_arrh(rxn):
    """
    Do a linear least square regression by writing the modified Arrhenius formula
    in linear form (natural log of both sides)
    """
    temps = sorted(rxn[2].keys())
    temps_R_inv = [1./(8.314*ti) for ti in temps]
    temps_ln = [np.log(ti) for ti in temps]
    rates = [rxn[2][ti] for ti in temps]
    rates_ln = np.array([np.log(ki) for ki in rates])
    
    B = []
    for i in range(len(temps)):
        B.append([1, temps_ln[i], temps_R_inv[i]])
        
    B = np.array(B)
    lnA,n,Ea = np.linalg.lstsq(B, rates_ln)[0]
    
    A = np.power(np.e,lnA)
    Ea = - Ea / 4.184
    
    return A, n, Ea
    
    
def make_fig(dir,species_name,index,kb_rxn,kinbot_arrh,hmigr, react_smiles, prod_smiles):
    """
    For each reaction that is found by KinBot and is in the original mechanism,
    a figure is made containing the reaction depiction and the kinetics both
    from KinBot as well as from the original mechanism.
    """
    #folder where all the plots are saved
    fig_dir = dir + 'figures/'
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)
    #create a depiction of the reaction based on the smiles
    if not os.path.exists(fig_dir + kb_rxn[3] + '.png'):
        create_rxn_depiction(react_smiles, prod_smiles,fig_dir,kb_rxn[3])
    #x-axis contains 1000K/T
    temps = sorted(kb_rxn[2].keys())
    temps_inv = [1000.0/ti for ti in temps]
    
    #The kinbot rates are plotted twice, ones with the output of MESS and ones with the regressed Arrhenius parameters
    kinbot_rates = [kb_rxn[2][ti] for ti in temps]
    kb_rates_log = [np.log10(ki) for ki in kinbot_rates]
    kb_rates_arrh = [kinbot_arrh[0]*np.power(1000./ti,kinbot_arrh[1])*np.exp(-kinbot_arrh[2]*4.184/8.314/(1000./ti)) for ti in temps_inv]
    kb_rates_arrh = [np.log10(ki) for ki in kb_rates_arrh]
    
    #The mechanism rate coefficients are plotted based on their Arrhenius parameters
    mech_rates = hmigr[2]
    mech_rates_arrh = [mech_rates[0]*np.power(1000./ti,mech_rates[1])*np.exp(-mech_rates[2]*4.184/8.314/(1000./ti)) for ti in temps_inv]
    mech_rates_arrh = [np.log10(ki) for ki in mech_rates_arrh]
    
    #position the depiction of the reaction in the bottom left corner
    img = mpimg.imread(fig_dir + kb_rxn[3] + '.png')
    imy = len(img) + 0.
    imx = len(img[0]) + 0.
    x_min = temps_inv[-1]
    if kb_rates_log[0] > -60 and kb_rates_log[0] < 60:
        y_min = min(mech_rates_arrh[-1],kb_rates_log[0])
    else:
        y_min = mech_rates_arrh[-1]
    if kb_rates_log[-1] > -60 and kb_rates_log[-1] < 60:
        y_max = max(mech_rates_arrh[0],kb_rates_log[-1])
    else:
        y_max = mech_rates_arrh[0]
    imw = (temps_inv[0]-x_min+0.)/(8.0)*imx/120
    imh = (y_max-y_min)/(6.0)*imy/120
    fig, ax = plt.subplots()
    extent=(x_min , x_min+imw, y_min, y_min+imh)
    ax.imshow(img, aspect='auto', extent=extent, zorder=-1)
        
    #plot the kinbot rate coefficients calculated through MESS
    plt.plot(temps_inv,kb_rates_log,'bo', label = 'KinBot')
    #plot the kinbot rate cofficients in their Arrhenius form
    plt.plot(temps_inv,kb_rates_arrh,color='blue', ls = 'solid')
    #plot the mechanism rate coefficients
    plt.plot(temps_inv,mech_rates_arrh,color='green', ls = 'dashed', label = 'Mechanism')
    
    plt.ylabel('Log(k / 1s)')
    plt.xlabel('1000K/T')
    
    plt.legend()
    
    plt.savefig(fig_dir + species_name + '_nr' + str(index) + '_plot.png',bbox_inches='tight')
    plt.close('all')


def equal_rxns(rxn1,rxn2):
    """
    Compare two reactions based on the inchis of the reactants and products
    Only do this in the forward direction
    """
    return equal_species(rxn1[0],rxn2[0]) and equal_species(rxn1[1],rxn2[1])

def equal_species(list1,list2):
    """
    Compare two species lists, this method only works if the length of the list is maximum 2
    """
    if len(list1) == len(list2):
        if len(list1) == 1:
            return list1[0] in list2[0] or list2[0] in list1[0]
        elif len(list1) == 2: 
            bool1 = list1[0] in list2[0] or list2[0] in list1[0]
            bool2 = list1[1] in list2[1] or list2[1] in list1[1]
            bool3 = list1[0] in list2[1] or list2[0] in list1[1]
            bool4 = list1[1] in list2[0] or list2[1] in list1[0]

            return (bool1 and bool2) or (bool3 and bool4)
        else:
            #cannot handle 3 or more products at the moment
            return 0
    else:
        return 0

def generate_name(smi,names):
    """
    Method to generate a name based on the molecular formula of the molecule, 
    which is unique, i.e. not in the names values
    """
    molform = get_molecular_formula(smi)
    
    name = molform
    it = 1
    while name in names:
        name = '{}_{}'.format(molform, it)
        it += 1
    return name

def preprocess(rxns):
    """
    Select the H-migration reactions for the reaction list
    and return the list of reactants and the list of H-migration reactions
    """
    hmigrs = []
    species = [] #list of smiles for which H migrations were found
    for rxn in rxns:
        if len(rxn[0]) == 1 and len(rxn[1]) == 1:
            obmol = create_ob_mol(rxn[0][0])
            if is_H_migration(rxn[0][0], rxn[0][0],obmol.spin):
                hmigrs.append(rxn)
                if not rxn[0][0] in species:
                    species.append(rxn[0][0])
    return species, hmigrs

def is_H_migration(reactant,product,mult):
    """
    Method to check if a reaction is a H-migration
    """
    # make the parent chemids of all the resonance isomers of the reacant
    react_ids = saturate(reactant,mult)
    
    # make the parent chemids of all the resonance isomers of the product
    prod_ids = saturate(product,mult)
    
    #check if at least one parent molecule (considering all resonance isomers) is identical
    for rid in react_ids:
        if rid in prod_ids:
            return 1

    return 0

def saturate(smi,mult):
    """
    Saturate each of the radical resonance isomer to its corresponding parent molecule
    by adding a H atom at the radical site
    
    Return the chemids of the parent molecules
    """
    par.mult = mult
    par.charge = 0
    obmol, structure = generate_3d_structure(smi)
    par.natom = len(obmol.atoms)
    
    structure = np.reshape(structure, (par.natom,4))
    par.atom = structure[:,0]
    geom = structure[:,1:4].astype(float)
    
    st_pt = stationary_pt('temp')
    st_pt.geom = geom
    st_pt.bond_mx(par.natom, par.atom)
    st_pt.characterize(par.natom, par.atom, par.mult, par.charge)
    
    par.natom += 1
    par.atom = np.append(par.atom,np.array(['H']))

    chemids = []
    for i, res in enumerate(st_pt.rads):
        if not 1 in res:
            return []
        b = np.zeros((par.natom,par.natom), dtype=int)
        for j in range(len(st_pt.bonds[i])):
            for k in range(len(st_pt.bonds[i])):
                b[j][k] = st_pt.bonds[i][j][k]
        r = np.where(res == 1)[0][0]
        b[par.natom-1][r] = 1
        b[r][par.natom-1] = 1
        st_pt.bond = b
        st_pt.calc_chemid(par.natom, par.atom, par.mult)
        chemids.append(st_pt.chemid)
    return chemids

def write_hmigrs(file, rxns):
    """
    This method writes a set of rxns to a separate file as such
    that the reactions do not have to be read and filtered each time
    this script is run
    """
    f_out = open(file,'w')
    for rxn in rxns:
        s = '+'.join(rxn[0])
        s += '<=>'
        s += '+'.join(rxn[1])
        s += ' '
        s += ' '.join([str(ki) for ki in rxn[2]])
        s += '\n'
        f_out.write(s)
    f_out.close()

def read_hmigrs(file):
    """
    This method reads the file and returns the reactions in the file as such
    that the reactions do not have to be read and filtered each time
    this script is run
    """
    rxns = []
    reactants = []
    f = open(file,'r')
    lines = f.read().split('\n')
    for line in lines:
        if len(line) > 0:
            pieces = line.split()
            reacts = pieces[0].split('<=>')[0].split('+')
            prods = pieces[0].split('<=>')[1].split('+')
            kin = [float(ki) for ki in pieces[1:]]
            rxns.append([reacts,prods,kin])
            for re in reacts:
                if not re in reactants:
                    reactants.append(re)
    return reactants, rxns

def read_sens(dir,sp,hmigrs):
    """
    Read the sensitivities of the mechanism in the ~/hept/1percent.txt file
    And return the reactions that correspond to the H-migrations
    """
    #convert the H-migration reactions to their inchi representation
    hmigrs_inchi = []
    for hmigr in hmigrs:
        re_inchi = [create_inchi_from_smi(re_smi) for re_smi in hmigr[0]]
        pr_inchi = [create_inchi_from_smi(pr_smi) for pr_smi in hmigr[1]]
        hmigrs_inchi.append([re_inchi,pr_inchi,hmigr[2]])
    
    #read the sensitivity file
    sens_file = '1percent.txt'
    sens_lines = open(dir + sens_file,'r').readlines()
    
    #for each reaciton in the sensitivity file, check if it is a H-migration
    #and add it to the new list of H-migration reactions
    sp_2 = []
    hmigrs_2 = []
    for line in sens_lines:
        rxn = line.split()[1]
        re_smi = rxn.split('>>')[0].split('.')
        re_smi = [smi.split('_')[0] for smi in re_smi]
        pr_smi = rxn.split('>>')[1].split('.')
        pr_smi = [smi.split('_')[0] for smi in pr_smi]
        if len(re_smi) ==1 and len(pr_smi) == 1:
            re_inchi = [create_inchi_from_smi(smi) for smi in re_smi]
            pr_inchi = [create_inchi_from_smi(smi) for smi in pr_smi]
            for index,hmigr in enumerate(hmigrs_inchi):
                if equal_rxns([re_inchi,pr_inchi],hmigr[:2]):
                    re_smiles = hmigrs[index][0][0]
                    hmigrs_2.append(hmigrs[index])
                    if not re_smiles in sp_2:
                        sp_2.append(re_smiles)
    return sp_2, hmigrs_2


def main(rmg_dir,use_sensitivity = 1):
    """
    main driver for the calculations of the heptane mechanism H-migration reactions
    The rmg_dir is the absolute path ~/hept/
    
    This script:
    
    1. Read the heptane mechanism
    2. Filters the mechanism to find all H-migration reactions
    3. Optionally filters these reaction depending on their sensitivity
    4. For each reactant, start a KinBot calculation only searching for H-migration reactions
    5. Postprocess the results by comparing the KinBot rate coefficients to the ones in the original mechanism
    """
    #add the KinBot code tot he path
    sys.path.append(os.path.expanduser('~/KinBot/code'))
    
    #directory where all the calculations are done (absolute path)
    dir = os.path.expanduser(rmg_dir)
    
    #read the mechanism (or the summary file if present)
    if not os.path.exists(dir + 'hmigr_rxns.txt'):
        rxns = read_mech(dir + 'nc7_16mech.dat.txt', dir + 'smi.txt')
        sp, hmigrs = preprocess(rxns)
        write_hmigrs(dir + 'hmigr_rxns.txt', hmigrs)
    else:
        sp, hmigrs = read_hmigrs(dir + 'hmigr_rxns.txt')
    
    
    
    # whether to use the sensitivity information
    if use_sensitivity:
        sp_2,hmigrs_2 = read_sens(dir,sp,hmigrs)
    else: 
        sp_2 = sp
        hmigrs_2  = hmigrs
    
    names = [] #all the names of the reactants 
    jobs = [] #all the jobs to be done
    job_names = {} #dictionary coupling the names to the jobs, the smiles and the inchis
    
    #iterate the reactants and write the job input file
    for smi in sp:
        name = generate_name(smi,names) #generate the name of the reactant
        names.append(name)
        obmol = create_ob_mol(smi)
        mult = obmol.spin
        
        #make the job directory
        sp_dir = '%s/%s/'%(dir,name)
        if not os.path.exists(sp_dir):
            os.makedirs(sp_dir)
        
        syms = {1:'H',6:'C',8:'O'}
        
        # create a kinbot input file
        header = 'title=\'%s\'\n\n'%name
        header += "charge = 0\nmult = %i\nmethod = 'b3lyp'\nbasis = '6-31g'\n"%mult
        header += "qc = 'gauss'\nconformer_search = 1\nreaction_search = 1\nbarrier_threshold = 200.\n"
        header += "families = ['intra_H_migration']\nme = 1\n" #only consider H-migration reactions
        header += "rotor_scan = 1\nhigh_level = 1\nhigh_level_method = 'WB97XD'\nhigh_level_basis = '6-31+G*'\n"
        header += "ga = 0\nngen = 0\nppn = 8\nscratch = '/scratch/rvandev'\n\n"
        s = header
        s += 'natom = %i\n\n'%len(obmol.atoms)
        s += "smiles = '%s'"%smi

        open('%sinput.inp'%(sp_dir),'w+').write(s)

        #name of the job
        job = '%s%s/'%(rmg_dir,name)
        
        #optionally filter jobs based on the sensitivity
        #(sp_2 and sp are equal if no sensitivity is used)
        if smi in sp_2:
            jobs.append(job)
            inchi = create_inchi_from_smi(smi)
            job_names[name] = [job,smi,inchi]
        
        #depict the heptane mechanism reactions for this species 
        create_depics = 0
        if create_depics:
            rxn_names = []
            for hmigr in hmigrs:
                re_smi = hmigr[0][0]
                if re_smi == si:
                    pr_smi = hmigr[1][0]
                    i = 1
                    rxn_name = 'rmg_' + name + '_nr' + str(i)
                    while rxn_name in rxn_names:
                        rxn_name = 'rmg_' + name + '_nr' + str(i)
                        i += 1
                    rxn_names.append(rxn_name)
                    create_rxn_depiction(re_smi,pr_smi,sp_dir, rxn_name)
    
    #run the KinBot jobs with the thread_kinbot script,
    #with a maximum of 5 jobs running simultaneously
    #run_threads(jobs, 'heptane', max_running = 5)
    
    #postprocess the results
    postprocess(dir,jobs,hmigrs_2,job_names)
    
    print 'Done!'


if __name__ == "__main__":
    rmg_dir = '~/hept/'
    main(rmg_dir)
    
    