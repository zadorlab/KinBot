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
import sys,os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

sys.dont_write_bytecode = True

from thread_kinbot import *
from cheminfo import *
from stationary_pt import *
import par

def postprocess(dir,jobs,job_names):
    # dictionary with the inchis of the reactants and products (for comparison to rmg)
    reactions = {}
    
    for name in job_names: 
        job = job_names[name][0]
        reactions[name] = []

        reac_chemid = ''
        for f in os.listdir(os.path.expanduser(job)):
            if 'summary' in f:
                reac_chemid = f.split('_')[1].split('.')[0]
        for f in os.listdir(os.path.expanduser(job)):
            if 'rates.out' in f:
                rate_lines = open(os.path.expanduser(job) + f).read().split('\n')
                temps = [300+100*i for i in range(18)]
                if 'temperatures' in rate_lines[0]:
                    temps = [float(ti) for ti in rate_lines[0].split()[1:]]
                for line in rate_lines:
                    if len(line) > 0:
                        if not 'temperatures' in line:
                            pieces = line.split()
                            react_name = pieces[0]
                            reactant_chemids = react_name.split('_')
                            if len(reactant_chemids ) == 1 and reac_chemid in reactant_chemids:
                                reactant = [create_inchi(job,chemid) for chemid in react_name.split('_')]
                                prod_name = pieces[1]
                                products = [create_inchi(job,chemid) for chemid in prod_name.split('_')]
                                rxn_name = react_name + '_' + prod_name
                                #rates = [float(xi) for xi in pieces[2:]]
                                rates = {}
                                for i,ti in enumerate(temps):
                                    rates[ti] = float(pieces[2+i])
                                reactions[name].append([reactant,products,rates,rxn_name])
    
    compare_reactions(dir,reactions,job_names)
    
def compare_reactions(dir,kb_reactions,job_names):
    rmg_reactions = {}
    # get the inchis of the reactants and products of each reaction
    with open('%s/rxnqueue.json'%dir) as json_file:
        json_data = json.load(json_file)
        for rxn in json_data:
            if rxn['Family'] == 'intra_H_migration':
                name = rxn['name']
                reactants = [sp['InChi'] for sp in rxn['Reactants']]
                products = [sp['InChi'] for sp in rxn['Products']]
                kfor = {}
                krev = {}
                for i in range(18):
                    ti = 300+100*i
                    kfor[ti] = rxn['k'+str(ti)]
                    krev[ti] = rxn['krev'+str(ti)]
                rmg_reactions[name] = [reactants,products,kfor,krev]
    
    #output file to write the summary to
    f_out = open('%s/identical_rxns.out'%dir,'w+')
    f_out_chemkin = open('%s/chemkin.inp'%dir,'w+')
    
    for i, hmigr_name in enumerate(rmg_reactions):
        hmigr = rmg_reactions[hmigr_name]
        rmg_arrh = regress_arrh(hmigr)
        react_smiles = [create_smiles(inchi) for inchi in hmigr[0]]
        prod_smiles = [create_smiles(inchi) for inchi in hmigr[1]]
        new = 1
        for name in job_names:
            if job_names[name][1] == hmigr[0][0]: #compare the inchis
                for index,kb_rxn in enumerate(kb_reactions[name]):
                    if equal_rxns(kb_rxn[:2],hmigr[:2]):
                        new = 0
                        kinbot_arrh = regress_arrh(kb_rxn)
                        make_fig(dir,name,index,kb_rxn,kinbot_arrh,hmigr,rmg_arrh, '.'.join(react_smiles), '.'.join(prod_smiles))
                        break
                break
        if new:
            line = '{}\t{}\t{}\t'.format(name,''.join(react_smiles),''.join(prod_smiles))
            line += '{}\t{}\t'.format(0,1)
            line += '{:.2E}\t{:.2f}\t{:.2f}\t{:.2E}\t{:.2f}\t{:.2f}\n'.format(rmg_arrh[0],rmg_arrh[1],rmg_arrh[2],0.,0.,0.)
        else:
            line = '{}\t{}\t{}\t'.format(name,''.join(react_smiles),''.join(prod_smiles))
            line += '{}\t{}\t'.format(1,0)
            line += '{:.2E}\t{:.2f}\t{:.2f}\t{:.2E}\t{:.2f}\t{:.2f}\n'.format(rmg_arrh[0],rmg_arrh[1],rmg_arrh[2],kinbot_arrh[0],kinbot_arrh[1],kinbot_arrh[2])
            
            re_names = get_names(dir,hmigr[0])
            pr_names = get_names(dir,hmigr[1])
            chemkin_line = '{}<=>{}\t'.format('+'.join(re_names),'+'.join(pr_names))
            chemkin_line += '{:.2E}\t{:.2f}\t{:.2f}\n'.format(kinbot_arrh[0],kinbot_arrh[1],kinbot_arrh[2])
            f_out_chemkin.write(chemkin_line)
        f_out.write(line)

    f_out_chemkin.close()
    f_out.close()

def get_names(dir,inchis):
    """
    Get the names of the species from the original mechanism
    """
    names = ['' for inchi in inchis]
    
    with open('%s/queue.json'%dir) as json_file:
        json_data = json.load(json_file)
        for molecule in json_data:
            for i,inchi in enumerate(inchis):
                if inchi == molecule['InChi']:
                    names[i] = molecule['name']
    return names

def regress_arrh(rxn):
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
    
    
def make_fig(dir,name,index,kb_rxn,kinbot_arrh,hmigr,rmg_arrh, react_smiles, prod_smiles):
    fig_dir = dir + 'figures/'
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)
    
    if not os.path.exists(fig_dir + kb_rxn[3] + '.png'):
        create_rxn_depiction(react_smiles, prod_smiles,fig_dir,kb_rxn[3])

    temps = sorted(kb_rxn[2].keys())
    temps_inv = [1000.0/ti for ti in temps]
    
    kb_rates = [kb_rxn[2][ti] for ti in temps]
    kb_rates_log = [np.log10(ki) for ki in kb_rates]
    kb_rates_arrh = [kinbot_arrh[0]*np.power(1000./ti,kinbot_arrh[1])*np.exp(-kinbot_arrh[2]*4.184/8.314/(1000./ti)) for ti in temps_inv]
    kb_rates_arrh = [np.log10(ki) for ki in kb_rates_arrh]
    
    rmg_rates = [hmigr[2][ti] for ti in temps]
    rmg_rates_log = [np.log10(ki) for ki in rmg_rates]
    rmg_rates_arrh = [rmg_arrh[0]*np.power(1000./ti,rmg_arrh[1])*np.exp(-rmg_arrh[2]*4.184/8.314/(1000./ti)) for ti in temps_inv]
    rmg_rates_arrh = [np.log10(ki) for ki in rmg_rates_arrh]
    
    img = mpimg.imread(fig_dir + kb_rxn[3] + '.png')
    imy = len(img) + 0.
    imx = len(img[0]) + 0.
    x_min = temps_inv[-1]
    if kb_rates_log[0] > -60 and kb_rates_log[0] < 60:
        y_min = min(rmg_rates_log[0],kb_rates_log[0])
    else:
        y_min = rmg_rates_log[0]
    if kb_rates_log[-1] > -60 and kb_rates_log[-1] < 60:
        y_max = max(rmg_rates_log[-1],kb_rates_log[-1])
    else:
        y_max = rmg_rates_log[-1]
    figw = 8. #default: 8
    figh = 6. #default: 6
    #dpi_fig = 120 #default: 100
    dpi_im = 120
    #font = 6
    imw = (temps_inv[0]-x_min+0.)/(figw)*imx/dpi_im
    imh = (y_max-y_min)/(figh)*imy/dpi_im
    
    #plt.rcParams["figure.figsize"]=[figw,figh]
    #plt.rcParams["figure.dpi"]=dpi_fig
    #plt.rcParams["font.size"]=font
    fig, ax = plt.subplots()
    
    extent=(x_min , x_min+imw, y_min, y_min+imh)
    ax.imshow(img, aspect='auto', extent=extent, zorder=-1)
    
    plt.plot(temps_inv,kb_rates_log,'bo', label = 'KinBot')
    plt.plot(temps_inv,rmg_rates_log,'gx', label = 'RMG')
    
    
    plt.plot(temps_inv,kb_rates_arrh,color='blue', ls = 'solid')
    plt.plot(temps_inv,rmg_rates_arrh,color='green', ls = 'dashed')
    
    plt.ylabel('Log(k / 1s)')
    plt.xlabel('1000K/T')
    
    plt.legend()
    
    plt.savefig(fig_dir + name + '_nr' + str(index) + '_plot.png',bbox_inches='tight')
    plt.close('all')
    
  

def equal_rxns(rxn1,rxn2):
    bool1 = equal_species(rxn1[0],rxn2[0]) and equal_species(rxn1[1],rxn2[1])
    #bool2 = equal_species(rxn1[0],rxn2[1]) and equal_species(rxn1[1],rxn2[0])
    
    #return bool1 or bool2
    return bool1

def equal_species(list1,list2):
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
    while name in names.values():
        name = '{}_{}'.format(molform, it)
        it += 1
    return name

def preprocess(dir):
    """
    Read the reaction queue file and extract the number of reactions per family. 
    For the h-migration reactions, search for the ones that include a peroxide
    """
    species = []
    names = []
    with open('%s/rxnqueue.json'%dir) as json_file:
        json_data = json.load(json_file)
        for rxn in json_data:
            if rxn['Family'] == 'intra_H_migration':
                names.append(rxn['name'])
                species.append(rxn['Reactants'][0]['InChi'])
    return species

def main(rmg_dir):
    dir = os.path.expanduser(rmg_dir)
    sp = preprocess(dir)
    inchis = []
    species = {}
    names = {}
    multi = {}
    syms = {1:'H',6:'C',8:'O'}
    jobs = []
    sys.path.append(os.path.expanduser('~/KinBot/code'))
    dir = os.path.expanduser(rmg_dir)
    with open('%s/queue.json'%dir) as json_file:
        json_data = json.load(json_file)
        for molecule in json_data:
            smi = ''
            inchi = ''
            name = ''
            mult = -1
            for mol_attrib in molecule: 
                if mol_attrib == 'SMILES':
                    smi = molecule[mol_attrib][0] # take the first one of the resonance isomers
                """
                do NOT use the rmg names: 
                KinBot uses the names as directory names for the calculations, but they
                contain brackets and other symbols which are not well handled by Linux
                when used in file and directory names
                if mol_attrib == 'name':
                    name = molecule[mol_attrib]
                """
                if mol_attrib == 'multiplicity':
                    mult = molecule[mol_attrib]
                if mol_attrib == 'InChi':
                    inchi = molecule[mol_attrib]
            if inchi in sp:
                inchis.append(inchi)
                species[inchi] = smi
                names[inchi] = generate_name(smi,names)
                multi[inchi] = mult
                
    job_names = {} #dictionary coupling the names to the jobs and the inchis
    
    for inchi in inchis:
        sp_dir = '%s/%s/'%(dir,names[inchi])
        if not os.path.exists(sp_dir):
            os.makedirs(sp_dir)
        
        syms = {1:'H',6:'C',8:'O'}
        
        m = create_ob_mol(str(species[inchi]))
        elem = 1 # checks if all the elements are either C, O or H
        for at in m.atoms:
            if not at.atomicnum in syms:
                elem = 0
        if elem:
            # create a kinbot input file
            header = 'title=\'%s\'\n\n'%names[inchi]
            header += "charge = 0\nmult = %i\nmethod = 'b3lyp'\nbasis = '6-31g'\n"%multi[inchi]
            header += "qc = 'gauss'\nconformer_search = 1\nreaction_search = 1\nbarrier_threshold = 200.\n"
            header += "families = ['intra_H_migration']\nme = 1\n"
            header += "rotor_scan = 1\nhigh_level = 1\nhigh_level_method = 'WB97XD'\nhigh_level_basis = '6-31+G*'\n"
            header += "ga = 0\nngen = 0\nppn = 8\nscratch = '/scratch/rvandev'\n\n"
            s = header
            s += 'natom = %i\n\n'%len(m.atoms)
            s += "smiles = '%s'"%str(species[inchi])

            open('%sinput.inp'%(sp_dir),'w+').write(s)

            job = '%s%s/'%(rmg_dir,names[inchi])
            
            jobs.append(job)
            job_names[names[inchi]] = [job,inchi]
            
            #depict the rmg rxn for this species 
            create_rmg_depics = 1
            if create_rmg_depics:
                rxn_names = []
                with open('%s/rxnqueue.json'%dir) as json_file:
                    json_data = json.load(json_file)
                    for rxn in json_data:
                        if rxn['Family'] == 'intra_H_migration':
                            if rxn['Reactants'][0]['InChi'] == inchi:
                                re_smi = rxn['Reactants'][0]['SMILES'][0]
                                pr_smi = rxn['Products'][0]['SMILES'][0]
                                #rxn_name = 'rmg_' + names[inchi]
                                i = 1
                                rxn_name = 'rmg_' + names[inchi] + '_nr' + str(i)
                                while rxn_name in rxn_names:
                                    i += 1
                                    rxn_name = 'rmg_' + names[inchi] + '_nr' + str(i)
                                rxn_names.append(rxn_name)
                                create_rxn_depiction(re_smi,pr_smi,dir + names[inchi], rxn_name)
    run_threads(jobs, 'butane', max_running = 10)

    postprocess(dir,jobs,job_names)
    
    print 'Done!'


if __name__ == "__main__":
    rmg_dir = '~/naturalgas/'
    main(rmg_dir)
    
    