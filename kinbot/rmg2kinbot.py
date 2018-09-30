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

sys.dont_write_bytecode = True

from thread_kinbot import *
from cheminfo import *


def postprocess(dir,jobs):

    # file with summary of all succeeded and failed reaction searches
    f_out = open('%s/summary.txt'%dir,'w+')
    f_out.write('Job name\tTotal rxns\tSuccessful rxns\n')
    
    # dictionary with the inchis of the reactants and products (for comparison to rmg)
    reactions = {}
    
    for job in jobs: 
        rxn_count = 0
        success = {}
        failed = []
        for f in os.listdir(os.path.expanduser(job)):
            if 'summary' in f:
                summary_lines = open(os.path.expanduser(job) + f).read().split('\n')[5:]
                for line in summary_lines:
                    if len(line) > 0:
                        rxn_count += 1
                        if 'SUCCESS' == line.split()[0]:
                            pieces = line.split()
                            name = pieces[2]
                            success[name] = float(pieces[1])
                            reactant = [create_inchi(job,name.split('_')[0])]
                            products = [create_inchi(job,chemid) for chemid in pieces[3:]]
                            reactions[name] = [reactant,products,float(pieces[1])]
                        else:
                            failed.append(line.split()[1])
        f_out.write('%s\t%i\t%i\n'%(job,rxn_count,len(success)))
    
    compare_reactions(dir,reactions)
    
def compare_reactions(dir,kb_reactions):
    rmg_reactions = {}
    # get the inchis of the reactants and products of each reaction
    with open('%s/rxnqueue.json'%dir) as json_file:
        json_data = json.load(json_file)
        for rxn in json_data:
            name = rxn['name']
            reactants = [sp['InChi'] for sp in rxn['Reactants']]
            products = [sp['InChi'] for sp in rxn['Products']]
            if len(reactants) == 1 or len(products) == 1: #only keep monomolecular reactions (at least in one direction)
                rmg_reactions[name] = [reactants,products]
    
    #compare the rmg and kb reactions
    f_iden_out = open('%s/identical_rxns.out'%dir,'w+')
    f_new_out = open('%s/new_rxns.out'%dir,'w+')
    
    for kb_rxn in kb_reactions:
        react_smiles = [create_smiles(smi) for smi in kb_reactions[kb_rxn][0]]
        prod_smiles = [create_smiles(smi) for smi in kb_reactions[kb_rxn][1]]
        line = '%.2f\t%s\t%s\t%s\n'%(kb_reactions[kb_rxn][2],kb_rxn,'\t'.join(react_smiles),'\t'.join(prod_smiles))
        new = 1
        for rmg_rxn in rmg_reactions:
            if equal_rxns(kb_reactions[kb_rxn][:2],rmg_reactions[rmg_rxn]):
                new = 0
        if new:
            f_new_out.write(line)
        else:
            f_iden_out.write(line)
    
    f_iden_out.close()
    f_new_out.close()
    

def equal_rxns(rxn1,rxn2):
    bool1 = equal_species(rxn1[0],rxn2[0]) and equal_species(rxn1[1],rxn2[1])
    bool2 = equal_species(rxn1[0],rxn2[1]) and equal_species(rxn1[1],rxn2[0])
    
    return bool1 or bool2

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
    
def main(rmg_dir):
    inchis = []
    species = {}
    names = {}
    multi = {}
    syms = {1:'H',6:'C',8:'O'}
    jobs = []
    inc = ''
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
            inchis.append(inchi)
            species[inchi] = smi
            names[inchi] = generate_name(smi,names)
            multi[inchi] = mult
    
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
            header += "qc = 'gauss'\nconformer_search = 0\nreaction_search = 1\nbarrier_threshold = 200.\n"
            header += "ga = 0\nngen = 0\nppn = 8\nscratch = '/scratch/rvandev'\n\n"
            s = header
            s += 'natom = %i\n\n'%len(m.atoms)
            s += "smiles = '%s'"%str(species[inchi])

            open('%sinput.inp'%(sp_dir),'w+').write(s)

            job = '%s%s/'%(rmg_dir,names[inchi])

            jobs.append(job)
    
    run_threads(jobs, 'rmg', max_running = 20)

    postprocess(dir,jobs)


if __name__ == "__main__":
    rmg_dir = '~/rmg/'
    main(rmg_dir)
    #postprocess(os.path.expanduser(rmg_dir),['~/rmg/CH3CHOOCHO/'])
    
    