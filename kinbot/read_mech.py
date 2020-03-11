import os,sys

def read_mech(mech_file,smi_file):
    """
    Read a chemkin input file and return a list of reactions
    Eeach reactions contains of [[reactant1, reactant2, ...],[product1, product2, ...], [A, n, Ea]]
    The reactants and products are defined based on their smiles
    """
    #make a dictionary {name:smi}
    smis = {}
    smi_lines = open(smi_file,'r').read().split('\n')
    for line in smi_lines:
        if len(line) > 0:
            smis[line.split()[0]] = line.split()[1]
    #read the mech file
    rxns = []
    arrows = ['<=>','=>']
    mech_lines = open(mech_file,'r').read().split('\n')
    for line in mech_lines:
        line = line.split('!')[0]
        pieces = line.split()
        if len(pieces) > 0:
            if not pieces[0].startswith('!'):
                for arrow in arrows:
                    if arrow in pieces[0]:
                        pieces[0] = pieces[0].replace('(+M)','')
                        reactants = get_species(pieces[0].split(arrow)[0].split('+'),smis)
                        products = get_species(pieces[0].split(arrow)[1].split('+'),smis)

                        kinetics = [float(ki) for ki in pieces[1:4]]
                        rxns.append([reactants,products,kinetics])
                        break
    return rxns
    
def get_species(sp,smis):
    """
    from a list of species, take the relevant ones and add them several times
    in case there is a 2X for species C
    """
    third_bodies = ['M']
    species = []
    for s in sp:
        try:
            n = int(s[0])
            s = s[1:]
        except ValueError:
            n = 1
        if not s in third_bodies:
            for ni in range(n):
                species.append(smis[s])
    return species
    
    
if __name__ == '__main__':
    read_mech('/home/rvandev/hept/nc7_16mech.dat.txt','/home/rvandev/hept/smi.txt')
    
