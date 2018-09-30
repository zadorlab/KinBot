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
from constants import *
from qc import *
from stationary_pt import *
from par import *
import license_message

def createSummaryFile(species,barriers,products):
    fname = 'summary_%s.out'%species.chemid
    f = open(fname,'w+')
    
    f.write(license_message.message)
    f.write('These calculations are done without IRCs,\n')
    f.write('The success thus means that a ts has been found,\n')
    f.write('but does not imply this ts is the correct one!!\n\n')
    f.write('Status\tEnergy\tName\n')
    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] == -1:
            if species.reac_type[index] == 'R_Addition_MultipleBond' and not par.high_level:
                mp2_energy = get_qc_energy(str(species.chemid) + '_well_mp2')[1]
                energy = (barriers[index].energy - mp2_energy) * AUtoKCAL
            else:
                energy = (barriers[index].energy - species.energy) * AUtoKCAL
            prod_name = ''
            name = []
            for prod in products[index]:
                name.append(str(prod.chemid))
            prod_name = ' '.join(sorted(name))
            
            f.write('SUCCESS\t%.2f\t%s\t%s\n'%(energy,species.reac_name[index],prod_name))
        else:
            f.write('FAILED\t\t%s\n'%species.reac_name[index])
            

def createPESViewerInput(species,barriers,products):
    fname = 'pesviewer.inp'
    dir_xyz = 'xyz/'
    if not os.path.exists(dir_xyz):
        os.mkdir(dir_xyz)
    f = open(fname,'w+')
    f.write("> <comments>")
    f.write(license_message.message)
    
    f.write("""This comment is not interpreted, so store any extra info here.
Keywords are case insensitive. Look at the help below.
IMPORTANT: avoid the use of '2d' and '3d' in the names of species, transition states and reactions
(these strings are employed when generating the 2d and 3d files of the molecules)
If you want to use 3D coordinates, store them in a xyz/ directory in the same directory as the python script""")
    f.write('\n\n')
    
    f.write('> <id> %s\n\n'%species.chemid)
    f.write("""> <options> 
units              kcal/mol  #energy units
use_xyz            1         # use xyz, put 0  to switch off
rescale            0         # no rescale , put the well or bimolecular name here to rescale to that value
fh                 9.        # figure height
fw                 18.       # figure width
margin             0.2       # margin fraction on the x and y axis
dpi                120       # dpi of the molecule figures
save               0         # does the plot need to be saved (1) or displayed (0)
write_ts_values    1         # booleans tell if the ts energy values should be written
write_well_values  1         # booleans tell if the well and bimolecular energy values should be written
bimol_color        red       # color of the energy values for the bimolecular products
well_color         blue      # color of the energy values of the wells
ts_color           green     # color or the energy values of the ts, put to 'none' to use same color as line
show_images        1         # boolean tells whether the molecule images should be shown on the graph
rdkit4depict       1         # boolean that specifies which code was used for the 2D depiction""")
    f.write('\n\n')
    
    f.write('> <wells> \n')
    make_xyz(par.atom,species.geom,str(species.chemid),dir_xyz)
    f.write('%s 0.0\n'%(species.chemid)) #use the well as point zero for the energy
    well_energy = species.energy
    wells = [str(species.chemid)]
    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] == -1:
            if len(products[index]) == 1:
                name = str(products[index][0].chemid)
                if not name in wells:
                    make_xyz(par.atom,products[index][0].geom,str(products[index][0].chemid),dir_xyz)
                    energy = (products[index][0].energy - well_energy) * AUtoKCAL
                    f.write('%s %.2f\n'%(products[index][0].chemid,energy))
                    wells.append(name)
    
    f.write('\n')
    
    f.write('> <bimolec> \n')
    bimolecs = []
    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] == -1:
            if len(products[index]) > 1:
                energy = 0. - well_energy
                names = []
                for i in range(len(products[index])):
                    prod = products[index][i]
                    energy += prod.energy
                    names.append(str(prod.chemid))
                    
                name = '_'.join(sorted(names))

                for i in range(len(products[index])):
                    prod = products[index][i]
                    # make twice the same file but with adifferent name ( TODO: is there no better way?)
                    make_xyz(prod.atom,prod.geom,name + str(i+1),dir_xyz) #this is for the pes viewer
                    make_xyz(prod.atom,prod.geom,str(prod.chemid),dir_xyz) #this is for the rmg postprocessing
                energy = energy * AUtoKCAL
                if not name in bimolecs:
                    f.write('%s %.2f\n'%(name,energy))
                    bimolecs.append(name)
    f.write('\n')
    
    f.write('> <ts> \n')
    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] == -1:
            if species.reac_type[index] == 'R_Addition_MultipleBond' and not par.high_level:
                we_energy = get_qc_energy(str(species.chemid) + '_well_mp2')[1]
                energy = (barriers[index].energy - we_energy) * AUtoKCAL
            else:
                energy = (barriers[index].energy - well_energy) * AUtoKCAL
            prod_name = ''
            name = []
            for prod in products[index]:
                name.append(str(prod.chemid))
            prod_name = '_'.join(sorted(name))
            f.write('%s %.2f %s %s\n'%(species.reac_name[index],energy,species.chemid,prod_name))
    f.write('\n')
    
    f.write('> <barrierless> \n\n')
    
    f.write("""> <help>
File follows the rules of SD file format for keywords. Keywords are case
insensitive when parsed.
Keywords:
units: units of the energies supplied above

usexyz: use the xyz coordinates of all the species and render a 2D/3D depiction

rescale: energies are rescaled relative to the energy of the species given here 

wells: all the wells of the PES, separated by lines
each line contains the name, the energy, and optionally the smiles

bimolec: all the bimolecular products of the PES, separated by lines
each line contains the name, the energy, and optionally the smiles of both bimolecular products

ts: all the transition states of the PES, separated by lines
each line contains the name, the energy, and the names of the reactant and product

barrierless: all the barrierless reactions of the PES, separated by lines
each line contains the name and the names of the reactant and product""")
    f.close()

def make_xyz(atoms,geom,name,dir):
    f = open(dir + name + '.xyz', 'w+')
    f.write('%i\n\n'%len(geom))
    
    for index in range(len(geom)):
        f.write('%s %.6f %.6f %.6f\n'%(atoms[index],geom[index][0],geom[index][1],geom[index][2]))
    
    f.write('\n\n')
    f.close()

