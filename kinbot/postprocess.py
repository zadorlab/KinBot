"""
This file contains the postprocessing of the KinBot run
It includes
1. Writing a summary file with all successful and failed
reactions, including their barrier heights and which products are formed

2. Writing an input file for the PES viewer.
"""
import os
import numpy as np

from kinbot import kb_path
from kinbot import license_message
from kinbot import constants


def creatMLInput(species, qc, par):
    """
    Create the files for the input of the Machine Learning tool
    of Ghent University.
    This input consists of the atom vectors and bond matrices of
    the reactant and all the products, and the transition state
    bond distances
    """
    directory = 'ml_input/'
    if not os.path.exists(directory):
        os.mkdir(directory)

    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] == -1:
            obj = species.reac_obj[index]
            name = species.reac_name[index]
            if not os.path.exists(directory + name):
                os.mkdir(directory + name)
            # make the reactant file
            s = ['{}'.format(species.natom)]
            s.append(' '.join(species.atom))
            s.append('\n')
            for bi in species.bond:
                s.append(' '.join([str(bij) for bij in bi]))
            s.append('\n')
            for i in range(species.natom):
                row = []
                for j in range(species.natom):
                    d = np.linalg.norm(species.geom[i] - species.geom[j])
                    row.append('{:.2f}'.format(d))
                s.append(' '.join(row))
            s.append('\n')
            with open(directory + name + '/reactant.txt', 'w') as f:
                f.write('\n'.join(s))
            # make the product file
            s = ['{}'.format(species.natom)]
            s.append(' '.join(species.atom))
            s.append('\n')
            for bi in obj.product_bonds:
                s.append(' '.join([str(bij) for bij in bi]))
            s.append('\n')
            with open(directory + name + '/product.txt', 'w') as f:
                f.write('\n'.join(s))
            # write the ts key distances
            s = []
            for i in range(species.natom):
                row = []
                for j in range(species.natom):
                    d = 0.0
                    if species.bond[i][j] != obj.product_bonds[i][j]:
                        d = np.linalg.norm(obj.ts.geom[i] - obj.ts.geom[j])
                    row.append('{:.2f}'.format(d))
                s.append(' '.join(row))
            s.append('\n')
            with open(directory + name + '/ts.txt', 'w') as f:
                f.write('\n'.join(s))


def create_summary_file(species, qc, par):
    """
    Create a summary file listing for each reaction
    1. whether its search was successful
    2. the barrier height
    3. the reaction name
    4. the product identifiers
    """
    # list of strings which will be put together for the output
    s = []
    # add the license message to the file
    s.append(license_message.message)
    # list of the products
    products = []
    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] == -1:
            ts = species.reac_obj[index].ts
            if species.reac_type[index] == 'R_Addition_MultipleBond' and not par['high_level']:
                mp2_energy = qc.get_qc_energy(str(species.chemid) + '_well_mp2')[1]
                mp2_zpe = qc.get_qc_zpe(str(species.chemid) + '_well_mp2')[1]
                energy = (ts.energy + ts.zpe - mp2_energy - mp2_zpe) * constants.AUtoKCAL
            else:
                energy = (ts.energy + ts.zpe - species.energy - species.zpe) * constants.AUtoKCAL
            prod_name = ''
            name = []
            for prod in species.reac_obj[index].products:
                name.append(str(prod.chemid))
            prod_name = ' '.join(sorted(name))
            products.append(prod_name)
            s.append('SUCCESS\t{energy:.2f}\t{name}\t{prod}'.format(energy=energy,
                                                                    name=species.reac_name[index],
                                                                    prod=prod_name))
        else:
            s.append('FAILED\t\t{name}'.format(name=species.reac_name[index]))

    # make a string out of all the lines
    s = '\n'.join(s)
    # write the string to a file
    fname = 'summary_{chemid}.out'.format(chemid=species.chemid)
    with open(fname, 'w') as f:
        f.write(s)


def createPESViewerInput(species, qc, par):
    """
    Write an input file for the PESViewer code
    """
    # make the directory for the well and bimolecular xyz files
    dir_xyz = 'xyz/'
    if not os.path.exists(dir_xyz):
        os.mkdir(dir_xyz)

    # list of the lines for the pesviewer input file
    wells = []
    # list of the names of the wells
    well_names = [str(species.chemid)]
    # make an xyz file for the initial well
    make_xyz(species.atom, species.geom, str(species.chemid), dir_xyz)
    # add the initial well to the wells list
    # use this well as point zero for the energy
    wells.append('{} 0.0'.format(species.chemid))
    well_energy = species.energy + species.zpe

    # iterate the reactions and search for single products
    # i.e. other wells on the pes
    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] == -1:
            if len(species.reac_obj[index].prod_opt) == 1:
                st_pt = species.reac_obj[index].prod_opt[0].species
                name = str(st_pt.chemid)
                if name not in well_names:
                    make_xyz(species.atom, st_pt.geom, str(st_pt.chemid), dir_xyz)
                    energy = (st_pt.energy + st_pt.zpe - well_energy) * constants.AUtoKCAL
                    wells.append('{name} {energy:.2f}'.format(name=st_pt.chemid, energy=energy))
                    well_names.append(name)

    # list of the lines for the pesviewer input file
    bimolecs = []
    # list of the names of the bimolecular products
    bimolec_names = []
    # add the bimolecular products from the regular reactions
    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] == -1:
            if len(species.reac_obj[index].prod_opt) > 1:
                energy = 0. - well_energy
                names = []
                for prod_opt in species.reac_obj[index].prod_opt:
                    st_pt = prod_opt.species
                    energy += st_pt.energy + st_pt.zpe
                    names.append(str(st_pt.chemid))
                name = '_'.join(sorted(names))

                for i, prod_opt in enumerate(species.reac_obj[index].prod_opt):
                    st_pt = prod_opt.species
                    with open("pesviewer_data.txt", 'a') as pesdata:
                        pesdata.write("Species: {}\n\tEnergy: {}\n\tZPE: {}\n".format(st_pt.chemid, st_pt.energy, st_pt.zpe))
                    # make twice the same file but with a different name
                    # TODO: is there no better way?
                    # this is for the pes viewer
                    make_xyz(st_pt.atom, st_pt.geom, name + str(i + 1), dir_xyz)
                    # this is for the rmg postprocessing
                    make_xyz(st_pt.atom, st_pt.geom, str(st_pt.chemid), dir_xyz)
                energy = energy * constants.AUtoKCAL
                if name not in bimolec_names:
                    bimolecs.append('{name} {energy:.2f}'.format(name=name, energy=energy))
                    bimolec_names.append(name)

    # list of the lines of the ts's
    tss = []
    # dict keeping track of the ts's
    # key: ts name
    # value: [energy,prod_names]
    ts_list = {}
    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] == -1:
            ts = species.reac_obj[index].ts
            if species.reac_type[index] == 'R_Addition_MultipleBond' and not par['high_level']:
                we_energy = qc.get_qc_energy(str(species.chemid) + '_well_mp2')[1]
                we_zpe = qc.get_qc_zpe(str(species.chemid) + '_well_mp2')[1]
                energy = (ts.energy + ts.zpe - we_energy - we_zpe) * constants.AUtoKCAL
            else:
                energy = (ts.energy + ts.zpe - well_energy) * constants.AUtoKCAL
            name = []
            for st_pt in species.reac_obj[index].products:
                name.append(str(st_pt.chemid))
            prod_name = '_'.join(sorted(name))
            add = 1
            for t in ts_list:
                if ts_list[t][1] == prod_name and np.abs(ts_list[t][0] - energy) < 1.0:
                    add = 0
            if add:
                ts_list[species.reac_name[index]] = [energy, prod_name]
                tss.append('{ts} {energy:.2f} {react} {prod}'.format(ts=species.reac_name[index],
                                                                     energy=energy,
                                                                     react=species.chemid,
                                                                     prod=prod_name))

    # make strings from the different lists
    wells = '\n'.join(wells)
    bimolecs = '\n'.join(bimolecs)
    tss = '\n'.join(tss)
    barrierless = ''

    # write everything to a file
    fname = 'pesviewer.inp'
    template_file_path = f'{kb_path}/tpl/{fname}.tpl'
    with open(template_file_path) as template_file:
        template = template_file.read()
    template = template.format(id=species.chemid, wells=wells, bimolecs=bimolecs, ts=tss, barrierless=barrierless)
    with open(fname, 'w') as f:
        f.write(template)


def make_xyz(atoms, geom, name, directory):
    s = []
    s.append('%i\n' % len(geom))
    for index in range(len(geom)):
        s.append('%s %.6f %.6f %.6f' % (atoms[index], geom[index][0], geom[index][1], geom[index][2]))
    with open(directory + '/' + name + '.xyz', 'w') as f:
        f.write('\n'.join(s))
