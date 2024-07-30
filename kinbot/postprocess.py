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
            for bi in obj.irc_prod.bonds[0]:
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
                    if species.bond[i][j] != obj.irc_prod.bonds[0][i][j]:
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
    optionally:
    5. the depth of vdW well
    6. the IRC direction to go to the vdW well
    """
    # list of strings which will be put together for the output
    s = []
    # add the license message to the file
    s.append(license_message.message)

    l3status = hasattr(species, 'l3energy')
    max_len = 0
    for index in range(len(species.reac_inst)):
        if len(species.reac_name[index]) > max_len:
            max_len = len(species.reac_name[index])
        if not l3status or species.reac_ts_done[index] != -1:
            continue
        l3status *= all([hasattr(p, 'l3energy')
                         for p in species.reac_obj[index].products])
        l3status *= hasattr(species.reac_obj[index].ts, 'l3energy')

    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] == -1:
            ts = species.reac_obj[index].ts
            if l3status and species.reac_type[index] != 'hom_sci':
                energy = (ts.l3energy + ts.zpe
                          - species.l3energy - species.zpe)
                energy *= constants.AUtoKCAL
            elif species.reac_type[index] == 'R_Addition_MultipleBond' \
                    and not par['high_level'] \
                    and qc.qc != 'nn_pes':
                mp2_energy = qc.get_qc_energy(str(species.chemid)
                                              + '_well_mp2')[1]
                mp2_zpe = qc.get_qc_zpe(str(species.chemid) + '_well_mp2')[1]
                energy = (ts.energy + ts.zpe
                          - mp2_energy - mp2_zpe) * constants.AUtoKCAL
            elif species.reac_type[index] == 'hom_sci':
                if l3status:
                    energy = (sum([pr.l3energy + pr.zpe
                                   for pr in species.reac_obj[index].products])
                              - species.energy - species.zpe)
                    energy *= constants.AUtoKCAL
                else:
                    energy = (sum([pr.energy + pr.zpe
                                   for pr in species.reac_obj[index].products])
                              - species.energy - species.zpe)
                    energy *= constants.AUtoKCAL
            else:
                energy = (ts.energy + ts.zpe - species.energy - species.zpe)
                energy *= constants.AUtoKCAL

            name = []
            for prod in species.reac_obj[index].products:
                name.append(str(prod.chemid))
            prod_name = ' '.join(sorted(name))
            status = "SUCCESS"
            if species.reac_obj[index].do_vdW:   
                vdW_energy = (species.reac_obj[index].irc_prod.energy +\
                              species.reac_obj[index].irc_prod.zpe -\
                              (species.energy + species.zpe))*constants.AUtoKCAL
                direction ="vdW{}".format(species.reac_obj[index]\
                                          .irc_prod.name.split(species.reac_obj[index]\
                                                               .instance_name)[1])        
                s.append('{status:7s}{energy:> 9.2f}  {name:{max_len}s} {prod} {vdW_energy:> 7.2f}  {direction}'.format(status=status,
                                                                    energy=energy,
                                                                    max_len=max_len+1,
                                                                    name=species.reac_name[index],
                                                                    prod=prod_name,
                                                                    vdW_energy=vdW_energy,
                                                                    direction=direction))
            else:
                s.append('{status:7s}{energy:> 9.2f}  {name:{max_len}s} {prod}'.format(status=status,
                                                                    energy=energy,
                                                                    name=species.reac_name[index],
                                                                    prod=prod_name,
                                                                    max_len=max_len+1))
        else:
            status = "FAILED"
            s.append('{status:16s}  {name:{max_len}s}'.format(status=status,
                                                 name=species.reac_name[index],
                                                 max_len=max_len+1))

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
        if species.reac_ts_done[index] != -1 \
                or len(species.reac_obj[index].prod_opt) != 1:
            continue
        st_pt = species.reac_obj[index].prod_opt[0].species
        name = str(st_pt.chemid)
        if name in well_names:
            continue
        make_xyz(species.atom, st_pt.geom, str(st_pt.chemid), dir_xyz)
        energy = (st_pt.energy + st_pt.zpe - well_energy) * constants.AUtoKCAL
        wells.append(f'{st_pt.chemid} {energy:.2f}')
        well_names.append(name)

    # list of the lines for the pesviewer input file
    bimolecs = []
    # list of the names of the bimolecular products
    bimolec_names = []
    # add the bimolecular products from the regular reactions
    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] != -1 \
                or len(species.reac_obj[index].prod_opt) <= 1:
            continue
        energy = 0. - well_energy
        names = []
        for prod_opt in species.reac_obj[index].prod_opt:
            st_pt = prod_opt.species
            energy += st_pt.energy + st_pt.zpe
            names.append(str(st_pt.chemid))
        name = '_'.join(sorted(names))

        for i, prod_opt in enumerate(species.reac_obj[index].prod_opt):
            st_pt = prod_opt.species
            with open('pesviewer_data.txt', 'a') as pesdata:
                pesdata.write(f'Species: {st_pt.chemid}\n'
                              f'\tEnergy: {st_pt.energy}\n'
                              f'\tZPE: {st_pt.zpe}\n')
            # make twice the same file but with a different name
            # TODO: is there no better way?
            # this is for the pes viewer
            make_xyz(st_pt.atom, st_pt.geom, name + str(i + 1), dir_xyz)
            # this is for the rmg postprocessing
            make_xyz(st_pt.atom, st_pt.geom, str(st_pt.chemid), dir_xyz)
        energy = energy * constants.AUtoKCAL
        if name not in bimolec_names:
            bimolecs.append(f'{name} {energy:.2f}')
            bimolec_names.append(name)
        if species.reac_obj[index].do_vdW:
            irc_prod = species.reac_obj[index].irc_prod_opt.species
            name = str(irc_prod.name)
            if name in well_names:
                continue
            make_xyz(species.atom, irc_prod.geom, str(irc_prod.name), dir_xyz)
            energy = (irc_prod.energy + irc_prod.zpe - well_energy) * constants.AUtoKCAL
            wells.append(f'{irc_prod.name} {energy:.2f}')
            well_names.append(name)
        

    # list of the lines of the ts's
    tss = []
    # list of the lines of the barrierless
    bless = []
    # dict keeping track of the ts's
    # key: ts name
    # value: [energy,prod_names]
    ts_list = {}
    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] != -1:
            continue
        ts = species.reac_obj[index].ts
        if species.reac_type[index] == 'R_Addition_MultipleBond' \
                and not par['high_level'] \
                and qc.qc != 'nn_pes':
            we_energy = qc.get_qc_energy(str(species.chemid) + '_well_mp2')[1]
            we_zpe = qc.get_qc_zpe(str(species.chemid) + '_well_mp2')[1]
            energy = (ts.energy + ts.zpe - we_energy - we_zpe) * constants.AUtoKCAL
        elif species.reac_type[index] == 'hom_sci':
            continue
        else:
            energy = (ts.energy + ts.zpe - well_energy) * constants.AUtoKCAL
        name = []
        if species.reac_obj[index].do_vdW:
            irc_prod = species.reac_obj[index].irc_prod_opt.species
            name.append(irc_prod.name)
            bimol_names = []
            for st_pt in species.reac_obj[index].products:
                bimol_names.append(str(st_pt.chemid))
            bimol_prod_name = '_'.join(sorted(bimol_names))
            bless.append(f'vdW_{index} {irc_prod.name} {bimol_prod_name}')
        else:
            for st_pt in species.reac_obj[index].products:
                name.append(str(st_pt.chemid))
        prod_name = '_'.join(sorted(name))
        add = 1
        for t in ts_list:
            if ts_list[t][1] == prod_name and np.abs(ts_list[t][0] - energy) < 1.0:
                add = 0
        if add:
            ts_list[species.reac_name[index]] = [energy, prod_name]
            tss.append(f'{species.reac_name[index]} {energy:.2f} '
                       f'{species.chemid} {prod_name}')
    
    # Barrierless reactions
    for index in range(len(species.reac_inst)):
        if species.reac_ts_done[index] != -1 or species.reac_type[index] != 'hom_sci':
            continue
        name = []
        for st_pt in species.reac_obj[index].products:
            name.append(str(st_pt.chemid))
        prod_name = '_'.join(sorted(name))
        if prod_name not in [v[1] for v in ts_list.values()]:
            bless.append(f'{species.reac_name[index]} {species.chemid} {prod_name}')

    # make strings from the different lists
    wells = '\n'.join(wells)
    bimolecs = '\n'.join(bimolecs)
    tss = '\n'.join(tss)
    barrierless = '\n'.join(bless)

    # write everything to a file
    fname = 'pesviewer.inp'
    template_file_path = f'{kb_path}/tpl/{fname}.tpl'
    with open(template_file_path) as template_file:
        template = template_file.read()
    template = template.format(id=species.chemid, wells=wells, ts=tss, 
                               bimolecs=bimolecs, barrierless=barrierless)
    with open(fname, 'w') as f:
        f.write(template)


def make_xyz(atoms, geom, name, directory):
    s = []
    s.append('%i\n' % len(geom))
    for index in range(len(geom)):
        s.append('%s %.6f %.6f %.6f' % (atoms[index], geom[index][0], 
                                        geom[index][1], geom[index][2]))
    with open(directory + '/' + name + '.xyz', 'w') as f:
        f.write('\n'.join(s))
