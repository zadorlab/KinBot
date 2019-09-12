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
from shutil import copyfile

from ase.db import connect

def copy_from_database_folder(well0_chemid, chemid, qc):
    # wait for this well to finish
    wait = 1
    # directory of the pes run
    dir = os.path.dirname(os.getcwd()) 
    # dir for this well
    dir_name = '{}/{}_db/'.format(dir, chemid)
    # check if the directory is there
    if os.path.exists(dir_name):
        # check for the running tag, and if it contains the chemid
        # of this well, continue from here
        temp_chemid = -1
        try:
            with open(dir_name + 'running') as f:
                temp_chemid = int(f.read().split()[0])
        except IOError:
            logging.error('Could not find running tag from well db directory.')
        if temp_chemid == well0_chemid:
            wait = 0
        elif os.path.exists(dir_name + 'done'):
            file_list = os.listdir(dir_name)
            for file in file_list:
                if '.com' in file or '.log' in file or '.fchk' in file:
                    copyfile(dir_name + file, os.getcwd() + '/' + file)
            try:
				hir_file_list = os.listdir(dir_name + 'hir/')
				for file in hir_file_list:
					if '.com' in file or '.log' in file:
						copyfile(dir_name + 'hir/' + file, os.getcwd() + '/hir/' + file)
            except IOError:
                logging.warning("hir dir/file not found for " + dir_name)
 
            try:
				conf_file_list = os.listdir(dir_name + 'conf/')
				for file in conf_file_list:
					if '.com' in file or '.log' in file:
						copyfile(dir_name + 'conf/' + file, os.getcwd() + '/conf/' + file)
            except IOError:
                logging.warning("conf_file dir/file not found for " + dir_name)
 
            # read the database and populate the current database
            data = connect(dir_name + '{}.db'.format(chemid))
            for row in data.select():
                ase_atoms = row.toatoms()
                ase_name = row.name
                ase_data = row.data
                qc.db.write(ase_atoms, name=ase_name, data=ase_data)
            wait = 0
    else:
        # directory is not yet made, make it now
        os.makedirs(dir_name)
        # make the running tag
        with open(dir_name + 'running', 'w') as f:
            f.write('{}'.format(well0_chemid))
        wait = 0
    
    return wait

def copy_to_database_folder(well0_chemid, chemid, qc):
    # directory of the pes run
    dir = os.path.dirname(os.getcwd()) 
    # dir for this well
    dir_name = '{}/{}_db/'.format(dir, chemid)
    # check if the directory is there
    if os.path.exists(dir_name):
        # check for the running tag, and if it contains the chemid
        # of this well, continue from here
        if not os.path.exists(dir_name + 'done'):
            if os.path.exists(dir_name + 'running'):
                temp_chemid = -1
                try:
                    with open(dir_name + 'running') as f:
                        temp_chemid = int(f.read().split()[0])
                except IOError:
                    logging.error('Could not find running tag from well db directory.')
                if temp_chemid == well0_chemid:
                    # copy the files
                    file_list = os.listdir(os.getcwd())
                    for file in file_list:
                        if '{}_well'.format(chemid) in file and '.log' in file:
                            copyfile(os.getcwd() + '/' + file, dir_name + file)
                        if '{}_well'.format(chemid) in file and '.com' in file:
                            copyfile(os.getcwd() + '/' + file, dir_name + file)
                        if '{}_well'.format(chemid) in file and '.fchk' in file:
                            copyfile(os.getcwd() + '/' + file, dir_name + file)

                    hir_file_list = os.listdir(os.getcwd() + '/hir/')
                    if not os.path.exists(dir_name + 'hir/'):
                        os.makedirs(dir_name + 'hir/')
                    for file in hir_file_list:
                        if '{}_hir'.format(chemid) in file and '.log' in file:
                            copyfile(os.getcwd() + '/hir/' + file, dir_name + 'hir/' + file)
                        if '{}_hir'.format(chemid) in file and '.com' in file:
                            copyfile(os.getcwd() + '/hir/' + file, dir_name + 'hir/' + file)

                    conf_file_list = os.listdir(os.getcwd() + '/conf/')
                    if not os.path.exists(dir_name + 'conf/'):
                        os.makedirs(dir_name + 'conf/')
                    for file in conf_file_list:
                        if '{}'.format(chemid) in file and '.log' in file:
                            copyfile(os.getcwd() + '/conf/' + file, dir_name + 'conf/' + file)
                        if '{}'.format(chemid) in file and '.com' in file:
                            copyfile(os.getcwd() + '/conf/' + file, dir_name + 'conf/' + file)
                    # read the database and populate the current database
                    ase_db = connect(dir_name + '{}.db'.format(chemid))
                    for row in qc.db.select():
                        if '{}'.format(chemid) in row.name:
                            ase_atoms = row.toatoms()
                            ase_name = row.name
                            ase_data = row.data
                            ase_db.write(ase_atoms, name=ase_name, data=ase_data)
                # make a done tag
                with open(dir_name + 'done', 'w') as f:
                    f.write('')
    
    
