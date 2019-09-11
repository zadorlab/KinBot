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


def copy_from_database_folder(dir_name, frag, qc):

    file_list = os.listdir(dir_name)
    for file in file_list:
        if '.com' in file or '.log' in file:
            copyfile(dir_name + '/' + file, os.getcwd() + file)

    hir_file_list = os.listdir(dir_name + '/hir/')
    for file in hir_file_list:
        if '.com' in file or '.log' in file:
            copyfile(dir_name + '/hir/' + file, os.getcwd() + '/hir/' + file)

    conf_file_list = os.listdir(dir_name + '/conf/')
    for file in conf_file_list:
        if '.com' in file or '.log' in file:
            copyfile(dir_name + '/conf/' + file, os.getcwd() + '/conf/' + file)

    # read the database and populate the current database
    data = connect(dir_name + '/{}.db'.format(frag.chemid)).select()
    for row in data:
        ase_atoms = row.toatoms()
        ase_name = row.name
        ase_data = row.data
        qc.db.write(ase_atoms, name=ase_name, data=ase_data)


def copy_to_database_folder(dir_name, frag, qc):
    file_list = os.listdir(os.getcwd())
    for file in file_list:
        if '{}_well'.format(frag.chemid) in file and '.log' in file:
            copyfile(os.getcwd() + '/' + file, dir_name + '/' + file)
        if '{}_well'.format(frag.chemid) in file and '.com' in file:
            copyfile(os.getcwd() + '/' + file, dir_name + '/' + file)
    
    hir_file_list = os.listdir(os.getcwd() + '/hir/')
    if not os.path.exists(dir_name + '/hir/'):
        os.makedirs(dir_name + '/hir/')
    for file in hir_file_list:
        if '{}_hir'.format(frag.chemid) in file and '.log' in file:
            copyfile(os.getcwd() + '/hir/' + file, dir_name + '/hir/' + file)
        if '{}_hir'.format(frag.chemid) in file and '.com' in file:
            copyfile(os.getcwd() + '/hir/' + file, dir_name + '/hir/' + file)

    conf_file_list = os.listdir(os.getcwd() + '/conf/')
    if not os.path.exists(dir_name + '/conf/'):
        os.makedirs(dir_name + '/conf/')
    for file in conf_file_list:
        if '{}'.format(frag.chemid) in file and '.log' in file:
            copyfile(os.getcwd() + '/conf/' + file, dir_name + '/conf/' + file)
        if '{}'.format(frag.chemid) in file and '.com' in file:
            copyfile(os.getcwd() + '/conf/' + file, dir_name + '/conf/' + file)
    # read the database and populate the current database
    ase_db = connect(dir_name + '/{}.db'.format(frag.chemid))
    for row in qc.db.select():
        if '{}'.format(frag.chemid) in row.name:
            ase_atoms = row.toatoms()
            ase_name = row.name
            ase_data = row.data
            ase_db.write(ase_atoms, name=ase_name, data=ase_data)
