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
        if file != 'done' and file != 'running':
            if '_hir_' in file:
                copyfile(dir_name + '/' + file, os.getcwd() + '/hir/' + file)
            elif '_conf_' in file:
                copyfile(dir_name + '/' + file, os.getcwd() + '/conf/' + file)
            else:
                copyfile(dir_name + '/' + file, os.getcwd() + '/' + file)
        # read the database and populate the current database
    data = connect(dir_name + '/{}.db'.format(frag.chemid)).select()
    for row in data:
        qc.db.write(row)

def copy_to_database_folder(dir_name, frag, qc):
    file_list = os.listdir(os.getcwd())
    for file in file_list:
        if '{}'.format(frag.chemid) in file and '_well' in file and '.log' in file:
            copyfile(os.getcwd() + '/' + file, dir_name + '/' + file)
        if '{}'.format(frag.chemid) in file and '_well' in file and '.com' in file:
            copyfile(os.getcwd() + '/' + file, dir_name + '/' + file)
        if '{}'.format(frag.chemid) in file and '_hir_' in file and '.log' in file:
            copyfile(os.getcwd() + '/' + file, dir_name + '/hir/' + file)
        if '{}'.format(frag.chemid) in file and '_hir_' in file and '.com' in file:
            copyfile(os.getcwd() + '/' + file, dir_name + '/hir/' + file)
        if '{}'.format(frag.chemid) in file and '_conf_' in file and '.log' in file:
            copyfile(os.getcwd() + '/' + file, dir_name + '/conf/' + file)
        if '{}'.format(frag.chemid) in file and '_conf_' in file and '.com' in file:
            copyfile(os.getcwd() + '/' + file, dir_name + '/conf/' + file)
        # read the database and populate the current database
    data = connect(dir_name + '/{}.db'.format(frag.chemid)).select()
    for row in self.qc.db.select():
        if '{}'.format(frag.chemid) in row.name:
            data.write(row)