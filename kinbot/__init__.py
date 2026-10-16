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
##   Judit Zádor                                 ##
##   Ruben Van de Vijver                         ##
##   Austin Dewyer                               ##
##   Carles Martí                                ##
##   Clément Soulié                              ##
##   Daniel T. Kendall                           ##
##   Luka Dockx                                  ##
##   Hrishikesh Ram                              ##
###################################################

kb_path = __path__[0]

from kinbot.rdkit_config import configure_rdkit

configure_rdkit()
