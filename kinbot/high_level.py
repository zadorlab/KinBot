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
import os, sys

from qc import *



def high_level(sp, ts = 0):
    """
    Reoptimize the stationary point and recalculate all the frequencies at a high level of theory
    
    sp: stationary point
    
    ts: boolean that tells if stationary point is a transition state
    """
    
    #reoptimize the stationary point
    qc_opt(well0, well0.geom, 0, par.natom, par.atom, par.mult, par.charge, high_level = 1)

    