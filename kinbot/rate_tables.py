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
import os,sys

class Rates:
    """
    Rates belonging to a paricular reactant/product pair
    """
    def __init__(self,reactant,product):
        self.reactant = reactant
        self.product = product
        self.temperatures = [] #list of temperatures
        self.pressures = [] #list of pressures
        self.rates = [] # 2D list, first corresponding to temp, than to p
