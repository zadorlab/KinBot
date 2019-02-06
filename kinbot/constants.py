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

AUtoKCAL = 627.5091809
AUtoCM = 219474.63068
CMtoKCAL = 0.0028591
GHZtoCM = 0.0333564
EVtoHARTREE = 0.03674932247495664
MEtoAMU = 1822.8885  # electron mass to atomic mass unit
BOHRtoCM = 5.2917720859E-9  # value taken from the Gaussian website
SPEEDofLIGHTcms = 2.99792458E10  # in cm per s
SPEEDofLIGHT = 137.0359996  # in atomic units
AUtoS = 2.418884326505E-17

# standard bond lengths, cutoffs, and oxidation numbers
st_bond = {'CC': 1.5*1.2}
# st_bond['CO'] = 1.4*1.2 # 1.4*1.54
st_bond['CO'] = 1.4*1.4
st_bond['CH'] = 1.09*1.2
st_bond['OO'] = 1.48*1.2
st_bond['HO'] = 0.95*1.2
st_bond['HH'] = 1.1  # 0.74*1.2
st_bond['OS'] = 1.82*1.2
st_bond['HS'] = 1.34*1.2
st_bond['SS'] = 2.07*1.2
st_bond['HN'] = 0.99*1.2
st_bond['CN'] = 1.47*1.2
st_bond['NO'] = 1.49*1.2
st_bond['NS'] = 1.6*1.2
st_bond['NN'] = 1.35*1.2

st_bond['ts'] = 1.15
st_bond['S'] = 6  # maximum valence, can also be 2 or 4
st_bond['C'] = 4
st_bond['O'] = 2
st_bond['H'] = 1
st_bond['N'] = 5  # maximum valence, can also be 3 

mass = {'H': 1}
mass['C'] = 12
mass['N'] = 14
mass['O'] = 16
mass['S'] = 32

exact_mass = {'H': 1.00782503207}
exact_mass['C'] = 12.00000000000
exact_mass['N'] = 14.003074
exact_mass['O'] = 15.99491461956
exact_mass['S'] = 31.97207100

znumber = {'H': 1}
znumber['C'] = 6
znumber['N'] = 7
znumber['O'] = 8
znumber['S'] = 16


def main():
    """
    This module contains constants and run-specific settings.
    """


if __name__ == "__main__":
    main()
