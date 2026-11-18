"""
File that contains exceptions for molecules that need
to have different multiplicity/spin or symmetry configurations
than the defaults set by Gaussian or other Quantum Chemistry programs
    Ex: CH2 radical, O2, OH radical
"""

# species
# 320320000000000000001 = O2
# 140260020000000000001 = CH2
# 170170000000000000002 = OH
# 320000000000000000001 = S


def is_proton(atoms, charge, multiplicity):
    """An isolated H nucleus has no electrons or vibrational modes."""
    return len(atoms) == 1 and atoms[0] == 'H' and charge == 1 and multiplicity == 1


def get_spin(chemid, spin):  # spin as defined in molpro = 2S
    if chemid == '320320000000000000001':
        spin = 2
    elif chemid == '140260020000000000001':
        spin = 2
    elif chemid == '170170000000000000002':
        spin == 2
    else:
        spin = spin
    return spin


def get_multiplicity(chemid, mult):  # 2S+1
    if chemid == '320320000000000000001':
        mult = 3
    elif chemid == '140260020000000000001':
        mult = 3
    elif chemid == '170170000000000000002':
        mult = 3
    else:
        mult = mult
    return mult


def get_symm(chemid, symm):  # molpro symm defined by point group
    if chemid == '320320000000000000001':
        symm = 1  # used Ag state, label = 1 for molpro
    elif chemid == '140260020000000000001':
        symm = 2  # B1 state, C2V symmetry, label = 2 for molpro
    elif chemid == '170170000000000000002':
        symm = 2  # 2Pi state ~ using B state = 2 for molpro
    else:
        symm = symm
    return symm
