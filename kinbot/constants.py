from ase import units
from ase import data

AUtoKCAL = 627.5091809
AUtoCM = 219474.63068
CMtoKCAL = 0.0028591
GHZtoCM = 0.0333564
EVtoHARTREE = 0.03674932247495664
MEtoAMU = 1822.8885  # electron mass to atomic mass unit
BOHRtoCM = 5.2917720859E-9  # value taken from the Gaussian website
BOHRtoANGSTROM = 0.529177
SPEEDofLIGHTcms = 2.99792458E10  # in cm per s
SPEEDofLIGHT = 137.0359996  # in atomic units
AUtoS = 2.418884326505E-17
AUtoMDYNE = 8.2387234983  # au of force to mDyne = 1e-8 N
KCALtoHARTREE = 0.0015936010974213599
R = 8.31446261815324  # J K−1 mol−1
CALtoJ = 4.184

# elements currently in KinBot
elements = ['C', 'H', 'O', 'N', 'S', 'F', 'Cl', 'Br', 'I']

znumber = {'H': 1}
znumber['C'] = 6
znumber['N'] = 7
znumber['O'] = 8
znumber['S'] = 16
znumber['F'] = 9
znumber['Cl'] = 17
znumber['Br'] = 35
znumber['I'] = 53

# standard bond lengths, cutoffs, and oxidation numbers

st_bond = {'CC': 1.5*1.2}
# st_bond['CO'] = 1.4*1.2 # 1.4*1.54
st_bond['CO'] = 1.4*1.4
st_bond['CH'] = 1.09*1.2
st_bond['CF'] = 1.37*1.2  # typical bond length = 1.37
st_bond['CCl'] = 1.76*1.2  # typical bond length = 1.76
st_bond['CBr'] = 1.94*1.2  # typical bond length = 1.94
st_bond['CI'] = 2.14*1.2  # typical bond length = 2.14
st_bond['OO'] = 1.48*1.2
st_bond['HO'] = 0.95*1.2
st_bond['HH'] = 1.1  # 0.74*1.2
st_bond['OS'] = 1.82*1.2
st_bond['HS'] = 1.34*1.2
st_bond['CS'] = 1.82*1.2
st_bond['SS'] = 2.07*1.2
st_bond['HN'] = 0.99*1.2
st_bond['CN'] = 1.47*1.2
st_bond['NO'] = 1.49*1.2
st_bond['NS'] = 1.6*1.2
st_bond['NN'] = 1.35*1.2
st_bond['FF'] = 1.55*1.2
st_bond['FO'] = 1.41*1.2
st_bond['FH'] = 0.91*1.2
st_bond['ClCl'] = 1.99*1.2
st_bond['ClO'] = 1.70*1.2
st_bond['ClH'] = 1.275*1.2
st_bond['BrBr'] = 2.29*1.2
st_bond['II'] = 1.48*1.2

st_bond['ts'] = 1.15
st_bond['S'] = 6  # maximum valence, can also be 2 or 4
st_bond['C'] = 4
st_bond['O'] = 2
st_bond['H'] = 1
st_bond['N'] = 5  # maximum valence, can also be 3
st_bond['F'] = 1
st_bond['Cl'] = 1
st_bond['Br'] = 1
st_bond['I'] = 1

for el1 in elements:
    for el2 in elements:
        try:
            st_bond[f'{el1}{el2}']
        except KeyError:
            try:
                st_bond[f'{el2}{el1}']
            except KeyError:
                z1 = znumber[el1]
                r1 = data.covalent_radii[z1]
                z2 = znumber[el2]
                r2 = data.covalent_radii[z2]
                st_bond[f'{el1}{el2}'] = (r1 + r2) * 1.2 
#print(st_bond)

mass = {'H': 1}
mass['C'] = 12
mass['N'] = 14
mass['O'] = 16
mass['S'] = 32
mass['F'] = 19
mass['Cl'] = 35  # isotopes = Cl-35, Cl-37, Cl-36 (3 most stable/common)
mass['Br'] = 79  # isotopes = Br-79, Br-81, Br-77 (3 most stable/common)
mass['I'] = 127  

exact_mass = {'H': 1.00782503207}
exact_mass['C'] = 12.00000000000
exact_mass['N'] = 14.003074
exact_mass['O'] = 15.99491461956
exact_mass['S'] = 31.97207100
exact_mass['F'] = 18.9984 
exact_mass['Cl'] = 34.9689
exact_mass['Br'] = 78.9183
exact_mass['I'] = 126.904477

# collision parameters
# Jasper & Miller, C&F 161, 101-110 (2014)
# data is expected in CHEMKIN format, where
# sigma is in Angstrom (no need to convert for MESS)
# epsilon is in epsilon/kB (K), need to be converted to cm-1
mass['He'] = 4.0
sigma = {'He': 2.715}  # Angstrom
epsilon = {'He': 11.442}  # e/kB

mass['N2'] = mass['N'] * 2.
sigma['N2'] = 3.610
epsilon['N2'] = 97.839

mass['Ar'] = 40.
sigma['Ar'] = 3.462
epsilon['Ar'] = 127.697

for e in epsilon:
    epsilon[e] = epsilon[e] * units.kB / units.invcm

# submission keywords
qsubmit = {'pbs': 'qsub'}
qsubmit['slurm'] = 'sbatch'
# extensions
qext = {'pbs': '.pbs',
        'slurm': '.sbatch',
        }


def main():
    """
    This module contains constants and run-specific settings.
    """


if __name__ == "__main__":
    main()
