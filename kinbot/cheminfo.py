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
import sys
import os
import numpy as np
import pkg_resources
from PIL import Image

# try to import pybel
try:
    import pybel
except ImportError:
    pass

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    pass

num_to_syms = {1: 'H', 6: 'C', 8: 'O', 16: 'S'}
syms_to_num = {'H': 1, 'C': 6, 'O': 8, 'S': 16}


def get_molecular_formula(smi):
    """
    Return the molecular formula of the molecule corresponding to the smiles
    """
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    return rdMolDescriptors.CalcMolFormula(mol)


def create_rxn_depiction(react_smiles, prod_smiles, dir, name):
    """
    Create a 2D depiction of a chemical reaction,
    react smiles: smiles of the reactants
    prod_smiles: smiles of the products
    dir: directory where to save the reaction depiction
    name: name of the depiction files
    """
    react_png = '{dir}/react.png'.format(dir=dir)
    prod_png = '{dir}/prod.png'.format(dir=dir)

    obmol = pybel.readstring("smi", react_smiles)
    obmol.draw(show=False, filename=react_png)

    obmol = pybel.readstring("smi", prod_smiles)
    obmol.draw(show=False, filename=prod_png)

    arrow = pkg_resources.resource_filename('tpl', 'arrow.png')
    images = map(Image.open, [react_png, arrow, prod_png])
    widths, heights = zip(*(i.size for i in images))

    total_width = sum(widths)
    total_height = max(heights)

    new_im = Image.new('RGB', (total_width, total_height), (255, 255, 255))

    x = 0
    for im in images:
        y = total_height / 2 - im.size[1] / 2
        new_im.paste(im, (x, y))
        x += im.size[0]

    new_im.save('{dir}/{name}.png'.format(dir=dir, name=name))


def generate_3d_structure(smi, obabel=1):
    """
    Method to generate the 3D coordinates of a molecule from its smiles
    The default code is OpenBabel, RDKit can also be used.
    """
    structure = []
    if obabel:  # use OpenBabel
        try:
            obmol = pybel.readstring('smi', smi)
        except NameError:
            message = '\nPybel is required to use the smiles input format.\n'
            message += 'If pybel is unavailable, use the geometry as input.\n'
            message += 'Else install OpenBabel with python bindings.\nExiting...\n'
            sys.exit(message)
        obmol.OBMol.AddHydrogens()
        obmol.make3D()
        bond = np.zeros((len(obmol.atoms), len(obmol.atoms)), dtype=int)
        for i in range(len(obmol.atoms)):
            for j in range(len(obmol.atoms)):
                if not obmol.OBMol.GetBond(i+1, j+1) is None:
                    order = obmol.OBMol.GetBond(i+1, j+1).GetBO()
                    bond[i][j] = order
        for at in obmol.atoms:
            pos = at.coords
            sym = num_to_syms[at.atomicnum]
            structure += [sym, pos[0], pos[1], pos[2]]
        return obmol, structure, bond
    else:  # use RDKit
        rdmol = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(rdmol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(rdmol)
        atoms = rdmol.GetAtoms()
        bond = np.zeros((len(atoms), len(atoms)), dtype=int)
        for i in range(len(rdmol.GetAtoms())):
            for j in range(len(rdmol.GetAtoms())):
                if not rdmol.GetBondBetweenAtoms(i, j) is None:
                    b = rdmol.GetBondBetweenAtoms(i, j)
                    order = int(b.GetBondTypeAsDouble())
                    bond[i][j] = order
        for i, atom in enumerate(rdmol.GetAtoms()):
            pos = rdmol.GetConformer(0).GetAtomPosition(i)
            sym = atom.GetSymbol()
            structure += [sym, pos.x, pos.y, pos.z]
        return rdmol, structure, bond


def create_ob_mol(smi):
    """
    Method to create a Molecule Object from ObenBabel
    """
    obmol = pybel.readstring('smi', smi)
    obmol.OBMol.AddHydrogens()
    return obmol


def create_rdkit_mol(bond, atom):
    """
    Method to create a RDKit Molecule object from a KinBot stationary_pt object
    """
    m = Chem.MolFromSmiles('[' + atom[0] + ']')
    mw = Chem.RWMol(m)
    for i in range(1, len(atom)):
        dummy = Chem.MolFromSmiles('[' + atom[i] + ']')
        at = dummy.GetAtoms()[0]
        #at = Chem.Atom(syms_to_num[atom[i]])
        #at.SetNoImplicit(True)
        mw.AddAtom(at)
    for i in range(len(atom)-1):
        for j in range(i, len(atom)):
            if bond[i][j] == 1:
                mw.AddBond(i, j, Chem.BondType.SINGLE)
            if bond[i][j] == 2:
                mw.AddBond(i, j, Chem.BondType.DOUBLE)
            if bond[i][j] == 3:
                mw.AddBond(i, j, Chem.BondType.TRIPLE)
    smi = Chem.MolToSmiles(mw)
    inchi = Chem.MolToInchi(mw)
    return mw, smi, inchi


def create_inchi_from_geom(atom, geom):
    xyz_file = 'temp.xyz'
    f = open(xyz_file, 'w')
    f.write(str(len(atom)) + '\n\n')
    for i, at in enumerate(atom):
        x, y, z = geom[i]
        f.write('{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z))
    f.write('\n\n')
    f.close()
    inchi = create_inchi('', '', xyz_file=xyz_file)
    # remove temp file
    os.remove(xyz_file)
    return inchi


def create_inchi(job, chemid, xyz_file=''):
    if xyz_file == '':
        xyz_file = os.path.expanduser(job) + 'xyz/' + chemid + '.xyz'
    obmol = pybel.readfile('xyz', xyz_file).next()
    return obmol.write("inchi", opt={'T': 'nostereo'}).split()[0]


def create_inchi_from_smi(smi):
    """
    Method to create the InChI of a structure given its smiles.
    OpenBabel is used for this.
    """
    obmol = pybel.readstring('smi', smi)
    return obmol.write("inchi").split()[0]


def create_smiles(inchi):
    """
    Method to create the smiles of a structure given its InChI.
    OpenBabel is used for this.
    """
    obmol = pybel.readstring('inchi', inchi)
    return obmol.write("smi").split()[0]
