import sys
import os
import logging

import numpy as np
from PIL import Image

from kinbot import kb_path

# try to import pybel
try:
    import pybel
    pybel.ob.obErrorLog.SetOutputLevel(0)
except ImportError:
    try:
        from openbabel import pybel
        pybel.ob.obErrorLog.SetOutputLevel(0)
    except:
        print('Warning: Pybel could not be imported.')
        print('Certain features or the whole code might not run properly.')
        pass

logger = logging.getLogger('KinBot')

num_to_syms = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 16: 'S'}
syms_to_num = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'S': 16}


def get_molecular_formula(smi):
    """
    Return the molecular formula of the molecule corresponding to the smiles
    """
    try:
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    except NameError:
        logger.error('RDKit is not installed or loaded correctly.')
        sys.exit()
    return rdMolDescriptors.CalcMolFormula(mol)


def create_rxn_depiction(react_smiles, prod_smiles, cdir, name):
    """
    Create a 2D depiction of a chemical reaction,
    react smiles: smiles of the reactants
    prod_smiles: smiles of the products
    dir: directory where to save the reaction depiction
    name: name of the depiction files
    """
    react_png = f'{cdir}/react.png'
    prod_png = f'{cdir}/prod.png'

    try:
        obmol = pybel.readstring("smi", react_smiles)
    except NameError:
        logger.error('Cannot create 2D structures, Pybel is not loaded or installed properly.')
        sys.exit()

    obmol.draw(show=False, filename=react_png)

    obmol = pybel.readstring("smi", prod_smiles)
    obmol.draw(show=False, filename=prod_png)

    arrow = f'{kb_path}/tpl/arrow.png'
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

    new_im.save(f'{cdir}/{name}.png')


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
                    try:
                        order = obmol.OBMol.GetBond(i+1, j+1).GetBO()
                    except:
                        try:
                            order = obmol.OBMol.GetBond(i+1, j+1).GetBondOrder()
                        except:
                            logger.error('Something went wrong with OpenBabel')
                            sys.exit()
                    bond[i][j] = order
        for at in obmol.atoms:
            pos = at.coords
            sym = num_to_syms[at.atomicnum]
            structure += [sym, pos[0], pos[1], pos[2]]
        return obmol, structure, bond
    else:  # use RDKit
        try:
            rdmol = Chem.AddHs(Chem.MolFromSmiles(smi))
        except NameError:
            logger.error('RDKit is not installed or loaded correctly.')
            sys.exit()
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
    try:
        obmol = pybel.readstring('smi', smi)
    except NameError:
        logger.error('Pybel is not installed or loaded correctly.')
        sys.exit()
    obmol.OBMol.AddHydrogens()
    return obmol


def create_rdkit_mol(bond, atom):
    """
    Method to create a RDKit Molecule object from a KinBot stationary_pt object
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem import rdMolDescriptors
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
    except ImportError:
        logger.warning('RDKit could not be imported.')
        pass

    try:
        m = Chem.MolFromSmiles('[' + atom[0] + ']')
    except NameError:
        logger.error('RDKit is not installed or loaded correctly.')
        sys.exit()

    mw = Chem.RWMol(m)
    for i in range(1, len(atom)):
        dummy = Chem.MolFromSmiles('[' + atom[i] + ']')
        at = dummy.GetAtoms()[0]
        # at = Chem.Atom(syms_to_num[atom[i]])
        # at.SetNoImplicit(True)
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
    return mw, smi


def create_inchi_from_geom(atom, geom):
    xyz_file = 'temp.xyz'
    with open(xyz_file, 'w') as f:
        f.write(str(len(atom)) + '\n\n')
        for i, at in enumerate(atom):
            x, y, z = geom[i]
            f.write('{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z))
        f.write('\n\n')
    inchi = create_inchi('', '', xyz_file=xyz_file)
    # remove temp file
    os.remove(xyz_file)
    return inchi


def create_inchi(job, chemid, xyz_file=''):
    if xyz_file == '':
        xyz_file = os.path.expanduser(job) + 'xyz/' + chemid + '.xyz'
    obmol = list(pybel.readfile('xyz', xyz_file))[0]
    try:
        obmol = list(pybel.readfile('xyz', xyz_file))[0]
    except NameError:
        logger.error('Pybel is not installed or loaded correctly.')
        sys.exit()

    # return obmol.write("inchi", opt={'T': 'nostereo'}).split()[0]
    return obmol.write("inchi").split()[0]


def create_inchi_from_smi(smi):
    """
    Method to create the InChI of a structure given its smiles.
    OpenBabel is used for this.
    """
    try:
        obmol = pybel.readstring('smi', smi)
    except NameError:
        logger.error('Pybel is not installed or loaded correctly.')
        sys.exit()

    return obmol.write("inchi").split()[0]


def create_smiles(inchi):
    """
    Method to create the smiles of a structure given its InChI.
    OpenBabel is used for this.
    """
    try:
        obmol = pybel.readstring('inchi', inchi)
    except NameError:
        logger.error('Pybel is not installed or loaded correctly.')
        sys.exit()

    return obmol.write("smi").split()[0]


def create_smi_from_geom(atom, geom):
    inchi = create_inchi_from_geom(atom, geom)
    return create_smiles(inchi)
