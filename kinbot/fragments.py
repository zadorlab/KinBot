from ase.db import connect
from ase import Atoms

class Fragment(Atoms):
    """
    Class that creates fragments in the form of Atom objects from the ASE package.

    """

    #def __init__(self, symbols=None, positions=None, numbers=None,tags=None,
    #             momenta=None, masses=None, magmoms=None, charges=None,
    #             scaled_positions=None, cell=None, pbc=None, celldisp=None,
    #             constraint=None, calculator=None, info=None):

       # Atoms.__init__(self, symbols, positions, numbers,
       #                tags, momenta, masses, magmoms, charges,
       #                scaled_positions, cell, pbc, celldisp,
       #                constraint, calculator, info)

    def __init__(self, symbols=None, positions=None, masses=None)

        Atoms.__init__(self, symbols, positions, masses)

    def create(self,par,chemid,parents_chemid)):
        #product_list: list of the fragments' chemid
        #parents: dictionary
        #par: json file parameters

        self.chemid = chemid

        #Identify the parent of the reaction from parent dictionary
        #self.parent_chemid = parent.get('_'.join(sorted(products_list)))
        self.parent_chemid = parent_chemid

        if par['high_level']:
            #Will read info from L2 structure
            self.basename = '{self.parent_chemid}/{self.chemid}_well_high'
        else
            #Will read info from L1 structure
            self.basename = '{self.parent_chemid}/{self.chemid}_well'

        #Recover Symbols, positions and masses from the database
        db = '{self.parent_chemid}/kinbot.db'
        for row in db.select(name='{self.basename}'):
            self.symbols = row.toatoms().get_chemical_symbols()
            self.positions = row.toatoms().get_positions()
            self.masses = row.toatoms().get_masses()

        return self
