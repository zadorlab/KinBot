from ase.db import connect
from stationary_pt import StationaryPoint
from fragments import Fragment

class vdW_Well(StationaryPoint):
    frag_number = 0 #Reset after creation of a vdW_Well object
    def __init__(self, parent, fragA, fragB, par):
        """
        The parent should be a StationaryPoint object, and fragA and fragB should be chemids.
        Create 2 Fragment objects for the fragments.
        The fragments need to be Fragments object to be linked with the VRC-TST treatment of barrierless reactions.
        """
        self.parent = parent #Stationary point of irc_prod
        self.fragA = self.create_Fragment(fragA, par)
        self.fragB = self.create_Fragment(fragB, par)

    def create_Fragment(self, chemid, par):
       if par['high_level']:
           #Will read info from L2 structure
           basename = f"{chemid}_well_high"
       else:
           #Will read info from L1 structure
           basename = f"{chemid}_well"

        #Create ase.atoms objects for each fragments
        #db = connect(f"{self.parent.chemid}/kinbot.db")
        #*_, last_row = db.select(name=f"{basename}", sort="-1")
        #atoms = last_row.toatoms() #This is an ase.atoms object
        #fragment = Fragment.from_ase_atoms(atoms=atoms,
        #                                   frag_number=frag_number,
        #                                   max_frag=2,
        #                                   chemid=chemid,
        #                                   parent_chemid=self.parent.chemid,
        #                                       ))
        #frag_number +=1
        #if frag_number == 2:
        #    frag_number = 0

        #return fragment
