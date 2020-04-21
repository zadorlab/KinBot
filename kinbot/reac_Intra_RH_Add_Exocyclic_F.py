from kinbot.reac_General import GeneralReac
from kinbot import geometry

class IntraRHAddExoF(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step)

        elif step < self.max_step:
            self.release_dihedrals(release)

            final_dist = [2.0,1.9,1.4,1.45,1.3]
            val = geometry.new_bond_length(self.species,self.instance[1],self.instance[-2],step - 11,10,final_dist[0],geom)
            self.set_bond(1, -2, val, change)
            
            val = geometry.new_bond_length(self.species,self.instance[-1],self.instance[-2],step - 11,10,final_dist[1],geom)
            self.set_bond(-1, -2, val, change)
            
            val = geometry.new_bond_length(self.species,self.instance[0],self.instance[-1],step - 11,10,final_dist[2],geom)
            self.set_bond(0, -1, val, change)
            
            val = geometry.new_bond_length(self.species,self.instance[0],self.instance[1],step - 11,10,final_dist[3],geom)
            self.set_bond(0, 1, val, change)
            
            val = geometry.new_bond_length(self.species,self.instance[-1],self.instance[1],step - 11,10,final_dist[4],geom)
            self.set_bond(-1, 1, val, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
