from kinbot.reac_General import GeneralReac
from kinbot import geometry
from kinbot.stationary_pt import StationaryPoint
import numpy as np

class HS(GeneralReac):
    scan = 0
    skip = 0
    max_step = 0
    family_name = 'hs'

    def product_bond_changes(self):
        return [(self.instance[0], self.instance[1], 0)]

    def get_barrierless_product(self):
        """Prepare the disconnected product graph before fragment optimization."""
        product = StationaryPoint(f'{self.instance_name}_prod', self.species.charge,
            self.species.mult, atom=self.species.atom, geom=self.species.geom.copy())
        product.characterize()
        for first, second, order in self.product_bond_changes():
            for matrix in (product.bond, product.bonds[0]):
                matrix[first][second] = matrix[second][first] = order
        return product
    
    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
