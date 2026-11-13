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
        changes = self.product_bond_changes()
        if (getattr(self.species, 'hydrogen_bridges', ())
                and self.family_name == 'h2_elim_barrierless'):
            # Localize only the bridge contacts broken by this reaction.
            # Otherwise its second contact leaves the intended H2 product
            # attached to the residue. Unaffected bridges retain both views.
            product.bonds = [matrix.copy() for matrix in self.species.bonds
                if all(matrix[first, second] > 0 for first, second, order in changes if order == 0)]
            if not product.bonds:
                raise ValueError('No bridge attachment supports the requested bond cleavage.')
            product.bond = np.maximum.reduce(product.bonds)
        for first, second, order in changes:
            matrices = ([product.bond, *product.bonds] if getattr(self.species, 'hydrogen_bridges', ())
                        else (product.bond, product.bonds[0]))
            for matrix in matrices:
                matrix[first][second] = matrix[second][first] = order
        if getattr(self.species, 'hydrogen_bridges', ()):
            product.hydrogen_bridges = [bridge for bridge in product.hydrogen_bridges
                if product.bond[bridge[0], bridge[1]] and product.bond[bridge[1], bridge[2]]]
            product.bond01 = (product.bond > 0).astype(int)
        return product
    
    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
