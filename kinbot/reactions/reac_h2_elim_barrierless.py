"""Declared 1,1-H2 elimination using ordinary dissociation product handling."""

import numpy as np

from kinbot.reactions.reac_homolytic_scission import HS


class H2ElimBarrierless(HS):
    family_name = 'h2_elim_barrierless'

    def product_bond_changes(self):
        first, carbon, second = self.instance
        return [(first, carbon, 0), (carbon, second, 0), (first, second, 1)]

    def get_barrierless_product(self):
        product = super().get_barrierless_product()
        first, _, second = self.instance
        for carbon1, hydrogen, carbon2 in getattr(self.species, 'hydrogen_bridges', ()):
            if hydrogen in (first, second):
                # Removing a shared proton in H2 leaves the three-centre
                # electron pair on C--C. Without this edit, protonated ethane
                # incorrectly becomes CH2 + CH3 + H2 instead of C2H5 + H2.
                for matrix in (product.bond, *product.bonds):
                    matrix[carbon1, carbon2] += 1
                    matrix[carbon2, carbon1] += 1
        if getattr(self.species, 'hydrogen_bridges', ()):
            if any(np.sum(matrix[i]) > 4 for matrix in product.bonds
                   for i, element in enumerate(product.atom) if element == 'C'):
                raise ValueError('H2 loss from this bridge exceeds carbon valence.')
        # The parent H-H distance is unsuitable for an isolated H2 minimum
        # calculation. Only the separated fragments will be optimized.
        axis = product.geom[second] - product.geom[first]
        axis /= np.linalg.norm(axis)
        midpoint = (product.geom[first] + product.geom[second]) / 2.
        product.geom[first] = midpoint - .37 * axis
        product.geom[second] = midpoint + .37 * axis
        product.bond01 = (product.bond > 0).astype(int)
        return product
