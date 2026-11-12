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
        # The parent H-H distance is unsuitable for an isolated H2 minimum
        # calculation. Only the separated fragments will be optimized.
        axis = product.geom[second] - product.geom[first]
        axis /= np.linalg.norm(axis)
        midpoint = (product.geom[first] + product.geom[second]) / 2.
        product.geom[first] = midpoint - .37 * axis
        product.geom[second] = midpoint + .37 * axis
        product.bond01 = (product.bond > 0).astype(int)
        return product
