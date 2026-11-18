"""Cation 1,1-H2 elimination using ordinary dissociation product handling."""

import numpy as np

from kinbot.product_states import InvalidDissociation
from kinbot.reactions.reac_homolytic_scission import HS


class H2ElimBarrierless(HS):
    family_name = 'h2_elim_barrierless'

    def product_bond_changes(self):
        first, carbon, second = self.instance
        return [(first, carbon, 0), (carbon, second, 0), (first, second, 1)]

    def edit_product_bonds(self, product):
        bridges = getattr(self.species, 'hydrogen_bridges', ())
        if not bridges:
            return super().edit_product_bonds(product)
        first, carbon, second = self.instance
        # Only localized bond assignments with both H attached to this carbon
        # can lose this H2. Removing a bridge returns its electron pair to C–C.
        carbons = np.asarray(product.atom) == 'C'
        bonds = []
        for original in self.species.bonds:
            if not original[first, carbon] or not original[carbon, second]:
                continue
            matrix = original.copy()
            for a, b, order in self.product_bond_changes():
                matrix[a, b] = matrix[b, a] = order
            for a, hydrogen, b in bridges:
                if hydrogen in (first, second):
                    matrix[a, b] += 1
                    matrix[b, a] += 1
            if all(matrix[carbons].sum(axis=1) <= 4):
                bonds.append(matrix)
        if not bonds:
            raise InvalidDissociation('No localized carbon bonding permits loss of this H pair.')
        product.bonds = bonds
        product.bond = np.maximum.reduce(bonds)

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
