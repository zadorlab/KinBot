"""Bond-order eligibility shared by rotor detection and rotor symmetry."""


def single_bond_candidate(atoms, bonds, left, right):
    """Keep ordinary all-single axes plus sulfur mixed-valence candidates.

    A sulfur axis still needs a single-bond form and unchanged connectivity.
    This is permission to scan, not evidence of a low torsional barrier. The
    caller separately checks rings, terminal atoms and linear references.
    """
    orders = [bond[left][right] for bond in bonds]
    return (bool(orders) and all(order > 0 for order in orders) and 1 in orders
            and (all(order == 1 for order in orders)
                 or 'S' in (atoms[left], atoms[right])))
