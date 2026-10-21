"""Bimolecular radical abstraction reaction family."""

import math

import numpy as np
from ase import Atoms
from ase.data import atomic_numbers, covalent_radii

from kinbot import geometry
from kinbot.bimol import (GOLDEN_ANGLE, approach_directions, open_directions,
                          rotate_between, unit_vector)
from kinbot.reactions.reac_bimol import BimolecularReaction


def _roles(species, instance):
    """Map combined abstraction indices onto donor and abstractor fragments."""
    size_a = species.fragA.natom

    def locate(index):
        return (0, index) if index < size_a else (1, index - size_a)

    abstractor = locate(instance[0])
    hydrogen = locate(instance[1])
    donor = locate(instance[2])
    if hydrogen[0] != donor[0] or abstractor[0] == donor[0]:
        raise ValueError('Invalid bimolecular abstraction atom assignment.')
    return abstractor, hydrogen, donor


def generate_encounters(species, geom_a, geom_b, instance, settings):
    """Place the abstractor around the donor-H approach axis."""
    fragments = [species.fragA, species.fragB]
    geometries = [np.asarray(geom_a), np.asarray(geom_b)]
    abstractor, hydrogen, donor = _roles(species, instance)
    donor_geom = geometries[donor[0]]
    abstractor_geom = geometries[abstractor[0]]
    dh_axis = unit_vector(
        donor_geom[hydrogen[1]] - donor_geom[donor[1]])
    element = fragments[abstractor[0]].atom[abstractor[1]]
    distance = (covalent_radii[atomic_numbers[element]]
                + covalent_radii[atomic_numbers['H']]
                + settings['encounter_extra_distance'])
    pose = 0
    for open_direction in open_directions(
            fragments[abstractor[0]], abstractor_geom, abstractor[1]):
        for approach in approach_directions(dh_axis, settings):
            radical_position = donor_geom[hydrogen[1]] + distance * approach
            incoming = -approach
            centered = abstractor_geom - abstractor_geom[abstractor[1]]
            centered = centered @ rotate_between(open_direction, incoming).T
            centered = centered @ geometry.rotation_matrix(
                incoming, pose * GOLDEN_ANGLE).T
            combined = [None, None]
            combined[donor[0]] = donor_geom
            combined[abstractor[0]] = centered + radical_position
            yield Atoms(species.atom, positions=np.vstack(combined)), {
                'pose': pose,
                'angle': math.degrees(math.acos(np.clip(
                    np.dot(-dh_axis, approach), -1., 1.))),
                'distance': float(distance),
            }
            pose += 1


class Abstraction(BimolecularReaction):
    """A-H + B(radical) -> A(radical) + B-H.

    ``instance`` is ``[abstractor, H, donor]`` in the combined KinBot atom
    order. The generic bimolecular engine clusters L0 saddles, validates each
    unique representative by L0 IRC, then hands it to L1 refinement.
    """

    family_name = 'abstraction'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.do_reac_vdW = False

    def get_bimol_search_data(self):
        """Return abstraction-specific placement and connectivity data."""
        return {
            'family_module': __name__,
            'encounter_function': 'generate_encounters',
            'formed_bonds': [self.instance[:2]],
            'broken_bonds': [self.instance[1:3]],
        }

    def get_bimol_site(self, instance):
        """Use abstractor and transferred-H classes for degeneracy."""
        return tuple(instance[:2])

    def prepare_hir_rotors(self):
        """Replace transferred-H dihedrals using the observed IRC reaction."""
        reactant = np.asarray(self.irc_reactant.bond01)
        product = np.asarray(self.irc_prod.bond01)
        formed = []
        broken = []
        for first in range(self.ts.natom):
            for second in range(first + 1, self.ts.natom):
                change = product[first, second] - reactant[first, second]
                if change > 0:
                    formed.append((first, second))
                elif change < 0:
                    broken.append((first, second))

        if len(formed) != 1 or len(broken) != 1:
            return None
        transferred = [atom for atom in set(formed[0]) & set(broken[0])
                       if self.ts.atom[atom] == 'H']
        if len(transferred) != 1:
            return None
        hydrogen = transferred[0]
        donor = next(atom for atom in broken[0] if atom != hydrogen)
        abstractor = next(atom for atom in formed[0] if atom != hydrogen)
        if donor == abstractor:
            return None
        donor_neighbors = sorted(
            int(atom) for atom in np.flatnonzero(reactant[donor])
            if atom != hydrogen)
        abstractor_neighbors = sorted(
            int(atom) for atom in np.flatnonzero(reactant[abstractor])
            if atom != hydrogen)
        replacement = None
        for donor_neighbor in donor_neighbors:
            for abstractor_neighbor in abstractor_neighbors:
                rotor = [donor_neighbor, donor, abstractor,
                         abstractor_neighbor]
                if len(set(rotor)) < 4:
                    continue
                _, warning = geometry.calc_dihedral(
                    *(self.ts.geom[atom] for atom in rotor))
                if not warning:
                    replacement = rotor
                    break
            if replacement is not None:
                break

        def replace_terminal_hydrogen(rotor):
            rotor = list(rotor)
            if hydrogen not in rotor:
                return rotor
            position = rotor.index(hydrogen)
            if position not in (0, 3):
                return None
            center = rotor[1] if position == 0 else rotor[2]
            other_center = rotor[2] if position == 0 else rotor[1]
            candidates = [int(atom) for atom in np.flatnonzero(
                reactant[center])
                if atom not in {hydrogen, other_center}]
            candidates.sort(key=lambda atom: (
                self.ts.atom[atom] != 'H', atom))
            for candidate in candidates:
                trial = rotor[:]
                trial[position] = candidate
                if len(set(trial)) < 4:
                    continue
                _, warning = geometry.calc_dihedral(
                    *(self.ts.geom[atom] for atom in trial))
                if not warning:
                    return trial
            return None

        def clean_rotors(rotors):
            cleaned = []
            axes = set()
            for rotor in rotors:
                rotor = replace_terminal_hydrogen(rotor)
                if rotor is None:
                    continue
                axis = tuple(sorted(rotor[1:3]))
                if axis not in axes:
                    cleaned.append(rotor)
                    axes.add(axis)
            return cleaned

        original_dihed = [list(rotor) for rotor in self.ts.dihed]
        self.ts.dihed = clean_rotors(self.ts.dihed)
        self.ts.conf_dihed = clean_rotors(self.ts.conf_dihed)
        axes = {tuple(sorted(rotor[1:3])) for rotor in self.ts.dihed}
        if (replacement is not None
                and tuple(sorted(replacement[1:3])) not in axes):
            self.ts.dihed.append(replacement)
        terminal_replacements = [rotor for rotor in self.ts.dihed
                                 if rotor not in original_dihed]
        self.hir_transfer = [abstractor, hydrogen, donor]
        self.hir_rotors = terminal_replacements[:]
        if replacement is not None and replacement not in self.hir_rotors:
            self.hir_rotors.append(replacement)
        self.hir_rotor = (replacement if replacement is not None else
                          (terminal_replacements[0]
                           if terminal_replacements else None))
        return self.hir_rotor
