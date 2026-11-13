"""Bounded fragment charge and multiplicity choices for dissociation products."""

import copy
import itertools

import numpy as np

from kinbot.reaction_path import endpoint_snapshot
from kinbot.species_routing import routing_name
from kinbot.stereo_identity import require_supported_identity


def spin_allowed(multiplicities, total):
    """Add integer 2S values, including all couplings of three or more spins."""
    spins = {0}
    for mult in multiplicities:
        spin = int(mult) - 1
        spins = {value for previous in spins
                 for value in range(abs(previous - spin), previous + spin + 1, 2)}
    return int(total) - 1 in spins


def multiplicities(fragment, excite=False):
    """Try valence-supported spins, optionally promoting one bonded pair.

    Valence remainders include carbene/biradical alternatives but miss excited
    triplet ethene. A multiple-bond promotion is considered only when the first
    set of fragment states cannot carry the required total spin.
    """
    parity = int(fragment.nel) % 2
    counts = []
    for bond, rad in zip(fragment.bonds, fragment.rads):
        values = np.maximum(rad, 0).copy()
        for atom, element in enumerate(fragment.atom):
            if element == 'N':
                # Characterization removes N valence remainders; atomic N and
                # NH must still admit quartet and triplet states respectively.
                values[atom] = max(values[atom], 3 - np.sum(bond[atom]))
        counts.append(int(values.sum()))
    unpaired = max(counts, default=parity)
    if fragment.natom == 1 and fragment.atom[0] == 'S':
        unpaired = 2  # The graph's sulfur valence of six is not six unpaired electrons.
    # Charge changes the electron count without changing the valence graph.
    # For example O+ can carry three unpaired electrons, whereas O has two.
    maximum = max(unpaired + abs(int(fragment.charge)), int(fragment.mult) - 1)
    if excite and any(np.any(np.asarray(bond) > 1) for bond in fragment.bonds):
        maximum += 2
    return list(range(parity + 1, min(maximum, int(fragment.nel)) + 2, 2))


def spin_copy(fragment, mult):
    point = copy.copy(fragment)
    for field in ('atom', 'geom', 'bond', 'bond01', 'bonds', 'rads'):
        setattr(point, field, copy.deepcopy(getattr(fragment, field)))
    if hasattr(fragment, 'hydrogen_bridges'):
        point.hydrogen_bridges = copy.deepcopy(fragment.hydrogen_bridges)
    point.mult = int(mult)
    point.__dict__.pop('optical_reference', None)
    point.calc_chemid()
    point.name = routing_name(point)
    require_supported_identity(point)
    return point


def fragment_state_combinations(endpoint):
    """Enumerate fragment states that can carry the parent's charge and spin.

    Atom-index groups keep alternative charge or spin assignments of one
    fragment from replacing a different physical fragment. A multiple-bond
    excitation is added only when the first bounded state set cannot carry
    the parent's spin, matching ordinary bimolecular product selection.
    """
    fragments, maps = endpoint_snapshot(endpoint).start_multi_molecular(vary_charge=True)
    if len(fragments) == 1:
        return fragments, [tuple(range(endpoint.natom))], [(0,)], False
    slots = {}
    for fragment, atoms in zip(fragments, maps):
        slots.setdefault(tuple(map(int, atoms)), []).append(fragment)
    for excite in (False, True):
        states, groups = [], []
        for variants in slots.values():
            group = []
            for fragment in variants:
                for mult in multiplicities(fragment, excite):
                    group.append(len(states))
                    states.append(spin_copy(fragment, mult))
            groups.append(group)
        combinations = [comb for comb in itertools.product(*groups)
            if sum(states[i].charge for i in comb) == endpoint.charge
            and spin_allowed([states[i].mult for i in comb], endpoint.mult)]
        if combinations:
            return states, list(slots), combinations, excite
    raise ValueError('No tested fragment multiplicity combination can carry '
                     f'total multiplicity {endpoint.mult} and charge {endpoint.charge}.')


def participating_fragment_states(endpoint):
    """Return only states occurring in a charge- and spin-compatible choice."""
    states, _, combinations, _ = fragment_state_combinations(endpoint)
    active = sorted(set(itertools.chain.from_iterable(combinations)))
    return [states[i] for i in active]
