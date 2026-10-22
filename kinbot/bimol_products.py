"""Compare fragment spin states without changing the IRC's total spin."""

import copy
import itertools
import logging

import numpy as np

from kinbot.bimol import load_selected_species_properties
from kinbot.calculation import load_calculation_record
from kinbot.optimize import Optimize
from kinbot.reaction_path import endpoint_snapshot
from kinbot.species_routing import routing_name, same_species, reusable_cached_product
from kinbot.stereo_identity import require_supported_identity

logger = logging.getLogger('KinBot')


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
    point.mult = int(mult)
    point.__dict__.pop('optical_reference', None)
    point.calc_chemid()
    point.name = routing_name(point)
    require_supported_identity(point)
    return point


def check_product_spins(products, endpoint):
    """Refuse products that cannot belong to the IRC's charge and spin."""
    if (sum(point.charge for point in products) != endpoint.charge
            or not spin_allowed([point.mult for point in products], endpoint.mult)):
        raise ValueError('Final product charges or multiplicities cannot reproduce '
                         f'IRC charge {endpoint.charge} and multiplicity {endpoint.mult}.')


class ProductSpinSelection:
    """Poll ordinary QC jobs for the bounded set of allowed product states.

    Only this in-memory calculation progress is new. Results and restart
    checks use the existing ASE database, theory routes, and product jobs.
    """

    def __init__(self, endpoint):
        self.total = int(endpoint.mult)
        self.charge = int(endpoint.charge)
        fragments, maps = endpoint_snapshot(endpoint).start_multi_molecular(vary_charge=True)
        self.needed = False
        if len(fragments) == 1:
            return  # A connected product already has the IRC multiplicity.
        slots = {}
        for fragment, atoms in zip(fragments, maps):
            slots.setdefault(tuple(map(int, atoms)), []).append(fragment)
        self.maps = list(slots)
        self.states = []
        self.combinations = []
        for excite in (False, True):
            self.states = []
            groups = []
            for variants in slots.values():
                group = []
                for fragment in variants:
                    for mult in multiplicities(fragment, excite):
                        group.append(len(self.states))
                        self.states.append(spin_copy(fragment, mult))
                groups.append(group)
            self.combinations = [comb for comb in itertools.product(*groups)
                if sum(self.states[i].charge for i in comb) == self.charge
                and spin_allowed([self.states[i].mult for i in comb], self.total)]
            if self.combinations:
                break
        if not self.combinations:
            raise ValueError('No tested fragment multiplicity combination can carry '
                             f'total multiplicity {self.total} and charge {self.charge}.')
        chosen = self.combinations[0]
        # A single allowed charge placement still needs this selection: the
        # ordinary charge ranking also considers spin-incompatible placements.
        self.needed = (len(self.combinations) > 1
            or any(len(originals) != 1 for originals in slots.values()) or any(
            not any(same_species(self.states[index], original) for original in originals)
            for index, originals in zip(chosen, slots.values())))
        self.active = sorted(set(itertools.chain.from_iterable(self.combinations)))
        names = {}
        for index in self.active:
            point = self.states[index]
            key = routing_name(point)
            if key in names and not same_species(point, names[key]):
                raise ValueError(f'{key}: distinct product charge/spin states have the '
                                 'same ordinary job name; cannot safely compare them.')
            names[key] = point
        self.sp_jobs = {}
        self.sp_energies = {}
        self.optimizations = {}
        self.selected = None
        self.evidence = {
            'policy': 'lowest E+ZPE among tested spin-compatible fragment minima',
            'total_multiplicity': self.total,
            'fragment_atom_indices': [list(atoms) for atoms in self.maps],
            'multiple_bond_excitation_added': excite,
            'electronic_state_connection': 'not verified; spin coupling and connectivity only',
        }

    def advance(self, qc, par, other_optimizations=()):
        """Return True only after all admissible contenders have been compared."""
        if self.selected is not None:
            return True
        for level in ('l1', 'l2') if par['high_level'] else ('l1',):
            backend = qc.for_level(level)
            if any(key.lower() in ('mult', 'multiplicity', 'charge')
                   for key in backend.par.get('calc_kwargs', {})):
                raise ValueError('Product spin comparison cannot use calculator_kwargs '
                                 'that override fragment charge or multiplicity.')
        pending = False
        for index in self.active:
            point = self.states[index]
            if index not in self.sp_jobs:
                self.sp_jobs[index] = qc.qc_product_spin(point, point.geom)
            job = self.sp_jobs[index]
            status = qc.check_qc(job)
            if status == 'error':
                raise ValueError(f'{job}: product spin single point failed; '
                                 'the state comparison remains unresolved.')
            if status != 'normal':
                pending = True
                continue
            error, energy = qc.get_qc_energy(job)
            if error or not np.isfinite(energy):
                raise ValueError(f'{job}: no finite product spin energy.')
            self.sp_energies[index] = float(energy)
        if pending:
            return False

        # Do not discard a state solely because its frozen-geometry energy is
        # high: relaxation can reverse the single-point ordering.
        for index in sorted(self.active, key=self.sp_energies.get):
            reference = self.states[index]
            key = routing_name(reference)
            if key not in self.optimizations:
                for opt in [*self.optimizations.values(), *other_optimizations]:
                    if same_species(reference, opt.species):
                        self.optimizations[key] = opt
                        break
            if key not in self.optimizations:
                point = reusable_cached_product(qc, spin_copy(reference, reference.mult))
                qc.qc_opt(point, point.geom)
                job = routing_name(point) + '_well'
                status = qc.check_qc(job)
                if status == 'error':
                    raise ValueError(f'{job}: product spin optimization failed; '
                                     'the state comparison remains unresolved.')
                if status != 'normal':
                    pending = True
                    continue
                load_calculation_record(point, qc, job)
                self._check_product(reference, point)
                self.optimizations[key] = Optimize(point, par, qc)
            opt = self.optimizations[key]
            if opt.shir != 1 and opt.shigh != -999:
                opt.do_optimization()
            if opt.shigh == -999:
                raise ValueError(f'{key}: product spin refinement failed; '
                                 'the state comparison remains unresolved.')
            if opt.shir != 1:
                pending = True
                continue
            load_selected_species_properties(opt, qc)
            self._check_product(reference, opt.species)
        if pending:
            return False
        energies = {index: self.optimizations[routing_name(self.states[index])].species.energy
                    + self.optimizations[routing_name(self.states[index])].species.zpe
                    for index in self.active}
        if not all(np.isfinite(value) for value in energies.values()):
            raise ValueError('Product spin comparison has a nonfinite optimized energy.')
        self.selected = min(self.combinations, key=lambda comb: sum(energies[i] for i in comb))
        self.evidence.update(
            level='l2' if par['high_level'] else 'l1',
            states=[{'index': i, 'species': routing_name(self.states[i]),
                     'charge': int(self.states[i].charge), 'multiplicity': self.states[i].mult,
                     'single_point_job': self.sp_jobs[i],
                     'single_point_energy_hartree': self.sp_energies[i],
                     'selected_job': self.optimizations[routing_name(self.states[i])].selected_job,
                     'energy_with_zpe_hartree': energies[i]} for i in self.active],
            allowed_combinations=[list(comb) for comb in self.combinations],
            selected_combination=list(self.selected))
        logger.info('Selected product multiplicities %s at total multiplicity %s; '
                    'lowest E+ZPE among tested spin-compatible states.',
                    [self.states[i].mult for i in self.selected], self.total)
        return True

    @staticmethod
    def _check_product(reference, product):
        observed = spin_copy(product, product.mult)
        observed.rads = []
        observed.characterize()
        observed.calc_chemid()
        if not same_species(reference, observed):
            raise ValueError(f'{routing_name(reference)}: product optimization changed '
                             'connectivity or configuration; spin comparison unresolved.')

    def selected_optimizations(self):
        return [self.optimizations[routing_name(self.states[i])] for i in self.selected]
