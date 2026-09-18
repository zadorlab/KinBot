"""Connectivity-based hierarchy reactions and 0 K formation enthalpies.

This implementation currently handles neutral, closed-shell molecules with
ordinary integral bond orders. It constructs capped graph fragments from the
same explicit-hydrogen SMILES connectivity that KinBot uses for its species.
Every generated reaction is checked for elemental balance. Radical/ionic and
aromatic-state extensions require state-aware fragment rules and fail closed.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
import math
from typing import Mapping

import networkx as nx
from ase.units import Hartree, kJ, mol

from kinbot.anl.atct import (ATcTRecord, ATcTTable, _canonical_smiles,
                            _preferred_formula_atoms)
from kinbot.energy import ZeroKEnergy


HARTREE_TO_KJ_MOL = Hartree * mol / kJ
_VALENCE = {6: 4, 7: 3, 8: 2, 9: 1, 17: 1}
_H2 = '[H][H]'
DEFAULT_METHOD_LADDER = ('ANL1-F12', 'ANL1', 'ANL0-F12', 'ANL0',
                         'L3', 'L2')


@dataclass(frozen=True)
class CBHReaction:
    rung: int
    target_smiles: str
    # Signed coefficients: negative for reactants, positive for products.
    stoichiometry: Mapping[str, int]
    formulas: Mapping[str, Mapping[str, int]]
    # Optional state constraints for externally specified radical reactions.
    # The graph generator currently creates neutral singlet references only.
    states: Mapping[str, tuple[int, int]] = field(default_factory=dict)

    @property
    def reference_smiles(self) -> tuple[str, ...]:
        return tuple(sorted(item for item in self.stoichiometry
                            if item != self.target_smiles))


@dataclass(frozen=True)
class FormationEnthalpy:
    target_smiles: str
    rung: int
    method: str
    reaction_energy_0k_kj_mol: float
    formation_0k_kj_mol: float
    references: Mapping[str, ATcTRecord]
    atct_version: str
    atct_source_sha256: str
    energy_sources: Mapping[str, str]


@dataclass(frozen=True)
class LadderSelection:
    formation: FormationEnthalpy
    reaction: CBHReaction
    skipped: tuple[str, ...]


def _graph_from_smiles(smiles: str, charge: int, multiplicity: int):
    if charge != 0 or multiplicity != 1:
        raise NotImplementedError('CBH graph fragmentation currently requires '
                                  'a neutral closed-shell species.')
    try:
        import pybel
    except ImportError:
        try:
            from openbabel import pybel
        except ImportError as exc:
            raise ImportError('Open Babel Python bindings are required for CBH SMILES.') from exc
    ob = pybel.ob
    try:
        molecule = pybel.readstring('smi', smiles)
    except (OSError, ValueError) as exc:
        raise ValueError(f'Invalid CBH SMILES {smiles!r}.') from exc
    if molecule.OBMol.GetTotalCharge() != 0:
        raise NotImplementedError('Charged CBH species need state-aware fragments.')
    canonical = molecule.write('can').split()[0]
    molecule.OBMol.AddHydrogens()
    graph = nx.Graph()
    for atom in molecule.atoms:
        native = atom.OBAtom
        z = native.GetAtomicNum()
        if z == 1:
            continue
        if (z not in _VALENCE or native.GetFormalCharge() != 0
                or native.IsAromatic()):
            raise NotImplementedError('CBH currently supports neutral, '
                                      'nonaromatic C/N/O/F/Cl graphs.')
        hydrogens = sum(neighbor.GetAtomicNum() == 1
                        for neighbor in ob.OBAtomAtomIter(native))
        graph.add_node(atom.idx, z=z, hydrogens=hydrogens)
    if not graph:
        raise NotImplementedError('CBH needs at least one heavy atom.')
    for bond in ob.OBMolBondIter(molecule.OBMol):
        left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if left not in graph or right not in graph:
            continue
        order = bond.GetBondOrder()
        if bond.IsAromatic() or order not in (1, 2, 3):
            raise NotImplementedError('Aromatic or nonintegral CBH bonds need '
                                      'a state-aware fragmentation rule.')
        graph.add_edge(left, right, order=order)
    if not nx.is_connected(graph):
        raise ValueError('CBH target SMILES must contain one connected molecule.')
    for vertex, data in graph.nodes(data=True):
        valence = data['hydrogens'] + sum(graph[vertex][other]['order']
                                           for other in graph[vertex])
        if valence != _VALENCE[data['z']]:
            raise NotImplementedError('CBH needs explicit radical, ion, or '
                                      'unusual-valence rules for this structure.')
    return graph, canonical


def _formula(graph: nx.Graph, vertices) -> Counter:
    from ase.data import chemical_symbols

    selected = set(vertices)
    result = Counter()
    for vertex in selected:
        data = graph.nodes[vertex]
        result[chemical_symbols[data['z']]] += 1
        result['H'] += data['hydrogens']
        result['H'] += sum(graph[vertex][other]['order']
                           for other in graph[vertex] if other not in selected)
    return +result


def _fragment(graph: nx.Graph, vertices):
    try:
        import pybel
    except ImportError:
        from openbabel import pybel
    ob = pybel.ob

    selected = set(vertices)
    if not selected or not nx.is_connected(graph.subgraph(selected)):
        raise NotImplementedError('CBH overlapping fragment is disconnected.')
    native = ob.OBMol()
    indices = {}
    for vertex in sorted(selected):
        source = graph.nodes[vertex]
        atom = native.NewAtom()
        atom.SetAtomicNum(source['z'])
        cap_h = source['hydrogens'] + sum(
            graph[vertex][other]['order'] for other in graph[vertex]
            if other not in selected)
        atom.SetImplicitHCount(cap_h)
        indices[vertex] = atom.GetIdx()
    for left, right, data in graph.subgraph(selected).edges(data=True):
        native.AddBond(indices[left], indices[right], data['order'])
    smiles = pybel.Molecule(native).write('can').split()[0]
    if not smiles:
        raise ValueError('Open Babel could not canonicalize a CBH fragment.')
    return smiles, _formula(graph, selected)


def _cores(graph: nx.Graph, rung: int) -> list[frozenset[int]]:
    if rung == 0:
        return [frozenset((vertex,)) for vertex in graph]
    if rung == 1:
        return [frozenset((left, right)) for left, right in graph.edges]
    if rung == 2:
        return [frozenset((vertex, *graph.neighbors(vertex)))
                for vertex in graph if graph.degree(vertex) >= 2]
    if rung == 3:
        return [frozenset((left, right, *graph.neighbors(left),
                           *graph.neighbors(right)))
                for left, right in graph.edges
                if graph.degree(left) >= 2 and graph.degree(right) >= 2]
    raise ValueError('CBH rung must be 0, 1, 2, or 3.')


def _inclusion_exclusion(cores: list[frozenset[int]]):
    coefficients: dict[frozenset[int], int] = {}
    for core in cores:
        changes = {core: 1}
        for existing, coefficient in coefficients.items():
            overlap = core & existing
            if overlap:
                changes[overlap] = changes.get(overlap, 0) - coefficient
        for fragment, coefficient in changes.items():
            coefficients[fragment] = coefficients.get(fragment, 0) + coefficient
            if not coefficients[fragment]:
                del coefficients[fragment]
        if len(coefficients) > 10000:
            raise ValueError('CBH overlap expansion exceeded 10000 fragments.')
    return coefficients


def _check_balance(stoichiometry, formulas):
    balance = Counter()
    for smiles, coefficient in stoichiometry.items():
        for element, count in formulas[smiles].items():
            balance[element] += coefficient * count
    bad = {element: count for element, count in balance.items() if count}
    if bad:
        raise ValueError(f'CBH reaction is not element balanced: {bad}.')


def generate_cbh_reaction(smiles: str, rung: int, *, charge: int = 0,
                          multiplicity: int = 1) -> CBHReaction | None:
    """Generate a balanced CBH rung from a connected closed-shell SMILES.

    Higher rungs use inclusion-exclusion of overlapping capped graph
    neighborhoods. A rung that contains the complete target as a fragment is
    self-referential and returns ``None`` so callers can try a lower rung.
    """
    graph, target = _graph_from_smiles(smiles, charge, multiplicity)
    cores = _cores(graph, rung)
    if not cores or set().union(*cores) != set(graph):
        return None
    if rung > 0 and any(not any({left, right} <= core for core in cores)
                        for left, right in graph.edges):
        return None
    fragments = _inclusion_exclusion(cores) if rung > 0 else {core: 1 for core in cores}
    if any(not nx.is_connected(graph.subgraph(vertices)) for vertices in fragments):
        return None
    stoichiometry = Counter({target: -1})
    formulas = {target: _formula(graph, graph)}
    for vertices, coefficient in fragments.items():
        fragment, formula = _fragment(graph, vertices)
        if fragment in formulas and formulas[fragment] != formula:
            raise ValueError(f'CBH canonical fragment {fragment} has inconsistent formulas.')
        formulas[fragment] = formula
        stoichiometry[fragment] += coefficient
    stoichiometry = {key: value for key, value in stoichiometry.items() if value}
    if target not in stoichiometry or stoichiometry[target] != -1:
        return None
    if rung == 0:
        hydrogen_excess = sum(coefficient * formulas[key].get('H', 0)
                              for key, coefficient in stoichiometry.items())
        if hydrogen_excess % 2:
            raise ValueError('CBH-0 saturation has an odd hydrogen imbalance.')
        if hydrogen_excess:
            stoichiometry[_H2] = -hydrogen_excess // 2
            formulas[_H2] = Counter({'H': 2})
    _check_balance(stoichiometry, formulas)
    return CBHReaction(rung, target, dict(sorted(stoichiometry.items())),
                       {key: dict(value) for key, value in formulas.items()
                        if key in stoichiometry})


def generate_for_stationary_point(species, *, max_rung: int = 3
                                  ) -> dict[int, CBHReaction | None]:
    """Use KinBot's accepted species SMILES and verify its atom inventory."""
    if getattr(species, 'wellorts', 0):
        raise ValueError('A transition state is not a CBH formation reference.')
    smiles = getattr(species, 'smiles', '')
    if not smiles:
        raise ValueError('KinBot species has no SMILES connectivity for CBH.')
    charge = getattr(species, 'charge', None)
    multiplicity = getattr(species, 'mult', None)
    graph, _ = _graph_from_smiles(smiles, charge, multiplicity)
    atoms = list(getattr(species, 'atom', []))
    if atoms and Counter(atoms) != _formula(graph, graph):
        raise ValueError('KinBot species atoms disagree with its CBH SMILES.')
    if max_rung not in range(4):
        raise ValueError('max_rung must be 0, 1, 2, or 3.')
    return {rung: generate_cbh_reaction(smiles, rung, charge=charge,
                                        multiplicity=multiplicity)
            for rung in range(max_rung + 1)}


def solve_formation_enthalpy(reaction: CBHReaction,
                             energies: Mapping[str, ZeroKEnergy],
                             atct: ATcTTable, *,
                             reference_ids: Mapping[str, str] | None = None
                             ) -> FormationEnthalpy:
    """Apply the signed CBH reaction to ANL electronic-plus-ZPE 0 K energies."""
    if reaction.stoichiometry.get(reaction.target_smiles) != -1:
        raise ValueError('CBH target must have stoichiometric coefficient -1.')
    _check_balance(reaction.stoichiometry, reaction.formulas)
    if set(energies) < set(reaction.stoichiometry):
        missing = sorted(set(reaction.stoichiometry) - set(energies))
        raise ValueError('Missing 0 K energies for: ' + ', '.join(missing))
    methods = set()
    for smiles in reaction.stoichiometry:
        energy = energies[smiles]
        if (not isinstance(energy, ZeroKEnergy) or not math.isfinite(energy.hartree)
                or not energy.method or not energy.source
                or _canonical_smiles(energy.smiles) != _canonical_smiles(smiles)
                or (energy.charge, energy.multiplicity) !=
                   reaction.states.get(smiles, (0, 1))):
            raise ValueError(f'{smiles}: invalid accepted 0 K energy.')
        methods.add(energy.method)
    if len(methods) != 1:
        raise ValueError('All species in one CBH reaction require the same '
                         '0 K electronic-structure method.')
    references = {}
    reference_ids = reference_ids or {}
    for smiles in reaction.reference_smiles:
        if reaction.states.get(smiles, (0, 1)) != (0, 1) and smiles not in reference_ids:
            raise ValueError(f'{smiles}: a non-singlet ATcT reference needs an explicit ID.')
        record = (atct.by_id(reference_ids[smiles]) if smiles in reference_ids
                  else atct.gas_by_smiles(smiles, formula=reaction.formulas[smiles]))
        if (record.phase != 'g' or record.formation_0k_kj_mol is None
                or record.version != atct.version
                or record.source_sha256 != atct.source_sha256
                or not record.atct_id.endswith('*0')
                or not record.preferred_formula.endswith('(g)')
                or _canonical_smiles(record.smiles) != _canonical_smiles(smiles)
                or _preferred_formula_atoms(record.preferred_formula)
                    != reaction.formulas[smiles]):
            raise ValueError(f'{smiles}: ATcT reference is not a pinned 0 K gas value.')
        references[smiles] = record
    delta = math.fsum(coefficient * energies[smiles].hartree
                      for smiles, coefficient in reaction.stoichiometry.items())
    delta *= HARTREE_TO_KJ_MOL
    reference_sum = math.fsum(reaction.stoichiometry[smiles] *
                              record.formation_0k_kj_mol
                              for smiles, record in references.items())
    # delta_r H(0) = sum(nu_i Hf_i(0)); nu_target is -1.
    target_hf = reference_sum - delta
    return FormationEnthalpy(
        reaction.target_smiles, reaction.rung, methods.pop(), delta,
        target_hf, references, atct.version, atct.source_sha256,
        {smiles: energies[smiles].source for smiles in reaction.stoichiometry})


def _method_family(label: str) -> str:
    if label.startswith('profiled:'):
        return label.split(':', 2)[1]
    if label.startswith('L3:'):
        return 'L3'
    return label


def select_cbh_ladder(
    species, energies_by_smiles: Mapping[str, Mapping[str, ZeroKEnergy]],
    atct: ATcTTable, *, max_rung: int = 3,
    method_ladder: tuple[str, ...] = DEFAULT_METHOD_LADDER,
    reference_ids: Mapping[str, str] | None = None,
) -> LadderSelection:
    """Choose the highest available rung, then highest common energy tier.

    The table's reference values are independent of the computed tier. Every
    molecule in the chosen CBH reaction must have an accepted 0 K energy at
    the *same exact method label*. A missing/failed tier falls to the next
    configured tier; inconsistent accepted data raises rather than falling.
    """
    if len(method_ladder) != len(set(method_ladder)) or not method_ladder:
        raise ValueError('CBH method ladder needs unique, ordered tiers.')
    reactions = generate_for_stationary_point(species, max_rung=max_rung)
    skipped = []
    for rung in range(max_rung, -1, -1):
        reaction = reactions[rung]
        if reaction is None:
            skipped.append(f'CBH-{rung}: self-referential or unavailable graph fragments')
            continue
        references = {}
        try:
            for smiles in reaction.reference_smiles:
                references[smiles] = (
                    atct.by_id(reference_ids[smiles])
                    if reference_ids is not None and smiles in reference_ids
                    else atct.gas_by_smiles(smiles, formula=reaction.formulas[smiles]))
        except LookupError as exc:
            if reference_ids is not None and smiles in reference_ids:
                raise
            skipped.append(f'CBH-{rung}: {exc}')
            continue
        needed = set(reaction.stoichiometry)
        common = set.intersection(*(set(energies_by_smiles.get(smiles, {}))
                                    for smiles in needed))
        for tier in method_ladder:
            candidates = sorted(label for label in common
                                if label == tier or _method_family(label) == tier)
            if len(candidates) > 1:
                raise ValueError(f'CBH-{rung}: multiple accepted {tier} profiles; '
                                 'select one exact method label.')
            if not candidates:
                skipped.append(f'CBH-{rung} {tier}: missing accepted 0 K energy')
                continue
            label = candidates[0]
            selected = {smiles: energies_by_smiles[smiles][label]
                        for smiles in needed}
            formation = solve_formation_enthalpy(
                reaction, selected, atct, reference_ids=reference_ids)
            return LadderSelection(formation, reaction, tuple(skipped))
    raise ValueError('No complete CBH reaction and energy tier is available. '
                     + '; '.join(skipped))
