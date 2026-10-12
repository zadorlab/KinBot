"""Connectivity-based hierarchy reactions and 0 K formation enthalpies.

This implementation handles neutral, closed-shell molecules with ordinary
integral bond orders and the neutral-doublet CBH-0 hydrogen-cap rule. It
constructs capped graph fragments from the same explicit-hydrogen SMILES
connectivity that KinBot uses for its species. Every generated reaction is
checked for elemental balance. Higher radical rungs, ions, and aromatic-state
extensions require additional state-aware fragment rules and fail closed.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass, field, fields
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Mapping

import networkx as nx
from ase.units import Hartree, kJ, mol

from kinbot.anl.atct import (ATcTRecord, ATcTTable, _canonical_smiles,
                            _preferred_formula_atoms, load_atct_references)
from kinbot.energy import ZeroKEnergy


HARTREE_TO_KJ_MOL = Hartree * mol / kJ
_VALENCE = {6: 4, 7: 3, 8: 2, 9: 1, 17: 1}
_H2 = '[H][H]'
_H = '[H]'
DEFAULT_METHOD_LADDER = ('ANL1-F12', 'ANL1', 'ANL0-F12', 'ANL0',
                         'L3', 'L2')
_SHA256 = re.compile(r'[0-9a-f]{64}\Z')


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
    # Immediate reference values.  CBH-0 terminates in pinned ATcT records;
    # higher rungs contain formation enthalpies solved from lower CBH rungs.
    references: Mapping[str, ATcTRecord | 'FormationEnthalpy']
    atct_version: str
    atct_source_sha256: str
    energy_sources: Mapping[str, str]
    # Preserve the actual reaction used for both the 0 K and finite-
    # temperature solves.  Older CBH-0 records are upgraded from their
    # top-level ``reaction`` object by the handoff reader.
    stoichiometry: Mapping[str, int] = field(default_factory=dict)
    formulas: Mapping[str, Mapping[str, int]] = field(default_factory=dict)
    states: Mapping[str, tuple[int, int]] = field(default_factory=dict)


@dataclass(frozen=True)
class LadderSelection:
    formation: FormationEnthalpy
    reaction: CBHReaction
    skipped: tuple[str, ...]


def formation_from_dict(value, reaction=None) -> FormationEnthalpy:
    """Reconstruct and validate a recursive CBH formation record."""
    if not isinstance(value, dict):
        raise ValueError('CBH formation record must be an object.')
    required = {
        'target_smiles', 'rung', 'method', 'reaction_energy_0k_kj_mol',
        'formation_0k_kj_mol', 'references', 'atct_version',
        'atct_source_sha256', 'energy_sources'}
    if not required <= set(value) or not isinstance(value['references'], dict):
        raise ValueError('CBH formation record is incomplete.')
    atct_fields = {item.name for item in fields(ATcTRecord)}
    references = {}
    for smiles, reference in value['references'].items():
        if not isinstance(reference, dict):
            raise ValueError(f'{smiles}: invalid CBH reference record.')
        if set(reference) == atct_fields:
            references[smiles] = ATcTRecord(**reference)
        else:
            references[smiles] = formation_from_dict(reference)
    stoichiometry = value.get('stoichiometry')
    formulas = value.get('formulas')
    states = value.get('states')
    if reaction is not None:
        stoichiometry = (reaction.stoichiometry if stoichiometry is None
                         else stoichiometry)
        formulas = reaction.formulas if formulas is None else formulas
        states = reaction.states if states is None else states
    stoichiometry = dict(stoichiometry or {})
    formulas = {smiles: dict(formula)
                for smiles, formula in (formulas or {}).items()}
    states = {smiles: tuple(state)
              for smiles, state in (states or {}).items()}
    result = FormationEnthalpy(
        target_smiles=value['target_smiles'], rung=value['rung'],
        method=value['method'],
        reaction_energy_0k_kj_mol=value['reaction_energy_0k_kj_mol'],
        formation_0k_kj_mol=value['formation_0k_kj_mol'],
        references=references, atct_version=value['atct_version'],
        atct_source_sha256=value['atct_source_sha256'],
        energy_sources=dict(value['energy_sources']),
        stoichiometry=stoichiometry, formulas=formulas, states=states)
    if (isinstance(result.rung, bool) or result.rung not in range(4)
            or not isinstance(result.method, str) or not result.method
            or not isinstance(result.atct_version, str)
            or not result.atct_version
            or not _SHA256.fullmatch(result.atct_source_sha256)
            or any(not isinstance(source, str) or not source
                   for source in result.energy_sources.values())
            or not all(math.isfinite(number) for number in (
                result.reaction_energy_0k_kj_mol,
                result.formation_0k_kj_mol))):
        raise ValueError('CBH formation identity or values are invalid.')
    state_species = (set(result.stoichiometry) if result.stoichiometry
                     else {result.target_smiles})
    for smiles, state in result.states.items():
        if (smiles not in state_species
                or len(state) != 2
                or any(isinstance(number, bool) or not isinstance(number, int)
                       for number in state)
                or state[1] < 1):
            raise ValueError('CBH electronic-state metadata are invalid.')
    if not result.stoichiometry:
        # A terminal wrapper is allowed only around the matching pinned ATcT
        # species.  It exists so higher rungs retain an explicit CBH-0 leaf.
        reference = result.references.get(result.target_smiles)
        if (len(result.references) != 1
                or not isinstance(reference, ATcTRecord)
                or result.rung != 0
                or reference.version != result.atct_version
                or reference.source_sha256 != result.atct_source_sha256
                or reference.phase != 'g'
                or _canonical_smiles(reference.smiles) !=
                   _canonical_smiles(result.target_smiles)
                or reference.formation_0k_kj_mol is None
                or reference.formation_298k_kj_mol is None
                or not math.isclose(
                    result.formation_0k_kj_mol,
                    reference.formation_0k_kj_mol,
                    rel_tol=0., abs_tol=1e-10)):
            raise ValueError('Invalid terminal CBH-0/ATcT reference.')
        return result
    if result.stoichiometry.get(result.target_smiles) != -1:
        raise ValueError('CBH formation target coefficient must be -1.')
    if (set(result.formulas) != set(result.stoichiometry)
            or set(result.energy_sources) != set(result.stoichiometry)):
        raise ValueError('CBH formation provenance differs from its reaction.')
    _check_balance(result.stoichiometry, result.formulas)
    expected_references = set(result.stoichiometry) - {result.target_smiles}
    if set(result.references) != expected_references:
        raise ValueError('CBH formation references disagree with its reaction.')
    for smiles, reference in result.references.items():
        if _canonical_smiles(reference.smiles if isinstance(reference, ATcTRecord)
                             else reference.target_smiles) != \
                _canonical_smiles(smiles):
            raise ValueError(f'{smiles}: CBH reference identity differs.')
        if isinstance(reference, ATcTRecord):
            if (result.rung != 0
                    or reference.version != result.atct_version
                    or reference.source_sha256 != result.atct_source_sha256
                    or reference.phase != 'g'
                    or reference.formation_0k_kj_mol is None
                    or reference.formation_298k_kj_mol is None):
                raise ValueError(f'{smiles}: invalid direct ATcT reference.')
        elif (result.rung == 0 or reference.rung >= result.rung
              or reference.atct_version != result.atct_version
              or reference.atct_source_sha256 != result.atct_source_sha256):
            raise ValueError(f'{smiles}: invalid lower-rung CBH reference.')
    reference_sum = math.fsum(
        result.stoichiometry[smiles] * reference.formation_0k_kj_mol
        for smiles, reference in result.references.items())
    expected = reference_sum - result.reaction_energy_0k_kj_mol
    if not math.isclose(result.formation_0k_kj_mol, expected,
                        rel_tol=0., abs_tol=1e-8):
        raise ValueError('CBH formation 0 K identity is inconsistent.')
    return result


def formation_species(formation: FormationEnthalpy) -> tuple[str, ...]:
    """List every calculated species required by a stored CBH tree."""
    if not isinstance(formation, FormationEnthalpy):
        raise TypeError('Expected a CBH formation record.')
    result = set(formation.stoichiometry)
    for reference in formation.references.values():
        if isinstance(reference, FormationEnthalpy):
            result.update(formation_species(reference))
    return tuple(sorted(result))


def _graph_from_smiles(smiles: str, charge: int, multiplicity: int):
    if charge != 0 or multiplicity not in (1, 2):
        raise NotImplementedError('CBH graph fragmentation currently requires '
                                  'a neutral singlet or doublet species.')
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
    deficits = []
    for vertex, data in graph.nodes(data=True):
        valence = data['hydrogens'] + sum(graph[vertex][other]['order']
                                           for other in graph[vertex])
        deficit = _VALENCE[data['z']] - valence
        data['radical_caps'] = deficit
        if deficit:
            deficits.append((vertex, deficit))
        if deficit < 0:
            raise NotImplementedError('CBH needs explicit radical, ion, or '
                                      'unusual-valence rules for this structure.')
    if multiplicity == 1 and deficits:
        raise NotImplementedError('CBH needs explicit radical, ion, or '
                                  'unusual-valence rules for this structure.')
    if multiplicity == 2 and (len(deficits) != 1 or deficits[0][1] != 1):
        raise NotImplementedError('Neutral-doublet CBH-0 requires exactly one '
                                  'one-electron valence deficit.')
    return graph, canonical


def _formula(graph: nx.Graph, vertices, *, cap_radicals=False) -> Counter:
    from ase.data import chemical_symbols

    selected = set(vertices)
    result = Counter()
    for vertex in selected:
        data = graph.nodes[vertex]
        result[chemical_symbols[data['z']]] += 1
        result['H'] += data['hydrogens']
        if cap_radicals:
            result['H'] += data['radical_caps']
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
        cap_h += source['radical_caps']
        atom.SetImplicitHCount(cap_h)
        indices[vertex] = atom.GetIdx()
    for left, right, data in graph.subgraph(selected).edges(data=True):
        native.AddBond(indices[left], indices[right], data['order'])
    smiles = pybel.Molecule(native).write('can').split()[0]
    if not smiles:
        raise ValueError('Open Babel could not canonicalize a CBH fragment.')
    return smiles, _formula(graph, selected, cap_radicals=True)


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
    """Generate a balanced CBH rung from connected molecular connectivity.

    Higher rungs use inclusion-exclusion of overlapping capped graph
    neighborhoods. A rung that contains the complete target as a fragment is
    self-referential and returns ``None`` so callers can try a lower rung.
    A neutral doublet is supported only at CBH-0, where its single radical
    valence is capped with H and atomic H carries the doublet state.
    """
    graph, target = _graph_from_smiles(smiles, charge, multiplicity)
    if multiplicity == 2 and rung > 0:
        return None
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
        states = {target: (0, multiplicity)} if multiplicity != 1 else {}
        if multiplicity == 2:
            stoichiometry[_H] = 1
            formulas[_H] = Counter({'H': 1})
            states[_H] = (0, 2)
            hydrogen_excess += 1
        if hydrogen_excess % 2:
            raise ValueError('CBH-0 saturation has an odd hydrogen imbalance.')
        if hydrogen_excess:
            stoichiometry[_H2] = -hydrogen_excess // 2
            formulas[_H2] = Counter({'H': 2})
    else:
        states = {}
    _check_balance(stoichiometry, formulas)
    return CBHReaction(rung, target, dict(sorted(stoichiometry.items())),
                       {key: dict(value) for key, value in formulas.items()
                        if key in stoichiometry}, states)


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


def solve_formation_enthalpy(
        reaction: CBHReaction,
        energies: Mapping[str, ZeroKEnergy],
        atct: ATcTTable | None = None, *,
        reference_ids: Mapping[str, str] | None = None,
        reference_formations: Mapping[str, FormationEnthalpy] | None = None,
        ) -> FormationEnthalpy:
    """Apply one signed CBH reaction at 0 K.

    CBH-0 is the only rung allowed to terminate in ATcT values.  A higher
    rung must be supplied formation enthalpies that were solved from lower
    CBH rungs.  This prevents a nominal CBH-1/2/3 result from silently using
    direct ATcT values for its immediate fragments.
    """
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
    references: dict[str, ATcTRecord | FormationEnthalpy] = {}
    reference_ids = reference_ids or {}
    if reaction.rung == 0:
        if atct is None:
            raise ValueError('CBH-0 requires a pinned ATcT table.')
        if reference_formations:
            raise ValueError('CBH-0 must terminate in ATcT references.')
        for smiles in reaction.reference_smiles:
            references[smiles] = _atct_reference(
                smiles, reaction.formulas[smiles],
                reaction.states.get(smiles, (0, 1)), atct, reference_ids)
        atct_version = atct.version
        atct_source_sha256 = atct.source_sha256
    else:
        supplied = dict(reference_formations or {})
        if set(supplied) != set(reaction.reference_smiles):
            missing = sorted(set(reaction.reference_smiles) - set(supplied))
            extra = sorted(set(supplied) - set(reaction.reference_smiles))
            raise ValueError(
                f'CBH-{reaction.rung} lower-rung references differ; '
                f'missing={missing}, extra={extra}.')
        for smiles in reaction.reference_smiles:
            value = supplied[smiles]
            if (not isinstance(value, FormationEnthalpy)
                    or _canonical_smiles(value.target_smiles) !=
                       _canonical_smiles(smiles)
                    or value.rung >= reaction.rung
                    or not math.isfinite(value.formation_0k_kj_mol)):
                raise ValueError(
                    f'{smiles}: invalid lower-rung CBH formation reference.')
            references[smiles] = value
        versions = {value.atct_version for value in supplied.values()}
        hashes = {value.atct_source_sha256 for value in supplied.values()}
        if len(versions) != 1 or len(hashes) != 1:
            raise ValueError('Lower-rung CBH references use different ATcT releases.')
        atct_version = versions.pop()
        atct_source_sha256 = hashes.pop()
    delta = math.fsum(coefficient * energies[smiles].hartree
                      for smiles, coefficient in reaction.stoichiometry.items())
    delta *= HARTREE_TO_KJ_MOL
    reference_sum = math.fsum(
        reaction.stoichiometry[smiles] * record.formation_0k_kj_mol
        for smiles, record in references.items())
    # delta_r H(0) = sum(nu_i Hf_i(0)); nu_target is -1.
    target_hf = reference_sum - delta
    return FormationEnthalpy(
        reaction.target_smiles, reaction.rung, methods.pop(), delta,
        target_hf, references, atct_version, atct_source_sha256,
        {smiles: energies[smiles].source for smiles in reaction.stoichiometry},
        dict(reaction.stoichiometry),
        {smiles: dict(formula) for smiles, formula in reaction.formulas.items()},
        dict(reaction.states))


def _lower_reaction(smiles: str, max_rung: int, *, charge: int,
                    multiplicity: int):
    """Return the highest structurally available lower-rung reaction."""
    for rung in range(max_rung, -1, -1):
        try:
            reaction = generate_cbh_reaction(
                smiles, rung, charge=charge, multiplicity=multiplicity)
        except NotImplementedError:
            return None
        if reaction is not None:
            return reaction
    return None


def _atct_reference(smiles, formula, state, atct, reference_ids):
    """Return one pinned gas-phase CBH-0 reference at 0 and 298 K."""
    if state != (0, 1) and smiles not in reference_ids:
        raise ValueError(
            f'{smiles}: a non-singlet ATcT reference needs an explicit ID.')
    record = (atct.by_id(reference_ids[smiles]) if smiles in reference_ids
              else atct.gas_by_smiles(smiles, formula=formula))
    if (record.phase != 'g' or record.formation_0k_kj_mol is None
            or record.formation_298k_kj_mol is None
            or record.version != atct.version
            or record.source_sha256 != atct.source_sha256
            or not record.atct_id.endswith('*0')
            or not record.preferred_formula.endswith('(g)')
            or _canonical_smiles(record.smiles) != _canonical_smiles(smiles)
            or _preferred_formula_atoms(record.preferred_formula) != formula):
        raise ValueError(
            f'{smiles}: ATcT reference is not a pinned 0/298 K gas value.')
    return record


def solve_cbh_network(reaction: CBHReaction,
                      energies: Mapping[str, ZeroKEnergy | Mapping[str, ZeroKEnergy]],
                      atct: ATcTTable, *,
                      reference_ids: Mapping[str, str] | None = None,
                      method_ladder: tuple[str, ...] = DEFAULT_METHOD_LADDER,
                      root_method: str | None = None,
                      ) -> FormationEnthalpy:
    """Solve a CBH reaction through lower rungs to verified CBH-0 leaves.

    Each reaction node selects one common electronic-structure method for all
    species in that reaction.  A lower-rung reference node selects its method
    independently, which is the defining ladder behavior: for example, a
    large target reaction may use F12 energies while its smaller reference
    species carry ANL0 formation enthalpies.
    """
    reference_ids = dict(reference_ids or {})
    if len(method_ladder) != len(set(method_ladder)) or not method_ladder:
        raise ValueError('CBH method ladder needs unique, ordered tiers.')
    profiles = {}
    for smiles, value in energies.items():
        if isinstance(value, ZeroKEnergy):
            profiles[smiles] = {value.method: value}
        elif isinstance(value, Mapping):
            profiles[smiles] = dict(value)
            if any(not isinstance(item, ZeroKEnergy)
                   or label != item.method
                   for label, item in profiles[smiles].items()):
                raise ValueError(f'{smiles}: invalid CBH energy profile.')
        else:
            raise ValueError(f'{smiles}: invalid CBH energy profile.')
    cache: dict[tuple[str, int, str], FormationEnthalpy] = {}

    def candidates(node, forced=None):
        missing = sorted(set(node.stoichiometry) - set(profiles))
        if missing:
            raise ValueError('Missing 0 K energies for: ' + ', '.join(missing))
        common = set.intersection(*(
            set(profiles[smiles]) for smiles in node.stoichiometry))
        if forced is not None:
            return [forced] if forced in common else []
        ordered = []
        for tier in method_ladder:
            matched = sorted(label for label in common
                             if label == tier or _method_family(label) == tier)
            if len(matched) > 1:
                raise ValueError(
                    f'CBH-{node.rung}: multiple accepted {tier} profiles; '
                    'select one exact method label.')
            ordered.extend(matched)
        ordered.extend(sorted(common - set(ordered)))
        return ordered

    def solve_node(node: CBHReaction, forced=None) -> FormationEnthalpy:
        labels = candidates(node, forced)
        if not labels:
            detail = f' at {forced}' if forced is not None else ''
            raise ValueError(
                f'CBH-{node.rung} has no common 0 K energy method{detail}.')
        failures = []
        for label in labels:
            try:
                return solve_label(node, label)
            except (LookupError, ValueError, NotImplementedError) as error:
                failures.append(f'{label}: {error}')
        raise ValueError(
            f'CBH-{node.rung} has no complete energy/reference tier: '
            + '; '.join(failures))

    def solve_label(node: CBHReaction, label: str) -> FormationEnthalpy:
        key = (node.target_smiles, node.rung, label)
        if key in cache:
            return cache[key]
        selected = {smiles: profiles[smiles][label]
                    for smiles in node.stoichiometry}
        if node.rung == 0:
            result = solve_formation_enthalpy(
                node, selected, atct, reference_ids=reference_ids)
        else:
            lower = {}
            for smiles in node.reference_smiles:
                state = node.states.get(smiles, (0, 1))
                child = _lower_reaction(
                    smiles, node.rung - 1, charge=state[0],
                    multiplicity=state[1])
                if child is None:
                    # Wrap a terminal CBH-0 datum in a zero-reaction node so
                    # every higher-rung input demonstrably descends through
                    # the CBH-0/ATcT boundary.
                    record = _atct_reference(
                        smiles, node.formulas[smiles], state, atct,
                        reference_ids)
                    lower[smiles] = FormationEnthalpy(
                        target_smiles=smiles, rung=0, method=label,
                        reaction_energy_0k_kj_mol=0.,
                        formation_0k_kj_mol=record.formation_0k_kj_mol,
                        references={smiles: record},
                        atct_version=atct.version,
                        atct_source_sha256=atct.source_sha256,
                        energy_sources={},
                        stoichiometry={},
                        formulas={smiles: dict(node.formulas[smiles])},
                        states={smiles: state})
                else:
                    lower[smiles] = solve_node(child)
            result = solve_formation_enthalpy(
                node, selected, reference_formations=lower)
        cache[key] = result
        return result

    return solve_node(reaction, root_method)


def cbh_network_requirements(reaction: CBHReaction):
    """Return all calculated species and terminal ATcT leaves in a CBH tree."""
    calculated = set()
    terminals = {}
    visiting = set()

    def visit(node):
        key = (node.target_smiles, node.rung)
        if key in visiting:
            raise ValueError('CBH reference network contains a cycle.')
        visiting.add(key)
        calculated.update(node.stoichiometry)
        for smiles in node.reference_smiles:
            state = node.states.get(smiles, (0, 1))
            if node.rung == 0:
                terminals[smiles] = (dict(node.formulas[smiles]), state)
                continue
            child = _lower_reaction(
                smiles, node.rung - 1, charge=state[0],
                multiplicity=state[1])
            if child is None:
                terminals[smiles] = (dict(node.formulas[smiles]), state)
            else:
                visit(child)
        visiting.remove(key)

    visit(reaction)
    return calculated, terminals


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
        for tier in method_ladder:
            candidates = sorted(
                label for label in energies_by_smiles.get(
                    reaction.target_smiles, {})
                                if label == tier or _method_family(label) == tier)
            if len(candidates) > 1:
                raise ValueError(f'CBH-{rung}: multiple accepted {tier} profiles; '
                                 'select one exact method label.')
            if not candidates:
                skipped.append(f'CBH-{rung} {tier}: missing accepted 0 K energy')
                continue
            label = candidates[0]
            try:
                formation = solve_cbh_network(
                    reaction, energies_by_smiles, atct,
                    reference_ids=reference_ids, method_ladder=method_ladder,
                    root_method=label)
            except (LookupError, ValueError) as exc:
                skipped.append(f'CBH-{rung} {tier}: {exc}')
                continue
            return LadderSelection(formation, reaction, tuple(skipped))
    raise ValueError('No complete CBH reaction and energy tier is available. '
                     + '; '.join(skipped))


def zero_k_from_composite_record(smiles: str, path: str | Path) -> ZeroKEnergy:
    """Load one complete, internally consistent composite JSON result."""
    path = Path(path).resolve()
    try:
        raw = path.read_bytes()
        record = json.loads(raw)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f'{path}: invalid composite result.') from exc
    if (not isinstance(record, dict) or record.get('schema') != 1
            or record.get('status') != 'complete'
            or not isinstance(record.get('recipe'), str)
            or not record['recipe'].strip()
            or isinstance(record.get('charge'), bool)
            or not isinstance(record.get('charge'), int)
            or isinstance(record.get('multiplicity'), bool)
            or not isinstance(record.get('multiplicity'), int)
            or record['multiplicity'] < 1):
        raise ValueError(f'{path}: composite result identity is incomplete.')
    components = record.get('components')
    if not isinstance(components, dict) or not components:
        raise ValueError(f'{path}: composite components are missing.')
    for key, component in components.items():
        if (not isinstance(key, str) or not isinstance(component, dict)
                or component.get('key') != key
                or component.get('review_required') is not False
                or not isinstance(component.get('value_hartree'), (int, float))
                or isinstance(component.get('value_hartree'), bool)
                or not math.isfinite(component['value_hartree'])
                or not _SHA256.fullmatch(str(component.get('source_sha256', '')))):
            raise ValueError(f'{path}: invalid accepted component {key!r}.')

    def total(terms, label):
        if not isinstance(terms, list) or not terms:
            raise ValueError(f'{path}: {label} expression is missing.')
        values = []
        seen = set()
        for term in terms:
            if (not isinstance(term, dict)
                    or set(term) != {'component', 'coefficient'}
                    or term['component'] not in components
                    or term['component'] in seen
                    or term['coefficient'] not in (-1, 1)):
                raise ValueError(f'{path}: invalid {label} expression.')
            seen.add(term['component'])
            values.append(term['coefficient'] *
                          components[term['component']]['value_hartree'])
        return math.fsum(values), seen

    electronic, electronic_keys = total(
        record.get('electronic_terms'), 'electronic')
    zero_point, zero_point_keys = total(
        record.get('zero_point_terms'), 'zero-point')
    if (electronic_keys & zero_point_keys
            or electronic_keys | zero_point_keys != set(components)):
        raise ValueError(f'{path}: expression does not use each component once.')
    expected = {
        'electronic_hartree': electronic,
        'zero_point_hartree': zero_point,
        'zero_k_hartree': math.fsum((electronic, zero_point)),
    }
    for key, value in expected.items():
        stored = record.get(key)
        if (not isinstance(stored, (int, float)) or isinstance(stored, bool)
                or not math.isfinite(stored) or abs(stored - value) > 5e-12):
            raise ValueError(f'{path}: stored {key} disagrees with components.')
    source = f'{path}#sha256={sha256(raw).hexdigest()}'
    return ZeroKEnergy(
        smiles, expected['zero_k_hartree'], record['recipe'], source,
        record['charge'], record['multiplicity'])


def solve_composite_records(
        target_smiles: str, rung: int,
        energy_files: Mapping[
            str, str | Path | Mapping[str, str | Path]], *, atct_version: str,
        atct_cache: str | Path,
        reference_ids: Mapping[str, str] | None = None,
        refresh_atct: bool = False) -> tuple[CBHReaction, FormationEnthalpy]:
    """Generate and solve a CBH network from composite JSON records.

    A species value may be one composite path or a ``method: path`` mapping.
    The latter lets each rung select its own highest complete method while
    keeping one common method within each individual reaction.
    """
    if target_smiles not in energy_files:
        raise ValueError('The target composite result is missing.')

    def load_profile(smiles, value):
        if isinstance(value, Mapping):
            if not value:
                raise ValueError(f'{smiles}: composite energy profile is empty.')
            profile = {}
            for declared, path in value.items():
                energy = zero_k_from_composite_record(smiles, path)
                if declared != energy.method or declared in profile:
                    raise ValueError(
                        f'{smiles}: composite profile label differs from its recipe.')
                profile[declared] = energy
            return profile
        energy = zero_k_from_composite_record(smiles, value)
        return {energy.method: energy}

    target_profile = load_profile(target_smiles, energy_files[target_smiles])
    target = next(iter(target_profile.values()))
    if any((item.charge, item.multiplicity) !=
           (target.charge, target.multiplicity)
           for item in target_profile.values()):
        raise ValueError('Target composite profiles use different states.')
    reaction = generate_cbh_reaction(
        target_smiles, rung, charge=target.charge,
        multiplicity=target.multiplicity)
    if reaction is None:
        raise ValueError(f'CBH-{rung} is unavailable for {target_smiles!r}.')
    required, terminals = cbh_network_requirements(reaction)
    missing = sorted(required - set(energy_files))
    extra = sorted(set(energy_files) - required)
    if missing or extra:
        raise ValueError(f'Composite energy mapping differs from the CBH-{rung} '
                         f'reference network; '
                         f'missing={missing}, extra={extra}.')
    energies = {smiles: (target_profile if smiles == target_smiles else
                         load_profile(smiles, energy_files[smiles]))
                for smiles in required}
    reference_ids = dict(reference_ids or {})
    table = load_atct_references(
        terminals, atct_version, atct_cache,
        reference_ids=reference_ids, refresh=refresh_atct)
    return reaction, solve_cbh_network(
        reaction, energies, table, reference_ids=reference_ids)


def _assignments(values, label):
    result = {}
    for value in values:
        if '=' not in value:
            raise ValueError(f'{label} must use SMILES=value syntax.')
        smiles, assigned = value.split('=', 1)
        if not smiles or not assigned or smiles in result:
            raise ValueError(f'Invalid or repeated {label} for {smiles!r}.')
        result[smiles] = assigned
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Generate and solve provenance-checked CBH reactions')
    commands = parser.add_subparsers(dest='action', required=True)
    generate = commands.add_parser('generate')
    generate.add_argument('target_smiles')
    generate.add_argument('--rung', type=int, required=True, choices=range(4))
    generate.add_argument('--charge', type=int, default=0)
    generate.add_argument('--multiplicity', type=int, default=1)
    solve = commands.add_parser('solve-records')
    solve.add_argument('target_smiles')
    solve.add_argument('output', type=Path)
    solve.add_argument('--rung', type=int, required=True, choices=range(4))
    solve.add_argument('--energy', action='append', required=True,
                       help='SMILES=complete_composite_result.json')
    solve.add_argument('--reference-id', action='append', default=[],
                       help='SMILES=ATcT_ID for an ambiguous state')
    solve.add_argument('--atct-version', required=True)
    solve.add_argument('--atct-cache', type=Path, required=True)
    solve.add_argument('--refresh-atct', action='store_true')
    args = parser.parse_args(argv)
    if args.action == 'generate':
        reaction = generate_cbh_reaction(
            args.target_smiles, args.rung, charge=args.charge,
            multiplicity=args.multiplicity)
        if reaction is None:
            raise ValueError(
                f'CBH-{args.rung} is unavailable for {args.target_smiles!r}.')
        print(json.dumps({
            'rung': reaction.rung,
            'target_smiles': reaction.target_smiles,
            'stoichiometry': dict(reaction.stoichiometry),
            'formulas': dict(reaction.formulas),
            'states': {key: list(value)
                       for key, value in reaction.states.items()},
        }, indent=2, sort_keys=True))
        return 0
    energy_files = _assignments(args.energy, 'energy')
    reference_ids = _assignments(args.reference_id, 'reference ID')
    reaction, formation = solve_composite_records(
        args.target_smiles, args.rung, energy_files,
        atct_version=args.atct_version, atct_cache=args.atct_cache,
        reference_ids=reference_ids, refresh_atct=args.refresh_atct)
    payload = {
        'schema': 1, 'status': 'complete',
        'reaction': {
            'rung': reaction.rung,
            'target_smiles': reaction.target_smiles,
            'stoichiometry': dict(reaction.stoichiometry),
            'formulas': dict(reaction.formulas),
            'states': {key: list(value)
                       for key, value in reaction.states.items()},
        },
        'formation': {
            **asdict(formation),
            'formation_0k_kcal_mol': formation.formation_0k_kj_mol / 4.184,
            'reaction_energy_0k_kcal_mol':
                formation.reaction_energy_0k_kj_mol / 4.184,
        },
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    print(args.output.resolve())
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
