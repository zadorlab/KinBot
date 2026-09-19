"""One 0 K energy convention for legacy and composite KinBot results."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math

from ase.units import Hartree, kJ, mol


KJ_PER_KCAL = 4.184
HARTREE_TO_KJ_MOL = Hartree * mol / kJ


@dataclass(frozen=True)
class ZeroKEnergy:
    smiles: str
    hartree: float
    method: str
    source: str
    charge: int = 0
    multiplicity: int = 1

    @classmethod
    def from_composite(cls, smiles: str, result):
        """Keep the accepted composite expression and native sources linked."""
        from kinbot.anl.model import CompositeResult

        if not isinstance(result, CompositeResult):
            raise TypeError('Expected an evaluated ANL CompositeResult.')
        payload = {'recipe': result.recipe, 'version': result.recipe_version,
                   'zero_k_hartree': result.zero_k_hartree,
                   'components': {key: value.source_sha256
                                  for key, value in sorted(result.components.items())}}
        identity = {(value.charge, value.multiplicity)
                    for value in result.components.values()}
        if len(identity) != 1:
            raise ValueError('Composite components disagree on electronic state.')
        charge, multiplicity = identity.pop()
        return cls(smiles, result.zero_k_hartree, result.recipe,
                   sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest(),
                   charge, multiplicity)


@dataclass(frozen=True)
class ZeroKBarrier:
    """A provenance-checked 0 K barrier relative to one reactant well."""

    reactant_smiles: str
    transition_state_smiles: str
    barrier_0k_kj_mol: float
    method: str
    reactant_source: str
    transition_state_source: str


def species_zero_k_hartree(species) -> float:
    """Use accepted ANL E0 once, otherwise preserve electronic + L2 ZPE."""
    final = getattr(species, 'final_zero_k_energy', None)
    if final is not None:
        if (not isinstance(final, ZeroKEnergy) or not math.isfinite(final.hartree)
                or (final.charge, final.multiplicity) !=
                   (species.charge, species.mult)
                or final.smiles != species.smiles
                or not final.source):
            raise ValueError(f'{species.name}: invalid final 0 K energy.')
        return final.hartree
    electronic = float(species.energy)
    zpe = float(species.zpe)
    if not math.isfinite(electronic) or not math.isfinite(zpe):
        raise ValueError(f'{species.name}: nonfinite legacy energy or ZPE.')
    return math.fsum((electronic, zpe))


def attach_composite_energy(species, result) -> ZeroKEnergy:
    """Publish an evaluated ANL result on a KinBot well or transition state."""
    if not getattr(species, 'smiles', ''):
        raise ValueError('KinBot species needs a SMILES identity for ANL handoff.')
    final = ZeroKEnergy.from_composite(species.smiles, result)
    if (final.charge, final.multiplicity) != (species.charge, species.mult):
        raise ValueError('ANL result and KinBot species have different states.')
    species.final_zero_k_energy = final
    return final


def attach_anharmonic_frequencies(species, parsed_vpt2, *, mode_map=None,
                                  match_tolerance_cm_inverse=100.):
    """Apply mode-resolved VPT2 shifts to KinBot's harmonic MESS modes.

    ``species.reduced_freqs`` is the Hessian after internal rotations have
    been projected out, so it may contain fewer modes than Gaussian's VPT2
    table.  Match the retained harmonics globally and add ``E(anh)-E(harm)``.
    An explicit one-based ``mode_map`` can replace automatic matching.
    """
    if (not isinstance(parsed_vpt2, dict)
            or parsed_vpt2.get('kind') != 'gaussian_vpt2'):
        raise ValueError('Expected a parsed Gaussian VPT2 result.')
    if parsed_vpt2.get('review_required'):
        raise ValueError('VPT2 fundamentals require explicit quality review.')
    harmonic = parsed_vpt2.get('harmonic_fundamentals_cm_inverse')
    anharmonic = parsed_vpt2.get('anharmonic_fundamentals_cm_inverse')
    if (not isinstance(harmonic, list) or not isinstance(anharmonic, list)
            or not harmonic or len(harmonic) != len(anharmonic)
            or any(not isinstance(value, (int, float))
                   or not math.isfinite(value)
                   for value in (*harmonic, *anharmonic))):
        raise ValueError('VPT2 fundamentals are missing or invalid.')
    target = [float(value) for value in getattr(species, 'reduced_freqs', ())]
    if not target or any(not math.isfinite(value) or value <= 0. for value in target):
        raise ValueError('KinBot harmonic MESS frequencies are missing or invalid.')
    if len(target) > len(harmonic):
        raise ValueError(f'VPT2 has {len(harmonic)} fundamentals; KinBot needs '
                         f'{len(target)}.')
    if mode_map is None:
        from scipy.optimize import linear_sum_assignment
        costs = [[abs(target_value - source_value) for source_value in harmonic]
                 for target_value in target]
        rows, columns = linear_sum_assignment(costs)
        mapping = [None] * len(target)
        for row, column in zip(rows, columns):
            mapping[row] = int(column)
    else:
        if (not isinstance(mode_map, (list, tuple))
                or len(mode_map) != len(target)
                or any(not isinstance(value, int) for value in mode_map)):
            raise ValueError('VPT2 mode map must give one integer mode per '
                             'retained KinBot frequency.')
        mapping = [value - 1 for value in mode_map]
        if (min(mapping) < 0 or max(mapping) >= len(harmonic)
                or len(set(mapping)) != len(mapping)):
            raise ValueError('VPT2 mode map is out of range or repeats a mode.')
    differences = [abs(target[index] - harmonic[source])
                   for index, source in enumerate(mapping)]
    if max(differences) > match_tolerance_cm_inverse:
        raise ValueError('VPT2-to-harmonic mode match exceeds '
                         f'{match_tolerance_cm_inverse:g} cm-1; provide a '
                         'reviewed explicit mode map or tolerance.')
    corrected = [target[index] + anharmonic[source] - harmonic[source]
                 for index, source in enumerate(mapping)]
    if any(not math.isfinite(value) or value <= 0. for value in corrected):
        raise ValueError('VPT2 correction gives a nonpositive MESS frequency.')
    species.anl_thermochemistry_frequencies = tuple(corrected)
    species.anl_thermochemistry_frequency_source = {
        'model': 'harmonic_plus_vpt2_mode_shifts',
        'method': parsed_vpt2['method'],
        'basis': parsed_vpt2['basis'],
        'dispersion': parsed_vpt2['dispersion'],
        'mode_map': [value + 1 for value in mapping],
        'maximum_harmonic_mismatch_cm_inverse': max(differences),
    }
    return species.anl_thermochemistry_frequencies


def attach_formation_enthalpy(species, formation):
    """Associate a CBH 0 K formation value with its accepted species energy."""
    from kinbot.anl.cbh import FormationEnthalpy

    final = getattr(species, 'final_zero_k_energy', None)
    if not isinstance(formation, FormationEnthalpy) or final is None:
        raise ValueError('A CBH formation value requires an accepted 0 K energy.')
    from kinbot.anl.atct import _canonical_smiles

    if (_canonical_smiles(formation.target_smiles) != _canonical_smiles(final.smiles)
            or formation.method != final.method
            or formation.energy_sources.get(formation.target_smiles) != final.source):
        raise ValueError('CBH formation value does not match this species energy.')
    species.formation_enthalpy_0k = formation
    return formation


def formation_enthalpy_0k_kj_mol(species) -> float:
    """Return a validated CBH formation enthalpy for a stable species."""
    from kinbot.anl.cbh import FormationEnthalpy
    from kinbot.anl.atct import _canonical_smiles

    formation = getattr(species, 'formation_enthalpy_0k', None)
    final = getattr(species, 'final_zero_k_energy', None)
    if not isinstance(formation, FormationEnthalpy) or not isinstance(final, ZeroKEnergy):
        raise ValueError(f'{species.name}: CBH/ANL MESS energy needs an accepted '
                         '0 K formation enthalpy.')
    if (_canonical_smiles(formation.target_smiles) != _canonical_smiles(final.smiles)
            or formation.method != final.method
            or formation.energy_sources.get(formation.target_smiles) != final.source
            or not math.isfinite(formation.formation_0k_kj_mol)):
        raise ValueError(f'{species.name}: invalid or mismatched CBH formation enthalpy.')
    species_zero_k_hartree(species)
    return formation.formation_0k_kj_mol


def attach_zero_k_barrier(reaction, reactant, transition_state=None) -> ZeroKBarrier:
    """Freeze a same-method ANL E0 difference on a KinBot reaction."""
    transition_state = transition_state or reaction.ts
    reactant_energy = getattr(reactant, 'final_zero_k_energy', None)
    transition_energy = getattr(transition_state, 'final_zero_k_energy', None)
    if not isinstance(reactant_energy, ZeroKEnergy) or not isinstance(
            transition_energy, ZeroKEnergy):
        raise ValueError('A CBH/ANL transition-state barrier needs accepted 0 K '
                         'energies for the reactant and transition state.')
    species_zero_k_hartree(reactant)
    species_zero_k_hartree(transition_state)
    if reactant_energy.method != transition_energy.method:
        raise ValueError('A CBH/ANL transition-state barrier must use one accepted '
                         'method for the reactant and transition state.')
    barrier = ZeroKBarrier(
        reactant_energy.smiles, transition_energy.smiles,
        (transition_energy.hartree - reactant_energy.hartree)
        * HARTREE_TO_KJ_MOL,
        reactant_energy.method, reactant_energy.source, transition_energy.source)
    reaction.zero_k_barrier = barrier
    return barrier


def reaction_barrier_0k_kj_mol(reaction, reactant) -> float:
    """Read an explicit barrier, or derive one from a consistent ANL pair."""
    barrier = getattr(reaction, 'zero_k_barrier', None)
    if barrier is None:
        barrier = attach_zero_k_barrier(reaction, reactant)
    if not isinstance(barrier, ZeroKBarrier) or not math.isfinite(barrier.barrier_0k_kj_mol):
        raise ValueError(f'{reaction.instance_name}: invalid 0 K barrier.')
    reactant_energy = getattr(reactant, 'final_zero_k_energy', None)
    transition_energy = getattr(reaction.ts, 'final_zero_k_energy', None)
    if (not isinstance(reactant_energy, ZeroKEnergy)
            or not isinstance(transition_energy, ZeroKEnergy)
            or barrier.reactant_smiles != reactant_energy.smiles
            or barrier.transition_state_smiles != transition_energy.smiles
            or barrier.method != reactant_energy.method
            or barrier.method != transition_energy.method
            or barrier.reactant_source != reactant_energy.source
            or barrier.transition_state_source != transition_energy.source):
        raise ValueError(f'{reaction.instance_name}: 0 K barrier provenance no longer '
                         'matches the accepted reactant and transition-state energies.')
    return barrier.barrier_0k_kj_mol
