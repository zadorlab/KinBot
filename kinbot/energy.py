"""One 0 K energy convention for legacy and composite KinBot results."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math


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
