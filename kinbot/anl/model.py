"""Provenance-checked arithmetic for complete 0 K composite energies.

The dispatcher supplies raw QC results. A later workflow layer must build
components from verified task records (including any CBS extrapolations)
before passing them to this evaluator.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import re
from typing import Mapping


_SHA256 = re.compile(r'[0-9a-f]{64}\Z')


class IncompleteRecipeError(ValueError):
    """A composite expression is missing a required validated component."""


@dataclass(frozen=True)
class ComponentResult:
    key: str
    value_hartree: float
    quantity: str
    method: str
    basis: str
    backend: str
    state_id: str
    charge: int
    multiplicity: int
    geometry_sha256: str | None
    source_sha256: str
    source: str
    settings: Mapping[str, object] = field(default_factory=dict)
    review_required: bool = False


@dataclass(frozen=True)
class ComponentRequirement:
    key: str
    quantity: str
    method: str
    basis: str
    geometry_role: str
    backends: tuple[str, ...]
    settings: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class ExpressionTerm:
    component: str
    coefficient: int = 1

    def __post_init__(self):
        if self.coefficient not in (-1, 1):
            raise ValueError('Composite coefficients must be +1 or -1.')


@dataclass(frozen=True)
class CompositeResult:
    recipe: str
    recipe_version: int
    electronic_hartree: float
    zero_point_hartree: float
    zero_k_hartree: float
    components: Mapping[str, ComponentResult]
    electronic_terms: tuple[ExpressionTerm, ...]
    zero_point_terms: tuple[ExpressionTerm, ...]


@dataclass(frozen=True)
class CompositeRecipe:
    name: str
    version: int
    requirements: tuple[ComponentRequirement, ...]
    electronic_terms: tuple[ExpressionTerm, ...]
    zero_point_terms: tuple[ExpressionTerm, ...]

    def __post_init__(self):
        keys = [item.key for item in self.requirements]
        if len(set(keys)) != len(keys):
            raise ValueError('Recipe component keys must be unique.')
        terms = self.electronic_terms + self.zero_point_terms
        if len({item.component for item in terms}) != len(terms):
            raise ValueError('A composite component may appear only once.')
        if {item.component for item in terms} != set(keys):
            raise ValueError('Every recipe requirement needs exactly one term.')

    def evaluate(self, components: Mapping[str, ComponentResult], *,
                 state_id: str, charge: int, multiplicity: int,
                 geometry_hashes: Mapping[str, str]) -> CompositeResult:
        """Assemble a complete expression, rejecting stale or mismatched data."""
        if not state_id or isinstance(charge, bool) or not isinstance(charge, int):
            raise ValueError('A valid electronic state and integer charge are required.')
        if isinstance(multiplicity, bool) or not isinstance(multiplicity, int) \
                or multiplicity < 1:
            raise ValueError('Multiplicity must be a positive integer.')
        missing = sorted({item.key for item in self.requirements} - components.keys())
        if missing:
            raise IncompleteRecipeError('Missing required components: ' + ', '.join(missing))
        selected = {}
        for expected in self.requirements:
            found = components[expected.key]
            if not isinstance(found, ComponentResult) or found.key != expected.key:
                raise ValueError(f'{expected.key}: component identity differs.')
            if not math.isfinite(found.value_hartree):
                raise ValueError(f'{expected.key}: nonfinite energy.')
            if found.review_required:
                raise IncompleteRecipeError(f'{expected.key}: native quality review is required.')
            if (found.state_id, found.charge, found.multiplicity) != (
                    state_id, charge, multiplicity):
                raise ValueError(f'{expected.key}: electronic state differs.')
            if (found.quantity, found.method, found.basis) != (
                    expected.quantity, expected.method, expected.basis):
                raise ValueError(f'{expected.key}: method, basis, or quantity differs.')
            backends = expected.backends
            if expected.key in ('hoe_high', 'hoe_tz_high', 'hoe_dz_low'):
                backends = ('cfour',) if multiplicity == 1 else ('mrcc',)
            if found.backend not in backends:
                raise ValueError(f'{expected.key}: backend differs.')
            if expected.geometry_role == 'state':
                if found.geometry_sha256 is not None:
                    raise ValueError(f'{expected.key}: state-only correction has a geometry.')
            else:
                geometry = geometry_hashes.get(expected.geometry_role)
                if not isinstance(geometry, str) or not _SHA256.fullmatch(geometry):
                    raise ValueError(f'Missing valid {expected.geometry_role} geometry hash.')
                if found.geometry_sha256 != geometry:
                    raise ValueError(f'{expected.key}: geometry differs.')
            if any(found.settings.get(key) != value
                   for key, value in expected.settings.items()):
                raise ValueError(f'{expected.key}: calculation settings differ.')
            if (not found.source or not isinstance(found.source_sha256, str)
                    or not _SHA256.fullmatch(found.source_sha256)):
                raise ValueError(f'{expected.key}: verified source provenance is missing.')
            selected[expected.key] = found
        electronic = math.fsum(term.coefficient * selected[term.component].value_hartree
                               for term in self.electronic_terms)
        zero_point = math.fsum(term.coefficient * selected[term.component].value_hartree
                               for term in self.zero_point_terms)
        if not all(map(math.isfinite, (electronic, zero_point,
                                       electronic + zero_point))):
            raise ValueError('Composite energy is nonfinite.')
        return CompositeResult(
            recipe=self.name, recipe_version=self.version,
            electronic_hartree=electronic, zero_point_hartree=zero_point,
            zero_k_hartree=math.fsum((electronic, zero_point)),
            components=selected, electronic_terms=self.electronic_terms,
            zero_point_terms=self.zero_point_terms,
        )
