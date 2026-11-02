"""First-passage events in a supplied kinetic network (no kinetic model fitting).

Rates must use a common clock within each source well. A common positive
well-specific multiplier cancels, so native MESS N(E) can be used for isolated
microcanonical branching probabilities, but not for lifetimes.
"""
from decimal import Decimal, localcontext
import re


def channel_selection_groups(contents, *, population_families=None):
    """Group work by explicit mirror or supplied population-family membership.

    This is only a selection grouping, never a lumped kinetic model. Explicit
    mirror metadata takes precedence over differing legacy display SMILES.
    Optional families can conservatively retain related conformers together;
    supplying them does not establish rapid conformer equilibration.
    """
    from kinbot.mess_networks import split_model, _HEADER
    blocks = split_model(contents)[1]
    parents, population, edges = {}, {}, {}
    for block in blocks:
        header = _HEADER.search(block); words = header[2].split('!')[0].split()
        if header[1] == 'Barrier':
            edges[words[0]] = words[1:3]
        else:
            parents[words[0]] = words[0]
            marker = re.search(r'^! kinbot_population (\S+) (\S+)$', block, re.M)
            if marker:
                population[marker[1]] = (words[0], marker[2])
    def root(name):
        while parents[name] != name:
            parents[name] = parents[parents[name]]
            name = parents[name]
        return name
    def union(a, b):
        a, b = sorted((root(a), root(b)))
        parents[b] = a
    for name, mirror in population.values():
        if mirror in population:
            union(name, population[mirror][0])
    families = {}
    for name, family in (population_families or {}).items():
        if name not in parents:
            raise ValueError(f'Unknown population in selection family: {name}')
        if family in families:
            union(name, families[family])
        else:
            families[family] = name
    groups = {}
    for name, endpoints in edges.items():
        key = tuple(sorted(root(n) for n in endpoints))
        groups.setdefault(key, []).append(name)
    return groups


def event_hit_probability(wells, rates, initial, event_rates, *, precision=80):
    """Probability of traversing ANY marked channel before leaving the wells.

    ``event_rates`` contains the marked contribution to each directed total
    rate. Marking a group simultaneously computes its union, not the sum of
    individual ever-hit probabilities (which double counts recrossing).
    Unreachable events and closed communicating classes correctly contribute
    zero. Multiple parallel channels can mark only part of a total rate.
    """
    wells = list(wells)
    if len(wells) != len(set(wells)) or initial not in wells:
        raise ValueError('Supply unique wells including the initial well.')
    if precision < 20:
        raise ValueError('At least 20 decimal digits are required.')
    with localcontext() as ctx:
        ctx.prec = precision
        raw = {key: Decimal(str(value)) for key, value in rates.items()}
        marked = {key: Decimal(str(value)) for key, value in event_rates.items()}
        for (source, target), value in raw.items():
            if source not in wells or source == target or not value.is_finite() or value < 0:
                raise ValueError('Rates must be finite nonnegative transitions out of wells.')
        for key, value in marked.items():
            if not value.is_finite() or not 0 <= value <= raw.get(key, Decimal(0)):
                raise ValueError('Event rate must be a nonnegative part of the total rate.')
        # Only wells from which a marked edge can still be reached are unknowns.
        # The complement has probability zero, including closed traps.
        active = {source for (source, _), value in marked.items() if value > 0}
        while True:
            added = {source for (source, target), value in raw.items()
                     if target in active and value > marked.get((source, target), 0)} - active
            if not added:
                break
            active.update(added)
        if initial not in active:
            return 0.
        active = [well for well in wells if well in active]
        index = {name: i for i, name in enumerate(active)}
        n = len(active)
        matrix = [[Decimal(0) for _ in range(n+1)] for _ in active]
        for (source, target), value in raw.items():
            if source not in index:
                continue
            i = index[source]; event = marked.get((source, target), 0)
            matrix[i][i] += value
            matrix[i][-1] += event
            if target in index:
                matrix[i][index[target]] -= value-event
        for j in range(n):
            pivot = max(range(j, n), key=lambda i: abs(matrix[i][j]))
            matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
            divisor = matrix[j][j]
            if not divisor:
                raise ArithmeticError('Singular first-passage system at requested precision.')
            matrix[j] = [value/divisor for value in matrix[j]]
            for i in range(n):
                if i != j:
                    factor = matrix[i][j]
                    matrix[i] = [a-factor*b for a, b in zip(matrix[i], matrix[j])]
        result = matrix[index[initial]][-1]
        if not Decimal('-1e-15') <= result <= Decimal('1.000000000000001'):
            raise ArithmeticError('First-passage probability outside [0, 1].')
        return float(max(Decimal(0), min(Decimal(1), result)))
