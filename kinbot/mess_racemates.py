"""Fold configured mirror observations into already weighted MESS populations.

QC objects and intermediate files keep their configured names. This last-step
operation changes network names only; it never adds an optical/site multiplier.
"""
import hashlib
import logging
import re

import numpy as np
from ase import Atoms

from kinbot.mess_mirrors import (_HEADER, _POP, _EDGE, _PATH, _insert,
                                _endpoint_references, _shift_eckart_depths)

_KEYS = re.compile(r'^! kinbot_racemic_keys (\S+) (\S+)$', re.M)
_ROUTE = re.compile(r'^! kinbot_racemic_route (\S+) (\S+) (\S+) (\S+)$', re.M)
_ORBIT = re.compile(r'^! kinbot_racemic_path (\S+)$', re.M)
_RRHO = re.compile(r'^\s*RRHO\s*\n.*?^\s*End ! RRHO[^\n]*', re.M | re.S)
_FAMILY = re.compile(r'^! kinbot_racemic_population (\S+) (\S+)$', re.M)
logger = logging.getLogger('KinBot')


def _harmonize_family_models(populations):
    """Reuse a complete intrinsic ensemble; keep each occurrence's energy zero.

    A fragment's absolute electronic energy is not present in assembled PES
    input. Consequently this is an explicit common-intrinsic-model approximation,
    not a reconstruction of endpoint energies or tunneling depths.
    """
    from kinbot.mess import apply_conformer_shifts
    from kinbot.mess_mirrors import _reflect
    energy = re.compile(r'ZeroEnergy\[kcal/mol\]\s+([-+\d.eE]+)')
    models, output = {}, []
    for original in populations:
        block = apply_conformer_shifts(original)
        # Re-find offsets after each replacement, since ensemble sizes may differ.
        markers = [m.groups() for m in _FAMILY.finditer(block)]
        for ordinal, (own, family) in enumerate(markers):
            marker = list(_FAMILY.finditer(block))[ordinal]
            start, end = marker.end(), len(block)
            if _HEADER.search(block)[1] == 'Bimolecular':
                fragments = list(re.finditer(r'^\s*Fragment[^\n]*\n', block[start:], re.M))
                if not fragments:
                    continue
                end = start + fragments[1].start() if len(fragments) > 1 else end
                start += fragments[0].end()
            segment = block[start:end]
            rrhos = list(_RRHO.finditer(segment))
            if not rrhos:
                continue
            first, last = rrhos[0].start(), rrhos[-1].end()
            union = re.search(r'^[ \t]*Union\b[^\n]*\n', segment[:first], re.M)
            union_end = re.search(r'^[ \t]*End ! Union[^\n]*', segment[last:], re.M)
            if union and union_end:
                first, last = union.start(), last + union_end.end()
            kernel = segment[first:last]
            values = [float(m[1]) for m in energy.finditer(kernel)]
            if not values:
                continue
            ground = min(values)
            normalized = energy.sub(lambda m: 'ZeroEnergy[kcal/mol] ' +
                                    str(round(float(m[1]) - ground, 10)), kernel)
            if family not in models:
                models[family] = (own, normalized)
                continue
            source, selected = models[family]
            if _same_ensemble(selected, normalized):
                continue
            if source != own:
                selected = _reflect(selected)
            selected = energy.sub(lambda m: 'ZeroEnergy[kcal/mol] ' +
                                  str(round(float(m[1]) + ground, 10)), selected)
            message = (f'Racemic species {own}: using the complete intrinsic model from {source}; '
                       'the existing endpoint energy reference is retained. Independent '
                       'mirror calculations disagree; common-model approximation.')
            logger.warning(message)
            replacement = '! WARNING: ' + message + '\n' + selected
            block = block[:start+first] + replacement + block[start+last:]
        output.append(block)
    return output


def _tokens(text):
    """Compare the supplied model, ignoring formatting and provenance comments."""
    return ' '.join(line.split('!')[0] for line in text.splitlines()).split()


def _same_tokens(left, right):
    if len(left) != len(right):
        return False
    for i, (a, b) in enumerate(zip(left, right)):
        try:
            a, b = float(a), float(b)
        except ValueError:
            if a != b:
                return False
        else:
            # Serialized references are rounded to .01 kcal/mol. Do not reuse
            # the 1 kcal/mol geometry-matching tolerance for endpoint energies.
            energy = i and left[i-1].startswith(('ZeroEnergy[', 'GroundEnergy['))
            if not np.isfinite([a, b]).all() or not np.isclose(
                    a, b, atol=.005 if energy else 1.e-6, rtol=0 if energy else 1.e-5):
                return False
    return True


def _rrho_model(block):
    """One RRHO kernel and its statistical weight, independent of orientation."""
    geometries = list(re.finditer(r'Geometry\[angstrom\]\s+(\d+)\s*\n', block))
    for geometry in reversed(geometries):
        n = int(geometry[1])
        lines = block[geometry.end():].splitlines(True)
        rows = [line.split() for line in lines[:n]]
        atoms = Atoms([r[0] for r in rows], positions=[[float(x) for x in r[1:4]] for r in rows])
        # A rigid rotor depends on principal moments, not lab orientation.
        moments = atoms.get_moments_of_inertia()
        marker = 'MassAndMoments ' + ' '.join(map(str, [atoms.get_masses().sum(), *moments])) + '\n'
        if re.search(r'\bRotor\b', block):
            # Rotor kinetic parameters depend on the indexed groups/axes as
            # well as the whole-molecule moments. Equal overall inertia alone
            # is not sufficient. Preserve indexed internal geometry for HIR.
            marker += 'RotorGeometry ' + ' '.join(atoms.get_chemical_symbols()) + ' '
            marker += ' '.join(map(str, atoms.get_all_distances().ravel())) + '\n'
        block = block[:geometry.start()] + marker + ''.join(lines[n:])
    factor = re.search(r'\bSymmetryFactor\s+([-+\d.eE]+)', block)
    if not factor or float(factor[1]) <= 0:
        raise ValueError('Racemic population comparison requires a positive RRHO symmetry factor')
    weight = 1. / float(factor[1])
    block = block[:factor.start()] + block[factor.end():]
    # Mirroring reverses a torsion's positive angular direction. The potential
    # [v0,v1,...] and [v0,...,v1] give the same one-dimensional rotor model.
    # Keep the reference angle fixed; do not fit/shift an unrelated profile.
    potentials = list(re.finditer(r'Potential\[kcal/mol\]\s+(\d+)\s*\n', block))
    for potential in reversed(potentials):
        number = int(potential[1])
        values = list(re.finditer(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?',
                                 block[potential.end():]))[:number]
        if len(values) != number or number == 0:
            raise ValueError('Incomplete rotor potential in racemic population model')
        samples = [float(v[0]) for v in values]
        samples = min(samples, samples[:1] + samples[:0:-1])
        end = potential.end() + values[-1].end()
        block = block[:potential.end()] + ' '.join(map(str, samples)) + block[end:]
    # Reflection may reverse reaction direction. Eckart's pair of well depths
    # describes the same potential in either order; leave selected output intact.
    depths = list(re.finditer(r'WellDepth\[kcal/mol\]\s+([-+\d.eE]+)', block))
    if len(depths) == 2:
        values = sorted(float(m[1]) for m in depths)
        for match, value in reversed(list(zip(depths, values))):
            block = block[:match.start()] + f'WellDepth[kcal/mol] {value}' + block[match.end():]
    return _tokens(block), weight


def _ensemble(block):
    models = []
    for match in _RRHO.finditer(block):
        tokens, weight = _rrho_model(match.group())
        for entry in models:
            if _same_tokens(entry[0], tokens):
                entry[1] += weight
                break
        else:
            models.append([tokens, weight])
    return models


def _same_ensemble(left, right):
    from kinbot.mess import apply_conformer_shifts
    left, right = apply_conformer_shifts(left), apply_conformer_shifts(right)
    a, b = _ensemble(left), _ensemble(right)
    if not a or len(a) != len(b):
        return False
    remaining = list(b)
    for tokens, weight in a:
        for i, (other, w) in enumerate(remaining):
            if _same_tokens(tokens, other) and np.isclose(weight, w, rtol=1.e-8, atol=1.e-8):
                remaining.pop(i)
                break
        else:
            return False
    return True


def _compatible_populations(left, right):
    """Compare whole well/fragment ensembles, not just their lowest energies."""
    from kinbot.mess import apply_conformer_shifts
    left, right = apply_conformer_shifts(left), apply_conformer_shifts(right)
    if _HEADER.search(left)[1] == 'Well':
        return _same_ensemble(left, right)
    if re.search(r'^\s*Dummy\s*$', left, re.M) or re.search(r'^\s*Dummy\s*$', right, re.M):
        # These sinks deliberately have no reverse partition function.
        def dummy_body(text):
            return _tokens(text[_HEADER.search(text).end():]) == ['Dummy']
        return dummy_body(left) and dummy_body(right)
    grounds = lambda s: re.findall(r'GroundEnergy\[kcal/mol\]\s+([-+\d.eE]+)', s)
    if not _same_tokens(['GroundEnergy[kcal/mol]', *grounds(left)],
                        ['GroundEnergy[kcal/mol]', *grounds(right)]):
        return False
    def fragments(text):
        return re.split(r'^\s*Fragment[^\n]*\n', text, flags=re.M)[1:]
    a, remaining = fragments(left), fragments(right)
    if not a or len(a) != len(remaining):
        return False
    for fragment in a:
        for i, other in enumerate(remaining):
            if _RRHO.search(fragment):
                same = _same_ensemble(fragment, other)
            else:
                atom = re.search(r'\bAtom\s*\n.*?End ! Atom', fragment, re.S)
                other_atom = re.search(r'\bAtom\s*\n.*?End ! Atom', other, re.S)
                same = bool(atom and other_atom and _same_tokens(_tokens(atom[0]), _tokens(other_atom[0])))
            if same:
                remaining.pop(i)
                break
        else:
            return False
    return True


def _validate_family_models(populations):
    """The same racemate must use the same intrinsic model in R+H and S+OH.

    Well grounds are global; fragment grounds are relative to the fragment.
    Compare excitation energies and full models here. Complete-population
    comparisons separately check the absolute network energy references.
    """
    from kinbot.mess import apply_conformer_shifts
    models = {}
    energy = re.compile(r'ZeroEnergy\[kcal/mol\]\s+([-+\d.eE]+)')
    for block in populations:
        kind = _HEADER.search(block)[1]
        for marker in _FAMILY.finditer(block):
            model = block
            if kind == 'Bimolecular':
                fragments = list(re.finditer(r'^\s*Fragment[^\n]*\n', block[marker.end():], re.M))
                if not fragments:
                    continue  # A Dummy sink supplies no reverse model.
                start = marker.end() + fragments[0].end()
                end = marker.end() + fragments[1].start() if len(fragments) > 1 else len(block)
                model = block[start:end]
            model = apply_conformer_shifts(model)
            energies = [float(m[1]) for m in energy.finditer(model)]
            if not energies:
                raise ValueError('Missing intrinsic model for racemic population ' + marker[1])
            ground = min(energies)
            model = energy.sub(lambda m: 'ZeroEnergy[kcal/mol] ' + str(round(float(m[1])-ground, 10)), model)
            if marker[2] in models and not _same_ensemble(models[marker[2]], model):
                raise ValueError('Conflicting intrinsic models of the same racemic species across wells/fragments: ' + marker[1])
            models[marker[2]] = model


def fold_racemic_network(prefix, blocks):
    """Collapse global mirror pairs; preserve relative-configuration pathways.

    A_R→B_R and A_R→B_S have different endpoint-pair orbits even when they
    share a reaction-family/site label. Only simultaneous global reflection
    (or reversal of the same reaction) identifies two route observations.
    """
    if not any(_KEYS.search(b) for b in blocks):
        return prefix, blocks
    references = _endpoint_references(blocks)
    blocks = _harmonize_family_models(blocks)
    # A substituted mirror block may be folded away; keep its warning visible.
    warning = re.compile(r'^! WARNING: Racemic species [^\n]*\n', re.M)
    for block in blocks:
        for message in warning.findall(block):
            if message not in prefix:
                prefix += message
    blocks = [warning.sub('', block) for block in blocks]
    populations, barriers, groups, aliases, key_aliases = [], [], {}, {}, {}
    for block in blocks:
        header = _HEADER.search(block)
        name = header[2].split()[0]
        if header[1] == 'Barrier':
            barriers.append(block)
            continue
        marker = _KEYS.search(block)
        if not marker:
            populations.append(block)
            continue
        key = (header[1], tuple(sorted(marker.groups())))
        if key in groups:
            prior = groups[key]
            if not _compatible_populations(prior, block):
                message = (f'Conflicting racemic population models: retaining the complete '
                           f'model {_HEADER.search(prior)[2].split()[0]} instead of {name}. '
                           'The first occurrence (including the requested entrance) sets the model.')
                logger.warning(message)
                position = populations.index(prior)
                prior = '! WARNING: ' + message + '\n' + prior
                populations[position] = groups[key] = prior
            aliases[name] = _HEADER.search(prior)[2].split()[0]
        else:
            groups[key] = block
            populations.append(block)
            aliases[name] = name
        canonical = _POP.search(groups[key])
        if not canonical:
            raise ValueError('Racemic population is missing its configured population metadata')
        for configured in marker.groups():
            key_aliases[configured] = canonical[1]
    # Fragments can recur in different product combinations. Their bookkeeping
    # comments must use one family representative even when the complete
    # bimolecular populations are different (R+H and R+OH, for example).
    _validate_family_models(populations)
    family_keys = {}
    def family_comment(match):
        key, family = match.groups()
        return '! kinbot_racemic_population ' + family_keys.setdefault(family, key) + ' ' + family
    populations = [_FAMILY.sub(family_comment, b) for b in populations]
    def redirect(match):
        return match[1] + aliases.get(match[2], match[2])
    prefix = re.sub(r'(^\s*Reactant\s+)(\S+)', redirect, prefix, flags=re.M)
    retained = {}
    for block in barriers:
        header = _HEADER.search(block)
        words = header[2].split('!')[0].split()
        if len(words) != 3:
            raise ValueError('Racemic network requires a barrier with two named endpoints')
        name, left, right = words
        new_left, new_right = aliases.get(left, left), aliases.get(right, right)
        if new_left == new_right and left != right:
            continue  # R↔S becomes internal motion of the same racemic well.
        original = _ROUTE.search(block)
        labels = _PATH.findall(block)
        if len(labels) > 1:
            raise ValueError('Ambiguous stereochemical path in racemic network')
        old_orbit = _ORBIT.search(block)
        if old_orbit:
            orbit = old_orbit[1]
        elif original:
            a, ma, b, mb = original.groups()
            pair = min(tuple(sorted((a, b))), tuple(sorted((ma, mb))))
            orbit = hashlib.sha256(repr((pair, labels)).encode()).hexdigest()[:24]
            block = _insert(block, '! kinbot_racemic_path ' + orbit)
            if labels:
                block = _PATH.sub('! original stereochemical path ' + labels[0], block)
            block = _insert(block, '! kinbot_stereopath racemic-' + orbit)
        else:
            if (left, right) != (new_left, new_right):
                raise ValueError('Cannot redirect a racemic barrier without original configured endpoint evidence')
            orbit = labels[0] if labels else None
        changes = [references[old] - references[new]
                   if old != new and old in references and new in references else 0.
                   for old, new in ((left, new_left), (right, new_right))]
        block = _shift_eckart_depths(block, changes)
        block = _HEADER.sub(lambda m: f'  Barrier {name} {new_left} {new_right}\n', block, count=1)
        block = _EDGE.sub(lambda m: '! kinbot_mirror_endpoints ' +
                          ' '.join(key_aliases.get(k, k) for k in m.groups()), block)
        # Unclassified/PST routes retain the pre-existing parallel-route guard.
        # Folding is not permission to silently deduplicate those observations.
        key = (tuple(sorted((new_left, new_right))), orbit,
               None if original or old_orbit else len(retained))
        if key in retained:
            # Keep the established lowest-barrier choice within a path class.
            # Keep a whole selected model, including its own Eckart depths.
            def ground(text):
                from kinbot.mess import apply_conformer_shifts
                values = re.findall(r'ZeroEnergy\[kcal/mol\]\s+([-+\d.eE]+)', apply_conformer_shifts(text))
                if not values:
                    raise ValueError('Mirror-duplicate barrier lacks a resolved energy reference')
                return min(map(float, values))
            prior = retained[key]
            if len(_RRHO.findall(prior)) > 1 or len(_RRHO.findall(block)) > 1:
                if not _same_ensemble(prior, block):
                    message = ('Conflicting racemic TS conformer ensembles: retaining the complete '
                               'lowest-barrier ensemble for this same mirror pathway; no extra optical factor.')
                    logger.warning(message)
                    chosen = block if ground(block) < ground(prior) else prior
                    retained[key] = '! WARNING: ' + message + '\n' + chosen
                    continue
            if ground(block) < ground(prior):
                retained[key] = block
        else:
            retained[key] = block
    return prefix, populations + list(retained.values())
