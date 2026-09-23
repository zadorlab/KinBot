"""Complete known global-mirror channels, without new QC jobs or site factors.

Field-free molecular models have identical mirror energies and spectra. Copy
only from an endpoint population already in the assembled network; a specified
chiral entrance does not by itself license its absent mirror population.
"""
import copy
import hashlib
import logging
import re
from kinbot.stereo_identity import canonical_identity
from kinbot.species_routing import routing_name
from kinbot.mess_networks import _HEADER, split_model

_POP = re.compile(r'^! kinbot_population (\S+) (\S+)$', re.M)
_EDGE = re.compile(r'^! kinbot_mirror_endpoints (\S+) (\S+) (\S+) (\S+)$', re.M)
_PATH = re.compile(r'^! kinbot_stereopath (\S+)$', re.M)
_XYZ = re.compile(r'^(\s*[A-Z][a-z]?\s+)([-+\d.eE]+)(\s+[-+\d.eE]+\s+[-+\d.eE]+\s*)$')
logger = logging.getLogger('KinBot')


def population_keys(species, population='specified'):
    own, mirror = [], []
    for point in species:
        identity = getattr(point, 'optical_reference', None) or canonical_identity(point)
        if identity['status'] != 'assigned':
            return None
        key = routing_name(point)
        own.append(key)
        if population == 'racemic' or identity['id'] == identity['mirror_id']:
            mirror.append(key)
        else:
            view = copy.copy(point)
            view.optical_reference = dict(identity, id=identity['mirror_id'], mirror_id=identity['id'])
            mirror.append(routing_name(view))
    return '_'.join(sorted(own)), '_'.join(sorted(mirror))


def _insert(block, comment):
    match = _HEADER.search(block)
    if match is None:
        raise ValueError('Missing MESS population/barrier header')
    return block[:match.end()] + comment + '\n' + block[match.end():]


def annotate_population(block, species, population, *, energy_reference=None):
    keys = population_keys(species, population)
    if keys:
        block = _insert(block, '! kinbot_population ' + ' '.join(keys))
        if population == 'racemic':
            configured = population_keys(species)
            if configured:
                block = _insert(block, '! kinbot_racemic_keys ' + ' '.join(configured))
    if energy_reference is not None:
        block = _insert(block, f'! kinbot_eckart_reference[kcal/mol] {energy_reference}')
    return block


def annotate_endpoints(block, reactants, products, population):
    left, right = population_keys(reactants, population), population_keys(products, population)
    if left and right:
        block = _insert(block, '! kinbot_mirror_endpoints ' + ' '.join(left + right))
        if population == 'racemic':
            configured_left, configured_right = population_keys(reactants), population_keys(products)
            if configured_left and configured_right:
                block = _insert(block, '! kinbot_racemic_route ' +
                                ' '.join(configured_left + configured_right))
    return block


def _reflect(block):
    lines = []
    for line in block.splitlines(True):
        match = _XYZ.match(line.rstrip('\n'))
        lines.append((match[1] + f'{-float(match[2]):.12f}' + match[3] + '\n') if match else line)
    block = _POP.sub(lambda m: f'! kinbot_population {m[2]} {m[1]}', ''.join(lines))
    block = _EDGE.sub(lambda m: f'! kinbot_mirror_endpoints {m[2]} {m[1]} {m[4]} {m[3]}', block)
    # Reflection reverses the positive torsional angle. Keep point zero and
    # reverse the remaining samples; energies and rotor levels are unchanged.
    potentials = list(re.finditer(r'Potential\[kcal/mol\]\s+(\d+)[ \t]*\n', block))
    for potential in reversed(potentials):
        values = list(re.finditer(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?',
                                 block[potential.end():]))[:int(potential[1])]
        if len(values) != int(potential[1]) or not values:
            continue
        samples = [value[0] for value in values]
        finish = potential.end() + values[-1].end()
        block = block[:potential.end()] + ' '.join(samples[:1] + samples[:0:-1]) + block[finish:]
    return block


def _population_energies(block):
    from kinbot.mess import apply_conformer_shifts
    resolved = apply_conformer_shifts(block)
    return sorted((match[1], float(match[2])) for match in re.finditer(
        r'^[ \t]*((?:Zero|Ground)Energy\[[^\]]+\])[ \t]+([-+\d.eE]+)', resolved, re.M))


def _endpoint_references(blocks):
    """Read Eckart references, then saved parents, with legacy grounds as fallback."""
    from kinbot.mess import apply_conformer_shifts
    references = {}
    # Eckart depths use the selected parent, which need not be the lowest
    # retained MC conformer. Read that existing convention, not an MC minimum.
    for block in blocks:
        header = _HEADER.search(block)
        if header[1] != 'Barrier':
            continue
        endpoints = header[2].split()[1:3]
        resolved = apply_conformer_shifts(block)
        depths = list(re.finditer(r'WellDepth\[kcal/mol\]\s+([-+\d.eE]+)', resolved))
        energies = (list(re.finditer(r'ZeroEnergy\[kcal/mol\]\s+([-+\d.eE]+)',
                                   resolved[:depths[0].start()])) if depths else [])
        if energies and len(depths) >= 2:
            energy = energies[-1]  # Earlier submerged MC conformers may have no Eckart block.
            for name, depth in zip(endpoints, depths):
                references.setdefault(name, float(energy[1]) - float(depth[1]))
    for block in blocks:
        header = _HEADER.search(block)
        if header[1] == 'Barrier':
            continue
        # A selected MC parent need not be its ensemble minimum and may have
        # no attached tunneling saddle. Existing Eckart data take precedence:
        # some channels deliberately use an inner complex or another QC level.
        saved = re.search(r'^! kinbot_eckart_reference\[kcal/mol\]\s+([-+\d.eE]+)$', block, re.M)
        if saved:
            references.setdefault(header[2].split()[0], float(saved[1]))
        kind = 'ZeroEnergy' if header[1] == 'Well' else 'GroundEnergy'
        values = [value for key, value in _population_energies(block) if key.startswith(kind)]
        if values:
            references.setdefault(header[2].split()[0], min(values))
    return references


def _shift_eckart_depths(block, changes):
    """Change endpoint references without changing the selected TS model."""
    from kinbot.mess import apply_conformer_shifts
    if not any(changes):
        return block
    block = apply_conformer_shifts(block)
    def update_tunneling(match):
        tunnel = match[0]
        depths = list(re.finditer(r'(WellDepth\[kcal/mol\]\s+)([-+\d.eE]+)', tunnel))
        if len(depths) != 2:
            return tunnel
        values = [round(float(depth[2]) + shift, 10) for depth, shift in zip(depths, changes)]
        if min(values) <= 0.:
            return '! barrier is submerged or has zero serialized tunneling depth'
        for depth, value in reversed(list(zip(depths, values))):
            tunnel = tunnel[:depth.start(2)] + str(value) + tunnel[depth.end(2):]
        return re.sub(r'(CutoffEnergy\[kcal/mol\]\s+)[-+\d.eE]+',
                      lambda m: m[1] + str(round(min(values), 10)), tunnel)
    return re.sub(r'Tunneling\s+Eckart\b.*?^[ \t]*End\b[^\n]*',
                  update_tunneling, block, flags=re.M | re.S)


def _reuse_mirror_populations(blocks):
    """Use one complete calculated model for each explicitly supplied mirror pair.

    The first supplied population is the reference. Its entire RRHO/HIR or MC
    model is reflected, never matched conformer-by-conformer to another search.
    Other pathways remain present. Their Eckart depths follow the substituted
    endpoint reference; the TS calculations themselves remain unchanged.
    """
    populations = {}
    for block in blocks:
        header, marker = _HEADER.search(block), _POP.search(block)
        if header[1] != 'Barrier' and marker:
            populations[marker[1]] = (header[2].split()[0], block)
    references = _endpoint_references(blocks)
    replacements, shifts, seen = {}, {}, set()
    for own, (source_name, source) in populations.items():
        mirror = _POP.search(source)[2]
        if own in seen or own == mirror or mirror not in populations:
            continue
        seen.update((own, mirror))
        target_name, target = populations[mirror]
        # Minimal/legacy blocks without a statistical model have nothing to copy.
        if not _population_energies(source) or not _population_energies(target):
            continue
        marker = f'! reflected population model from {source_name}; complete calculation reused'
        if marker in target:
            continue
        reflected = _HEADER.sub(lambda m: _HEADER.search(target)[0], _reflect(source), count=1)
        replacements[target_name] = _insert(reflected, marker)
        shifts[target_name] = references.get(target_name, 0.) - references.get(source_name, 0.)
        logger.info('Mirror population %s uses the complete model from %s; '
                    'independent calculation files are retained.', target_name, source_name)
    output = []
    for block in blocks:
        header = _HEADER.search(block)
        words = header[2].split()
        if words[0] in replacements and header[1] != 'Barrier':
            block = replacements[words[0]]
        elif header[1] == 'Barrier' and any(name in shifts for name in words[1:3]):
            block = _shift_eckart_depths(block, [shifts.get(name, 0.) for name in words[1:3]])
        output.append(block)
    return output


def _separate_routes(block):
    """Unpack only our outer path Union; preserve nested MC RRHO Unions."""
    if not re.search(r'Union ! \d+ stereochemical pathways', block):
        return [block]
    parts = re.split(r'^! pathway source: (Barrier[^\n]+)\n', block, flags=re.M)
    routes = []
    for i in range(1, len(parts), 2):
        content = parts[i+1].split('    End ! stereochemical pathways')[0].rstrip('\n') + '\n'
        routes.append('  ' + parts[i] + '\n' + content)
    return routes


def complete_mirror_channels(contents):
    """Complete missing classified mirror routes, then combine parallel paths.

    The endpoint-pair/path key retains ordinary lowest-barrier selection while
    preserving each demonstrated pathway, including a partly explored mirror
    side. No independent fragment racemates are invented.
    """
    from kinbot.mess import union_stereochemical_barriers
    body, separator, tail = contents.rpartition('End ! end kinetics')
    if not separator:
        return contents
    headers = list(_HEADER.finditer(body))
    if not headers:
        return contents
    prefix, blocks, ending = split_model(contents)
    blocks = [route for block in blocks for route in _separate_routes(block)]
    from kinbot.mess_racemates import fold_racemic_network
    prefix, blocks = fold_racemic_network(prefix, blocks)
    blocks = _reuse_mirror_populations(blocks)
    def finish():
        populations = [b for b in blocks if _HEADER.search(b)[1] != 'Barrier']
        barriers = [b for b in blocks if _HEADER.search(b)[1] == 'Barrier']
        return prefix + ''.join(populations + union_stereochemical_barriers(barriers)) + ending
    if not any(m[1] != m[2] or m[3] != m[4] for m in _EDGE.finditer(contents)):
        return finish()
    populations, edges, templates, names = {}, set(), [], set()
    def route_id(block):
        labels = _PATH.findall(block)
        if len(labels) > 1:
            raise ValueError('A mirror route has ambiguous stereochemical path metadata')
        return labels[0] if labels else None
    for block in blocks:
        header = _HEADER.search(block)
        words = header[2].split('!')[0].split()
        names.add(words[0])
        if header[1] == 'Barrier':
            if len(words) >= 3:
                edges.add((frozenset(words[1:3]), route_id(block)))
            if _EDGE.search(block):
                templates.append(block)
        else:
            marker = _POP.search(block)
            if marker:
                if marker[1] in populations and populations[marker[1]][0] != words[0]:
                    raise ValueError('A configured MESS population has two different names')
                populations[marker[1]] = (words[0], block)
    by_name = {name: _POP.search(block).groups() for name, block in populations.values()}
    for block in blocks:
        header = _HEADER.search(block)
        words = header[2].split('!')[0].split()
        if header[1] == 'Barrier' and not _EDGE.search(block) and len(words) == 3:
            if words[1] in by_name and words[2] in by_name:
                templates.append(_insert(block, '! kinbot_mirror_endpoints ' +
                                         ' '.join(by_name[words[1]] + by_name[words[2]])))
    def name_for(key, kind):
        name = kind + '_mirror_' + hashlib.sha256(key.encode()).hexdigest()[:16]
        if name in names:
            raise ValueError('Derived mirror name collides with a supplied MESS object')
        names.add(name)
        return name
    changed = True
    while changed:
        changed = False
        for block in templates:
            left, ml, right, mr = _EDGE.search(block).groups()
            if left not in populations or right not in populations:
                raise ValueError('Mirror endpoint metadata lacks its declared population block')
            # Separated products are exits, never a route into another well.
            if not any(key in populations and _HEADER.search(populations[key][1])[1] == 'Well'
                       for key in (ml, mr)):
                continue
            endpoints = frozenset(populations[key][0] for key in (ml, mr) if key in populations)
            if ml in populations and mr in populations and (endpoints, route_id(block)) in edges:
                continue
            for own, mirror in ((left, ml), (right, mr)):
                if mirror in populations:
                    continue
                source_name, source = populations[own]
                name = name_for(mirror, _HEADER.search(source)[1].lower())
                derived = _HEADER.sub(lambda m: f'  {m[1]} {name}\n', _reflect(source), count=1)
                derived = _insert(derived, f'! derived by global reflection from {source_name}; no additional QC calculation')
                derived = _insert(derived, f'! reflected population model from {source_name}; complete calculation reused')
                populations[mirror] = (name, derived)
                blocks.append(derived)
                changed = True
            endpoints = frozenset((populations[ml][0], populations[mr][0]))
            source_name = _HEADER.search(block)[2].split()[0]
            name = name_for(source_name + '|' + ml + '|' + mr, 'ts')
            derived = _HEADER.sub(lambda m: f'  Barrier {name} {populations[ml][0]} {populations[mr][0]}\n', _reflect(block), count=1)
            derived = _insert(derived, f'! derived by global reflection from {source_name}; one specified channel')
            blocks.append(derived)
            edges.add((endpoints, route_id(block)))
            changed = True
    return finish()
