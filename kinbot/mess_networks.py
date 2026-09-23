"""Keep independent groups of wells in independent MESS calculations.

Separated products are exits, not connections between wells. This module only
groups the already selected models; it does not select reactions or change
energies, symmetry factors, or conformer sums.
"""
import json
import logging
from pathlib import Path
import re

logger = logging.getLogger('KinBot')
_HEADER = re.compile(r'^[ \t]*(Well|Bimolecular|Barrier)[ \t]+([^\n]+)\n', re.M)
_REACTANT = re.compile(r'^([ \t]*Reactant[ \t]+)(\S+)([^\n]*)$', re.M)


def split_model(contents):
    """Separate KinBot's top-level model blocks, keeping their own comments."""
    body, separator, tail = contents.rpartition('End ! end kinetics')
    headers = list(_HEADER.finditer(body))
    if not separator or not headers:
        raise ValueError('Missing completed MESS model or population blocks')
    starts = []
    for match in headers:
        start = offset = match.start()
        for line in reversed(body[:start].splitlines(keepends=True)):
            if line.strip() and not line.lstrip().startswith('!'):
                break
            offset -= len(line)
            if (line.startswith('! kinbot_stereopath ')
                    or line.startswith('! kinbot_racemic_population ')):
                start = offset
        starts.append(start)
    return (body[:starts[0]],
            [body[start:starts[i+1] if i+1 < len(starts) else len(body)]
             for i, start in enumerate(starts)], separator + tail)


def connected_models(contents):
    """Return all supplied well groups, primary Reactant's group first.

The caller has already selected the wells to include. Preserve every such well;
choose the first supplied well as Reactant in each additional group. Products
shared by independent groups are copied into each, without combining rates.
"""
    prefix, blocks, ending = split_model(contents)
    reactants = _REACTANT.findall(prefix)
    if len(reactants) != 1:
        raise ValueError('Expected one Reactant in the completed MESS header')
    root = reactants[0][1]
    populations, barriers, names = {}, [], set()
    for block in blocks:
        header = _HEADER.search(block)
        kind = header[1]
        words = header[2].split('!')[0].split()
        if not words or words[0] in names:
            raise ValueError('Duplicate or missing MESS object name')
        names.add(words[0])
        if kind == 'Barrier':
            if len(words) != 3:
                raise ValueError('MESS barrier must declare two endpoints: ' + header[2])
            barriers.append((words, block))
        else:
            populations[words[0]] = (kind, block)
    wells = {name: set() for name, (kind, _) in populations.items() if kind == 'Well'}
    if root not in wells:
        raise ValueError('MESS Reactant is not a supplied well: ' + root)
    for (_, left, right), _ in barriers:
        if left not in populations or right not in populations:
            raise ValueError(f'MESS barrier refers to missing endpoint: {left} {right}')
        if left not in wells and right not in wells:
            raise ValueError('MESS barrier has no well endpoint: ' + left + ' ' + right)
        if left in wells and right in wells:
            wells[left].add(right)
            wells[right].add(left)

    remaining, models = set(wells), []
    for first in [root] + [name for name in wells if name != root]:
        if first not in remaining:
            continue
        group, pending = set(), [first]
        while pending:
            well = pending.pop()
            if well in group:
                continue
            group.add(well)
            pending.extend(wells[well] - group)
        remaining -= group
        routes = [(words, block) for words, block in barriers
                  if words[1] in group or words[2] in group]
        exits = {name for words, _ in routes for name in words[1:]
                 if name not in group}
        selected = group | exits
        model_prefix = _REACTANT.sub(lambda m: m[1] + first + m[3], prefix)
        text = (model_prefix + ''.join(block for name, (_, block) in populations.items()
                                      if name in selected)
                + ''.join(block for _, block in routes) + ending)
        models.append(dict(contents=text, reactant=first,
                           wells=[name for name in wells if name in group],
                           products=[name for name in populations if name in exits],
                           barriers=[words[0] for words, _ in routes]))
    # Preserve the existing complete input exactly when no grouping is needed.
    if len(models) == 1 and len(populations) == len(models[0]['wells']) + len(models[0]['products']):
        models[0]['contents'] = contents
    return models


def write_network_inputs(writer, contents, uq_index):
    """Publish inputs plus a list of current jobs; never discover jobs by glob."""
    models = connected_models(contents)
    folder = Path('me')
    folder.mkdir(exist_ok=True)
    if not hasattr(writer, 'mess_jobs'):
        writer.mess_jobs = []
    # Rewriting a UQ sample replaces its job list, not its saved solver outputs.
    writer.mess_jobs = [job for job in writer.mess_jobs if job['uq_index'] != uq_index]
    for number, model in enumerate(models, 1):
        stem = f'mess_{uq_index:04d}' + (f'_group_{number:04d}' if number > 1 else '')
        text = model.pop('contents')
        # All supported launchers run inside me; keep auxiliary outputs distinct
        # between both network groups and uncertainty samples, and relocatable.
        text = re.sub(r'^([ \t]*MicroRateOutput[ \t]+)\S+',
                      lambda m: m[1] + stem + '.micro', text, flags=re.M)
        (folder / (stem + '.inp')).write_text(text)
        status = 'ready' if model['barriers'] else 'no_reactions'
        job = dict(stem=stem, uq_index=uq_index, status=status, **model)
        writer.mess_jobs.append(job)
        if status == 'no_reactions':
            logger.warning('MESS input %s has no reactions; retained for inspection, '
                           'but no rate calculation will be submitted.', stem)
        elif len(models) > 1:
            logger.info('Independent MESS calculation %s: Reactant %s, wells %s',
                        stem, model['reactant'], ', '.join(model['wells']))
    (folder / 'mess_networks.json').write_text(json.dumps(writer.mess_jobs, indent=2) + '\n')
