"""Explicit alternative guesses for a bounded reaction-search revision."""
from copy import deepcopy
import hashlib
import json
import math
import re

import numpy as np


def validate(guesses):
    if not isinstance(guesses, dict):
        raise ValueError('reaction_search_guesses must be an object.')
    for name, selection in guesses.items():
        if not isinstance(name, str) or not re.fullmatch(
                r'[^/]+_R_Addition_MultipleBond_\d+_\d+_\d+', name):
            raise ValueError('Alternative guesses require an exact three-atom R_Addition_MultipleBond search.')
        variants = selection if isinstance(selection, list) else [selection]
        if not variants:
            raise ValueError('An alternative guess list cannot be empty.')
        for options in variants:
            if (not isinstance(options, dict) or set(options) != {'bond_scale', 'angle'}
                    or options['angle'] != 'parent'):
                raise ValueError('Alternative guesses require bond_scale and angle="parent".')
            scale = options['bond_scale']
            if isinstance(scale, bool) or not isinstance(scale, (int, float)) or not math.isfinite(scale) or scale <= 1:
                raise ValueError('The reaction-guess bond_scale must be finite and greater than one.')
        if len({alternative_name(name, v) for v in variants}) != len(variants):
            raise ValueError('Alternative guesses must be distinct.')
    return deepcopy(guesses)


def targets(reaction):
    source = reaction.instance_name.rsplit('_guess_', 1)[0]
    selection = reaction.par.get('reaction_search_guesses', {}).get(source)
    if selection is None:
        return None
    validate({source: selection})
    variants = selection if isinstance(selection, list) else [selection]
    options = next((v for v in variants if reaction.instance_name == alternative_name(source, v)), None)
    if options is None:
        return None
    if reaction.family_name != 'radditionmultiplebond' or len(reaction.instance) != 3:
        raise ValueError('Alternative guess does not apply to this reaction family/instance.')
    a, b, c = reaction.instance
    geom = np.asarray(reaction.species.geom, dtype=float)
    left, right = geom[a] - geom[b], geom[c] - geom[b]
    lengths = np.linalg.norm(left), np.linalg.norm(right)
    if (not np.isfinite(geom).all() or min(lengths) <= 0
            or reaction.species.bond[b][c] <= 0):
        raise ValueError('Alternative guess requires a finite parent geometry and an existing breaking bond.')
    angle = math.degrees(math.acos(float(np.clip(np.dot(left, right) / math.prod(lengths), -1., 1.))))
    return float(lengths[1] * options['bond_scale']), angle


def alternative_name(source, options):
    digest = hashlib.sha256(json.dumps(options, sort_keys=True).encode()).hexdigest()[:12]
    return source + '_guess_' + digest


def add_alternatives(finder):
    """Give an existing discovered motif an independent, durable attempt."""
    from kinbot.reactions.reac_R_Addition_MultipleBond import RAdditionMultipleBond
    guesses = validate(finder.par.get('reaction_search_guesses', {}))
    species = finder.species
    for index, source in enumerate(list(species.reac_name)):
        if source not in guesses:
            continue
        selection = guesses[source]
        variants = selection if isinstance(selection, list) else [selection]
        for options in variants:
            name = alternative_name(source, options)
            if name in species.reac_name:
                continue
            obj = RAdditionMultipleBond(species, finder.qc, finder.par,
                                       list(species.reac_inst[index]), name)
            targets(obj)  # Validate the native motif before any submission.
            species.reac_obj.append(obj)
            species.reac_name.append(name)
            species.reac_type.append('R_Addition_MultipleBond')
            species.reac_inst.append(list(species.reac_inst[index]))
            species.reac_step.append(0)
            species.reac_ts_done.append(0)
            species.reac_scan_energy.append([])
            species.reac_ts_geom.append(0)
            species.reac_ts_freq.append(0)
