"""Acknowledge already reviewed, installed models before reaction replay.

Called only by the worker owning its input. This does not approve new orbital
choices, scan geometries, or rates. Invalid/unavailable models keep the hold.
"""
import json
import logging
from pathlib import Path
from types import SimpleNamespace

from kinbot.execution import atomic_json


def acknowledge_installed(parameters, species, directory):
    from kinbot.species_routing import routing_name
    from kinbot.vrc_multistate import channel_for_reaction, state_parameters
    from kinbot.vrc_correction import (correction_model, resolved_state_model,
        sampling_settings, spin_splitting_settings, trusted_settings, validate_vrc_levels)
    from kinbot.vrc_reference_review import validate, digest
    from kinbot.electronic_states import state_average_roots
    directory = Path(directory)
    paths = sorted((directory/'vrctst/orbital_review').glob('*/request.json'))
    if not paths:
        return []
    route = routing_name(species)
    reports = []
    coordinate_requests = {}
    for path in paths:
        try:
            saved = json.loads(path.read_text())
            bond = saved.get('bond')
            if (saved.get('parent') == route
                    and saved.get('status') in ('needs_review', 'model_accepted')
                    and saved.get('reaction') == path.parent.name
                    and isinstance(bond, list) and len(bond) == 2
                    and all(type(i) is int for i in bond)):
                coordinate_requests.setdefault(tuple(sorted(bond)), set()).add(saved['reaction'])
        except (OSError, ValueError, AttributeError):
            pass  # The individual request below retains its validation error.
    for path in paths:
        try:
            record = json.loads(path.read_text())
            if not isinstance(record, dict):
                raise ValueError('Invalid saved orbital review request.')
            if record.get('status') != 'needs_review':
                continue
            if record.get('review_kind') == 'fragment_identity':
                raise ValueError('Owning VRC driver must revalidate the fragment identity; '
                                 'an installed orbital model cannot clear this hold.')
            if (record.get('parent') != route or record.get('reaction') != path.parent.name
                    or record.get('symbols') != list(map(str, species.atom))
                    or record.get('charge') != int(species.charge)
                    or record.get('multiplicity') != int(species.mult)):
                raise ValueError('Review request does not match the current parent identity.')
            bond = record.get('bond')
            if (not isinstance(bond, list) or len(bond) != 2 or bond[0] == bond[1]
                    or any(type(i) is not int or not 1 <= i <= len(species.atom) for i in bond)):
                raise ValueError('Invalid saved review coordinate.')
            products = record.get('products')
            if not isinstance(products, list) or len(products) != 2:
                raise ValueError('Review request lacks the separated product identities.')
            fragments = []
            for product in products:
                if (not isinstance(product, dict)
                        or type(product.get('multiplicity')) is not int or product['multiplicity'] < 1
                        or type(product.get('charge')) is not int or not product.get('smiles')):
                    raise ValueError('Invalid saved fragment identity.')
                fragments.append(SimpleNamespace(smiles=product['smiles'], charge=product['charge'],
                                                 mult=product['multiplicity']))
            selection = channel_for_reaction(parameters, SimpleNamespace(
                species=species, instance_name=record['reaction'],
                scan_coo=[i-1 for i in bond], products=fragments))
            if selection is None:
                raise ValueError('No reviewed surface channel installed.')
            channel, source = selection
            if (channel.get('reaction') is None
                    and len(coordinate_requests.get(tuple(sorted(bond)), ())) > 1):
                raise ValueError('Different reaction requests share this coordinate; '
                                 'an exact reaction-specific surface model is required.')
            for stage in channel['bottlenecks']:
                for state in stage['states']:
                    leaf = state_parameters(parameters, stage, state, [p.mult for p in fragments])
                    validate_vrc_levels(leaf)
                    validate(leaf, directory)  # New spin references need their native review too.
                    spin = (spin_splitting_settings(leaf) if correction_model(leaf) == 'spin_splitting'
                            else None)
                    multiplicities = ([spin['target_multiplicity'], spin['reference_multiplicity']]
                                      if spin else [int(leaf['mult'])])
                    roots = state_average_roots(fragments, multiplicities, leaf)
                    levels = [sampling_settings(leaf)]
                    levels += [spin['reference'], *spin['gap_levels']] if spin else [trusted_settings(leaf)]
                    bases = {leaf['vrc_tst_scan_basis']}
                    for level in levels:
                        bases.update(level['basis'] if isinstance(level['basis'], list) else [level['basis']])
                    for basis in bases:
                        resolved_state_model(leaf, species.atom, int(species.charge), spin,
                                             roots=roots, basis=basis)
            record.update(status='model_accepted', surface_model=source,
                          installed_channel_sha256=digest(channel),
                          acknowledgement='Installed reviewed channel validated before reaction replay; '
                                          'scan and rate validation remain separate.')
            atomic_json(path, record)
            reports.append(dict(reaction=record['reaction'], status='model_accepted'))
        except (OSError, ValueError, KeyError, TypeError) as error:
            reports.append(dict(reaction=path.parent.name, status='needs_review', error=str(error)))
            logging.getLogger('KinBot').warning('VRC orbital hold retained for %s: %s', path.parent.name, error)
    return reports
