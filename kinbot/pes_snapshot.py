"""Accepted L2 evidence published before downstream planning and CBH expansion.

This diagnostic is deliberately independent of successful CBH construction.
It does not select reactions, launch calculations, or change their acceptance.
"""
from datetime import datetime, timezone
from pathlib import Path

from kinbot.execution import atomic_json
from kinbot.species_routing import routing_name


def snapshot(species, parameters):
    from kinbot.anl.automatic import _l2_species_record, _resolved_l2_theory
    from kinbot.reaction_path import reaction_path_id
    theory = _resolved_l2_theory(parameters)
    network = dict(initial=routing_name(species), route_to_smiles={},
                   stable_species=[], transition_states={}, reactions=[])
    errors = []
    stable = {}
    def point(obj, **kwargs):
        return _l2_species_record(obj, theory, **kwargs)
    def add(obj):
        name = routing_name(obj)
        network['route_to_smiles'][name] = str(obj.smiles)
        try:
            stable[name] = point(obj)
        except (ValueError, TypeError, AttributeError) as error:
            errors.append(dict(name=name, error=str(error)))
        return name
    add(species)
    for i, reaction in enumerate(species.reac_obj):
        if species.reac_ts_done[i] != -1:
            continue
        name = str(reaction.instance_name)
        family = species.reac_type[i]
        row = dict(name=name, family=family, reactant=routing_name(species),
                   products=sorted(add(p) for p in reaction.products),
                   barrierless=family == 'hom_sci',
                   transition_state=None if family == 'hom_sci' else name)
        if not row['barrierless']:
            try:
                row['path_class'] = reaction_path_id(reaction)
                network['transition_states'][name] = point(
                    reaction.ts, identity_smiles=species.smiles, routing_identity=name)
            except (ValueError, TypeError, AttributeError) as error:
                errors.append(dict(name=name, error=str(error)))
        network['reactions'].append(row)
    network['stable_species'] = list(stable.values())
    return dict(schema=1, source='accepted_l2_before_downstream', network=network,
                captured_utc=datetime.now(timezone.utc).isoformat(), errors=errors)


def publish(generator):
    """One small atomic write when exploration finishes, including live revisits."""
    if generator is not None and getattr(generator, 'par', {}).get('composite_method'):
        atomic_json(Path('accepted_l2_network.json'),
                    snapshot(generator.species, generator.par))
