"""Attach completed ANL/CBH records to a reconstructed KinBot network."""

from __future__ import annotations

from dataclasses import fields
from hashlib import sha256
import json
import math
from pathlib import Path

from kinbot.anl.atct import ATcTRecord, _canonical_smiles
from kinbot.anl.cbh import FormationEnthalpy, zero_k_from_composite_record
from kinbot.anl.workflow import attach_task_vpt2_frequencies
from kinbot.energy import attach_formation_enthalpy
from kinbot.species_routing import routing_name


def _digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _path(base, value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{label} must be a nonempty path.')
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def _network_species(root):
    stable = {routing_name(root): root}
    transition_states = {}
    for index, reaction in enumerate(root.reac_obj):
        if root.reac_ts_done[index] != -1:
            continue
        if reaction.do_vdW:
            point = reaction.irc_prod_opt.species
            stable.setdefault(routing_name(point), point)
        for optimization in reaction.prod_opt:
            point = optimization.species
            key = routing_name(point)
            previous = stable.setdefault(key, point)
            if ((previous.charge, previous.mult, previous.smiles) !=
                    (point.charge, point.mult, point.smiles)):
                raise ValueError(f'{key}: conflicting stable species in network.')
        channel = ('hom_sci' in reaction.instance_name and not reaction.do_vdW
                   and len(reaction.prod_opt) == 2)
        if not channel:
            transition_states[reaction.instance_name] = reaction.ts
    return stable, transition_states


def _formation(path, final):
    try:
        record = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f'{path}: invalid CBH result.') from error
    if (not isinstance(record, dict) or record.get('schema') != 1
            or record.get('status') != 'complete'
            or not isinstance(record.get('formation'), dict)):
        raise ValueError(f'{path}: CBH result is incomplete.')
    raw = record['formation']
    required = {
        'target_smiles', 'rung', 'method', 'reaction_energy_0k_kj_mol',
        'formation_0k_kj_mol', 'references', 'atct_version',
        'atct_source_sha256', 'energy_sources'}
    if not required <= set(raw):
        raise ValueError(f'{path}: CBH formation fields are incomplete.')
    reference_fields = {field.name for field in fields(ATcTRecord)}
    references = {}
    if not isinstance(raw['references'], dict):
        raise ValueError(f'{path}: CBH references must be an object.')
    for smiles, value in raw['references'].items():
        if not isinstance(value, dict) or set(value) != reference_fields:
            raise ValueError(f'{path}: invalid ATcT record for {smiles!r}.')
        references[smiles] = ATcTRecord(**value)
    formation = FormationEnthalpy(
        target_smiles=raw['target_smiles'], rung=raw['rung'],
        method=raw['method'],
        reaction_energy_0k_kj_mol=raw['reaction_energy_0k_kj_mol'],
        formation_0k_kj_mol=raw['formation_0k_kj_mol'],
        references=references, atct_version=raw['atct_version'],
        atct_source_sha256=raw['atct_source_sha256'],
        energy_sources=raw['energy_sources'])
    if (isinstance(formation.rung, bool)
            or formation.rung not in range(4)
            or not all(math.isfinite(value) for value in (
                formation.reaction_energy_0k_kj_mol,
                formation.formation_0k_kj_mol))
            or _canonical_smiles(formation.target_smiles) !=
               _canonical_smiles(final.smiles)):
        raise ValueError(f'{path}: invalid CBH formation identity or values.')
    return formation


def _entry_map(value, label):
    if value is None:
        return {}
    if not isinstance(value, dict) or any(
            not isinstance(key, str) or not isinstance(entry, dict)
            for key, entry in value.items()):
        raise ValueError(f'{label} must be an object of keyed records.')
    return value


def _attach_vpt2(point, value, base, label):
    """Attach one explicitly accepted, mode-resolved VPT2 result."""
    if value is None:
        return None
    if not isinstance(value, dict) or set(value) - {
            'run', 'task_id', 'review', 'mode_map',
            'match_tolerance_cm_inverse'}:
        raise ValueError(f'{label}: invalid VPT2 handoff options.')
    options = {}
    if 'mode_map' in value:
        options['mode_map'] = value['mode_map']
    if 'match_tolerance_cm_inverse' in value:
        options['match_tolerance_cm_inverse'] = (
            value['match_tolerance_cm_inverse'])
    run_path = _path(base, value.get('run'), f'{label}.vpt2.run')
    review = (_path(base, value['review'], f'{label}.vpt2.review')
              if value.get('review') is not None else None)
    task_id = value.get('task_id', 'gaussian_vpt2')
    attach_task_vpt2_frequencies(
        point, run_path, task_id, review_file=review, **options)
    return dict(point.anl_thermochemistry_frequency_source)


def apply_handoff(root, manifest_file, *, report_file='me/anl_handoff.json'):
    """Apply hash-verifiable composite, CBH, and VPT2 records.

    The live network is reconstructed through KinBot's normal database path.
    This function only decorates exact routing keys and reaction names; it
    never creates species or silently matches by formula.
    """
    manifest_path = Path(manifest_file).expanduser().resolve()
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f'{manifest_path}: invalid ANL handoff manifest.') from error
    if (not isinstance(manifest, dict) or manifest.get('schema') != 1
            or manifest.get('mode') not in ('anl', 'cbh-anl')):
        raise ValueError('ANL handoff requires schema 1 and mode anl or cbh-anl.')
    strict = manifest.get('strict', True)
    if not isinstance(strict, bool):
        raise ValueError('ANL handoff strict must be boolean.')
    species_entries = _entry_map(manifest.get('species'), 'species')
    ts_entries = _entry_map(
        manifest.get('transition_states'), 'transition_states')
    stable, transition_states = _network_species(root)
    unknown = sorted(set(species_entries) - set(stable))
    unknown_ts = sorted(set(ts_entries) - set(transition_states))
    if unknown or unknown_ts:
        raise ValueError('ANL handoff contains unknown network identities: '
                         f'species={unknown}, transition_states={unknown_ts}.')
    if strict:
        missing = sorted(set(stable) - set(species_entries))
        missing_ts = sorted(set(transition_states) - set(ts_entries))
        if missing or missing_ts:
            raise ValueError('ANL handoff does not cover the complete network: '
                             f'species={missing}, transition_states={missing_ts}.')

    base = manifest_path.parent
    applied = {'species': {}, 'transition_states': {}}
    for key, entry in species_entries.items():
        point = stable[key]
        composite_path = _path(base, entry.get('composite'),
                               f'{key}.composite')
        final = zero_k_from_composite_record(point.smiles, composite_path)
        if (final.charge, final.multiplicity) != (point.charge, point.mult):
            raise ValueError(f'{key}: composite electronic state disagrees with KinBot.')
        point.final_zero_k_energy = final
        record = {
            'composite': str(composite_path),
            'composite_sha256': _digest(composite_path),
            'method': final.method,
            'zero_k_hartree': final.hartree,
        }
        if manifest['mode'] == 'cbh-anl':
            formation_path = _path(base, entry.get('formation'),
                                   f'{key}.formation')
            formation = _formation(formation_path, final)
            attach_formation_enthalpy(point, formation)
            record.update(
                formation=str(formation_path),
                formation_sha256=_digest(formation_path),
                formation_0k_kj_mol=formation.formation_0k_kj_mol)
        elif entry.get('formation') is not None:
            raise ValueError(f'{key}: formation is only valid in cbh-anl mode.')
        vpt2 = _attach_vpt2(point, entry.get('vpt2'), base, key)
        if vpt2 is not None:
            record['vpt2'] = vpt2
        applied['species'][key] = record

    for reaction_name, entry in ts_entries.items():
        point = transition_states[reaction_name]
        composite_path = _path(
            base, entry.get('composite'), f'{reaction_name}.composite')
        if not point.smiles:
            raise ValueError(
                f'{reaction_name}: transition state has no SMILES identity.')
        final = zero_k_from_composite_record(point.smiles, composite_path)
        if (final.charge, final.multiplicity) != (point.charge, point.mult):
            raise ValueError(
                f'{reaction_name}: composite state disagrees with KinBot TS.')
        point.final_zero_k_energy = final
        record = {
            'composite': str(composite_path),
            'composite_sha256': _digest(composite_path),
            'method': final.method, 'zero_k_hartree': final.hartree,
        }
        vpt2 = _attach_vpt2(
            point, entry.get('vpt2'), base, reaction_name)
        if vpt2 is not None:
            record['vpt2'] = vpt2
        applied['transition_states'][reaction_name] = record

    result = {
        'schema': 1, 'status': 'complete', 'mode': manifest['mode'],
        'strict': strict, 'manifest': str(manifest_path),
        'manifest_sha256': _digest(manifest_path), **applied,
    }
    destination = Path(report_file)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    return result
