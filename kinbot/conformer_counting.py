"""Member-specific RRHO counting; site degeneracy is deliberately separate."""
from dataclasses import replace
import math
import logging
from pathlib import Path
from kinbot import constants, symmetry
from kinbot.conformer_records import ConformerRecord, retain
from kinbot.molecular_symmetry import (equivalent_geometry,
                                      OPTICAL_RMSD_TOLERANCE)
from kinbot.optical import evaluate_optical, compare_rigid, unresolved_optical_default
from kinbot.stereo_identity import canonical_identity, optical_scope
from kinbot.reaction_path import path_geometry_allowed
from kinbot.optical_harmonic import (apply_midpoint_heuristic, conformer_midpoint_diagnostic,
                                    conformer_pair_midpoint)

logger = logging.getLogger('KinBot')


class CountingError(ValueError):
    """The observations do not support a definitive automatic RRHO sum."""


def preserve_counting_error(species, error):
    from kinbot.stereo_routing import preserve_observations
    inventory = getattr(species, 'conformer_inventory', ())
    path = getattr(species, 'conformer_counting_error_evidence', None)
    if not (path and Path(path).exists()
            and getattr(species, 'conformer_counting_error', None) == str(error)
            and getattr(species, '_counting_error_inventory', None) == inventory):
        path = preserve_observations(str(error), [species])
    species.conformer_counting_error = str(error)
    species.conformer_counting_error_evidence = path
    species._counting_error_inventory = inventory
    error.evidence_path = path


def evaluate_members(species, records, population='specified', tolerance=.05, *, strict=True,
                     optical_tolerance=OPTICAL_RMSD_TOLERANCE):
    """Count rigid members and check coverage of analytically known mirrors.

    A self-mirror test establishes the possible pair size. Comparing members
    only checks whether that known mirror is explicitly represented; it never
    estimates missing reaction paths. Mirror partners share a pruning group.
    ``tolerance`` remains the stricter proper-duplicate threshold; optical
    matching has its own tolerance for numerical distortions of a mirror.
    """
    scope = optical_scope(species, population)
    species.optical_population = population
    species.optical_counting_scope = scope
    legacy_scope = scope['identity'].get('status') != 'assigned'
    if legacy_scope:
        from kinbot.stereo_identity import legacy_stereo_warning
        legacy_stereo_warning(species, scope['identity'].get('reason', 'unknown reference'))
    result = list(records)
    for record in result:
        if record.status == 'valid' and (record.zero_energy_hartree is None
                                        or not math.isfinite(record.zero_energy_hartree)):
            raise CountingError(f'MC member {record.member_id} lacks a finite E + ZPE.')
    unique = []
    order = sorted(range(len(result)), key=lambda i: (
        result[i].zero_energy_hartree if result[i].status == 'valid' else float('inf'),
        result[i].member_id))
    for offset in order:
        record = result[offset]
        if record.status != 'valid':
            continue
        if (getattr(species, 'wellorts', 0)
                and not path_geometry_allowed(species, record.geometry, population)):
            result[offset] = replace(record, exclusion_reason='different stereochemical pathway')
            continue
        identity = canonical_identity(species, record.geometry)
        if (not legacy_scope and not getattr(species, 'wellorts', 0)
                and identity.get('status') != 'assigned'):
            logger.warning('Excluding conformer %s: its configuration cannot be checked against '
                           'the specified stereoisomer (%s).', record.member_id,
                           identity.get('reason', 'unknown observation'))
            result[offset] = replace(record, exclusion_reason='unassigned configured stereoisomer')
            continue
        reference = scope['identity']
        allowed = {reference.get('id')}
        if population == 'racemic':
            allowed.add(reference.get('mirror_id'))
        if (reference.get('status') == identity.get('status') == 'assigned'
                and identity['id'] not in allowed):
            result[offset] = replace(record, exclusion_reason='different configured stereoisomer')
            continue
        numbers = symmetry.conformer_symmetry_numbers(species, record.geometry)
        optical = evaluate_optical(species, geometry=record.geometry, population=population,
                                   tolerance=optical_tolerance)
        previous_warnings = (record.optical_evidence or {}).get('warnings', ())
        if previous_warnings:
            optical['warnings'] = list(dict.fromkeys([*previous_warnings, *optical.get('warnings', ())]))
        optical = unresolved_optical_default(optical)
        if optical.get('fallback'):
            optical = apply_midpoint_heuristic(optical, conformer_midpoint_diagnostic(species, record))
        if strict and optical['remaining_multiplier'] is None:
            raise CountingError(f'Conformer {record.member_id}: {optical["reason"]}')
        record = replace(record, sigma_ext=numbers['sigma_ext'],
                         mirror_states=optical['total_optical_states'], optical_evidence=optical,
                         stereo_identity=identity.get('id'),
                         optical_population=population)
        result[offset] = record
        for previous in unique:
            other = result[previous]
            difference = abs(record.zero_energy_hartree-other.zero_energy_hartree) * constants.AUtoKCAL
            if equivalent_geometry(species, record.geometry, other.geometry, tolerance):
                if difference >= .2:
                    message = (f'Geometrically duplicate conformers {record.member_id} and '
                               f'{other.member_id} disagree by {difference:g} kcal/mol; '
                               f'retaining the complete lower-energy calculation {other.member_id}.')
                    logger.warning(message)
                    evidence = dict(other.optical_evidence)
                    evidence['warnings'] = list(dict.fromkeys([*evidence.get('warnings', ()), message]))
                    result[previous] = replace(other, optical_evidence=evidence)
                result[offset] = replace(record, duplicate_of=other.index,
                                         exclusion_reason='proper spatial duplicate')
                break
        else:
            unique.append(offset)
    groups = []
    remaining = set(unique)
    for offset in unique:
        if offset not in remaining:
            continue
        record = result[offset]
        partners = []
        if (scope['mirror_allowed'] and (record.mirror_states == 2
                or record.optical_evidence.get('remaining_multiplier') == 2.)):
            for other in unique:
                if other == offset or other not in remaining:
                    continue
                match = compare_rigid(species, record.geometry, result[other].geometry,
                                      optical_tolerance)['status']
                if match == 'undetermined':
                    pair = conformer_pair_midpoint(species, record, result[other])
                    for i, partner in ((offset, other), (other, offset)):
                        evidence = dict(result[i].optical_evidence)
                        comparisons = dict(evidence.get('explicit_mirror_comparisons', {}))
                        comparisons[result[partner].member_id] = pair
                        evidence['explicit_mirror_comparisons'] = comparisons
                        result[i] = replace(result[i], optical_evidence=evidence)
                    match = pair['status']
                if match == 'undetermined':
                    for i in (offset, other):
                        evidence = dict(result[i].optical_evidence,
                            status='unresolved', remaining_multiplier=None,
                            reason='Uncertain explicit mirror coverage.')
                        result[i] = replace(result[i], optical_evidence=unresolved_optical_default(evidence))
                if match == 'match':
                    partners.append(other)
        group = [offset] + partners
        warnings = []
        if len(partners) > 1:
            warnings.append('Several distinct conformers match the same mirror: '
                            + ', '.join(result[i].member_id for i in group)
                            + '. Mirror grouping is uncertain; retaining each calculation '
                              'once with optical factor 1, without an extra mirror multiplier.')
        for other in partners:
            difference = abs(record.zero_energy_hartree-result[other].zero_energy_hartree) * constants.AUtoKCAL
            if difference >= 1.:
                warnings.append(f'Explicit mirror conformers {record.member_id} and '
                                f'{result[other].member_id} disagree by {difference:g} kcal/mol; '
                                'retaining their own energies and properties with optical factor 1 each.')
        for message in warnings:
            logger.warning(message)
        for member in group:
            item = result[member]
            evidence = dict(item.optical_evidence)
            if warnings:
                evidence['warnings'] = list(dict.fromkeys([*evidence.get('warnings', ()), *warnings]))
            weight = 2. if item.mirror_states == 2 and scope['mirror_allowed'] and not partners else 1.
            if partners:
                # Explicit coverage wins over a single-conformer approximation.
                weight = 1.
                evidence.pop('fallback', None)
                if evidence.get('heuristic') == 'harmonic_midpoint':
                    evidence['single_conformer_heuristic'] = evidence.pop('heuristic')
                    evidence['harmonic_midpoint'] = dict(evidence['harmonic_midpoint'],
                        used_for_optical_counting=False,
                        superseded_by='explicit mirror coverage')
                approximate_pair = any(evidence.get('explicit_mirror_comparisons', {}).get(
                    result[i].member_id, {}).get('status') == 'match' for i in group if i != member)
                evidence['status'] = 'assumed' if approximate_pair or warnings else 'resolved'
            elif item.optical_evidence['status'] != 'resolved':
                weight = item.optical_evidence['remaining_multiplier']
            coverage = ('uncertain explicit mirrors; optical weight one' if len(partners) > 1 else
                        'explicit mirror' if partners else
                        'unresolved symmetry; optical weight one' if item.optical_evidence.get('fallback') else
                        'harmonic midpoint approximation' if evidence.get('heuristic') else
                        'undetermined; not used for statistical pruning' if weight is None else
                        'analytic missing mirror' if weight == 2 else
                        'self-mirror' if item.mirror_states == 1 else
                        'mirror outside specified population' if scope['identity'].get('status') == 'assigned'
                        else 'unknown configuration scope; no added mirror')
            result[member] = replace(item, remaining_optical_weight=weight,
                                     optical_evidence=dict(evidence,
                                         remaining_multiplier=weight,
                                         representation='MC-RRHO',
                                         reason=(evidence['reason']
                                                 if not partners and (evidence.get('fallback')
                                                     or evidence.get('heuristic')) else coverage),
                                         states_covered_by_mc_tst=(None if len(partners) > 1 else len(group)),
                                         explicit_mirror_ids=([result[i].member_id for i in group if i != member]
                                                              if len(partners) <= 1 else [])),
                                     mirror_coverage=coverage,
                                     mirror_partner_ids=(tuple(result[i].member_id
                                                              for i in group if i != member)
                                                         if len(partners) <= 1 else ()))
        groups.append(group)
        remaining.difference_update(group)
    return tuple(result), groups


def writer_members(species, population='specified', *, preserve_errors=True):
    """Refresh accepted final arrays and return records keyed by array offset.

    This also supports old restarts with no inventory. Electronic energy is
    not reconstructed from an ambiguous legacy conformer_energy array.
    """
    records = []
    offsets = {}
    # A previous writer can exclude a duplicate without removing its legacy
    # array entry. Its source still belongs to the complete inventory.
    existing = {record.index: record for record in getattr(species, 'conformer_inventory', ())}
    existing.update(getattr(species, 'conformer_records', {}))
    for offset, index in enumerate(getattr(species, 'conformer_index', [])):
        if index < 0:
            continue
        previous = existing.get(index)
        record = ConformerRecord(previous.member_id if previous else f'{species.name}:conformer:{index}', index,
                                 previous.source_job if previous else None, 'valid',
                                 geometry=tuple(tuple(map(float, atom)) for atom in species.conformer_geom[offset]),
                                 zero_energy_hartree=float(species.conformer_zeroenergy[offset]),
                                 frequencies_cm1=tuple(map(float, species.conformer_freq[offset])))
        if previous is not None:
            record = replace(record, electronic_energy_hartree=previous.electronic_energy_hartree,
                             zpe_hartree=previous.zpe_hartree,
                             hessian=previous.hessian, hessian_reference=previous.hessian_reference,
                             optical_evidence=({'warnings': previous.optical_evidence['warnings']}
                                               if (previous.optical_evidence or {}).get('warnings') else None),
                             attempted_source_job=previous.attempted_source_job)
        offsets[len(records)] = offset
        records.append(record)
    # Keep failed/filtered observations from the original calculation inventory.
    refreshed = {record.index: record for record in records}
    original = getattr(species, 'conformer_inventory', ())
    inventory = tuple(refreshed.pop(record.index, record) for record in original)
    inventory += tuple(refreshed.values())
    retain(species, inventory, [])
    try:
        records, groups = evaluate_members(species, records, population)
    except CountingError as error:
        if preserve_errors:
            preserve_counting_error(species, error)
        raise
    accepted = [member for group in groups for member in group]
    refreshed = {record.index: record for record in records}
    inventory = tuple(refreshed.get(record.index, record) for record in inventory)
    species.conformer_counting_error = None
    retain(species, inventory, [records[member].index for member in accepted])
    return {offsets[member]: records[member] for member in accepted}


def representative_record(species, population='specified', *, preserve_errors=True):
    """Count an accepted MC fallback without replacing its failed inventory."""
    record = ConformerRecord(f'{species.name}:selected-parent', 0,
        getattr(species, 'source_job', None), 'valid',
        geometry=tuple(tuple(map(float, xyz)) for xyz in species.geom),
        electronic_energy_hartree=float(species.energy), zpe_hartree=float(species.zpe),
        zero_energy_hartree=float(species.energy + species.zpe),
        hessian=tuple(tuple(map(float, row)) for row in getattr(species, 'hess', ())) or None,
        hessian_reference=getattr(species, 'optical_hessian_reference', None),
        frequencies_cm1=tuple(map(float, getattr(species, 'freq', species.reduced_freqs))))
    try:
        records, groups = evaluate_members(species, [record], population)
        if not groups:
            raise CountingError('Selected MC fallback is outside the configured population.')
    except CountingError as error:
        if preserve_errors:
            preserve_counting_error(species, error)
        raise
    species.rrho_representative_counting = records[0].as_dict()
    return records[0]
