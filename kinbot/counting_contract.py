"""Counting for the represented stationary states, independent of a rate adapter."""
import copy
import numpy as np

from kinbot.calculation import array_fingerprint
from kinbot.conformer_counting import writer_members, representative_record, CountingError
from kinbot.molecular_symmetry import OPTICAL_RMSD_TOLERANCE


def counting_view(species):
    """Evaluate current MC arrays without changing selection or writing files."""
    view = copy.copy(species)
    if getattr(view, 'conformer_representation', None) == 'MC-RRHO':
        try:
            indices = getattr(view, 'conformer_index', [])
            lengths = [len(getattr(view, key, [])) for key in (
                'conformer_geom', 'conformer_zeroenergy', 'conformer_freq')]
            if any(not isinstance(i, (int, np.integer)) or isinstance(i, bool) for i in indices):
                raise ValueError('Conformer calculation indices must be integers.')
            active = [i for i in indices if i >= 0]
            if (len(set(active)) != len(active)
                    or (active and any(n != len(indices) for n in lengths))
                    or (not active and any(n not in (0, len(indices)) for n in lengths))):
                raise ValueError('Conformer arrays have ambiguous lengths or duplicate calculation indices.')
            members = writer_members(view, getattr(view, 'optical_population', 'specified'),
                                     preserve_errors=False)
            if not members:
                representative_record(view, getattr(view, 'optical_population', 'specified'),
                                      preserve_errors=False)
        except (CountingError, ValueError, TypeError, IndexError, AttributeError) as error:
            view.conformer_counting_error = str(error)
            if not isinstance(error, CountingError):
                view.unassociated_conformer_arrays = {
                    key: getattr(species, key, None) for key in (
                        'conformer_index', 'conformer_geom', 'conformer_zeroenergy', 'conformer_freq')}
    return view


def _reference_issue(species, hir, projection):
    """A mutable species pointer is not evidence of a scan's original input."""
    for label, reference in [('scan', hir.get('scan_reference')),
                             ('projection', projection.get('reference'))]:
        if reference is None:
            return f'{label} input provenance was not recorded'
        if reference.get('geometry_sha256') != array_fingerprint(species.geom):
            return f'{label} input differs from the selected geometry'
        if reference.get('atoms') != list(map(str, species.atom)):
            return f'{label} atom order differs from the selected calculation'
        for key in ('source_job', 'source_row_id'):
            selected = getattr(species, key, None)
            if selected is None or reference.get(key) is None:
                return f'{label} {key} is unknown'
            if reference[key] != selected:
                return f'{label} {key} belongs to another selected calculation'
        if reference.get('dihedrals') != [list(map(int, r)) for r in species.dihed]:
            return f'{label} rotor definitions differ from the current rotor axes'
    reference = projection['reference']
    if not reference.get('hessian_sha256'):
        return 'projection Hessian input was not recorded'
    hessian = getattr(species, 'hess', [])
    if len(hessian) and reference['hessian_sha256'] != array_fingerprint(hessian):
        return 'projection used a different Hessian from the selected record'
    return None


def _projection_decisions(species, hir, projection):
    """Check that the exported frequency set represents the declared rotors."""
    decisions = {}
    usable = {r['index'] for r in hir['rotors'] if r['usable']}
    if not projection:
        issue = 'rotor projection decisions were not recorded' if usable else None
        if len(species.freq) != len(species.reduced_freqs):
            issue = 'frequency mode counts differ without a recorded projection'
        return decisions, issue
    for entry in projection.get('rotors', []):
        index, projected = entry.get('rotor_index'), entry.get('projected')
        if (not isinstance(index, (int, np.integer)) or isinstance(index, bool)
                or index in decisions or not isinstance(projected, (bool, np.bool_))):
            return decisions, 'rotor projection decisions are duplicated or malformed'
        decisions[index] = bool(projected)
    if set(decisions) != set(range(len(species.dihed))):
        return decisions, 'individual rotor projection decisions are missing or inconsistent'
    rank = projection.get('internal_rank')
    if (not isinstance(rank, (int, np.integer)) or isinstance(rank, bool)
            or rank != sum(decisions.values())):
        return decisions, 'internal projection rank disagrees with the rotor decisions'
    if any(projected and index not in usable for index, projected in decisions.items()):
        return decisions, 'a projected rotor is now unavailable; harmonic frequencies must be restored'
    if len(species.freq) - len(species.reduced_freqs) != rank:
        return decisions, 'frequency mode counts disagree with the internal projection rank'
    return decisions, None


def optical_counting(species, hir, tolerance=OPTICAL_RMSD_TOLERANCE, energy_tolerance=1.):
    """Validate the calculation model, then use the shared optical evaluator."""
    from kinbot.optical import (evaluate_optical, unresolved_optical_default, METHOD,
                               SCAN_MIRROR_RMSD_TOLERANCE)
    data = dict(method=METHOD, status='unresolved', remaining_multiplier=None,
                population=getattr(species, 'optical_population', 'specified'),
                scope='selected configuration and its global mirror; not independent fragment racemates',
                total_optical_states=None, total_configurational_states=None,
                allowed_global_mirror_states=None, states_covered_by_hir=None,
                states_covered_by_mc_tst=None, site_degeneracy_included=False,
                members=[], local_rmsd_tolerance_angstrom=tolerance,
                measured_mirror_rmsd_tolerance_angstrom=min(tolerance, SCAN_MIRROR_RMSD_TOLERANCE),
                mirror_energy_tolerance_kcal_mol=energy_tolerance)
    projection = getattr(species, 'rotor_projection', None) or {}
    decisions, issue = _projection_decisions(species, hir, projection)
    active = [r for r in hir['rotors'] if r['usable'] and decisions.get(r['index']) is not False]
    if len({tuple(sorted(r['axis'])) for r in active}) != len(active):
        data['reason'] = 'multiple rotor records describe the same axis.'
        return data
    if getattr(species, 'conformer_representation', None) == 'MC-RRHO':
        from kinbot.stereo_identity import optical_scope
        data['representation'] = 'MC-RRHO'
        data['population_scope'] = optical_scope(species, data['population'])
        if active or issue or getattr(species, 'conformer_counting_error', None):
            data['reason'] = ('MC-RRHO and active HIR cannot describe the same modes'
                              if active else issue or species.conformer_counting_error)
            return data
        members = [r.as_dict() for r in getattr(species, 'conformer_records', {}).values()]
        if not members and getattr(species, 'rrho_representative_counting', None):
            members = [species.rrho_representative_counting]
        data['members'] = [{key: record.get(key) for key in (
            'member_id', 'sigma_ext', 'mirror_states', 'remaining_optical_weight',
            'mirror_partner_ids', 'mirror_coverage', 'optical_evidence')} for record in members]
        data.update(status='per_member' if members else 'unresolved',
                    reason='Apply each conformer optical weight once; no ensemble-wide multiplier.')
        return data
    data['representation'] = 'HIR' if active else 'RRHO'
    if issue or (active and (issue := _reference_issue(species, hir, projection))):
        data['reason'] = issue
        # The shape classification is still useful raw evidence when the
        # statistical model cannot be used (for example, stale projection).
        from kinbot.optical import rigid_mirror
        try:
            data['rigid_mirror'] = rigid_mirror(species, tolerance=tolerance)
            data['total_optical_states'] = data['rigid_mirror']['mirror_states']
        except (ValueError, TypeError, IndexError):
            pass
        data['reported_optical_states'] = getattr(species, 'optical_isomers', None)
        data['reported_count_disagrees'] = (data['total_optical_states'] is not None
            and data['reported_optical_states'] not in (None, -1, data['total_optical_states']))
        return data
    try:
        result = evaluate_optical(species, rotors=active, tolerance=tolerance,
                                  energy_tolerance=energy_tolerance)
    except (ValueError, TypeError, IndexError) as error:
        data['reason'] = str(error)
        return data
    data.update(result)
    if (data.get('population_scope', {}).get('identity', {}).get('status') == 'unavailable'
            and isinstance(getattr(species, 'nopt', None), (int, float, np.number))
            and np.isfinite(species.nopt) and species.nopt > 0):
        data.update(status='legacy_unverified', remaining_multiplier=float(species.nopt),
                    reason='RDKit unavailable; retained legacy optical convention.')
    data['reported_optical_states'] = getattr(species, 'optical_isomers', None)
    data['reported_count_disagrees'] = (data['total_optical_states'] is not None
        and data['reported_optical_states'] not in (None, -1, data['total_optical_states']))
    if active:
        data.update(hir_coverage_convention='full torsional circle divided by assigned internal symmetry',
                    hir_symmetry_periodicity_assumed=True)
    assumption = getattr(species, 'optical_factor_assumption', None)
    # Only geometric model ambiguity is overridable. Calculation provenance,
    # population restrictions and contradictory measured data are not.
    if (assumption and data['status'] == 'unresolved'
            and data.get('population_scope', {}).get('mirror_allowed')
            and data.get('reason') in (
                'The represented motions do not establish mirror coverage or retained handedness.',
                'Rigid mirror comparison is numerically undetermined.')):
        data.update(status='assumed', remaining_multiplier=assumption['multiplier'],
                    assumption=dict(assumption), reason=assumption['reason'])
    data = unresolved_optical_default(data)
    if data.get('fallback') == 'unresolved_symmetry':
        from kinbot.optical_harmonic import selected_midpoint_diagnostic, apply_midpoint_heuristic
        data = apply_midpoint_heuristic(data, selected_midpoint_diagnostic(species, active))
    return data


def site_counting(degeneracy, members, relations, stereopath_id, *, conflict=False):
    """Convert labelled replicas of ONE declared unlabelled contribution.

    This is not a test of physical path or conformer completeness. Rotational,
    internal-rotor and optical factors belong to the stationary-state records.
    """
    result = {'degeneracy': degeneracy, 'site_members': members, 'site_relations': relations,
              'stereopath_id': stereopath_id, 'physical_path_completeness': 'not_established',
              'unlabelled_model_conversion': {'status': 'unresolved',
                  'convention': 'canonical unlabelled stationary-state partition functions',
                  'additional_site_multiplier': None, 'labelled_replica_normalization': None}}
    conversion = result['unlabelled_model_conversion']
    if conflict:
        conversion['reason'] = 'endpoint, energy or site-class observations conflict'
    elif (not isinstance(stereopath_id, str) or not stereopath_id
          or not isinstance(degeneracy, (int, np.integer))
          or isinstance(degeneracy, bool) or degeneracy < 1
          or not isinstance(members, dict) or not members or not isinstance(relations, dict)
          or set(members) != set(relations)):
        conversion['reason'] = 'site-class evidence is missing or inconsistent'
    elif (any(not isinstance(atoms, (list, tuple, np.ndarray))
              or (isinstance(atoms, np.ndarray) and atoms.ndim != 1)
              or not len(atoms)
              or any(not isinstance(i, (int, np.integer)) or isinstance(i, (bool, np.bool_))
                     or i < 0 for i in atoms)
              or len(atoms) != len(set(atoms)) for atoms in members.values())
          or int(np.prod([len(atoms) for atoms in members.values()])) != degeneracy
          or any(not isinstance(r, str) or r not in {'homotopic', 'enantiotopic', 'diastereotopic', 'singleton'}
                 for r in relations.values())):
        conversion['reason'] = 'site-class members, relations and degeneracy disagree'
    else:
        conversion.update(status='defined_for_declared_unlabelled_model',
            additional_site_multiplier=1., labelled_replica_normalization=1. / degeneracy,
            reason='No additional raw site factor. The 1/d conversion only splits this same canonical contribution into labelled replicas; it does not equate independent TS basins.')
    return result
