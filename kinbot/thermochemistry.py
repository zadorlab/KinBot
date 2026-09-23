"""Backend-neutral thermochemical evidence, not a rate-model specification.

These records preserve observations and unknowns. They neither infer missing
reaction paths nor assert that a one-dimensional rotor represents all states.
"""
import numpy as np

from kinbot import constants, frequencies


def json_record(record):
    """Return native JSON values and report nonfinite numerical observations."""
    issues = []

    def native(value, path):
        if isinstance(value, np.ndarray):
            value = value.tolist()
        elif isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, dict):
            return {str(key): native(item, f'{path}.{key}') for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [native(item, f'{path}[{index}]') for index, item in enumerate(value)]
        if isinstance(value, float) and not np.isfinite(value):
            issues.append({'path': path, 'reason': 'nonfinite numerical observation'})
            return None
        if value is None or isinstance(value, (str, int, float, bool)):
            return value
        raise TypeError(f'Unsupported record value at {path}: {type(value).__name__}')

    result = native(record, '$')
    if issues:
        result.setdefault('data_issues', []).extend(issues)
    return result


def _number(value):
    return float(value) if value is not None and np.isfinite(value) else None


def hir_evidence(species):
    """Serialize completed/incomplete HIR objects without reading job files."""
    hir = getattr(species, 'hir', None)
    if hir is None:
        return {'status': 'not_recorded', 'rotors': []}
    records = []
    names = {-1: 'pending', 0: 'successful', 1: 'failed', 2: 'skipped'}
    for index, rotor in enumerate(species.dihed):
        statuses = hir.hir_status[index] if index < len(hir.hir_status) else []
        energies = hir.hir_energies[index] if index < len(hir.hir_energies) else []
        raw_arrays = getattr(hir, 'hir_raw_energies', [])
        raw = raw_arrays[index] if index < len(raw_arrays) else None
        raw = energies if raw is None else raw
        reference = (_number(raw[0]) if len(raw) and len(statuses)
                     and statuses[0] == 0 else None)
        coefficients = (hir.hir_fourier[index]
                        if index < len(hir.hir_fourier) else None)
        geometries = hir.hir_geoms[index] if index < len(hir.hir_geoms) else []
        points = []
        jobs = hir.scan_jobs[index] if index < len(hir.scan_jobs) else []
        observations = getattr(hir, 'point_observations', [])
        observations = observations[index] if index < len(observations) else []
        prior = getattr(hir, 'statuses_before_demotion', None)
        prior = prior[index] if prior is not None and index < len(prior) else statuses
        for point in range(hir.nrotation):
            status = statuses[point] if point < len(statuses) else -1
            energy = (_number(raw[point]) if status == 0 and point < len(raw) else None)
            geom = np.asarray(geometries[point]) if point < len(geometries) else np.array([])
            saved_observation = (observations[point] if point < len(observations) else None)
            measured_geom = (saved_observation.get('geometry_angstrom') if saved_observation is not None
                             else geom.tolist() if geom.shape == (species.natom, 3)
                             and np.all(np.isfinite(geom)) and (status == 0 or np.any(geom)) else None)
            angle = point * 2 * np.pi / hir.nrotation
            fit = (_number(hir.get_fit_value(angle, rotor=index))
                   if coefficients is not None else None)
            points.append({
                'index': point, 'angle_offset_degrees': point * 360. / hir.nrotation,
                'source_job': jobs[point] if point < len(jobs) else None,
                'status': names.get(status, 'unknown'), 'status_code': int(status),
                'observation': (saved_observation if saved_observation is not None else {
                                    'source_job': jobs[point] if point < len(jobs) else None,
                                    'qc_geometry_status': 'not_recorded',
                                    'electronic_energy_hartree': (
                                        _number(raw[point]) if point < len(raw)
                                        and point < len(prior) and prior[point] == 0 else None),
                                    'geometry_angstrom': measured_geom,
                                    'provenance': 'available legacy arrays; QC rejection details not recorded',
                                }),
                'data_issue': ('successful point lacks a finite measured energy'
                               if status == 0 and energy is None else None),
                'electronic_energy_hartree': energy,
                'relative_energy_kcal_mol': (
                    (energy - reference) * constants.AUtoKCAL
                    if energy is not None and reference is not None else None),
                'fitted_relative_energy_kcal_mol': fit,
                'geometry_angstrom': measured_geom,
            })
        first, second = frequencies.partition(species, rotor, species.natom)
        matrix = np.asarray(getattr(species, 'sigma_int', []))
        sigma = (int(species.sigma_int[rotor[1]][rotor[2]])
                 if matrix.shape == (species.natom, species.natom) else None)
        diagnostics = getattr(hir, 'hir_fit_diagnostics', [])
        records.append({
            'index': index, 'atom_index_base': 0,
            'dihedral': list(map(int, rotor)), 'axis': list(map(int, rotor[1:3])),
            'partition': [list(map(int, first)), list(map(int, second))],
            'sigma_int': sigma, 'sigma_int_source': 'KinBot graph heuristic',
            'potential_periodicity_verified': None,
            'scan_domain_degrees': [0., 360.], 'endpoint_included': False,
            'represented_domain_degrees': [0., 360. / sigma] if sigma and sigma > 0 else None,
            'reference_electronic_energy_hartree': reference,
            'usable': hir.is_valid_rotor(index),
            'exclusion_reason': hir.invalid_rotor_reason(index),
            'points': points,
            'fourier': {'energy_unit': 'kcal/mol', 'angle_unit': 'radian',
                        'basis': 'sum a[k]*(1-cos(k*angle)) + b[k]*sin(k*angle), k=1..N',
                        'coefficient_order': 'a[1..N], b[1..N]',
                        'coefficients': coefficients,
                        'diagnostics': diagnostics[index] if index < len(diagnostics) else None},
            'mirror_coverage': {'status': 'unknown', 'covered': None},
        })
    return json_record({'status': 'recorded',
            'backend': (getattr(hir, 'scan_reference', None) or {}).get('backend'),
            'backend_provenance': 'configured at scan generation; unknown if not recorded',
            'rigid_scan': hir.rigid_scan,
            'scan_reference': getattr(hir, 'scan_reference', None),
            'demotion_reason': getattr(hir, 'demotion_reason', None),
            'energy_kind': 'electronic; no ZPE added at scan points',
            'provenance': 'in-memory HIR scan results; failed points retain their status',
            'rotors': records})


def thermochemistry_evidence(species):
    """Explicit units and representation semantics for a stationary point."""
    from kinbot.counting_contract import counting_view, optical_counting
    species = counting_view(species)
    hir = hir_evidence(species)
    optical = optical_counting(species, hir)
    return json_record({
        'schema_version': 2,
        'rate_model_ready': False,
        'source_job': getattr(species, 'source_job', None),
        'source_row_id': getattr(species, 'source_row_id', None),
        'calculation_provenance': getattr(species, 'calculation_provenance', None),
        'atoms': list(map(str, species.atom)),
        'geometry_angstrom': np.asarray(species.geom).tolist(),
        'charge': getattr(species, 'charge', None),
        'multiplicity': getattr(species, 'mult', None),
        'electronic_energy_hartree': species.energy,
        'zpe_hartree': species.zpe,
        'sigma_ext': getattr(species, 'sigma_ext', None),
        'sigma_ext_source': 'KinBot graph symmetry rules',
        'legacy_nopt': getattr(species, 'nopt', None),
        'zpe_semantics': 'harmonic ZPE of the selected calculation; not rotor corrected',
        'raw_harmonic_frequencies_cm-1': [float(f) for f in getattr(species, 'freq', [])],
        'thermochemical_frequencies_cm-1': [float(f) for f in getattr(species, 'reduced_freqs', [])],
        'frequency_projection': getattr(species, 'rotor_projection', None),
        'conformer_representation': getattr(species, 'conformer_representation', None),
        'conformer_counting_error': getattr(species, 'conformer_counting_error', None),
        'conformer_counting_error_evidence': getattr(species, 'conformer_counting_error_evidence', None),
        'unassociated_conformer_arrays': getattr(species, 'unassociated_conformer_arrays', None),
        'rrho_representative_counting': getattr(species, 'rrho_representative_counting', None),
        'mess_tunneling_reference': getattr(species, 'mess_tunneling_reference', None),
        'conformer_inventory': [record.as_dict() for record in
                                getattr(species, 'conformer_inventory', ())],
        'optical_population_scope': optical.get('population_scope'),
        'hir': hir,
        'optical_counting': optical,
    })
