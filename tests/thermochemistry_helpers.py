"""Inspection helpers used by counting and HIR regressions."""
import numpy as np
from kinbot.thermochemistry import hir_evidence, json_record


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
