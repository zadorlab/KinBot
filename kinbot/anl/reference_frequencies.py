"""Recover harmonic reference spectra without changing saved QC or models."""
from copy import deepcopy
from pathlib import Path

import numpy as np
from ase.db import connect

from kinbot import constants
from kinbot.anl.ts_reuse import digest
from kinbot.execution import file_sha256


POLICY = 'selected-row-harmonic-spectrum-v1'


def recover(record):
    """Derive frequencies from the exact accepted row, keeping original hashes.

    Harmonic spectra, rather than rotor-reduced spectra, establish L2 reuse.
    For a rotor-free thermal model the same native spectrum also repairs old
    exports made before Gaussian Hessians were aligned to the selected frame.
    No geometry, energy, model file or dispatcher input is rewritten.
    """
    model = record['thermochemistry_model']
    source = model['source']
    job, row_id = source.get('source_job'), source.get('source_row_id')
    root = Path(record['thermochemistry_model_file']).parent
    database = root/'kinbot.db'
    if not job or row_id is None or not database.is_file():
        return record  # Older/external models need separate source evidence.
    with connect(str(database)) as db:
        row = db.get(id=int(row_id))
    data = row.data
    if (row.name != job or data.get('status') != 'normal'
            or list(row.symbols) != model['symbols']
            or not np.allclose(row.positions, model['positions'], atol=1e-8, rtol=0)
            or not np.isclose(data.get('energy', np.nan)*constants.EVtoHARTREE,
                              record['energy_hartree'], atol=1e-10, rtol=0)
            or not np.isclose(data.get('zpe', np.nan), record['harmonic_zpe_hartree'],
                              atol=1e-10, rtol=0)):
        raise ValueError(f'{root}: reference frequencies do not match the selected L2 row.')
    native = data.get('frequencies')
    # ASE database round-trips preserve the arrays written by QC templates,
    # including a length-zero array for an atom. JSON-derived rows use lists.
    if not isinstance(native, (list, tuple, np.ndarray)):
        raise ValueError(f'{root}: selected L2 row lacks finite harmonic frequencies.')
    try:
        native = np.asarray(native)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{root}: selected L2 row has malformed harmonic frequencies.') from exc
    if (native.ndim != 1 or native.dtype.kind not in 'iuf'
            or not np.isfinite(native).all()):
        raise ValueError(f'{root}: selected L2 row lacks finite harmonic frequencies.')
    native = native.astype(float).tolist()
    n = len(model['symbols'])
    counts = {0} if n == 1 else {1} if n == 2 else {3*n-6, 3*n-5}
    if len(native) not in counts or any(f < -50. for f in native):
        raise ValueError(f'{root}: selected L2 reference has invalid minimum frequencies.')
    thermal = [abs(f) for f in native]
    rotor_free = (not model.get('hindered_rotor_blocks', '').strip()
                  and model.get('hindered_rotor_source') in
                  ('none-no-rotatable-bonds', 'harmonic_fallback'))
    reexport = (bool(model.get('hindered_rotor_blocks', '').strip())
                and record.get('l2_theory', {}).get('calculator') == 'gaussian'
                and source.get('hessian_frame_policy') != 'selected-geometry-v1')
    changes = record['frequencies_cm_inverse'] != native or (
        rotor_free and model['frequencies_cm_inverse'] != thermal) or reexport
    if not changes:
        return record
    result = deepcopy(record)
    result['frequencies_cm_inverse'] = native
    if rotor_free:
        result['thermochemistry_model']['frequencies_cm_inverse'] = thermal
    # A legacy model with rotor blocks cannot be repaired by replacing its
    # reduced spectrum with the full harmonic list (that double counts modes).
    # Leave that model intact; only the harmonic comparison is corrected here.
    evidence = dict(policy=POLICY, previous_record_sha256=digest(record),
        source_job=job, source_row_id=row.id, harmonic_frequencies_cm_inverse=native,
        thermal_spectrum_recovered=rotor_free, model_sha256=file_sha256(
            record['thermochemistry_model_file']))
    evidence['thermal_spectrum_requires_reexport'] = reexport
    result['frequency_recovery'] = evidence
    if rotor_free:
        result['thermochemistry_model']['frequency_source'] = evidence
    return result


def verified_migration(old, new):
    """Only accept the exact read-time derivation from an unchanged source."""
    evidence = new.get('frequency_recovery', {})
    if (evidence.get('policy') != POLICY
            or evidence.get('previous_record_sha256') != digest(old)
            or old.get('provider') != 'automatic-kinbot-reference-l2'):
        return False
    return recover(old) == new
