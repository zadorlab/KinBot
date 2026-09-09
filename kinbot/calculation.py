"""Load the coherent record selected by a completed optimization."""
import numpy as np

from kinbot import constants


def load_calculation_record(species, qc, job):
    """Load geometry and molecular properties from one calculation record."""
    rows = list(qc.db.select(name=job))
    if not rows:
        raise ValueError(f'No database record found for selected job {job}.')
    row = rows[-1]
    data = row.data
    missing = [key for key in ('energy', 'zpe', 'frequencies')
               if data.get(key) is None]
    if missing:
        raise ValueError(
            f'Selected job {job} is missing properties: {missing}.')
    species.geom = np.asarray(row.positions, dtype=float)
    species.energy = float(data['energy']) * constants.EVtoHARTREE
    species.zpe = float(data['zpe'])
    species.freq = [float(value) for value in data['frequencies']]
    species.kinbot_freqs = list(species.freq)
    species.reduced_freqs = list(species.freq)
    hessian = data.get('hess')
    # Native Hessians can require formchk or a frequency calculation. Defer
    # that work until projection actually needs it; raw thermochemistry does not.
    if hessian is None:
        hessian = []
    species.hess = np.asarray(hessian).tolist() if len(hessian) else []
    species.source_job = job
    return job


def selected_calculation_job(optimization):
    """Return the one job selected by an Optimize instance."""
    if getattr(optimization, 'selected_job', None):
        return optimization.selected_job
    if optimization.par['high_level'] and optimization.shigh == 1:
        return optimization.log_name(1)
    conformers = getattr(optimization.species, 'confs', None)
    selected = getattr(conformers, 'selected_job', None)
    if optimization.par['conformer_search'] and selected:
        return selected
    return optimization.log_name(0)


def record_optimization_selection(optimization):
    """Publish an accepted source reference, without claiming IRC validation.

    Ordinary unchanged calculations keep their legacy lookup. Conformer or
    refinement selections have an explicit source_job and need provenance.
    """
    job = getattr(optimization.species, 'source_job', None)
    if not job:
        return
    if (optimization.shigh != 1 or optimization.shir != 1
            or getattr(optimization, 'defer_hir', False)):
        return
    db = optimization.qc.db
    sources = list(db.select(name=job))
    if not sources:
        return
    source = sources[-1]
    point = optimization.species
    # A later database row is not automatically a newly accepted geometry.
    # Publish only a record matching the properties accepted by Optimize.
    if (source.data.get('status') != 'normal'
            or source.data.get('energy') is None
            or not np.isfinite(point.energy) or not np.isfinite(point.zpe)
            or point.energy != source.data['energy'] * constants.EVtoHARTREE
            or point.zpe != source.data.get('zpe')
            or not np.array_equal(point.geom, source.positions)
            or not np.array_equal(point.freq, source.data.get('frequencies'))):
        return
    data = {
        'schema_version': 1, 'status': 'accepted',
        'source_job': job, 'source_row_id': source.id,
        'high_level': bool(optimization.par['high_level']),
        'conformer_search': bool(optimization.par['conformer_search']),
        'rotor_scan': bool(optimization.par['rotor_scan']),
    }
    name = f'optimization/{optimization.log_name(0)}'
    # Deduplicate polls of this optimization, not a new acceptance in a new
    # run that happens to select the same cached source again.
    if getattr(optimization, '_recorded_selection', None) != (name, data):
        db.write(source.toatoms(), name=name, data=data)
        optimization._recorded_selection = (name, data)


def optimization_selection_row(db, base_job, high_level, conformer_search,
                               rotor_scan=None, newer_than=0):
    """Resolve a current accepted selection, or use the legacy PES lookup.

    Never guess a selection from job names or choose a lower-energy row.
    Existing QC cache invalidation still governs method/basis changes.
    """
    selections = list(db.select(name=f'optimization/{base_job}'))
    if not selections:
        return None
    selection = selections[-1]
    data = selection.data
    if (selection.id <= newer_than
            or data.get('schema_version') != 1 or data.get('status') != 'accepted'
            or data.get('high_level') != bool(high_level)
            or data.get('conformer_search') != bool(conformer_search)
            or (rotor_scan is not None
                and data.get('rotor_scan') != bool(rotor_scan))):
        return None
    job = data.get('source_job')
    if not job:
        return None
    rows = list(db.select(name=job))
    if not rows or rows[-1].id != data.get('source_row_id'):
        return None
    row = rows[-1]
    values = [row.data.get('energy'), row.data.get('zpe')]
    if (row.data.get('status') != 'normal'
            or any(value is None or not np.isfinite(value) for value in values)):
        return None
    return row
