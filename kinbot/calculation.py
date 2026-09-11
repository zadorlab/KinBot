"""Load the coherent record selected by a completed optimization."""
import logging
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


def publish_optimization_result(optimization):
    """Publish the complete accepted calculation under its conventional name.

    Keep the original calculation as evidence. PES and legacy restart readers
    continue to use the latest result under the usual well/high/low job name.
    """
    logger = logging.getLogger('KinBot')
    point = optimization.species
    job = getattr(point, 'source_job', None)
    if (not job or optimization.shigh != 1 or optimization.shir != 1
            or getattr(optimization, 'defer_hir', False)):
        logger.debug('Not publishing %s: optimization is not final or has no selected source.',
                     optimization.name)
        return
    sources = list(optimization.qc.db.select(name=job))
    if not sources:
        logger.warning('Cannot publish accepted result for %s: source %s is missing.',
                       optimization.name, job)
        return
    source = sources[-1]
    # Do not substitute a later observation for the result accepted by Optimize.
    if (source.data.get('status') != 'normal'
            or source.data.get('energy') is None
            or not np.isfinite(point.energy) or not np.isfinite(point.zpe)
            or point.energy != source.data['energy'] * constants.EVtoHARTREE
            or point.zpe != source.data.get('zpe')
            or not np.array_equal(point.geom, source.positions)
            or not np.array_equal(point.freq, source.data.get('frequencies'))):
        logger.warning('Cannot publish accepted result for %s: latest source %s '
                       'does not match the accepted calculation.', optimization.name, job)
        return
    if optimization.par['high_level']:
        target = optimization.log_name(1)
    elif optimization.par['conformer_search'] and not optimization.just_high:
        target = f'conf/{optimization.name}_low'
    else:
        target = optimization.log_name(0)
    optimization.qc.publish_result(source, target)
