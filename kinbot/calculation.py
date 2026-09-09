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
