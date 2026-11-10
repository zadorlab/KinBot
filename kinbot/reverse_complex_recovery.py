"""Owner-scoped numerical recovery for complexes attached to accepted reverse TSs."""
from pathlib import Path

from kinbot.calculation import selected_calculation_job


def drive(optimization, transition_state=None):
    """Poll/recover one minimum, never reinterpret a failed saddle search.

    Existing SCF/ODE recoveries remain in the QC driver. This context allows
    one Cartesian fallback after an internal-coordinate minimum exhausts its
    step budget. Normal minimum/Hessian and connectivity validation still run.
    """
    if (transition_state is not None and transition_state.shigh == 1
            and transition_state.species.wellorts
            and optimization.just_high and not optimization.species.wellorts
            and optimization.par.get('high_level')):
        job = optimization.log_name(1)
        ts_job = selected_calculation_job(transition_state)
        qc = optimization.qc
        backend = qc._backend_for(job) if hasattr(qc, '_backend_for') else qc
        ts_backend = qc._backend_for(ts_job) if hasattr(qc, '_backend_for') else qc
        rows = list(ts_backend.db.select(name=ts_job))
        if rows and rows[-1].data.get('status') == 'normal':
            context = dict(transition_state=ts_job, row_id=rows[-1].id,
                           directory=str(Path.cwd()), complex_job=job)
            if not hasattr(backend, '_reverse_complex_recovery'):
                backend._reverse_complex_recovery = {}
            backend._reverse_complex_recovery[job] = context
            if optimization.shigh == -999:
                status = qc.check_qc(job)
                if status in ('running', 'normal'):
                    # Adopt the replacement through the owning optimizer so
                    # its geometry, modes and published result stay coherent.
                    optimization.shigh = 0
    if optimization.shigh != -999:
        optimization.do_optimization()
