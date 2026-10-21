"""One migration retry for the old unsymmetrized open-shell diatomic failure."""
import ast
import json
from pathlib import Path

from ase.io import read
from kinbot.execution import file_sha256


def evidence(qc, job, job_id):
    from kinbot.ase_modules.calculators.gaussian import (
        diatomic_recovery_parameters, density_scf_failure, qc_scf_failure)
    from kinbot.timeout_retry import _matching
    if getattr(qc, 'qc', None) != 'gauss':
        return None
    report, worker = Path(job+'.scf_recovery.json'), Path(job+'.py')
    if not report.is_file() or not worker.is_file():
        return None
    outcome = json.loads(report.read_text())
    if outcome.get('status') != 'failed' or outcome.get('recovery_kind') == 'diatomic_symmetry':
        return None
    log = Path(job+'.log')
    if not (density_scf_failure(log) or qc_scf_failure(log)):
        return None
    # Literal input only, never execute a saved Python worker to inspect it.
    kwargs = None
    for node in ast.parse(worker.read_text()).body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id=='kwargs' for t in node.targets):
            try:
                kwargs = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                return None
            break
    if not isinstance(kwargs, dict) or not kwargs.get('scf_recovery', True) or not kwargs.get('diatomic_scf_symmetry', True):
        return None
    atoms = read(job+'.com', format='gaussian-in')
    if diatomic_recovery_parameters(kwargs, atoms) is None:
        return None
    allocation = Path('.kinbot_timeout_allocations')/(job_id+'.json')
    if not allocation.is_file():
        return None  # Do not guess an expired allocation's owner or resources.
    native = json.loads(allocation.read_text())
    if not _matching(native, Path.cwd(), job_id) or not native.get('TimeLimit') or not native.get('Partition'):
        return None
    return dict(job_id=job_id, source='diatomic_symmetry',
                previous_walltime=native['TimeLimit'], partition=native['Partition'],
                report_sha256=file_sha256(report), native_sha256=file_sha256(log))
