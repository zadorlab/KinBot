import ast
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from kinbot.qc import QuantumChemistry


def load_helper():
    path = (Path(__file__).parents[1] / 'examples/anl/tools/check_minimum_recovery.py')
    spec = importlib.util.spec_from_file_location('minimum_check_helper', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize('level', ['l1', 'l2'])
def test_check_uses_original_geometry_and_production_worker_without_touching_old_run(tmp_path, monkeypatch, level):
    helper = load_helper()
    old = tmp_path/'old'
    old.mkdir()
    settings = old/'input.json'
    settings.write_text(json.dumps({'theory_preset': 'b3lyp-b2plyp-anl', 'mult': 3,
        'barrier_threshold': 1000., 'slurm_exclude_nodes': ['node16']}))
    worker = old/'well_high.py'
    geom = [[0., 0., 0.], [.9, 0., 0.], [-.2, .8, 0.]]
    worker.write_text("raise RuntimeError('old worker must never execute')\n"
        f"mol = Atoms(symbols=['O', 'H', 'H'], positions={geom!r})\n"
        "kwargs = {'charge': 0, 'mult': 1}\n")
    before = {p: p.read_bytes() for p in old.iterdir()}
    events = []
    def submit(qc, job, nproc, **kwargs):
        assert qc.par['slurm_exclude_nodes'] == ['node16']
        script = Path(job+'.py')
        ast.parse(script.read_text())
        np.testing.assert_allclose(helper.starting_geometry(script)[1], geom)
        assert helper.worker_state(script) == (0, 1)
        assert "frequency_mode = 'native_hessian'" in script.read_text()
        assert 'gaussian_minimum_restart' in script.read_text()
        events.append((Path.cwd(), job, nproc))
        qc.job_ids[job] = 123
        return 1
    monkeypatch.setattr(QuantumChemistry, 'submit_qc', submit)
    destination = tmp_path/'check'
    cwd = Path.cwd()
    record = helper.submit_check(settings, worker, destination, level)
    assert Path.cwd() == cwd
    assert len(events) == 1 and events[0][0] == destination
    assert record['status'] == 'submitted' and record['job_ids']
    assert record['multiplicity'] == 1
    assert before == {p: p.read_bytes() for p in old.iterdir()}
    with pytest.raises(FileExistsError, match='duplicate'):
        helper.submit_check(settings, worker, destination)
    assert len(events) == 1


def test_worker_source_cannot_execute_expressions(tmp_path):
    worker = tmp_path/'untrusted.py'
    worker.write_text("mol = Atoms(symbols=['H'], positions=__import__('os').system('false'))")
    with pytest.raises(ValueError):
        load_helper().starting_geometry(worker)
