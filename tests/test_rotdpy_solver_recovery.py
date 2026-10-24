from copy import deepcopy
import json
from pathlib import Path

import pytest
from ase import Atoms
from tests.test_rotdpy_profiles import sulfur_correction
from kinbot.pes import _rotdpy_validation_metadata, _rotdpy_sampling_calculator

NATIVE = ('Number of p-space configurations: 4\nCatastrophic failure in diagonalization(hsdel)\n'
          'The expansion set has become singular\n')


@pytest.fixture
def sample(tmp_path, monkeypatch):
    pytest.importorskip('rotd_py')
    from kinbot.rotdpy_molpro import MolproSample
    par,data,path = sulfur_correction(tmp_path)
    model = _rotdpy_sampling_calculator(par,data,_rotdpy_validation_metadata(par,data,path,False,4))
    symbols = [s for f in data['frags_atom'] for s in f]
    atoms = Atoms(symbols, positions=[[i,0,0] for i in range(len(symbols))])
    calc = dict(scratch=str(tmp_path),processors=1,mem=300,method=model['method'],basis=model['basis'],kinbot_vrc=model)
    monkeypatch.chdir(tmp_path)
    sample = MolproSample('sample',atoms,calc); sample.create_input()
    return sample


def test_exact_sample_two_retries_preserve_science_and_archive(sample, monkeypatch):
    calls = []; child = sample.basis_jobs[0]
    old = deepcopy(child.kinbot_retry)
    def native(self, job):
        calls.append(Path(job.name+'.inp').read_text())
        if len(calls) <= 2:
            Path(job.name+'.out').write_text(NATIVE)
            error = RuntimeError(f'Molpro sample {job.name} exited with status 1.')
            error.native_returncode = 1
            error.native_output_current = True
            raise error
        Path(job.name+'.out').write_text('normal native output')
    monkeypatch.setattr(type(sample),'_run_basis_job',native)
    sample._run_with_solver_recovery(child)
    assert len(calls) == 3
    assert 'pspace,0,8;' not in calls[0] and 'pspace,0,8;' in calls[1]
    for deck in calls:
        assert deck.split('geometry={')[1].split('}')[0] == calls[0].split('geometry={')[1].split('}')[0]
        assert 'shift=0.1,ipea=0.25' in deck
    assert child.kinbot_retry == old
    journal = next(Path('.kinbot_sample_attempts').rglob('recovery.json'))
    record = json.loads(journal.read_text())
    assert len(record['attempts']) == 2 and record['status'] == 'native_executed'
    assert all((Path(a['archive'])/(child.name+'.out')).read_text() == NATIVE for a in record['attempts'])


def test_exhaustion_persists_across_worker_restart(sample, monkeypatch):
    calls=[]
    def native(self, job):
        calls.append(job.name); Path(job.name+'.out').write_text(NATIVE)
        error=RuntimeError(f'Molpro sample {job.name} exited with status 1.');error.native_returncode=1
        error.native_output_current=True
        raise error
    monkeypatch.setattr(type(sample),'_run_basis_job',native)
    with pytest.raises(RuntimeError):sample._run_with_solver_recovery(sample.basis_jobs[0])
    assert len(calls) == 3
    sample.create_input()
    with pytest.raises(RuntimeError,match='exhausted'):sample._run_with_solver_recovery(sample.basis_jobs[0])
    assert len(calls) == 3


def test_stale_output_is_not_used_as_native_retry_evidence(sample, monkeypatch):
    child=sample.basis_jobs[0]
    Path(child.name+'.out').write_text(NATIVE)
    def launcher_failed(self, job):
        error=RuntimeError(f'Molpro sample {job.name} exited with status 1.')
        error.native_returncode=1;error.native_output_current=False
        raise error
    monkeypatch.setattr(type(sample),'_run_basis_job',launcher_failed)
    with pytest.raises(RuntimeError):
        sample._run_with_solver_recovery(child)
    assert not Path('.kinbot_sample_attempts').exists()
