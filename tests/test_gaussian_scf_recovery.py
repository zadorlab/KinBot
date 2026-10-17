from pathlib import Path
import json

import pytest

from kinbot.ase_modules.calculators.gaussian import (
    Gaussian, FileIOCalculator, qc_scf_failure, scf_recovery_parameters,
)


def test_solver_change_preserves_physics_and_threshold():
    original = dict(method='B3LYP', basis='cc-pVDZ', charge=0, mult=3,
                    scf='XQC,Conver=10', integral='UltraFine', EmpiricalDispersion='GD3BJ')
    result = scf_recovery_parameters(original)
    assert original['scf'] == 'XQC,Conver=10'
    assert result.pop('scf') == 'Conver=10,CDIIS,NoVarAcc,NoIncFock,MaxCycle=256'
    assert result == {k: v for k, v in original.items() if k != 'scf'}
    assert 'MaxCycle=512' in scf_recovery_parameters({'SCF': '(QC,Conver=12,MaxCycle=512)'})['scf']
    from kinbot.gaussian_input import merge_gaussian_kwargs
    assert merge_gaussian_kwargs({}, {'scf_recovery': True}) == {'scf_recovery': True}
    assert merge_gaussian_kwargs({}, {'scf_recovery': False}) == {'scf_recovery': False}
    from io import StringIO
    from ase import Atoms
    from kinbot.ase_modules.io.gaussian import write_gaussian_in
    stream = StringIO()
    write_gaussian_in(stream, Atoms('S'), method='B3LYP', basis='cc-pVDZ',
                      mult=3, scf_recovery=False)
    assert 'scf_recovery' not in stream.getvalue()


@pytest.mark.parametrize('text', ['QPErr\nError termination via l1.exe',
    'No space left on device', 'Convergence failure\nError termination via l502.exe',
    'No lower point found\nNormal termination of Gaussian'])
def test_does_not_retry_unrelated_or_successful_output(tmp_path, text):
    path = tmp_path/'test.log'
    path.write_text(text)
    assert not qc_scf_failure(path)


@pytest.mark.parametrize('second_fails,enabled', [(False, True), (True, True), (False, False)])
def test_one_bounded_retry_and_native_archive(tmp_path, monkeypatch, second_fails, enabled):
    monkeypatch.chdir(tmp_path)
    calls = []
    def fake_calculate(self, *args, **kwargs):
        calls.append(dict(self.parameters))
        Path(self.label+'.com').write_text(str(self.parameters))
        Path('test.chk').write_bytes(b'original wavefunction')
        if len(calls) == 1 or second_fails:
            Path(self.label+'.log').write_text('No lower point found -- run aborted.\n'
                'Error termination via Lnk1e in /opt/gaussian/g16/l508.exe')
            raise RuntimeError('Gaussian exited with error code 1')
        Path(self.label+'.log').write_text('SCF Done\nNormal termination of Gaussian 16')
        self.results = {'energy': -10.}
    monkeypatch.setattr(FileIOCalculator, 'calculate', fake_calculate)
    monkeypatch.setattr(Gaussian, 'command', 'true')
    calc = Gaussian(label='test', command='true', method='B3LYP', mult=3,
                    scf='xqc,conver=10', chk='test', scf_recovery=enabled)
    if second_fails or not enabled:
        with pytest.raises(RuntimeError, match='Gaussian exited'):
            calc.calculate()
    else:
        calc.calculate()
        assert calc.results['energy'] == -10.
    assert len(calls) == (2 if enabled else 1)
    assert calc.parameters['scf'] == 'xqc,conver=10'
    assert 'scf_recovery' not in calc.parameters  # never a Gaussian route keyword
    if enabled:
        assert calls[1]['scf'] == 'conver=10,CDIIS,NoVarAcc,NoIncFock,MaxCycle=256'
        report = json.loads(Path('test.scf_recovery.json').read_text())
        assert report['status'] == ('failed' if second_fails else 'complete')
        archived = Path(report['preserved_files']['test.log']['archive'])
        assert 'No lower point found' in archived.read_text()
        assert list(tmp_path.glob('restart_gaussian_scf/*/test.chk'))[0].read_bytes() == b'original wavefunction'


def test_successful_calculation_is_not_retried(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calls = []
    monkeypatch.setattr(Gaussian, 'command', 'true')
    monkeypatch.setattr(FileIOCalculator, 'calculate', lambda *a, **k: calls.append(1))
    Gaussian(command='true').calculate()
    assert calls == [1]
    assert not list(tmp_path.glob('*.scf_recovery.json'))


def test_retry_renders_native_route_through_real_calculation_path(tmp_path, monkeypatch):
    from ase import Atoms
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(Gaussian, 'command', 'true')
    decks = []
    def execute(self):
        decks.append(Path(self.label+'.com').read_text())
        if len(decks) == 1:
            Path(self.label+'.log').write_text('No lower point found -- run aborted.\n'
                'Error termination via Lnk1e in /opt/gaussian/g16/l508.exe')
            raise RuntimeError('Gaussian exited with error code 1')
        Path(self.label+'.log').write_text('Normal termination of Gaussian 16')
    monkeypatch.setattr(Gaussian, 'execute', execute)
    monkeypatch.setattr(Gaussian, 'read_results', lambda self: self.results.update(energy=-1.))
    atom = Atoms('S')
    atom.calc = Gaussian(label='atom', method='B3LYP', basis='cc-pVDZ',
                         mult=3, charge=0, scf='xqc,conver=10',
                         integral='UltraFine', EmpiricalDispersion='GD3BJ')
    assert atom.get_potential_energy() == -1.
    assert len(decks) == 2
    assert 'scf(xqc,conver=10)' in decks[0]
    assert 'scf(conver=10,CDIIS,NoVarAcc,NoIncFock,MaxCycle=256)' in decks[1]
    for text in decks:
        assert 'scf_recovery' not in text and 'B3LYP/cc-pVDZ' in text
        assert 'integral(UltraFine)' in text and '\n0 3\n' in text
