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
    assert 'scf(xqc,conver=10,Symm,NoVarAcc,NoIncFock,MaxCycle=256)' in decks[0]
    assert 'scf(conver=10,CDIIS,MaxCycle=256,Symm,NoVarAcc,NoIncFock)' in decks[1]
    for text in decks:
        assert 'scf_recovery' not in text and 'B3LYP/cc-pVDZ' in text
        assert 'integral(UltraFine)' in text and '\n0 3\n' in text


def test_atom_policy_is_general_and_can_be_disabled():
    from io import StringIO
    from ase import Atoms
    from kinbot.ase_modules.io.gaussian import write_gaussian_in
    from kinbot.gaussian_input import atomic_scf_parameters, merge_gaussian_kwargs
    original = {'Symm': 'None', 'SCF': 'XQC,Conver=10,NoSymm,MaxCycle=512',
                'integral': 'UltraFine', 'guess': 'Read', 'chk': 'example'}
    effective = atomic_scf_parameters(original)
    assert effective['scf'] == 'XQC,Conver=10,MaxCycle=512,Symm,NoVarAcc,NoIncFock'
    assert effective['guess'] == 'Read' and effective['integral'] == 'UltraFine'
    assert 'SCF' in original and 'Symm' in original
    assert atomic_scf_parameters(effective) == effective
    assert merge_gaussian_kwargs({}, {'atomic_scf_symmetry': False}) == {'atomic_scf_symmetry': False}
    for symbol, mult in [('S', 3), ('O', 3), ('H', 2), ('Ne', 1)]:
        for enabled in (True, False):
            stream = StringIO()
            write_gaussian_in(stream, Atoms(symbol, positions=[[1, 2, 3]]),
                method='B3LYP', basis='cc-pVDZ', charge=0, mult=mult,
                atomic_scf_symmetry=enabled, **original)
            native = stream.getvalue()
            assert 'atomic_scf_symmetry' not in native
            assert ('Symm(None)' in native) == (not enabled)
            assert ('NoIncFock' in native) == enabled
            assert f'\n0 {mult}\n' in native
            assert 'Conver=10' in native and 'PG=' not in native
    stream = StringIO()
    write_gaussian_in(stream, Atoms('OH'), method='B3LYP', basis='cc-pVDZ', mult=2, **original)
    assert 'Symm(None)' in stream.getvalue() and 'NoIncFock' not in stream.getvalue()


# Condensed native sulfur evidence, not a mocked read_results implementation.
ORIENTATION = ''' Input orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1         16           0        0.000000    0.000000    0.000000
 ---------------------------------------------------------------------
'''
UNCONVERGED = ''' RMSDP=3.86D-06 MaxDP=3.22D-05 DE=4.48D-11 OVMax=9.68D-05
 >>>>>>>>>> Convergence criterion not met.
 SCF Done:  E(UB3LYP) = -398.125080444 A.U. after 257 cycles
 NFock=256 Conv=0.39D-05 -V/T=2.0022
'''
CONVERGED = ''' RMSDP=1.26D-11 MaxDP=1.32D-10 DE=0.00D+00 OVMax=2.11D-13
 SCF Done: E(UB3LYP) = -398.125080451 A.U. after 9 cycles
 NFock=9 Conv=0.13D-10 -V/T=2.0022
'''
QC_CONVERGED = ''' Quadratic Convergence SCF Method.
 Rotation gradient small -- convergence achieved.
 SCF Done: E(UB3LYP) = -398.125080452 a.u. after 1 cycles
'''


@pytest.mark.parametrize('tail', ['', 'Normal termination of Gaussian 16\n',
    r'1\1\GINC'+'\nNormal termination of Gaussian 16\n',
    ORIENTATION + QC_CONVERGED + 'Normal termination of Gaussian 16\n'])
def test_reader_rejects_native_nonconvergence_even_with_zero_exit(tail):
    from io import StringIO
    from kinbot.ase_modules.io.gaussian import read_gaussian_out
    with pytest.raises(RuntimeError, match='Gaussian SCF did not converge'):
        read_gaussian_out(StringIO(ORIENTATION + UNCONVERGED + tail))


@pytest.mark.parametrize('diagnostic', [False, True])
def test_startup_cleanup_retains_rejected_scf_for_recovery(
        tmp_path, monkeypatch, caplog, diagnostic):
    from kinbot.utils import clean_files
    from kinbot.ase_modules.io.formats import read
    monkeypatch.chdir(tmp_path)
    path = tmp_path/'330330000000000000002_well.log'
    native = (' Entering Gaussian System\n'+ORIENTATION+UNCONVERGED+
              'Error termination via Lnk1e in /opt/gaussian/g16/l701.exe\n')
    path.write_text(native)
    checkpoint = tmp_path/'330330000000000000002_well.chk'
    checkpoint.write_bytes(b'failed-attempt evidence')
    with pytest.raises(RuntimeError, match='Gaussian SCF did not converge'):
        read(path)
    clean_files(qc='gauss', diagnostic=diagnostic)
    assert 'Startup cleanup retained unreadable QC output' in caplog.text
    assert path.read_text() == native
    assert checkpoint.read_bytes() == b'failed-attempt evidence'
    # Surviving cleanup must never turn the failed native output into a result.
    with pytest.raises(RuntimeError, match='Gaussian SCF did not converge'):
        read(path)


@pytest.mark.parametrize('scf', [CONVERGED, UNCONVERGED + QC_CONVERGED])
def test_reader_accepts_converged_diis_and_successful_xqc(scf):
    from io import StringIO
    from ase.units import Hartree
    from kinbot.ase_modules.io.gaussian import read_gaussian_out
    atom = read_gaussian_out(StringIO(ORIENTATION + scf + 'Normal termination of Gaussian 16\n'))
    assert atom.get_potential_energy()/Hartree == pytest.approx(-398.125080451, abs=2e-9)


def test_zero_exit_unconverged_retry_fails_and_clears_results(tmp_path, monkeypatch):
    from ase import Atoms
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(Gaussian, 'command', 'true')
    calls = []
    def execute(self):
        calls.append(1)
        if len(calls) == 1:
            Path(self.label+'.log').write_text('No lower point found -- run aborted.\n'
                'Error termination via Lnk1e in /opt/gaussian/g16/l508.exe')
            raise RuntimeError('Gaussian exited with error code 1')
        Path(self.label+'.log').write_text(ORIENTATION + UNCONVERGED +
                                         'Normal termination of Gaussian 16\n')
        self.results['energy'] = -1.  # must not survive a rejected native result
    monkeypatch.setattr(Gaussian, 'execute', execute)
    atom = Atoms('S')
    atom.calc = Gaussian(label='atom', method='B3LYP', basis='cc-pVDZ', mult=3,
                         scf='xqc,conver=10')
    with pytest.raises(RuntimeError, match='Gaussian SCF did not converge'):
        atom.get_potential_energy()
    assert calls == [1, 1]
    assert atom.calc.results == {}
    assert atom.calc.parameters['scf'] == 'xqc,conver=10'
    report = json.loads(Path('atom.scf_recovery.json').read_text())
    assert report['status'] == 'failed'
    assert 'SCF did not converge' in report['error']


def test_diatomic_symmetry_fallback_keeps_requested_physics_and_avoids_custom_states():
    from ase import Atoms
    from kinbot.ase_modules.calculators.gaussian import diatomic_recovery_parameters
    original=dict(method='B3LYP',basis='cc-pVDZ',mult=2,charge=0,
        scf='xqc,conver=8',Symm='None',guess='Read',integral='UltraFine',EmpiricalDispersion='GD3BJ')
    result=diatomic_recovery_parameters(original,Atoms('SH'))
    assert result['guess']=='Huckel' and 'Symm' not in result
    assert 'conver=8' in result['scf'] and result['scf'].endswith(',Symm')
    for k in ('method','basis','mult','charge','integral','EmpiricalDispersion'):
        assert result[k]==original[k]
    assert diatomic_recovery_parameters(original,Atoms('CH3S')) is None
    assert diatomic_recovery_parameters(dict(original,guess='Alter'),Atoms('SH')) is None
    assert diatomic_recovery_parameters(dict(original,mult=1),Atoms('SH')) is None


def test_link701_density_failure_uses_symmetry_once(tmp_path,monkeypatch):
    from ase import Atoms
    monkeypatch.chdir(tmp_path);monkeypatch.setattr(Gaussian,'command','true')
    decks=[]
    def execute(self):
        decks.append(Path(self.label+'.com').read_text())
        if len(decks)==1:
            Path(self.label+'.log').write_text('Convergence criterion not met.\nSCFChk: convergence required\nError termination in l701.exe')
            raise RuntimeError('Gaussian failed')
        Path(self.label+'.log').write_text('Normal termination of Gaussian')
    monkeypatch.setattr(Gaussian,'execute',execute)
    monkeypatch.setattr(Gaussian,'read_results',lambda self:self.results.update(energy=-1.))
    atoms=Atoms('SH',positions=[[.06,-2.55,.77],[-.52,-1.97,1.85]])
    atoms.calc=Gaussian(label='sh',method='B3LYP',basis='cc-pVDZ',mult=2,charge=0,
                        scf='xqc,conver=8',Symm='None',guess='Read')
    assert atoms.get_potential_energy()==-1.
    assert len(decks)==2 and 'guess(Huckel)' in decks[1] and 'Symm(None)' not in decks[1]
    assert 'conver=8' in decks[1] and 'CDIIS' in decks[1]
    assert atoms.calc.parameters['Symm']=='None'
    assert json.loads(Path('sh.scf_recovery.json').read_text())['recovery_kind']=='diatomic_symmetry'


def test_symmetry_force_table_is_already_in_input_frame(tmp_path,monkeypatch):
    import numpy as np
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator
    from kinbot.ase_modules.calculators import gaussian
    original=Atoms('SH',positions=[[.066735,-2.557438,.770077],[-.518863,-1.974158,1.854251]])
    native=Atoms('SH',positions=[[0,0,.080194],[0,0,-1.283101]])
    forces=np.array([[-.000012558,.000012509,.000023250],[.000012558,-.000012509,-.000023250]])
    native.calc=SinglePointCalculator(native,energy=-1.,forces=forces)
    monkeypatch.setattr(gaussian,'read',lambda *a,**k:native)
    calc=Gaussian(command='true');calc.atoms=original.copy();calc._diatomic_symmetry_retry_active=True
    calc.read_results()
    assert calc.results['forces']==pytest.approx(forces)
    native.positions[1,2]-=.01
    with pytest.raises(RuntimeError,match='changed geometry'):
        calc.read_results()
    assert calc.results=={}
