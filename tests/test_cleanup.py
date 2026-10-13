"""Cleanup must not discard completed calculations or mistake paths for methods."""
import logging

import pytest

from kinbot.utils import clean_files


GAUSSIAN_LOG = """ Entering Gaussian System
 Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          1           0        0.000000    0.000000    0.000000
      2          1           0        0.000000    0.000000    0.740000
 ---------------------------------------------------------------------
 SCF Done:  E(RHF) =  -1.0000000000     A.U. after 1 cycles
"""
NORMAL = ' Normal termination of Gaussian 16\n'
ERROR = ' Error termination via Lnk1e\n'


def calculation(tmp_path, route, ending, name='hydrogen'):
    # The checkpoint path and title deliberately mention AM1 even for DFT.
    (tmp_path / f'{name}.com').write_text(
        f'%chk={tmp_path}/methanol_am1_uma_s_m062x/{name}\n'
        f'{route}\n\nam1 comparison\n\n0 1\nH 0 0 0\nH 0 0 .74\n\n')
    logfile = tmp_path / f'{name}.log'
    logfile.write_text(GAUSSIAN_LOG + ending)
    return logfile


@pytest.mark.parametrize('diagnostic', [False, True])
@pytest.mark.parametrize('qc', ['gauss', 'fc'])
@pytest.mark.parametrize('route, am1', [
    ('#P M062X/6-31G(d) ! am1 comparison\nopt\nfreq', False),
    ('#P AM1D\nopt\nfreq', False),
    ('#p am1\nopt\nfreq', True),
    ('#P opt\n AM1\nfreq', True),
    ('#P AM1! route comment\nopt\nfreq', True),
])
@pytest.mark.parametrize('complete', [False, True])
def test_cleanup_uses_method_and_completion(tmp_path, monkeypatch, caplog,
                                            diagnostic, qc, route, am1, complete):
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO, logger='KinBot')
    logfile = calculation(tmp_path, route,
                          (NORMAL if complete else ERROR) + 'done\n')
    original = logfile.read_bytes()

    clean_files(diagnostic=diagnostic, qc=qc)

    unsuccessful_am1 = am1 and not complete
    assert logfile.exists() == (diagnostic or not unsuccessful_am1)
    if logfile.exists():
        assert logfile.read_bytes() == original
    messages = [r.message for r in caplog.records if 'AM1' in r.message]
    assert bool(messages) == unsuccessful_am1


@pytest.mark.parametrize('diagnostic', [False, True])
@pytest.mark.parametrize('ending, complete', [
    (NORMAL, True),
    (NORMAL + 'done\n\n', True),
    (NORMAL + ERROR + 'done\n', False),
    (NORMAL + ' Entering Link 1\n SCF iteration in progress\n', False),
])
def test_cleanup_checks_the_final_gaussian_calculation(tmp_path, monkeypatch,
                                                       caplog, diagnostic,
                                                       ending, complete):
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO, logger='KinBot')
    logfile = calculation(tmp_path, '#P AM1 opt freq', ending)

    clean_files(diagnostic=diagnostic, qc='gauss')

    assert logfile.exists() == (diagnostic or complete)
    assert any('AM1' in r.message for r in caplog.records) == (not complete)


def test_missing_output_is_not_reported_as_a_stuck_calculation(tmp_path,
                                                             monkeypatch,
                                                             caplog):
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO, logger='KinBot')
    logfile = calculation(tmp_path, '#p am1', ERROR)
    logfile.unlink()

    clean_files(diagnostic=True, qc='gauss')

    assert not any('AM1' in r.message for r in caplog.records)
