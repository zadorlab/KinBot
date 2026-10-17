"""Execute generated workers through native-frequency parsing and publication.

The optimizer and QC executable are substituted; this exercises their Python
handoff, not Gaussian's electronic-structure results or Sella convergence.
"""
import json
import pickle
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase.calculators.calculator import Calculator, all_changes
from ase.db import connect
from ase.io import write

from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.frequencies import frequency_validation_error


@pytest.mark.parametrize('order,frequencies,rejected', [
    (0, [100., 200., 300.], False), (0, [-20., 200., 300.], False),
    (0, [-587.334, 200., 300.], True), (0, [-10., -20., 300.], True),
    (1, [-587., 200., 300.], False), (1, [-587., -20., 300.], False),
    (1, [-587., -80., 300.], True), (1, [10., 200., 300.], True),
    (0, None, True), (0, [], True), (0, [float('nan')], True),
    (1, [float('inf')], True),
])
def test_frequency_acceptance_preserves_signed_modes(order, frequencies, rejected):
    before = None if frequencies is None else list(frequencies)
    assert bool(frequency_validation_error(frequencies, order)) == rejected
    if before is not None:
        np.testing.assert_equal(frequencies, before)


@pytest.mark.parametrize('ts', [False, True])
@pytest.mark.parametrize('scenario', ['valid', 'wrong_order', 'optimizer_exception',
                                     'optimizer_limit', 'frequency_exception'])
def test_generated_worker_reports_actual_rejection(tmp_path, monkeypatch, capsys, ts, scenario):
    monkeypatch.chdir(tmp_path)
    from kinbot.ase_modules.calculators import gaussian
    import sella
    frames = []
    frequencies = [-587.334, 200., 300.] if ts else [100., 200., 300.]
    if scenario == 'wrong_order':
        frequencies = [100., 200., 300.] if ts else [-587.334, 200., 300.]

    class ReplayGaussian(Calculator):
        implemented_properties = ['energy', 'forces']

        def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            self.results = {'energy': -10., 'forces': np.zeros((len(atoms), 3))}
            frames.append(dict(self.parameters))
            if 'freq' in self.parameters:
                if scenario == 'frequency_exception':
                    raise RuntimeError('native Hessian evaluation failed')
                Path(self.label+'.log').write_text(
                    ' Frequencies --  '+' '.join(map(str, frequencies))+'\n'
                    ' Zero-point correction=  0.034729 (Hartree/Particle)\n'
                    ' Normal termination of Gaussian 16\n')

    class ConvergedOptimizer:
        def __init__(self, atoms, trajectory, logfile, **kwargs):
            self.atoms, self.trajectory, self.logfile = atoms, trajectory, logfile

        def run(self, **kwargs):
            if scenario == 'optimizer_exception':
                raise RuntimeError('force calculation failed')
            self.atoms.get_forces()
            write(self.trajectory, self.atoms)
            Path(self.logfile).write_text('Sella   0 19:11:56 -10.0 0.001 0.0\n')
            return scenario != 'optimizer_limit'

    monkeypatch.setattr(gaussian, 'Gaussian', ReplayGaussian)
    monkeypatch.setattr(sella, 'Sella', ConvergedOptimizer)
    Path('input.json').write_text(json.dumps({'theory_preset': 'uma-b2plyp-anl',
                                            'minimum_recovery_attempts': 0,
                                            'barrier_threshold': 1000}))
    qc = QuantumChemistry(Parameters('input.json').par)
    monkeypatch.setattr(qc.l2, 'submit_qc', lambda *a, **kw: 1)
    species = SimpleNamespace(chemid=123, name='probe', mult=1, charge=0, nel=10,
        natom=3, wellorts=int(ts), atom=['O', 'H', 'H'],
        geom=np.array([[0., 0., 0.], [.9, 0., 0.], [-.2, .8, 0.]]))
    if ts:
        qc.qc_opt_ts(species, species.geom, high_level=1)
        label = 'probe_high'
    else:
        qc.qc_opt(species, species.geom, high_level=1)
        label = '123_well_high'
    path = Path(label+'.py')
    exec(compile(path.read_text(), str(path), 'exec'), {'__name__': '__main__'})
    if ts:
        data = connect('kinbot.db').get(name=label).data
    else:
        with Path(label+'.pkl').open('rb') as stream:
            data = pickle.load(stream)['data']
    log = Path(label+'_sella.log').read_text()
    stderr = capsys.readouterr().err
    if scenario == 'valid':
        assert data['status'] == 'normal'
        assert list(data['frequencies']) == frequencies
        assert data['zpe'] == pytest.approx(.034729)
        assert not stderr
    else:
        assert data['status'] == 'error'
        assert data['error'] in stderr and data['error'] in log
        assert data['error']
        assert 'energy' not in data
        if scenario == 'wrong_order':
            assert 'Frequency validation failed' in data['error']
            assert list(data['frequencies']) == frequencies
    native = [p for p in frames if 'freq' in p]
    if scenario in ('optimizer_exception', 'optimizer_limit'):
        assert not native
    else:
        assert len(native) == 1
        assert native[0]['symmetry'] == 'PG=C1'
        assert 'opt' not in native[0] and 'force' not in native[0]
