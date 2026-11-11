import json
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from sella import Sella

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot.optimization_budget import run_optimizer
from kinbot.optimizer_cycle import FloorCycle, OptimizerCycleError, attach_cycle_monitor


def test_gaussian_energy_and_forces_share_native_request(tmp_path, monkeypatch):
    requests = []
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(Gaussian, 'command', 'true')
    def execute(calc):
        requests.append((tmp_path / 'test.com').read_text())
    def read_results(calc):
        calc.results = dict(energy=float(len(requests)), forces=np.zeros((2, 3)))
    monkeypatch.setattr(Gaussian, 'execute', execute)
    monkeypatch.setattr(Gaussian, 'read_results', read_results)
    atoms = Atoms('HH', positions=[[0, 0, 0], [0, 0, 1]])
    atoms.calc = Gaussian(label='test', method='B3LYP', basis='cc-pVDZ', scf='conver=8')
    def run(**kwargs):
        first = atoms.get_potential_energy()
        atoms.get_forces()
        assert atoms.get_potential_energy() == first
        return True
    assert run_optimizer(SimpleNamespace(atoms=atoms, run=run), fmax=.001, steps=10, chunks=1)
    assert len(requests) == 1 and '\nforce\n' in requests[0].lower()
    assert not atoms.calc._coupled_derivatives
    # Ordinary single-point jobs retain their original route.
    atoms.positions[1, 2] += .1
    atoms.get_potential_energy()
    assert '\nforce\n' not in requests[1].lower()


def test_derivative_mode_restored_on_exception_and_native_tasks(tmp_path, monkeypatch):
    monkeypatch.setattr(Gaussian, 'command', 'true')
    atoms = Atoms('H')
    atoms.calc = Gaussian(label=str(tmp_path/'job'), method='B3LYP', basis='cc-pVDZ')
    def fail(**kwargs):
        assert atoms.calc._coupled_derivatives
        raise RuntimeError('failed QC')
    with pytest.raises(RuntimeError, match='failed QC'):
        run_optimizer(SimpleNamespace(atoms=atoms, run=fail), fmax=.001, steps=10)
    assert not atoms.calc._coupled_derivatives
    atoms.calc.set(freq='')
    with atoms.calc.optimizer_derivatives():
        assert not atoms.calc._coupled_derivatives


def cycle_optimizer():
    return SimpleNamespace(nsteps=0, delta=1e-4, delta_min=1e-4, eta=1e-4,
                           rho=-1., pes=SimpleNamespace(curr={}))


def feed(monitor, count=12, drift=0., force=.024):
    opt = monitor.optimizer
    for _ in range(count):
        opt.nsteps += 1
        x = (opt.nsteps % 4) * opt.delta_min / 4 + drift * opt.nsteps
        signed = force if drift else force * (-1 if opt.nsteps % 4 < 2 else 1)
        opt.pes.curr = dict(x=np.array([x, 0, 0]), f=x*x, g=np.array([signed, 0, 0]))
        monitor()


def test_floor_response_is_bounded_and_preserves_acceptance(tmp_path):
    opt = cycle_optimizer()
    monitor = FloorCycle(opt, .0005, tmp_path/'cycle.json')
    feed(monitor)
    assert opt.delta_min == pytest.approx(1e-5)
    assert opt.eta == 1e-4 and monitor.fmax == .0005
    feed(monitor)
    assert opt.delta_min == pytest.approx(1e-6)
    with pytest.raises(OptimizerCycleError):
        feed(monitor)
    record = json.loads((tmp_path/'cycle.json').read_text())
    assert len(record['events']) == 3 and record['status'] == 'exhausted'
    resumed = cycle_optimizer()
    restored = FloorCycle(resumed, .0005, tmp_path/'cycle.json')
    resumed.delta = resumed.delta_min
    with pytest.raises(OptimizerCycleError):
        feed(restored)
    assert resumed.delta_min == pytest.approx(1e-6)


@pytest.mark.parametrize('setting', ['drift', 'converged', 'large_radius', 'good_model'])
def test_slow_progress_or_convergence_does_not_trigger(setting):
    opt = cycle_optimizer()
    monitor = FloorCycle(opt, .0005, None)
    if setting == 'large_radius': opt.delta = .1
    if setting == 'good_model': opt.rho = 1.
    feed(monitor, 30, drift=1e-4 if setting == 'drift' else 0.,
         force=.0001 if setting == 'converged' else .024)
    assert not monitor.events


class HarmonicBond(Calculator):
    implemented_properties = ['energy', 'forces']
    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        displacement = atoms.positions[1] - atoms.positions[0]
        length = np.linalg.norm(displacement)
        force = 100 * (length - 1.) * displacement / length
        self.results = dict(energy=50*(length-1.)**2, forces=np.array([force, -force]))


def test_normal_sella_entry_point_and_repeated_invocation(tmp_path):
    atoms = Atoms('HH', positions=[[0,0,0], [0,0,1.2]])
    atoms.calc = HarmonicBond()
    opt = Sella(atoms, order=0, internal=False, logfile=None)
    for _ in range(2):
        assert run_optimizer(opt, fmax=.0005, steps=40, chunks=1, report=tmp_path/'budget.json')
    assert atoms.get_distance(0,1) == pytest.approx(1., abs=1e-6)
    assert not opt._kinbot_floor_cycle.events
    assert sum(isinstance(o[0], FloorCycle) for o in opt.observers) == 1
    ts = Sella(atoms, order=1, internal=False, logfile=None)
    attach_cycle_monitor(ts, .0005, None)
    assert not hasattr(ts, '_kinbot_floor_cycle')
