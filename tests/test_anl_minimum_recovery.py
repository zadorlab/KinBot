"""Physical model surfaces test ANL curvature recovery without native QC."""
import json

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.optimize import BFGS

from kinbot.minimum_recovery import (optimize_anl_minimum, projected_force_modes,
                                     validate_anl_recovery)


class BondSurface(Calculator):
    implemented_properties = ['energy', 'forces']

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        delta = atoms.positions[1] - atoms.positions[0]
        distance = np.linalg.norm(delta)
        q = distance - 1.2
        energy = q**4 - .1*q*q
        gradient = (4*q**3 - .2*q) * delta / distance
        self.results = {'energy': energy, 'forces': np.array([gradient, -gradient])}


def molecule():
    atoms = Atoms('HH', positions=[[1., 2., 3.], [1.72, 2.96, 3.]])
    atoms.calc = BondSurface()
    return atoms


def test_projected_force_hessian_excludes_five_diatomic_rigid_modes(tmp_path):
    atoms = molecule()
    original = atoms.positions.copy()
    frequencies, modes = projected_force_modes(atoms, .001, tmp_path/'curvature.json')
    assert len(frequencies) == 1 and frequencies[0] < -100
    assert np.linalg.norm(modes[0]) == pytest.approx(1.)
    assert np.linalg.norm(modes[0].sum(axis=0)) < 1.e-10
    np.testing.assert_array_equal(atoms.positions, original)
    data = json.loads((tmp_path/'curvature.json').read_text())
    assert data['rigid_rank'] == 5 and len(data['probes']) == 12
    assert data['status'] == 'complete'


def test_force_converged_saddle_recovers_then_checks_fresh_curvature(tmp_path):
    atoms = molecule()
    assert np.linalg.norm(atoms.get_forces()) < 1.e-12
    starts = []

    def optimize():
        starts.append(atoms.positions.copy())
        return BFGS(atoms, logfile=str(tmp_path/'optimization.log')).run(fmax=1.e-6)

    generated = optimize_anl_minimum(atoms, optimize, tmp_path, {})
    assert len(starts) == 2
    assert abs(atoms.get_distance(0, 1)-1.2) == pytest.approx(np.sqrt(.05), abs=1.e-4)
    report = json.loads((tmp_path/'minimum_recovery.json').read_text())
    assert report['status'] == 'minimum_curvature_accepted'
    assert report['attempts'][0]['lowest_frequency_cm_inverse'] < 0
    assert report['attempts'][1]['lowest_frequency_cm_inverse'] > 0
    assert len(report['attempts'][0]['trials']) == 2
    assert all((tmp_path/name).exists() for name in generated)


def test_bounded_recovery_never_releases_persistent_saddle(tmp_path):
    atoms = molecule()
    original = atoms.positions.copy()
    count = []

    def falsely_converged():
        count.append(1)
        atoms.positions = original.copy()
        return True

    with pytest.raises(RuntimeError, match='after bounded recovery'):
        optimize_anl_minimum(atoms, falsely_converged, tmp_path, {'attempts': 1})
    assert len(count) == 2
    assert json.loads((tmp_path/'minimum_recovery.json').read_text())['status'] == 'failed'


@pytest.mark.parametrize('invalid', [{'attempts': True}, {'attempts': -1},
    {'displacement_angstrom': 0}, {'difference_step_angstrom': float('nan')},
    {'imaginary_threshold_cm_inverse': False}, {'unexpected': 1}])
def test_recovery_controls_reject_bad_values(invalid):
    with pytest.raises(ValueError):
        validate_anl_recovery(invalid)


def test_anl_dispatch_checks_curvature_before_publishing_final_geometry(tmp_path, monkeypatch):
    from ase.io import read, write
    from kinbot.anl.dispatch import _run_ase_optimize
    atoms = molecule()
    write(tmp_path/'geometry.xyz', atoms)
    monkeypatch.setattr('kinbot.ase_modules.calculators.factory.build_calculator',
                        lambda *a, **kw: BondSurface())
    record = {'molecule': {'charge': 0, 'multiplicity': 1}, 'task': {
        'id': 'general_geometry', 'profile': {'calculator': 'gaussian',
         'method': 'HF', 'basis': 'sto-3g', 'command': 'g16', 'optimizer': 'sella'},
        'resources': {'cores': 1, 'memory_mb': 2000},
        'optimizer': {'fmax': 1.e-6, 'minimum_recovery': {}}}}
    details = _run_ase_optimize(tmp_path, record)
    final = read(tmp_path/'final.xyz')
    assert abs(final.get_distance(0, 1)-1.2) == pytest.approx(np.sqrt(.05), abs=1.e-4)
    assert 'minimum_recovery.json' in details['generated_files']


def test_disabled_extra_curvature_and_optimizer_failure_are_distinct(tmp_path):
    atoms = molecule()
    optimize_anl_minimum(atoms, lambda: True, tmp_path, {'attempts': 0})
    assert json.loads((tmp_path/'minimum_recovery.json').read_text())['status'] == 'disabled'
    assert not list(tmp_path.glob('minimum_curvature*'))
    with pytest.raises(RuntimeError, match='did not converge'):
        optimize_anl_minimum(atoms, lambda: False, tmp_path, {})


def test_curvature_probe_error_preserves_candidate_and_records_failure(tmp_path):
    class Failed(Calculator):
        implemented_properties = ['energy', 'forces']

        def calculate(self, *args, **kwargs):
            raise RuntimeError('native force calculation failed')

    atoms = molecule()
    atoms.calc = Failed()
    original = atoms.positions.copy()
    with pytest.raises(RuntimeError, match='native force'):
        optimize_anl_minimum(atoms, lambda: True, tmp_path, {})
    np.testing.assert_array_equal(atoms.positions, original)
    assert json.loads((tmp_path/'minimum_curvature_0.json').read_text())['status'] == 'failed'
    assert json.loads((tmp_path/'minimum_recovery.json').read_text())['status'] == 'failed'


def test_nonlinear_rigid_projection_leaves_three_internal_water_modes(tmp_path):
    # A central-force triangle has three positive internal stretch curvatures;
    # all six rigid motions have zero curvature at its equilibrium geometry.
    class Springs(Calculator):
        implemented_properties = ['energy', 'forces']

        def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            forces, energy = np.zeros((3, 3)), 0.
            for i, j in ((0, 1), (0, 2), (1, 2)):
                v = atoms.positions[j]-atoms.positions[i]
                r = np.linalg.norm(v)
                dr = r-np.linalg.norm(base[j]-base[i])
                energy += .5*dr*dr
                forces[i] += dr*v/r
                forces[j] -= dr*v/r
            self.results = {'energy': energy, 'forces': forces}

    base = np.array([[0., 0., 0.], [1., 0., 0.], [-.3, .9, 0.]])
    atoms = Atoms('OHH', positions=base)
    atoms.calc = Springs()
    frequencies, modes = projected_force_modes(atoms, .001, tmp_path/'curvature.json')
    assert len(frequencies) == 3 and np.all(frequencies > 0.)
    assert modes.shape == (3, 3, 3)
    assert json.loads((tmp_path/'curvature.json').read_text())['rigid_rank'] == 6
