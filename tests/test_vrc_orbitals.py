from types import SimpleNamespace

from ase import Atoms
from ase.io import write
import numpy as np
import pytest

from kinbot.anl.dispatch import _render_input, validate_spec
from kinbot.vrc_orbitals import (
    molpro_pivot_spec, molpro_pivot_settings, validate_molpro_pivot)
from kinbot.vrc_tst_scan import VTS


def parameters(**changes):
    return {
        'vrc_tst_pivot_backend': 'molpro',
        'vrc_tst_pivot_method': 'UKS/B3LYP3',
        'vrc_tst_pivot_basis': 'cc-pVDZ', 'pp_orient': 'homo',
        'vrc_tst_walltime': '04:00:00', 'single_point_ppn': 8,
        'vrc_tst_min_stack_mw': 1024, 'queue_name': 'compute',
        **changes,
    }


def product(symbols=('C', 'H', 'H', 'H'), mult=2):
    return SimpleNamespace(atom=list(symbols), charge=0, mult=mult,
                           geom=np.zeros((len(symbols), 3)))


def test_molpro_pivot_keeps_charge_spin_alpha_selection_and_native_artifacts():
    spec = molpro_pivot_spec(product(), parameters())
    task = spec['tasks'][0]
    task['resources'] = {'cores': 2, 'memory_mb': 24000,
                         'walltime': '04:00:00', 'partition': 'compute'}
    validate_spec(spec)
    rendered = _render_input(task['input_template'],
        Atoms('CH3', positions=np.zeros((4, 3))), spec['molecule'],
        task['resources'], backend='molpro')
    assert '{uks,b3lyp3,energy=1.d-8,maxit=200;wf,9,1,1;' in rendered
    assert 'orbital,5.1,record=2200.2,type=alpha' in rendered
    assert 'optg,maxit=100,gradient=' in rendered
    assert 'forces\nput,xyzgrad,final.xyzgrad' in rendered
    assert task['geometry_output'] == 'final.xyz'
    assert {'homo.cube', 'final.xyzgrad', 'final.xyz'} <= set(task['required_outputs'])


def test_atomic_pivot_has_no_optimization_or_unneeded_cube():
    spec = molpro_pivot_spec(product(('S',), mult=3), parameters())
    task = spec['tasks'][0]
    assert 'optg' not in task['input_template']
    assert 'cube,' not in task['input_template']
    assert task['required_outputs'] == ['pivot.out', 'final.xyz']


@pytest.mark.parametrize('kwargs', [
    {'vrc_tst_pivot_method': 'ub3lyp'},
    {'vrc_tst_pivot_method': 'UKS/B3LYP3;hf'},
    {'vrc_tst_pivot_basis': 'cc-pVDZ\nhf'},
    {'vrc_tst_pivot_grid_step': float('nan')},
    {'vrc_tst_pivot_fmax': 0},
    {'vrc_tst_pivot_maxit': 101},
])
def test_molpro_pivot_validates_native_options(kwargs):
    with pytest.raises(ValueError):
        molpro_pivot_settings(parameters(**kwargs))


def test_molpro_pivot_rejects_restricted_open_shell():
    with pytest.raises(ValueError, match='Open-shell'):
        molpro_pivot_spec(product(), parameters(vrc_tst_pivot_method='RKS/B3LYP3'))


def test_molpro_vrc_does_not_probe_gaussian_tools(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    def no_gaussian(*args):
        pytest.fail('Molpro pivots must not look for formchk or cubegen')
    monkeypatch.setattr('kinbot.anl.site._executable', no_gaussian)
    vts = VTS(SimpleNamespace(reac_obj=[], wellorts=1, name='parent', chemid='parent'), parameters(
        rotdpy_run=0, vrc_tst_scan={}, vrc_tst_noscan={}), None)
    vts.calculate_correction_potentials()


def test_fragment_handoff_uses_molpro_geometry_once_per_fragment(monkeypatch):
    calls = []
    def run(fragment, par, path):
        calls.append(path)
        return Atoms('CH3', positions=np.ones((4, 3)))
    monkeypatch.setattr('kinbot.vrc_orbitals.prepare_molpro_fragment', run)
    frag = product()
    frag.name, frag.wellorts = 'methyl', 1
    vts = VTS(SimpleNamespace(), parameters(), None)
    vts.scan_reac = {'a': SimpleNamespace(products=[frag]),
                     'b': SimpleNamespace(products=[frag])}
    vts.save_products(['a', 'b'])
    assert calls == ['vrctst/methyl_vts']
    assert np.array_equal(frag.geom, np.ones((4, 3)))


def test_pivot_final_gradient_is_required_even_after_normal_termination(tmp_path):
    par = parameters(pp_orient='geometric')
    spec = molpro_pivot_spec(product(('H', 'H'), mult=1), par)
    (tmp_path/'pivot.out').write_text('Molpro calculation terminated\n')
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, .74]])
    write(tmp_path/'final.xyz', atoms)
    (tmp_path/'final.xyzgrad').write_text(
        '2\nFinal gradient\nH 0 0 0 0 0 0 .1\nH 0 0 .74 0 0 0 -.1\n')
    with pytest.raises(ValueError, match='not converged'):
        validate_molpro_pivot(tmp_path, spec)
    (tmp_path/'final.xyzgrad').write_text(
        '2\nFinal gradient\nH 0 0 0 0 0 0 .00001\nH 0 0 .74 0 0 0 -.00001\n')
    assert np.allclose(validate_molpro_pivot(tmp_path, spec).positions, atoms.positions)
    (tmp_path/'pivot.out').write_text('NO CONVERGENCE!\nMolpro calculation terminated\n')
    with pytest.raises(ValueError, match='unconverged'):
        validate_molpro_pivot(tmp_path, spec)
