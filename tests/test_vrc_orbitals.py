from types import SimpleNamespace
from pathlib import Path

from ase import Atoms
from ase.io import write
import numpy as np
import pytest

from kinbot.anl.dispatch import _render_input, validate_spec
from kinbot.vrc_orbitals import (
    molpro_pivot_spec, molpro_pivot_settings, validate_molpro_pivot)
from kinbot.vrc_tst_scan import VTS


def test_gaussian_fixed_cube_conversion_runs_in_local_workspace(tmp_path, monkeypatch):
    import os
    import sys
    from kinbot.anl import dispatch
    from kinbot.vrc_orbitals import fixed_pivot_spec
    from tests.test_orbital_cube import write_cube
    from kinbot.parameters import Parameters
    inp = tmp_path/'input.json'
    inp.write_text('{"barrier_threshold": 100, "theory_preset": "b3lyp-b2plyp-anl"}')
    par = Parameters(str(inp)).par
    atoms = Atoms('CO', positions=[[0, 0, 0], [1.13, 0, 0]])
    frag = SimpleNamespace(atom=atoms.get_chemical_symbols(), geom=atoms.positions, mult=1, charge=0)
    spec = fixed_pivot_spec(frag, par)
    task = spec['tasks'][0]
    task['resources'] = {'cores': 1, 'memory_mb': 512, 'walltime': '00:10:00'}
    dispatch.validate_spec(spec)
    task['command'] = [sys.executable, '-c',
        "from pathlib import Path; Path('pivot.chk').write_text('native checkpoint'); "
        "print('Normal termination')"]
    cube = write_cube(tmp_path/'source.cube', atoms=[(6, atoms.positions[0]), (8, atoms.positions[1])])
    binaries = tmp_path/'bin'; binaries.mkdir()
    work = tmp_path/'task'; work.mkdir()
    (work/'pivot.com').write_text('fixed geometry input\n')
    write(work/'geometry.xyz', atoms)
    scratch = tmp_path/'node-local'; scratch.mkdir()
    runtime = {'scratch_directory': str(scratch)}
    for name, action in (
        ('formchk', "Path('pivot.fchk').write_text(str(Path.cwd()))"),
        ('cubegen', f"Path('homo.cube').write_text({cube.read_text()!r})")):
        executable = binaries/name
        executable.write_text(f'#!{sys.executable}\nfrom pathlib import Path\n{action}\n')
        executable.chmod(0o755)
    env = dict(os.environ, PATH=str(binaries)+os.pathsep+os.environ['PATH'])
    monkeypatch.setattr(dispatch, 'qc_runtime_environment', lambda *a, **kw: (env, runtime))
    monkeypatch.setattr(dispatch, 'cleanup_qc_runtime', lambda runtime: None)
    dispatch._run_external(work, {'task': task})
    assert Path((work/'pivot.fchk').read_text()) == scratch/'work'
    assert (work/'homo.cube').read_text() == cube.read_text()
    assert runtime['native_work_directory'] == str(scratch/'work')
    (work/'pivot.log').write_text('Convergence criterion not met.\nNormal termination\n')
    with pytest.raises(RuntimeError, match='did not converge'):
        dispatch._validate_external_outputs(work, task)


def test_restart_controls_validate_reviewed_examples_and_only_change_geometry(tmp_path, monkeypatch):
    import json
    from kinbot import vrc_restart
    from kinbot.parameters import Parameters
    monkeypatch.setattr('kinbot.anl.site._executable', lambda *a: '/fake/tool')
    repo = Path(__file__).resolve().parents[1]
    for case in ('ch3s', 'ch3ss'):
        path = repo/f'examples/anl/{case}_anl1_f12/{case}.json'
        updates = vrc_restart.settings('avdz-f12', 'asymptotic_scan', 'L2')
        vrc_restart.validate_input(path, updates)
        par = Parameters(str(path)).par
        assert par['vrc_tst_scan_basis'] == 'avdz-f12'
        assert par['vrc_tst_high_basis'] == ['avdz-f12', 'avtz-f12']
        assert par['vrc_tst_high_extrapolation_power'] == 3.7
        spin = par['vrc_tst_spin_splitting']
        assert spin['reference']['basis'] == ['avdz-f12', 'avtz-f12']
        assert spin['reference']['extrapolation_power'] == 3.7
        assert all(level['basis'] == ['avdz-f12', 'avtz-f12']
                   and level['extrapolation_power'] == 3.7 for level in spin['gap_levels'])
        assert par['vrc_tst_scan_ipea'] == .25 and par['vrc_tst_scan_shift'] == .1
        assert par['vrc_tst_fragment_geometry'] == 'asymptotic_scan'
        assert par['vrc_tst_pivot_level'] == 'L2'
        previous = json.loads(path.read_text())
        previous['vrc_tst_scan_basis'] = 'avtz-f12'
        stage = previous['vrc_tst_surface_model']['channels'][0]['bottlenecks'][0]
        stage['parameters']['vrc_tst_scan_basis'] = 'avtz-f12'
        result = vrc_restart.revised_input(previous, updates)
        assert result['vrc_tst_surface_model']['channels'][0]['bottlenecks'][0]['parameters']['vrc_tst_scan_basis'] == 'avdz-f12'
        assert previous['vrc_tst_scan_basis'] == 'avtz-f12'
        with pytest.raises(ValueError, match='no orbital record'):
            vrc_restart.validate_input(path, vrc_restart.settings(
                reference_bases=['avtz-f12', 'avqz-f12'], reference_power=3.7))


def test_reference_cbs_migration_preserves_reviewed_spaces_and_gap_controls():
    import json
    from copy import deepcopy
    from kinbot import vrc_restart
    repo = Path(__file__).resolve().parents[1]
    original = json.loads((repo/'examples/anl/ch3ss_anl1_f12/ch3ss.json').read_text())
    stage = original['vrc_tst_surface_model']['channels'][0]['bottlenecks'][0]['parameters']
    stage['vrc_tst_spin_splitting'] = deepcopy(original['vrc_tst_spin_splitting'])
    before = deepcopy(original)
    updates = vrc_restart.settings(reference_bases=['avtz-f12', 'avqz-f12'], reference_power=4.0)
    result = vrc_restart.revised_input(original, updates)
    assert original == before
    targets = [result, result['vrc_tst_surface_model']['channels'][0]['bottlenecks'][0]['parameters']]
    for target in targets:
        ref = target['vrc_tst_spin_splitting']['reference']
        assert ref == {'method': 'uccsd(t)-f12b', 'basis': ['avtz-f12', 'avqz-f12'],
                       'extrapolation_power': 4.0}
        # Restore only this requested change; everything else must be identical.
        target['vrc_tst_spin_splitting']['reference'] = deepcopy(before['vrc_tst_spin_splitting']['reference'])
    assert result == before
    single = vrc_restart.revised_input(original, vrc_restart.settings(reference_bases=['avtz-f12']))
    assert single['vrc_tst_spin_splitting']['reference'] == {
        'method': 'uccsd(t)-f12b', 'basis': 'avtz-f12'}


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


def test_l2_pivot_is_fixed_geometry_and_inherits_full_profile(tmp_path):
    import json
    from kinbot.parameters import Parameters
    from kinbot.vrc_orbitals import fixed_pivot_spec
    path = tmp_path/'input.json'
    path.write_text(json.dumps(dict(barrier_threshold=100,
        theory_preset='b3lyp-b2plyp-anl', vrc_tst_fragment_geometry='asymptotic_scan',
        vrc_tst_pivot_level='L2', anl_scf_convergence=1.e-8)))
    par = Parameters(str(path)).par
    spec = fixed_pivot_spec(product(), par)
    task = spec['tasks'][0]
    task['resources'] = {'cores': 2, 'memory_mb': 4000, 'walltime': '01:00:00'}
    validate_spec(spec)
    deck = task['input_template'].lower()
    assert 'b2plyp/cc-pvtz' in deck and 'empiricaldispersion(gd3bj)' in deck
    assert 'integral(ultrafine)' in deck and 'scf(xqc,conver=8)' in deck
    assert all(word not in deck for word in ('opt(', 'freq', 'force', 'geom=check', 'guess(read)'))
    assert task['backend'] == 'gaussian'
    assert task['gaussian_orbital_export']['orbital'] == 'AMO=HOMO'
    assert spec['intent']['fixed_geometry'] and spec['intent']['orbital_only']


def test_custom_molpro_pivot_does_not_reoptimize():
    from kinbot.vrc_orbitals import fixed_pivot_spec
    spec = fixed_pivot_spec(product(), parameters(vrc_tst_pivot_level='custom'))
    task = spec['tasks'][0]
    assert 'optg' not in task['input_template'] and 'forces' not in task['input_template']
    assert 'type=alpha' in task['input_template']
    assert 'final.xyzgrad' not in task['required_outputs']
    assert spec['intent']['fixed_geometry']


def test_asymptotic_mode_never_launches_legacy_fragment_optimization():
    vts = VTS(SimpleNamespace(), parameters(vrc_tst_fragment_geometry='asymptotic_scan'), None)
    # No legacy QC accessor is available; save_products must defer to the MR path.
    vts.save_products(['not-yet-resolved'])
