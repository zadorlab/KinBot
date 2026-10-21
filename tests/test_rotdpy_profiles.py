"""Exercise the sulfur VRC-to-ROTD boundary with native-format fixtures."""

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

from ase import Atoms
from ase.units import Hartree
import numpy as np
import pytest

from kinbot.anl.extrapolation import two_point_cbs
from kinbot.parameters import Parameters
from kinbot.pes import (_rotdpy_sampling_calculator,
                        _rotdpy_validation_metadata, create_rotdpy_inputs)
from kinbot.vrc_correction import expected_levels, spin_splitting_settings


def sulfur_correction(tmp_path, name='ch3s', roots=3, sample_basis=None):
    source = (Path(__file__).parents[1] / 'examples' / 'anl' /
              f'{name}_anl1_f12' / f'{name}.json')
    par = Parameters(str(source), show_warnings=False).par
    if sample_basis is not None:
        par['vrc_tst_sample_basis'] = sample_basis
    par.update(input_file=str(tmp_path / 'input.json'), rotdpy_run=0,
               vrc_tst_state_model_file=f'{name}_state_model.json',
               pp_orient='geometric', rotdpy_walltime='01:00:00',
               vrc_tst_scan={'parent': ['parent_hom_sci_1_2']})
    models = {}
    for basis in ('avdz-f12', 'avtz-f12', 'avqz-f12'):
        wfu = tmp_path / f'{basis}.wfu'
        wfu.write_bytes(f'test reviewed orbitals {basis}'.encode())
        models[basis] = {
            'orbital_selection': {
                'mode': 'molpro_record', 'path': str(wfu), 'record': '2140.2',
                'sha256': hashlib.sha256(wfu.read_bytes()).hexdigest()},
            'state_average': [{'multiplicity': mult, 'roots': 3} for mult in (2, 4)],
        }
    (tmp_path / f'{name}_state_model.json').write_text(json.dumps({
        'schema': 2, 'active_space': par['vrc_tst_active_space'],
        'basis_models': models}))
    symbols = ['C', 'S', 'H', 'H', 'H']
    if name == 'ch3ss':
        symbols.insert(2, 'S')
    order = [0, *range(len(symbols) - 3, len(symbols)), *range(1, len(symbols) - 3)]
    split = 4
    if name == 'ch3ss' and roots == 6:
        order, split = [0, 1, 3, 4, 5, 2], 5
    fragment_symbols = [symbols[index] for index in order]
    data = {
        'vrc_molecule': {'symbols': symbols, 'charge': 0, 'multiplicity': 2},
        'sampling_atom_order': order, 'state_average_roots': {2: roots, 4: roots},
        'frags_atom': [fragment_symbols[:split], fragment_symbols[split:]],
        'frags_mult': [2, 3],
        'frags_geom': [[[float(i), 0., 0.] for i in range(split)],
                       [[0., 0., float(i)] for i in range(len(symbols) - split)]],
        'dist': [2.5, 4., 8., 30.],
        'e_samp': [-1., -.3, -.1, 0.], 'e_high': [-1.2, -.4, -.15, 0.],
        'e_inf_samp': -438., 'ra': [[0], [0]], 'scan_ref': [[0, 4]],
        'unique': [[[0]], [[0]]],
    }
    settings = spin_splitting_settings(par, parent_multiplicity=2,
                                        fragment_multiplicities=[2, 3])
    data['levels'] = expected_levels(
        par, symbols=symbols, charge=0, settings=settings,
        state_roots=data['state_average_roots'])
    path = tmp_path / 'vrctst' / 'corr_parent_hom_sci_1_2.json'
    path.parent.mkdir()
    path.write_text(json.dumps(data))
    return par, data, path


@pytest.mark.parametrize('name,roots', [('ch3s', 3), ('ch3ss', 1), ('ch3ss', 6)])
def test_production_provenance_uses_the_actual_sulfur_channel(tmp_path, name, roots):
    par, data, path = sulfur_correction(tmp_path, name, roots)
    metadata = _rotdpy_validation_metadata(par, data, path, False, 4)
    assert metadata['trusted_correction_level']['target_multiplicity'] == 2
    assert metadata['trusted_correction_level']['reference_multiplicity'] == 4
    calculator = _rotdpy_sampling_calculator(par, data, metadata)
    assert calculator['multiplicity'] == 2
    assert calculator['basis'] == 'avdz-f12'
    for model in calculator['state_models_by_basis'].values():
        assert [entry['roots'] for entry in model['state_average']] == [roots, roots]
    # Resolving per-basis orbital models must not append correction bases to
    # the sampling basis list in the immutable provenance declaration.
    assert data['levels']['sampling']['basis'] == 'avdz-f12'


@pytest.mark.parametrize('flux_target', [1, 5])
def test_generated_cbs_input_is_valid_python_and_selects_worker_bridge(tmp_path, monkeypatch, flux_target):
    import ast
    par, data, path = sulfur_correction(tmp_path, sample_basis=['avtz-f12', 'avqz-f12'])
    par['rotdpy_flux_parameters']['flux_rel_err'] = flux_target
    monkeypatch.chdir(tmp_path)
    par['rotdpy_directory'] = 'custom_rotd'
    class SurfaceFixture:
        def __repr__(self):
            return 'Surface()'
    monkeypatch.setattr('kinbot.pp_settings.create_all_surf_for_dist',
                        lambda **kw: ([[1]], [[0]], [SurfaceFixture()]))
    generated = create_rotdpy_inputs(
        par, [['parent', 'parent_hom_sci_1_2', ['methyl', 'sulfur'], 0.]], [],
        correction_root=tmp_path)
    assert generated == ['custom_rotd/parent_hom_sci_1_2.py']
    text = Path(generated[0]).read_text()
    compile(text, generated[0], 'exec')
    assert "'basis': ['avtz-f12', 'avqz-f12']" in text
    assert 'from kinbot.rotdpy_molpro import KinBotMultiSample as MultiSample' in text
    assert 'calc["kinbot_vrc"]' in text
    assert "'sampling_atom_order': [0, 2, 3, 4, 1]" in text
    assignments = {node.targets[0].id: node.value for node in ast.parse(text).body
                   if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)}
    assert ast.literal_eval(assignments['flux_parameter']) == par['rotdpy_flux_parameters']
    metadata = assignments['validation_metadata']
    recorded_flux = next(value for key, value in zip(metadata.keys, metadata.values)
                         if ast.literal_eval(key) == 'flux_parameters')
    assert ast.literal_eval(recorded_flux) == par['rotdpy_flux_parameters']
    assert (tmp_path / 'custom_rotd' / 'qu.tpl').is_file()


@pytest.mark.parametrize('coordinate', ['reactive', 'legacy'])
def test_real_surface_builder_covers_r7_grid_and_custom_output_directory(tmp_path, monkeypatch, coordinate):
    import re
    par, data, path = sulfur_correction(tmp_path)
    grid = ([i / 10 for i in range(20, 81, 2)]
            + [i / 2 for i in range(17, 25)] + list(range(13, 21)))
    par.update(rotdpy_dist=grid, vrc_tst_scan_points=grid,
               _vrc_dividing_coordinate=coordinate, rotdpy_directory='custom_rotd')
    data['dist'] = [2., 8., 20., 30.]
    # Non-collinear methyl geometry: the real geometric pivot construction
    # runs here; no QC, surface builder or Surface constructor is mocked.
    data['frags_geom'] = [[[0., 0., 0.], [1.08, 0., 0.],
                          [-.54, .935, 0.], [-.54, -.935, 0.]], [[0., 0., 0.]]]
    path.write_text(json.dumps(data))
    monkeypatch.chdir(tmp_path)
    generated = create_rotdpy_inputs(par,
        [['parent', 'parent_hom_sci_1_2', ['methyl', 'sulfur'], 0.]], [], correction_root=tmp_path)
    text = Path(generated[0]).read_text()
    compile(text, generated[0], 'exec')
    assert set(map(float, re.findall(r'#Distance: ([\d.]+) Angstrom', text))) == set(grid)
    far = text[:re.search(r'#Distance: 20(?:\.0)? Angstrom', text).start()].split('#Surface id:')[-1]
    assert ('#Fragment 0: on atom' in far) == (coordinate == 'reactive')
    assert ('#Fragment 0: COM' in far) == (coordinate == 'legacy')
    assert list(Path('custom_rotd').glob('S_*_F_*.xyz'))
    assert not Path('rotdPy').exists()
    assert not re.search(r'\b(?:nan|inf)\b', text, re.I)


@pytest.mark.parametrize('height', [0., .3, -.3])
def test_geometric_pivots_are_finite_and_translation_invariant(height):
    from kinbot.fragments import Fragment
    geom = np.array([[0., 0., height], [1., 0., 0.],
                     [-.5, .866, 0.], [-.5, -.866, 0.]])
    bonds = np.zeros((1, 4, 4), dtype=int)
    bonds[0, 0, 1:] = bonds[0, 1:, 0] = 1
    def orientations(positions):
        return Fragment.pp_bipyramide_triangle_base(
            SimpleNamespace(geom=positions, bonds=bonds), 0)
    with np.errstate(all='raise'):
        actual = orientations(geom)
        translated = orientations(geom + [3., -2., 5.])
    assert np.all(np.isfinite(actual))
    assert np.allclose(actual, translated)
    assert np.allclose(np.linalg.norm(actual, axis=1), 1.)
    if height:
        assert actual.shape == (1, 3)
        assert actual[0, 2] * height > 0
    else:
        assert actual.shape == (2, 3)
        assert np.allclose(actual[0], -actual[1])


def test_explicit_reactive_grid_rejects_uncovered_pivot_range(tmp_path, monkeypatch):
    par, _, _ = sulfur_correction(tmp_path)
    par.update(rotdpy_dist=[3., 4., 5.], _vrc_dividing_coordinate='reactive',
               pp_oriented=[1., 2.], pp_on_atom=[5., 12.])
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match='No reactive-coordinate ROTD surface at 3.0'):
        create_rotdpy_inputs(par,
            [['parent', 'parent_hom_sci_1_2', ['methyl', 'sulfur'], 0.]], [], correction_root=tmp_path)


@pytest.mark.parametrize('name,roots,electrons,partition,sulfur_x', [
    ('ch3s', 3, 25, 'occ,14\nclosed,10', 12.),
    ('ch3ss', 1, 41, 'occ,22\nclosed,17', 12.),
    ('ch3ss', 6, 41, 'occ,22\nclosed,17', 3.),
])
@pytest.mark.parametrize('sample_basis', ['avtz-f12', ['avtz-f12', 'avqz-f12']])
def test_worker_evaluates_both_bases_with_correct_spin_and_parent_atom_order(
        tmp_path, monkeypatch, name, roots, electrons, partition, sulfur_x, sample_basis):
    pytest.importorskip('rotd_py')
    from kinbot.rotdpy_molpro import MolproSample, LegacyMolpro
    par, data, path = sulfur_correction(tmp_path, name, roots, sample_basis)
    metadata = _rotdpy_validation_metadata(par, data, path, False, 4)
    model = _rotdpy_sampling_calculator(par, data, metadata)
    symbols = [symbol for fragment in data['frags_atom'] for symbol in fragment]
    atoms = Atoms(symbols, positions=np.arange(3 * len(symbols)).reshape(len(symbols), 3))
    calc = {'scratch': str(tmp_path), 'processors': 1, 'mem': 300,
            'method': model['method'], 'basis': model['basis'], 'kinbot_vrc': model}
    monkeypatch.chdir(tmp_path)
    sample = MolproSample('sample', atoms, calc)
    sample.create_input()
    bases = model['basis'] if isinstance(model['basis'], list) else [model['basis']]
    for index, basis in enumerate(bases):
        deck = (tmp_path / f'sample_basis_{index}.inp').read_text()
        assert f'basis={basis}\n' in deck
        assert deck.index('file,2,') < deck.index('geometry={') < deck.index('basis=')
        assert 'gprint,orbitals,orben,civector,ref' in deck
        assert 'gthresh,thrprint=0.05,printci=0.05' in deck
        assert 'orbprint,2' in deck
        assert f'wf,{electrons},1,1' in deck and f'wf,{electrons},1,3' in deck
        assert partition in deck
        assert f'state,{roots}' in deck
        assert 'start,2140.2' in deck and 'diab,2140.2' in deck
        assert 'KBVRCSEL=energy(1)' in deck
        assert 'set,spin=1' in deck
        sulfur_line = [line for line in deck.splitlines() if line.startswith('S ')][0]
        assert float(sulfur_line.split()[1]) == sulfur_x
    calls = []
    def native_run(self, child):
        calls.append(child.name)
        value = (-438., -438.1)[len(calls) - 1]
        Path(f'{child.name}.out').write_text(
            'PROGRAM * RS2-F12\n'
            f' SETTING KB_GEOM_ENERGY = {value} AU\n'
            'Molpro calculation terminated\n')
    monkeypatch.setattr(MolproSample, '_run_basis_job', native_run)
    sample.run()
    assert calls == [f'sample_basis_{index}' for index in range(len(bases))]
    expected = (-438. if len(bases) == 1 else
                two_point_cbs(-438., -438.1, upper_cardinal=4, power=3.7))
    assert sample.read_energy() / Hartree == pytest.approx(expected)
    # A missing normal termination cannot create a valid sample energy.
    Path(f'sample_basis_{len(bases) - 1}.out').write_text('SETTING KB_GEOM_ENERGY = -438.1 AU\n')
    with pytest.raises(ValueError, match='normally'):
        sample.read_energy()


def test_sampling_rejects_changed_atom_mapping_and_changed_reviewed_wfu(tmp_path):
    par, data, path = sulfur_correction(tmp_path)
    metadata = _rotdpy_validation_metadata(par, data, path, False, 4)
    data['sampling_atom_order'] = [0, 1, 2, 3, 4]
    with pytest.raises(ValueError, match='atom order'):
        _rotdpy_sampling_calculator(par, data, metadata)
    (tmp_path / 'avtz-f12.wfu').write_bytes(b'changed reviewed orbitals')
    with pytest.raises(ValueError, match='SHA-256 mismatch'):
        _rotdpy_validation_metadata(par, data, path, False, 4)
