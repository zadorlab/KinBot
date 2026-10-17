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
               pp_orient='geometric',
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


def test_generated_cbs_input_is_valid_python_and_selects_worker_bridge(tmp_path, monkeypatch):
    par, data, path = sulfur_correction(tmp_path, sample_basis=['avtz-f12', 'avqz-f12'])
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
    assert (tmp_path / 'custom_rotd' / 'qu.tpl').is_file()


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
