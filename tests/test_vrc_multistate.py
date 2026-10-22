"""Physical multistate fluxes stay separate from correction variants."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from kinbot.execution import file_sha256
from kinbot.vrc_multistate import (validate_surface_model, channel_for_reaction,
    state_parameters, calculate_channel, assemble_states, stage_channel, render_mess)


def model():
    def state(mult, root=1, offset=0.):
        return {'multiplicity': mult, 'root': root, 'energy_offset_cm_inverse': offset,
                'energy_reference_note': 'fixture common spin-free asymptote'}
    return {'schema': 1, 'reviewed': True, 'review_note': 'fixture state assignment',
        'channels': [{'bond': [1, 2], 'combination': 'statistical', 'equilibration': 'rapid',
            'bottlenecks': [
                {'name': 'inner', 'coordinate': 'reactive', 'distances_angstrom': [3., 4., 5.],
                 'states': [state(2), state(2, 2)]},
                {'name': 'outer', 'coordinate': 'com', 'distances_angstrom': [6., 8., 10.],
                 'states': [state(2), state(4, offset=100.)]},
            ]}]}


@pytest.mark.parametrize('change,match', [
    (lambda data: data.update(reviewed=False), 'reviewed=true'),
    (lambda data: data['channels'][0].pop('equilibration'), 'rapid equilibration'),
    (lambda data: data['channels'][0]['bottlenecks'][0]['states'][0].pop('energy_offset_cm_inverse'), 'energy offset'),
    (lambda data: data['channels'][0]['bottlenecks'][0]['states'].append(
        deepcopy(data['channels'][0]['bottlenecks'][0]['states'][0])), 'Duplicate'),
])
def test_model_rejects_unreviewed_or_ambiguous_chemistry(change, match):
    data = model()
    change(data)
    with pytest.raises(ValueError, match=match):
        validate_surface_model(data)


def test_channel_matching_uses_bond_and_input_terms_not_molecule_names(tmp_path):
    from test_rotdpy_profiles import sulfur_correction
    from kinbot.stationary_pt import StationaryPoint
    par, _, _ = sulfur_correction(tmp_path)
    path = tmp_path / 'surfaces.json'
    path.write_text(json.dumps(model()))
    par['vrc_tst_surface_model'] = None
    par['vrc_tst_surface_model_file'] = str(path)
    products = [StationaryPoint('a', 0, 2, smiles='[CH3]'),
                StationaryPoint('b', 0, 3, smiles='[S]')]
    for product in products:
        product.characterize()
    reaction = SimpleNamespace(scan_coo=[1, 0], products=products)
    selected, provenance = channel_for_reaction(par, reaction)
    assert selected['bottlenecks'][0]['states'][1]['root'] == 2
    assert provenance['sha256'] == file_sha256(path)
    data = model()
    data['channels'][0]['bottlenecks'][0]['states'][1]['root'] = 4
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='outside'):
        channel_for_reaction(par, reaction)


def test_scalar_parameters_preserve_parent_and_select_each_native_root(tmp_path):
    from test_rotdpy_profiles import sulfur_correction
    from kinbot.vrc_correction import resolved_state_model, spin_splitting_settings
    par, _, _ = sulfur_correction(tmp_path)
    before = deepcopy(par)
    channel = validate_surface_model(model())['channels'][0]
    inner = channel['bottlenecks'][0]
    leaf = state_parameters(par, inner, inner['states'][1], [2, 3])
    settings = spin_splitting_settings(leaf, parent_multiplicity=2, fragment_multiplicities=[2, 3])
    resolved = resolved_state_model(leaf, ['C', 'S', 'H', 'H', 'H'], 0, settings,
                                    roots={2: 3, 4: 3}, basis='avtz-f12')
    assert resolved['selected_roots'] == {2: 2, 4: 1}
    outer = channel['bottlenecks'][1]
    high = state_parameters(par, outer, outer['states'][1], [2, 3])
    assert high['mult'] == 4 and high['vrc_tst_correction_model'] == 'direct'
    assert high['_vrc_dividing_coordinate'] == 'com'
    assert par == before


def test_new_pes_well_cannot_inherit_other_wells_orbitals(tmp_path,monkeypatch):
    from test_rotdpy_profiles import sulfur_correction
    from kinbot import species_routing
    from kinbot.stationary_pt import StationaryPoint
    par,_,_=sulfur_correction(tmp_path)
    products=[StationaryPoint('a',0,2,smiles='[CH3]'),StationaryPoint('b',0,3,smiles='[S]')]
    for p in products: p.characterize()
    monkeypatch.setattr(species_routing,'routing_name',lambda parent:parent.route)
    reaction=SimpleNamespace(species=SimpleNamespace(route='new',smiles='C[S]'),scan_coo=[0,1],products=products)
    par.update(vrc_tst_surface_model=model(),vrc_tst_surface_model_file='',vrc_tst_surface_model_scope='original')
    with pytest.raises(ValueError,match='newly explored well'):
        channel_for_reaction(par,reaction)
    par['vrc_tst_surface_models']={'new':model()}
    with pytest.raises(ValueError,match='each new-parent bottleneck'):
        channel_for_reaction(par,reaction)
    for stage in par['vrc_tst_surface_models']['new']['channels'][0]['bottlenecks']:
        stage['parameters']={'vrc_tst_state_model_file':'new-reviewed-state.json'}
    assert channel_for_reaction(par,reaction)[0]['bond']==[1,2]


def test_correction_orchestrator_uses_independent_paths_and_original_species(tmp_path, monkeypatch):
    from kinbot.vrc_tst_scan import VTS
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.conformers import Conformers
    from kinbot.parameters import Parameters
    from kinbot.qc import QuantumChemistry
    from test_rotdpy_profiles import sulfur_correction
    monkeypatch.chdir(tmp_path)
    Path('profiles.json').write_text(json.dumps({
        'theory_preset': 'b3lyp-b2plyp-anl', 'barrier_threshold': 50.}))
    profiles = Parameters('profiles.json', show_warnings=False).par
    qc = QuantumChemistry(profiles)
    par, _, _ = sulfur_correction(tmp_path)
    products = [StationaryPoint('h', 0, 2, natom=1, atom=['H'], geom=np.zeros((1, 3))),
                StationaryPoint('s', 0, 3, natom=1, atom=['S'], geom=np.ones((1, 3)))]
    for product in products:
        product.confs = Conformers(product, profiles, qc)
    # Live conformer/QC objects must remain shared, even when copy support on
    # the router changes. The regression uses the actual production objects.
    def reject_service_copy(*args):
        pytest.fail('VRC copied a live conformer/QC controller')
    monkeypatch.setattr(Conformers, '__deepcopy__', reject_service_copy, raising=False)
    original_geometries = [p.geom.copy() for p in products]
    parent = SimpleNamespace(mult=2)
    reaction = SimpleNamespace(species=parent, products=products)
    vts = VTS(SimpleNamespace(), par, qc)
    vts.scan_reac['reaction'] = reaction
    calls, leaves = [], []
    def scan(child, names):
        leaf = child.scan_reac[names[0]]
        leaves.append(leaf)
        for index, product in enumerate(leaf.products):
            assert product is not products[index]
            assert product.confs is products[index].confs
            assert product.confs.qc is qc
            assert np.array_equal(product.geom, original_geometries[index])
            # Include in-place changes, not just assignment of new arrays.
            product.geom[:] += len(leaves)
        calls.append((names[0], leaf.species.mult,
                      child.par['vrc_tst_selected_roots']))
    def energies(child, names):
        Path(f'vrctst/corr_{names[0]}.json').write_text('{}')
    monkeypatch.setattr(VTS, 'do_scan', scan)
    monkeypatch.setattr(VTS, 'energies', energies)
    path = calculate_channel(vts, 'reaction', (validate_surface_model(model())['channels'][0],
                                             {'source': 'review.json', 'sha256': 'fixture'}))
    payload = json.loads(path.read_text())
    assert len(calls) == 4 and len(set(name for name, _, _ in calls)) == 4
    assert calls[-1][1] == 4 and calls[1][2]['2'] == 2
    assert parent.mult == 2
    for index, product in enumerate(products):
        assert np.array_equal(product.geom, original_geometries[index])
        for offset, leaf in enumerate(leaves, 1):
            assert np.array_equal(leaf.products[index].geom,
                                  original_geometries[index] + offset)
    assert payload['kind'] == 'multistate'
    assert all(file_sha256(state['correction']) == state['correction_sha256']
               for stage in payload['bottlenecks'] for state in stage['states'])


def test_flux_sum_aligns_thresholds_counts_spin_once_and_does_not_extrapolate(tmp_path):
    table = tmp_path / 'Ne.out'
    table.write_text('0 0\n100 1\n200 2\n300 3\n')
    entries = [{'source': table, 'multiplicity': 2, 'energy_offset_cm_inverse': 0.},
               {'source': table, 'multiplicity': 4, 'energy_offset_cm_inverse': 100.}]
    result = tmp_path / 'sum.out'
    info = assemble_states(entries, result, baseline_cm_inverse=0.)
    data = np.loadtxt(result)
    assert data[:, 1] == pytest.approx([0., 2., 8., 14.])
    assert info['upper_energy_cm_inverse'] == 300.
    assert data[-1, 0] == 300.  # no extrapolation of the doublet beyond its grid


def test_doublet_components_have_union_temperature_weights_not_an_extra_degeneracy(tmp_path):
    """An analytic shifted ramp has Laplace weight T^2 exp(-threshold/T).

    Different thresholds must produce a temperature-dependent contribution;
    equal thresholds must sum two doublets, never average or multiply twice.
    Integration uses cm^-1 energy and kT in the same units.
    """
    energy = np.arange(0., 40001., 1.)
    source = tmp_path/'ramp.out'
    np.savetxt(source, np.column_stack((energy,energy)))
    target = tmp_path/'sum.out'
    entries = [dict(source=source,multiplicity=2,energy_offset_cm_inverse=offset)
               for offset in (0.,500.)]
    info = assemble_states(entries,target,baseline_cm_inverse=0.)
    values = np.loadtxt(target)
    for temperature in (100.,300.,1000.,2000.):
        kt = .69503476*temperature
        # Common fragment vibration/rotor factors multiply both sides and
        # cancel in this comparison; the two state averages are irrelevant.
        measured = np.trapezoid(values[:,1]*np.exp(-values[:,0]/kt),values[:,0])
        expected = 2*kt**2*(1+np.exp(-500./kt))
        assert measured == pytest.approx(expected,rel=3e-5)
    entries[1]['energy_offset_cm_inverse']=0.
    assemble_states(entries,target,baseline_cm_inverse=0.)
    assert np.loadtxt(target)[100,1] == pytest.approx(4*100.)
    assert info['electronic_assembly']=='sum_of_shifted_microcanonical_fluxes'


def test_si_level_convention_counts_resolved_degeneracy_instead_of_extra_spin(tmp_path):
    contract = json.loads((Path(__file__).parents[1] / 'docs/NH2O_SI_verification.json').read_text())
    assert contract['status'] == 'passed'
    components = contract['entrance_barriers']['B1v']['components']
    source = tmp_path / 'Ne.out'
    source.write_text('-186 0\n0 100\n100 200\n200 300\n300 400\n')
    entries = [{'source': source, 'multiplicity': 2, 'energy_offset_cm_inverse': 0.,
                'electronic_levels': component['electronic_levels']} for component in components]
    target = tmp_path / 'sum.out'
    info = assemble_states(entries, target, baseline_cm_inverse=0.)
    table = np.loadtxt(target)
    # At the negative bound every component is zero; at E=0 only the
    # ten unsplit ground-level contributions and shifted excited tails enter.
    assert table[0, 0] == -186. and table[0, 1] == 0.
    assert info['ground_cm_inverse'] == -186.
    assert table[np.where(table[:, 0] == 0.)[0][0], 1] < 1800.
    assert len(contract['files_sha256']) == 10
    assert max(contract['canonical_level_counting_relative_errors']) < .001


def test_single_bottleneck_never_requires_or_creates_a_complex():
    data = model()
    channel = data['channels'][0]
    channel['combination'] = 'single'
    channel.pop('equilibration')
    channel['bottlenecks'] = [channel['bottlenecks'][0]]
    channel['bottlenecks'][0]['name'] = 'single'
    assert validate_surface_model(data)['channels'][0]['combination'] == 'single'
    result = render_mess({'combination': 'single', 'staged_bottlenecks': [{'name': 'single'}]},
                         lambda stage: 'RRHO physical state flux\n')
    assert 'Well' not in result and 'Variational' not in result


def test_unbounded_negative_state_flux_cannot_be_staged(tmp_path):
    source = tmp_path / 'Ne.out'
    source.write_text('-186 100\n0 200\n100 300\n')
    with pytest.raises(ValueError, match='zero-flux lower boundary'):
        assemble_states([{'source': source, 'multiplicity': 2,
                         'energy_offset_cm_inverse': 0.}],
                        tmp_path / 'sum.out', baseline_cm_inverse=0.)


def _completed_channel(tmp_path):
    model_file = tmp_path / 'surfaces.json'
    model_file.write_text(json.dumps(model()))
    correction = tmp_path / 'corr.json'
    correction.write_text('{}')
    payload = {'schema': 1, 'kind': 'multistate', 'reaction': 'r', 'combination': 'statistical',
        'surface_model': {'source': str(model_file), 'sha256': file_sha256(model_file)},
        'correction': str(correction), 'correction_sha256': file_sha256(correction), 'bottlenecks': []}
    for stage in validate_surface_model(model())['channels'][0]['bottlenecks']:
        entries = []
        for state in stage['states']:
            name = f'{stage["name"]}_m{state["multiplicity"]}_r{state["root"]}'
            inp = tmp_path / f'{name}.py'
            inp.write_text('dynamical_correction=0.85\n')
            out = tmp_path / f'kb_{name}'
            (out / 'output').mkdir(parents=True)
            files = [out / 'Ne_0.out', out / 'Ne_1.out', out / 'output' / 'surface_0.dat']
            files[0].write_text('0 0\n100 100\n200 200\n300 300\n')
            files[1].write_text('0 0\n100 1\n200 2\n300 3\n')
            files[2].write_text('Canonical 1 2\n')
            relative = [str(path.relative_to(tmp_path)) for path in files]
            inp.with_suffix('.rotdpy.json').write_text(json.dumps({
                'status': 'complete', 'surface_count': 1, 'result_files': relative,
                'result_sha256': {name: file_sha256(path) for name, path in zip(relative, files)}}))
            entries.append({**state, 'input': inp.name, 'input_sha256': file_sha256(inp),
                            'correction': str(correction), 'correction_sha256': file_sha256(correction)})
        payload['bottlenecks'].append({**stage, 'states': entries})
    manifest = tmp_path / 'r.multistate.json'
    manifest.write_text(json.dumps(payload))
    return manifest


def test_staging_selects_corrected_variant_and_revalidates_all_native_sources(tmp_path):
    manifest = _completed_channel(tmp_path)
    provenance = stage_channel(manifest, tmp_path / 'me')
    inner = np.loadtxt(tmp_path / 'me' / 'rotd_r_inner.out')
    assert inner[-1, 1] == 12.  # two doublets, corrected Ne=3 each
    block = render_mess(provenance,
        lambda stage: '    RRHO\n      Core Rotd\n        File ' + stage['file'] + '\n      End\n    End\n')
    assert '2TSMethod statistical' in block and 'OuterRRHO' in block
    assert block.count('Core Rotd') == 2
    assert len(provenance['staged_bottlenecks'][1]['states']) == 2
    # An alternate correction index was not counted as another physical state.
    table = tmp_path / 'kb_inner_m2_r2' / 'Ne_1.out'
    table.write_text('0 0\n100 999\n')
    with pytest.raises(RuntimeError, match='hash changed'):
        stage_channel(manifest, tmp_path / 'me')


@pytest.mark.parametrize('duplicate_input',[False,True])
def test_flux_manifest_cannot_count_a_surface_twice(tmp_path,duplicate_input):
    manifest=_completed_channel(tmp_path)
    payload=json.loads(manifest.read_text())
    states=payload['bottlenecks'][0]['states']
    if duplicate_input:
        states[1]['input']=states[0]['input']
        states[1]['input_sha256']=states[0]['input_sha256']
    else:
        states.append(deepcopy(states[0]))
    manifest.write_text(json.dumps(payload))
    with pytest.raises(ValueError,match='Duplicate electronic surface'):
        stage_channel(manifest,tmp_path/'me')


def test_mess_handoff_uses_unit_electronic_degeneracy_for_summed_flux(tmp_path, monkeypatch):
    from kinbot.mess import MESS
    from kinbot import kb_path
    manifest = _completed_channel(tmp_path)
    monkeypatch.chdir(tmp_path)
    writer = MESS.__new__(MESS)
    writer.par = {'rotdpy_mess_mode': 'required', 'pes': 0, 'rotdpy_directory': str(tmp_path)}
    writer.species = SimpleNamespace(mult=2)
    writer.rotdpy_sources = {}
    writer.rotdrrhotpl = Path(kb_path, 'tpl/mess_rotd_rrho.tpl').read_text()
    writer._stable_frequencies = lambda species: []
    writer.make_rotors = lambda *args, **kwargs: ''
    writer.make_freq = lambda *args: ''
    writer._parent_symmetry = lambda species: 1
    products = [SimpleNamespace(mult=2, chemid=1), SimpleNamespace(mult=3, chemid=2)]
    monkeypatch.setattr('kinbot.mess.routing_key', lambda species: str(species.chemid))
    block = writer._rotd_model(products, 'CH3S', 10., 0., 1., 1., SimpleNamespace(instance_name='r'))
    assert block.count('ElectronicLevels[1/cm]      1') == 2
    assert block.count('            0 1') == 2
    assert '2TSMethod statistical' in block
    assert writer.rotdpy_sources['r']['kind'] == 'multistate'


def test_requested_multistate_cannot_reuse_an_old_scalar_result(tmp_path):
    from kinbot.mess import MESS
    writer = MESS.__new__(MESS)
    writer.par = {'rotdpy_mess_mode': 'auto', 'pes': 0,
                  'rotdpy_directory': str(tmp_path),
                  'vrc_tst_surface_model_file': 'reviewed.json'}
    (tmp_path / 'r.py').write_text('# old scalar input\n')
    (tmp_path / 'r.rotdpy.json').write_text('{}')
    with pytest.raises(FileNotFoundError, match='state-resolved'):
        writer._rotd_model([], '', 0., 0., 1., 1., SimpleNamespace(instance_name='r'))


def test_full_pes_requires_worker_context_for_selected_rotd(tmp_path):
    from kinbot.rotdpy_mess import resolve_models
    with pytest.raises(ValueError, match='regenerate the PES worker'):
        resolve_models('RRHO\nCore PhaseSpaceTheory\n', tmp_path,
            {'rotdpy_mess_mode': 'production'}, 10., {}, 'r')


def test_multistate_leaf_preserves_a_reviewed_reference_root(tmp_path):
    from test_rotdpy_profiles import sulfur_correction
    from kinbot.vrc_correction import resolved_state_model, spin_splitting_settings
    par, _, _ = sulfur_correction(tmp_path)
    stage = validate_surface_model(model())['channels'][0]['bottlenecks'][0]
    stage['parameters']['vrc_tst_selected_roots'] = {'4': 3}
    leaf = state_parameters(par, stage, stage['states'][1], [2, 3])
    settings = spin_splitting_settings(leaf, parent_multiplicity=2, fragment_multiplicities=[2, 3])
    resolved = resolved_state_model(leaf, ['C', 'S', 'H', 'H', 'H'], 0,
        settings, roots={2: 3, 4: 3}, basis='avtz-f12')
    assert resolved['selected_roots'] == {2: 2, 4: 3}


def test_actual_scalar_templates_are_generated_for_each_state_and_coordinate(tmp_path, monkeypatch):
    from test_rotdpy_profiles import sulfur_correction
    from kinbot.pes import create_rotdpy_inputs
    from kinbot.vrc_correction import expected_levels, spin_splitting_settings
    par, data, original = sulfur_correction(tmp_path)
    monkeypatch.chdir(tmp_path)
    channel = validate_surface_model(model())['channels'][0]
    stages = []
    for stage in channel['bottlenecks']:
        states = []
        for state in stage['states']:
            leaf_par = state_parameters(par, stage, state, [2, 3])
            name = f'r__{stage["name"]}_m{state["multiplicity"]}_r{state["root"]}'
            leaf = deepcopy(data)
            leaf['vrc_molecule']['multiplicity'] = state['multiplicity']
            settings = (spin_splitting_settings(leaf_par, parent_multiplicity=state['multiplicity'],
                fragment_multiplicities=[2, 3]) if leaf_par['vrc_tst_correction_model'] == 'spin_splitting' else None)
            roots = {2: 3, 4: 3} if settings else {state['multiplicity']: 3}
            leaf['state_average_roots'] = roots
            leaf['levels'] = expected_levels(leaf_par, symbols=data['vrc_molecule']['symbols'],
                                            charge=0, settings=settings, state_roots=roots)
            from kinbot.vrc_multistate import _OVERRIDES
            leaf['surface_parameters'] = {key: value for key, value in leaf_par.items()
                if key in _OVERRIDES or key in {'vrc_tst_selected_roots', 'mult', 'rotdpy_dist', '_vrc_dividing_coordinate'}}
            path = original.with_name(f'corr_{name}.json')
            path.write_text(json.dumps(leaf))
            states.append({**state, 'reaction': name, 'correction': str(path),
                           'correction_sha256': file_sha256(path)})
        stages.append({**stage, 'states': states})
    master = {'schema': 1, 'kind': 'multistate', 'reaction': 'parent_hom_sci_1_2',
              'combination': 'statistical', 'bottlenecks': stages,
              'surface_model': {'source': 'fixture', 'sha256': 'fixture'}}
    original.write_text(json.dumps(master))
    class SurfaceFixture:
        def __repr__(self):
            return 'Surface()'
    calls = []
    def com(*args, **kwargs):
        calls.append('com')
        return [1], [0], SurfaceFixture()
    monkeypatch.setattr('kinbot.pp_settings.create_surface', com)
    monkeypatch.setattr('kinbot.pp_settings.create_all_surf_for_dist',
        lambda **kwargs: ([[1]], [[0]], [SurfaceFixture()]))
    files = create_rotdpy_inputs(par,
        [['parent', 'parent_hom_sci_1_2', ['methyl', 'sulfur'], 0.]], [], correction_root=tmp_path)
    assert len(files) == 4 and len(set(files)) == 4
    assert len(calls) == 6  # two outer-state runs, three COM surfaces each
    texts = [Path(path).read_text() for path in files]
    for path, text in zip(files, texts):
        compile(text, path, 'exec')
    assert "'selected_roots': {2: 2, 4: 1}" in texts[1]
    assert "'multiplicity': 4" in texts[-1]
    assert "'state_average': [{'multiplicity': 4" in texts[-1]
    assert Path('rotdPy/parent_hom_sci_1_2.multistate.json').is_file()
