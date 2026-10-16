"""Chemical comparisons preserve evidence and explicit user/stage choices."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from kinbot import orbital_review as review
from kinbot import molpro_reconnaissance as recon
from kinbot.execution import atomic_json, file_sha256


SOURCE = Path(__file__).parents[1] / 'examples/anl/ch3s_anl1_f12/orbital_reconnaissance.json'


def test_correlated_order_diagnostic_keeps_method_controls_and_root_labels():
    row = {'state_model': {'selected_roots': {'2': 1}},
           'correlated_energies_hartree': [
               {'method': 'caspt2(3,3)', 'shift': shift, 'multiplicity': 2,
                'root': root, 'energy_hartree': energy}
               for shift, values in ((0.0, ((1, -10.), (2, -10.2))),
                                     (0.3, ((2, -9.9), (1, -10.1))))
               for root, energy in values]}
    original = deepcopy(row)
    first, second = review._state_order_diagnostics(row)
    assert first['roots_by_increasing_correlated_energy'] == [2, 1]
    assert first['selected_reference_root'] == 1
    assert first['selected_root_above_lowest_computed_hartree'] == pytest.approx(.2)
    assert second['roots_by_increasing_correlated_energy'] == [1, 2]
    assert second['selected_root_above_lowest_computed_hartree'] == 0.
    assert row == original


def test_correlated_order_diagnostic_does_not_invent_missing_selected_energy():
    row = {'state_model': {'selected_roots': {3: 2}},
           'correlated_energies_hartree': [
               {'method': 'mrci+q(2,2)', 'multiplicity': 3,
                'root': 1, 'energy_hartree': -1.} ]}
    assert review._state_order_diagnostics(row)[0]['selected_root_above_lowest_computed_hartree'] is None
    assert review._state_order_diagnostics({}) == []
    assert review._state_order_diagnostics({'state_model': None}) == []


@pytest.mark.parametrize('symbol,bond_length,fragment_mult,counts,multiplicities', [
    ('H', .74, 2, (2, 2), [1, 3]),
    ('N', 1.10, 4, (10, 8), [1, 3, 5, 7]),
])
def test_non_sulfur_channels_derive_spaces_and_spin_blocks(
        tmp_path, symbol, bond_length, fragment_mult, counts, multiplicities):
    source = tmp_path / 'recon.json'
    atomic_json(source, {'schema': 1, 'title': 'unrelated_reaction',
        'symbols': [symbol, symbol],
        'geometry_angstrom': [[0., 0., 0.], [0., 0., bond_length]],
        'charge': 0, 'multiplicity': 1, 'probe_active_space': True,
        'bases': [{'name': 'test_basis', 'basis': 'vdz'}],
        'coordinates': [{'name': 'dissociation', 'bond': [1, 2],
            'distances_angstrom': [bond_length, 2., 4., 20.],
            'fragment_states': [{'atom': atom, 'charge': 0,
                'multiplicity': fragment_mult, 'smiles': f'[{symbol}]'}
                for atom in (1, 2)]}]})
    plan = review.plan(source, tmp_path / 'plan.json')
    definition = json.loads(plan.read_text())
    model = definition['models'][0]
    assert (model['active_electrons'], model['active_orbitals']) == counts
    assert model['fragment_context']['multiplicities'] == multiplicities
    assert model['channel'] == 'dissociation'
    definition['variants'] = [definition['variants'][1]]  # native CASSCF preparation
    atomic_json(plan, definition)
    campaign = json.loads(review.prepare(plan, tmp_path / 'campaign').read_text())
    manifest = json.loads(Path(campaign['systems'][0]['manifest']).read_text())
    assert manifest['jobs'][0]['root_states'] == [
        {'multiplicity': mult, 'root': 1} for mult in multiplicities]


def test_general_plan_and_shared_campaign_prepare(tmp_path):
    first = review.plan(SOURCE, tmp_path / 'plan.json')
    definition = json.loads(first.read_text())
    model = definition['models'][0]
    assert model['fragment_context']['multiplicities'] == [2, 4]
    assert model['region'] == 'diagnostic'
    assert 4 <= len(model['anchors']) <= 8
    assert len(definition['available_anchors']) > len(model['anchors'])
    candidate = deepcopy(model)
    candidate.update(name='chemical_candidate', region='inner', active_electrons=5, active_orbitals=4,
                     chemical_note='Radical orbital and complete atomic p shell')
    candidate['state_model']['orbital_selection']['centers'] = [
        {'atom': 1, 'orbitals': ['2pz']}, {'atom': 2, 'orbitals': ['3p']}]
    candidate['state_model']['state_average'] = [{'multiplicity': 2, 'roots': 2}]
    candidate['state_model']['selected_roots'] = {2: 1}
    definition['models'].append(candidate)
    atomic_json(first, definition)
    second = review.plan(SOURCE, tmp_path / 'second.json')
    path = review.prepare(first, tmp_path / 'campaign', additional_plans=[second])
    campaign = json.loads(path.read_text())
    assert len(campaign['systems']) == 9
    for system in campaign['systems']:
        manifest = json.loads(Path(system['manifest']).read_text())
        assert manifest['scan_policy'] is None  # finite panel, no adaptive expansion
        assert Path(manifest['source']).is_file()  # not a deleted temporary config
        assert Path(json.loads(Path(manifest['source']).read_text())['kinbot_input']).is_file()
        for job in manifest['jobs']:
            deck = (Path(job['directory']) / (job['name'] + '.inp')).read_text()
            if not job['orbital_predecessor']:
                assert 'start,2100.2' in deck
            if job['probe_method'] == 'casscf':
                assert 'cas-ci' not in deck
                assert 'natorb,2140.2' in deck
                if job['orbital_predecessor']:
                    assert 'start,2140.2,type=natural' in deck
            if 'chemical_candidate' in system['title']:
                assert job['root_states'] == [{'multiplicity': 2, 'root': 1},
                                               {'multiplicity': 2, 'root': 2}]
    with pytest.raises(FileExistsError):
        review.prepare(first, tmp_path / 'campaign')


def test_plan_does_not_accept_stale_chemistry(tmp_path):
    source = tmp_path / 'source.json'
    data = json.loads(SOURCE.read_text())
    data['kinbot_input'] = str(SOURCE.parent / data['kinbot_input'])
    atomic_json(source, data)
    path = review.plan(source, tmp_path / 'plan.json')
    data['charge'] = 1
    atomic_json(source, data)
    with pytest.raises(ValueError, match='source changed'):
        review.prepare(path, tmp_path / 'campaign')


def test_plan_exposes_requested_rotd_controls_and_unreviewed_stages(tmp_path):
    context = json.loads((SOURCE.parent / 'ch3s.json').read_text())
    context.update(title='arbitrary_title', vrc_tst_sample_basis='avdz-f12',
                   vrc_tst_scan_points=[4., 4.3, 4.6], rotdpy_dist=[6., 7.],
                   rotdpy_flux_parameters={'flux_rel_err': 2.5, 'pot_smp_max': 1500,
                       'pot_smp_min': 100, 'tot_smp_max': 20000, 'tot_smp_min': 100,
                       'smp_len': 1},
                   vrc_tst_surface_model_file='surfaces.json')
    atomic_json(tmp_path / 'context.json', context)
    stages = {'schema': 1, 'reviewed': False, 'channels': [{'bond': [1, 2],
        'bottlenecks': [{'name': 'outer', 'coordinate': 'com',
                        'parameters': {'vrc_tst_sample_basis': 'avtz-f12'}}]}]}
    atomic_json(tmp_path / 'surfaces.json', stages)
    data = json.loads(SOURCE.read_text())
    data['kinbot_input'] = 'context.json'
    atomic_json(tmp_path / 'source.json', data)
    result = json.loads(review.plan(tmp_path / 'source.json', tmp_path / 'plan.json').read_text())
    target = result['rotd_target']
    assert target['parameters']['vrc_tst_sample_basis'] == 'avdz-f12'
    assert target['parameters']['vrc_tst_scan_points'] == [4., 4.3, 4.6]
    assert target['parameters']['rotdpy_dist'] == [6., 7.]
    assert target['parameters']['rotdpy_flux_parameters'] == context['rotdpy_flux_parameters']
    assert target['surface_model']['model'] == stages  # never approves it
    assert target['surface_model']['sha256'] == file_sha256(tmp_path / 'surfaces.json')
    assert {'accuracy', 'continuity', 'chemistry', 'cost', 'kinetics'} <= set(target['acceptance_checks'])


def test_correlation_checks_share_production_energy_contract(tmp_path):
    data = recon._load(SOURCE)
    data.update(probe_method='casscf', avas_reference='rhf', correlation_checks=[
        {'method': 'caspt2-f12b(13,8)', 'shift': 0.2, 'ipea': 0.25},
        {'method': 'mrci-f12b+q(13,8)'}])
    deck = recon._deck(data, 'avtz-f12', data['anchors'][0], 'trial')
    assert deck.count('{multi,') == 1
    assert 'cas-ci' not in deck
    assert 'canonical,2141.2\n' in deck and 'natorb,2140.2\n' in deck
    assert 'rs2-f12,shift=0.2,ipea=0.25,thrden=1.d-8,thrvar=1.d-10;orbital,2141.2,type=canonical;wf,25,1,1;state,1,3;' in deck
    assert 'mrci-f12,scalf12=1,cabs_singles=1,thrden=1.d-8,thrvar=1.d-10;orbital,2140.2,type=natural;wf,25,1,3;state,1,3;' in deck
    assert 'KBVRCSEL=energd(1)' in deck
    assert 'reconc12=KBVRCSEL' in deck
    recon.validate_input_blocks(deck)


def test_shared_reference_checks_and_explicit_high_spin_anchor(tmp_path):
    from kinbot.molpro import _vrc_multireference_method
    data = recon._load(SOURCE)
    data.update(probe_method='casscf', correlation_checks=[
        {'method': 'caspt2-f12(13,8)', 'shift': .3, 'ipea': 0., 'reference_policy': 'state_average'},
        {'method': 'mrci-f12+q(13,8)', 'reference_policy': 'state_average'},
        {'method': 'uccsd(t)-f12b', 'multiplicity': 4}])
    deck = recon._deck(data, 'avtz-f12', data['anchors'][0], 'trial')
    assert 'state,1,3;refstate,3,1,2,3;weight,' in deck
    assert 'state,1,3;refstate,3,1,2,3;option,cluster=0;' in deck
    assert deck.count('{multi,') == 1
    assert 'wf,25,1,3;orbital,2150.2}\n{uccsd(t)-f12}\nreconc13=energy(2)' in deck
    states = [{'multiplicity': m, 'root': r} for m in (2, 4) for r in (1, 2, 3)]
    labels = recon._correlation_labels(data['correlation_checks'], states)
    assert len(labels) == 13 and labels[-1]['multiplicity'] == 4 and labels[-1]['root'] == 1
    model = recon._probe_model(data, data['anchors'][0])
    model['correlation_options'] = {'reference_policy': 'state_average', 'max_reference_iterations': 200}
    model['selected_roots'][2] = 2
    production, _ = _vrc_multireference_method('caspt2-f12(13,8)', 25, 1, state_model=model)
    assert 'state,1,2;refstate,3,1,2,3;weight,' in production
    assert 'maxiti=200' in production
    recon.validate_input_blocks(deck)


def test_cc_anchor_rejects_failed_scf_and_wrong_f12_selector():
    job = {'probe_method': 'casscf', 'correlation_checks': [{'method': 'uccsd(t)-f12b', 'multiplicity': 4}]}
    output = ('PROGRAM * MULTI\n CONVERGENCE REACHED!\n'
              'PROGRAM * Restricted Hartree-Fock\n NO CONVERGENCE!\n'
              'PROGRAM * CCSD\n RHF-UCCSD-F12b energy -20.0\n Triples (T) contribution -0.1\n')
    assert any('Single-reference' in error for error in recon.convergence_diagnostics(output, job)['errors'])
    labels = recon._correlation_labels(job['correlation_checks'], [])
    recon._validate_cc_check_energies(output, labels, {1: -20.1})
    with pytest.raises(ValueError, match='differs'):
        recon._validate_cc_check_energies(output, labels, {1: -20.2})


def test_followup_preparation_keeps_new_basis_fresh_and_uses_basis_specific_seeds(
        tmp_path, partition_source, monkeypatch):
    import importlib.util
    folder = SOURCE.parents[1] / 'orbital_reviews'
    monkeypatch.syspath_prepend(str(folder))
    spec = importlib.util.spec_from_file_location('example_followup', folder / 'prepare_followup.py')
    followup = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(followup)
    row, template = partition_source
    original_campaign = json.loads((Path(row['campaign']) / 'campaign.json').read_text())
    record = original_campaign['systems'][0]
    def all_done(*args, **kwargs):
        raise ValueError('All validation pilots are complete; no new campaign needed.')
    monkeypatch.setattr(review, 'prepare_incomplete', all_done)
    recipe = json.loads((folder / 'followup_recipe.json').read_text())
    case = {'label': 'unrelated_channel', 'spaces': [
        {'electrons': 5, 'orbitals': 4, 'close_lowest': 4, 'radial_neighbors': True},
        {'electrons': 13, 'orbitals': 8, 'close_lowest': 0, 'radial_neighbors': False}]}
    recipe['inner_average_cases'] = [case['label']]
    original = []
    for space in case['spaces']:
        group = f'unrelated_channel_cas{space["electrons"]}_{space["orbitals"]}'
        for point, distance in [('center', 3.), ('far', 30.), ('tilted', 3.), ('left', 2.8), ('right', 3.2)]:
            if point in ('left', 'right') and not space['radial_neighbors']:
                continue
            original.append({**deepcopy(template), 'name': group + '_' + point, 'model': group,
                'target_distance_angstrom': distance, 'close_lowest': space['close_lowest'],
                'active_electrons': space['electrons'], 'active_orbitals': space['orbitals'],
                'correlation_checks': [{'method': f'caspt2-f12({space["electrons"]},{space["orbitals"]})',
                                        'shift': .3, 'ipea': .25}]})
    prepared = followup.prepare_references(tmp_path / 'followup', tmp_path / 'old', original, [case], recipe)
    saved = json.loads(prepared.read_text())
    metadata = {r['title']: r for r in json.loads((prepared.parent / 'review_index.json').read_text())['systems']}
    assert len(saved['systems']) == 14
    basis_rows = []
    for system in saved['systems']:
        loaded = recon._load(system['configuration'])  # All temporary input/context paths survived copying.
        manifest = json.loads(Path(system['pilot_manifest']).read_text())
        job = manifest['jobs'][0]
        if metadata[system['title']].get('purpose') == 'basis_reference':
            assert job['basis'] == 'avdz-f12' and not job.get('initial_record')
            assert '{avas' in Path(job['directory'], job['name'] + '.inp').read_text()
            wfu = Path(job['directory']) / 'completed.wfu'
            wfu.write_bytes(b'explicit synthetic DZ seed')
            job.update(status='complete', phase='pilot', campaign=str(prepared.parent),
                review=metadata[system['title']], electron_count=25, output_sha256='synthetic-output',
                wfu=str(wfu), wfu_sha256=file_sha256(wfu), native_active_space={
                    'active_electrons': 13, 'active_orbitals': 8})
            basis_rows.append(job)
    monkeypatch.setattr(review, '_review_rows', lambda root:
                        (saved, basis_rows) if root == prepared.parent else (original_campaign, [row]))
    second = followup.prepare_potentials(tmp_path / 'potentials', prepared.parent, original, [case], recipe)
    systems = json.loads(second.read_text())['systems']
    assert len(systems) == 20  # 8 DZ + 6 TZ + 6 reduced-average probes.
    for system in systems:
        manifest = json.loads(Path(system['pilot_manifest']).read_text())
        job = manifest['jobs'][0]
        if system['title'].endswith('_dz'):
            assert job['basis'] == 'avdz-f12'
            assert Path(job['initial_record']['path']).is_relative_to(prepared.parent)
        else:
            assert job['basis'] == 'avtz-f12' and job['initial_record']['path'] == row['wfu']


@pytest.fixture
def partition_source(tmp_path, monkeypatch):
    plan_path = review.plan(SOURCE, tmp_path / 'source_plan.json')
    plan = json.loads(plan_path.read_text())
    plan['variants'] = [plan['variants'][1]]
    plan['models'][0]['anchors'] = plan['models'][0]['anchors'][:1]
    atomic_json(plan_path, plan)
    campaign_path = review.prepare(plan_path, tmp_path / 'original')
    campaign = json.loads(campaign_path.read_text())
    manifest = json.loads(Path(campaign['systems'][0]['manifest']).read_text())
    row = manifest['jobs'][0]
    wfu = Path(row['directory']) / 'reviewed.wfu'
    wfu.write_bytes(b'verified source natural orbitals')
    metadata = json.loads((campaign_path.parent / 'review_index.json').read_text())['systems'][0]
    row.update(status='complete', campaign=str(campaign_path.parent), review=metadata,
               wfu=str(wfu), wfu_sha256=file_sha256(wfu), output_sha256='source-output-hash',
               electron_count=manifest['electron_count'],
               native_active_space={key: row['probe_active_space'][key]
                                    for key in ('active_electrons', 'active_orbitals')})
    monkeypatch.setattr(review, '_review_rows', lambda root: (campaign, [row]))
    item = {'name': 'trial_near', 'model': 'trial', 'source_directory': row['directory'],
            'region': 'diagnostic', 'chemical_note': 'Radical and complete sulfur p manifold',
            'active_electrons': 5, 'active_orbitals': 4, 'close_lowest': 4, 'discard_highest': 0,
            'correlation_checks': [{'method': 'caspt2-f12b(5,4)', 'shift': .3, 'ipea': .25}]}
    return row, item


def test_partition_reoptimizes_reviewed_natural_orbitals_without_avas_or_basis_change(tmp_path, partition_source):
    row, item = partition_source
    source_bytes = Path(row['wfu']).read_bytes()
    specification = tmp_path / 'partition.json'
    atomic_json(specification, {'schema': 1, 'models': [item]})
    result = review.prepare_partitions(specification, tmp_path / 'candidate')
    campaign = json.loads(result.read_text())
    for key in ('pilot_manifest', 'manifest'):
        manifest = json.loads(Path(campaign['systems'][0][key]).read_text())
        assert len(manifest['jobs']) == 1
        job = manifest['jobs'][0]
        assert job['basis'] == row['basis'] and job['geometry_angstrom'] == row['geometry_angstrom']
        assert job['root_states'] == row['root_states']
        assert job['initial_record'] == {'path': row['wfu'], 'sha256': row['wfu_sha256']}
        assert job['orbital_predecessor'] is None
        deck = Path(job['directory'], job['name'] + '.inp').read_text()
        assert deck.count('file,2,') == 1 and deck.count('{multi,') == 1
        assert 'occ,14\nclosed,10\n' in deck
        assert 'start,2140.2,type=natural' in deck and '{avas' not in deck and 'diab,' not in deck
        assert 'accuracy,energy=1.d-8' in deck
        assert 'orbital,2141.2,type=canonical' in deck
        assert 'cas-ci' not in deck and job.get('status') != 'complete'
    assert Path(row['wfu']).read_bytes() == source_bytes
    config = json.loads(Path(campaign['systems'][0]['configuration']).read_text())
    config.pop('initial_records')
    atomic_json(tmp_path / 'unseeded.json', config)
    with pytest.raises(ValueError, match='own verified'):
        recon.prepare(tmp_path / 'unseeded.json', tmp_path / 'unseeded')


def test_partition_iteration_controls_reach_every_correlated_root(tmp_path, partition_source):
    _, item = partition_source
    item['casscf_options'] = {'max_iterations': 100}
    item['correlation_options'] = {'max_reference_iterations': 200}
    specification = tmp_path / 'partition.json'
    atomic_json(specification, {'schema': 1, 'models': [item]})
    result = review.prepare_partitions(specification, tmp_path / 'candidate')
    system = json.loads(result.read_text())['systems'][0]
    job = json.loads(Path(system['pilot_manifest']).read_text())['jobs'][0]
    deck = Path(job['directory'], job['name'] + '.inp').read_text()
    assert 'maxit,100;' in deck and 'accuracy,energy=1.d-8' in deck
    assert deck.count('option,maxiti=200;') == len(job['root_states'])
    assert 'rs2-f12,shift=0.3,ipea=0.25,thrden=1.d-8,thrvar=1.d-10;' in deck
    assert job['state_model']['correlation_options'] == item['correlation_options']


def test_partition_orbital_exchange_is_source_relative_and_hash_bound(tmp_path, partition_source):
    row, item = partition_source
    item.update(orbital_swaps=[[2, 6]], source_wfu_sha256=row['wfu_sha256'])
    source = tmp_path / 'swap.json'
    atomic_json(source, {'schema': 1, 'models': [item]})
    result = review.prepare_partitions(source, tmp_path / 'candidate')
    system = json.loads(result.read_text())['systems'][0]
    job = json.loads(Path(system['pilot_manifest']).read_text())['jobs'][0]
    deck = Path(job['directory'], job['name'] + '.inp').read_text()
    assert 'start,2140.2,type=natural\nrotate,8.1,12.1,0' in deck
    assert 'occ,14\nclosed,10' in deck
    item['source_wfu_sha256'] = 'different reviewed checkpoint'
    atomic_json(source, {'schema': 1, 'models': [item]})
    with pytest.raises(ValueError, match='orbital identity'):
        review.prepare_partitions(source, tmp_path / 'changed')


@pytest.mark.parametrize('swaps', [[[1, 1]], [[0, 2]], [[1, 9]], [[True, 2]], [[1]], '1,2'])
def test_partition_rejects_invalid_orbital_exchanges(tmp_path, partition_source, swaps):
    _, item = partition_source
    item['orbital_swaps'] = swaps
    source = tmp_path / 'swap.json'
    atomic_json(source, {'schema': 1, 'models': [item]})
    with pytest.raises(ValueError, match='orbital_swaps'):
        review.prepare_partitions(source, tmp_path / 'candidate')


def test_prepare_incomplete_preserves_complete_points_and_only_raises_limits(tmp_path, partition_source, monkeypatch):
    from types import SimpleNamespace
    _, item = partition_source
    second = deepcopy(item)
    second['name'] = 'other_point'
    specification = tmp_path / 'partition.json'
    atomic_json(specification, {'schema': 1, 'models': [item, second]})
    original = review.prepare_partitions(specification, tmp_path / 'candidate')
    campaign = json.loads(original.read_text())
    originals = {Path(s['configuration']): Path(s['configuration']).read_bytes()
                 for s in campaign['systems']}
    def audit(path):
        manifest = json.loads(Path(path).read_text())
        job = manifest['jobs'][0]
        if Path(path).parent.name == '01_pilot':
            job.update(status='complete', output_sha256='validated-output')
        else:
            job.update(status='incomplete', validation_errors=['Reference CI did not converge'])
        return manifest
    monkeypatch.setattr(recon, '_audit', audit)
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **kw: SimpleNamespace(stdout=''))
    review.main(['prepare-incomplete', str(original.parent), str(tmp_path / 'remaining'),
                 '--casscf-max-iterations', '150', '--reference-max-iterations', '200'])
    fixed = json.loads((tmp_path / 'remaining/campaign.json').read_text())
    assert len(fixed['systems']) == 1 and fixed['systems'][0]['title'] == second['name']
    previous = json.loads(originals[Path(campaign['systems'][1]['configuration'])])
    updated = json.loads(Path(fixed['systems'][0]['configuration']).read_text())
    for key in ('anchors', 'bases', 'scf_options', 'scf_convergence', 'initial_records', 'correlation_checks'):
        assert updated.get(key) == previous.get(key)
    expected = deepcopy(previous['casscf_options'])
    expected['max_iterations'] = 150
    assert updated['casscf_options'] == expected
    previous['probe_definition']['state_model']['correlation_options'] = {'max_reference_iterations': 200}
    assert updated['probe_definition'] == previous['probe_definition']
    assert all(p.read_bytes() == content for p, content in originals.items())
    provenance = json.loads((tmp_path / 'remaining/review_plan.json').read_text())
    assert len(provenance['completed_preserved']) == 1
    job = json.loads(Path(fixed['systems'][0]['pilot_manifest']).read_text())['jobs'][0]
    deck = Path(job['directory'], job['name'] + '.inp').read_text()
    assert 'maxit,150;' in deck and 'option,maxiti=200;' in deck


def test_prepare_incomplete_refuses_multigeometry_dependencies(tmp_path, monkeypatch):
    from types import SimpleNamespace
    path = review.prepare(review.plan(SOURCE, tmp_path / 'plan.json'), tmp_path / 'original')
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **kw: SimpleNamespace(stdout=''))
    with pytest.raises(ValueError, match='single-geometry'):
        review.prepare_incomplete(path.parent, tmp_path / 'remaining')
    assert not (tmp_path / 'remaining').exists()


@pytest.mark.parametrize('options', [None, {'max_iterations': True},
    {'max_reference_iterations': 0}, {'max_reference_iterations': 2.5},
    {'max_reference_vectors': True}, {'max_reference_vectors': 0}, {'ignore_error': True}])
def test_invalid_correlation_controls_cannot_enter_a_selection(tmp_path, options):
    from kinbot.molpro import validate_mr_correlation_options
    with pytest.raises(ValueError, match='correlation_options'):
        validate_mr_correlation_options(options)


@pytest.mark.parametrize('invalid', ['count', 'native_count', 'basis', 'incomplete', 'casci'])
def test_partition_refuses_incompatible_or_unvalidated_sources(tmp_path, partition_source, invalid):
    row, item = partition_source
    if invalid == 'count':
        item['close_lowest'] = 3
    elif invalid == 'native_count':
        row['native_active_space']['active_orbitals'] -= 1
    elif invalid == 'basis':
        item['basis'] = 'avdz-f12'
    elif invalid == 'casci':
        row['probe_method'] = 'casci'
    else:
        row['status'] = 'incomplete'
    specification = tmp_path / 'partition.json'
    atomic_json(specification, {'schema': 1, 'models': [item]})
    with pytest.raises(ValueError):
        review.prepare_partitions(specification, tmp_path / 'candidate')


def test_partition_can_seed_a_new_radial_distance_without_relabeling_source(tmp_path, partition_source):
    import numpy as np
    row, item = partition_source
    original_geometry = deepcopy(row['geometry_angstrom'])
    item['target_distance_angstrom'] = 3.3
    specification = tmp_path / 'partition.json'
    atomic_json(specification, {'schema': 1, 'models': [item]})
    result = review.prepare_partitions(specification, tmp_path / 'candidate')
    system = json.loads(result.read_text())['systems'][0]
    job = json.loads(Path(system['pilot_manifest']).read_text())['jobs'][0]
    xyz = np.array(job['geometry_angstrom'])
    assert np.linalg.norm(xyz[1] - xyz[0]) == pytest.approx(3.3)
    assert job['initial_record']['path'] == row['wfu']
    assert job.get('status') != 'complete' and row['geometry_angstrom'] == original_geometry
    metadata = json.loads((result.parent / 'review_index.json').read_text())['systems'][0]
    assert metadata['source_geometry_angstrom'] == original_geometry
    item['target_distance_angstrom'] = float('nan')
    atomic_json(specification, {'schema': 1, 'models': [item]})
    with pytest.raises(ValueError, match='positive and finite'):
        review.prepare_partitions(specification, tmp_path / 'invalid')


def test_partition_explicit_geometry_is_preserved_with_original_seed_provenance(tmp_path, partition_source):
    row, item = partition_source
    original_geometry = deepcopy(row['geometry_angstrom'])
    geometry = deepcopy(original_geometry)
    geometry[0][0] += .023  # Change an internal coordinate, not only separation.
    geometry[-1][2] -= .017
    item['target_geometry_angstrom'] = geometry
    source = tmp_path / 'fixed_geometry.json'
    atomic_json(source, {'schema': 1, 'models': [item]})
    result = review.prepare_partitions(source, tmp_path / 'candidate')
    system = json.loads(result.read_text())['systems'][0]
    for phase in ('pilot_manifest', 'manifest'):
        job = json.loads(Path(system[phase]).read_text())['jobs'][0]
        assert job['geometry_angstrom'] == geometry
        assert job['initial_record'] == {'path': row['wfu'], 'sha256': row['wfu_sha256']}
        assert job['root_states'] == row['root_states']
        assert job.get('status') != 'complete'
    metadata = json.loads((result.parent / 'review_index.json').read_text())['systems'][0]
    assert metadata['source_geometry_angstrom'] == original_geometry
    assert metadata['target_geometry_angstrom'] == geometry
    assert metadata['channel'] == row['review']['channel']
    assert row['geometry_angstrom'] == original_geometry


@pytest.mark.parametrize('invalid', ['shape', 'count', 'nonfinite', 'boolean', 'conflict'])
def test_partition_refuses_ambiguous_or_invalid_target_geometry(tmp_path, partition_source, invalid):
    row, item = partition_source
    geometry = deepcopy(row['geometry_angstrom'])
    if invalid == 'shape':
        geometry[0].pop()
    elif invalid == 'count':
        geometry.pop()
    elif invalid == 'nonfinite':
        geometry[0][0] = float('nan')
    elif invalid == 'boolean':
        geometry[0][0] = True
    else:
        item['target_distance_angstrom'] = 3.3
    item['target_geometry_angstrom'] = geometry
    source = tmp_path / 'bad_geometry.json'
    atomic_json(source, {'schema': 1, 'models': [item]})
    with pytest.raises(ValueError, match='[Tt]arget (geometry|distance)'):
        review.prepare_partitions(source, tmp_path / 'candidate')
    assert not (tmp_path / 'candidate').exists()


def test_potential_profiles_match_root_labels_and_keep_unlike_curves_separate():
    base = {'review': {'model': 'candidate', 'region': 'inner', 'channel': 'reaction'},
            'basis': 'basis', 'campaign': 'campaign', 'phase': 'pilot', 'symbols': ['H', 'H'],
            'electron_count': 2, 'probe_method': 'casscf', 'state_model': {},
            'probe_active_space': {'active_electrons': 2, 'active_orbitals': 2}}
    rows = []
    for distance, values in [(2., [-10., -9.]), (3., [-9.5, -8.9]), (30., [-9.4, -8.8])]:
        states = [{'multiplicity': 1, 'root': r + 1, 'energy_hartree': e} for r, e in enumerate(values)]
        if distance == 30.:
            states.reverse()  # A list order cannot determine the energy pairing.
        rows.append({**deepcopy(base), 'directory': str(distance), 'distance_angstrom': distance,
                     'geometry_angstrom': [[0, 0, 0], [0, 0, distance]], 'state_energies_hartree': states})
    tilted = deepcopy(rows[0])
    tilted['orientation'] = {'axis': 'bond', 'angle_degrees': 60.}
    profiles = review._potential_profiles(rows + [tilted])
    radial = profiles[0]
    assert radial['points'][0]['states'][0]['relative_to_far_hartree'] == pytest.approx(-.6)
    assert radial['points'][1]['states'][0]['secant_from_previous_hartree_per_angstrom'] == pytest.approx(.5)
    assert profiles[1]['reference_directory'] is None  # no fake zero interaction for a lone orientation
    assert all(p['accuracy_status'].startswith('unreviewed') for p in profiles)
    duplicate = review._potential_profiles(rows + [rows[0]])[0]
    assert duplicate['reference_directory'] is None

    # Additional far-point diagnostics do not split an otherwise matching
    # curve, and a different IPEA/reference treatment cannot become its zero.
    for index, row in enumerate(rows):
        row['correlation_checks'] = [{'method': 'caspt2(2,2)', 'ipea': .25}]
        row['correlated_energies_hartree'] = [{
            **row['correlation_checks'][0], 'multiplicity': 1, 'root': 1,
            'energy_hartree': -11. + index * .2}]
    rows[-1]['correlation_checks'].append({'method': 'caspt2(2,2)', 'ipea': 0.})
    rows[-1]['correlated_energies_hartree'].append({
        **rows[-1]['correlation_checks'][-1], 'multiplicity': 1, 'root': 1,
        'energy_hartree': -12.})
    radial = review._potential_profiles(rows)[0]
    assert len(radial['points']) == 3
    assert radial['points'][0]['states'][-1]['relative_to_far_hartree'] == pytest.approx(-.4)
    assert len(radial['points'][0]['states']) == 3  # No invented IPEA=0 energy.


def _selected(tmp_path, electrons=5, orbitals=4, averages=None, casscf_options=None,
              correlation_options=None):
    wfu = tmp_path / f'space{electrons}.wfu'
    wfu.write_bytes(b'reviewed checkpoint ' + str(electrons).encode())
    return recon.select_record(wfu, tmp_path / f'space{electrons}.json', record='2140.2',
        active_electrons=electrons, active_orbitals=orbitals,
        state_average=averages or [{'multiplicity': 2, 'roots': 2, 'weights': [0.7, 0.3]},
                                   {'multiplicity': 4, 'roots': 1}],
        selected_roots={2: 1, 4: 1}, review_note='Reviewed candidate interaction energies and characters',
        basis='avtz-f12', casscf_options=casscf_options, correlation_options=correlation_options)


def test_reviewed_averages_survive_production_resolution(tmp_path):
    from kinbot.vrc_correction import resolved_state_model
    from kinbot.molpro import stage_vrc_state_model, _vrc_multireference_method
    controls = {'solver': 'so_sci', 'energy_convergence': 1e-8,
                'gradient_convergence': 1e-5, 'quasi_newton': True,
                'exact_diagonalization': True, 'max_iterations': 70}
    correlation = {'max_iterations': 80, 'max_reference_iterations': 200, 'max_reference_vectors': 60}
    selection = _selected(tmp_path, casscf_options=controls, correlation_options=correlation)
    parameters = {'vrc_tst_active_space': json.loads(selection.read_text())['active_space'],
                  'vrc_tst_state_model_file': str(selection), 'mult': 2}
    settings = {'target_multiplicity': 2, 'reference_multiplicity': 4}
    model = resolved_state_model(parameters, ['C', 'S', 'H', 'H', 'H'], 0,
                                settings, roots={2: 3, 4: 3}, basis='avtz-f12')
    assert model['state_average'][0] == {'multiplicity': 2, 'roots': 2, 'weights': [0.7, 0.3]}
    assert model['state_average'][1]['roots'] == 1
    assert model['casscf_options'] == controls
    assert model['correlation_options'] == correlation
    staged = stage_vrc_state_model(model, tmp_path / 'native', 'job')
    for method in ('caspt2-f12b(5,4)', 'mrci-f12b+q(5,4)'):
        deck, _ = _vrc_multireference_method(method, 25, 1, state_model=staged)
        assert '{multi,so-sci,exact_diag=1,qn=1;' in deck
        assert 'accuracy,energy=1.d-8,gradient=1.d-5' in deck
        assert 'maxit,70;' in deck
        assert 'option,maxit=80,maxiti=200,maxvi=60;' in deck
    parameters['vrc_tst_state_average_roots'] = {'2': 3}
    with pytest.raises(ValueError, match='conflict'):
        resolved_state_model(parameters, ['C', 'S', 'H', 'H', 'H'], 0,
                             settings, basis='avtz-f12')
    parameters.pop('vrc_tst_state_average_roots')
    parameters['vrc_tst_selected_roots'] = {'2': 3}
    with pytest.raises(ValueError):
        resolved_state_model(parameters, ['C', 'S', 'H', 'H', 'H'], 0,
                             settings, basis='avtz-f12')


def test_different_inner_outer_spaces_install_without_selecting_physics(tmp_path):
    from test_vrc_multistate import model
    from kinbot.vrc_multistate import load_surface_model, state_parameters
    from kinbot.vrc_correction import resolved_state_model
    source = tmp_path / 'input.json'
    atomic_json(source, {'title': 'workflow', 'mult': 2, 'composite_method': 'ANL1-F12',
        'rotdpy_run': 1, 'me': 1, 'run_me': 1, 'vrc_tst_state_average_roots': {'2': 9},
        'vrc_tst_selected_roots': {'2': 8}})
    surface = model()
    surface['reviewed'] = False  # installing orbitals must not approve kinetic assumptions
    surface_file = tmp_path / 'surfaces.json'
    atomic_json(surface_file, surface)
    first = review.install_selection(source, _selected(tmp_path), tmp_path / 'inner.json',
                                      surface_file=surface_file, bond=[1, 2], stage='inner')
    larger = _selected(tmp_path, 9, 6)
    installed = review.install_selection(first, larger, tmp_path / 'final.json',
        surface_file=tmp_path / 'inner.surfaces.json', bond=[2, 1], stage='outer')
    data = json.loads(installed.read_text())
    assert data['composite_method'] == 'ANL1-F12' and data['run_me'] == 1
    with pytest.raises(ValueError, match='reviewed=true'):
        load_surface_model(data)
    # Emulate the user's separate physical-surface review, then exercise the
    # exact production loader/leaf path with two independently selected spaces.
    final_surface = json.loads(Path(data['vrc_tst_surface_model_file']).read_text())
    final_surface['reviewed'] = True
    atomic_json(Path(data['vrc_tst_surface_model_file']), final_surface)
    loaded, _ = load_surface_model(data)
    for stage, electrons in zip(loaded['channels'][0]['bottlenecks'], [5, 9]):
        leaf = state_parameters(data, stage, stage['states'][0], [2, 3])
        assert leaf['vrc_tst_active_space']['active_electrons'] == electrons
        resolved = resolved_state_model(leaf, ['C', 'S', 'H', 'H', 'H'], 0,
            None, roots={2: 3, 4: 3}, basis='avtz-f12')
        assert resolved['state_average'][0]['weights'] == [0.7, 0.3]
    assert json.loads(source.read_text())['vrc_tst_selected_roots'] == {'2': 8}


def test_report_checks_equivalent_geometry_and_preserves_native_character_data(tmp_path):
    path = review.plan(SOURCE, tmp_path / 'plan.json')
    definition = json.loads(path.read_text())
    definition['variants'] = definition['variants'][1:]  # two optimized initializations
    definition['models'][0]['anchors'] = definition['models'][0]['anchors'][:1]
    atomic_json(path, definition)
    campaign_file = review.prepare(path, tmp_path / 'campaign')
    campaign = json.loads(campaign_file.read_text())
    for index, system in enumerate(campaign['systems']):
        manifest_file = Path(system['manifest'])
        manifest = json.loads(manifest_file.read_text())
        for job in manifest['jobs']:
            directory = Path(job['directory'])
            out = directory / (job['name'] + '.out')
            out.write_text('PROGRAM * MULTI\nCONVERGENCE REACHED!\n'
                           'NATURAL ORBITALS\nMO character and CI coefficients\n')
            wfu = directory / 'test.wfu'
            wfu.write_bytes(b'finished')
            job.update(status='complete', output_sha256=file_sha256(out),
                       wfu=str(wfu), wfu_sha256=file_sha256(wfu),
                       state_energies_hartree=[{'multiplicity': 2, 'root': 1,
                                               'energy_hartree': -400. + index * 0.1}])
        atomic_json(manifest_file, manifest)
    result = json.loads(review.report(campaign_file.parent, tmp_path / 'report.json').read_text())
    assert result['complete'] == 2
    assert result['equivalent_geometry_comparisons'][0]['state_energy_differences_hartree'] == pytest.approx([0.1])
    assert 'CI coefficients' in result['jobs'][0]['native_orbital_and_ci_text']
    assert result['selection_status'] == 'user_review_required'
    wfu.write_bytes(b'changed')
    with pytest.raises(ValueError, match='evidence changed'):
        review.report(campaign_file.parent, tmp_path / 'changed.json')


@pytest.mark.parametrize('different_average', [False, True])
def test_report_combines_repair_pilots_without_accepting_failed_attempts(tmp_path, different_average):
    plan = review.plan(SOURCE, tmp_path / 'plan.json')
    definition = json.loads(plan.read_text())
    definition['variants'] = definition['variants'][1:2]
    definition['models'][0]['anchors'] = definition['models'][0]['anchors'][:1]
    atomic_json(plan, definition)
    roots = [tmp_path / 'original', tmp_path / 'repair']
    for index, root in enumerate(roots):
        campaign = json.loads(review.prepare(plan, root).read_text())
        record = campaign['systems'][0]
        # An earlier full-panel success and a later pilot success are both
        # usable, even though neither campaign has all attempts complete.
        path = Path(record['pilot_manifest' if index else 'manifest'])
        manifest = json.loads(path.read_text())
        job = manifest['jobs'][0]
        directory = Path(job['directory'])
        out, wfu = directory / (job['name'] + '.out'), directory / 'test.wfu'
        out.write_text('PROGRAM * MULTI\nCONVERGENCE REACHED!\nNATURAL ORBITALS\n')
        wfu.write_bytes(b'valid orbitals')
        job.update(status='complete', output_sha256=file_sha256(out), wfu=str(wfu),
                   wfu_sha256=file_sha256(wfu), state_energies_hartree=[
                       {'multiplicity': 2, 'root': 1, 'energy_hartree': -400. + index * 0.01}])
        if index:
            job['initial_record'] = {'path': '/old/casci.wfu', 'sha256': 'old-checkpoint-hash'}
            if different_average:
                job['state_model']['state_average'][0]['weights'] = [1., 0., 0.]
        atomic_json(path, manifest)
        other_path = Path(record['manifest' if index else 'pilot_manifest'])
        other = json.loads(other_path.read_text())
        other['jobs'][0].update(status='incomplete', validation_errors=['Native optimizer did not converge.'])
        atomic_json(other_path, other)
    output = review.report(roots[0], tmp_path / 'combined.json', additional_campaigns=[roots[1], roots[0]])
    result = json.loads(output.read_text())
    assert result['complete'] == 2 and result['total'] == 4
    assert len(result['campaigns']) == 2
    assert result['selection_status'] == 'user_review_required'
    assert len([row for row in result['jobs'] if row.get('validation_errors')]) == 2
    pilot = next(row for row in result['jobs'] if row['phase'] == 'pilot' and row['status'] == 'complete')
    assert pilot['initialization'] == {'kind': 'external_checkpoint',
                                     'path': '/old/casci.wfu', 'sha256': 'old-checkpoint-hash'}
    assert len(result['equivalent_geometry_comparisons']) == (0 if different_average else 1)
    assert {profile['phase'] for profile in result['radial_profiles']} == {'pilot', 'panel'}
    assert all(len(profile['points']) == 1 for profile in result['radial_profiles'])


def test_radial_report_never_subtracts_unlike_correlated_controls(tmp_path, monkeypatch):
    common = {'status': 'complete', 'basis': 'basis', 'probe_method': 'casscf',
        'symbols': ['H', 'H'], 'electron_count': 2, 'campaign': str(tmp_path), 'phase': 'pilot',
        'review': {'title': 'test', 'model': 'test'}, 'channel': 'bond',
        'state_energies_hartree': []}
    selected = {'method': 'caspt2-f12(2,2)', 'ipea': 0., 'multiplicity': 1, 'root': 1}
    shifted = {**selected, 'ipea': .25}
    missing = {**selected, 'shift': .1}
    rows = [{**deepcopy(common), 'name': 'near', 'directory': '/near', 'distance_angstrom': 3.,
             'geometry_angstrom': [[0, 0, 0], [3, 0, 0]],
             'correlated_energies_hartree': [{**selected, 'energy_hartree': -10.1},
                                             {**missing, 'energy_hartree': -10.2}]},
            {**deepcopy(common), 'name': 'far', 'directory': '/far', 'distance_angstrom': 30.,
             'geometry_angstrom': [[0, 0, 0], [30, 0, 0]],
             'correlated_energies_hartree': [{**shifted, 'energy_hartree': -9.},
                                             {**selected, 'energy_hartree': -10.}]}]
    monkeypatch.setattr(review, '_review_rows', lambda _: ({}, rows))
    result = json.loads(review.report(tmp_path, tmp_path / 'report.json').read_text())
    states = result['radial_profiles'][0]['points'][0]['correlated']
    assert states[0]['relative_to_far_hartree'] == pytest.approx(-.1)
    assert 'relative_to_far_hartree' not in states[1]


@pytest.mark.parametrize('solver', ['second_order', 'so_sci'])
def test_convergence_repair_preserves_completed_work_and_chemistry(tmp_path, monkeypatch, solver):
    from types import SimpleNamespace
    plan = review.plan(SOURCE, tmp_path / 'plan.json')
    definition = json.loads(plan.read_text())
    definition['variants'] = definition['variants'][1:]
    definition['models'][0]['anchors'] = definition['models'][0]['anchors'][:2]
    atomic_json(plan, definition)
    original = review.prepare(plan, tmp_path / 'old')
    campaign = json.loads(original.read_text())
    before = {}
    for i, system in enumerate(campaign['systems']):
        for key in ('pilot_manifest', 'manifest'):
            path = Path(system[key])
            manifest = json.loads(path.read_text())
            for j, job in enumerate(manifest['jobs']):
                text = ('NO CONVERGENCE! MAXIMAL NUMBER OF ITERATIONS REACHED!'
                        if i == 0 and j == 0 else 'CONVERGENCE REACHED!')
                Path(job['directory'], job['name'] + '.out').write_text('PROGRAM * MULTI\n' + text + '\n')
                before[Path(job['directory'], job['name'] + '.inp')] = job['input_sha256']
            before[path] = file_sha256(path)
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **kw: SimpleNamespace(stdout=''))
    repaired = review.repair_convergence(original.parent, tmp_path / 'repair', solver=solver)
    records = json.loads(repaired.read_text())['systems']
    assert len(records) == 1
    manifest = json.loads(Path(records[0]['manifest']).read_text())
    assert len(manifest['jobs']) == 1
    job = manifest['jobs'][0]
    deck = Path(job['directory'], job['name'] + '.inp').read_text()
    hf = 'so' if solver == 'second_order' else 'so-sci'
    multi = 'exact_diag=1,qn=0' if solver == 'second_order' else 'so-sci,exact_diag=1,qn=1'
    assert '{rhf,' + hf + ',maxit=100;' in deck
    assert '{multi,' + multi + ';' in deck
    assert 'gthresh,energy=1.d-10,coeff=1.d-5' in deck
    assert 'accuracy,energy=1.d-10\n' in deck
    assert job['root_states'] == json.loads(Path(campaign['systems'][0]['manifest']).read_text())['jobs'][0]['root_states']
    assert all(file_sha256(path) == digest for path, digest in before.items())
    data = recon._load(records[0]['configuration'])
    seeded = recon._deck(data, 'avtz-f12', data['anchors'][0], 'trial', previous='previous')
    assert '{multi,' + multi + ';\nstart,2140.2,type=natural\ndiab,2140.2\n' in seeded
    recon.validate_input_blocks(seeded)
    # Repeating a failed repair must not spend more jobs with identical controls.
    for key in ('pilot_manifest', 'manifest'):
        repeated = json.loads(Path(records[0][key]).read_text())
        for job in repeated['jobs']:
            Path(job['directory'], job['name'] + '.out').write_text(
                'PROGRAM * MULTI\nNO CONVERGENCE! MAXIMAL NUMBER OF ITERATIONS REACHED!\n')
    with pytest.raises(ValueError, match='repeats the failed solver'):
        review.repair_convergence(repaired.parent, tmp_path / 'duplicate', solver=solver)
    assert not (tmp_path / 'duplicate').exists()


def test_so_sci_repair_can_exclude_optional_casci_controls(tmp_path, monkeypatch):
    from types import SimpleNamespace
    plan = review.plan(SOURCE, tmp_path / 'plan.json')
    original = review.prepare(plan, tmp_path / 'original')
    for system in json.loads(original.read_text())['systems']:
        manifest = json.loads(Path(system['manifest']).read_text())
        job = manifest['jobs'][0]
        Path(job['directory'], job['name'] + '.out').write_text(
            'No convergence in rhfpr\nPROGRAM * MULTI\nNO CONVERGENCE!\n')
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **kw: SimpleNamespace(stdout=''))
    review.main(['repair-convergence', str(original.parent), str(tmp_path / 'repair'),
                 '--solver', 'so_sci', '--probe-method', 'casscf'])
    campaign = json.loads((tmp_path / 'repair' / 'campaign.json').read_text())
    assert len(campaign['systems']) == 2
    for system in campaign['systems']:
        config = recon._load(system['configuration'])
        assert config['probe_method'] == 'casscf'
        assert config['scf_options']['solver'] == config['casscf_options']['solver'] == 'so_sci'
        assert len(config['anchors']) == 1


@pytest.mark.parametrize('explicit_override', [False, True])
def test_so_sci_repair_can_change_only_local_accuracy_once(tmp_path, monkeypatch, explicit_override):
    from types import SimpleNamespace
    from kinbot import orbital_campaign
    source = tmp_path / 'input.json'
    data = json.loads(SOURCE.read_text())
    data.pop('kinbot_input', None)
    data.update(probe_method='casscf', scf_options={'solver': 'so_sci', 'max_iterations': 100},
                casscf_options={'solver': 'so_sci', 'exact_diagonalization': True, 'quasi_newton': True})
    if explicit_override:
        data['casscf_options'].update(energy_convergence=1.e-10, gradient_convergence=1.e-5)
    atomic_json(source, data)
    loaded = recon._load(source)
    data.pop('coordinates', None)
    data['anchors'] = loaded['anchors']
    data['probe_definition'] = {
        'active_electrons': loaded['_probe_space']['active_electrons'],
        'active_orbitals': loaded['_probe_space']['active_orbitals'],
        'state_model': recon._probe_model(loaded, loaded['anchors'][0]),
        'chemical_note': 'Valence reference for a threshold repair test.'}
    atomic_json(source, data)
    renderer = recon._deck
    with monkeypatch.context() as old:
        if not explicit_override:
            old.setattr(recon, '_deck', lambda *args, **kw:
                        renderer(*args, **kw).replace('accuracy,energy=1.d-10\n', ''))
        original = orbital_campaign.prepare(tmp_path / 'old', [source])
    campaign = json.loads(original.read_text())
    system = campaign['systems'][0]
    atomic_json(original.parent / 'review_index.json', {'systems': [{'title': system['title']}]})
    pilot = json.loads(Path(system['pilot_manifest']).read_text())['jobs'][0]
    inp = Path(pilot['directory'], pilot['name'] + '.inp')
    before = inp.read_bytes()
    Path(pilot['directory'], pilot['name'] + '.out').write_text(
        'PROGRAM * MULTI\nConvergence thresholds 0.10E-04 (gradient) 0.10E-11 (energy)\n'
        'NO CONVERGENCE! MAXIMAL NUMBER OF ITERATIONS REACHED!\n')
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **kw: SimpleNamespace(stdout=''))
    extra = ['--energy-convergence', '1e-8', '--verbosity', '4'] if explicit_override else []
    review.main(['repair-convergence', str(original.parent), str(tmp_path / 'fixed'),
                 '--solver', 'so_sci', '--probe-method', 'casscf', *extra])
    repaired = tmp_path / 'fixed' / 'campaign.json'
    record = json.loads(repaired.read_text())['systems'][0]
    config = json.loads(Path(record['configuration']).read_text())
    assert config['casscf_options']['energy_convergence'] == (1.e-8 if explicit_override else 1.e-10)
    assert config['scf_options'] == data['scf_options']
    assert config['probe_definition'] == json.loads(Path(system['configuration']).read_text())['probe_definition']
    job = json.loads(Path(record['pilot_manifest']).read_text())['jobs'][0]
    deck = Path(job['directory'], job['name'] + '.inp').read_text()
    expected = 'accuracy,energy=1.d-8,gradient=1.d-5' if explicit_override else 'accuracy,energy=1.d-10'
    assert expected + '\n' in deck
    assert 'gthresh,energy=1.d-10,coeff=1.d-5' in deck
    if explicit_override:
        assert 'so-sci,exact_diag=1,qn=1,verbosity=4;' in deck
        Path(job['directory'], job['name'] + '.out').write_text(
            'PROGRAM * MULTI\nNO CONVERGENCE! MAXIMAL NUMBER OF ITERATIONS REACHED!\n')
        with pytest.raises(ValueError, match='repeats the failed solver'):
            review.repair_convergence(repaired.parent, tmp_path / 'identical',
                                     solver='so_sci', verbosity=3)
    assert inp.read_bytes() == before


def test_reuses_checked_old_wfu_as_casscf_seed_without_relabeling_results(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot import orbital_campaign
    old_campaign = orbital_campaign.prepare(tmp_path / 'old', [SOURCE])
    old = json.loads(old_campaign.read_text())
    old_manifest = Path(old['systems'][0]['manifest'])
    manifest = json.loads(old_manifest.read_text())
    prior = manifest['jobs'][0]
    wfu = Path(prior['directory']) / 'old.wfu'
    wfu.write_bytes(b'original CASCI checkpoint')
    prior.update(status='complete', wfu=str(wfu), wfu_sha256=file_sha256(wfu))
    atomic_json(old_manifest, manifest)
    before = {p: p.read_bytes() for p in (old_manifest, wfu)}
    plan = review.plan(SOURCE, tmp_path / 'plan.json')
    path = review.prepare(plan, tmp_path / 'new', reuse_from=old_campaign.parent)
    campaign = json.loads(path.read_text())
    systems = [s for s in campaign['systems'] if s['title'].endswith('casscf_propagated')]
    prepared = Path(systems[0]['manifest'])
    jobs = json.loads(prepared.read_text())['jobs']
    reused = next(j for j in jobs if j['initial_record'])
    assert reused['initial_record']['sha256'] == file_sha256(wfu)
    assert reused['orbital_predecessor'] is None
    deck = Path(reused['directory'], reused['name'] + '.inp').read_text()
    assert '{{PREVIOUS_WFU_SHA256}}' not in deck
    assert 'cas-ci' not in deck and 'start,2140.2,type=natural' in deck
    assert reused.get('status') != 'complete'
    monkeypatch.setattr(recon, 'qc_runtime_environment', lambda *a, **k: ({}, {'scratch_directory': str(tmp_path)}))
    monkeypatch.setattr(recon, 'cleanup_qc_runtime', lambda *a: None)
    def native(command, **kwargs):
        assert Path(kwargs['cwd'], reused['name'] + '.wfu').read_bytes() == wfu.read_bytes()
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', native)
    recon.execute(prepared, reused['name'])
    assert all(p.read_bytes() == content for p, content in before.items())
    wfu.write_bytes(b'changed')
    with pytest.raises(ValueError, match='changed before execution'):
        recon.execute(prepared, reused['name'])
