"""Chemical comparisons preserve evidence and explicit user/stage choices."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from kinbot import orbital_review as review
from kinbot import molpro_reconnaissance as recon
from kinbot.execution import atomic_json, file_sha256


SOURCE = Path(__file__).parents[1] / 'examples/anl/ch3s_anl1_f12/orbital_reconnaissance.json'


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


def test_correlation_checks_share_production_energy_contract(tmp_path):
    data = recon._load(SOURCE)
    data.update(probe_method='casscf', avas_reference='rhf', correlation_checks=[
        {'method': 'caspt2-f12b(13,8)', 'shift': 0.2, 'ipea': 0.25},
        {'method': 'mrci-f12b+q(13,8)'}])
    deck = recon._deck(data, 'avtz-f12', data['anchors'][0], 'trial')
    assert deck.count('{multi,') == 1
    assert 'cas-ci' not in deck
    assert 'rs2-f12,shift=0.2,ipea=0.25;wf,25,1,1;state,1,3;' in deck
    assert 'mrci-f12,scalf12=1,cabs_singles=1;wf,25,1,3;state,1,3;' in deck
    assert 'KBVRCSEL=energd(2)' in deck
    assert 'reconc12=KBVRCSEL' in deck
    recon.validate_input_blocks(deck)


def _selected(tmp_path, electrons=5, orbitals=4, averages=None):
    wfu = tmp_path / f'space{electrons}.wfu'
    wfu.write_bytes(b'reviewed checkpoint ' + str(electrons).encode())
    return recon.select_record(wfu, tmp_path / f'space{electrons}.json', record='2140.2',
        active_electrons=electrons, active_orbitals=orbitals,
        state_average=averages or [{'multiplicity': 2, 'roots': 2, 'weights': [0.7, 0.3]},
                                   {'multiplicity': 4, 'roots': 1}],
        selected_roots={2: 1, 4: 1}, review_note='Reviewed candidate interaction energies and characters',
        basis='avtz-f12')


def test_reviewed_averages_survive_production_resolution(tmp_path):
    from kinbot.vrc_correction import resolved_state_model
    selection = _selected(tmp_path)
    parameters = {'vrc_tst_active_space': json.loads(selection.read_text())['active_space'],
                  'vrc_tst_state_model_file': str(selection), 'mult': 2}
    settings = {'target_multiplicity': 2, 'reference_multiplicity': 4}
    model = resolved_state_model(parameters, ['C', 'S', 'H', 'H', 'H'], 0,
                                settings, roots={2: 3, 4: 3}, basis='avtz-f12')
    assert model['state_average'][0] == {'multiplicity': 2, 'roots': 2, 'weights': [0.7, 0.3]}
    assert model['state_average'][1]['roots'] == 1
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
    assert 'gthresh,energy=1.d-10,coeff=1.d-10' in deck
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
