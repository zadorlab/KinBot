"""Theory-profile integration checks without calculations or queue submission."""

import ast
import copy
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from ase.build import molecule

from kinbot.irc import IRC
from kinbot.parameters import Parameters
from kinbot.profiled_qc import ProfiledQuantumChemistry, MissingL1Result
from kinbot.qc import QuantumChemistry
from kinbot.reac_family import carry_out_reaction
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint
from kinbot.theory import resolve_profiles


def fairchem(model, device='cpu'):
    return dict(calculator='fairchem', model_path=model, task_name='omol',
                device=device, optimizer='sella', frequency_mode='ase_forces')


def gaussian(method='M062X', command='/l2/g16'):
    return dict(calculator='gaussian', method=method, basis='6-31G(d)',
                command=command, optimizer='native',
                frequency_mode='native_hessian',
                calculator_kwargs={'integral': 'UltraFine'})


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('conf').mkdir()

    def submit(backend, job, nproc, **kwargs):
        backend.job_ids[job] = '314'
        return 1

    monkeypatch.setattr(QuantumChemistry, 'submit_qc', submit)
    return tmp_path


def parameters(**options):
    Path('input.json').write_text(json.dumps(dict(
        barrier_threshold=100., qc='gauss', qc_command='/global/g09',
        method='HF', basis='sto-3g', fc_model_path='/models/global.pt',
        fc_device='cpu', ppn=4, queue_job_limit=5, **options)))
    return Parameters('input.json', show_warnings=False).par


def profile_parameters(**options):
    values = dict(l0_profile=fairchem('/models/small.pt'),
                  l1_profile=fairchem('/models/medium.pt', 'cuda'),
                  l2_profile=gaussian())
    values.update(options)
    return parameters(**values)


def methanol(wellorts=0, name='methanol'):
    atoms = molecule('CH3OH')
    point = StationaryPoint(name, 0, 1, wellorts=wellorts,
                            atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    point.characterize()
    return point


def script(job):
    text = Path(job + '.py').read_text()
    ast.parse(text)
    return text


def literal(text, name):
    for node in ast.walk(ast.parse(text)):
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id == name):
            return ast.literal_eval(node.value)
    raise AssertionError(f'Missing literal assignment: {name}')


def reaction(qc, par, name='methanol_ts'):
    point = methanol(wellorts=1, name=name)
    return SimpleNamespace(
        qc=qc, par=par, species=point, instance_name=name, max_step=1,
        scan=False, skip=False, instance=[0, 1, 2],
        get_constraints=lambda step, geom: (step, [], [], []))


def assert_model(text, model, device):
    assert model in text
    assert device in text
    assert '/models/global.pt' not in text


def test_ordinary_am1_prescreen_keeps_requested_g09(workspace):
    par = parameters(semi_emp_conformer_search=1, semi_emp_method='am1')
    assert resolve_profiles(par) == {}
    qc = QuantumChemistry(par)
    assert type(qc) is QuantumChemistry
    point = methanol()
    qc.qc_conf(point, point.geom, 0, semi_emp=1)
    text = script(f'conf/{routing_name(point)}_semi_emp_0000')
    assert literal(text, 'kwargs')['method'].lower() == 'am1'
    assert '/global/g09' in text
    assert '/l2/g16' not in text
    assert not Path('.kinbot_theory_jobs.json').exists()


def test_l0_legacy_alias_is_fallback_beneath_explicit_profile(workspace):
    par = profile_parameters(semi_emp_method='pm6',
                             l0_profile=gaussian('AM1', '/l0/g09'))
    profiles = resolve_profiles(par)
    assert par['l0_method'] == 'pm6'
    assert profiles['l0'].method == 'AM1'
    assert profiles['l0'].command == '/l0/g09'
    assert profiles['l2'].command == '/l2/g16'


def test_calculator_aliases_and_explicit_precedence(workspace):
    par = parameters(l0='fc', l1='uma', l2='gauss',
                     l0_profile={'model_path': '/models/small.pt'},
                     l1_profile={'model_path': '/models/medium.pt'},
                     l2_profile=gaussian(), calc_kwargs={'global_only': 'yes'})
    profiles = resolve_profiles(par)
    assert [profiles[level].calculator for level in ('l0', 'l1', 'l2')] == [
        'fairchem', 'fairchem', 'gaussian']
    assert 'global_only' not in profiles['l0'].calculator_kwargs
    assert 'global_only' not in profiles['l1'].calculator_kwargs
    assert profiles['l2'].method == 'M062X'


def test_each_model_and_device_stays_in_its_backend(workspace):
    par = profile_parameters(l2_profile=fairchem('/models/large.pt', 'cuda:1'))
    original = copy.deepcopy(par)
    qc = QuantumChemistry(par)
    assert isinstance(qc, ProfiledQuantumChemistry)
    for level, model, device in (
            ('l0', '/models/small.pt', 'cpu'),
            ('l1', '/models/medium.pt', 'cuda'),
            ('l2', '/models/large.pt', 'cuda:1')):
        backend = qc.for_level(level)
        assert backend.par['fc_model_path'] == model
        assert backend.par['fc_device'] == device
        assert backend.par['fc_task_name'] == 'omol'
        assert backend.ppn == 4 and backend.queue_job_limit == 5
    assert par == original


def test_explicit_l2_overrides_legacy_theory_and_merges_keywords(workspace):
    par = parameters(l1_profile=fairchem('/models/medium.pt'),
                     calc_kwargs={'integral': 'UltraFine', 'scf': 'xqc'},
                     l2_profile={'method': 'M062X', 'basis': '6-31G(d)',
                                 'calculator_kwargs': {'scf': 'tight'}})
    profiles = resolve_profiles(par)
    assert profiles['l1'].calculator == 'fairchem'
    assert profiles['l2'].method == 'M062X'
    assert profiles['l2'].basis == '6-31G(d)'
    assert profiles['l2'].calculator_kwargs['integral'] == 'UltraFine'
    assert profiles['l2'].calculator_kwargs['scf'] == 'tight'


@pytest.mark.parametrize('preset', ['uma-b2plyp-anl', 'uma-b3lyp-anl-low'])
@pytest.mark.parametrize('explicit', [False, True])
def test_named_presets_fail_before_creating_calculations(workspace, preset, explicit):
    options = {'l1_profile': fairchem('/models/medium.pt')} if explicit else {}
    with pytest.raises(ValueError, match='theory_preset.*not supported.*l1_profile'):
        parameters(theory_preset=preset, **options)
    assert not Path('.kinbot_theory_jobs.json').exists()
    assert not Path('kinbot.db').exists()


@pytest.mark.parametrize('option, value', [
    ('composite_method', 'ANL0'), ('l3_overrides', {'method': 'ANL1'}),
    ('l3_resource_overrides', {'ppn': 8}),
])
def test_composite_recipe_inputs_remain_unsupported(workspace, option, value):
    with pytest.raises(IOError, match=f'does not recognize option {option}'):
        parameters(**{option: value})


@pytest.mark.parametrize('explicit', [False, True])
def test_saved_inputs_with_empty_preset_keep_their_calculators_and_records(workspace, explicit):
    par = profile_parameters() if explicit else parameters()
    qc = QuantumChemistry(par)
    point = methanol()
    qc.qc_opt(point, point.geom)
    qc.db.write(molecule('CH3OH'), name=routing_name(point) + '_well',
                data=dict(energy=-100., status='normal'))
    # PES writes all parameter defaults, including the former preset field.
    assert par['theory_preset'] == ''
    Path('saved_input.json').write_text(json.dumps(par))
    before = file_contents(workspace)
    loaded = Parameters('saved_input.json', show_warnings=False).par
    restarted = QuantumChemistry(loaded)
    assert resolve_profiles(loaded) == resolve_profiles(par)
    assert type(restarted) is type(qc)
    assert restarted.for_level('l1').qc == qc.for_level('l1').qc
    assert file_contents(workspace) == before


@pytest.mark.parametrize('option', ['ppn', 'queue_job_limit', 'fc_compile'])
def test_per_level_execution_options_are_not_silently_accepted(workspace, option):
    profile = fairchem('/models/small.pt')
    profile[option] = 8
    with pytest.raises(ValueError, match='Unknown l0_profile options'):
        profile_parameters(l0_profile=profile)


@pytest.mark.parametrize('enabled', [False, True])
@pytest.mark.parametrize('profiled', [False, True])
def test_global_compile_option_reaches_every_fairchem_worker(workspace, enabled, profiled):
    if profiled:
        par = profile_parameters(fc_compile=enabled,
                                 l2_profile=fairchem('/models/large.pt'))
    else:
        par = dict(parameters(fc_compile=enabled), qc='fc')
    qc = QuantumChemistry(par)
    point = methanol()
    qc.qc_opt(point, point.geom)
    qc.qc_opt(point, point.geom, high_level=1)
    qc.qc_conf(point, point.geom, 0, semi_emp=1)
    qc.qc_conf(point, point.geom, 1)
    qc.qc_ring_conf(point, point.geom, [[2, 1, 0, 3]], [60.], 0, 0)
    Path('hir').mkdir()
    qc.qc_hir(point, point.geom, 0, 0,
              [[int(atom) + 1 for atom in point.dihed[0]]], False)

    def check(text):
        calls = [node for node in ast.walk(ast.parse(text))
                 if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                 and node.func.id == 'load_predictor']
        assert len(calls) == 1
        options = {arg.arg: ast.literal_eval(arg.value) for arg in calls[0].keywords}
        assert options['compile_model'] is enabled

    rxn = reaction(qc, par)
    backend = qc.for_level('l1')
    backend.check_qc = Mock(return_value='normal')
    backend.is_in_database = Mock(return_value=False)
    backend.get_qc_geom = Mock(return_value=(0, rxn.species.geom))
    for step in (0, 1):
        carry_out_reaction(rxn, step, 'unused')
        check(script(rxn.instance_name))
    IRC(rxn, par).do_irc_calculations()
    for path in workspace.rglob('*.py'):
        check(path.read_text())

    if profiled:
        restarted = QuantumChemistry(dict(par, fc_compile=not enabled))
        assert restarted._fingerprint == qc._fingerprint
        assert restarted._routes == qc._routes
        for level in ('l0', 'l1', 'l2'):
            assert restarted.for_level(level).par['fc_compile'] is not enabled


def test_well_and_refined_well_use_their_own_uma_models(workspace):
    qc = QuantumChemistry(profile_parameters(
        l2_profile=fairchem('/models/large.pt', 'cuda:1')))
    point = methanol()
    qc.qc_opt(point, point.geom)
    qc.qc_opt(point, point.geom, high_level=1)
    name = routing_name(point)
    assert_model(script(name + '_well'), '/models/medium.pt', 'cuda')
    assert_model(script(name + '_well_high'), '/models/large.pt', 'cuda:1')
    assert qc.for_job(name + '_well') is qc.l1
    assert qc.for_job(name + '_well_high') is qc.l2


def test_prescreen_ring_and_l1_conformer_workers_use_assigned_models(workspace):
    qc = QuantumChemistry(profile_parameters())
    point = methanol()
    qc.qc_conf(point, point.geom, 0, semi_emp=1)
    qc.qc_conf(point, point.geom, 1)
    qc.qc_ring_conf(point, point.geom, [[2, 1, 0, 3]], [60.], 0, 0)
    name = 'conf/' + routing_name(point)
    for job in (name + '_semi_emp_0000', name + '_r0000_0000'):
        assert_model(script(job), '/models/small.pt', 'cpu')
        assert qc.for_job(job) is qc.l0
    assert_model(script(name + '_0001'), '/models/medium.pt', 'cuda')
    assert qc.for_job(name + '_0001') is qc.l1


def test_search_and_final_ts_workers_use_l1_model(workspace):
    par = profile_parameters()
    qc = QuantumChemistry(par)
    rxn = reaction(qc, par)
    qc.l1.check_qc = Mock(return_value='normal')
    qc.l1.is_in_database = Mock(return_value=False)
    qc.l1.get_qc_geom = Mock(return_value=(0, rxn.species.geom))
    for step in (0, 1):
        carry_out_reaction(rxn, step, '/incorrect/global/command')
        assert_model(script(rxn.instance_name), '/models/medium.pt', 'cuda')
        assert qc.for_job(rxn.instance_name) is qc.l1


def test_gaussian_ts_search_and_end_keep_profile_method_command_and_keywords(workspace):
    par = profile_parameters(l1_profile=gaussian('B3LYP', '/l1/g16'))
    qc = QuantumChemistry(par)
    rxn = reaction(qc, par)
    qc.l1.check_qc = Mock(return_value='normal')
    qc.l1.is_in_database = Mock(return_value=False)
    qc.l1.get_qc_geom = Mock(return_value=(0, rxn.species.geom))
    for step in (0, 1):
        carry_out_reaction(rxn, step, '/incorrect/global/command')
        text = script(rxn.instance_name)
        kwargs = literal(text, 'kwargs')
        assert kwargs['method'] == ('am1' if step == 0 else 'B3LYP')
        assert kwargs['basis'] == ('' if step == 0 else '6-31G(d)')
        assert kwargs['integral'] == 'UltraFine'
        assert '/l1/g16' in text
        assert '/incorrect/global/command' not in text


def test_ordinary_irc_uses_l1_model_for_path_and_endpoint(workspace):
    par = profile_parameters()
    qc = QuantumChemistry(par)
    rxn = reaction(qc, par)
    qc.l1.get_qc_geom = Mock(return_value=(0, rxn.species.geom))
    IRC(rxn, par).do_irc_calculations()
    for direction in ('F', 'R'):
        job = rxn.instance_name + '_IRC_' + direction
        assert_model(script(job), '/models/medium.pt', 'cuda')
        assert qc.for_job(job) is qc.l1
        assert qc.for_job(job + '_prod') is qc.l1






def test_restart_restores_backend_and_job_id_for_nonstandard_job_name(workspace):
    par = profile_parameters()
    qc = QuantumChemistry(par)
    qc.l0.submit_qc('candidate', 4)
    qc.l2.submit_qc('accepted_IRC_F', 4)
    restarted = QuantumChemistry(par)
    for name, level in (('candidate', 'l0'), ('accepted_IRC_F', 'l2'),
                        ('accepted_IRC_F_prod', 'l2')):
        assert restarted.for_job(name) is restarted.for_level(level)
    assert restarted.l0.job_ids['candidate'] == '314'
    assert restarted.l2.job_ids['accepted_IRC_F'] == '314'
    restarted.l2.get_qc_energy = Mock(return_value=(0, -115.))
    restarted.l1.get_qc_energy = Mock(side_effect=AssertionError('wrong surface'))
    assert restarted.get_qc_energy('accepted_IRC_F_prod') == (0, -115.)


def test_restart_refuses_to_mix_results_after_model_change(workspace):
    par = profile_parameters()
    QuantumChemistry(par).l1.submit_qc('methanol_well', 4)
    changed = copy.deepcopy(par)
    changed['l1_profile']['model_path'] = '/models/different.pt'
    before = Path('.kinbot_theory_jobs.json').read_bytes()
    with pytest.raises(ValueError, match='fresh directory') as error:
        QuantumChemistry(changed)
    assert 'Keep .kinbot_theory_jobs.json with the existing results' in str(error.value)
    assert 'archive or remove' not in str(error.value)
    assert Path('.kinbot_theory_jobs.json').read_bytes() == before


@pytest.mark.parametrize('calculator, option, original, replacement', [
    ('orca', 'orcasimpleinput', 'D3BJ', 'D4'),
    ('orca', 'orcablocks', '%method Grid 4 end', '%method Grid 5 end'),
    ('gaussian', 'integral', 'UltraFine', 'FineGrid'),
])
def test_restart_checks_backend_settings_before_reusing_methanol(
        workspace, calculator, option, original, replacement):
    profile = dict(calculator=calculator, method='B3LYP', basis='6-31G(d)',
                   optimizer='sella', frequency_mode='ase_forces')
    high_level = int(calculator == 'gaussian')
    level = 'l2' if high_level else 'l1'
    par = profile_parameters(high_level=high_level, conformer_search=1,
        queuing='local', **{f'{level}_profile': profile, option: original})
    qc = QuantumChemistry(par)
    point = methanol()
    qc.qc_opt(point, point.geom, high_level=high_level)
    job = routing_name(point) + ('_well_high' if high_level else '_well')
    used = literal(script(job), 'kwargs')
    assert original in used[option]
    qc.db.write(molecule('CH3OH'), name=job, data=dict(energy=-100., status='normal'))
    Path(job + '_sella.log').write_text('done\n')
    before = file_contents(workspace)
    restarted = QuantumChemistry(par)
    assert restarted.get_qc_energy(job)[0] == 0
    assert file_contents(workspace) == before

    changed = dict(par, **{option: replacement})
    with pytest.raises(ValueError, match='fresh directory'):
        QuantumChemistry(changed)
    assert file_contents(workspace) == before


@pytest.mark.parametrize('calculator, level, option, replacement, arguments, keyword', [
    ('gaussian', 'l1', 'guessmix', 1, {}, 'guess'),
    ('gaussian', 'l1', 'scan_method', 'B3LYP', {'scan': 1}, 'method'),
    ('qchem', 'l1', 'scan_basis', '6-31G(d)', {'scan': 1}, 'basis'),
    ('gaussian', 'l2', 'vrc_tst_scan_method', 'HF', {'vts': 1}, 'method'),
    ('gaussian', 'l2', 'vrc_tst_scan_basis', 'STO-3G', {'vts': 1}, 'basis'),
])
def test_restart_checks_method_choices_outside_profiles(
        workspace, calculator, level, option, replacement, arguments, keyword):
    from kinbot.profiled_qc import _backend_parameters

    par = parameters(l0_profile=fairchem('/models/small.pt'), queuing='local')
    par['qc'] = calculator
    qc = QuantumChemistry(par)
    args = ('methanol', 1, 0, 18)
    original = qc.for_level(level).get_qc_arguments(*args, **arguments)
    changed = dict(par, **{option: replacement})
    # Confirm this setting changes the actual requested calculation, not just
    # the recorded input. Construct only a backend in this temporary test.
    backend = QuantumChemistry(_backend_parameters(changed, resolve_profiles(changed)[level]),
                               profile_backend=True)
    requested = backend.get_qc_arguments(*args, **arguments)
    assert requested.get(keyword) != original.get(keyword)
    before = file_contents(workspace)
    with pytest.raises(ValueError, match='fresh directory'):
        QuantumChemistry(changed)
    assert file_contents(workspace) == before


@pytest.mark.parametrize('calculator', ['gaussian', 'orca'])
def test_old_records_without_backend_settings_refuse_restart(workspace, calculator):
    import hashlib

    par = profile_parameters(l1_profile=dict(calculator=calculator, optimizer='sella'),
                             l2_profile=fairchem('/models/large.pt'))
    QuantumChemistry(par)
    path = Path('.kinbot_theory_jobs.json')
    record = json.loads(path.read_text())
    for profile in record['profiles'].values():
        profile.pop('calculation_settings', None)
    record['profile_sha256'] = hashlib.sha256(json.dumps(
        record['profiles'], sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    path.write_text(json.dumps(record))
    before = file_contents(workspace)
    with pytest.raises(ValueError, match='saved theory record lacks required settings'):
        QuantumChemistry(par)
    assert file_contents(workspace) == before


def test_unused_calculator_settings_do_not_change_uma_records(workspace):
    par = profile_parameters(l2_profile=fairchem('/models/large.pt'))
    qc = QuantumChemistry(par)
    before = file_contents(workspace)
    changed = dict(par, orcasimpleinput='D4', orcablocks='%method Grid 5 end',
                   integral='FineGrid', guessmix=1)
    assert QuantumChemistry(changed)._fingerprint == qc._fingerprint
    assert file_contents(workspace) == before


@pytest.mark.parametrize('explicit', [{'l1': 'gauss'},
    {'l1_profile': {'command': '/new/g16'}}, {'profiled_theory': 1}])
def test_explicit_profile_keeps_existing_reaction_search_defaults(workspace, explicit):
    par = parameters(l0_profile=fairchem('/models/small.pt'), queuing='local')
    qc = QuantumChemistry(par)
    changed = QuantumChemistry(dict(par, **explicit))
    for job, options in (
            ('methanol_ts', dict(ts=1, step=0, max_step=10)),
            ('methanol_ts', dict(ts=1, step=10, max_step=10)),
            ('methanol_scan', dict(scan=1))):
        before = qc.l1.get_qc_arguments(job, 1, 0, 18, **options)
        after = changed.l1.get_qc_arguments(job, 1, 0, 18, **options)
        assert (before['method'], before['basis']) == (after['method'], after['basis'])
    assert changed._fingerprint == qc._fingerprint


def test_restart_allows_executable_and_device_changes(workspace):
    par = profile_parameters()
    qc = QuantumChemistry(par)
    qc.l1.submit_qc('methanol_well', 4)
    changed = copy.deepcopy(par)
    changed['l1_profile']['device'] = 'cpu'
    changed['l2_profile']['command'] = '/other/g16'
    restarted = QuantumChemistry(changed)
    assert restarted._routes == qc._routes
    assert restarted.l1.par['fc_device'] == 'cpu'
    assert restarted.l2.qc_command == '/other/g16'
    assert restarted._fingerprint == qc._fingerprint
    for field, value in [('method', 'HF'), ('optimizer', 'sella'),
                         ('frequency_mode', 'auto')]:
        invalid = copy.deepcopy(changed)
        invalid['l2_profile'][field] = value
        with pytest.raises(ValueError, match='fresh directory'):
            QuantumChemistry(invalid)


def test_old_explicit_search_record_requires_new_directory(workspace):
    par = profile_parameters(l1_profile=gaussian())
    QuantumChemistry(par)
    path = Path('.kinbot_theory_jobs.json')
    record = json.loads(path.read_text())
    record['schema'] = 1
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='saved theory record lacks required settings'):
        QuantumChemistry(par)


@pytest.mark.parametrize('high_level', [0, 1])
@pytest.mark.parametrize('replacement_scan', [False, True])
def test_restart_refuses_rotor_points_from_the_other_level(
        workspace, high_level, replacement_scan):
    par = parameters(l1_profile=dict(gaussian('HF'), basis='STO-3G'),
                     l2_profile=gaussian('B3LYP'), high_level=high_level,
                     conformer_search=1, queuing='local')
    qc = QuantumChemistry(par)
    point = methanol()
    Path('hir').mkdir()
    fix = [[int(atom) + 1 for atom in point.dihed[0]]]
    geom = point.geom.copy()
    job = qc.qc_hir(point, geom, 0, 0, fix, False)
    if replacement_scan:
        geom[0, 0] += .1
        job = qc.qc_hir(point, geom, 0, 0, fix, False)
        assert '_recovery_' in job
    qc.db.write(molecule('CH3OH'), name=job, data=dict(energy=-100., status='normal'))
    Path(job + '.log').write_text('done\n')
    before = file_contents(workspace)
    restarted = QuantumChemistry(par)
    assert restarted.qc_hir(point, geom, 0, 0, fix, False) == job
    assert file_contents(workspace) == before

    changed = dict(par, high_level=1 - high_level)
    with pytest.raises(ValueError, match='rotor.*L[12].*L[12].*fresh directory'):
        QuantumChemistry(changed)
    assert file_contents(workspace) == before


def test_enabling_l2_before_rotor_scans_keeps_l1_results(workspace):
    par = profile_parameters(high_level=0, conformer_search=1)
    qc = QuantumChemistry(par)
    qc.qc_opt(methanol(), methanol().geom)
    before = file_contents(workspace)
    restarted = QuantumChemistry(dict(par, high_level=1))
    assert restarted._routes == qc._routes
    assert file_contents(workspace) == before


@pytest.mark.parametrize('saved_result', ['database', 'worker', 'rotor'])
def test_adding_profiles_refuses_unverified_old_calculations(workspace, saved_result):
    par = parameters(high_level=0)
    ordinary = QuantumChemistry(par)
    point = methanol()
    job = routing_name(point) + '_well'
    if saved_result == 'database':
        ordinary.db.write(molecule('CH3OH'), name=job, data=dict(
            energy=-100., zpe=.01, frequencies=[100.] * 12, status='normal'))
        old_file = Path('kinbot.db')
    elif saved_result == 'worker':
        old_file = Path(job + '.py')
        old_file.write_text('# previously submitted HF calculation\n')
    else:
        Path('hir').mkdir()
        old_file = Path('hir') / (job + '_hir_0_00.pkl')
        old_file.write_bytes(b'previous rotor result')
    before = old_file.read_bytes()
    changed = dict(par, l1_profile=gaussian('B3LYP'))

    with pytest.raises(ValueError, match='theory.*fresh directory'):
        QuantumChemistry(changed)

    assert old_file.read_bytes() == before
    assert not Path('.kinbot_theory_jobs.json').exists()


@pytest.mark.parametrize('route_file', ['.kinbot_theory_jobs.json', 'routes/jobs.json'])
@pytest.mark.parametrize('remove_path', [False, True])
def test_removing_profiles_refuses_existing_profile_results(workspace, route_file, remove_path):
    par = parameters(l1_profile=gaussian('B3LYP'), profiled_route_file=route_file)
    qc = QuantumChemistry(par)
    point = methanol()
    qc.qc_opt(point, point.geom)
    qc.db.write(molecule('CH3OH'), name=routing_name(point) + '_well', data=dict(
        energy=-100., zpe=.01, frequencies=[100.] * 12, status='normal'))
    before = file_contents(workspace)
    changed = dict(par, l1_profile={})
    if remove_path:
        changed.pop('profiled_route_file')

    with pytest.raises(ValueError, match='profiles.*removed'):
        QuantumChemistry(changed)

    assert file_contents(workspace) == before


def file_contents(directory):
    return {str(path.relative_to(directory)): path.read_bytes()
            for path in directory.rglob('*') if path.is_file()}


def test_existing_custom_record_gains_flag_only_on_calculation_restart(workspace):
    from kinbot.profiled_qc import validate_theory_restart
    from kinbot.run_format import ensure_current_run

    par = profile_parameters(profiled_route_file='routes/jobs.json')
    QuantumChemistry(par)
    marker = Path('.kinbot_run.json')
    previous = json.loads(marker.read_text())
    previous.pop('profiled_theory')
    marker.write_text(json.dumps(previous))
    before = file_contents(workspace)
    validate_theory_restart(par)
    assert file_contents(workspace) == before
    QuantumChemistry(par)
    assert ensure_current_run()['profiled_theory'] is True


@pytest.mark.parametrize('change', [
    'l1_model', 'l2_method', 'remove_profiles', 'orcasimpleinput', 'orcablocks',
    'integral', 'rotor_l1_to_l2', 'rotor_l2_to_l1',
])
@pytest.mark.parametrize('entry', ['postprocess', 'write_input', 'write_input_keep', 'cli'])
def test_pes_rejects_changed_theory_before_reading_or_overwriting_results(
        workspace, monkeypatch, change, entry):
    import logging
    import sys
    from kinbot import pes
    from kinbot.run_format import ensure_current_run

    options = {}
    if change in ('orcasimpleinput', 'orcablocks'):
        options = dict(l1_profile=dict(calculator='orca', method='B3LYP',
            basis='6-31G(d)'), orcasimpleinput='D3BJ', orcablocks='%method Grid 4 end')
    elif change == 'integral':
        options = dict(integral='UltraFine')
    elif change.startswith('rotor_'):
        options = dict(high_level=int(change == 'rotor_l2_to_l1'), conformer_search=1)
    par = profile_parameters(profiled_route_file='routes/jobs.json', queuing='local',
                             mult=1, **options)
    point = methanol()
    name = routing_name(point)
    ensure_current_run(create=True)
    Path(name).mkdir()
    with monkeypatch.context() as context:
        context.chdir(name)
        qc = QuantumChemistry(par)
        qc.db.write(molecule('CH3OH'), name=name + '_well', data=dict(
            energy=-100., zpe=.01, frequencies=[100.] * 12, status='normal'))
        if change.startswith('rotor_'):
            Path('hir').mkdir()
            job = qc.qc_hir(point, point.geom, 0, 0,
                           [[int(atom) + 1 for atom in point.dihed[0]]], False)
            qc.db.write(molecule('CH3OH'), name=job,
                        data=dict(energy=-100., status='normal'))
        Path(name + '.json').write_text(json.dumps(par))
    changed = copy.deepcopy(par)
    if change == 'l1_model':
        changed['l1_profile']['model_path'] = '/models/different.pt'
    elif change == 'l2_method':
        changed['l2_profile']['method'] = 'B3LYP'
    elif change == 'remove_profiles':
        for key in ('l0_profile', 'l1_profile', 'l2_profile', 'profiled_route_file'):
            changed.pop(key)
    elif change.startswith('rotor_'):
        changed['high_level'] = 1 - par['high_level']
    else:
        changed[change] = dict(orcasimpleinput='D4', orcablocks='%method Grid 5 end',
                              integral='FineGrid')[change]
    changed['smiles'] = 'CO'
    Path('input.json').write_text(json.dumps(changed))
    before = file_contents(workspace)
    monkeypatch.setattr(pes, 'get_energy', Mock(side_effect=AssertionError('read old energy')))
    monkeypatch.setattr(pes, 'config_log', lambda *a, **k: logging.getLogger('KinBot'))
    monkeypatch.setattr(QuantumChemistry, '__init__',
                        Mock(side_effect=AssertionError('created a calculator')))

    with pytest.raises(ValueError, match='Theory profiles|saved rotor scan'):
        if entry == 'postprocess':
            pes.postprocess(changed, [name], 'all', [], point.mass)
        elif entry == 'write_input':
            pes.write_input('input.json', point, 100., None, str(workspace), 0)
        elif entry == 'write_input_keep':
            pes.write_input_keep('input.json', name, str(workspace))
        else:
            monkeypatch.setattr(sys, 'argv', ['pes', 'input.json', 'no-kinbot'])
            pes._main()
    assert file_contents(workspace) == before


def test_pes_refuses_profiles_added_to_ordinary_results(workspace, monkeypatch):
    from kinbot import pes
    from kinbot.run_format import ensure_current_run

    par = profile_parameters()
    ensure_current_run('methanol', create=True)
    Path('methanol/kinbot.db').write_bytes(b'previous HF results')
    before = file_contents(workspace)
    monkeypatch.setattr(pes, 'get_energy', Mock(side_effect=AssertionError('read old energy')))
    with pytest.raises(ValueError, match='Cannot verify the theory'):
        pes.postprocess(par, ['methanol'], 'all', [], 32.)
    assert file_contents(workspace) == before


def test_profile_record_precedes_first_job_and_allows_startup_directories(workspace):
    from kinbot.run_format import ensure_current_run
    from kinbot.utils import make_dirs

    par = profile_parameters(conformer_search=1, rotor_scan=1)
    ensure_current_run(create=True)
    make_dirs(par)
    Path('kinbot.log').write_text('Starting KinBot\n')
    QuantumChemistry(par)
    recorded = Path('.kinbot_theory_jobs.json').read_bytes()
    assert json.loads(recorded)['jobs'] == {}
    # A restart before the first calculation must still use these profiles.
    QuantumChemistry(par)
    assert Path('.kinbot_theory_jobs.json').read_bytes() == recorded


def test_qchem_l0_keeps_ordinary_conformers_but_refuses_ring_jobs(workspace):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from kinbot.conformers import Conformers

    par = profile_parameters(l0_profile=dict(
        calculator='qchem', method='HF', basis='sto-3g', optimizer='native'),
        conformer_search=1)
    qc = QuantumChemistry(par)
    qc.l0.submit_qc = Mock(return_value=1)
    point = methanol()
    qc.qc_conf(point, point.geom, 0, semi_emp=1)
    qc.l0.submit_qc.assert_called_once()
    assert 'QChem' in script(f'conf/{routing_name(point)}_semi_emp_0000')
    qc.l0.submit_qc.reset_mock()

    mol = Chem.AddHs(Chem.MolFromSmiles('C1CCCCC1'))
    assert AllChem.EmbedMolecule(mol, randomSeed=81) == 0
    ring = StationaryPoint('cyclohexane', 0, 1,
        atom=[atom.GetSymbol() for atom in mol.GetAtoms()],
        geom=mol.GetConformer().GetPositions())
    ring.characterize()
    with pytest.raises(ValueError, match='Q-Chem.*ring.*Gaussian.*FairChem'):
        Conformers(ring, par, qc).generate_ring_conformers(ring.geom)
    qc.l0.submit_qc.assert_not_called()
    assert not list(Path('conf').glob('*_r*.py'))


@pytest.mark.parametrize('level', ['l0', 'l1', 'l2'])
@pytest.mark.parametrize('optimizer', ['sella', 'bfgs'])
def test_qchem_profiles_reject_force_workers_before_submission(workspace, level, optimizer):
    with pytest.raises(ValueError, match='Q-Chem requires optimizer=native'):
        profile_parameters(**{f'{level}_profile': dict(
            calculator='qchem', method='B3LYP', basis='6-31G(d)', optimizer=optimizer)})
    assert not Path('.kinbot_theory_jobs.json').exists()
    assert not list(Path('.').glob('*.py'))


@pytest.mark.parametrize('alias', [False, True])
def test_switching_to_qchem_defaults_to_native_workers(workspace, alias):
    profile = dict(method='B3LYP', basis='6-31G(d)', command='/l2/qchem')
    options = {'l2': 'qchem'} if alias else {}
    if not alias:
        profile['calculator'] = 'qchem'
    qc = QuantumChemistry(profile_parameters(l2_profile=profile, **options))
    assert qc.l2.optimizer == 'native' and not qc.l2.use_sella
    point = methanol()
    qc.qc_opt(point, point.geom, high_level=1)
    text = script(routing_name(point) + '_well_high')
    assert 'from sella import' not in text
    assert "kwargs['jobtype'] = 'freq'" in text
    assert '/l2/qchem' in text


def test_ordinary_qchem_optimizer_setting_is_unchanged(workspace):
    par = parameters(use_sella=1)
    par['qc'] = 'qchem'
    assert resolve_profiles(par) == {}
    assert QuantumChemistry(par).use_sella


@pytest.mark.parametrize('override, expected', [
    ({}, 'B3LYP/cc-pVTZ (EmpiricalDispersion=GD3BJ)'),
    ({'command': '/different/g16'}, 'B3LYP/cc-pVTZ (EmpiricalDispersion=GD3BJ)'),
    ({'method': 'M062X', 'basis': '6-31G(d)'},
     'M062X/6-31G(d) (EmpiricalDispersion=GD3BJ)'),
    ({'basis': 'cc-pVDZ'}, 'B3LYP/cc-pVDZ (EmpiricalDispersion=GD3BJ)'),
    ({'calculator_kwargs': {'EmpiricalDispersion': ''}}, 'B3LYP/cc-pVTZ'),
    ({'method': 'M062X', 'label': 'My selected theory'}, 'My selected theory'),
])
def test_explicit_profile_labels_describe_resolved_theory(workspace, override, expected):
    from kinbot.mess import MESS
    from kinbot.theory import calculation_label

    profile = dict(gaussian('B3LYP'), basis='cc-pVTZ')
    profile.update(override)
    par = parameters(l2_profile=profile,
                     calc_kwargs={'EmpiricalDispersion': 'GD3BJ'},
                     high_level=1, conformer_search=1)
    assert calculation_label(par) == expected
    assert MESS(par, methanol()).calculation_label() == expected


@pytest.mark.parametrize('high_level, expected', [
    (0, 'L1 (FairChem /models/medium.pt (omol))'),
    (1, 'L2 (FairChem /models/large.pt (omol))'),
])
def test_pes_logs_the_selected_profile(workspace, monkeypatch, caplog, high_level, expected):
    import logging
    from kinbot import pes
    par = profile_parameters(l2_profile=fairchem('/models/large.pt'),
                             high_level=high_level, conformer_search=1, queuing='local',
                             profiled_route_file='routes/jobs.json')
    Path('methanol').mkdir()
    with monkeypatch.context() as context:
        context.chdir('methanol')
        QuantumChemistry(par)
    Path('methanol/summary_methanol.out').write_text('# no reactions\n')
    before = file_contents(Path('methanol'))
    monkeypatch.setattr(QuantumChemistry, '__init__',
                        Mock(side_effect=AssertionError('created a calculator')))
    monkeypatch.setattr(pes, 'logger', logging.getLogger('KinBot'), raising=False)
    monkeypatch.setattr(pes, 'get_energy', Mock(return_value=(-1., .01)))
    monkeypatch.setattr(pes, 'get_l3energy', Mock(return_value=(0, -1.)))
    for name in ('copy_from_kinbot', 'create_pesviewer_input', 'create_rotdpy_inputs', 't1_analysis'):
        monkeypatch.setattr(pes, name, Mock())
    with caplog.at_level(logging.INFO):
        pes.postprocess(par, ['methanol'], 'all', [], 32.)
    assert f'Energies used are at the {expected} level of theory.' in caplog.text
    assert 'M062X' not in caplog.text
    assert file_contents(Path('methanol')) == before


def test_hessian_weighting_follows_source_job_backend(workspace):
    par = profile_parameters(l1_profile=dict(
        calculator='qchem', method='B3LYP', basis='6-31G(d)',
        command='/l1/qchem', optimizer='native', frequency_mode='native_hessian'))
    qc = QuantumChemistry(par)
    qc.l1.submit_qc('native_qchem_selected', 4)
    qc.l2.submit_qc('native_gaussian_selected', 4)
    qc.l0.submit_qc('uma_selected', 4)
    assert qc.hessian_is_massweighted('native_qchem_selected')
    assert not qc.hessian_is_massweighted('native_gaussian_selected')
    assert not qc.hessian_is_massweighted('uma_selected')
    expected = np.eye(3)
    qc.l1.read_qc_hess = Mock(return_value=expected)
    np.testing.assert_array_equal(qc.read_qc_hess('native_qchem_selected', 1), expected)


def test_frequency_recovery_uses_source_profile_despite_caller_flag(workspace):
    from kinbot.ase_modules.calculators.gaussian import Gaussian

    qc = QuantumChemistry(profile_parameters())
    point = methanol()
    source = 'conf/selected_with_nonstandard_name'
    qc.l2.submit_qc(source, 4)
    qc.l2.db.write(molecule('CH3OH'), name=source, data={'status': 'normal'})
    job = qc.qc_freq(point, source, high_level=0)
    kwargs = literal(script(job), 'kwargs')
    assert kwargs['method'] == 'M062X'
    assert kwargs['basis'] == '6-31G(d)'
    assert 'opt' not in kwargs
    calc = Gaussian(**kwargs)
    calc.write_input(molecule('CH3OH'), properties=['energy'])
    assert (Path(calc.directory) / kwargs['chk']).resolve() == Path(job).resolve()
    assert f"%chk={kwargs['chk']}" in Path(job + '.com').read_text()
    assert qc.for_job(job) is qc.l2


def test_published_result_keeps_source_profile(workspace):
    qc = QuantumChemistry(profile_parameters())
    qc.l0.submit_qc('initial_candidate', 4)
    qc.l0.publish_result = Mock(return_value=17)
    source = SimpleNamespace(name='initial_candidate')
    assert qc.publish_result(source, 'accepted_candidate') == 17
    assert qc.for_job('accepted_candidate') is qc.l0


def test_l3_arguments_keep_ordinary_global_settings(workspace):
    qc = QuantumChemistry(profile_parameters())
    kwargs = qc.get_qc_arguments('single_point', 1, 0, 18, L3=True)
    ordinary = qc.legacy.get_qc_arguments('single_point', 1, 0, 18, L3=True)
    assert kwargs == ordinary
    assert kwargs['method'] == 'HF'
    assert kwargs['basis'] == 'sto-3g'
    assert qc.legacy.qc_command == '/global/g09'




def test_sella_gaussian_workers_receive_profile_frequency_policy(workspace):
    profile = gaussian()
    profile['optimizer'] = 'sella'
    qc = QuantumChemistry(profile_parameters(l2_profile=profile))
    point = methanol()
    qc.qc_opt(point, point.geom, high_level=1)
    text = script(routing_name(point) + '_well_high')
    assert 'profile_vibrations' in text
    assert 'native_hessian' in text
    ts = methanol(wellorts=1, name='selected_ts')
    qc.qc_opt_ts(ts, ts.geom, high_level=1)
    text = script('selected_ts_high')
    assert 'profile_vibrations' in text
    assert 'native_hessian' in text


def test_imaginary_mode_reader_follows_profile_frequency_policy(workspace, monkeypatch):
    from kinbot import reader_gauss, reader_sella

    profile = gaussian()
    profile['optimizer'] = 'sella'
    qc = QuantumChemistry(profile_parameters(l2_profile=profile))
    native = Mock(return_value=np.array([1., 0., 0.]))
    finite_difference = Mock(return_value=np.array([0., 1., 0.]))
    monkeypatch.setattr(reader_gauss, 'read_imag_mode', native)
    monkeypatch.setattr(reader_sella, 'read_imag_mode', finite_difference)
    qc.l2.submit_qc('native_selected', 4)
    qc.l1.submit_qc('uma_selected', 4)
    np.testing.assert_array_equal(qc.read_imag_mode('native_selected', 1), [1., 0., 0.])
    np.testing.assert_array_equal(qc.read_imag_mode('uma_selected', 1), [0., 1., 0.])
    native.assert_called_once_with('native_selected', 1)
    finite_difference.assert_called_once_with('uma_selected', 1)


@pytest.mark.parametrize('case', [1, 2])
def test_requested_three_level_methanol_settings_are_accepted(workspace, case):
    if case == 1:
        profiles = dict(
            l0_profile=dict(calculator='gaussian', method='AM1', basis='',
                            command='g16', optimizer='sella', frequency_mode='ase_forces'),
            l1_profile=fairchem('uma-s-1p1', 'cuda'),
            l2_profile=fairchem('uma-m-1p1', 'cuda'))
    else:
        profiles = dict(l0_profile=fairchem('uma-s-1p1', 'cuda'),
                        l1_profile=fairchem('uma-m-1p1', 'cuda'),
                        l2_profile=gaussian())
    point = methanol()
    structure = [item for atom, xyz in zip(point.atom, point.geom)
                 for item in (str(atom), *map(float, xyz))]
    par = parameters(structure=structure, charge=0, mult=1, high_level=1,
                     conformer_search=1, l0_conformer_search=1, pes=0, L3_calc=0,
                     me=0, run_me=0, **profiles)
    assert set(resolve_profiles(par)) == {'l0', 'l1', 'l2'}


def test_qchem_profile_controls_early_ts_method_and_command(workspace):
    par = profile_parameters(l1_profile=dict(calculator='qchem', method='wb97x-d',
                            basis='def2-svp', command='/l1/qchem', optimizer='native'))
    qc = QuantumChemistry(par)
    rxn = reaction(qc, par)
    qc.l1.is_in_database = Mock(return_value=False)
    carry_out_reaction(rxn, 0, '/global/qchem')
    kwargs = literal(script(rxn.instance_name), 'kwargs')
    assert kwargs['method'] == 'b3lyp'
    assert kwargs['basis'] == 'sto-3g'
    assert kwargs['command'] == '/l1/qchem'
    legacy = QuantumChemistry({**par, 'qc': 'qchem'},
                              profile_backend=True)
    kwargs = legacy.get_qc_arguments('ts', 1, 0, 18, ts=1, step=0, max_step=2)
    assert kwargs['method'] == 'b3lyp' and kwargs['basis'] == 'sto-3g'


def test_calculator_instances_use_the_profile_command(workspace):
    from kinbot.ase_modules.calculators.gaussian import Gaussian
    from kinbot.ase_modules.calculators.qchem import QChem
    qc = QuantumChemistry(profile_parameters())
    arguments = qc.l2.merge_kwargs(qc.l2.get_qc_arguments('methanol', 1, 0, 18))
    calculator = Gaussian(**arguments)
    assert calculator.command == '/l2/g16 < PREFIX.com > PREFIX.log'
    calculator = QChem(command='/site/qchem', nt=4)
    assert calculator.command == '/site/qchem -nt 4 PREFIX.inp PREFIX.out'


def test_native_frequency_after_sella_keeps_geometry_and_executable(workspace, monkeypatch):
    from kinbot.ase_modules.calculators.gaussian import Gaussian
    from kinbot.frequencies import profile_vibrations
    from kinbot import reader_gauss, utils
    mol = molecule('CH3OH')
    before = mol.positions.copy()
    mol.calc = Gaussian(method='M062X', basis='6-31G(d)', force='',
                        opt='TS', guess='Read', geom='Check',
                        command='/profile/g16 < PREFIX.com > PREFIX.log')
    monkeypatch.setattr(Gaussian, 'get_potential_energy', lambda *a, **k: -1.)
    monkeypatch.setattr(utils, 'iowait', lambda *a: None)
    monkeypatch.setattr(reader_gauss, 'read_freq', lambda *a: [100., 200.])
    monkeypatch.setattr(reader_gauss, 'read_zpe', lambda *a: .01)
    freqs, zpe, hess = profile_vibrations(mol, 'conf/methanol', 'native_hessian')
    assert (freqs, zpe, hess) == ([100., 200.], .01, None)
    assert np.array_equal(mol.positions, before)
    assert mol.calc.command.startswith('/profile/g16 ')
    assert mol.calc.parameters['chk'] == 'methanol'
    assert not {'force', 'opt', 'guess', 'geom'} & set(mol.calc.parameters)


def test_saved_hir_without_a_route_is_not_assigned_a_level(workspace):
    qc = QuantumChemistry(profile_parameters(high_level=0))
    job = 'hir/methanol_hir_0_01'
    Path('hir').mkdir()
    Path(job + '_sella.log').write_text('done\n')
    qc.db.write(molecule('CH3OH'), name=job, data={'energy': -1., 'status': 'normal'})
    for high_level in (0, 1):
        qc.par['high_level'] = high_level
        with pytest.raises(ValueError, match='no recorded theory level'):
            qc.get_qc_energy(job)
        with pytest.raises(ValueError, match='no recorded theory level'):
            qc.for_level('l2' if high_level else 'l1').check_qc(job)
    qc._record(job, 'l1')
    assert qc.for_job(job) is qc.l1
    assert qc.get_qc_energy(job)[0] == 0


def test_profiled_sella_conformer_handoff_preserves_logs_modes_and_checkpoint(workspace, monkeypatch):
    import runpy
    import sella
    from ase.calculators.calculator import Calculator, all_changes
    from ase.io import write
    from kinbot import frequencies
    from kinbot.ase_modules.calculators import gaussian as gaussian_module

    par = profile_parameters(l1_profile=dict(calculator='gaussian', method='B3LYP',
        basis='6-31G(d)', command='/never/executed', optimizer='sella',
        frequency_mode='ase_forces'), queuing='local')
    qc = QuantumChemistry(par)
    point = methanol(wellorts=1, name='methanol_ts')
    qc.qc_conf(point, point.geom, 0)
    job = 'conf/methanol_ts_0000'
    kwargs = literal(script(job), 'kwargs')
    calc = gaussian_module.Gaussian(**kwargs)
    assert (Path(calc.directory) / kwargs['chk']).resolve() == Path(job).resolve()
    initial_guesses = []

    class FakeGaussian(Calculator):
        implemented_properties = ['energy', 'forces']

        def calculate(self, atoms=None, properties=None, system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            initial_guesses.append(self.parameters.get('guess'))
            self.results = {'energy': -10., 'forces': np.zeros((len(atoms), 3))}

    class FakeSella:
        def __init__(self, mol, trajectory, logfile, **kwargs):
            assert initial_guesses == [None]
            assert mol.calc.parameters['guess'] == 'Read'
            self.mol, self.trajectory, self.logfile = mol, trajectory, logfile

        def run(self, **kwargs):
            write(self.trajectory, self.mol)
            Path(self.logfile).write_text('mock converged\n')
            return True

    def vibrate(mol, label, mode):
        path = Path(Path(label).name + '_vib')
        path.mkdir()
        path.joinpath('vib.xyz').write_text(str(len(mol)) + '\nmode\n' + ''.join(
            f'{symbol} 0 0 0 1 2 3\n' for symbol in mol.get_chemical_symbols()))
        return [-100.] + [100.] * (3 * len(mol) - 7), .01, np.eye(3 * len(mol))

    monkeypatch.setattr(gaussian_module, 'Gaussian', FakeGaussian)
    monkeypatch.setattr(sella, 'Sella', FakeSella)
    monkeypatch.setattr(frequencies, 'profile_vibrations', vibrate)
    runpy.run_path(job + '.py')
    assert Path(job + '_sella.log').read_text().endswith('done\n')
    assert Path(job + '_vib.xyz').is_file()
    np.testing.assert_array_equal(qc.read_imag_mode(job, point.natom), [[1, 2, 3]] * point.natom)
    assert qc.check_qc(job) == 'normal'


@pytest.mark.parametrize('profiled_sella, semi_emp, multi_conf_tst', [
    (True, 1, 0), (True, 0, 0), (False, 0, 1), (False, 1, 0)])
def test_gaussian_conformer_checkpoint_matches_execution_directory(
        workspace, profiled_sella, semi_emp, multi_conf_tst):
    from kinbot.ase_modules.calculators.gaussian import Gaussian

    if profiled_sella:
        profile = dict(calculator='gaussian', method='AM1', basis='',
                       command='/never/executed', optimizer='sella',
                       frequency_mode='ase_forces')
        par = profile_parameters(l0_profile=profile, l1_profile=profile)
    else:
        par = parameters(use_sella=0, multi_conf_tst=multi_conf_tst,
                         conformer_search=1)
    qc = QuantumChemistry(par)
    point = methanol()
    qc.qc_conf(point, point.geom, 0, semi_emp=semi_emp)
    suffix = 'semi_emp_0000' if semi_emp else '0000'
    job = f'conf/{routing_name(point)}_{suffix}'
    kwargs = literal(script(job), 'kwargs')
    calc = Gaussian(**kwargs)
    calc.write_input(molecule('CH3OH'), properties=['energy'])
    assert Path(calc.directory).resolve() == workspace / 'conf'
    text = Path(job + '.com').read_text()
    if profiled_sella or multi_conf_tst:
        checkpoint = (Path(calc.directory) / kwargs['chk']).resolve()
        assert checkpoint == workspace / job
        assert checkpoint.parent.is_dir()
        assert f"%chk={kwargs['chk']}" in text
    else:
        assert 'chk' not in kwargs
        assert '%chk' not in text


@pytest.mark.parametrize('formula', ['H2CO', 'CO', 'H2', 'H'])
def test_gaussian_startup_failure_finishes_profiled_conformer_as_error(
        workspace, monkeypatch, capsys, formula):
    import runpy
    import sella
    from ase import optimize
    from kinbot import frequencies
    from kinbot.ase_modules.calculators.gaussian import Gaussian

    profile = dict(calculator='gaussian', method='AM1', basis='',
                   command='/never/executed', optimizer='sella',
                   frequency_mode='ase_forces')
    qc = QuantumChemistry(profile_parameters(l0_profile=profile))
    atoms = molecule(formula)
    point = StationaryPoint(formula, 0, 2 if formula == 'H' else 1,
                            atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    point.characterize()
    qc.qc_conf(point, point.geom, 0, semi_emp=1)
    job = f'conf/{routing_name(point)}_semi_emp_0000'
    calculate = Mock(side_effect=RuntimeError('Gaussian startup failed'))
    optimizer = Mock(side_effect=AssertionError('must not optimize after startup failure'))
    vibrate = Mock(side_effect=AssertionError('must not evaluate frequencies after startup failure'))
    monkeypatch.setattr(Gaussian, 'calculate', calculate)
    monkeypatch.setattr(sella, 'Sella', optimizer)
    monkeypatch.setattr(optimize, 'BFGS', optimizer)
    monkeypatch.setattr(frequencies, 'profile_vibrations', vibrate)

    runpy.run_path(job + '.py')

    calculate.assert_called_once()
    optimizer.assert_not_called()
    vibrate.assert_not_called()
    rows = list(qc.l0.db.select(name=job))
    assert len(rows) == 1 and rows[0].data.status == 'error'
    assert 'energy' not in rows[0].data
    assert Path(job + '_sella.log').read_text().endswith('done\n')
    assert 'Gaussian startup failed' in capsys.readouterr().err
    # Simulate a finished Slurm job; missing files would start the grace wait.
    import kinbot.qc as qc_module
    process = Mock()
    process.communicate.return_value = (b'JOBID NAME\n', b'')
    monkeypatch.setattr(qc_module.subprocess, 'Popen', Mock(return_value=process))
    qc.l0.queuing = 'slurm'
    assert qc.check_qc(job) == 'error'


@pytest.mark.parametrize('calculator', ['orca', 'fairchem'])
def test_force_calculators_reject_unsupported_optimizer(workspace, calculator):
    with pytest.raises(ValueError, match='requires optimizer=sella or bfgs'):
        profile_parameters(l1_profile=dict(calculator=calculator, optimizer='native'))
    assert not Path('kinbot.db').exists()


def test_error_names_inherited_level_instead_of_unsupplied_profile(workspace):
    par = parameters(l0_profile=fairchem('/models/small.pt'))
    par.update(qc='qchem', use_sella=1)
    with pytest.raises(ValueError, match='L1 settings inherited from qc/use_sella'):
        resolve_profiles(par)


@pytest.mark.parametrize('status', ['normal', 'error', 'running', 0])
def test_reuse_only_never_submits_l1(workspace, monkeypatch, status):
    submit = Mock(return_value=1)
    monkeypatch.setattr(QuantumChemistry, 'submit_qc', submit)
    par = profile_parameters(profiled_l1_reuse_only=1)
    qc = QuantumChemistry(par)
    qc.l1.check_qc = Mock(return_value=status)
    if status == 0:
        with pytest.raises(MissingL1Result, match='no replacement job was submitted'):
            qc.l1.submit_qc('missing_well', 4)
    else:
        assert qc.l1.submit_qc('existing_well', 4, singlejob=0) == 0
    submit.assert_not_called()
    assert qc.l2.submit_qc('well_high', 4) == 1
    submit.assert_called_once()


@pytest.mark.parametrize('value', [2, '1'])
def test_reuse_only_input_is_validated(workspace, value):
    with pytest.raises(ValueError, match='profiled_l1_reuse_only must be 0 or 1'):
        profile_parameters(profiled_l1_reuse_only=value)


def test_reuse_only_cannot_be_ignored_without_profiles(workspace):
    with pytest.raises(ValueError, match='requires theory profiles'):
        parameters(profiled_l1_reuse_only=1)


def test_job_level_is_recorded_before_submission_returns(workspace, monkeypatch):
    def launch(backend, job, *args, **kwargs):
        routes = json.loads(Path('.kinbot_theory_jobs.json').read_text())['jobs']
        assert routes[job]['level'] == backend.theory_level
        assert routes[job]['job_id'] is None
        assert routes[job + '_prod']['level'] == backend.theory_level
        backend.job_ids[job] = '42'
        return 1
    monkeypatch.setattr(QuantumChemistry, 'submit_qc', launch)
    qc = QuantumChemistry(profile_parameters())
    qc.l2.submit_qc('test_IRC_F', 4)
    assert qc._routes['test_IRC_F']['job_id'] == '42'
    save = Mock(wraps=qc._save_routes)
    qc._save_routes = save
    qc._record('test_IRC_F', 'l2', '42')
    save.assert_not_called()
    with pytest.raises(ValueError, match='cannot replace recorded L2 with L1'):
        qc.l1.submit_qc('test_IRC_F', 4)


def test_missing_router_attribute_raises_without_recursion():
    qc = object.__new__(ProfiledQuantumChemistry)
    with pytest.raises(AttributeError, match='l1'):
        getattr(qc, 'l1')
    with pytest.raises(AttributeError, match='__setstate__'):
        getattr(qc, '__setstate__')


def test_profiled_product_registration_preserves_chemical_context(workspace):
    from ase import Atoms
    from kinbot.calculation import load_calculation_record
    from kinbot.species_routing import register_optimized_product
    from kinbot.stereo_routing import guard_well_job
    par = profile_parameters(queuing='local')
    qc = QuantumChemistry(par)
    initial = StationaryPoint('initial', 0, 1, smiles='C[C@H](O)CC')
    initial.characterize()
    qc.qc_opt(initial, initial.geom, initial_product=True)
    job = routing_name(initial) + '_product'
    product = StationaryPoint('optimized', 0, 1, atom=initial.atom,
                              geom=initial.geom * [-1., 1., 1.])
    product.characterize()
    data = dict(status='normal', energy=-100., zpe=.01,
                frequencies=[100.] * (3 * product.natom-6),
                hess=np.eye(3 * product.natom), charge=0, multiplicity=1)
    row = qc.db.write(Atoms(product.atom, positions=product.geom), name=job, data=data)
    Path(job + '_sella.log').write_text('done\n')
    accepted = register_optimized_product(qc, initial, product, job)
    target = routing_name(accepted) + '_well'
    assert qc.for_job(target) is qc.l1
    saved = next(qc.db.select(name=target, sort='-id', limit=1))
    assert saved.data['copied_from_row_id'] == row
    from kinbot.stereo_routing import _chemical_context
    assert saved.data['chemical_context'] == _chemical_context(product)
    restarted = QuantumChemistry(par)
    guard_well_job(restarted, accepted, accepted.geom, target)
    load_calculation_record(accepted, restarted, target)
    np.testing.assert_array_equal(accepted.hess, data['hess'])


@pytest.mark.parametrize('job', ['well', 'conf/well_0000', 'hir/well_hir_0_01'])
def test_gaussian_checkpoint_remains_relative_after_directory_move(workspace, job):
    from kinbot.ase_modules.calculators.gaussian import Gaussian
    qc = QuantumChemistry(parameters())
    kwargs = qc.get_qc_arguments(job, 1, 0, 18)
    assert kwargs['chk'] == Path(job).name
    # ASE runs the executable in the label directory; the same checkpoint
    # filename must work under any relocated calculation root.
    calculator = Gaussian(**kwargs)
    relocated = Path('/another/run')
    assert relocated / calculator.directory / kwargs['chk'] == relocated / job


@pytest.mark.parametrize('calculator', ['fc', 'orca', 'qchem'])
def test_flat_l0_backend_switch_has_same_optimizer_as_alias(workspace, calculator):
    flat = resolve_profiles(parameters(l0_qc=calculator, calc_kwargs={'gaussian_only': 1}))
    alias = resolve_profiles(parameters(l0=calculator, calc_kwargs={'gaussian_only': 1}))
    assert flat['l0'].optimizer == alias['l0'].optimizer
    assert flat['l0'].calculator == alias['l0'].calculator
    assert flat['l0'].calculator_kwargs == alias['l0'].calculator_kwargs == {}


@pytest.mark.parametrize('status', ['normal', 'error', 'running', 'incomplete', 'missing'])
@pytest.mark.parametrize('step', [0, 1])
def test_reuse_only_reaction_search_finishes_or_reports_missing_ts(
        workspace, monkeypatch, status, step):
    submit = Mock(side_effect=AssertionError('reuse-only must not submit'))
    monkeypatch.setattr(QuantumChemistry, 'submit_qc', submit)
    par = profile_parameters(profiled_l1_reuse_only=1, queuing='local')
    qc = QuantumChemistry(par)
    rxn = reaction(qc, par)
    job = rxn.instance_name
    qc._record(job, 'l1')
    if status != 'missing':
        frequencies = [-100., 100.] if status == 'normal' else []
        qc.db.write(molecule('CH3OH'), name=job, data=dict(
            status='error' if status == 'error' else 'normal',
            energy=-10., frequencies=frequencies))
        Path(job + '_sella.log').write_text('done\n')
    if status == 'running':
        qc.l1.check_qc = Mock(return_value='running')
    if status in ('missing', 'incomplete'):
        with pytest.raises(MissingL1Result, match='completed L1 TS|L1 TS search is incomplete'):
            carry_out_reaction(rxn, step, '')
    else:
        result = carry_out_reaction(rxn, step, '')
        assert result == {'normal': rxn.max_step+1, 'error': -1, 'running': step}[status]
    assert not Path(job + '.py').exists()
    submit.assert_not_called()


@pytest.mark.parametrize('missing', ['ts', 'frequencies', 'irc'])
def test_reuse_only_skips_one_reaction_and_keeps_its_files(
        workspace, monkeypatch, caplog, missing):
    from kinbot.reaction_generator import ReactionGenerator
    submit = Mock(side_effect=AssertionError('reuse-only must not submit'))
    monkeypatch.setattr(QuantumChemistry, 'submit_qc', submit)
    par = profile_parameters(profiled_l1_reuse_only=1, queuing='local',
                             delete_intermediate_files=1)
    qc = QuantumChemistry(par)
    parent = methanol()
    bad, good = [reaction(qc, par, name) for name in ('missing_ts', 'completed_ts')]
    for item in (bad, good):
        item.species = parent
        item.products = []
        qc._record(item.instance_name, 'l1')
    parent.reac_obj = [bad, good]
    parent.reac_inst = [bad.instance, good.instance]
    parent.reac_name = [bad.instance_name, good.instance_name]
    parent.reac_type = ['intra_H_migration'] * 2
    parent.reac_step = [0, 0]
    parent.reac_ts_done = [int(missing == 'irc'), 0]
    for job in (routing_name(parent) + '_well', bad.instance_name, good.instance_name):
        if job == bad.instance_name and missing == 'ts':
            continue
        frequencies = [] if job == bad.instance_name and missing == 'frequencies' else [-100., 100.]
        qc.db.write(molecule('CH3OH'), name=job, data=dict(
            status='normal', energy=-10., zpe=.02, frequencies=frequencies))
        Path(job + '_sella.log').write_text('done\n')
    checkpoint = Path(bad.instance_name + '.chk')
    checkpoint.write_bytes(b'keep the earlier calculation')
    original_rows = qc.db.count()

    # Stop after the unaffected TS has passed its frequency check, before IRC.
    class CheckedBothReactions(Exception):
        pass
    def after_poll(_):
        if parent.reac_ts_done[1] == 1:
            raise CheckedBothReactions
    monkeypatch.setattr('kinbot.reaction_generator.time.sleep', after_poll)
    with pytest.raises(CheckedBothReactions):
        ReactionGenerator(parent, par, qc, 'input.json').generate()
    assert parent.reac_ts_done == [-999, 1]
    assert parent.reac_step[1] == good.max_step + 1
    assert bad.missing_l1_result in caplog.text
    assert 'other reactions can continue' in caplog.text
    assert checkpoint.read_bytes() == b'keep the earlier calculation'
    assert qc.db.count() == original_rows
    submit.assert_not_called()


@pytest.mark.parametrize('high_level', [0, 1])
def test_hir_description_names_the_calculator_used_for_scans(workspace, high_level):
    from kinbot.hindered_rotors import HIR, recover_hir_model
    par = profile_parameters(high_level=high_level, conformer_search=1, rotor_scan=1)
    qc = QuantumChemistry(par)
    point = methanol()
    point.hir = HIR(point, qc, par)
    Path('hir').mkdir()
    point.hir.generate_hir_geoms(point.geom, False)
    expected = 'gauss' if high_level else 'fc'
    assert point.hir.scan_reference['backend'] == expected
    for job in point.hir.scan_jobs[0]:
        assert qc.for_job(job).qc == expected
    point.hir.scan_reference['source_job'] = 'outdated_reference'
    recover_hir_model(point, qc, par)
    assert point.hir.scan_reference['backend'] == expected


def test_failed_ts_frequency_calculation_writes_error_and_completion(workspace, monkeypatch):
    import runpy
    import sella
    from ase.io import write
    from kinbot import frequencies
    from kinbot.ase_modules.calculators.gaussian import Gaussian
    par = profile_parameters(l1_profile=dict(gaussian(), optimizer='sella'))
    qc = QuantumChemistry(par)
    rxn = reaction(qc, par)
    qc.l1.check_qc = Mock(return_value='normal')
    qc.l1.get_qc_geom = Mock(return_value=(0, rxn.species.geom))
    carry_out_reaction(rxn, rxn.max_step, '')
    class Converged:
        def __init__(self, mol, trajectory, **kwargs):
            self.mol, self.trajectory = mol, trajectory
        def run(self, **kwargs):
            write(self.trajectory, self.mol)
            return True
    monkeypatch.setattr(sella, 'Sella', Converged)
    monkeypatch.setattr(Gaussian, 'get_potential_energy', lambda *a, **k: -1.)
    monkeypatch.setattr(frequencies, 'profile_vibrations',
                        Mock(side_effect=RuntimeError('frequency calculation failed')))
    runpy.run_path(rxn.instance_name + '.py')
    assert qc.db.get(name=rxn.instance_name).data['status'] == 'error'
    assert Path(rxn.instance_name + '_sella.log').read_text().endswith('done\n')


def test_record_with_retired_barrierless_settings_still_restarts(workspace):
    import hashlib
    # A record written before the barrierless_saddle family was removed
    # carries its settings; they do not identify the method.
    par = profile_parameters(l1_profile=gaussian())
    qc = QuantumChemistry(par)
    qc.l1.submit_qc('methanol_well', 4)
    path = Path('.kinbot_theory_jobs.json')
    record = json.loads(path.read_text())
    for level in record['profiles'].values():
        level.setdefault('calculation_settings', {}).update({
            'barrierless_saddle_method': 'b3lyp', 'barrierless_saddle_basis': '6-31G',
            'barrierless_saddle_method_high': 'b3lyp', 'barrierless_saddle_basis_high': '6-31G'})
    record['profile_sha256'] = hashlib.sha256(json.dumps(
        record['profiles'], sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    path.write_text(json.dumps(record))
    assert QuantumChemistry(par)._routes == qc._routes
    # A real method change in such a record is still refused.
    record['profiles']['l1']['calculation_settings']['scan_method'] = 'hf'
    record['profile_sha256'] = hashlib.sha256(json.dumps(
        record['profiles'], sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='fresh directory'):
        QuantumChemistry(par)
