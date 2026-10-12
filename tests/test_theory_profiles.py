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
from kinbot.profiled_qc import ProfiledQuantumChemistry
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


@pytest.mark.parametrize('option', ['ppn', 'queue_job_limit'])
def test_per_level_execution_options_are_not_silently_accepted(workspace, option):
    profile = fairchem('/models/small.pt')
    profile[option] = 8
    with pytest.raises(ValueError, match='Unknown l0_profile options'):
        profile_parameters(l0_profile=profile)


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
        assert kwargs['method'] == 'B3LYP'
        assert kwargs['basis'] == '6-31G(d)'
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


@pytest.mark.parametrize('change', ['l1_model', 'l2_method', 'remove_profiles'])
@pytest.mark.parametrize('entry', ['postprocess', 'write_input', 'write_input_keep', 'cli'])
def test_pes_rejects_changed_theory_before_reading_or_overwriting_results(
        workspace, monkeypatch, change, entry):
    import logging
    import sys
    from kinbot import pes
    from kinbot.run_format import ensure_current_run

    par = profile_parameters(profiled_route_file='routes/jobs.json', queuing='local', mult=1)
    point = methanol()
    name = routing_name(point)
    ensure_current_run(create=True)
    Path(name).mkdir()
    with monkeypatch.context() as context:
        context.chdir(name)
        qc = QuantumChemistry(par)
        qc.db.write(molecule('CH3OH'), name=name + '_well', data=dict(
            energy=-100., zpe=.01, frequencies=[100.] * 12, status='normal'))
        Path(name + '.json').write_text(json.dumps(par))
    changed = copy.deepcopy(par)
    if change == 'l1_model':
        changed['l1_profile']['model_path'] = '/models/different.pt'
    elif change == 'l2_method':
        changed['l2_profile']['method'] = 'B3LYP'
    else:
        for key in ('l0_profile', 'l1_profile', 'l2_profile', 'profiled_route_file'):
            changed.pop(key)
    changed['smiles'] = 'CO'
    Path('input.json').write_text(json.dumps(changed))
    before = file_contents(workspace)
    monkeypatch.setattr(pes, 'get_energy', Mock(side_effect=AssertionError('read old energy')))
    monkeypatch.setattr(pes, 'config_log', lambda *a, **k: logging.getLogger('KinBot'))
    monkeypatch.setattr(QuantumChemistry, '__init__',
                        Mock(side_effect=AssertionError('created a calculator')))

    with pytest.raises(ValueError, match='Theory profiles'):
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


@pytest.mark.parametrize('level', ['l0', 'l1', 'l2'])
@pytest.mark.parametrize('optimizer', ['sella', 'native'])
def test_nwchem_profiles_fail_before_generating_jobs(workspace, level, optimizer):
    with pytest.raises(ValueError, match='NWChem.*Gaussian.*FairChem'):
        profile_parameters(**{f'{level}_profile': dict(
            calculator='nwchem', method='HF', basis='sto-3g', optimizer=optimizer)})
    assert not Path('.kinbot_theory_jobs.json').exists()
    assert not list(Path('.').glob('*.py'))


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
    qc = QuantumChemistry(profile_parameters())
    point = methanol()
    source = 'selected_with_nonstandard_name'
    qc.l2.submit_qc(source, 4)
    qc.l2.db.write(molecule('CH3OH'), name=source, data={'status': 'normal'})
    job = qc.qc_freq(point, source, high_level=0)
    kwargs = literal(script(job), 'kwargs')
    assert kwargs['method'] == 'M062X'
    assert kwargs['basis'] == '6-31G(d)'
    assert 'opt' not in kwargs
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
    assert kwargs['method'] == 'wb97x-d'
    assert kwargs['basis'] == 'def2-svp'
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


def test_saved_hir_without_a_route_follows_enabled_refinement(workspace):
    qc = QuantumChemistry(profile_parameters(high_level=0))
    assert qc.for_job('hir/methanol_hir_0_01') is qc.l1
    qc.par['high_level'] = 1
    assert qc.for_job('hir/methanol_hir_0_01') is qc.l2
    qc._record('hir/methanol_hir_0_01', 'l1')
    assert qc.for_job('hir/methanol_hir_0_01') is qc.l1


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
    assert Path(literal(script(job), 'kwargs')['chk']).resolve() == Path(job).resolve()

    class FakeGaussian(Calculator):
        implemented_properties = ['energy', 'forces']

        def calculate(self, atoms=None, properties=None, system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            self.results = {'energy': -10., 'forces': np.zeros((len(atoms), 3))}

    class FakeSella:
        def __init__(self, mol, trajectory, logfile, **kwargs):
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
