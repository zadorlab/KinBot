"""Network selection tests; simple blocks test grouping, not molecular physics."""
import json
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot.mess_mirrors import complete_mirror_channels
from kinbot.mess_networks import connected_models, write_network_inputs


def population(name, own, mirror, kind='Well'):
    return (f' {kind} {name}\n! kinbot_population {own} {mirror}\n'
            ' ZeroEnergy[kcal/mol] 0\n End\n')


def route(name, left, right, keys, path='ordinary'):
    return (f'! kinbot_stereopath {path}\n Barrier {name} {left} {right}\n'
            f'! kinbot_mirror_endpoints {" ".join(keys)}\n'
            ' RRHO\n ZeroEnergy[kcal/mol] 30\n End\n End\n')


def network(*blocks, root='r'):
    return f'Reactant {root}\nMicroRateOutput micro.out\nModel\n' + ''.join(blocks) + 'End ! end kinetics\n'


def independent_wells():
    return network(population('r', 'R', 'S'), population('s', 'S', 'R'),
                   population('b', 'B', 'B', 'Bimolecular'),
                   route('tr', 'r', 'b', ('R', 'S', 'B', 'B')),
                   route('ts', 's', 'b', ('S', 'R', 'B', 'B')))


def test_specified_reactant_does_not_acquire_mirror_from_product_sink():
    text = network(population('r', 'R', 'S'), population('b', 'B', 'B', 'Bimolecular'),
                   route('tr', 'r', 'b', ('R', 'S', 'B', 'B')))
    completed = complete_mirror_channels(text)
    assert 'derived by global reflection' not in completed
    models = connected_models(completed)
    assert len(models) == 1
    assert models[0]['wells'] == ['r']
    assert models[0]['products'] == ['b']


def test_bound_achiral_intermediate_does_allow_other_entrance_mirror():
    text = network(population('r', 'R', 'S'), population('a', 'A', 'A'),
                   route('tr', 'r', 'a', ('R', 'S', 'A', 'A')))
    completed = complete_mirror_channels(text)
    models = connected_models(completed)
    assert len(models) == 1
    assert len(models[0]['wells']) == 3
    assert len(models[0]['barriers']) == 2
    assert '! kinbot_population S R' in completed
    assert complete_mirror_channels(completed) == completed


def test_achiral_well_still_produces_both_configured_product_mirrors():
    text = network(population('a', 'A', 'A'), population('b', 'R', 'S', 'Bimolecular'),
                   route('tr', 'a', 'b', ('A', 'A', 'R', 'S')), root='a')
    models = connected_models(complete_mirror_channels(text))
    assert len(models) == 1
    assert models[0]['wells'] == ['a']
    assert len(models[0]['products']) == len(models[0]['barriers']) == 2


def test_separate_requested_wells_sharing_products_become_separate_calculations():
    text = complete_mirror_channels(independent_wells())
    models = connected_models(text)
    assert [model['wells'] for model in models] == [['r'], ['s']]
    assert [model['products'] for model in models] == [['b'], ['b']]
    assert [model['barriers'] for model in models] == [['tr'], ['ts']]
    for model in models:
        assert 'Reactant ' + model['reactant'] + '\n' in model['contents']
        assert 'ZeroEnergy[kcal/mol] 30' in model['contents']
    swapped = connected_models(text.replace('Reactant r', 'Reactant s'))
    assert swapped[0]['wells'] == ['s']


def test_inner_barrier_connects_requested_wells_and_input_is_unchanged():
    text = independent_wells().replace('End ! end kinetics',
        route('inner', 'r', 's', ('R', 'S', 'S', 'R')) + 'End ! end kinetics')
    assert connected_models(text)[0]['contents'] == text
    assert len(connected_models(text)) == 1


def test_parallel_paths_and_nested_mc_unions_are_kept_with_their_group():
    text = independent_wells().replace('End ! end kinetics',
        route('extra', 'r', 'b', ('R', 'S', 'B', 'B'), 'another') + 'End ! end kinetics')
    text = text.replace(' RRHO\n ZeroEnergy', ' Union\n RRHO\n ZeroEnergy').replace(
        ' 30\n End\n End\n', ' 30\n End\n End ! MC Union\n End\n')
    completed = complete_mirror_channels(text)
    models = connected_models(completed)
    assert len(models) == 2
    for model in models:
        assert 'Union ! 2 stereochemical pathways' in model['contents']
        assert model['contents'].count('End ! MC Union') == 2
        assert '! kinbot_stereopath another' in model['contents']


def test_racemic_fold_precedes_grouping_and_redirects_primary_reactant():
    from test_mess_racemates import well, barrier, network as racemic_network
    text = racemic_network(well('r', 'R', 'S'), well('s', 'S', 'R'), well('a', 'A', 'A'),
        barrier('tr', 'r', 'a', ('R', 'S', 'A', 'A')),
        barrier('ts', 's', 'a', ('S', 'R', 'A', 'A')), reactant='s')
    completed = complete_mirror_channels(text)
    models = connected_models(completed)
    assert len(models) == 1
    assert models[0]['reactant'] == 'r'
    assert models[0]['wells'] == ['r', 'a']
    assert len(models[0]['barriers']) == 1
    assert models[0]['contents'] == completed


@pytest.mark.parametrize('change,message', [
    (lambda s: s.replace('Reactant r', 'Reactant missing'), 'Reactant is not'),
    (lambda s: s.replace('Barrier tr r b', 'Barrier tr r missing'), 'missing endpoint'),
    (lambda s: s.replace('Barrier tr r b', 'Barrier tr b b'), 'no well endpoint'),
    (lambda s: s.replace('Well s', 'Well r'), 'Duplicate'),
])
def test_bad_network_is_not_silently_truncated(change, message):
    with pytest.raises(ValueError, match=message):
        connected_models(change(independent_wells()))


def test_uq_and_groups_have_distinct_files_and_current_manifest(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    writer = SimpleNamespace()
    for index in range(2):
        write_network_inputs(writer, independent_wells(), index)
    jobs = json.loads(Path('me/mess_networks.json').read_text())
    assert len(jobs) == 4
    assert len({job['stem'] for job in jobs}) == 4
    for job in jobs:
        assert job['status'] == 'ready'
        text = Path('me', job['stem'] + '.inp').read_text()
        assert f"MicroRateOutput {job['stem']}.micro" in text
    # Old outputs and inputs are preserved, but no longer scheduled by a glob.
    Path('me/mess_0000_group_0002.out').write_text('saved output')
    write_network_inputs(writer, connected_models(independent_wells())[0]['contents'], 0)
    assert len(writer.mess_jobs) == 3
    assert Path('me/mess_0000_group_0002.out').read_text() == 'saved output'


def test_empty_group_is_preserved_but_not_scheduled(tmp_path, monkeypatch, caplog):
    from kinbot.mess_execution import run_mess
    monkeypatch.chdir(tmp_path)
    writer = SimpleNamespace(par=dict(queuing='local', uq_n=1, run_me=True, uq_max_runs=1))
    write_network_inputs(writer, network(population('r', 'R', 'S')), 0)
    assert writer.mess_jobs[0]['status'] == 'no_reactions'
    assert 'no rate calculation will be submitted' in caplog.text
    assert run_mess(writer) == 0
    assert Path('batch_me.sub').read_text() == ''
    assert Path('me/mess_0000.inp').exists()


def test_actual_pes_assembly_preserves_disconnected_selected_wells(tmp_path, monkeypatch):
    from kinbot import pes
    from kinbot.parameters import Parameters
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(pes, 'logger', logging.getLogger('KinBot'), raising=False)
    Path('input.json').write_text(json.dumps(dict(barrier_threshold=100., smiles='O',
        high_level=0, me=0, epsilon=100., sigma=3., rotor_scan=0)))
    par = Parameters('input.json', show_warnings=False).par
    for name in ('123', '456'):
        Path(name).mkdir()
        Path(name, name + '_0000.mess').write_text(
            ' Well {name}\n RRHO\n ZeroEnergy[kcal/mol] {zeroenergy}\n End\n End\n')
    pes.create_mess_input(par, ['123', '456'], [], [], [], [],
        {'123': 0., '456': 5.}, {}, {'123': '123', '456': '456'}, 18., False)
    jobs = json.loads(Path('me/mess_networks.json').read_text())
    assert len(jobs) == 2
    assert all(job['status'] == 'no_reactions' for job in jobs)
    assert 'ZeroEnergy[kcal/mol] 0.0' in Path('me/mess_0000.inp').read_text()
    assert 'ZeroEnergy[kcal/mol] 5.0' in Path('me/mess_0000_group_0002.inp').read_text()


@pytest.mark.parametrize('mc', [0, 1])
def test_pes_passes_all_active_groups_and_uq_samples_to_runner(tmp_path, monkeypatch, mc):
    from kinbot import pes
    from kinbot.mess import MESS
    from kinbot.parameters import Parameters
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(pes, 'logger', logging.getLogger('KinBot'), raising=False)
    Path('input.json').write_text(json.dumps(dict(barrier_threshold=100., smiles='O',
        high_level=0, me=0, epsilon=100., sigma=3., rotor_scan=0,
        multi_conf_tst=mc, conformer_search=1)))
    par = Parameters('input.json', show_warnings=False).par
    par.update(uq_n=2, me=1)
    wells = ['123', '456', '789', '901']
    for name in wells:
        Path(name).mkdir()
        for index in range(2):
            Path(name, name + f'_{index:04d}.mess').write_text(
                ' Well {name}\n RRHO\n ZeroEnergy[kcal/mol] {zeroenergy}\n End\n End\n')
    for name, well in [('ta', '123'), ('tb', '789')]:
        for index in range(2):
            Path(well, name + f'_{index:04d}.mess').write_text(
                ' Barrier {name}\n RRHO\n ZeroEnergy[kcal/mol] {zeroenergy}\n End\n End\n')
    calls = []
    def run(writer):
        calls.append(writer.mess_jobs)
        assert len(writer.mess_jobs) == 4
        assert all(job['status'] == 'ready' for job in writer.mess_jobs)
        for job in writer.mess_jobs:
            text = Path('me', job['stem'] + '.inp').read_text()
            assert text.count(' Well ') == 2
            assert text.count(' Barrier ') == 1
            assert '{zeroenergy}' not in text
    monkeypatch.setattr(MESS, 'run', run)
    pes.create_mess_input(par, wells, [],
        [['123', 'ta', ['456'], 30.], ['789', 'tb', ['901'], 32.]], [], [],
        dict.fromkeys(wells, 0.), {}, {name: name for name in wells}, 18., False)
    assert len(calls) == 1
