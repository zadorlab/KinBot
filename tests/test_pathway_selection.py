"""Placeholder dissociation must not erase or block stereochemical saddles."""
import itertools

import pytest
from unittest.mock import patch
from pathlib import Path
import copy

from kinbot.reaction_path import compare_pathways
from kinbot.pes import select_summary_reaction, remove_unused_complexes


def test_saddle_preference_precedes_stereo_comparison():
    assert compare_pathways('hom_sci_2_10', None, 30., 'saddle_H10', 'H10', 35.) == 'replace'
    assert compare_pathways('saddle_H10', 'H10', 35., 'hom_sci_2_10', None, 30.) == 'keep'
    assert compare_pathways('saddle_H10', 'H10', 35., 'saddle_H11', 'H11', 32.) == 'distinct'
    assert compare_pathways('saddle_H10', 'H10', 35., 'ordinary_unknown', None, 30.) == 'keep'
    assert compare_pathways('ordinary_unknown', None, 30., 'saddle_H10', 'H10', 35.) == 'replace'


def test_invalid_selected_path_is_rejected_locally(caplog):
    from types import SimpleNamespace
    from kinbot.reaction_path import reject_invalid_pathway
    reactions = [SimpleNamespace(instance_name='bad'), SimpleNamespace(instance_name='good')]
    species = SimpleNamespace(reac_obj=reactions, reac_ts_done=[-1, -1], reac_type=['test', 'test'])
    with patch('kinbot.reaction_path.reaction_path_id', side_effect=ValueError('changed configuration')):
        assert reject_invalid_pathway(species, 0)
    assert species.reac_ts_done == [-999, -1]
    assert 'network is incomplete' in caplog.text


def test_unclassified_route_cannot_merge_two_known_classes():
    rows = [['R', 'a', ['P'], 32.], ['R', 'b', ['P'], 30.], ['R', 'unknown', ['P'], 1.]]
    for order in itertools.permutations(rows):
        retained = []
        for row in order:
            select_summary_reaction(retained, row, {'a': 'stereo-a', 'b': 'stereo-b'})
        assert sorted(row[1] for row in retained) == ['a', 'b']


def test_pes_retains_both_stereochemical_classes_in_every_order():
    rows = [['R', 'hom_sci_2_10', ['P', 'H'], 30.],
            ['R', 'hom_sci_2_11', ['P', 'H'], 31.],
            ['R', 'saddle_H10', ['P', 'H'], 35.],
            ['R', 'saddle_H11', ['P', 'H'], 32.],
            ['R', 'saddle_H10_higher', ['P', 'H'], 36.]]
    paths = {'saddle_H10': 'H10', 'saddle_H10_higher': 'H10', 'saddle_H11': 'H11'}
    for order in itertools.permutations(rows):
        reactions = []
        for row in order:
            select_summary_reaction(reactions, row, paths)
        assert sorted(row[1] for row in reactions) == ['saddle_H10', 'saddle_H11']


def test_ordinary_minimum_reverse_direction_and_complex_tie():
    reactions = [['R', 'higher', ['P'], 35.]]
    candidate = ['P', 'lower', ['R'], 30.]
    select_summary_reaction(reactions, candidate, {})
    assert reactions == [candidate]
    complex_route = ['P', 'equal_with_complex', ['R'], 30., '12', 'vdW_IRC_F_prod']
    select_summary_reaction(reactions, complex_route, {})
    assert reactions == [complex_route]
    select_summary_reaction(reactions, ['R', 'equal_without_complex', ['P'], 30.], {})
    assert reactions == [complex_route]


def test_rejected_complex_cleanup_removes_exact_name_not_last_well():
    old = ['R', 'high', ['P', 'H'], 35., '10', 'vdW_IRC_F_prod']
    new = ['R', 'low', ['P', 'H'], 30.]
    reactions = [old]
    discarded = select_summary_reaction(reactions, new, {})
    wells = ['R', 'high_IRC_F_prod', 'unrelated']
    flags = [False, True, False]
    parents = {'R': 'R', 'high_IRC_F_prod': 'R', 'H_P': 'high_IRC_F_prod', 'unrelated': 'R'}
    remove_unused_complexes(discarded, reactions, wells, flags, parents)
    assert reactions == [new]
    assert wells == ['R', 'unrelated']
    assert flags == [False, False]
    assert parents['H_P'] == 'R'
    assert 'high_IRC_F_prod' not in parents


def test_complex_still_used_by_another_route_is_preserved():
    old = ['R', 'saddle', ['P', 'H'], 35., '10', 'vdW_IRC_F_prod']
    reactions = [old]
    wells = ['R', 'saddle_IRC_F_prod']
    flags = [False, True]
    parents = {'saddle_IRC_F_prod': 'R'}
    remove_unused_complexes([old], reactions, wells, flags, parents)
    assert wells == ['R', 'saddle_IRC_F_prod']
    assert flags == [False, True]
    assert parents == {'saddle_IRC_F_prod': 'R'}


def test_direct_mess_keeps_both_saddle_classes_and_no_placeholder(tmp_path, monkeypatch):
    from test_complex_mess import TestComplexMESS
    from kinbot.mess import MESS
    from kinbot.parameters import Parameters
    fixture = TestComplexMESS()
    parent, saddles = fixture.reactions()
    for reaction in saddles:
        reaction.do_vdW = False
    placeholder = copy.copy(saddles[0])
    placeholder.instance_name = 'hom_sci_2_10'
    placeholder.ts = copy.copy(parent)
    for index, order in enumerate(itertools.permutations([placeholder, *saddles])):
        directory = tmp_path/str(index)
        directory.mkdir()
        (directory/'me').mkdir()
        monkeypatch.chdir(directory)
        parent.reac_obj = order
        parent.reac_ts_done = [-1]*3
        parent.reac_type = ['hom_sci' if r is placeholder else 'test' for r in order]
        Path('input.json').write_text('{"barrier_threshold": 100.0}')
        par = Parameters('input.json', show_warnings=False).par
        par.update(pes=0, multi_conf_tst=0, rotor_scan=0, epsilon=100., sigma=3.)
        with patch('kinbot.mess.reaction_path_id', side_effect=lambda r: None if r is placeholder else r.instance_name):
            MESS(par, parent).write_input(None)
        text = Path('me/mess_0000.inp').read_text()
        assert 'Union ! 2 stereochemical pathways' in text
        assert 'Core PhaseSpaceTheory' not in text


@pytest.mark.parametrize('ordinary_energy,stereo_energy', [(79.50, 92.29), (92.29, 79.50)])
def test_examined_ordinary_and_stereo_keep_classes_in_every_order(ordinary_energy, stereo_energy):
    rows = [['R', 'ordinary', ['P'], ordinary_energy],
            ['R', 'ordinary_higher', ['P'], ordinary_energy + 1.],
            ['R', 'relay', ['P'], stereo_energy],
            ['R', 'hom_sci_1_2', ['P'], 0.]]
    paths = {'ordinary': 'ordinary', 'ordinary_higher': 'ordinary', 'relay': 'stereopath:relay'}
    for order in itertools.permutations(rows):
        selected = []
        for row in order:
            select_summary_reaction(selected, row, paths)
        assert sorted(row[1] for row in selected) == ['ordinary', 'relay']
    assert compare_pathways('ordinary', 'ordinary', ordinary_energy, 'legacy', None, 0.) == 'keep'


def test_saved_butanol_elimination_and_relay_are_examined_classes():
    import json
    import numpy as np
    from types import SimpleNamespace
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.reaction_path import prepare_stereopath, reaction_path_id, summary_path_line, read_summary_paths
    data = json.loads((Path(__file__).parent/'reference/butanol_pathways.json').read_text())
    parent = StationaryPoint('butanol', 0, 1, atom=data['atoms'], geom=np.array(data['reactant_geometry']))
    parent.characterize()
    paths = []
    for route in data['routes']:
        product = StationaryPoint('endpoint', 0, 1, atom=data['atoms'], geom=np.array(route['endpoint_geometry']))
        product.characterize()
        ts = StationaryPoint(route['name'], 0, 1, atom=data['atoms'], geom=np.array(route['ts_geometry']), wellorts=1)
        prepare_stereopath(ts, parent, product)
        reaction = SimpleNamespace(ts=ts, instance_name=ts.name)
        path = reaction_path_id(reaction)
        assert read_summary_paths([summary_path_line(reaction)])[ts.name] == path
        paths.append(path)
    assert paths[0] == 'ordinary'
    assert paths[1] == 'htransfer:5ea045c6ffc781763ef0e71d57580c598ed6bf6a5a8d64baeecbcd98155b47b9'
    assert compare_pathways('ordinary', paths[0], 79.50, 'relay', paths[1], 92.29) == 'distinct'
