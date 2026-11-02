import pytest
from kinbot.mess_branching import event_hit_probability as hit
from kinbot.mess_branching import channel_selection_groups


def test_explicit_mirrors_group_even_when_legacy_display_labels_differ():
    text = '''Model
 Well a ! starting
! kinbot_population A A
 End
 Well r ! resonance display one
! kinbot_population R S
 End
 Well s ! resonance display two
! kinbot_population S R
 End
 Barrier ar a r
 End
 Barrier as a s
 End
End ! end kinetics
'''
    assert list(channel_selection_groups(text).values()) == [['ar', 'as']]
    assert list(channel_selection_groups(text, population_families={
        'a': 'initial', 'r': 'display one', 's': 'display two'}).values()) == [['ar', 'as']]


def test_exact_union_does_not_sum_overlapping_ever_hit_events():
    rates = {('a', 'b'): 1, ('a', 'exit'): 1, ('b', 'a'): 1, ('b', 'exit'): 1}
    forward = {('a', 'b'): 1}; reverse = {('b', 'a'): 1}
    assert hit(['a', 'b'], rates, 'a', forward) == pytest.approx(.5)
    assert hit(['a', 'b'], rates, 'a', reverse) == pytest.approx(.25)
    assert hit(['a', 'b'], rates, 'a', forward | reverse) == pytest.approx(.5)


def test_mirror_events_group_before_threshold_and_keep_partial_parallel_flux():
    rates = {('a', 'r'): 1, ('a', 's'): 1, ('a', 'exit'): 198,
             ('r', 'exit'): 1, ('s', 'exit'): 1}
    assert hit(['a', 'r', 's'], rates, 'a', {('a', 'r'): 1}) == .005
    assert hit(['a', 'r', 's'], rates, 'a', {('a', 'r'): 1, ('a', 's'): 1}) == .01
    assert hit(['a', 'r', 's'], rates, 'a', {('a', 'r'): .5}) == .0025


def test_closed_traps_unreachable_events_and_well_clock_invariance():
    wells = ['a', 'b', 'c', 'unused']
    rates = {('a', 'exit'): 1, ('a', 'b'): 1, ('b', 'c'): 1,
             ('c', 'b'): 1, ('unused', 'exit'): 1}
    assert hit(wells, rates, 'a', {('a', 'exit'): 1}) == .5
    assert hit(wells, rates, 'a', {('unused', 'exit'): 1}) == 0
    scaled = {key: value*(1e20 if key[0]=='a' else .003) for key, value in rates.items()}
    assert hit(wells, scaled, 'a', {('a', 'exit'): 1e20}) == .5


def test_tiny_escape_and_invalid_event():
    rates = {('a', 'b'): 1, ('b', 'a'): 1, ('b', 'exit'): 1e-40}
    assert hit(['a', 'b'], rates, 'a', {('b', 'exit'): 1e-40}) == 1
    with pytest.raises(ValueError):
        hit(['a', 'b'], rates, 'a', {('a', 'b'): 2})
