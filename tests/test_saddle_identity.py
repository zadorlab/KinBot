"""Native distinct/duplicate saddles through local and outer-PES selection."""
from copy import deepcopy
from pathlib import Path
import itertools
import json
import numpy as np

from kinbot import saddle_identity as identity
from kinbot.pes import select_summary_reaction
from kinbot.mess import union_stereochemical_barriers
from kinbot.mess_mirrors import complete_mirror_channels, _reflect
from kinbot.mess_networks import split_model


def records():
    return json.loads((Path(__file__).parent/'reference/parallel_saddles.json').read_text())


def test_native_saddles_survive_all_discovery_orders_and_reverse_direction():
    saddles = records()
    assert identity.compare(saddles['loose'], saddles['repeat'])['status'] == 'equivalent'
    assert identity.compare(saddles['loose'], saddles['compact'])['status'] == 'distinct'
    rows = [['R', 'loose', ['P'], 42.41], ['P', 'repeat', ['R'], 42.41],
            ['R', 'compact', ['P'], 62.04], ['R', 'hom_sci_1_2', ['P'], 30.]]
    for order in itertools.permutations(rows):
        selected, audit = [], []
        for row in order:
            select_summary_reaction(selected, row, dict.fromkeys(saddles, 'ordinary'), saddles, audit)
        names = {r[1] for r in selected}
        assert 'compact' in names and len(names) == 2
        assert bool(names & {'loose', 'repeat'})
        assert any(r['status'] == 'equivalent' for r in audit)


def test_rotated_permuted_reverse_mode_is_duplicate():
    a = records()['loose']; b = deepcopy(a)
    order = [2, 1, 0, 5, 4, 3]
    rotation = np.array([[0.,-1,0],[1,0,0],[0,0,1]])
    for key in ('positions', 'mode'):
        b[key] = (np.asarray(b[key])[order] @ rotation).tolist()
    b['positions'] = (np.asarray(b['positions'])+[3,4,5]).tolist()
    b['mode'] = (-np.asarray(b['mode'])).tolist()
    b['symbols'] = [b['symbols'][i] for i in order]
    assert identity.compare(a,b)['status'] == 'equivalent'
    b['mode'] = np.roll(b['mode'],1,axis=0).tolist()
    assert identity.compare(a,b)['status'] == 'distinct'


def test_missing_mode_and_distinct_states_never_remove_a_saddle():
    a = records()['loose']; b = deepcopy(a); b.pop('mode')
    assert identity.compare(a,b)['status'] == 'unresolved'
    b['multiplicity'] = 4
    assert identity.compare(a,b)['status'] == 'distinct'
    b = deepcopy(a); b['isotopes'][0] = 13
    assert identity.compare(a,b)['status'] != 'equivalent'


def test_symmetric_geometry_tries_alternative_maps_for_the_reactive_mode():
    a = dict(symbols=['C','H','H','O'],isotopes=[0]*4,charge=0,multiplicity=1,
        positions=[[0.,0,0],[1,0,0],[-1,0,0],[0,0,2]],
        mode=[[.2,.3,.4],[1,.2,0],[0,1,.4],[.1,.1,-.5]],
        frequencies=[-100.,200.,300.,400.,500.,600.],energy_hartree=-100.)
    b = deepcopy(a); order=[0,2,1,3]
    for key in ('positions','mode'):
        b[key]=np.asarray(a[key])[order].tolist()
    result=identity.compare(a,b)
    assert result['status']=='equivalent'
    assert result['measurements']['mode_overlap'] > .999999


def test_duplicate_ensemble_requires_every_member_and_preserves_missing_evidence():
    a=records()['loose']; b=deepcopy(a)
    a['members']=[records()['loose'],records()['compact']]
    b['members']=list(reversed(deepcopy(a['members'])))
    assert identity.compare(a,b)['status']=='equivalent'
    b['members'][0].pop('mode')
    assert identity.compare(a,b)['status']=='unresolved'


def test_mapping_limit_is_recorded_as_unresolved(monkeypatch):
    def exhausted(*args,**kwargs):
        raise identity.MappingLimitReached('limit')
    monkeypatch.setattr(identity,'aligned_permutations',exhausted)
    record=records()['loose']
    result=identity.compare(record,record)
    assert result['status']=='unresolved' and result['reason']=='atom_mapping_limit'


def barrier(name, record, left='a', right='r', keys='A A R S'):
    block = (f'! kinbot_stereopath ordinary\n Barrier {name} {left} {right}\n'
             f'! kinbot_mirror_endpoints {keys}\n RRHO\n ZeroEnergy[kcal/mol] 42\n End\n')
    return identity.annotate(block,record)


def test_union_keeps_native_parallel_surfaces_and_deduplicates_repeat():
    r = records()
    blocks = [barrier(n,e) for n,e in r.items()]
    union = ''.join(union_stereochemical_barriers(blocks))
    assert 'Union ! 2 stereochemical pathways' in union
    assert 'Barrier compact ' in union and 'Barrier loose ' in union
    assert 'Barrier repeat ' not in union


def test_mirrors_preserve_two_surfaces_and_repeated_completion_is_stable():
    r = records()
    text = ('Model\n Well a\n! kinbot_population A A\n End\n'
            ' Well r\n! kinbot_population R S\n End\n'
            ' Well s\n! kinbot_population S R\n End\n')
    text += barrier('loose',r['loose']) + barrier('compact',r['compact'])
    text += _reflect(barrier('repeat',r['repeat'])).replace('Barrier repeat a r','Barrier repeat a s')
    result = complete_mirror_channels(text+'End ! end kinetics\n')
    assert result.count('Union ! 2 stereochemical pathways') == 2
    assert complete_mirror_channels(result) == result
    assert sum('pathway source: Barrier compact ' in b or 'compact' in b
               for b in split_model(result)[1]) >= 2


def test_outer_selection_uses_native_evidence_without_modifying_snapshots():
    from kinbot.pes_selection import choose
    r = records()
    data = dict(errors=[], network=dict(stable_species=[], reactions=[], transition_states={}))
    for name, evidence in r.items():
        data['network']['reactions'].append(dict(name=name, reactant='R', products=['P'],
            barrierless=False, transition_state=name, path_class='ordinary', saddle_evidence=evidence))
        data['network']['transition_states'][name] = dict(energy_hartree=evidence['energy_hartree'],
            harmonic_zpe_hartree=.03, l2_theory={'method':'B2PLYP', 'basis':'cc-pVTZ'})
    before = deepcopy(data)
    result = choose({'worker':data})
    assert result['workers']['worker']['selected'] == ['compact','loose']
    assert result['workers']['worker']['omitted']['repeat']['retained'] == ['loose']
    assert result == choose({'worker':data}) and data == before
