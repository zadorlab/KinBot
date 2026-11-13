"""Alternative beta-scission guesses through the stopped-owner retry path."""
from copy import deepcopy
import ast
import json
from pathlib import Path
from types import SimpleNamespace as NS
from unittest.mock import patch

from ase import Atoms
from ase.db import connect
import numpy as np
import pytest

from kinbot import network_revisit, reaction_guess, run_refresh, search_control
from kinbot.execution import atomic_json
from kinbot.reactions.reac_R_Addition_MultipleBond import RAdditionMultipleBond
from tests.test_run_refresh import live_tree, sampling_tree, scientific_files


NAME = '123_R_Addition_MultipleBond_1_2_3'
GUESS = dict(bond_scale=1.6, angle='parent')


def species():
    bond = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    return NS(atom=np.array(['C', 'S', 'H']), natom=3, linear=[],
              geom=np.array([[0., 0., 0.], [1.8, 0., 0.], [2.1, 1.3, 0.]]),
              bond=bond, bonds=[bond])


def test_only_selected_guess_changes_and_final_ts_is_released():
    s = species()
    par = dict(reaction_search_guesses={NAME: GUESS})
    reaction = RAdditionMultipleBond(s, None, par, [0, 1, 2], reaction_guess.alternative_name(NAME, GUESS))
    distance, angle = reaction_guess.targets(reaction)
    assert distance == pytest.approx(np.linalg.norm(s.geom[2]-s.geom[1])*1.6)
    assert angle == pytest.approx(102.9946167919)
    _, _, changes, _ = reaction.get_constraints(0, s.geom.copy())
    assert changes[0] == [1, 2, 3, angle]
    assert reaction.get_constraints(reaction.max_step, s.geom.copy())[2] == []
    other = RAdditionMultipleBond(s, None, par, [0, 1, 2], 'other_R_Addition_MultipleBond_1_2_3')
    assert reaction_guess.targets(other) is None
    assert other.get_constraints(0, s.geom.copy())[2][0][-1] == 90.


@pytest.mark.parametrize('value', [None, [], {'bad': GUESS}, {NAME: dict(GUESS, bond_scale=float('nan'))},
                                  {NAME: dict(GUESS, bond_scale=True)}, {NAME: dict(GUESS, angle=90)}])
def test_invalid_guess_rejected(value):
    with pytest.raises(ValueError):
        reaction_guess.validate(value)


def test_guess_requires_exact_retry_scope():
    with pytest.raises(ValueError, match='independent'):
        network_revisit.settings(rejected=True, reaction_guesses={NAME: GUESS})
    with pytest.raises(ValueError, match='independent'):
        network_revisit.settings(rejected=['another'], reaction_guesses={NAME: GUESS})


def test_multiple_guesses_have_distinct_names_and_parent_geometry_targets():
    variants = [dict(GUESS, bond_scale=x) for x in (1.3, 1.6, 1.9)]
    par = dict(reaction_search_guesses={NAME: variants})
    reaction_guess.validate(par['reaction_search_guesses'])
    s = species()
    for options in variants:
        r = RAdditionMultipleBond(s, None, par, [0, 1, 2], reaction_guess.alternative_name(NAME, options))
        assert reaction_guess.targets(r)[0] == pytest.approx(np.hypot(.3, 1.3)*options['bond_scale'])
    with pytest.raises(ValueError, match='distinct'):
        reaction_guess.validate({NAME: [GUESS, GUESS]})


@pytest.mark.parametrize('source_status', [-1, -999])
def test_refresh_owner_retry_and_repeat_preserve_other_work(sampling_tree, monkeypatch, source_status):
    t = sampling_tree
    worker = t.root/'ch3ss/123'
    for case in ('ch3s', 'ch3ss'):
        w = t.root/case/'123'
        atomic_json(w/'123.json', dict(high_level=1, profiled_l1_reuse_only=0))
        rejected = NAME if case == 'ch3ss' else 'unrelated'
        (w/'kinbot_monitor.out').write_text(f'-1 13 accepted\n{source_status} 13 {rejected}\n')
    original, jobs = scientific_files(t.root), deepcopy(t.jobs)
    options = network_revisit.settings(reaction_guesses={NAME: GUESS})
    result = run_refresh.refresh_when_idle(t.root, apply=True, wait_idle=30, network_options=options)
    assert result['phase'] == 'launched'
    assert t.jobs == jobs and scientific_files(t.root) == original
    assert {p['pid'] for p in t.processes} == {700, 701, 702}
    assert run_refresh.refresh(t.root, apply=True, network_options=options) == result
    assert len(t.launches) == 2
    par = json.loads((worker/'123.json').read_text())
    assert par['high_level'] == 1 and par['reaction_search_guesses'] == {NAME: GUESS}
    assert json.loads((t.root/'ch3s/123/.kinbot_revisit.json').read_text())['rejected'] == []
    monkeypatch.chdir(worker)
    qc = NS(db=connect('retry-evidence.db'), job_ids={})
    for name in (NAME, 'accepted', 'product_well'):
        qc.db.write(Atoms('H'), name=name, data={'status': 'normal'})
    (worker/(NAME+'.log')).write_text('original native rejected saddle')
    s = species()
    s.reac_obj = [RAdditionMultipleBond(s, qc, par, [0, 1, 2], n) for n in (NAME, 'accepted')]
    s.reac_name, s.reac_type = [NAME, 'accepted'], ['R_Addition_MultipleBond']*2
    s.reac_inst, s.reac_scan_energy = [[0, 1, 2]]*2, [[], []]
    s.reac_ts_geom, s.reac_ts_freq = [0, 0], [0, 0]
    s.reac_ts_done, s.reac_step = [source_status, -1], [13, 13]
    g = NS(species=s, qc=qc, par=par, inp='123.json')
    reaction_guess.add_alternatives(g)
    search_control.reactions(g)
    assert s.reac_ts_done == [source_status, -1, 0]
    assert reaction_guess.targets(s.reac_obj[0]) is None
    assert reaction_guess.targets(s.reac_obj[2])[0] == pytest.approx(1.6*np.hypot(.3, 1.3))
    assert (worker/(NAME+'.log')).read_text() == 'original native rejected saddle'
    assert len(list(qc.db.select(name='product_well'))) == 1
    s.reac_ts_done[2] = -999
    count = qc.db.count()
    reaction_guess.add_alternatives(g)
    search_control.reactions(g)
    assert s.reac_ts_done == [source_status, -1, -999] and qc.db.count() == count
    assert not network_revisit.reassesses_accepted(options)
    from kinbot.search_state import SearchState
    checkpoint = SearchState(g)
    s.reac_ts_done[2] = -999
    checkpoint.save(2)
    s.reac_ts_done[2] = 0
    SearchState(g)
    assert s.reac_ts_done[2] == -999


def test_cli_carries_hash_pinned_guess_payload(tmp_path, monkeypatch):
    path = tmp_path/'guesses.json'
    atomic_json(path, {NAME: GUESS})
    called = []
    monkeypatch.setattr(run_refresh, 'refresh_when_idle', lambda root, **kw:
                        called.append(kw) or dict(phase='reviewed', jobs=[]))
    run_refresh.main([str(tmp_path), '--reaction-guesses', str(path)])
    assert called[0]['network_options']['reaction_search_guesses'] == {NAME: GUESS}


def test_normal_entry_renders_changed_search_then_unconstrained_ts(tmp_path, monkeypatch):
    from kinbot.parameters import Parameters
    from kinbot.qc import QuantumChemistry
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.reac_family import carry_out_reaction
    monkeypatch.chdir(tmp_path)
    s = StationaryPoint('test', 0, 2, smiles='[CH2]S')
    s.characterize()
    sulfur = list(s.atom).index('S')
    carbon = list(s.atom).index('C')
    hydrogen = next(i for i, symbol in enumerate(s.atom) if symbol == 'H' and s.bond[sulfur][i])
    name = f'test_R_Addition_MultipleBond_{carbon+1}_{sulfur+1}_{hydrogen+1}'
    atomic_json('input.json', dict(theory_preset='b3lyp-b2plyp-anl', barrier_threshold=150,
                                reaction_search_guesses={name: GUESS}))
    par = Parameters('input.json', show_warnings=False).par
    qc = QuantumChemistry(par)
    name = reaction_guess.alternative_name(name, GUESS)
    r = RAdditionMultipleBond(s, qc, par, [carbon, sulfur, hydrogen], name)
    with patch.object(qc.l1, 'submit_qc', return_value=1):
        assert carry_out_reaction(r, 0, 'g16') == 1
        text = Path(name+'.py').read_text()
        assert 'order=0' in text
        assert 'B3LYP' in text and 'cc-pVDZ' in text
        with patch.object(qc, 'check_qc', return_value='running'):
            assert carry_out_reaction(r, 0, 'g16') == 0
            assert Path(name+'.py').read_text() == text
        with patch.object(qc, 'check_qc', return_value='normal'), \
                patch.object(qc, 'get_qc_geom', return_value=(0, s.geom.copy())):
            assert carry_out_reaction(r, r.max_step, 'g16') == r.max_step+1
    final = Path(name+'.py').read_text()
    ast.parse(final)
    assert 'order=1' in final and 'frequency_validation_error(freqs, 1' in final
    assert 'constraints=' not in final
