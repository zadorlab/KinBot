from copy import deepcopy
import json
from types import SimpleNamespace

import pytest
from kinbot import barrier_policy as policy, constants
from kinbot.execution import atomic_json


def reference():
    return dict(zero_k_hartree=-100.,ceiling_kcal_mol=400.,source='accepted initial well')


def test_absolute_ceiling_does_not_accumulate_child_offsets(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    p=dict(pes_barrier_reference=reference(),high_level=1,barrier_threshold_add=10.)
    well=SimpleNamespace(energy=-100.+50./constants.AUtoKCAL,zpe=0.)
    policy.apply(well,p)
    assert p['barrier_threshold_L2']==pytest.approx(350.)
    assert p['barrier_threshold']==pytest.approx(360.)
    policy.apply(well,p)
    assert p['barrier_threshold_L2']==pytest.approx(350.)
    # A child starts with the same absolute zero, not its parent's local limit.
    policy.apply(SimpleNamespace(energy=-100.+80./constants.AUtoKCAL,zpe=0.),p)
    assert p['barrier_threshold_L2']==pytest.approx(320.)
    assert json.loads((tmp_path/'barrier_policy.json').read_text())['reference']==reference()


def test_existing_accepted_network_over_new_ceiling_blocks_before_mutation(tmp_path):
    atomic_json(tmp_path/'pes_exploration.json',dict(initial_worker='111'))
    atomic_json(tmp_path/'111/111.json',{})
    path=tmp_path/'111/anl_automatic/network_request.json'
    plan=dict(network=dict(stable_species=[dict(routing_name='111',energy_hartree=-100.,harmonic_zpe_hartree=0.)],
        transition_states={'ts':dict(energy_hartree=-100.+450./constants.AUtoKCAL,harmonic_zpe_hartree=0.)}))
    atomic_json(path,plan);saved=path.read_bytes()
    with pytest.raises(ValueError,match='explicit network pruning'):
        policy.for_case(tmp_path,400.)
    assert path.read_bytes()==saved
    assert policy.for_case(tmp_path,500.)['zero_k_hartree']==-100.


def test_legacy_parameters_without_reference_are_unchanged():
    parameters=dict(barrier_threshold=100.,high_level=0)
    prior=deepcopy(parameters);policy.apply(object(),parameters)
    assert parameters==prior
