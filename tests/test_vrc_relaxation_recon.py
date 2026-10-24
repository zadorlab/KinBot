import importlib.util
from pathlib import Path
import json
import numpy as np
import pytest
from kinbot.execution import atomic_json
from tests.test_requested_orbital_recon import request_tree, completed_screen, RECIPE, SELECTION, helper as screening

ROOT=Path(__file__).resolve().parents[1]
from kinbot import vrc_relaxation_recon as helper


def test_relaxed_panel_uses_verified_seeds_and_one_constraint(request_tree):
    request,case,output=request_tree
    completed_screen(request_tree,RECIPE)
    par=dict(queue_name='test-partition',vrc_tst_scan_basis='avdz-f12',vrc_tst_scan_shift=.1,
        vrc_tst_scan_ipea=.25,vrc_tst_scan_gradient='numerical',vrc_tst_scf_convergence=1.e-8,
        vrc_tst_correlation_convergence=1.e-7,vrc_tst_scan_fmax=.0005,barrier_threshold=400.,
        smiles='C=[SH]',mult=2)
    atomic_json(case,par)
    before={p:p.read_bytes() for p in output.rglob('*') if p.is_file()}
    specs=helper.specifications(output,request,SELECTION,case)
    assert len(specs)==4
    assert {s['intent']['relaxation_recon']['distance_angstrom'] for s in specs}=={2.,3.5,6.,30.}
    for s in specs:
        proof=s['intent']['relaxation_recon'];t=s['tasks'][0];xyz=np.array(s['molecule']['positions']);i,j=proof['coordinate']
        assert np.linalg.norm(xyz[i]-xyz[j])==pytest.approx(proof['distance_angstrom'])
        assert t['optimizer']['constraints']==dict(fixed_bonds=[[i,j]])
        assert t['profile']['method']=='caspt2-f12(5,4)'
        assert s['intent']['production_model_accepted'] is False
        assert t['optimizer']['steps']==50 and t['optimizer']['step_retries']==0
    assert all(p.read_bytes()==b for p,b in before.items())
    # Rejected channels must not consume diagnostic or production allocations.
    atomic_json(request.parents[3]/'.kinbot_reaction_decisions.json',dict(schema=1,reactions={
        json.loads(request.read_text())['reaction']:dict(decision='reject')}))
    with pytest.raises(ValueError,match='active review'):
        helper.specifications(output,request,SELECTION,case)


@pytest.mark.parametrize('bond',['cs','ss'])
def test_actual_preliminary_channels_stage_relaxations_and_repeat_without_duplicates(tmp_path,bond):
    base=ROOT/'examples/anl/orbital_review_batches/r7_additional_20261021'
    recipe=base/(f'ch2shs_{bond}.json');selection=base/(f'ch2shs_{bond}_selected.json')
    match=json.loads(recipe.read_text())['request_match']
    positions=[[1.2000100119,-.4655510032,.1013882575],[.9455572501,.7518250541,1.3204837046],
        [.2070132629,-.3355235642,2.8547670677],[-.134624981,1.3246843886,.7341794664],
        [.3345123258,-.9791964859,-.2862124727],[2.1494724455,-.9711253636,.1676844543]]
    production=tmp_path/'production';case=production/'case/input.json'
    request=case.parent/match['parent']/'vrctst/orbital_review'/match['reaction']/'request.json'
    atomic_json(request,dict(match,schema=1,status='needs_review',positions=positions))
    atomic_json(case,dict(queue_name='test-partition'))
    atomic_json(production/'.kinbot_resource_policy.json',dict(schema=1,max_jobs=22,excluded_nodes=[],manage_pending=True,orbital_review_reserve=3))
    source=completed_screen((request,case,tmp_path/'screen'),recipe)
    par=dict(queue_name='test-partition',vrc_tst_scan_basis='avdz-f12',vrc_tst_scan_shift=.1,
        vrc_tst_scan_ipea=.25,vrc_tst_scan_gradient='numerical',vrc_tst_scf_convergence=1.e-8,
        vrc_tst_correlation_convergence=1.e-7,vrc_tst_scan_fmax=.0005,barrier_threshold=400.,
        smiles='C=[SH][S]',mult=2)
    atomic_json(case,par)
    plan=tmp_path/'plan.json';atomic_json(plan,dict(panels=[dict(source=str(source),request=str(request),
        selection=str(selection),case_input=str(case))]))
    first=helper.prepare(production,plan,tmp_path/'relaxations')
    second=helper.prepare(production,plan,tmp_path/'relaxations')
    assert first==second and len(first['graphs'])==4
    for g in first['graphs']:
        spec=helper.read(Path(g)/'workflow.json')
        assert spec['tasks'][0]['profile']['method']=='caspt2-f12(5,4)'
        assert spec['tasks'][0]['optimizer']['constraints']['fixed_bonds']==[[i-1 for i in match['bond']]]
