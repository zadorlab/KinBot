"""Exercise the complete review campaign without a scheduler or native QC."""
from copy import deepcopy
import importlib.util
from pathlib import Path

import pytest
import numpy as np

from kinbot import run_policy
from kinbot.execution import atomic_json, file_sha256

HERE=Path(__file__).resolve().parents[1]/'examples/anl/ch3ss_anl1_f12/orbital_reviews/r9_completion_20261112'
spec=importlib.util.spec_from_file_location('r9_completion',HERE/'campaign.py')
campaign=importlib.util.module_from_spec(spec);spec.loader.exec_module(campaign)


@pytest.fixture
def staged(tmp_path,monkeypatch):
    plan=campaign.read(HERE/'plan.json')
    # Both electron counts and all mandatory old-channel checks remain in the test.
    for check in plan['existing_checks']:
        for b,seed in check['seeds'].items():
            for key,sha in (('path','sha256'),('output','output_sha256')):
                p=tmp_path/f'{check["case_id"]}_{b}_{key}';p.write_text('Synthetic native prerequisite '+str(p))
                seed[key]=str(p);seed[sha]=file_sha256(p)
    recipe=tmp_path/'recipe';recipe.mkdir();atomic_json(recipe/'plan.json',plan)
    monkeypatch.setattr(campaign,'HERE',recipe)
    pool=tmp_path/'pool';pool.mkdir();run_policy.update(pool,max_jobs=20,exclude=['node9'],fair_submission=False)
    output=pool/'completion'
    record=campaign.stage(output,pool=pool)
    return plan,pool,output,record


def test_offline_stage_and_repeat_preserve_sources_and_science(staged,monkeypatch):
    plan,pool,out,record=staged
    monkeypatch.setattr(campaign.recon.subprocess,'run',lambda *a,**k:pytest.fail('Offline staging queried the scheduler.'))
    before={p:p.read_bytes() for p in out.rglob('*.inp')}
    assert campaign.stage(out,pool=pool)==record
    assert all(p.read_bytes()==value for p,value in before.items())
    jobs=[j for m in record['manifests'] for j in campaign.read(m)['jobs']]
    assert len(jobs)==92
    assert {c['reference_space']['active_electrons'] for c in plan['h_loss_cases']}=={13,19}
    for job in jobs:
        if job['natural_orbital_probe']:
            assert all(job['probe_active_space']['active_electrons']%2 == (s['multiplicity']-1)%2 for s in job['root_states'])
    seed=plan['existing_checks'][0]['seeds']['avdz-f12']
    Path(seed['path']).write_text('Changed')
    with pytest.raises(ValueError,match='Changed original native evidence'):
        campaign.stage(pool/'changed',pool=pool)


def fake_partition(case, reference, atomic):
    total=sum(campaign.atomic_numbers[s] for s in case['symbols'])
    e,o=case['reference_space']['active_electrons'],case['reference_space']['active_orbitals']
    closed=(total-e)//2;occupied=(e-1)//2
    return dict(radical=closed+occupied+1,
                frontier=[closed+occupied+1],
                occupied=list(range(closed+1,closed+occupied+1)),
                virtual=list(range(closed+occupied+2,closed+o+1)),
                hydrogen_projection=1.,occupations=[2.]*occupied+[1.]+[0.]*(o-occupied-1),
                source_sha256='synthetic')


def test_full_owner_partial_completion_and_repeat_adopt_finished_work(staged,monkeypatch):
    plan,pool,out,record=staged
    finished=[];first=True
    def fake_many(paths,**options):
        nonlocal first
        assert options['max_jobs']==20 and options['keep_going']
        for path in paths:
            data=campaign.read(path);title=campaign.read(Path(path).parents[1]/'input.json')['title']
            for job in data['jobs']:
                if job.get('status')=='complete':continue
                d=Path(job['directory']);output=d/(job['name']+'.out');wfu=d/(job['name']+'.wfu')
                output.write_text('Simulated native result');wfu.write_text('Simulated checkpoint '+str(wfu))
                # Deterministic, finite energies exercise own-asymptote ranking.
                level=-100.+(.001 if job['anchor']!='r30' else 0.)
                if '_pair_' in title:level+=int(title.split('_pair_')[1].split('_')[0])*1.e-7
                job.update(status='complete',output_sha256=file_sha256(output),wfu=str(wfu),wfu_sha256=file_sha256(wfu),
                           correlated_energies_hartree=[dict(method='caspt2-f12',multiplicity=2,root=1,energy_hartree=level)])
                finished.append((path,job['name']));atomic_json(path,data)
                if first:
                    first=False
                    raise RuntimeError('Simulated controller interruption after native completion')
            atomic_json(path,data)
    monkeypatch.setattr(campaign.recon,'run_many',fake_many)
    monkeypatch.setattr(campaign,'hydrogen_partition',fake_partition)
    monkeypatch.setattr(campaign,'distinct_pair_trials',lambda case,scores,jobs,bases:(
        sorted((r for r in scores if r['name']!='single'),key=lambda r:r['name'])[:3],{}))
    with pytest.raises(RuntimeError,match='interruption'):
        campaign.run(out)
    campaign.run(out)
    result=campaign.read(out/'completion.json')
    assert result['status']=='native_complete_review_required' and not result['production_model_accepted']
    assert not result['holds']
    before=list(finished)
    campaign.run(out)
    assert finished==before
    assert len(set(finished))==len(finished)
    for cid in (101,102):
        chosen=campaign.read(out/f'h{cid}_followup_trials.json')
        assert not chosen['production_model_accepted']
        assert set(chosen['selected']) >= {'single','union_2','union_3'}
    assert all(campaign.read(p)['shared_policy_file']==str(pool/run_policy.NAME) for p in result['manifests'])


def test_failed_native_prerequisite_does_not_create_guessed_partitions(staged,monkeypatch):
    plan,pool,out,record=staged
    monkeypatch.setattr(campaign,'run_manifests',lambda manifests,cap:'Native prerequisites incomplete')
    campaign.run(out)
    result=campaign.read(out/'completion.json')
    assert result['status']=='native_partial_review_required'
    assert len(result['holds'])==3
    assert not (out/'screen').exists()
    assert not result['production_model_accepted']


def test_equivalent_native_spaces_are_not_repeated_on_full_paths(tmp_path,monkeypatch):
    from kinbot.radical_pair_evidence import Orbitals
    case=campaign.read(HERE/'plan.json')['h_loss_cases'][0]
    bases=['avdz-f12','avtz-f12'];jobs={}
    for name in ('pair_0_0','pair_1_0','pair_2_0'):
        for basis in bases:
            for point in ('near','mid','r30'):
                title=f'h101_{name}_{basis.replace("-","_")}'
                jobs[(title,basis,point)]=dict(directory=str(tmp_path/name),name=point,
                    probe_active_space=dict(active_electrons=3,active_orbitals=3),
                    correlated_energies_hartree=[dict(method='caspt2-f12',multiplicity=2,root=1,energy_hartree=-10.)])
    def orbitals(path):
        matrix=np.eye(15)
        if 'pair_2_0' in str(path):
            angle=.1
            matrix[np.ix_([11,14],[11,14])]=[[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]]
        return Orbitals(matrix,np.ones(15),['H'],np.zeros((1,3)),['same basis'],[list(range(15))])
    monkeypatch.setattr(campaign,'read_orbitals',orbitals)
    scores=[dict(name=name,max_potential_difference=0.) for name in ('pair_0_0','pair_1_0','pair_2_0')]
    chosen,equivalent=campaign.distinct_pair_trials(case,scores,jobs,bases)
    assert [r['name'] for r in chosen]==['pair_0_0','pair_2_0']
    assert equivalent=={'pair_1_0':'pair_0_0'}


def test_hydrogen_identity_uses_native_orbital_projection(tmp_path,monkeypatch):
    from kinbot.radical_pair_evidence import Orbitals
    case=campaign.read(HERE/'plan.json')['h_loss_cases'][0]
    atom=case['coordinate']['moving_atoms'][0]-1
    positions=np.array(case['geometry']);columns=[[0,1,2,3,4],[5,6,7,8,9],[12],[10,11],[13,14]]
    assert atom==2
    occupations=np.array([2.]*12+[1.,0.,0.])
    dimer=Orbitals(np.eye(15),occupations,case['symbols'],positions,['heavy','heavy','H basis','H basis','H basis'],columns)
    hydrogen=Orbitals(np.ones((1,1)),np.ones(1),['H'],positions[atom:atom+1],['H basis'],[[0]])
    directory=tmp_path/'reference';directory.mkdir();(directory/'far.natural.orbitals.xml').write_text('Synthetic native matrix fixture')
    monkeypatch.setattr(campaign,'read_orbitals',lambda path:hydrogen if str(path).endswith('.scf.orbitals.xml') else dimer)
    result=campaign.hydrogen_partition(case,dict(directory=str(directory),name='far'),dict(directory=str(tmp_path),name='H'))
    assert result['radical']==13 and result['hydrogen_projection']==pytest.approx(1.)
    assert len(result['occupied'])==6 and len(result['virtual'])==2
    # A rotation in a degenerate singly occupied span must preserve the H test.
    dimer.occupations[10:13]=1.;dimer.occupations[13]=2.
    angle=np.pi/4
    dimer.coefficients[np.ix_([11,12],[11,12])]=[[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]]
    rotated=campaign.hydrogen_partition(case,dict(directory=str(directory),name='far'),dict(directory=str(tmp_path),name='H'))
    assert rotated['frontier']==[11,12,13]
    assert rotated['hydrogen_projection']==pytest.approx(1.)
    assert rotated['hydrogen_electrons']==pytest.approx(1.)
    hydrogen.occupations[:]=2.
    with pytest.raises(ValueError,match='one electron'):
        campaign.hydrogen_partition(case,dict(directory=str(directory),name='far'),dict(directory=str(tmp_path),name='H'))


def test_text_report_reads_waiting_campaign_without_changes(staged,capsys):
    plan,pool,out,record=staged
    spec=importlib.util.spec_from_file_location('r9_text_report',HERE/'report.py')
    report=importlib.util.module_from_spec(spec);spec.loader.exec_module(report)
    before={p:file_sha256(p) for p in out.rglob('*') if p.is_file()}
    assert report.report(out)==0
    text=capsys.readouterr().out
    assert '"waiting": 92' in text and '"production_model_accepted": false' in text
    assert {p:file_sha256(p) for p in out.rglob('*') if p.is_file()}==before
