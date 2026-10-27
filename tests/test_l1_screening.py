from copy import deepcopy
from dataclasses import asdict
import json
from types import SimpleNamespace

from ase import Atoms
from ase.db import connect
import pytest

from kinbot import constants, l1_screening as screen
from kinbot.execution import atomic_json
from kinbot.parameters import Parameters
from kinbot.theory import resolve_profiles


def case(tmp_path):
    worker = tmp_path/'100'; worker.mkdir()
    atomic_json(tmp_path/'pes_exploration.json', {'initial_worker': '100'})
    inp = worker/'100.json'
    atomic_json(inp, dict(structure=['H',0,0,0,'H',0,0,0.7], mult=1,
                         barrier_threshold=200, l1_pes_energy_ceiling=150,
                         l1_profile=dict(calculator='gaussian',method='HF',basis='sto-3g')))
    par = Parameters(str(inp)).par
    profiles = {k: asdict(v) for k,v in resolve_profiles(par).items()}
    routes = dict(profile_sha256='test', profiles=profiles,
                  jobs={n:dict(level='l1') for n in ['100_well', '100_test', '200_well', 'H_well']})
    atomic_json(worker/'.kinbot_theory_jobs.json', routes)
    db = connect(worker/'kinbot.db')
    for name, symbols, energy in [('100_well','H2',-1), ('100_test','H2',-1+160/constants.AUtoKCAL),
                                  ('200_well','H2',-1+150/constants.AUtoKCAL), ('H_well','H',-.2)]:
        db.write(Atoms(symbols), name=name,
                 data=dict(status='normal', energy=energy/constants.EVtoHARTREE, zpe=0.))
    return worker, par


def test_absolute_reference_and_threshold_boundary(tmp_path):
    worker, par = case(tmp_path)
    before = (worker/'kinbot.db').read_bytes()
    assert screen.screen(par,['100_test'],reaction='100_test',kind='barrier',directory=worker)
    assert not screen.screen(par,['200_well'],reaction='test',kind='products',directory=worker)
    result=json.loads((worker/'l1_screening.json').read_text())['decisions']['100_test:barrier']
    assert result['relative_kcal_mol']==pytest.approx(160)
    assert (worker/'kinbot.db').read_bytes()==before
    assert not list(worker.glob('*.lock'))


def test_profile_stoichiometry_status_zpe_and_special_theory_are_not_guessed(tmp_path):
    worker, par = case(tmp_path)
    for names, comparable in [(['H_well'],True), (['missing'],True), (['100_test'],False)]:
        assert not screen.screen(par,names,reaction='test',kind='barrier',directory=worker,comparable=comparable)
    route=worker/'.kinbot_theory_jobs.json'
    r=json.loads(route.read_text()); r['jobs']['100_test']['level']='l2'; atomic_json(route,r)
    assert not screen.screen(par,['100_test'],reaction='test',kind='barrier',directory=worker)
    db=connect(worker/'kinbot.db')
    db.write(Atoms('H2'),name='200_well',data=dict(status='normal',energy=-.2))
    assert not screen.screen(par,['200_well'],reaction='test',kind='products',directory=worker)


def test_sum_fragments_only_for_same_stoichiometry(tmp_path):
    worker, par = case(tmp_path)
    assert screen.screen(par,['H_well','H_well'],reaction='test',kind='products',directory=worker)
    other=deepcopy(par);other['charge']=1
    assert not screen.screen(other,['H_well','H_well'],reaction='test',kind='products',directory=worker)


def test_disabled_makes_no_files_or_database_reads(tmp_path):
    assert not screen.screen({},['no-db'],reaction='x',kind='barrier',directory=tmp_path)
    assert not list(tmp_path.iterdir())


def test_report_is_read_only_and_does_not_claim_vdw_product_filter(tmp_path):
    worker,par=case(tmp_path)
    atomic_json(worker/'accepted_l2_network.json',dict(network=dict(reactions=[
        dict(name='100_test',family='test',products=['200'],product_complex='complex') ])))
    before={str(p):p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    report=screen.report(tmp_path,150)
    assert report['counts']['accepted_reaction_decisions']['exclude']==1
    assert report['workers'][0]['channels'][0]['products']['decision']=='unknown'
    assert before=={str(p):p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}


def test_accepted_replay_is_not_removed(tmp_path, monkeypatch):
    from kinbot.reaction_generator import ReactionGenerator
    generator=SimpleNamespace(search_checkpoint=SimpleNamespace(replaying={'accepted'}),par={})
    def forbidden(*a,**k):
        pytest.fail('Accepted replay must not be re-screened')
    monkeypatch.setattr(screen,'screen',forbidden)
    assert not ReactionGenerator._screen_l1(generator,SimpleNamespace(instance_name='accepted'),[], 'products')


def test_driver_stops_high_l1_before_irc_and_preserves_rejection_on_restart(tmp_path, monkeypatch):
    from kinbot.reaction_generator import ReactionGenerator
    worker, par = case(tmp_path)
    monkeypatch.chdir(worker)
    monkeypatch.setattr('kinbot.reaction_generator.time.sleep',lambda *a:None)
    par['delete_intermediate_files']=1
    native=worker/'100_test.log';native.write_text('retained native evidence')
    def make(status):
        reaction=SimpleNamespace(instance_name='100_test',family_name='test',instance=[0,1],
            max_step=12,scan=False,skip=False,products=[])
        species=SimpleNamespace(reac_obj=[reaction],reac_inst=[[0,1]],reac_type=['test'],
            reac_ts_done=[status],reac_step=[13],reac_scan_energy=[[]])
        qc=SimpleNamespace(qc='gauss',job_ids={},check_qc=lambda *a:'normal',
            get_qc_energy=lambda *a:(0,-1+160/constants.AUtoKCAL),get_qc_zpe=lambda *a:(0,0.))
        return ReactionGenerator(species,par,qc,str(worker/'100.json'))
    generator=make(1)
    generator.generate()
    assert generator.species.reac_ts_done==[-999]
    restored=make(0)
    restored.qc.check_qc=lambda *a:pytest.fail('Restart reopened a screened reaction')
    restored.generate()
    assert restored.species.reac_ts_done==[-999]
    assert native.read_text()=='retained native evidence'


@pytest.mark.parametrize('value',[True,0,-1,float('nan'),float('inf'),'150'])
def test_invalid_limits(value):
    with pytest.raises(ValueError):screen.validate(value)
