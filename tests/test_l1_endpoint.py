from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect

from kinbot import constants, l1_endpoint, l1_screening
from kinbot.execution import atomic_json
from test_l1_screening import case


@pytest.fixture
def endpoint(tmp_path, monkeypatch):
    worker, par = case(tmp_path); monkeypatch.chdir(worker)
    par['l1_vdw_depth_check'] = True
    name = '100_test_IRC_F_prod'
    atoms = Atoms('H2', positions=[[0,0,0], [0,0,.7]])
    db = connect(worker/'kinbot.db')
    source_id = db.write(atoms, name=name, data=dict(status='normal', energy=-1/constants.EVtoHARTREE))
    job = name+'_freq_recovery_'+str(source_id)
    routes = l1_endpoint.json.loads((worker/'.kinbot_theory_jobs.json').read_text())
    routes['jobs'][name[:-5]] = dict(level='l1')
    routes['jobs'][job] = dict(level='l1')
    atomic_json(worker/'.kinbot_theory_jobs.json', routes)
    species = SimpleNamespace(name=name, atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    qc = SimpleNamespace(qc_freq=lambda *a, **k: job, check_qc=lambda j: 'running')
    return worker, par, db, species, qc, job, atoms


def test_wait_then_same_level_completed_frequency_and_restart_reuse(endpoint):
    worker, par, db, species, qc, job, atoms = endpoint
    before = (worker/'kinbot.db').read_bytes()
    result = l1_endpoint.ensure(qc, species, par, family='test')
    assert result['status'] == 'waiting' and result['source']['zpe_hartree'] is None
    assert (worker/'kinbot.db').read_bytes() == before
    db.write(atoms, name=job, data=dict(status='normal', energy=-1/constants.EVtoHARTREE,
                                       zpe=.001, frequencies=[500.]))
    qc.check_qc = lambda j: 'normal'
    result = l1_endpoint.ensure(qc, species, par, family='test')
    assert result['status'] == 'complete' and result['point']['zpe_hartree'] == .001
    journal = (worker/'l1_endpoint_checks.json').read_bytes()
    assert l1_endpoint.ensure(qc, species, par, family='test') == result
    assert (worker/'l1_endpoint_checks.json').read_bytes() == journal
    assert db.get(name=species.name).data.get('zpe') is None


@pytest.mark.parametrize('bad', ['energy', 'geometry', 'zpe', 'frequency', 'level', 'failed'])
def test_unusable_frequency_cannot_become_a_depth_or_trigger_l2(endpoint, bad):
    worker, par, db, species, qc, job, atoms = endpoint
    data = dict(status='normal', energy=-1/constants.EVtoHARTREE, zpe=.001, frequencies=[500.])
    if bad == 'energy': data['energy'] = float('nan')
    if bad == 'geometry': atoms.positions[1,2] += .01
    if bad == 'zpe': data.pop('zpe')
    if bad == 'frequency': data['frequencies'] = [-200.]
    if bad == 'level':
        routes = l1_endpoint.json.loads((worker/'.kinbot_theory_jobs.json').read_text())
        routes['jobs'][job]['level'] = 'l2'; atomic_json(worker/'.kinbot_theory_jobs.json', routes)
    db.write(atoms, name=job, data=data)
    qc.check_qc = lambda j: 'error' if bad == 'failed' else 'normal'
    result = l1_endpoint.ensure(qc, species, par, family='test')
    assert result['status'] == 'unresolved' and 'point' not in result
    from kinbot.reaction_generator import ReactionGenerator
    generator = SimpleNamespace(qc=qc, par=par, search_checkpoint=SimpleNamespace(replaying=set()))
    obj = SimpleNamespace(instance_name='100_test', family_name='test', irc_prod=species)
    assert ReactionGenerator._endpoint_l1_check(generator, obj)['status'] == 'unresolved'


def test_existing_valid_zpe_needs_no_new_job(endpoint):
    worker, par, db, species, qc, _, atoms = endpoint
    db.write(atoms, name=species.name, data=dict(status='normal', energy=-1/constants.EVtoHARTREE,
                                               zpe=.001, frequencies=[500.]))
    qc.qc_freq = lambda *a, **k: pytest.fail('Valid existing endpoint must be reused')
    assert l1_endpoint.ensure(qc, species, par, family='test')['status'] == 'complete'


def test_native_frequency_offsets_reuse_coherent_records_and_retain_old_holds(tmp_path):
    from tests.endpoint_fixture import stage_native_endpoints
    fixture = stage_native_endpoints(tmp_path)
    original = (tmp_path/'kinbot.db').read_bytes()
    results = []
    for species, example in zip(fixture.species, fixture.examples):
        result = l1_endpoint.ensure(fixture.qc, species, fixture.par,
                                   family='test', directory=tmp_path)
        assert result['status'] == 'complete'
        point = result['point']
        assert point['energy_hartree'] == pytest.approx(example['frequency']['energy_hartree'], abs=1e-10)
        assert point['zpe_hartree'] == example['frequency']['zpe_hartree']
        assert result['energy_zpe_source'] == point['name']
        assert result['energy_comparison']['diagnostic_only']
        assert 2e-6 < result['energy_comparison']['difference_hartree'] < 3e-6
        results.append(result)
    assert (tmp_path/'kinbot.db').read_bytes() == original
    journal = (tmp_path/'l1_endpoint_checks.json').read_bytes()
    saved = l1_endpoint.json.loads(journal)
    assert len(saved['history']) == 3
    assert all(r['status'] == 'unresolved' for r in saved['history'])
    for species, result in zip(fixture.species, results):
        assert l1_endpoint.ensure(fixture.qc, species, fixture.par,
                                  family='test', directory=tmp_path) == result
    assert (tmp_path/'l1_endpoint_checks.json').read_bytes() == journal


def test_disabled_check_preserves_existing_behavior(endpoint):
    from kinbot.reaction_generator import ReactionGenerator
    _, par, _, species, qc, _, _ = endpoint
    obj = SimpleNamespace(instance_name='accepted', irc_prod=species)
    qc.qc_freq = lambda *a, **k: pytest.fail('No calculation for disabled check')
    generator = SimpleNamespace(qc=qc, par=par, search_checkpoint=SimpleNamespace(replaying={'accepted'}))
    generator.par = dict(par, l1_vdw_depth_check=False)
    assert ReactionGenerator._endpoint_l1_check(generator, obj) is None


NATIVE_REPLAYS=l1_endpoint.json.loads((Path(__file__).parent/'reference/l1_accepted_replay.json').read_text())['records']


@pytest.mark.parametrize('example',NATIVE_REPLAYS,ids=lambda e:e['reaction'])
def test_accepted_native_replay_reads_completed_zpe_instead_of_inventing_complex(tmp_path,monkeypatch,example):
    from tests.endpoint_fixture import stage_native_replay
    from kinbot.reaction_generator import ReactionGenerator
    fixture=stage_native_replay(tmp_path,example);monkeypatch.chdir(tmp_path)
    owner=SimpleNamespace(qc=fixture.qc,par=fixture.par,
        search_checkpoint=SimpleNamespace(replaying={example['reaction']}))
    before=(tmp_path/'kinbot.db').read_bytes()
    for _ in range(2):
        checked=ReactionGenerator._endpoint_l1_check(owner,fixture.reaction)
        assert checked is not None and checked['status']=='complete'
        point=checked['point']
        separated=sum(p['energy_hartree']+p['zpe_hartree'] for p in checked['fragment_points'])
        depth=(separated-point['energy_hartree']-point['zpe_hartree'])*constants.AUtoKCAL
        assert depth==pytest.approx(example['matched_depth_kcal_mol'],abs=1e-7)
        assert depth<.5<example['missing_zpe_depth_kcal_mol']
        assert fixture.qc.check_qc(example['source']['name']+'_high')==0
    assert (tmp_path/'kinbot.db').read_bytes()==before
    assert not list(tmp_path.glob('*.py')) and not list(tmp_path.glob('*.sbatch'))


@pytest.mark.parametrize('depth', [.4, .8])
@pytest.mark.parametrize('energy_offset', [0., 2.529e-6])
def test_reaction_owner_waits_reconstructs_and_uses_matched_depth(tmp_path, monkeypatch, depth, energy_offset):
    """Run the real search loop/checkpoint with simulated native observations.

    Two identical fragments exercise shared calculation adoption. The owner
    dies while the endpoint is pending, then again while the fragment is
    pending; reconstruction must reach the same jobs without reaching L2.
    """
    from kinbot import reaction_generator as module
    from kinbot.species_routing import routing_name
    from kinbot.stationary_pt import StationaryPoint
    worker, par = case(tmp_path); monkeypatch.chdir(worker)
    par.update(l1_vdw_depth_check=True, vdW_detection=.5)
    atoms = Atoms('H4', positions=[[0,0,0], [0,0,.7], [3,0,0], [3,0,.7]])
    endpoint_name = '100_test_IRC_F_prod'
    def point(name):
        obj = StationaryPoint(name, 0, 1, atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        obj.characterize()
        return obj
    endpoint = point(endpoint_name)
    fragments, _ = endpoint.start_multi_molecular(vary_charge=True)
    assert len(fragments) == 2
    fragment_name = routing_name(fragments[0])+'_well'
    db = connect(worker/'kinbot.db')
    def record(name, geometry, energy, **extra):
        return db.write(geometry, name=name, data=dict(status='normal',
            energy=energy/constants.EVtoHARTREE, **extra))
    record('100_well', atoms, -1., zpe=0.)
    record(routing_name(endpoint)+'_well', atoms, -1., zpe=0.)
    record('100_test', atoms, -.99, zpe=.001, frequencies=[-500.])
    endpoint_id = record(endpoint_name, atoms, -1.01)
    fragment_atoms = Atoms(fragments[0].atom, positions=fragments[0].geom)
    fragment_id = record(fragment_name, fragment_atoms, -.5)
    jobs = {endpoint_name: endpoint_name+'_freq_recovery_'+str(endpoint_id),
            fragment_name: fragment_name+'_freq_recovery_'+str(fragment_id)}
    routes = l1_endpoint.json.loads((worker/'.kinbot_theory_jobs.json').read_text())
    for name in [endpoint_name[:-5], fragment_name, *jobs.values()]:
        routes['jobs'][name] = dict(level='l1')
    atomic_json(worker/'.kinbot_theory_jobs.json', routes)
    calls = []
    def freq(species, source_job, **kwargs):
        assert kwargs == dict(high_level=0, endpoint=True)
        calls.append(source_job)
        return jobs[source_job]
    def row(name):
        return list(db.select(name=name))[-1]
    def geom(name, natom, reorder=False):
        observed = row(name)
        result = (0, observed.positions.copy())
        return (*result, np.array(observed.toatoms().get_chemical_symbols())) if reorder else result
    qc = SimpleNamespace(qc='gauss', db=db, job_ids={}, qc_freq=freq,
        qc_opt=lambda *a, **k: None,
        check_qc=lambda name: 'normal' if db.count(name=name) else 'running',
        get_qc_geom=geom,
        get_qc_energy=lambda name: (0,row(name).data['energy']*constants.EVtoHARTREE),
        get_qc_zpe=lambda name: (0,row(name).data.get('zpe',0.)),
        get_qc_freq=lambda name, natom: (0,row(name).data.get('frequencies',[])))
    # Completed IRC work is reconstructed from saved evidence after each
    # owner restart. No new IRC/optimization is permitted in this test.
    monkeypatch.setattr(module, 'IRC', lambda *a: SimpleNamespace(
        check_irc=lambda: ('normal','normal'), irc2stationary_pt=lambda: point(endpoint_name)))
    monkeypatch.setattr(module, 'Optimize', lambda *a, **k: pytest.fail('Premature L2 work'))
    class StopPoll(Exception): pass
    def make_owner():
        obj = SimpleNamespace(instance_name='100_test', family_name='test', instance=[0,1],
            max_step=12, scan=False, skip=False, products=[], prod_opt=[],
            irc_prod=point(endpoint_name), do_vdW=False)
        species = point('parent')
        species.reac_obj=[obj]; species.reac_inst=[[0,1]]; species.reac_type=['test']
        species.reac_ts_done=[2]; species.reac_step=[13]; species.reac_scan_energy=[[]]
        return module.ReactionGenerator(species,par,qc,str(worker/'100.json'))
    def poll(owner, waiting):
        initial = len(calls)
        def sleep(_):
            if owner.species.reac_ts_done == [3] or len(calls) >= initial+2:
                raise StopPoll
        monkeypatch.setattr(module.time, 'sleep', sleep)
        with pytest.raises(StopPoll): owner.generate()
        assert owner.species.reac_ts_done == [2 if waiting else 3]
    poll(make_owner(), True)
    assert set(calls) == {endpoint_name}
    record(jobs[endpoint_name], atoms, -1.01+energy_offset,
           zpe=.012-depth/constants.AUtoKCAL-energy_offset,
           frequencies=[100.]*6)
    calls.clear()
    poll(make_owner(), True)
    assert set(calls) == {endpoint_name,fragment_name}
    record(jobs[fragment_name], fragment_atoms, -.5, zpe=.001, frequencies=[500.])
    owner = make_owner()
    poll(owner, False)
    reaction = owner.species.reac_obj[0]
    assert reaction.vdW_depth == pytest.approx(depth)
    assert reaction.do_vdW == (depth > .5)
    assert db.get(id=endpoint_id).data.get('zpe') is None
    assert db.get(id=fragment_id).data.get('zpe') is None
    decision = l1_endpoint.json.loads((worker/'l1_screening.json').read_text())['decisions']['100_test:products']
    expected = [jobs[endpoint_name]] if depth > .5 else [jobs[fragment_name]]*2
    assert [p['name'] for p in decision['points']] == expected
    assert all(fragment.freq == [500.] for fragment in reaction.irc_fragments)

    # Complete cached L2 observations represent an already accepted search.
    # The shallow endpoint deliberately has NO *_prod_high job. Restoring
    # it must not turn missing raw-IRC ZPE into an invented deep complex.
    from kinbot.qc import QuantumChemistry
    from types import MethodType
    record('100_test_high', atoms, -.99, zpe=.001, frequencies=[-500.])
    record(fragment_name+'_high', fragment_atoms, -.5, zpe=.001, frequencies=[500.])
    if depth > .5:
        record(endpoint_name+'_high', atoms, -1.01+energy_offset,
               zpe=.012-depth/constants.AUtoKCAL-energy_offset,
               frequencies=[100.]*6)
    qc.check_qc = lambda name: 'normal' if db.count(name=name) else 0
    qc.use_sella=True; qc.frequency_mode='native_hessian'
    # Use the real fixed-frequency cache adoption and QC submission guard.
    qc.qc_freq = MethodType(QuantumChemistry.qc_freq, qc)
    qc.qc_opt = lambda sp, geom, **kw: QuantumChemistry._submit_qc(
        qc, routing_name(sp)+'_well', 1)
    dispatched = []
    class CachedOptimize:
        def __init__(self, sp, par, backend, just_high=False):
            self.species=sp; self.just_high=just_high
            self.name=(sp.name if just_high or sp.wellorts else routing_name(sp)+'_well')+'_high'
            self.shigh=1; self.shir=1
        def do_optimization(self):
            dispatched.append(self.name)
            QuantumChemistry._submit_qc(qc, self.name, 1)
            observed = row(self.name)
            self.species.energy = observed.data['energy']*constants.EVtoHARTREE
            self.species.zpe = observed.data['zpe']
            self.species.freq = list(observed.data['frequencies'])
            self.species.reduced_freqs = list(self.species.freq)
    monkeypatch.setattr(module, 'Optimize', CachedOptimize)
    owner.species.reac_ts_done[0] = -1
    owner.search_checkpoint.save(0)
    db_before=(worker/'kinbot.db').read_bytes()
    screening_before=(worker/'l1_screening.json').read_bytes()
    for _ in range(2):
        replay = make_owner()
        def stop_after_dispatch(_):
            if replay.species.reac_ts_done == [4]: raise StopPoll
        monkeypatch.setattr(module.time, 'sleep', stop_after_dispatch)
        with pytest.raises(StopPoll): replay.generate()
        assert replay.species.reac_obj[0].do_vdW == (depth > .5)
        assert replay.species.reac_obj[0].vdW_depth == pytest.approx(depth)
        assert l1_endpoint.json.loads(replay.search_checkpoint.path.read_text())['reactions']['100_test']['status']==-1
    assert (endpoint_name+'_high' in dispatched) == (depth > .5)
    # Also finish final product-complex reassessment and re-publish the
    # accepted reaction, rather than stopping at the dispatch boundary only.
    published=[]
    monkeypatch.setattr(module.postprocess, 'createPESViewerInput',
                        lambda *a: published.append(True))
    polls=[]
    def bounded_poll(_):
        polls.append(True)
        assert len(polls)<30, 'Accepted replay failed to complete'
    monkeypatch.setattr(module.time, 'sleep', bounded_poll)
    final=make_owner(); final.generate()
    assert final.species.reac_ts_done==[-1] and published
    assert final.species.reac_obj[0].do_vdW == (depth>.5)
    assert (worker/'kinbot.db').read_bytes()==db_before
    assert (worker/'l1_screening.json').read_bytes()==screening_before


@pytest.mark.parametrize('code', ['gauss', 'qchem'])
@pytest.mark.parametrize('outcome', ['normal', 'missing_zpe', 'imaginary'])
def test_fixed_frequency_worker_parses_native_evidence_before_publishing(tmp_path, monkeypatch, code, outcome):
    from ase.calculators.calculator import Calculator, all_changes
    from kinbot import kb_path
    label = 'endpoint'
    marker = tmp_path/(label+('.log' if code == 'gauss' else '.out'))
    output = tmp_path/(label+('.log' if code == 'gauss' else '_freq.out'))
    frequencies = [-300., 300., 400.] if outcome == 'imaginary' else [200., 300., 400.]
    text = (' Frequencies -- '+ ' '.join(map(str,frequencies)) if code == 'gauss'
            else ' Frequency: '+ ' '.join(map(str,frequencies)))+'\n'
    if outcome != 'missing_zpe':
        text += (' Zero-point correction= 0.001 (Hartree/Particle)' if code == 'gauss'
                 else ' Zero point vibrational energy: 1.0 kcal/mol')+'\n'
    class NativeFixture(Calculator):
        implemented_properties = ['energy']
        def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
            super().calculate(atoms,properties,system_changes)
            output.write_text(text)
            self.results = dict(energy=-10.)
    monkeypatch.chdir(tmp_path)
    module = ('kinbot.ase_modules.calculators.gaussian.Gaussian' if code == 'gauss'
              else 'kinbot.ase_modules.calculators.qchem.QChem')
    monkeypatch.setattr(module, NativeFixture)
    monkeypatch.setattr('kinbot.utils.iowait', lambda *a: None)
    template = (Path(kb_path)/'tpl/ase_endpoint_freq.tpl.py').read_text()
    script = template.format(label=label, atom=['O','H','H'],
        geom=[[0,0,0],[0,0,1],[0,1,0]], kwargs={}, code=code,
        qc_command='unused-native-fixture', ppn=1, imagfreq_threshold=50., working_dir=tmp_path)
    exec(compile(script, '<endpoint fixture>', 'exec'), {})
    row = connect(tmp_path/'kinbot.db').get(name=label)
    assert row.data['status'] == ('normal' if outcome == 'normal' else 'error')
    assert marker.read_text().endswith('done\n')
    if outcome == 'normal':
        assert row.data['zpe'] > 0
        assert list(row.data['frequencies']) == frequencies
