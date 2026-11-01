"""Disconnected MESS populations must not be identified by ChemID sums alone."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from ase import Atoms
from ase.db import connect
import numpy as np
import pytest

from kinbot import constants, pes, symmetry
from kinbot.execution import file_sha256
from kinbot.mess_mirrors import population_keys, complete_mirror_channels, _POP, _EDGE
from kinbot.mess_networks import split_model, _HEADER, connected_models
from kinbot.parameters import Parameters
from kinbot.reaction_path import prepare_stereopath
from kinbot.run_format import ensure_current_run
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import require_supported_identity
from kinbot.stereo_routing import guard_well_job
from tests.test_mess_pes_artifacts import at, write_worker

FIXTURE = Path(__file__).parent/'reference/mess_complex_identity_collision.json'


def test_disconnected_population_identity_is_canonical_not_a_geometry_or_atom_order():
    row=json.loads(FIXTURE.read_text())['records'][0]['points'][2]
    points=[]
    for order in (np.arange(6), np.array([5,2,0,4,1,3])):
        p=StationaryPoint('encounter',row['charge'],row['multiplicity'],
            atom=np.array(row['symbols'])[order],
            geom=np.array(row['low']['positions'])[order][:,[2,0,1]]+[7.,-3.,1.])
        p.characterize()
        points.append(p)
    assert population_keys([points[0]])==population_keys([points[1]])
    key=population_keys([points[0]])
    for charge,spin in ((1,2),(0,4)):
        other=deepcopy(points[0]);other.charge=charge;other.mult=spin
        assert population_keys([other])!=key


def test_disconnected_mirror_uses_mirrored_identity_and_preserves_declared_namespace():
    # Structural metadata fixture, not an optimized molecular calculation.
    bond=np.zeros((6,6),int);bond[0,1:5]=1;bond[1:5,0]=1
    p=SimpleNamespace(name='encounter',chemid=123,charge=0,mult=2,
        atom=['C','H','F','Cl','Br','H'],bond=bond,bonds=[bond],
        geom=np.array([[0.,0.,0.],[1.,1.,1.],[1.,-1.,-1.],
                       [-1.,1.,-1.],[-1.,-1.,1.],[7.,0.,0.]]))
    identity=require_supported_identity(p)
    own,mirror=population_keys([p])
    assert identity['id']!=identity['mirror_id'] and own!=mirror
    reflected=deepcopy(p);reflected.geom *= [-1.,1.,1.]
    assert population_keys([reflected])==(mirror,own)
    assert population_keys([p],'racemic')==(own,own)
    p.optical_reference=identity
    p.geom=reflected.geom.copy()
    assert population_keys([p])==(own,mirror)
    # Existing connected species keep their public names unchanged.
    connected=deepcopy(reflected)
    connected.atom=connected.atom[:-1];connected.geom=connected.geom[:-1]
    connected.bond=connected.bond[:-1,:-1];connected.bonds=[connected.bond]
    assert population_keys([connected])[0]==routing_name(connected)


def native_points(record):
    points = []
    for row in record['points']:
        p = StationaryPoint(row['name'], row['charge'], row['multiplicity'],
            atom=row['symbols'], geom=np.array(row['low']['positions']))
        p.characterize()
        if row['role'] != 'ts':
            p.optical_reference = require_supported_identity(p)
        points.append(p)
    parent, ts, complex_, *products = points
    parent.name = record['worker']
    for product, name in zip(products, record['products']):
        product.name = name
    ts.wellorts = 1
    ts.optical_reference = parent.optical_reference
    prepare_stereopath(ts, parent, complex_)
    for p, row in zip(points, record['points']):
        high = row['high']
        p.geom = np.array(high['positions'])
        p.energy = high['energy_eV'] * constants.EVtoHARTREE
        p.zpe = high['zpe_hartree']
        p.freq = p.reduced_freqs = high['frequencies'].copy()
        symmetry.calculate_symmetry(p)
    return points


def stage_collision_network(root, *, mc=False, population='specified', reverse=False):
    root.mkdir(parents=True, exist_ok=True)
    records = json.loads(FIXTURE.read_text())['records']
    if reverse:
        records.reverse()
    settings = dict(barrier_threshold=150., high_level=1, rotor_scan=0,
        multi_conf_tst=int(mc), conformer_search=1, me=2, uq=0,
        epsilon=100., sigma=3., queuing='local', smiles='C', mult=2, pes=1,
        rotdpy_mess_mode='off', optical_population=population)
    source = root/'collision_input.json'
    source.write_text(json.dumps(settings))
    par = Parameters(str(source), show_warnings=False).par
    cases = []
    for record in records:
        parent, ts, complex_, *products = native_points(record)
        reaction = SimpleNamespace(instance_name=record['reaction'], ts=ts,
            products=products, prod_opt=[SimpleNamespace(species=p) for p in products],
            irc_prod=complex_, irc_prod_opt=SimpleNamespace(species=complex_), do_vdW=True, mp2=0)
        parent.reac_obj=[reaction];parent.reac_inst=[None];parent.reac_name=[reaction.instance_name]
        parent.reac_ts_done=[-1];parent.reac_type=['test']
        worker=root/record['worker'];worker.mkdir();(worker/'me').mkdir()
        ensure_current_run(worker, create=True)
        db=connect(str(worker/'kinbot.db'))
        for p, row in zip([parent, ts, complex_, *products], record['points']):
            name=row['name']+'_high'
            if row['role'] in ('parent','fragment'):
                guard_well_job(SimpleNamespace(db=db), p, p.geom, row['name'])
                guard_well_job(SimpleNamespace(db=db), p, p.geom, name)
            db.write(Atoms(p.atom, positions=p.geom),name=name,data=dict(status='normal',
                energy=row['high']['energy_eV'],zpe=p.zpe,frequencies=p.freq))
        cases.append(SimpleNamespace(root=root,worker=worker,par=par,parent=parent,
            reactions=[reaction],db_hash=file_sha256(worker/'kinbot.db')))
    return cases


def assemble_collision_network(cases, *, skip=False):
    t=cases[0]
    with at(t.root), patch('kinbot.pes.copy_from_kinbot'), \
            patch('kinbot.pes.create_pesviewer_input'), patch('kinbot.pes.create_rotdpy_inputs'), \
            patch('kinbot.pes.t1_analysis'), patch('kinbot.pes.get_l3energy',return_value=(0,-1)):
        pes.postprocess(dict(t.par,me_skip_vdW=int(skip)),[c.worker.name for c in cases],
                        'all',[],t.parent.mass)
    jobs=json.loads((t.root/'me/mess_networks.json').read_text())
    assert len(jobs)==2 and all(j['status']=='ready' for j in jobs)
    return {j['stem']:(t.root/'me'/(j['stem']+'.inp')).read_text() for j in jobs}


def assert_separate_complexes(cases, models, population='specified'):
    complexes=[c.reactions[0].irc_prod for c in cases]
    keys=[population_keys([p],population)[0] for p in complexes]
    assert len(set(keys))==2
    blocks=[b for text in models.values() for b in split_model(text)[1]]
    wells={_POP.search(b)[1]:_HEADER.search(b)[2].split()[0] for b in blocks
           if _HEADER.search(b)[1]=='Well'}
    assert set(keys) <= wells.keys() and wells[keys[0]] != wells[keys[1]]
    assert sum(text.count('Core PhaseSpaceTheory') for text in models.values())==2
    for text in models.values():
        local=split_model(text)[1]
        present={_POP.search(b)[1] for b in local if _HEADER.search(b)[1]=='Well'}
        selected=[i for i,key in enumerate(keys) if key in present]
        assert len(selected)==1
        index=selected[0]; key=keys[index]
        products=population_keys(cases[index].reactions[0].products,population)[0]
        assert {_POP.search(b)[1] for b in local if _HEADER.search(b)[1]=='Bimolecular'}=={products}
        assert any(_HEADER.search(b)[1]=='Barrier' and _EDGE.search(b)
                   and _EDGE.search(b)[3]==key and _HEADER.search(b)[2].split()[2]==wells[key]
                   for b in local)
        assert len(connected_models(text))==1
        assert complete_mirror_channels(text).split()==text.split()


@pytest.mark.parametrize('mc',[False,True])
@pytest.mark.parametrize('population',['specified','racemic'])
@pytest.mark.parametrize('reverse',[False,True])
def test_native_collision_survives_real_worker_and_pes_assembly(tmp_path,mc,population,reverse):
    cases=stage_collision_network(tmp_path,mc=mc,population=population,reverse=reverse)
    complexes=[c.reactions[0].irc_prod for c in cases]
    assert routing_name(complexes[0]) == routing_name(complexes[1])
    assert complexes[0].optical_reference['id'] != complexes[1].optical_reference['id']
    original=[deepcopy([(p.energy,p.zpe,p.geom.tolist(),p.freq) for p in
              [c.parent,c.reactions[0].ts,c.reactions[0].irc_prod,*c.reactions[0].products]]) for c in cases]
    previous=None
    for _ in range(2):
        for case in cases:write_worker(case)
        models=assemble_collision_network(cases)
        assert all(c.parent.reac_ts_done==[-1] for c in cases)
        keys=[population_keys([p],population)[0] for p in complexes]
        assert_separate_complexes(cases,models,population)
        omitted=assemble_collision_network(cases,skip=True)
        assert all('Core PhaseSpaceTheory' not in text for text in omitted.values())
        assert all(key not in {_POP.search(b)[1] for text in omitted.values()
                              for b in split_model(text)[1] if _HEADER.search(b)[1]=='Well'} for key in keys)
        for c, before in zip(cases,original):
            assert file_sha256(c.worker/'kinbot.db')==c.db_hash
            assert before==[(p.energy,p.zpe,p.geom.tolist(),p.freq) for p in
                           [c.parent,c.reactions[0].ts,c.reactions[0].irc_prod,*c.reactions[0].products]]
        if previous is not None:assert (models,omitted)==previous
        previous=models,omitted
