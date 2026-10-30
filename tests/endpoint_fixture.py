"""Exported endpoint observations for restart/acceptance regression tests."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace, MethodType

from ase import Atoms
from ase.db import connect

from kinbot import constants, l1_screening
from kinbot.execution import atomic_json


def stage_native_endpoints(worker):
    fixture = json.loads((Path(__file__).parent/'reference/l1_endpoint_frequency_offsets.json').read_text())
    profile = fixture['records'][0]['old_check']['theory']
    par = dict(qc='gauss', qc_command='g16', method=profile['method'], basis=profile['basis'],
               high_level_method='B2PLYP', high_level_basis='cc-pVTZ',
               l1_profile=profile, imagfreq_threshold=50., l1_vdw_depth_check=True)
    assert l1_screening.theory(par) == profile
    db = connect(worker/'kinbot.db')
    checks, jobs, species, routes = {}, {}, [], {}
    for example in fixture['records']:
        source, frequency = example['source'], example['frequency']
        for point in (source, frequency):
            db.write(Atoms(point['symbols'], positions=point['positions']), name=point['name'],
                     data=dict(status='normal', energy=point['energy_hartree']/constants.EVtoHARTREE,
                               zpe=point['zpe_hartree'], frequencies=point['frequencies']))
            routes[point['name']] = dict(level='l1')
        checks[source['name']] = deepcopy(example['old_check'])
        jobs[source['name']] = frequency['name']
        species.append(SimpleNamespace(name=source['name'], atom=source['symbols'],
                                       geom=source['positions'], mult=2, charge=0))
    atomic_json(worker/'.kinbot_theory_jobs.json', dict(schema=1, profiles={'l1': profile},
                profile_sha256='native-fixture', jobs=routes))
    atomic_json(worker/'l1_endpoint_checks.json', dict(schema=1, checks=checks))
    calls = []
    def completed_frequency(species, source_job, **kwargs):
        assert kwargs == dict(high_level=0, endpoint=True)
        calls.append(source_job)
        return jobs[source_job]  # Completed native job, no execution.
    qc = SimpleNamespace(db=db, qc_freq=completed_frequency, check_qc=lambda job: 'normal')
    return SimpleNamespace(par=par, qc=qc, species=species, examples=fixture['records'], calls=calls)


def stage_native_replay(worker, example):
    """Rebase exported row IDs into a disposable DB, preserving observations."""
    from kinbot.qc import QuantumChemistry
    profile=example['theory']
    par=dict(qc='gauss', qc_command='g16', method=profile['method'], basis=profile['basis'],
             high_level_method='B2PLYP', high_level_basis='cc-pVTZ', l1_profile=profile,
             scan_method=profile['method'], scan_basis=profile['basis'],
             imagfreq_threshold=50., l1_vdw_depth_check=True, vdW_detection=.5)
    db=connect(worker/'kinbot.db');routes={}
    def write(point):
        routes[point['name']]=dict(level='l1')
        return db.write(Atoms(point['symbols'],positions=point['positions']),name=point['name'],
                        data=dict(status='normal',energy=point['energy_hartree']/constants.EVtoHARTREE,
                                  zpe=point['zpe_hartree'],frequencies=point['frequencies']))
    source=example['source'];frequency=deepcopy(example['frequency'])
    source_id=write(source)
    frequency['name']=source['name']+'_freq_recovery_'+str(source_id)
    write(frequency)
    products=[]
    for point in example['fragments']:
        write(point)
        products.append(SimpleNamespace(name=point['name'].removesuffix('_well'),
            chemid=point['name'].removesuffix('_well'),atom=point['symbols'],geom=point['positions'],
            optical_reference={'status':'assigned','canonical_graphs':[]}))
    atomic_json(worker/'.kinbot_theory_jobs.json',dict(schema=1,profiles={'l1':profile},
                profile_sha256='native-replay-fixture',jobs=routes))
    def status(job):
        rows=list(db.select(name=job))
        return rows[-1].data['status'] if rows else 0
    qc=SimpleNamespace(db=db,qc='gauss',use_sella=True,frequency_mode='native_hessian',
                       check_qc=status,_search_replay=example['reaction'])
    qc.qc_freq=MethodType(QuantumChemistry.qc_freq,qc)
    qc.submit_qc=MethodType(QuantumChemistry.submit_qc,qc)
    obj=SimpleNamespace(instance_name=example['reaction'],
        family_name=example['search_entry']['signature']['family'],products=products,
        irc_prod=SimpleNamespace(name=source['name'],atom=source['symbols'],geom=source['positions']))
    return SimpleNamespace(qc=qc,par=par,reaction=obj,frequency=frequency)
