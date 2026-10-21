"""An optional absolute L2 ceiling, referenced to the initial PES well."""
import json
import math
from pathlib import Path

from kinbot import constants
from kinbot.execution import atomic_json


def validate(reference):
    if reference is None:
        return
    if (not isinstance(reference,dict) or set(reference)!= {'zero_k_hartree','ceiling_kcal_mol','source'}
            or not isinstance(reference['source'],str)):
        raise ValueError('pes_barrier_reference needs zero_k_hartree, ceiling_kcal_mol and source.')
    for key in ('zero_k_hartree','ceiling_kcal_mol'):
        value=reference[key]
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
            raise ValueError(f'pes_barrier_reference.{key} must be finite.')
    if reference['ceiling_kcal_mol']<=0:
        raise ValueError('PES barrier ceiling must be positive.')


def apply(species, parameters):
    ref=parameters.get('pes_barrier_reference')
    if ref is None:
        return
    validate(ref)
    if not parameters.get('high_level'):
        raise ValueError('An L2 PES barrier reference requires high_level.')
    energy=float(species.energy)+float(species.zpe)
    if not math.isfinite(energy):
        raise ValueError('PES barrier requires a finite accepted well energy.')
    relative=(energy-ref['zero_k_hartree'])*constants.AUtoKCAL
    parameters['barrier_threshold_L2']=ref['ceiling_kcal_mol']-relative
    parameters['barrier_threshold']=parameters['barrier_threshold_L2']+parameters['barrier_threshold_add']
    atomic_json(Path('barrier_policy.json'),dict(reference=ref,well_relative_kcal_mol=relative,
        barrier_threshold_L2=parameters['barrier_threshold_L2'],barrier_threshold=parameters['barrier_threshold']))


def for_case(case, ceiling):
    phase=json.loads((case/'pes_exploration.json').read_text())
    initial=phase.get('initial_worker')
    if not initial:
        initial=(case/'chemids').read_text().splitlines()[0]
    inp=case/initial/(initial+'.json')
    parameters=json.loads(inp.read_text())
    p=inp.parent/parameters.get('anl_automatic_directory','anl_automatic')/'network_request.json'
    plan=json.loads(p.read_text())
    records=[r for r in plan['network']['stable_species'] if r['routing_name']==initial]
    if len(records)!=1:
        raise ValueError('Initial PES well L2 energy is ambiguous.')
    record=records[0]
    ref=dict(zero_k_hartree=record['energy_hartree']+record['harmonic_zpe_hartree'],
             ceiling_kcal_mol=ceiling,source=str(p))
    validate(ref)
    # This migration changes search limits, not scientific acceptance of an
    # existing handoff. Reject an implicit removal before stopping controllers.
    for inp in case.glob('[0-9]*/*.json'):
        if inp.name!=inp.parent.name+'.json':
            continue
        par=json.loads(inp.read_text())
        q=inp.parent/par.get('anl_automatic_directory','anl_automatic')/'network_request.json'
        if not q.exists():
            continue
        network=json.loads(q.read_text())['network']
        for name,r in network['transition_states'].items():
            value=(r['energy_hartree']+r['harmonic_zpe_hartree']-ref['zero_k_hartree'])*constants.AUtoKCAL
            if value>ceiling:
                raise ValueError(f'{name}: saved accepted barrier {value:.2f} exceeds {ceiling}; explicit network pruning is required.')
    return ref
