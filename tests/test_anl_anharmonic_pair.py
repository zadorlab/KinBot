"""Native scalar CBS verification and independent per-basis mode projection."""
from copy import deepcopy
import json
import re

import numpy as np
import pytest
from ase.io import read

from kinbot import frequencies
from kinbot.anl import anharmonic, thermal, dispatch
from kinbot.anl.extrapolation import extrapolation_coefficient
from kinbot.anl.recipes import recipe
from kinbot.calculation import array_fingerprint
from kinbot.execution import atomic_json
from tests.test_anl_workflow import _completed_zpe_pair
from tests.test_anl_thermal import FIX, point
from tests.test_rotor_projection import rotor_fixture


@pytest.mark.parametrize('tier', ['ANL1', 'ANL1-F12'])
def test_native_pair_is_same_correction_expression_in_both_recipes(tmp_path,tier):
    graph,_,upper_hash=_completed_zpe_pair(tmp_path,kind='vpt2')
    _,spec,state=dispatch._load(graph)
    # This native fixture contains optimized TZ/QZ geometries and VPT2 results;
    # its QZ optimized geometry supplies the completed fixture gate.
    spec['intent']=dict(anharmonic_basis_extension=dict(lower_task='vpt2_tz',upper_task='vpt2_qz',
        geometry_task='geom_qz',harmonic_task='geom_qz',bases=['cc-pVTZ','cc-pVQZ'],extrapolation_power=3.7))
    from kinbot.anl.live_control import commit
    commit(graph,spec,state,set(),{})
    requirement=next(r for r in recipe(tier,vpt2_method='B2PLYP-D3BJ',vpt2_cbs=True).requirements
                     if r.key=='vpt2_correction')
    result=anharmonic.component(graph,state_id='state-A',requirement=requirement)
    from ase.units import invcm,Hartree
    assert result.value_hartree==pytest.approx((-136+extrapolation_coefficient(4,power=3.7))*invcm/Hartree)
    assert result.quantity=='correction' and result.geometry_sha256==upper_hash
    _,spec,state=dispatch._load(graph)
    spec['intent']['live_protocol']=dict(omitted_tasks=['vpt2_qz'])
    commit(graph,spec,state,set(),{})
    assert not anharmonic.available(spec,state)
    with pytest.raises(ValueError,match='not complete'):
        anharmonic.component(graph,state_id='state-A')


def test_each_basis_matches_its_own_vectors_at_its_own_geometry(tmp_path):
    p=point(); atoms=read(FIX/'minimum.l2.xyz'); p.geom=atoms.positions
    p.anl_thermal_model['l2_geometry_sha256']=array_fingerprint(p.geom)
    lower_native=FIX/'minimum.gaussian_modes.txt'
    f,_=thermal.gaussian_modes(lower_native.read_text(),list(p.atom),p.geom)
    lower_shift=np.array([-20.,-40.,-60.]); upper_shift=np.array([-22.,-43.,-64.])
    def parsed(order,basis,shifts):
        return dict(harmonic_fundamentals_cm_inverse=f[order],
            anharmonic_fundamentals_cm_inverse=(f+shifts)[order],
            method='B2PLYP',basis=basis,dispersion='GD3BJ')
    thermal.attach_vpt2(p,parsed([2,1,0],'cc-pVTZ',lower_shift),lower_native,atoms)
    upper_atoms=atoms.copy(); upper_atoms.positions*=1.002
    text=lower_native.read_text()
    rows=re.compile(r'(?m)^\s*(\d+)\s+(\d+)\s+0\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$')
    text=rows.sub(lambda m:f'{int(m[1]):7d}{int(m[2]):11d}           0 '+
                  ' '.join(f'{x:12.6f}' for x in upper_atoms.positions[int(m[1])-1]),text)
    upper_native=tmp_path/'qz.log';upper_native.write_text(text)
    result=thermal.attach_vpt2_upper(p,parsed([1,2,0],'cc-pVQZ',upper_shift),upper_native,
        upper_atoms,power=3.7,evidence={'source_sha256':'qz-native'})
    expected=upper_shift+extrapolation_coefficient(4,power=3.7)*(upper_shift-lower_shift)
    assert result==pytest.approx(np.asarray(p.freq)+expected,abs=.01)
    proof=p.anl_thermochemistry_frequency_source
    assert proof['mode_map']==[3,2,1]
    assert proof['anharmonic_extrapolation']['upper_mode_map']==[3,1,2]
    assert proof['anharmonic_extrapolation']['upper_source']['source_sha256']=='qz-native'
    thermal.update_model(p)
    assert p.anl_thermochemistry_frequencies==pytest.approx(result)
    # A wrong method is not a CBS pair just because the vectors overlap.
    wrong=parsed([1,2,0],'cc-pVQZ',upper_shift);wrong['method']='B3LYP'
    with pytest.raises(ValueError,match='methods'):
        thermal.attach_vpt2_upper(p,wrong,upper_native,upper_atoms,power=3.7,evidence={})


def test_rotor_demotion_reprojects_both_shift_sets_without_changing_l2():
    p,hessian,_=rotor_fixture('methanol_co')
    p.freq,p.reduced_freqs,q,_=frequencies.get_frequencies(p,hessian,p.geom,return_modes=True)
    l2=deepcopy(p.__dict__)
    def shifts(value):
        return dict(shifts=[value]*len(p.freq),harmonic=list(p.freq),mass_weighted_modes=q.tolist(),
                    mode_map=None,match_tolerance_cm_inverse=100.,method='B2PLYP',basis='cc-pVTZ')
    lower=shifts(-2.);lower.update(upper=shifts(-3.),extrapolation_power=3.7,upper_evidence={})
    lower['upper']['basis']='cc-pVQZ'
    p.anl_thermal_model=dict(run='fixture',atom_map=list(range(p.natom)),geometry=None,
        l2_geometry_sha256=array_fingerprint(p.geom),vpt2=lower,
        harmonic=dict(positions=p.geom.tolist(),frequencies=p.freq,mass_weighted_modes=q.tolist()))
    thermal.update_model(p)
    shift=-3.-extrapolation_coefficient(4,power=3.7)
    assert p.anl_thermochemistry_frequencies==pytest.approx(np.asarray(p.reduced_freqs)+shift)
    assert p.rotor_projection==l2['rotor_projection']
    p.rotor_projection['rotors'][0]['projected']=False
    thermal.update_model(p)
    assert p.anl_thermochemistry_frequencies==pytest.approx(np.asarray(p.freq)+shift)


def test_retained_l2_harmonics_use_vectors_for_the_same_paired_shifts():
    p=point(); atoms=read(FIX/'minimum.l2.xyz'); p.geom=atoms.positions
    native=FIX/'minimum.gaussian_modes.txt'
    f,_=thermal.gaussian_modes(native.read_text(),list(p.atom),p.geom)
    p.freq=f.tolist(); p.reduced_freqs=f.tolist()
    del p.anl_thermal_model
    parsed=dict(harmonic_fundamentals_cm_inverse=f,
        anharmonic_fundamentals_cm_inverse=f-20,
        method='B2PLYP',basis='cc-pVTZ',dispersion='GD3BJ')
    thermal.attach_l2_harmonic_reference(p,parsed,native,atoms,run='fixture')
    thermal.attach_vpt2(p,parsed,native,atoms)
    upper={**parsed,'basis':'cc-pVQZ','anharmonic_fundamentals_cm_inverse':f-22}
    result=thermal.attach_vpt2_upper(p,upper,native,atoms,power=3.7,evidence={})
    assert result==pytest.approx(f-22-2*extrapolation_coefficient(4,power=3.7),abs=.01)
    assert p.anl_thermochemistry_frequency_source['model'].startswith('l2_harmonic')
    p.freq[0]+=5
    with pytest.raises(ValueError,match='retained accepted L2'):
        thermal.attach_l2_harmonic_reference(p,parsed,native,atoms,run='fixture')


@pytest.mark.parametrize('tier', ['ANL1', 'ANL1-F12'])
def test_recipe_assembly_retains_tz_cc_but_uses_separate_qz_vpt2_geometry(tmp_path,monkeypatch,tier):
    from dataclasses import replace
    from kinbot.anl import protocol, validation as v
    from kinbot.anl.model import ComponentResult
    from tests.test_anl_live_protocol import specification, fixed_resources
    monkeypatch.setattr(dispatch,'assign_partitions',fixed_resources)
    spec=protocol.replacement(specification(tier=tier),
                              dict(geometry='tz',harmonic='tz',omit_components=['higher_order_tz']))
    state={'tasks':{t['id']:{'status':'complete'} for t in spec['tasks']}}
    monkeypatch.setattr(v,'_load',lambda run:(tmp_path,spec,state))
    monkeypatch.setattr(v,'molecule_from_completed_run',
                        lambda *a:{'source':{'geometry_sha256':'3'*64}})
    equation=v.versioned_recipe(tier,vpt2_method='B2PLYP-D3BJ',vpt2_cbs=True,
        component_version=spec['intent']['component_version'],multiplicity=1,atomic=False,
        augmented=spec['intent']['anl_augmented'],
        f12_augmented=True,omitted_families=('higher_order_tz',))
    equation=replace(equation,requirements=tuple(replace(r,geometry_role='vpt2')
        if r.key=='vpt2_correction' else r for r in equation.requirements))
    equation=protocol.equation(equation,spec)
    components={r.key:ComponentResult(key=r.key,value_hartree=.001,quantity=r.quantity,
        method=r.method,basis=r.basis,backend=r.backends[0],state_id='test',charge=0,multiplicity=1,
        geometry_sha256='4'*64 if r.key=='vpt2_correction' else '3'*64,
        source_sha256='5'*64,source='verified-native-fixture',settings=dict(r.settings))
        for r in equation.requirements}
    calls=[]
    def task(run,ident,*,key,state_id):
        calls.append(ident);return components[key]
    monkeypatch.setattr(v,'task_component',task)
    monkeypatch.setattr(v,'cbs_task_component',lambda *a,**k:components[k['requirement'].key])
    monkeypatch.setattr(v,'_dual_f12_reference_component',lambda *a,**k:components['reference_cbs'])
    monkeypatch.setattr(v,'core_valence_task_component',lambda *a,**k:components['core_valence_cbs'])
    monkeypatch.setattr(v,'scalar_relativistic_task_component',lambda *a,**k:components['scalar_relativistic'])
    monkeypatch.setattr(anharmonic,'component',lambda *a,**k:components['vpt2_correction'])
    result=v._assemble_profiled_anl(tier,tmp_path,tmp_path,tmp_path,state_id='test',
        spin_orbit_hartree=0.,spin_orbit_source='spin-free',omitted_families=('higher_order_tz',))
    assert result['zero_point_hartree']==pytest.approx(.002)
    assert result['geometry_sha256']['l3']==result['geometry_sha256']['harmonic']=='3'*64
    assert result['geometry_sha256']['vpt2']=='4'*64
    assert 'harmonic_tz' in calls and 'harmonic_qz' not in calls
    assert not {'ccsdt_tz','ccsdtq_tz'} & set(calls)
    assert {'ccsdtq_dz','ccsdtqp_dz'} <= set(calls)
    assert 'cc-pVTZ' in result['components']['harmonic_zpe']['basis']
    assert 'CBS(cc-pVTZ,cc-pVQZ)'==result['components']['vpt2_correction']['basis']
