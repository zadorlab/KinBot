"""Optical decisions under the actual represented motions, without QC calls."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.io import read

from kinbot.conformer_counting import evaluate_members
from kinbot.conformer_records import ConformerRecord
from kinbot.counting_contract import optical_counting
from kinbot.mess import MESS
from kinbot.optical import _Parts, rigid_mirror, evaluate_optical
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.thermochemistry import hir_evidence, thermochemistry_evidence
from tests.counting_fixtures import saved_point, methanol_data
from test_hir_optical_domains import pyramidal_product, puckered_ts, relaxed_radical, rigid_mmff_scan


@pytest.mark.parametrize('factory,factor', [(pyramidal_product,2),
    (lambda:puckered_ts()[0],2),(relaxed_radical,1),
    (lambda:saved_point(methanol_data()['saddles'][0]),1)])
def test_saved_chemistry_and_global_coordinate_invariances(factory,factor):
    p=factory(); rotors=hir_evidence(p)['rotors']
    assert evaluate_optical(p,rotors=rotors)['remaining_multiplier']==factor
    rng=np.random.default_rng(177)
    order=rng.permutation(p.natom); reverse=np.argsort(order)
    rotation=np.linalg.qr(rng.normal(size=(3,3)))[0]
    rotation[:,0]*=np.linalg.det(rotation)
    p.atom=np.asarray(p.atom)[order];p.geom=p.geom[order]@rotation+4.
    p.bond=p.bond[np.ix_(order,order)]
    p.bonds=[b[np.ix_(order,order)] for b in p.bonds]
    p.rads=[r[order] for r in p.rads]
    if hasattr(p,'reac_bond'):p.reac_bond=p.reac_bond[np.ix_(order,order)]
    for r in rotors:
        r['axis']=list(reverse[r['axis']])[::-1]
        r['dihedral']=list(reverse[r['dihedral']])[::-1]
        for point in r['points']:
            point['geometry_angstrom']=(np.asarray(point['geometry_angstrom'])[order]@rotation+4.).tolist()
    assert evaluate_optical(p,rotors=rotors)['remaining_multiplier']==factor


def test_active_rotors_together_cover_mirror_and_inactive_ones_do_not():
    p=relaxed_radical();e=hir_evidence(p)
    result=evaluate_optical(p,rotors=e['rotors'])
    assert result['remaining_multiplier']==1
    assert result['coordinate_coverage']['maximum_displacement_angstrom']<.1
    assert result['coordinate_coverage']['additional_qc_calculations']==0
    assert evaluate_optical(p)['remaining_multiplier']==2


def test_handed_core_cannot_disappear_inside_a_long_planar_tail():
    # Analytic tetrahedron attached to an arbitrarily long planar substituent.
    xyz=[[0,0,0],[0,1,0],[0,0,1],[.8,-.5,-.5],[0,-1,0]]
    xyz += [[0,-2-i,.1*(i%2)] for i in range(120)]
    atoms=['C','H','F','Cl']+['C']*(len(xyz)-4)
    bond=np.zeros((len(xyz),len(xyz)),int)
    for a,b in [(0,1),(0,2),(0,3),(0,4)]+[(i,i+1) for i in range(4,len(xyz)-1)]:
        bond[a,b]=bond[b,a]=1
    p=SimpleNamespace(atom=np.array(atoms),geom=np.array(xyz,float),bond=bond,bonds=[bond],rads=[])
    # No stereo assignment is assumed for this abstract local-geometry fixture.
    assert _Parts(p,[]).compare(p.geom,p.geom,.1)['status']=='distinct'


def ethanol_near_planar(angle=6.):
    atoms=read(Path(__file__).parents[1]/'tests/reference/ethanol_gauche_mirror.xyz')
    atoms.set_dihedral(0,1,2,8,angle,indices=[2,8])
    p=StationaryPoint('ethanol',0,1,atom=atoms.get_chemical_symbols(),geom=atoms.positions)
    p.characterize()
    return p


def record(p,index,geometry=None):
    return ConformerRecord(f'{p.name}:{index}',index,None,'valid',
        geometry=tuple(map(tuple,p.geom if geometry is None else geometry)),
        zero_energy_hartree=-100.,frequencies_cm1=(100.,)*max(0,3*p.natom-6))


@pytest.mark.parametrize('strict',[False,True])
def test_inconclusive_shape_keeps_weight_one_in_mc_with_unresolved_evidence(strict):
    p=ethanol_near_planar()
    assert rigid_mirror(p)['mirror_states'] is None
    records,groups=evaluate_members(p,[record(p,0)],strict=strict)
    assert groups==[[0]] and records[0].remaining_optical_weight == 1
    assert records[0].optical_evidence['status']=='unresolved'
    par=dict(multi_conf_tst=1,pes=0,freq_uq_ref=1000.,freq_uq_max_exp=1.)
    writer=MESS(par,p)
    text=writer._member_rrho(p,records[0],1.,0.)
    assert 'WARNING: unresolved symmetry number' in text
    assert 'using optical factor 1' in text


def test_local_threshold_is_01_and_failed_matches_remain_unresolved():
    from kinbot.molecular_symmetry import OPTICAL_RMSD_TOLERANCE
    from kinbot.optical import SCAN_MIRROR_RMSD_TOLERANCE
    assert OPTICAL_RMSD_TOLERANCE==.1 and SCAN_MIRROR_RMSD_TOLERANCE==.1
    assert rigid_mirror(ethanol_near_planar(4.))['status']=='match'
    result=rigid_mirror(ethanol_near_planar(6.))
    assert result['status']=='undetermined'
    assert result['largest_local_rmsd_angstrom']>.1


def test_mc_inconclusive_pair_is_not_doubled_twice():
    p=ethanol_near_planar()
    first=Atoms(p.atom,positions=p.geom)
    first.set_dihedral(0,1,2,8,60.,indices=[2,8])
    p.geom=first.positions.copy()
    second=Atoms(p.atom,positions=p.geom*[-1,1,1])
    second.rotate_dihedral(0,1,2,8,12.5,indices=[2,8])
    pair=[record(p,0),record(p,1,second.positions)]
    records,groups=evaluate_members(p,pair)
    assert len(groups)==2 and all(r.remaining_optical_weight==1 for r in records)
    assert all(r.optical_evidence['status']=='unresolved' for r in records)


def test_unknown_stereo_scope_uses_explicit_legacy_outcome():
    p=pyramidal_product();p.hir=None;p.reduced_freqs=p.freq;p.rotor_projection=None
    p.optical_reference={'status':'unavailable'}
    result=optical_counting(p,hir_evidence(p))
    assert result['status']=='legacy_unverified'
    p.optical_reference={'status':'unsupported','reason':'unrepresented configuration'}
    result=optical_counting(p,hir_evidence(p))
    assert result['status']=='legacy_unverified'
    assert result['remaining_multiplier'] > 0


@pytest.mark.parametrize('value',[{'x':1},{'x':{'multiplier':3,'reason':'x'}},
    {'x':{'multiplier':1,'reason':''}},{'x':{'multiplier':True,'reason':'x'}}])
def test_assumption_parameter_rejects_invalid_entries(tmp_path,value):
    path=tmp_path/'input.json';path.write_text(json.dumps({'optical_factor_assumptions':value}))
    with pytest.raises(ValueError,match='optical_factor_assumptions'):
        Parameters(str(path),show_warnings=False)


def test_selected_mc_mirror_weights_have_consistent_evidence():
    p=pyramidal_product()
    records,_=evaluate_members(p,[record(p,0),record(p,1,p.geom*[-1,1,1])])
    for r in records:
        assert r.remaining_optical_weight==1
        assert r.optical_evidence['remaining_multiplier']==1
        assert r.optical_evidence['reason']=='explicit mirror'
        assert r.optical_evidence['states_covered_by_mc_tst']==2


def loose_ts():
    data=json.loads((Path(__file__).parent/'reference/loose_ts_optical.json').read_text())
    p=saved_point(data)
    p.name='751992412091250510002_r12_insertion_R_3_2_1'
    return p


@pytest.mark.parametrize('pes',[0,1])
def test_explicit_assumption_remains_available_without_bypassing_bad_data(pes):
    from kinbot.species_routing import routing_name
    p=ethanol_near_planar();p.reduced_freqs=list(p.freq)
    assert thermochemistry_evidence(p)['optical_counting']['status']=='unresolved'
    par=dict(multi_conf_tst=0,pes=pes,optical_factor_assumptions={routing_name(p):{
        'multiplier':1,'reason':'Explicit optical approximation.'}})
    writer=MESS(par,p)
    assert writer._parent_symmetry(p)==p.sigma_ext
    result=thermochemistry_evidence(p)['optical_counting']
    assert result['status']=='assumed' and result['remaining_multiplier']==1
    assert result['states_covered_by_hir'] is None
    assert result['rigid_mirror']['status']=='undetermined'
    assert 'explicit assumption' in writer._optical_comment(p)
    # An assumption cannot conceal stale calculation data.
    p=loose_ts();p.hir.scan_reference['source_row_id']=-1
    result=thermochemistry_evidence(p)['optical_counting']
    assert result['status']=='unresolved' and result['remaining_multiplier'] is None


def test_fixed_configuration_controls_hir_population_and_assumption_scope():
    p=rigid_mmff_scan('C[C@H](F)Cl')
    p.optical_factor_assumption=dict(multiplier=2,reason='Not permission to add another population.')
    assert optical_counting(p,hir_evidence(p))['remaining_multiplier']==1
    p.optical_population='racemic'
    result=optical_counting(p,hir_evidence(p))
    assert result['remaining_multiplier']==2
    assert result['reason']=='The represented torsions retain the assigned fixed configuration.'


def test_torsional_coverage_does_not_invent_a_resolved_rigid_pair_size():
    p=rigid_mmff_scan('CCO')
    atoms=Atoms(p.atom,positions=p.geom)
    # RDKit's explicit atom ordering here has the hydroxyl H last.
    atoms.set_dihedral(0,1,2,8,8.,indices=[2,8])
    p.geom=atoms.positions
    result=evaluate_optical(p,rotors=hir_evidence(p)['rotors'])
    assert result['total_optical_states'] is None
    assert result['remaining_multiplier']==1
    assert result['states_covered_by_hir'] is None


def test_preflight_permits_unresolved_states_and_ignores_placeholder_ts(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    p=ethanol_near_planar();p.name='parent';p.reduced_freqs=list(p.freq)
    q=copy.deepcopy(p);q.name='product'
    placeholder=SimpleNamespace(name='hom_sci_placeholder')  # no physical TS properties
    p.reac_type=['hom_sci']
    p.reac_obj=[SimpleNamespace(instance_name='hom_sci',ts=placeholder,
                              prod_opt=[SimpleNamespace(species=q)],do_vdW=False)]
    writer=MESS(dict(multi_conf_tst=0,pes=0),p)
    writer._check_optical_models(['hom_sci'],['hom_sci'])
    for state in (p,q):
        assert state.mess_optical_counting['status']=='unresolved'
        assert state.mess_optical_counting['remaining_multiplier']==1
    assert not hasattr(placeholder,'mess_optical_counting')


@pytest.mark.parametrize('pes',[0,1])
def test_unresolved_well_is_written_with_weight_one_warning_and_export(tmp_path,monkeypatch,caplog,pes):
    from kinbot import symmetry
    from kinbot.species_routing import routing_key
    p=ethanol_near_planar();p.freq=p.reduced_freqs=[100.]*(3*p.natom-6)
    p.energy=-100.;p.zpe=.1
    symmetry.calculate_symmetry(p)
    path=tmp_path/'input.json';path.write_text(json.dumps({'barrier_threshold':100.}))
    par=Parameters(str(path),show_warnings=False).par
    par.update(multi_conf_tst=0,pes=pes)
    writer=MESS(par,p);writer.well_names={routing_key(p):'w1'}
    monkeypatch.chdir(tmp_path)
    text=writer.write_well(p,0.,1.,0)
    assert '! WARNING: unresolved symmetry number' in text
    assert 'using optical factor 1' in text
    assert 'unresolved symmetry number' in caplog.text
    assert writer._parent_symmetry(p)==p.sigma_ext
    result=thermochemistry_evidence(p)['optical_counting']
    assert result['status']=='unresolved' and result['remaining_multiplier']==1
    assert result['states_covered_by_hir'] is None
    assert result['local_rmsd_tolerance_angstrom']==.1
    assert result['measured_mirror_rmsd_tolerance_angstrom']==.1


def test_saved_loose_ts_uses_explicit_midpoint_approximation():
    result=thermochemistry_evidence(loose_ts())['optical_counting']
    assert result['status']=='assumed' and result['remaining_multiplier']==1
    assert .1<result['coordinate_coverage']['largest_local_rmsd_angstrom']<.25
    assert result['harmonic_midpoint']['status']=='complete'
    assert result['harmonic_midpoint']['cutoff_kcal_mol']==4.
    assert result['heuristic']=='harmonic_midpoint'
