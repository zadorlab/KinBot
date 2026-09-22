"""Analytic quadratic checks and saved-chemistry midpoint diagnostics; no QC."""
import copy
import json
from pathlib import Path

import numpy as np
import pytest
from ase import units

from kinbot import constants, frequencies, geometry
from kinbot.counting_contract import optical_counting
from kinbot.optical_harmonic import (harmonic_midpoint_energy,
    mirror_midpoint_diagnostic, selected_midpoint_diagnostic)
from kinbot.stationary_pt import StationaryPoint
from kinbot.thermochemistry import hir_evidence
from tests.counting_fixtures import saved_point
from test_optical_evaluator import loose_ts


def quadratic(values=(.1,.2,.3), displacement=(.2,.1,.3)):
    p=StationaryPoint('quadratic',0,1,atom=['O','H','H'],
                     geom=np.array([[0.,0.,0.],[.9,0.,0.],[-.2,.9,0.]]))
    p.characterize()
    mass=np.repeat([constants.exact_mass[a] for a in p.atom],3)
    centered=p.geom-geometry.get_center_of_mass(p.geom,p.atom)
    t,r=frequencies.rigid_body_vectors(centered,p.atom)
    _,_,v=np.linalg.svd(np.vstack((t,r)),full_matrices=True)
    basis=v[6:]
    mw=basis.T@np.diag(values)@basis
    hessian=mw*np.sqrt(np.outer(mass,mass))
    delta=(basis.T@np.asarray(displacement))*units.Bohr/np.sqrt(mass)
    return p,p.geom+delta.reshape((-1,3)),hessian,mw


def test_analytic_midpoint_has_one_eighth_not_one_half_of_full_displacement_cost():
    p,mirror,hessian,_=quadratic()
    result=harmonic_midpoint_energy(p,mirror,hessian,hessian_unit='hartree / bohr^2')
    expected=np.dot([.1,.2,.3],np.square([.2,.1,.3]))/8*constants.AUtoKCAL
    assert result['stable_midpoint_energy_kcal_mol']==pytest.approx(expected)
    assert result['signed_vibrational_energy_kcal_mol']==pytest.approx(expected)
    assert result['external_rank']==6


def test_cartesian_massweighted_and_ev_hessians_give_the_same_energy():
    p,mirror,cartesian,mw=quadratic()
    matrices=[(cartesian,'hartree / bohr^2'),(mw,'hartree / (bohr^2 * amu)'),
              (cartesian*units.Hartree/units.Bohr**2,'eV / angstrom^2')]
    results=[harmonic_midpoint_energy(p,mirror,h,hessian_unit=u) for h,u in matrices]
    assert all(r['stable_midpoint_energy_kcal_mol']==pytest.approx(
        results[0]['stable_midpoint_energy_kcal_mol']) for r in results)
    with pytest.raises(ValueError,match='Unknown Hessian units'):
        harmonic_midpoint_energy(p,mirror,cartesian,hessian_unit='unknown')


def test_negative_ts_mode_does_not_cancel_stable_displacement_energy():
    p,mirror,hessian,_=quadratic((-.2,.2,.3),(1.,1.,0.))
    p.wellorts=1
    result=harmonic_midpoint_energy(p,mirror,hessian,hessian_unit='hartree / bohr^2')
    expected=.2/8*constants.AUtoKCAL
    assert result['stable_midpoint_energy_kcal_mol']==pytest.approx(expected)
    assert result['negative_mode_energy_kcal_mol']==pytest.approx(-expected)
    assert abs(result['signed_vibrational_energy_kcal_mol'])<1.e-10
    assert result['negative_mode_displacement_fraction']==pytest.approx(.5)
    assert result['negative_mode_count']==1
    assert result['transition_surface_preserved']=='not established'


def test_coordinate_frame_and_atom_permutation_do_not_change_midpoint_energy():
    p,mirror,hessian,_=quadratic()
    original=harmonic_midpoint_energy(p,mirror,hessian,hessian_unit='hartree / bohr^2')
    rng=np.random.default_rng(162)
    rotation=np.linalg.qr(rng.normal(size=(3,3)))[0]
    rotation[:,0]*=np.linalg.det(rotation)
    transform=np.kron(np.eye(3),rotation)
    p.geom=p.geom@rotation.T+3.;mirror=mirror@rotation.T+3.
    hessian=transform@hessian@transform.T
    order=np.array([2,0,1]);indices=(3*order[:,None]+np.arange(3)).ravel()
    p.geom=p.geom[order];p.atom=np.array(p.atom)[order];mirror=mirror[order]
    p.bond=p.bond[np.ix_(order,order)]
    result=harmonic_midpoint_energy(p,mirror,hessian[np.ix_(indices,indices)],
                                  hessian_unit='hartree / bohr^2')
    assert result['stable_midpoint_energy_kcal_mol']==pytest.approx(
        original['stable_midpoint_energy_kcal_mol'])


def test_saved_rotor_construction_is_invariant_to_axis_direction_and_rotor_order():
    p=loose_ts();rotors=hir_evidence(p)['rotors']
    before=p.geom.copy();hessian=np.array(p.hess).copy()
    original=mirror_midpoint_diagnostic(p,p.hess,hessian_unit='hartree / bohr^2',rotors=rotors)
    reverse=copy.deepcopy(rotors[::-1])
    for rotor in reverse:
        rotor['dihedral']=rotor['dihedral'][::-1];rotor['axis']=rotor['axis'][::-1]
    result=mirror_midpoint_diagnostic(p,p.hess,hessian_unit='hartree / bohr^2',rotors=reverse)
    assert result['global_rmsd_angstrom']==pytest.approx(original['global_rmsd_angstrom'])
    assert result['stable_midpoint_energy_kcal_mol']==pytest.approx(original['stable_midpoint_energy_kcal_mol'])
    np.testing.assert_allclose(result['aligned_mirror_geometry_angstrom'],original['aligned_mirror_geometry_angstrom'])
    np.testing.assert_array_equal(p.geom,before)
    np.testing.assert_array_equal(p.hess,hessian)
    assert result['additional_qc_calculations']==0 and result['cutoff_kcal_mol'] is None


def test_capped_mapping_search_is_not_reported_as_exhaustive():
    p=loose_ts()
    result=mirror_midpoint_diagnostic(p,p.hess,hessian_unit='hartree / bohr^2',
                                     rotors=hir_evidence(p)['rotors'],max_mappings=1)
    assert result['status']=='limited' and not result['mapping_search_complete']
    assert result['mappings_tested']==1 and not result['optical_weight_changed']


@pytest.mark.parametrize('change',['absent','geometry','hessian','job','row','unit'])
def test_selected_diagnostic_never_combines_incompatible_calculation_properties(change):
    p=loose_ts();rotors=hir_evidence(p)['rotors']
    if change=='absent':p.hess=[]
    elif change=='geometry':p.geom[0,0]+=.01
    elif change=='hessian':p.hess[0][0]+=.01
    elif change=='job':p.source_job='another-job'
    elif change=='row':p.source_row_id+=1
    else:p.rotor_projection['reference']['hessian_unit']='unknown'
    result=selected_midpoint_diagnostic(p,rotors)
    assert result['status']=='unavailable' and not result['optical_weight_changed']


@pytest.mark.parametrize('fixture',['pyramidal_product_hir','puckered_ts_hir',
                                  'relaxed_radical_hir','loose_ts_optical'])
def test_corrected_fixture_units_reproduce_saved_spectrum_and_do_not_change_count(fixture):
    data=json.loads((Path(__file__).parent/'reference'/f'{fixture}.json').read_text())
    assert 'hessian_native_eV_angstrom-2' not in data
    p=saved_point(data);rotors=hir_evidence(p)['rotors']
    original=optical_counting(p,hir_evidence(p))
    result=mirror_midpoint_diagnostic(p,p.hess,hessian_unit='hartree / bohr^2',rotors=rotors)
    modes=[m['frequency_cm1'] for m in result['mode_contributions']]
    np.testing.assert_allclose(sorted(modes),sorted(p.freq),atol=1.e-7)
    assert optical_counting(p,hir_evidence(p))['remaining_multiplier']==original['remaining_multiplier']
    assert result['cutoff_kcal_mol'] is None and not result['optical_weight_changed']


def test_midpoint_quadratic_is_not_a_double_well_barrier():
    # V=a(x^2-d^2)^2 has curvature 8ad^2 at its two minima. Their
    # midpoint harmonic estimate is four times the actual barrier ad^4.
    a,d=.2,.7
    curvature=8*a*d*d
    true_barrier=a*d**4
    midpoint_estimate=.5*curvature*d*d
    assert midpoint_estimate==pytest.approx(4*true_barrier)


def test_linear_molecule_removes_five_external_motions():
    p=StationaryPoint('CO2',0,1,atom=['O','C','O'],
                     geom=np.array([[-1.16,0.,0.],[0.,0.,0.],[1.16,0.,0.]]))
    p.characterize()
    result=harmonic_midpoint_energy(p,p.geom,np.eye(9),hessian_unit='hartree / bohr^2')
    assert result['external_rank']==5
    assert len(result['mode_contributions'])==4
    assert result['stable_midpoint_energy_kcal_mol']==0.


@pytest.mark.parametrize('scale,weight',[(1.,1.),(3.,2.)])
def test_mess_midpoint_approximation_divisor_and_comment_agree(caplog,scale,weight):
    from kinbot.mess import MESS
    from kinbot.calculation import array_fingerprint
    p=loose_ts()
    # Scaling is an analytic writer test, not another physical calculation.
    p.hess=(np.asarray(p.hess)*scale).tolist()
    p.freq=(np.asarray(p.freq)*np.sqrt(scale)).tolist()
    p.reduced_freqs=(np.asarray(p.reduced_freqs)*np.sqrt(scale)).tolist()
    p.rotor_projection['reference']['hessian_sha256']=array_fingerprint(p.hess)
    writer=MESS(dict(multi_conf_tst=0,optical_population='specified'),p)
    assert writer._parent_symmetry(p)==p.sigma_ext/weight
    text=writer._optical_comment(p)
    for message in ('Harmonic midpoint approximation',f'optical factor {weight:g}',
                    '4 kcal/mol',
                    'counting heuristic','not an inversion barrier'):
        assert message in text and message in caplog.text
    assert 'unresolved symmetry number' not in text


@pytest.mark.parametrize('energy,weight',[(0.,1.),(3.999,1.),(4.,1.),(4.001,2.),(50.,2.)])
def test_heuristic_threshold_and_no_claim_of_measured_coverage(energy,weight):
    from kinbot.optical_harmonic import apply_midpoint_heuristic
    count=dict(status='unresolved',fallback='unresolved_symmetry',remaining_multiplier=1.,
        reason='Rigid mirror comparison is numerically undetermined.',
        population_scope={'mirror_allowed':True},states_covered_by_hir=None)
    diagnostic=dict(status='complete',mapping_search_complete=True,
                    stable_midpoint_energy_kcal_mol=energy)
    result=apply_midpoint_heuristic(count,diagnostic)
    assert result['status']=='assumed' and result['remaining_multiplier']==weight
    assert result['states_covered_by_hir'] is None and 'fallback' not in result
    assert count['status']=='unresolved' and 'cutoff_kcal_mol' not in diagnostic


@pytest.mark.parametrize('change',['fixed_population','resolved','assumption','energy_conflict',
    'overlapping_rotors','mc_pair','missing','limited','nan','negative','infinite','bool','text'])
def test_heuristic_cannot_override_reliable_counts_conflicts_or_missing_evidence(change):
    from kinbot.optical_harmonic import apply_midpoint_heuristic
    count=dict(status='unresolved',fallback='unresolved_symmetry',remaining_multiplier=1.,
        reason='Rigid mirror comparison is numerically undetermined.',
        population_scope={'mirror_allowed':True})
    diagnostic=dict(status='complete',mapping_search_complete=True,
                    stable_midpoint_energy_kcal_mol=50.)
    if change=='fixed_population':count['population_scope']['mirror_allowed']=False
    elif change in ('resolved','assumption'):count['status']='resolved' if change=='resolved' else 'assumed'
    elif change=='energy_conflict':count['reason']='Measured mirror geometries have contradictory energies.'
    elif change=='overlapping_rotors':count['reason']='Multiple relaxed rotor coordinates independently contain the same mirror; overlapping coverage is unresolved.'
    elif change=='mc_pair':count['reason']='Uncertain explicit mirror coverage.'
    elif change in ('missing','limited'):diagnostic['status']='unavailable' if change=='missing' else 'limited'
    else:diagnostic['stable_midpoint_energy_kcal_mol']={
        'nan':float('nan'),'negative':-1.,'infinite':float('inf'),'bool':True,'text':'4'}[change]
    result=apply_midpoint_heuristic(count,diagnostic)
    assert result['status']==count['status'] and result['remaining_multiplier']==1.
    assert 'heuristic' not in result


def test_mc_diagnostic_does_not_borrow_a_parent_hessian():
    from kinbot.conformer_counting import evaluate_members
    from test_optical_evaluator import ethanol_near_planar, record
    p=ethanol_near_planar()
    p.hess=np.eye(3*p.natom)
    records,_=evaluate_members(p,[record(p,0)])
    result=records[0].optical_evidence['harmonic_midpoint']
    assert result['status']=='unavailable'
    assert 'no own Hessian' in result['reason']
    assert records[0].remaining_optical_weight==1.
