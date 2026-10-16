from copy import deepcopy
import numpy as np
import pytest

from kinbot.constants import AUtoKCAL
from kinbot.orbital_geometry import audit_symmetry, project_symmetry
from kinbot.orbital_validation import assess_path, compare_potentials


def test_declared_mirror_projection_and_element_identity():
    symbols = ['O', 'H', 'H']
    xyz = [[0.,0.,0.],[.8,.001,.6],[-.8,0.,.6]]
    op = {'name':'mirror','matrix':[[-1,0,0],[0,1,0],[0,0,1]],
          'origin_angstrom':[0,0,0],'permutation_0based':[0,2,1]}
    assert not audit_symmetry(symbols,xyz,[op],1.e-8)['passed']
    fixed, report = project_symmetry(symbols,xyz,[op],maximum_displacement_angstrom=.01)
    assert report['passed'] and not report['orbital_continuity_verified']
    assert np.linalg.norm(fixed-np.array(xyz),axis=1).max() < .001
    bad=deepcopy(op);bad['permutation_0based']=[1,0,2]
    with pytest.raises(ValueError,match='element-preserving'):
        audit_symmetry(symbols,xyz,[bad],1.e-8)
    with pytest.raises(ValueError,match='displacement'):
        project_symmetry(symbols,xyz,[op],maximum_displacement_angstrom=1.e-6)


def fixture():
    points={p:{'distance_angstrom':r,'geometry_angstrom':[['H',0,0,0],['H',r,0,0]]}
            for p,r in [('far',30.),('axial',3.),('offaxis',3.)]}
    def model(shift, differences):
        return {p:{'geometry_angstrom':v['geometry_angstrom'],
                   'energies_hartree':{'a':shift+differences[p]/AUtoKCAL,
                                       'b':shift+(differences[p]+1)/AUtoKCAL}}
                for p,v in points.items()}
    hi=model(-100,{'far':0,'axial':-5,'offaxis':-2})
    lo=model(-99,{'far':0,'axial':-3,'offaxis':.05})
    return points,lo,hi


def test_radial_bias_can_be_corrected_but_angular_error_remains():
    p,lo,hi=fixture()
    r=compare_potentials(p,lo,hi,asymptote='far',angular_reference={'offaxis':'axial'},tolerance_kcal_mol=.1)
    assert r['status']=='angular_test_passed'
    assert r['angular_residuals'][0]['maximum_absolute_residual_kcal_mol']==pytest.approx(.05)
    assert not r['production_approved']
    lo['offaxis']['energies_hartree']['a']+=.2/AUtoKCAL
    assert compare_potentials(p,lo,hi,asymptote='far',angular_reference={'offaxis':'axial'},tolerance_kcal_mol=.1)['status']=='angular_test_failed'


def test_common_zero_retains_asymptotic_state_errors_and_missing_data():
    p,lo,hi=fixture();lo['far']['energies_hartree']['b']+=.5/AUtoKCAL
    r=compare_potentials(p,lo,hi,asymptote='far',angular_reference={'offaxis':'axial'},tolerance_kcal_mol=.1)
    far=next(x for x in r['raw_comparisons'] if x['point']=='far')
    assert far['candidate_minus_reference_kcal_mol']['b']==pytest.approx(.5)
    del hi['offaxis']
    r=compare_potentials(p,lo,hi,asymptote='far',angular_reference={'offaxis':'axial'},tolerance_kcal_mol=.1)
    assert r['status']=='incomplete' and r['missing']==[{'model':'reference','point':'offaxis'}]


def test_mismatched_geometry_state_and_distance_fail():
    p,lo,hi=fixture();hi['axial']=deepcopy(hi['axial']);hi['axial']['geometry_angstrom'][1][1]=3.1
    with pytest.raises(ValueError,match='geometries'):
        compare_potentials(p,lo,hi,asymptote='far',angular_reference={},tolerance_kcal_mol=.1)
    p,lo,hi=fixture();del hi['far']['energies_hartree']['b']
    with pytest.raises(ValueError,match='manifolds'):
        compare_potentials(p,lo,hi,asymptote='far',angular_reference={},tolerance_kcal_mol=.1)
    p,lo,hi=fixture();p['offaxis']['distance_angstrom']=4
    with pytest.raises(ValueError,match='same separation'):
        compare_potentials(p,lo,hi,asymptote='far',angular_reference={'offaxis':'axial'},tolerance_kcal_mol=.1)


def test_assessment_requires_unchanged_native_evidence(tmp_path):
    import json
    from kinbot.orbital_validation import assess
    from kinbot.execution import file_sha256
    p,lo,hi=fixture()
    native=tmp_path/'native.out';native.write_text('synthetic energy evidence')
    for model in (lo,hi):
        for row in model.values():
            row['evidence']=[{'path':'native.out','sha256':file_sha256(native)}]
    payload={'schema':1,'points':p,'models':{'small':lo,'comparator':hi},'comparisons':[
        {'id':'arbitrary_channel','question':'Does the extra orbital manifold change angular dependence?',
         'candidate':'small','reference':'comparator','points':list(p),'asymptote':'far',
         'angular_reference':{'offaxis':'axial'},'tolerance_kcal_mol':.1,
         'state_assignment':'Explicit synthetic state labels'}]}
    source=tmp_path/'input.json';source.write_text(json.dumps(payload))
    result=assess(source,tmp_path/'assessment.json')
    assert json.loads(result.read_text())['comparisons'][0]['status']=='angular_test_passed'
    native.write_text('changed evidence')
    with pytest.raises(ValueError,match='evidence changed'):
        assess(source,tmp_path/'changed.json')


def test_return_hysteresis_and_missing_intervals_are_not_hidden():
    points, model, _ = fixture()
    points['returned'] = deepcopy(points['axial'])
    model['returned'] = deepcopy(model['axial'])
    model['returned']['energies_hartree']['b'] += .02/AUtoKCAL
    options = dict(path=['far', 'axial'], asymptote='far',
        return_pairs=[['axial', 'returned']], return_tolerance_kcal_mol=.01,
        spin_pairs=[{'label':'explicit pairing', 'target':'a', 'reference':'b'}])
    result = assess_path(points, model, **options)
    assert result['status'] == 'return_test_failed'
    assert result['return_checks'][0]['maximum_absolute_difference_kcal_mol'] == pytest.approx(.02)
    assert result['spin_gaps']['axial']['explicit pairing']['change_from_far_kcal_mol'] == pytest.approx(0, abs=1.e-9)
    assert not result['state_tracking_verified']
    points['middle'] = {'distance_angstrom':8., 'geometry_angstrom':[['H',0,0,0],['H',8,0,0]]}
    options['path'] = ['far', 'middle', 'axial']
    result = assess_path(points, model, **options)
    assert result['status'] == 'incomplete' and result['missing'] == ['middle']
    assert result['radial_secants'] == []


def test_path_rejects_changed_return_geometry_and_inconsistent_states():
    points, model, _ = fixture()
    options = dict(path=['far', 'axial'], asymptote='far',
                   return_pairs=[['axial','offaxis']], return_tolerance_kcal_mol=.01)
    points['offaxis']['geometry_angstrom'][1][2] = .1
    with pytest.raises(ValueError, match='identical nuclear'):
        assess_path(points, model, **options)
    points, model, _ = fixture()
    del model['axial']['energies_hartree']['b']
    with pytest.raises(ValueError, match='manifolds differ'):
        assess_path(points, model, **options)
    points, model, _ = fixture()
    points['axial']['distance_angstrom'] = float('nan')
    with pytest.raises(ValueError, match='finite point separations'):
        assess_path(points, model, **options)
