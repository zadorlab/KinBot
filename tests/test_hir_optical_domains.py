"""Ordinary full 1D HIR models can include mirrors between grid points."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from kinbot import symmetry
from kinbot.calculation import geometry_reference, array_fingerprint
from kinbot.hindered_rotors import HIR
from kinbot.mess import MESS
from kinbot.molecular_symmetry import geometric_mirror_states
from kinbot.thermochemistry import hir_evidence
from kinbot.counting_contract import optical_counting
from test_configured_rotational_symmetry import optimized
from tests.counting_fixtures import saved_point

def rigid_mmff_scan(smi, angle=None):
 p=optimized(smi)
 mol=Chem.AddHs(Chem.MolFromSmiles(smi)); conf=Chem.Conformer(p.natom)
 for i,xyz in enumerate(p.geom): conf.SetAtomPosition(i,xyz)
 mol.AddConformer(conf)
 if angle is not None:
  rdMolTransforms.SetDihedralDeg(mol.GetConformer(),*p.dihed[0],angle)
  p.geom=mol.GetConformer().GetPositions()
 symmetry.calculate_symmetry(p)
 p.source_job='structural-rigid-scan'; p.source_row_id=1
 p.hess=np.eye(p.natom*3)
 p.freq=[500.]*(3*p.natom-6); p.reduced_freqs=p.freq[:-len(p.dihed)]
 h=p.hir=HIR(p,SimpleNamespace(qc='mmff'),{'nrotation':12,'plot_hir_profiles':False,'rotor_0_test':1})
 h.rigid_scan=True; h.scan_reference=geometry_reference(p)
 h.scan_reference.update(backend='mmff',dihedrals=copy.deepcopy(p.dihed))
 props=AllChem.MMFFGetMoleculeProperties(mol)
 for i,rotor in enumerate(p.dihed):
  base=rdMolTransforms.GetDihedralDeg(mol.GetConformer(),*rotor)
  geoms=[]; energies=[]
  for j in range(12):
   m=Chem.Mol(mol); rdMolTransforms.SetDihedralDeg(m.GetConformer(),*rotor,base+30*j)
   geoms.append(m.GetConformer().GetPositions())
   energies.append(AllChem.MMFFGetMoleculeForceField(m,props).CalcEnergy()/627.509)
  h.hir_status.append([0]*12); h.hir_energies.append(energies); h.hir_geoms.append(geoms)
  h.scan_jobs.append([f'rotor{i}-{j}' for j in range(12)]); h.point_observations.append([])
  assert h.fourier_fit('structural',np.arange(12)*np.pi/6,i)
 ref=geometry_reference(p); ref.update(dihedrals=copy.deepcopy(p.dihed),hessian_sha256=array_fingerprint(p.hess))
 p.rotor_projection={'reference':ref,'internal_rank':len(p.dihed),'rotors':[{'rotor_index':i,'projected':True} for i in range(len(p.dihed))]}
 return p


def test_two_ethanol_rotors_do_not_need_another_optical_factor():
 p=rigid_mmff_scan('CCO')
 evidence=hir_evidence(p)
 result=optical_counting(p,evidence)
 assert len(evidence['rotors']) == 2
 assert result['remaining_multiplier'] == 1.
 assert MESS({'multi_conf_tst':0},p)._parent_symmetry(p) == p.sigma_ext
 assert sum(r['mirror_coverage']['status']=='observed_pair' for r in evidence['rotors']) == 1


def test_peroxide_mirror_between_grid_points_is_in_full_torsional_domain():
 p=rigid_mmff_scan('OO',111.5)
 evidence=hir_evidence(p)
 result=optical_counting(p,evidence,tolerance=.05)
 assert result['remaining_multiplier'] == 1.
 assert not evidence['rotors'][0]['mirror_coverage']['witnesses']
 assert result['coordinate_coverage']['additional_qc_calculations'] == 0
 assert not result['coordinate_coverage']['coupled_potential_accuracy_established']
 assert MESS({'multi_conf_tst':0},p)._parent_symmetry(p) == p.sigma_ext


def test_coordinate_coverage_does_not_invent_a_missing_partial_scan_observation():
 p=rigid_mmff_scan('OO',111.5)
 p.hir.hir_status[0][5]=1
 assert p.hir.fourier_fit('partial',np.arange(12)*np.pi/6,0)
 evidence=hir_evidence(p)
 result=optical_counting(p,evidence,tolerance=.05)
 assert result['remaining_multiplier'] == 1.
 assert not evidence['rotors'][0]['mirror_coverage']['witnesses']
 assert 'coordinate_coverage' in result
 assert evidence['rotors'][0]['points'][5]['status'] != 'successful'


def pyramidal_product():
    data = json.loads((Path(__file__).parent / 'reference/pyramidal_product_hir.json').read_text())
    return saved_point(data)


@pytest.mark.parametrize('pes', [0, 1])
def test_saved_pyramidal_product_writes_hir_with_one_missing_mirror(tmp_path, monkeypatch, pes):
    from kinbot.parameters import Parameters
    from kinbot.species_routing import routing_key
    p = pyramidal_product()
    original_geom = p.geom.copy()
    original_sigma = copy.deepcopy(p.sigma_int)
    assert p.sigma_ext == 1
    assert len(p.freq) == 30 and len(p.reduced_freqs) == 27
    input_file = tmp_path / 'input.json'
    input_file.write_text(json.dumps({'barrier_threshold': 100.}))
    par = Parameters(str(input_file)).par
    par.update(multi_conf_tst=0, rotor_scan=1, pes=pes)
    writer = MESS(par, p)
    writer.well_names = {routing_key(p): 'w1'}
    monkeypatch.chdir(tmp_path)
    text = writer.write_well(p, 0., 1., 0)
    result = p.mess_optical_counting
    assert result['status'] == 'resolved'
    assert result['remaining_multiplier'] == 2.
    assert result['states_covered_by_hir'] == 1
    assert result['coordinate_coverage']['rejected_mapping_examples']
    assert result['coordinate_coverage']['successful_observations'] == 36
    assert sum(line.split() == ['Rotor', 'Hindered'] for line in text.splitlines()) == 3
    assert any(line.split() == ['SymmetryFactor', '0.5'] for line in text.splitlines())
    np.testing.assert_array_equal(p.geom, original_geom)
    np.testing.assert_array_equal(p.sigma_int, original_sigma)


@pytest.mark.parametrize('alteration', ['opposite_branch', 'planar', 'failed', 'missing_geometry', 'unknown_domain'])
def test_pyramidal_branch_evidence_must_be_complete_and_consistent(alteration):
    p = pyramidal_product()
    evidence = hir_evidence(p)
    point = evidence['rotors'][0]['points'][5]
    if alteration == 'opposite_branch':
        point['geometry_angstrom'] = (np.array(point['geometry_angstrom']) * [-1., 1., 1.]).tolist()
        # A branch change without compatible energy is not mirror coverage.
        point['electronic_energy_hartree'] += 1.
    elif alteration == 'planar':
        geom = np.array(point['geometry_angstrom'])
        neighbours = geom[[2, 4, 5]]
        normal = np.cross(neighbours[1]-neighbours[0], neighbours[2]-neighbours[0])
        normal /= np.linalg.norm(normal)
        geom[3] -= np.dot(geom[3]-neighbours[0], normal) * normal
        point['geometry_angstrom'] = geom.tolist()
    elif alteration == 'unknown_domain':
        evidence['rotors'][0]['represented_domain_degrees'] = None
    else:
        point['geometry_angstrom'] = None
    result = optical_counting(p, evidence)
    # A planar point or a different inverted scan geometry removes the
    # handedness proof. Neither is a measured mirror of the reference.
    # Missing data remain invalid; exact mirror energy conflicts have their
    # own regression in test_counting_contract.
    if alteration in ('planar','opposite_branch'):
        from kinbot.optical import evaluate_optical
        assert evaluate_optical(p,rotors=evidence['rotors'])['status']=='unresolved'
        # The new >4 rule is an approximation, not restored handedness proof.
        assert result['status']=='assumed' and result['remaining_multiplier']==2.
        assert result['heuristic']=='harmonic_midpoint'
        assert result['states_covered_by_hir'] is None
    else:
        assert result['status']=='unresolved' and result['remaining_multiplier'] is None


def test_unrepresented_pyramid_is_not_an_oxygen_specific_rule():
    p = rigid_mmff_scan('CCN(C)O')
    result = optical_counting(p, hir_evidence(p))
    assert result['remaining_multiplier'] == 2.
    assert any('N' in [p.atom[i] for i in w['atom_indices']] for w in result['coordinate_coverage']['rejected_mapping_examples'])
    # Two identical methyl substituents do not define two inversion states.
    p = rigid_mmff_scan('CCN(C)C')
    assert optical_counting(p, hir_evidence(p))['remaining_multiplier'] == 1.
    # A fixed carbon configuration excludes the global mirror by default,
    # even when an unscanned pyramidal centre is present elsewhere.
    p = rigid_mmff_scan('C[C@H](F)N(C)O')
    result = optical_counting(p, hir_evidence(p))
    assert result['allowed_global_mirror_states'] == 1
    assert result['remaining_multiplier'] == 1.


def puckered_ts():
    data = json.loads((Path(__file__).parent / 'reference/puckered_ts_hir.json').read_text())
    p = saved_point(data)
    p.name = '751992412091250510002_intra_H_migration_5_9'
    endpoints = []
    for observation in data['endpoints']:
        from kinbot.stationary_pt import StationaryPoint
        endpoint = StationaryPoint(observation['name'], observation['charge'],
            observation['multiplicity'], atom=observation['atoms'],
            geom=np.array(observation['geometry_angstrom']))
        endpoint.bond = np.array(observation['bond'])
        endpoint.bonds = [np.array(b) for b in observation['bonds']]
        endpoint.rads = [np.array(r) for r in observation['rads']]
        endpoint.characterize(bond_mx=endpoint.bond)
        endpoints.append(endpoint)
    return p, endpoints


@pytest.mark.parametrize('pes', [0, 1])
def test_saved_puckered_ts_writes_barrier_with_mirror_outside_methyl_rotor(tmp_path, monkeypatch, pes):
    from kinbot.parameters import Parameters
    from kinbot.species_routing import routing_key
    p, (reactant, product) = puckered_ts()
    original_geom = p.geom.copy()
    original_sigma = copy.deepcopy(p.sigma_int)
    original_freqs = list(p.reduced_freqs)
    assert p.dihed == [[5, 0, 1, 2]]
    assert p.sigma_int[0][1] == 3 and p.sigma_ext == 1
    assert len(p.freq) == 30 and len(p.reduced_freqs) == 29
    assert p.reduced_freqs[0] < 0 < p.reduced_freqs[1]
    input_file = tmp_path / 'input.json'
    input_file.write_text(json.dumps({'barrier_threshold': 100.}))
    par = Parameters(str(input_file)).par
    par.update(multi_conf_tst=0, rotor_scan=1, pes=pes)
    reaction = SimpleNamespace(instance_name=p.name, ts=p, products=[product], do_vdW=False)
    reactant.reac_type = ['intra_H_migration']
    writer = MESS(par, reactant)
    writer.ts_names = {p.name: 'ts1'}
    writer.well_names = {routing_key(reactant): 'w1', routing_key(product): 'w2'}
    monkeypatch.chdir(tmp_path)
    text, _ = writer.write_barrier(reaction, 0, 33., 25., 0., 1., 1., 0)
    result = p.mess_optical_counting
    assert result['remaining_multiplier'] == 2.
    assert result['allowed_global_mirror_states'] == 2
    assert result['states_covered_by_hir'] == 1
    assert result['coordinate_coverage']['successful_observations'] == 12
    assert result['coordinate_coverage']['rejected_mapping_examples']
    assert sum(line.split() == ['Rotor', 'Hindered'] for line in text.splitlines()) == 1
    assert any(line.split() == ['SymmetryFactor', '0.5'] for line in text.splitlines())
    np.testing.assert_array_equal(p.geom, original_geom)
    np.testing.assert_array_equal(p.sigma_int, original_sigma)
    np.testing.assert_array_equal(p.reduced_freqs, original_freqs)


def test_relaxed_scan_that_contains_ts_mirror_does_not_get_factor_two():
    p, _ = puckered_ts()
    evidence = hir_evidence(p)
    points = evidence['rotors'][0]['points']
    # A true measured mirror at 30 degrees is inside the methyl rotor's
    # 120-degree domain. Relaxation can couple the torsion to ring inversion.
    points[1]['geometry_angstrom'] = (p.geom * [-1., 1., 1.]).tolist()
    points[1]['electronic_energy_hartree'] = points[0]['electronic_energy_hartree']
    result = optical_counting(p, evidence)
    assert result['remaining_multiplier'] == 1.
    assert result['states_covered_by_hir'] == 2
    assert 'unrepresented_mirror' not in result


def relaxed_radical():
    data = json.loads((Path(__file__).parent / 'reference/relaxed_radical_hir.json').read_text())
    return saved_point(data)


def test_usable_partial_scan_retains_rigid_mirror_evidence_and_failed_points():
    p, _ = puckered_ts()
    evidence = hir_evidence(p)
    rotor = evidence['rotors'][0]
    for index in (1, 5, 9):
        rotor['points'][index]['status'] = 'failed'
        # A rejected geometry cannot change the adopted scan representation.
        rotor['points'][index]['geometry_angstrom'] = (p.geom * [-1., 1., 1.]).tolist()
    original = copy.deepcopy(evidence)
    result = optical_counting(p, evidence)
    assert result['remaining_multiplier'] == 2.
    witness = result['coordinate_coverage']
    assert witness['successful_observations'] == 9
    assert [point['point_index'] for point in witness['omitted_points']] == [1, 5, 9]
    assert witness['coupled_potential_accuracy_established'] is False
    assert evidence['rotors'][0]['points'] == original['rotors'][0]['points']


def test_ring_handedness_cannot_be_claimed_after_a_planar_scan_observation():
    p,_=puckered_ts()
    evidence=hir_evidence(p)
    geom=np.array(evidence['rotors'][0]['points'][5]['geometry_angstrom'])
    # Flattening the ring destroys every signed-volume obstruction. The
    # hypothetical point is not an observed full mirror. Only the explicitly
    # approximate energy rule may now choose a count; no handedness is proved.
    geom[:,2]=0.
    evidence['rotors'][0]['points'][5]['geometry_angstrom']=geom.tolist()
    result=optical_counting(p,evidence)
    from kinbot.optical import evaluate_optical
    assert evaluate_optical(p,rotors=evidence['rotors'])['status']=='unresolved'
    assert result['status']=='assumed' and result['remaining_multiplier']==2
    assert result['heuristic']=='harmonic_midpoint'
    assert result['states_covered_by_hir'] is None
