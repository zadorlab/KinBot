"""L3 RRHO/L2 rotor separation and vector-based anharmonic mode matching."""
from copy import deepcopy
from pathlib import Path
import json
import hashlib
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.io import read, write
from scipy.spatial.transform import Rotation

from kinbot import frequencies
from kinbot.anl import thermal
from kinbot.anl.results import parse_molpro_harmonic
from kinbot.anl.dispatch import prepare, advance, _geometry_hash
from kinbot.calculation import array_fingerprint
from kinbot.mess import MESS
from kinbot.stationary_pt import StationaryPoint
from tests.anl_fixture import dispatch_spec
from tests.test_anl_workflow import _complete_task, _complete_optimization
from tests.test_rotor_projection import rotor_fixture

FIX = Path(__file__).parent/'fixtures/anl_thermal'


def native(label='minimum'):
    text = (FIX/(label+'.molpro.txt')).read_text()
    symbols = ['S','H','H'] if label == 'minimum' else ['C','S','H','H','H']
    parsed = parse_molpro_harmonic(text, basis='cc-pVTZ',
        stationary_point='minimum' if label == 'minimum' else 'transition_state')
    coords, freq, modes = thermal.molpro_modes(text, symbols, parsed)
    return symbols, coords, freq, modes, parsed


def point(label='minimum'):
    symbols, coords, freq, modes, parsed = native(label)
    p = StationaryPoint(label, 0, 1 if label=='minimum' else 2,
        atom=symbols, geom=coords.copy())
    p.characterize(); p.wellorts = int(label != 'minimum')
    p.freq = freq.tolist(); p.reduced_freqs = freq.tolist(); p.dihed = []
    p.anl_thermal_model = dict(run='fixture', atom_map=list(range(p.natom)),
        l2_geometry_sha256=array_fingerprint(p.geom), geometry=None, vpt2=None,
        harmonic=dict(positions=coords.tolist(), frequencies=freq.tolist(),
            mass_weighted_modes=modes.tolist(), task_id='harmonic', method='CCSD(T)',
            basis='cc-pVTZ', native_output='fixture', source_sha256='a'*64))
    return p


@pytest.mark.parametrize('label', ['minimum','saddle'])
def test_native_spectral_hessian_reproduces_signed_modes(label):
    p = point(label)
    thermal.update_model(p)
    assert p.anl_thermochemistry_frequencies == pytest.approx(p.freq, abs=.001)
    assert sum(v<0 for v in p.anl_thermochemistry_frequencies) == p.wellorts
    assert not hasattr(p, 'rotor_projection') or not p.rotor_projection


def test_native_mode_tables_reject_missing_vectors_and_wrong_geometry():
    symbols, coords, f, q, parsed = native()
    t = (FIX/'minimum.molpro.txt').read_text()
    with pytest.raises(ValueError):
        thermal.molpro_modes(t.replace('SX1','SX2'), symbols, parsed)
    with pytest.raises(ValueError):
        thermal.molpro_modes(t.replace('2732.57','2832.57'), symbols, parsed)
    a = read(FIX/'minimum.l2.xyz')
    with pytest.raises(ValueError, match='accepted L2 geometry'):
        thermal.gaussian_modes((FIX/'minimum.gaussian_modes.txt').read_text(), symbols, a.positions*1.1)


def test_l2_shift_mapping_uses_vectors_despite_reversed_mode_table():
    p = point(); a = read(FIX/'minimum.l2.xyz')
    p.geom = a.positions
    p.anl_thermal_model['l2_geometry_sha256'] = array_fingerprint(p.geom)
    # Native Gaussian harmonic vectors; deliberately synthetic anharmonic shifts.
    f, q = thermal.gaussian_modes((FIX/'minimum.gaussian_modes.txt').read_text(), list(p.atom), a.positions)
    harmonic = f[::-1]; shifts = np.array([-60.,-40.,-20.])
    parsed = dict(harmonic_fundamentals_cm_inverse=harmonic,
        anharmonic_fundamentals_cm_inverse=harmonic+shifts,
        method='B2PLYP', basis='cc-pVTZ', dispersion='GD3BJ')
    out = thermal.attach_vpt2(p, parsed, FIX/'minimum.gaussian_modes.txt', a)
    assert out == pytest.approx(np.array(p.freq)+shifts[::-1], abs=.01)
    assert p.anl_thermochemistry_frequency_source['mode_map'] == [3,2,1]
    assert min(p.anl_thermochemistry_frequency_source['squared_overlaps']) > .98


def test_verified_vpt2_handoff_uses_overlap_and_keeps_evidence_on_reprojection(tmp_path):
    from kinbot.anl.workflow import attach_task_vpt2_frequencies
    p = point(); atoms = read(FIX/'minimum.l2.xyz')
    p.geom = atoms.positions
    p.anl_thermal_model['l2_geometry_sha256'] = array_fingerprint(p.geom)
    text = (FIX/'minimum.gaussian_modes.txt').read_text()
    f, _ = thermal.gaussian_modes(text, list(p.atom), atoms.positions)
    # Native vector excerpt plus deliberately synthetic, reversed VPT2 bands.
    bands = '\n'.join(f'{i}(1) active {v:.4f} {v-20:.4f} 1.0 1.0'
                      for i, v in enumerate(f[::-1], 1))
    output = ('#p B2PLYP/cc-pVTZ Freq=Anharmonic EmpiricalDispersion=GD3BJ\n'+text+
        '\nFundamental Bands\nMode(n) Status E(harm) E(anharm) I(harm) I(anharm)\n'+bands+
        '\nOvertones\nAnharmonic Zero Point Energy\n'+
        f'Harmonic       : cm-1 = {sum(f)/2:.5f} ;\n'+
        'Anharmonic Pot.: cm-1 = -30.00000 ;\nWatson+Coriolis: cm-1 = 0.00000 ;\n'+
        f'Total Anharm   : cm-1 = {sum(f)/2-30:.5f} ;\n'+
        'Normal termination of Gaussian 16\n')
    spec = dispatch_spec()
    spec['molecule'].update(symbols=list(p.atom), positions=p.geom.tolist())
    task = next(t for t in spec['tasks'] if t['id']=='gaussian_vpt2')
    geometry = next(t for t in spec['tasks'] if t['id']=='l2_geometry')
    task.update(depends_on=[]); spec['tasks']=[geometry,task]
    path=tmp_path/'vpt2.json'; path.write_text(json.dumps(spec))
    run=prepare(path,tmp_path/'vpt2'); _complete_optimization(run,geometry)
    advance(run,submit=False); _complete_task(run,task,output)
    attach_task_vpt2_frequencies(p,run,'gaussian_vpt2')
    assert p.anl_thermochemistry_frequency_source['mode_map']==[3,2,1]
    thermal.member_spectrum(p)
    thermal.update_model(p)
    proof=p.anl_thermochemistry_frequency_source
    assert proof['task_id']=='gaussian_vpt2'
    assert proof['source_sha256']==hashlib.sha256(output.encode()).hexdigest()
    assert p.anl_thermochemistry_frequencies==pytest.approx(np.array(p.freq)-20.,abs=.01)
    native_file=run/'tasks/gaussian_vpt2/vpt2.log'
    native_file.write_text(output+'changed')
    with pytest.raises(RuntimeError,match='hash|artifact|changed|differs'):
        attach_task_vpt2_frequencies(p,run,'gaussian_vpt2')


def test_ambiguous_modes_require_review_and_explicit_map_is_bounded():
    p = point(); q = np.array(p.anl_thermal_model['harmonic']['mass_weighted_modes'])
    mixed = q.copy(); mixed[1]=(q[1]+q[2])/np.sqrt(2); mixed[2]=(q[1]-q[2])/np.sqrt(2)
    p.anl_thermal_model['vpt2'] = dict(shifts=[-20,-40,-60], harmonic=p.freq,
        mass_weighted_modes=mixed.tolist(), mode_map=None, match_tolerance_cm_inverse=100.,
        method='synthetic', basis='fixture')
    with pytest.raises(ValueError, match='Ambiguous'):
        thermal.update_model(p)
    p.anl_thermal_model['vpt2']['mode_map'] = [1,2,3]
    thermal.update_model(p)
    assert p.anl_thermochemistry_frequencies == pytest.approx(np.array(p.freq)-[20,40,60],abs=.01)
    p.anl_thermal_model['vpt2']['harmonic'][0] += 500
    with pytest.raises(ValueError, match='mismatch exceeds'):
        thermal.update_model(p)


@pytest.mark.parametrize('saddle', [False,True])
def test_l3_projection_tracks_l2_rotor_demotion_and_preserves_l2(saddle):
    p,hess,eigen = rotor_fixture('methanol_co', saddle)
    p.freq,p.reduced_freqs,q,_ = frequencies.get_frequencies(p,hess,p.geom,return_modes=True)
    old = deepcopy(p.__dict__)
    p.anl_thermal_model = dict(run='fixture', atom_map=list(range(p.natom)),
        geometry=None, vpt2=None, l2_geometry_sha256=array_fingerprint(p.geom),
        harmonic=dict(positions=p.geom.tolist(), frequencies=(np.array(p.freq)*1.05).tolist(),
                      mass_weighted_modes=q.tolist()))
    thermal.update_model(p)
    assert p.anl_thermochemistry_frequencies == pytest.approx(np.array(p.reduced_freqs)*1.05)
    assert p.rotor_projection == old['rotor_projection']
    np.testing.assert_array_equal(p.geom,old['geom'])
    # The owning driver demotes the L2 rotor; its L3 mode returns to MESS.
    p.rotor_projection['rotors'][0]['projected'] = False
    thermal.update_model(p)
    assert p.anl_thermochemistry_frequencies == pytest.approx(np.array(p.freq)*1.05)


def completed_graph(tmp_path):
    symbols, xyz, f, q, parsed = native()
    spec = dispatch_spec(); spec['molecule'].update(symbols=symbols,positions=xyz.tolist())
    geom = next(t for t in spec['tasks'] if t['id']=='l3_geometry')
    geom['geometry_from']='initial'
    harmonic = next(t for t in spec['tasks'] if t['id']=='harmonic')
    spec['tasks']=[geom,harmonic]
    path=tmp_path/'spec.json'; path.write_text(json.dumps(spec))
    run=prepare(path,tmp_path/'run'); _complete_optimization(run,geom)
    advance(run,submit=False)
    _complete_task(run,harmonic,(FIX/'minimum.molpro.txt').read_text())
    state=json.loads((run/'state.json').read_text()); state['tasks']['harmonic']['status']='complete'
    (run/'state.json').write_text(json.dumps(state))
    return run


def test_verified_handoff_maps_rotated_permuted_atoms_and_detects_tampering(tmp_path):
    run=completed_graph(tmp_path)
    p=point(); p.__dict__.pop('anl_thermal_model')
    perm=[1,0,2]; rot=Rotation.from_rotvec([.2,.3,.4]).as_matrix()
    p.atom=np.array(p.atom)[perm]; p.geom=p.geom[perm]@rot+[3,4,5]
    p.characterize(); p.dihed=[]
    options=thermal.handoff_options(run)
    proof=thermal.attach_model(p,options)
    assert p.anl_thermochemistry_frequencies == pytest.approx([1204.06,2716.78,2732.57],abs=.01)
    np.testing.assert_allclose(p.anl_thermal_model['geometry']['positions'],p.geom,atol=1e-7)
    assert proof['atom_map'][1] == 0
    thermal.attach_model(p,options)  # repeated handoff regenerates the same model
    out=run/'tasks/harmonic/harmonic.out';out.write_text(out.read_text()+'tampered\n')
    with pytest.raises(RuntimeError, match='hash|artifact|changed|differs'):
        thermal.attach_model(p,options)


def test_handoff_selects_best_available_geometry_and_hessian_independently(tmp_path):
    run=completed_graph(tmp_path)
    s=json.loads((run/'state.json').read_text())
    s['tasks']['l3_geometry']['status']='failed'
    (run/'state.json').write_text(json.dumps(s))
    assert thermal.handoff_options(run)['geometry_task'] is None
    assert thermal.handoff_options(run)['harmonic_task']=='harmonic'


def test_mess_main_l3_geometry_keeps_explicit_l2_rotor_geometry():
    p,h,_=rotor_fixture('methanol_co')
    p.freq,p.reduced_freqs,q,_=frequencies.get_frequencies(p,h,p.geom,return_modes=True)
    p.sigma_ext=1; p.sigma_int=np.ones((p.natom,p.natom));p.nopt=1
    p.hir=SimpleNamespace(is_valid_rotor=lambda i:True, invalid_rotor_reason=lambda i:None, nrotation=6,
                         hir_energies=[[0,.002,.004,.002,.001,.002]])
    p.anl_thermal_model=dict(run='fixture', geometry={'positions':(p.geom*1.02).tolist(),'sigma_ext':1},
        harmonic=dict(positions=p.geom.tolist(), frequencies=p.freq,mass_weighted_modes=q.tolist()),
        vpt2=None, l2_geometry_sha256=array_fingerprint(p.geom))
    writer=MESS({'rotor_scan':1,'free_rotor_thrs':.1},p)
    main=writer.rotor_geom(p); rotor=writer.make_rotors(p,1.)
    assert writer.make_geom(p.geom*1.02,p.atom) in main
    assert writer.make_geom(p.geom,p.atom) in rotor
    assert 'Geometry[angstrom]' in rotor
    assert writer.make_geom(p.geom*1.02,p.atom) not in rotor


def test_selected_conformer_only_gets_l3_and_full_rrho_spectrum():
    p=point();p.sigma_ext=1
    p.anl_thermal_model['geometry']=dict(positions=(p.geom*1.01).tolist(),sigma_ext=2)
    w=MESS({},p)
    r=SimpleNamespace(geometry=p.geom,frequencies_cm1=[1,2,3],sigma_ext=1)
    g,f,sigma=w._member_data(p,r)
    assert f==pytest.approx(p.freq,abs=.01); assert sigma==2
    r.geometry=p.geom*1.2
    g,f,sigma=w._member_data(p,r)
    assert f==[1,2,3];assert sigma==1


def test_external_reference_gets_l3_geometry_and_spectrum(tmp_path):
    run=completed_graph(tmp_path);p=point()
    m=dict(name='hydrogen sulfide',symbols=list(p.atom),positions=p.geom.tolist(),
        charge=0,multiplicity=1,smiles='S',symmetry_factor=2.,frequencies_cm_inverse=[1000,2000,3000],
        hindered_rotor_blocks='',source={'provider':'kinbot-accepted-l2'})
    out=thermal.enrich_reference_model(m,run)
    assert out['frequencies_cm_inverse']==pytest.approx(p.freq,abs=.01)
    assert out['source']==m['source']
    assert out['thermal_source']['frequencies']['model']=='l3_harmonic'
    assert m['frequencies_cm_inverse']==[1000,2000,3000]


def test_complex_assessment_uses_accepted_energy_and_l3_harmonics():
    from kinbot.product_complex import reassess_product_complex
    from kinbot.energy import ZeroKEnergy
    a=point(); b=point(); c=point()
    for p,e in [(a,-398.),(b,-398.),(c,-796.01)]:
        p.energy=-100.; p.zpe=.01
        p.final_zero_k_energy=ZeroKEnergy(p.smiles,e,'matched-harmonic','fixture',0,1)
    # L2 would give an unbound complex and a spurious unstable mode.
    c.freq[0]=-80.
    r=SimpleNamespace(do_vdW=True,instance_name='outer-exit',
        irc_prod_opt=SimpleNamespace(species=c),prod_opt=[SimpleNamespace(species=p) for p in (a,b)])
    reassess_product_complex(r,dict(vdW_detection=0,imagfreq_threshold=50))
    assert r.do_vdW
    assert r.final_vdW_assessment['frequency_source']=='accepted_l3_harmonic'
    assert r.final_vdW_assessment['energy_source']=='anl_e0'
    assert r.vdW_depth==pytest.approx(.01*627.509474,rel=1e-5)


def test_vpt2_vectors_cannot_come_from_another_geometry():
    p = point(); a = read(FIX/'minimum.l2.xyz')
    p.geom = a.positions*1.1
    p.anl_thermal_model['l2_geometry_sha256'] = array_fingerprint(p.geom)
    f, _ = thermal.gaussian_modes((FIX/'minimum.gaussian_modes.txt').read_text(), list(p.atom), a.positions)
    parsed = dict(harmonic_fundamentals_cm_inverse=f,
        anharmonic_fundamentals_cm_inverse=f-20.,
        method='B2PLYP', basis='cc-pVTZ', dispersion='GD3BJ')
    with pytest.raises(ValueError, match='accepted L2 thermal model'):
        thermal.attach_vpt2(p, parsed, FIX/'minimum.gaussian_modes.txt', a)


def test_linear_molecule_keeps_one_stretch():
    symbols=['H','H'];g=np.array([[0.,0.,0.],[0.,0.,.74]])
    p=StationaryPoint('hydrogen',0,1,atom=symbols,geom=g);p.characterize()
    p.freq=p.reduced_freqs=[4400.]
    q=np.array([[0,0,-1,0,0,1]])/np.sqrt(2)
    p.anl_thermal_model=dict(run='fixture',geometry=None,vpt2=None,
        l2_geometry_sha256=array_fingerprint(g),harmonic=dict(positions=g.tolist(),
        frequencies=[4400.],mass_weighted_modes=q.tolist()))
    thermal.update_model(p)
    assert p.anl_thermochemistry_frequencies==pytest.approx([4400.])
