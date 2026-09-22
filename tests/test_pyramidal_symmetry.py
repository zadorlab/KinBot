"""Localized pyramidal XY3 has three rotations, without optical doubling."""
from pathlib import Path
import numpy as np
import pytest
from kinbot import symmetry
from kinbot.mess import MESS
from kinbot.molecular_symmetry import geometric_mirror_states
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint


def point(center, arm, height, radius, charge=0, mult=1):
    angles=np.arange(3)*2*np.pi/3
    xyz=np.vstack(([0.,0.,height],np.column_stack((radius*np.cos(angles),radius*np.sin(angles),np.zeros(3)))))
    p=StationaryPoint(center+arm+'3',charge,mult,atom=[center]+[arm]*3,geom=xyz)
    p.characterize()
    p.energy,p.zpe=-100.,.1
    p.freq=p.reduced_freqs=[500.]*6
    return p


@pytest.mark.parametrize('center,arm,height,radius,charge,mult,expected',[
    ('N','H',.388,.94,0,1,3),('O','H',.3,.92,1,1,3),
    ('C','H',0.,1.08,0,2,6),('N','H',0.,1.,0,1,6),
    ('C','H',.01,1.08,0,2,6)])
def test_planar_and_pyramidal_graph_equivalent_arms(center,arm,height,radius,charge,mult,expected):
    p=point(center,arm,height,radius,charge,mult)
    assert len(symmetry.get_neighbors(p,0))==3
    symmetry.calculate_symmetry(p)
    assert p.sigma_ext==expected
    assert p.nopt==1
    assert geometric_mirror_states(p)==1
    assert np.all(np.array(p.sigma_int)==1)
    assert MESS({'multi_conf_tst':0},p)._parent_symmetry(p)==expected


def test_local_planarity_tolerance_and_coordinate_invariances():
    p=point('N','H',.04,1.)
    assert symmetry._threefold_atom_symmetry(p.geom,0,[1,2,3])==6
    p.geom[0,2]=.06
    assert symmetry._threefold_atom_symmetry(p.geom,0,[1,2,3])==3
    p.geom[0,2]=.388
    rotation=np.linalg.qr(np.array([[1.,2.,3.],[4.,2.,1.],[2.,4.,3.]]))[0]
    for xyz in (p.geom,p.geom@rotation+12.,p.geom*10.,p.geom*[-1.,1.,1.]):
        assert symmetry._threefold_atom_symmetry(xyz,0,[3,1,2])==3


def test_harmonic_and_mc_mess_use_three_without_optical_workaround(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    p=point('N','H',.388,.94)
    symmetry.calculate_symmetry(p)
    Path("input.json").write_text('{"barrier_threshold":100}')
    par=Parameters("input.json",show_warnings=False).par
    par.update(pes=0,rotor_scan=0,multi_conf_tst=0)
    writer=MESS(par,p);writer.create_short_names()
    text=writer.write_well(p,0.,1.,0)
    assert float(text.split('SymmetryFactor')[1].split()[0])==3.
    assert p.mess_optical_counting['remaining_multiplier']==1
    p.conformer_index=[0];p.conformer_geom=[p.geom.copy()]
    p.conformer_zeroenergy=[p.energy+p.zpe];p.conformer_freq=[p.freq]
    par['multi_conf_tst']=1
    writer=MESS(par,p);writer.create_short_names()
    text=writer.write_well(p,0.,1.,0)
    assert float(text.split('SymmetryFactor')[1].split()[0])==3.
    assert p.nopt==1


def test_planar_inversion_ts_keeps_six():
    p=point('N','H',0.,1.)
    p.wellorts=1;p.freq=p.reduced_freqs=[-1000.]+[500.]*5
    symmetry.calculate_symmetry(p)
    assert p.sigma_ext==6
    assert p.nopt==1


def test_long_bond_pyramid_and_planar_frame_without_new_element_support():
    # Geometry-only PH3/BF3-sized frames. This does not imply P/B support in
    # KinBot's existing atom characterization or QC workflow.
    for height,radius,expected in ((.65,1.26,3),(0.,1.31,6),(.01,1.31,6)):
        angles=np.arange(3)*2*np.pi/3
        xyz=np.vstack(([0.,0.,height],np.column_stack((radius*np.cos(angles),radius*np.sin(angles),np.zeros(3)))))
        assert symmetry._threefold_atom_symmetry(xyz,0,[1,2,3])==expected
