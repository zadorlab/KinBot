"""Keep distinct stereochemical channels on the selected lowest network route."""
import logging
import sys
from types import SimpleNamespace
from pathlib import Path
import numpy as np
from kinbot import pes
from kinbot.parameters import Parameters


def test_lowestpath_retains_classified_selected_edges_only():
    wells=['a','b','c','d']
    reactions=[['a','x',['b'],30.],['b','y',['a'],32.],['b','z',['c'],40.],
               ['a','off1',['d'],35.],['d','off2',['c'],45.]]
    conn=np.array([[0,1,0,1],[1,0,1,0],[0,1,0,1],[1,0,1,0]])
    paths={r[1]:'path-'+r[1] for r in reactions}
    args=({},wells,[],reactions,conn,None,{},'lowestpath',['a','c'])
    _,_,selected,_=pes.filter_stat_points(*args,stereopaths=paths)
    assert [r[1] for r in selected]==['x','z','y']
    # No path evidence: retain legacy single-route filtering.
    _,_,ordinary,_=pes.filter_stat_points(*args)
    assert [r[1] for r in ordinary]==['x','z']
    # A different configured product is not the same edge.
    assert all(r[1] not in ('off1','off2') for r in selected)


def test_no_kinbot_cli_forwards_lowestpath_and_species_names(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path("input.json").write_text('{"barrier_threshold":100}')
    par=Parameters("input.json",show_warnings=False).par
    par.update(smiles='[H][H]',charge=0,mult=1,verbose=0)
    monkeypatch.setattr(sys,'argv',['pes','input.json','no-kinbot','lowestpath','a','b'])
    monkeypatch.setattr(pes,'Parameters',lambda *a,**k:SimpleNamespace(par=par))
    monkeypatch.setattr(pes,'config_log',lambda *a,**k:logging.getLogger('test-pes'))
    monkeypatch.setattr(pes,'write_input',lambda *a,**k:None)
    monkeypatch.setattr(pes,'get_wells',lambda *a:None)
    monkeypatch.setattr(pes,'check_status',lambda *a:False)
    monkeypatch.setattr(pes.time,'sleep',lambda *a:None)
    def no_job(*a,**k):
        raise AssertionError('no-kinbot must not submit a job')
    monkeypatch.setattr(pes,'submit_job',no_job)
    seen=[]
    monkeypatch.setattr(pes,'postprocess',lambda par,jobs,task,names,mass:seen.append((task,names)))
    pes.main()
    assert seen==[('lowestpath',['a','b'])]


def test_lowestpath_keeps_ordinary_and_stereo_on_the_selected_edge():
    wells = ['a', 'b']
    conn = np.array([[0, 1], [1, 0]])
    for ordinary, stereo in [(30., 32.), (32., 30.)]:
        reactions = [['a', 'ordinary', ['b'], ordinary], ['a', 'stereo', ['b'], stereo]]
        paths = {'ordinary': 'ordinary', 'stereo': 'stereopath:relay'}
        _, _, selected, _ = pes.filter_stat_points({}, wells, [], reactions,
            conn, None, {}, 'lowestpath', ['a', 'b'], stereopaths=paths)
        assert {row[1] for row in selected} == {'ordinary', 'stereo'}
