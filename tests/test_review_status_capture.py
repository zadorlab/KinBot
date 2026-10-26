"""The handoff capture retains requests even without a readable endpoint DB."""
import importlib.util
import json
from pathlib import Path

from ase import Atoms
from ase.db import connect


def test_capture_reads_endpoints_and_failures_without_changing_or_creating_databases(tmp_path):
    script = Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches/r7_root_spectrum_20261026/capture_status.py'
    spec = importlib.util.spec_from_file_location('capture_status', script)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    for parent in ('123', '456'):
        w=tmp_path/'ch3ss'/parent
        request=w/'vrctst/orbital_review/exit/request.json'
        request.parent.mkdir(parents=True)
        request.write_text(json.dumps(dict(status='needs_review', reaction=parent+'_exit')))
        (w/'accepted_l2_network.json').write_text(json.dumps(dict(network=dict(reactions=[
            dict(name=parent+'_exit', product_complex='complex')]))))
    db=connect(str(tmp_path/'ch3ss/123/kinbot.db'))
    db.write(Atoms('SH',positions=[[0,0,0],[0,0,1.3]]),name='123_exit_IRC_R_prod_high',data={'status':'normal'})
    (tmp_path/'.kinbot_live_control.json').write_text('{"revision": 1}')
    log=tmp_path/'ch3ss/123/kinbot.log'; log.write_text('A'*140000+'END')
    bad=tmp_path/'ch3ss/123/vrctst/orbital_review/bad/request.json'
    bad.parent.mkdir(); bad.write_text('{partial')
    before={p:p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    result=module.capture(tmp_path)
    after={p:p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    assert before==after
    assert len(result['requests'])==2
    assert result['requests'][0]['geometries'][0]['symbols']==['S','H']
    assert result['requests'][0]['network'][0]['product_complex']=='complex'
    assert result['requests'][1]['geometries']==[]
    assert len(result['errors'])==2
    assert not (tmp_path/'ch3ss/456/kinbot.db').exists()
    tail=next(f for f in result['files'] if f['path'].endswith('kinbot.log'))
    assert tail['tail_bytes']==131072 and tail['text'].endswith('END') and not tail['changed']
    assert any(f['path']=='.kinbot_live_control.json' for f in result['files'])
