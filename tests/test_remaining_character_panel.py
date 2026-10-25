"""Finite character stage: no inferred valence CI, publication or native execution."""
import json
from pathlib import Path

import numpy as np
import pytest

from kinbot import requested_orbital_recon as review
from kinbot.execution import atomic_json, file_sha256
from test_remaining_orbital_screens import tree, PRODUCTS  # common supplied geometry fixture
from test_remaining_orbital_screens import test_campaign_launch_and_resume_share_retained_production_capacity

PLAN=Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches/r7_remaining_character_20261025'


@pytest.fixture
def plan():
    return PLAN


def test_exact_panel_preflight_staging_and_repeat_preserve_live_work(tree, monkeypatch):
    production, output=tree
    before={p:p.read_bytes() for p in production.rglob('*') if p.is_file()}
    monkeypatch.setattr(review.subprocess,'run',lambda *a,**k:pytest.fail('No scheduler call during preparation'))
    result=review.batch(production, PLAN/'screens.json', output, shared_pool=True, preflight=True)
    assert result['status']=='preflight_passed' and result['tasks']==16
    assert not output.exists()
    assert before=={p:p.read_bytes() for p in production.rglob('*') if p.is_file()}
    result=review.batch(production, PLAN/'screens.json', output, shared_pool=True)
    assert (result['screens'],result['tasks'],result['max_jobs'])==(4,16,3)
    assert all(p.read_bytes()==v for p,v in before.items())
    for child in review._saved(output)['reviews']:
        child=Path(child); config=review.recon._load(child/'input.json')
        manifest=json.loads((child/'native/reconnaissance.json').read_text())
        assert len(manifest['jobs'])==4
        assert manifest['orbital_character_only']
        assert review.recon._coverage(manifest)=={}
        for job in manifest['jobs']:
            text=(Path(job['directory'])/(job['name']+'.inp')).read_text().lower()
            assert '{rhf' in text and '{multi' not in text and '{avas' not in text and 'rs2' not in text
            assert 'gprint,orbitals=8' in text and 'gthresh,thrprint=0.05' in text
            assert job['probe_active_space'] is None and not job.get('orbital_predecessor')
        if child.name in PRODUCTS:
            assert np.allclose(review.recon._anchor_geometry(config, config['anchors'][0]), PRODUCTS[child.name])
            assert config['anchors'][0]['name']=='exit_l2'
            assert config['anchors'][1]['distance_angstrom']==30.
    staged={p:p.read_bytes() for p in output.rglob('*') if p.is_file()}
    assert review.batch(production,PLAN/'screens.json',output,shared_pool=True,resume=True)['status']=='resume_ready'
    assert staged=={p:p.read_bytes() for p in output.rglob('*') if p.is_file()}


def test_resume_rejects_changed_recipe_without_mutation(tree,tmp_path):
    production,output=tree
    review.batch(production,PLAN/'screens.json',output,shared_pool=True)
    changed=tmp_path/'changed';changed.mkdir()
    for p in PLAN.glob('*.json'): (changed/p.name).write_bytes(p.read_bytes())
    recipe=json.loads((changed/'ch_h2ss_cs.json').read_text())
    recipe['reconnaissance']['multiplicity']=4
    atomic_json(changed/'ch_h2ss_cs.json',recipe)
    before={p:p.read_bytes() for p in production.rglob('*') if p.is_file()}
    with pytest.raises(ValueError,match='recipe changed'):
        review.batch(production,changed/'screens.json',output,shared_pool=True,resume=True)
    assert before=={p:p.read_bytes() for p in production.rglob('*') if p.is_file()}


def test_scf_character_report_is_hash_checked_read_only_and_no_cas_label(tree,monkeypatch):
    production,output=tree
    review.batch(production,PLAN/'screens.json',output,shared_pool=True)
    record=json.loads((output/'review_launch.json').read_text())
    manifest_path=Path(record['systems'][0]['manifest'])
    manifest=json.loads(manifest_path.read_text());job=manifest['jobs'][0]
    path=Path(job['directory'])/(job['name']+'.out')
    path.write_text('PROGRAM * RHF-SCF\n ELECTRON ORBITALS\n\n'
        ' 19.1   2.00000   -0.45000   2 1 pz 0.69000\n'
        ' 21.1   1.00000   -0.21000   1 2 px -0.71000\n'
        ' 22.1   0.00000    0.01000   3 1 py 0.82000\n\nEND\n')
    wfu=path.with_suffix('.wfu');wfu.write_bytes(b'checkpoint')
    xml=path.with_suffix('.xml');xml.write_text('<molpro/>')
    job.update(status='complete',energy_hartree=-835.,output_sha256=file_sha256(path),
               wfu=str(wfu),wfu_sha256=file_sha256(wfu),
               orbital_xml={'scf':dict(path=str(xml),sha256=file_sha256(xml))})
    atomic_json(manifest_path,manifest)
    before={p:p.read_bytes() for p in output.rglob('*') if p.is_file()}
    monkeypatch.setattr(review.subprocess,'run',lambda *a,**k:pytest.fail('Read-only report'))
    text=review.character(output)
    assert 'SCF Eh -835.' in text and 'SCF 21.1 occ=1.0 energy_Eh=-0.21' in text
    assert '(1, 2, \'px\', -0.71)' in text and 'XML scf' in text
    assert 'CASCI Eh' not in text and 'ACTIVE ' not in text
    assert before=={p:p.read_bytes() for p in output.rglob('*') if p.is_file()}
    xml.write_text('<changed/>')
    with pytest.raises(ValueError,match='XML changed'): review.character(output)


@pytest.mark.parametrize('key,value', [('orbital_print_virtuals',True),('orbital_print_virtuals',-1),
    ('orbital_print_virtuals',1.5),('orbital_print_threshold',float('nan')),('orbital_print_threshold',-1)])
def test_print_controls_fail_before_preparation(tree,tmp_path,key,value):
    production,output=tree
    review.batch(production,PLAN/'screens.json',output,shared_pool=True)
    source=Path(review._saved(output)['systems'][0]['configuration'])
    data=json.loads(source.read_text());data[key]=value
    invalid=tmp_path/'invalid.json';atomic_json(invalid,data)
    with pytest.raises(ValueError,match=key): review.recon._load(invalid)
