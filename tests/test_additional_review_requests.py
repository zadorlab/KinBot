from pathlib import Path
import json

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect

from kinbot import requested_orbital_recon as review, run_policy
from kinbot.execution import atomic_json
from test_state_support_panel import (
    support_tree, combined, spectrum, continuation_tree, resolution_tree,
    resolution_plan, tree, followup_plan, plan, original_tree, snapshot,
)

PLAN = Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches/r7_additional_requests_20261027'


def test_all_four_new_requests_native_decks_maps_and_immutable_source(support_tree, tmp_path, monkeypatch):
    production, _ = support_tree
    evidence = review._read(PLAN/'captured_geometry.json')
    worker = production/'ch3ss'/evidence['parent']; worker.mkdir(exist_ok=True)
    network = []
    for item in review._read(PLAN/'screens.json')['screens']:
        recipe = review._read(PLAN/item['recipe']); request = recipe['request_match']
        atomic_json(worker/'vrctst/orbital_review'/request['reaction']/'request.json',
                    dict(request, schema=1, status='needs_review'))
        if 'geometry_source' in recipe:
            network.append(dict(name=request['reaction'],
                                product_complex=recipe['geometry_source']['product_complex']))
    atomic_json(worker/'accepted_l2_network.json', dict(network=dict(reactions=network)))
    db = connect(worker/'kinbot.db')
    for row in evidence['endpoint_rows']:
        db.write(Atoms('CSSHHH', positions=row['positions']), name=row['name'], data=dict(status='normal'))
    pool = tmp_path/'new_pool'; pool.mkdir(); run_policy.update(pool, max_jobs=3)
    output = pool/'additional'; before = snapshot(production)
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **k: pytest.fail('Staging submits nothing'))
    opts = dict(shared_pool=True, resource_root=pool)
    assert review.batch(production, PLAN/'screens.json', output, preflight=True, **opts)['tasks'] == 16
    assert review.batch(production, PLAN/'screens.json', output, **opts)['tasks'] == 16
    assert review.batch(production, PLAN/'screens.json', output, resume=True, **opts)['tasks'] == 16
    for system in review._saved(output)['systems']:
        cfg = review.recon._load(system['configuration'])
        assert not cfg.get('_probe_space')
        assert cfg['multiplicity'] == (2 if system['title'] == 'stereo_h_exit_character' else 4)
        for job in review._read(system['manifest'])['jobs']:
            deck = (Path(job['directory'])/(job['name']+'.inp')).read_text()
            assert '{multi' not in deck and 'rs2' not in deck
        if 'exit' in system['title']:
            near = cfg['anchors'][0]; assert near['name'] == 'exit_l2'
            assert np.array_equal(review.recon._anchor_geometry(cfg, near), cfg['geometry_angstrom'])
            assert sorted(near['moving_atoms']) == ([5] if system['title'] == 'stereo_h_exit_character' else [2,3,6])
    assert snapshot(production) == before
