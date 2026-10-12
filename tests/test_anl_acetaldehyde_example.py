"""The published small-system example remains one hands-off KinBot input."""

import json
from pathlib import Path

from kinbot.parameters import Parameters


def test_acetaldehyde_example_uses_automatic_cbh_and_portable_resources():
    path = (Path(__file__).parents[1] / 'examples' / 'anl' /
            'acetaldehyde_cbh_anl0_f12' / 'acetaldehyde.json')
    raw = json.loads(path.read_text())
    assert 'cbh_rung' not in raw
    assert 'fc_model_path' not in raw
    assert 'vrc_tst_max_nodes' not in raw
    assert 'rotdpy_max_jobs' not in raw
    parsed = Parameters(str(path), show_warnings=False)
    assert parsed.par['cbh_rung'] == 'auto'
    assert parsed.par['queue_job_limit'] == 8
    assert parsed.par['vrc_tst_max_nodes'] == 8
    assert parsed.par['rotdpy_max_jobs'] == 8
    assert parsed.theory_profiles['l1'].model_path == 'uma-s-1p2'
    assert parsed.par['composite_method'] == 'ANL0-F12'
    assert parsed.par['run_messpf'] == parsed.par['run_me'] == 1
