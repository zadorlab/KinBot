"""The published small-system example remains one hands-off KinBot input."""

import json
from pathlib import Path

from kinbot.parameters import Parameters
from kinbot.reaction_finder import ReactionFinder
from kinbot.stationary_pt import StationaryPoint


def test_acetaldehyde_example_uses_automatic_cbh_and_portable_resources():
    path = (Path(__file__).parents[1] / 'examples' / 'anl' /
            'acetaldehyde_cbh_anl0_f12' / 'acetaldehyde.json')
    raw = json.loads(path.read_text())
    assert 'cbh_rung' not in raw
    assert 'fc_model_path' not in raw
    assert 'vrc_tst_max_nodes' not in raw
    assert 'rotdpy_max_jobs' not in raw
    assert 'reaction_search' not in raw
    assert 'families' not in raw
    assert 'selected_reactions' not in raw
    assert 'stop_product_reaction_search' not in raw
    parsed = Parameters(str(path), show_warnings=False)
    assert parsed.par['reaction_search'] == 1
    assert parsed.par['families'] == ['all']
    assert parsed.par['selected_reactions'] == {}
    assert parsed.par['stop_product_reaction_search'] == 0
    assert parsed.par['cbh_rung'] == 'auto'
    assert parsed.par['queue_job_limit'] == 8
    assert parsed.par['vrc_tst_max_nodes'] == 8
    assert parsed.par['rotdpy_max_jobs'] == 8
    assert parsed.theory_profiles['l1'].model_path == 'uma-s-1p2'
    assert parsed.par['composite_method'] == 'ANL0-F12'
    assert parsed.par['vrc_tst_scan_method'] == 'caspt2(2,2)'
    assert parsed.par['vrc_tst_scan_basis'] == 'aug-cc-pVTZ'
    assert parsed.par['run_messpf'] == parsed.par['run_me'] == 1


def test_acetaldehyde_example_discovers_validation_channels_without_selectors():
    path = (Path(__file__).parents[1] / 'examples' / 'anl' /
            'acetaldehyde_cbh_anl0_f12' / 'acetaldehyde.json')
    parameters = Parameters(str(path), show_warnings=False).par
    species = StationaryPoint('acetaldehyde', 0, 1, smiles='CC=O')
    species.characterize()
    finder = ReactionFinder(species, parameters, None)
    finder.find_reactions()

    assert [6, 1, 0] in finder.reactions['r12_insertion_R']
    assert [3, 0, 1, 6] in finder.reactions['h2_elim']
    assert [0, 1] in finder.reactions['hom_sci']
