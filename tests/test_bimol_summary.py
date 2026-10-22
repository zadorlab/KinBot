"""Report successful searches and grouped saddles without hiding real failures."""
from types import SimpleNamespace

import pytest

from kinbot import bimol, postprocess
from kinbot.reaction_path import read_summary_paths
from kinbot.species_routing import routing_key
from tests.test_bimol_stereo_integration import donor_and_hydrogen, abstraction


@pytest.mark.parametrize('bimolecular', [False, True])
def test_summary_distinguishes_searches_duplicates_and_failures(tmp_path, monkeypatch, bimolecular):
    monkeypatch.chdir(tmp_path)
    species = donor_and_hydrogen()
    reaction = abstraction(species)
    bimol.build_bimolecular_ts(reaction, reaction.qc)
    bimol.update_ts_topology_from_irc(reaction)
    reaction.products = reaction.irc_prod.start_multi_molecular()[0]
    reaction.do_vdW = False
    reaction.search_parent = 'methanol_search'
    parent = SimpleNamespace(instance_name=reaction.search_parent, pose_search_template=True)
    duplicate = SimpleNamespace(instance_name='same_saddle', final_conformer_role='duplicate',
                                final_representative=reaction.instance_name)
    conformer = SimpleNamespace(instance_name='other_conformer', final_conformer_role='archived conformer',
                                final_representative=reaction.instance_name)
    failed_parent = SimpleNamespace(instance_name='failed_search', pose_search_template=True)
    failed = SimpleNamespace(instance_name='failed_refinement', search_parent='failed_search')
    exhausted = SimpleNamespace(instance_name='exhausted_search', pose_search_template=True)
    equivalent = SimpleNamespace(instance_name='equivalent_search', pose_search_template=True)
    reaction.search_members = [dict(reaction='equivalent_search')]
    species.reac_obj = [parent, reaction, duplicate, conformer, failed_parent, failed, exhausted, equivalent]
    species.reac_ts_done = [-998, -1, -998, -998, -998, -999, -998, -998]
    species.reac_inst = [reaction.instance] * len(species.reac_obj)
    species.reac_name = [r.instance_name for r in species.reac_obj]
    species.reac_type = ['abstraction'] * len(species.reac_obj)
    original_statuses = species.reac_ts_done.copy()
    postprocess.create_summary_file(species, reaction.qc, dict(bimol=bimolecular, high_level=0))
    text = (tmp_path / f'summary_{routing_key(species)}.out').read_text()
    failed_names = {line.split()[1] for line in text.splitlines() if line.startswith('FAILED')}
    assert failed_names == ({'failed_refinement', 'exhausted_search'} if bimolecular else
                           {'methanol_search', 'same_saddle', 'other_conformer',
                            'failed_search', 'failed_refinement', 'exhausted_search', 'equivalent_search'})
    assert (f'# DUPLICATE same_saddle; represented by {reaction.instance_name}' in text) is bimolecular
    assert (f'# CONFORMER other_conformer; represented by {reaction.instance_name}' in text) is bimolecular
    assert read_summary_paths(text.splitlines()) == {reaction.instance_name: 'ordinary'}
    assert len([line for line in text.splitlines() if line.startswith('SUCCESS')]) == 1
    assert species.reac_ts_done == original_statuses
