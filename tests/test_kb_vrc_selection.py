"""VRC auto-selection uses public ReactionFinder identities."""

from types import SimpleNamespace

from kinbot.kb import _accepted_homolytic_channels


def test_completed_hom_sci_with_internal_hs_label_is_selected():
    scission = SimpleNamespace(
        instance_name='parent_hom_sci_1_2', family_name='hs',
        products=[object(), object()])
    tight = SimpleNamespace(
        instance_name='parent_h2_elim_4_7', family_name='h2elim',
        products=[object(), object()])
    species = SimpleNamespace(
        reac_obj=[scission, tight], reac_ts_done=[-1, -1],
        reac_type=['hom_sci', 'h2_elim'])

    assert _accepted_homolytic_channels(species) == [
        'parent_hom_sci_1_2']


def test_incomplete_or_nonbimolecular_scission_is_not_selected():
    incomplete = SimpleNamespace(
        instance_name='incomplete', family_name='hs',
        products=[object(), object()])
    one_product = SimpleNamespace(
        instance_name='one_product', family_name='hs', products=[object()])
    species = SimpleNamespace(
        reac_obj=[incomplete, one_product], reac_ts_done=[0, -1],
        reac_type=['hom_sci', 'hom_sci'])

    assert _accepted_homolytic_channels(species) == []
