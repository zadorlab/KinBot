"""Bimolecular IRC reactant identity must respect the requested population."""
import numpy as np
import pytest

from kinbot.irc import IRC, BimolecularIRCResult
from kinbot.stereo_identity import canonical_identity
from test_bimol_stereo_integration import donor_and_hydrogen, abstraction


@pytest.mark.parametrize('population,reflected,accepted', [
    ('specified', False, True), ('specified', True, False), ('racemic', True, True)])
def test_irc_checks_the_requested_reactant_configuration(population, reflected, accepted):
    species = donor_and_hydrogen('CC[C@H](O)C')
    reaction = abstraction(species)
    reaction.qc.for_level = lambda level: reaction.qc
    reflect = np.array([-1., 1., 1.]) if reflected else np.ones(3)
    reactant_geom = species.geom * reflect
    product_geom = reaction.irc_prod.geom * reflect
    reaction.qc.get_qc_geom = lambda name, *args, **kwargs: (
        0, reactant_geom if '_IRC_F_' in name else product_geom)
    result = IRC(reaction, {'bimol': 1, 'optical_population': population}).irc2stationary_pt()
    assert isinstance(result, BimolecularIRCResult) == accepted
    if accepted:
        assert result.product_matches_expected


def test_irc_allows_a_reacting_chiral_centre_to_lose_its_configuration():
    species = donor_and_hydrogen('CC[C@H](O)C')
    reaction = abstraction(species, donor=2)
    reaction.qc.for_level = lambda level: reaction.qc
    assert canonical_identity(species)['is_chiral_configuration']
    assert not canonical_identity(reaction.irc_prod)['is_chiral_configuration']
    reaction.qc.get_qc_geom = lambda name, *args, **kwargs: (
        0, species.geom if '_IRC_F_' in name else reaction.irc_prod.geom)
    result = IRC(reaction, {'bimol': 1, 'optical_population': 'specified'}).irc2stationary_pt()
    assert isinstance(result, BimolecularIRCResult)
    assert result.product_matches_expected
