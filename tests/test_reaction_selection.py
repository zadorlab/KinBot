"""General atom-motif reaction selection without connectivity IDs."""

from types import SimpleNamespace
import unittest

from kinbot.reaction_finder import ReactionFinder


class TestConfiguredReactionSelection(unittest.TestCase):
    def finder(self, selected, *, short_name='w1'):
        species = SimpleNamespace(
            name='generated-name', chemid=123, short_name=short_name)
        finder = object.__new__(ReactionFinder)
        finder.species = species
        finder.par = {'selected_reactions': selected}
        finder.reactions = {
            'hom_sci': [[0, 1], [1, 2]],
            'h2_elim': [[3, 0, 1, 6]],
            'other': [[0, 2]],
        }
        return finder

    def test_initial_selector_keeps_only_exact_requested_motifs(self):
        finder = self.finder({
            'initial': {
                'hom_sci': [[0, 1]],
                'h2_elim': [[3, 0, 1, 6]],
            }})
        finder._apply_configured_reaction_selection()
        self.assertEqual(finder.reactions, {
            'hom_sci': [[0, 1]],
            'h2_elim': [[3, 0, 1, 6]],
        })

    def test_explicit_species_key_precedes_initial(self):
        finder = self.finder({
            'initial': {'hom_sci': [[0, 1]]},
            '123': {'hom_sci': [[1, 2]]},
        })
        finder._apply_configured_reaction_selection()
        self.assertEqual(finder.reactions, {'hom_sci': [[1, 2]]})

    def test_unselected_product_is_unrestricted(self):
        finder = self.finder(
            {'initial': {'hom_sci': [[0, 1]]}}, short_name='p1')
        original = finder.reactions.copy()
        finder._apply_configured_reaction_selection()
        self.assertEqual(finder.reactions, original)

    def test_missing_requested_motif_fails_with_available_motifs(self):
        finder = self.finder({'initial': {'hom_sci': [[2, 3]]}})
        with self.assertRaisesRegex(ValueError, 'available motifs'):
            finder._apply_configured_reaction_selection()


if __name__ == '__main__':
    unittest.main()
