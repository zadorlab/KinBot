"""Identical reaction fragments share computation without losing stoichiometry."""

from types import SimpleNamespace

from kinbot.reaction_generator import ReactionGenerator


def test_identical_products_share_current_reaction_optimizer():
    product = SimpleNamespace(chemid=101)
    optimizer = SimpleNamespace(species=product)
    reaction = SimpleNamespace(products=[product, product], prod_opt=[])
    parent = SimpleNamespace(reac_obj=[reaction], reac_ts_done=[3])
    generator = ReactionGenerator(parent, {}, None, None)

    duplicate = SimpleNamespace(chemid=101)
    assert generator._existing_product_optimizer(
        0, duplicate, [optimizer]) is optimizer


def test_product_optimizer_is_reused_across_reactions():
    product = SimpleNamespace(chemid=101)
    optimizer = SimpleNamespace(species=product)
    previous = SimpleNamespace(products=[product], prod_opt=[optimizer])
    current = SimpleNamespace(products=[SimpleNamespace(chemid=101)],
                              prod_opt=[])
    parent = SimpleNamespace(reac_obj=[previous, current],
                             reac_ts_done=[4, 3])
    generator = ReactionGenerator(parent, {}, None, None)

    assert generator._existing_product_optimizer(
        1, current.products[0], []) is optimizer
