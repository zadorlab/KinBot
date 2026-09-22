import copy
import json
from unittest.mock import Mock
import numpy as np
from kinbot.parameters import Parameters
from kinbot.reaction_finder import ReactionFinder
from kinbot.species_routing import configured_selection, routing_name
from test_reaction_paths import peroxy


def test_configured_selector_precedes_legacy_for_two_enantiomers(tmp_path):
    p = peroxy()
    mirror = copy.copy(p)
    mirror.geom = p.geom * [-1., 1., 1.]
    legacy = str(p.chemid)
    names = [routing_name(p), routing_name(mirror)]
    assert names[0] != names[1]
    path = tmp_path / 'input.json'
    path.write_text(json.dumps({'barrier_threshold': 100.}))
    par = Parameters(path, show_warnings=False).par
    for mapping in ({legacy: [[0, 1]]}, {legacy: [[0, 1]], names[0]: [[1, 2]], names[1]: []}):
        for species in (p, mirror):
            expected = mapping.get(routing_name(species), mapping[legacy])
            params = dict(par, barrierless_saddle=mapping, homolytic_bonds=mapping)
            finder = ReactionFinder(species, params, None)
            finder.new_reaction = Mock()
            finder.search_hom_sci(species.natom, species.atom, species.bond, species.rads[0])
            assert finder.barrierless_saddle == expected
            assert finder.new_reaction.call_args.args[0] == expected
            assert configured_selection(mapping, species) == expected


def test_vrc_explicit_centres_accept_configured_keys():
    from kinbot.vrc_tst_scan import VTS
    p = peroxy()
    scanner = object.__new__(VTS)
    scanner.par = {'vrc_tst_scan_reac_cent': {routing_name(p): [p.atomid[0]]}}
    explicit = []
    scanner.explicit(p, -1, explicit, list(range(p.natom)))
    assert 0 in explicit
