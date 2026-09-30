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


def test_vrc_scan_overrides_agree_before_qc_and_pes_export(tmp_path, monkeypatch):
    import logging
    from types import SimpleNamespace
    import pytest
    from kinbot import pes
    from kinbot.vrc_tst_scan import VTS

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(pes, 'logger', logging.getLogger('KinBot'), raising=False)
    p = peroxy()
    mirror = copy.copy(p)
    mirror.geom = p.geom * [-1., 1., 1.]
    base, exact = str(p.chemid), routing_name(p)
    suffixes = ('_hom_sci_1_2', '_hom_sci_2_3')
    for option in ('vrc_tst_scan', 'vrc_tst_noscan'):
        for override in (None, [], [base + suffixes[1]]):
            settings = {base: [base + suffixes[0]]}
            if override is not None:
                settings[exact] = override
            par = dict(vrc_tst_scan={}, vrc_tst_noscan={})
            par[option] = settings
            for point in (p, mirror):
                name = routing_name(point)
                jobs = [name + suffix for suffix in suffixes]
                point.reac_obj = [SimpleNamespace(instance_name=job) for job in jobs]
                expected = ([] if override == [] else [jobs[1]]) if (
                    name == exact and override is not None) else [jobs[0]]
                qc = Mock()
                scanner = VTS(point, par, qc)
                for method in ('opt_products', 'save_products', 'find_scan_coos',
                               'find_equiv', 'do_scan', 'energies'):
                    setattr(scanner, method, Mock())
                scanner.calculate_correction_potentials()
                if expected:
                    scanner.opt_products.assert_called_once_with(expected)
                    scanner.do_scan.assert_called_once_with(expected, noscan=option == 'vrc_tst_noscan')
                else:
                    scanner.opt_products.assert_not_called()
                    scanner.do_scan.assert_not_called()
                assert qc.mock_calls == []
                for job in jobs:
                    # A selected export reaches its missing correction file;
                    # an excluded export must not try to read that file.
                    reaction = [name, job, ['fragment1', 'fragment2'], 0.]
                    if job in expected:
                        with pytest.raises(KeyError, match='Results of scan'):
                            pes.create_rotdpy_inputs(par, [reaction], [])
                    else:
                        pes.create_rotdpy_inputs(par, [reaction], [])
