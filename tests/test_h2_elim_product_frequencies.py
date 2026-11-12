"""Declared H2-loss exits still require optimized products to be minima."""

import copy
from types import SimpleNamespace

import numpy as np
import pytest

from kinbot.reaction_generator import ReactionGenerator
from kinbot.stationary_pt import StationaryPoint


@pytest.mark.parametrize('family', ['hom_sci', 'h2_elim_barrierless'])
@pytest.mark.parametrize('lowest, expected', [(100., -1), (-20., -1), (-100., -999)])
def test_declared_exit_retains_product_frequency_check(
        tmp_path, monkeypatch, family, lowest, expected):
    monkeypatch.chdir(tmp_path)
    parent = StationaryPoint('ethyl cation', 1, 1, smiles='C[CH2+]')
    parent.characterize()
    heavy = StationaryPoint('vinyl cation', 1, 1, smiles='C=[CH+]')
    heavy.characterize()
    heavy.reduced_freqs = [lowest, 100.]
    hydrogen = StationaryPoint('hydrogen', 0, 1, atom=['H', 'H'],
        geom=np.array([[0., 0., 0.], [.74, 0., 0.]]))
    hydrogen.characterize()
    hydrogen.reduced_freqs = [4400.]
    # Product properties have already been selected. A declared exit uses the
    # parent well in place of a saddle; positive parent modes must be accepted.
    parent.reduced_freqs = [100.]
    reaction = SimpleNamespace(instance_name='example_' + family,
        species=parent, products=[heavy, hydrogen],
        prod_opt=[SimpleNamespace(species=heavy), SimpleNamespace(species=hydrogen)],
        ts=copy.copy(parent), do_vdW=False)
    parent.reac_obj = [reaction]
    parent.reac_inst = [[0, 1, 2]]
    parent.reac_type = [family]
    parent.reac_name = [reaction.instance_name]
    parent.reac_step = [0]
    parent.reac_ts_done = [5]
    par = dict(bimol=0, skip_reactions=[], high_level=0,
        barrier_threshold_L2=None, pes=0, imagfreq_threshold=50.,
        delete_intermediate_files=0)
    monkeypatch.setattr('kinbot.reaction_generator.time.sleep', lambda _: None)
    monkeypatch.setattr('kinbot.reaction_generator.postprocess.createPESViewerInput',
                        lambda *args: None)

    ReactionGenerator(parent, par, None, 'unused.json').generate()

    assert parent.reac_ts_done == [expected]
    if -50. <= lowest < 0.:
        assert heavy.reduced_freqs[0] == -lowest
