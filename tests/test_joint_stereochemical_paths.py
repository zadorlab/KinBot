"""Joint reacting-site fixtures; structures are not optimized TS calculations."""
import copy
from types import SimpleNamespace
import numpy as np
from kinbot.reaction_path import prepare_stereopath, path_geometry_allowed, summary_path_line
from test_reaction_paths import peroxy


def elimination(second_hydrogen):
    reactant = peroxy('CCCCC')
    reactant.mult = 1
    product = copy.deepcopy(reactant)
    product.bond[1, 8] = product.bond[8, 1] = 0
    product.bond[2, second_hydrogen] = product.bond[second_hydrogen, 2] = 0
    product.bond[8, second_hydrogen] = product.bond[second_hydrogen, 8] = 1
    product.bond[1, 2] = product.bond[2, 1] = 2
    product.bonds = [product.bond.copy()]
    product.calc_chemid(); product.find_atom_eqv()
    ts = copy.deepcopy(reactant)
    ts.wellorts = 1
    ts.bond = np.maximum(reactant.bond, product.bond)
    ts.bonds = [ts.bond.copy()]
    ts.calc_chemid(); ts.find_atom_eqv()
    prepare_stereopath(ts, reactant, product)
    return reactant, product, ts


def test_joint_relative_configuration_survives_path_metadata_and_json():
    keys = []
    for hydrogen in (10, 11):
        reactant, product, ts = elimination(hydrogen)
        assert ts.stereopath_metadata['site_relation'] == 'joint diastereotopic'
        assert path_geometry_allowed(ts, ts.geom)
        assert summary_path_line(SimpleNamespace(ts=ts, instance_name='joint'))
        keys.append(ts.stereopath_id)
        reverse = copy.deepcopy(ts)
        prepare_stereopath(reverse, product, reactant)
        assert reverse.stereopath_id == ts.stereopath_id
    assert keys[0] != keys[1]


def test_joint_path_identity_survives_atom_renumbering_and_global_reflection():
    reactant, product, ts = elimination(10)
    key = ts.stereopath_id
    order = np.arange(reactant.natom)[::-1]
    inverse = {int(old): new for new, old in enumerate(order)}
    for point in (reactant, product, ts):
        point.atom = np.asarray(point.atom)[order]
        point.geom = point.geom[order] * [-1., 1., 1.]
        point.bond = point.bond[np.ix_(order, order)]
        point.bonds = [bond[np.ix_(order, order)] for bond in point.bonds]
        point.rads = [np.asarray(rad)[order] for rad in point.rads]
        point.atom_eqv = [[inverse[i] for i in group] for group in point.atom_eqv]
    prepare_stereopath(ts, reactant, product)
    assert ts.stereopath_id == key


def test_ephemeral_changing_bond_mask_keeps_other_alkene_configuration():
    from kinbot.stereo_identity import canonical_identity
    # Ignore the first alkene only, as if its order changes in a reaction.
    first = peroxy('C/C=C/C/C=C/C')
    second = peroxy('C/C=C/C/C=C\\C')
    first.stereo_ignored_bonds = second.stereo_ignored_bonds = [(1, 2)]
    left, right = canonical_identity(first), canonical_identity(second)
    assert left['status'] == right['status'] == 'assigned'
    assert left['id'] != right['id']  # spectator E/Z remains explicit
    assert any('/' in text or '\\' in text for text in left['canonical_graphs'])


def test_joint_paths_reach_mess_once_each_with_duplicate_lowest_barrier_selection(tmp_path, monkeypatch):
    import json
    from pathlib import Path
    from kinbot import constants, symmetry
    from kinbot.mess import MESS
    from kinbot.parameters import Parameters
    from kinbot.species_routing import routing_name
    from test_mess_conformers import values
    monkeypatch.chdir(tmp_path)
    Path('me').mkdir()
    Path('input.json').write_text(json.dumps(dict(barrier_threshold=100., high_level=0,
        rotor_scan=0, multi_conf_tst=0, conformer_search=0, me=0, uq=0,
        epsilon=100., sigma=3., queuing='local', pes=0)))
    par = Parameters('input.json', show_warnings=False).par
    routes = []
    for hydrogen, barrier in [(10, 30.), (11, 32.), (10, 35.)]:
        parent, product, ts = elimination(hydrogen)
        parent.name = routing_name(parent)
        product.name = routing_name(product)
        product.energy = parent.energy - 5./constants.AUtoKCAL
        symmetry.calculate_symmetry(product)
        ts.name = parent.name + f'_h2_elim_9_{hydrogen+1}_{int(barrier)}'
        ts.energy = parent.energy + barrier/constants.AUtoKCAL
        ts.freq = ts.reduced_freqs = [-1000.] + [500.]*(3*ts.natom-7)
        symmetry.calculate_symmetry(ts)
        routes.append(SimpleNamespace(instance_name=ts.name, ts=ts, products=[product],
            prod_opt=[SimpleNamespace(species=product)], do_vdW=False, mp2=0))
    # Deliberately structural endpoint graphs, not optimized H2-elimination
    # thermochemistry. This fixture isolates final pathway selection/output.
    assert routing_name(routes[0].products[0]) == routing_name(routes[1].products[0])
    parent.reac_obj, parent.reac_ts_done = routes, [-1]*3
    parent.reac_type, parent.reac_inst = ['h2_elim']*3, [None]*3
    parent.reac_name = [reaction.instance_name for reaction in routes]
    writer = MESS(par, parent)
    # These deliberately unoptimized product coordinates do not define an
    # optical model. Keep this test about pathway selection and energies.
    monkeypatch.setattr(writer, '_parent_symmetry', lambda p, **kwargs: p.sigma_ext)
    writer.write_input(None)
    text = Path('me/mess_0000.inp').read_text()
    barrier = text[text.index('  Barrier'):]
    assert 'Union ! 2 stereochemical pathways' in barrier
    assert values(barrier, 'ZeroEnergy') == [30., 32.]
    assert routes[2].instance_name not in barrier
