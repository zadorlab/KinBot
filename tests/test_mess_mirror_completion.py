import copy
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from kinbot import constants, symmetry
from kinbot.mess import MESS
from kinbot.mess_mirrors import complete_mirror_channels, population_keys
from kinbot.reaction_path import prepare_stereopath
from kinbot.stereo_identity import canonical_identity
from kinbot.species_routing import routing_name, routing_key
from test_reaction_paths import peroxy
import unittest
import test_stereopath_mess as fixtures


def structural_transfer():
    p = peroxy('[C](C)(CC)CCC')
    h = next(i for i in range(p.natom) if p.atom[i] == 'H' and p.bond[6,i])
    q = copy.deepcopy(p)
    q.bond[6,h] = q.bond[h,6] = 0
    q.bond[0,h] = q.bond[h,0] = 1
    q.bonds = [q.bond.copy()]
    axis = np.cross(p.geom[1]-p.geom[0], p.geom[2]-p.geom[0]); axis /= np.linalg.norm(axis)
    q.geom[h] = q.geom[0] + 1.1*axis
    q.calc_chemid(); q.find_atom_eqv(); symmetry.calculate_symmetry(q)
    ts = copy.deepcopy(p)
    ts.wellorts = 1; ts.geom = q.geom.copy(); ts.geom[h] = (q.geom[h]+p.geom[h])/2
    ts.bond = np.maximum(p.bond, q.bond); ts.bonds = [ts.bond.copy()]
    ts.optical_reference = canonical_identity(p)
    ts.calc_chemid(); ts.find_atom_eqv(); ts.find_cycle(); ts.find_conf_dihedral()
    prepare_stereopath(ts, p, q); symmetry.calculate_symmetry(ts)
    ts.freq = ts.reduced_freqs = [-1000.] + [500.]*(3*p.natom-7)
    p.name = routing_name(p); q.name = routing_name(q)
    ts.name = p.name + '_intra_H_migration_1_7'
    ts.energy = p.energy + 30/constants.AUtoKCAL
    q.energy = p.energy - 5/constants.AUtoKCAL
    reaction = SimpleNamespace(instance_name=ts.name, ts=ts, products=[q],
        prod_opt=[SimpleNamespace(species=q)], do_vdW=False, mp2=0)
    p.reac_obj=[reaction]; p.reac_ts_done=[-1]; p.reac_type=['intra_H_migration']
    p.reac_inst=[None]; p.reac_name=[ts.name]
    return p, q, ts, reaction


class TestMirrorOutput(unittest.TestCase):
    setUp = fixtures.TestStereopathMESS.setUp

    def test_real_writer_represents_both_configured_products_without_ts_factor_two(self):
        p, q, ts, reaction = structural_transfer()
        Path('me').mkdir()
        writer = MESS(dict(self.par, pes=0), p)
        writer.write_input(None)
        text = Path('me/mess_0000.inp').read_text()
        assert text.count('derived by global reflection from') == 2
        assert text.count('  Barrier ') == 2
        assert text.count('  Well ') == 3
        assert writer._parent_symmetry(ts) == ts.sigma_ext
        assert complete_mirror_channels(text) == text
        mirror_key = population_keys([q])[1]
        assert mirror_key in text
        assert text.index('Well well_mirror_') < text.index('  Barrier ')

    def test_explicit_mc_mirror_is_assigned_to_its_own_product_channel(self):
        p, q, ts, reaction = structural_transfer()
        ts.conformer_index = [0, 1]
        ts.conformer_geom = [ts.geom.copy(), ts.geom * [-1., 1., 1.]]
        ts.conformer_zeroenergy = [ts.energy + ts.zpe]*2
        ts.conformer_freq = [ts.freq]*2
        from kinbot.conformer_counting import writer_members
        records = writer_members(ts)
        assert len(records) == 1
        assert next(iter(records.values())).remaining_optical_weight == 1.

    def test_actual_pes_assembly_completes_the_same_mirror_channels(self):
        import shutil
        from kinbot.pes import create_mess_input
        p, q, ts, reaction = structural_transfer()
        Path('me').mkdir()
        par = dict(self.par, pes=1, smiles='[C](C)(CC)CCC', mult=2)
        writer = MESS(par, p)
        writer.write_input(None)
        Path(p.name).mkdir()
        for file in Path('.').glob('*.mess'):
            shutil.copyfile(file, Path(p.name)/file.name)
        create_mess_input(par, [p.name, q.name], [], [[p.name, ts.name, [q.name], 30.]], [], [],
                          {p.name: 0., q.name: -5.}, {}, {p.name: p.name, q.name: p.name}, p.mass, False)
        text = Path('me/mess_0000.inp').read_text()
        assert text.count('derived by global reflection from') == 2
        assert text.count('  Barrier ') == 2
        assert text.count('  Well ') == 3

    def test_termolecular_dummy_keeps_configured_population_annotation(self):
        p,q,ts,reaction=structural_transfer()
        fragments=[q,peroxy('[H]'),peroxy('[H]')]
        key='_'.join(sorted(routing_name(point) for point in fragments))
        writer=MESS(dict(self.par,pes=0),p)
        writer.termolec_names={key:'sink'}
        text=writer.write_termol(fragments,reaction,0)
        assert '! kinbot_population '+population_keys(fragments)[0] in text
        assert 'Dummy' in text
        assert 'ZeroEnergy' not in text
        assert 'Bimolecular sink' in text

        # PES assigns its own product name when combining several well searches.
        writer.par['pes'] = 1
        template = writer.write_termol(fragments, reaction, 0)
        from kinbot.species_routing import mess_filename
        assert Path(mess_filename(key, 0)).read_text() == template
        rendered = template.format(name='pr_5')
        assert rendered == text.replace('Bimolecular sink', 'Bimolecular pr_5')
        from kinbot.mess_networks import connected_models
        model = ('Reactant w_1\nModel\n  Well w_1\n Species\n End\n'
                 + rendered + '\n  Barrier rxn_7 w_1 pr_5\n RRHO\n End\n'
                 + 'End ! end kinetics\n')
        assert connected_models(model)[0]['products'] == ['pr_5']

    def test_same_ts_has_same_weight_from_either_discovery_direction(self):
        p, q, ts, reaction = structural_transfer()
        from kinbot.reaction_path import set_endpoint_populations
        writer = MESS(dict(self.par, pes=0), p)
        set_endpoint_populations(ts, [p], [q])
        forward = writer._parent_symmetry(ts)
        ts.optical_reference = canonical_identity(q)
        set_endpoint_populations(ts, [q], [p])
        assert writer._parent_symmetry(ts) == forward


def test_completion_works_after_pes_placeholders_and_is_not_a_racemate_default():
    text = '''Model\n  Well {a}\n! kinbot_population A A\n Species\n End\n  Well {r}\n! kinbot_population R S\n Species\n Geometry[angstrom] 1\n C 1.0 2.0 3.0\n End\n  Barrier {ts} {a} {r}\n! kinbot_mirror_endpoints A A R S\n RRHO\n ZeroEnergy[kcal/mol] {energy}\n End\nEnd ! end kinetics\n'''
    completed = complete_mirror_channels(text.format(a='w1', r='w2', ts='t1', energy=30))
    assert 'C -1.000000000000 2.0 3.0' in completed
    assert completed.count('ZeroEnergy[kcal/mol] 30') == 2
    assert complete_mirror_channels(completed) == completed
    # Two specified chiral endpoints do not authorize their absent mirrors.
    restricted = text.replace('A A', 'A a').format(a='w1', r='w2', ts='t1', energy=30)
    assert 'derived by global reflection' not in complete_mirror_channels(restricted)


def test_explicit_mirror_channel_is_not_added_twice():
    text = 'Model\n Well a\n! kinbot_population A A\n End\n Well r\n! kinbot_population R S\n End\n Well s\n! kinbot_population S R\n End\n Barrier tr a r\n! kinbot_mirror_endpoints A A R S\n End\n Barrier ts a s\n! kinbot_mirror_endpoints A A S R\n End\nEnd ! end kinetics\n'
    assert complete_mirror_channels(text) == text


def mirror_path_fixture():
    return """Model
 Well a
! kinbot_population A A
 ZeroEnergy[kcal/mol] 0
 End
 Well r
! kinbot_population R S
 ZeroEnergy[kcal/mol] -5
 End
 Well s
! kinbot_population S R
 ZeroEnergy[kcal/mol] -5
 End
! kinbot_stereopath path1
 Barrier ar1 a r
! kinbot_mirror_endpoints A A R S
 RRHO
 ZeroEnergy[kcal/mol] 30
 WellDepth[kcal/mol] 30
 WellDepth[kcal/mol] 35
 End
! kinbot_stereopath path2
 Barrier ar2 a r
! kinbot_mirror_endpoints A A R S
 Union
 RRHO
 ZeroEnergy[kcal/mol] 32
 End ! RRHO
 End ! Union
! kinbot_stereopath path1
 Barrier as1 a s
! kinbot_mirror_endpoints A A S R
 RRHO
 ZeroEnergy[kcal/mol] 30
 End
End ! end kinetics
"""


def test_partly_searched_mirror_side_keeps_every_path_once():
    text = complete_mirror_channels(mirror_path_fixture())
    assert text.count('Union ! 2 stereochemical pathways') == 2
    assert text.count('! kinbot_stereopath path1') == 2
    assert text.count('! kinbot_stereopath path2') == 2
    assert text.count('ZeroEnergy[kcal/mol] 32') == 2
    assert text.count('End ! Union') == 2  # nested MC content survives
    repeated = complete_mirror_channels(text)
    assert repeated.count('! kinbot_stereopath path2') == 2
    assert repeated.count('ZeroEnergy[kcal/mol] 32') == 2


def test_independent_mirror_energy_does_not_prevent_complete_model_reuse():
    text = mirror_path_fixture().replace('S R\n ZeroEnergy[kcal/mol] -5',
                                       'S R\n ZeroEnergy[kcal/mol] -4')
    result = complete_mirror_channels(text)
    assert 'S R\n ZeroEnergy[kcal/mol] -5' in result
    assert 'ZeroEnergy[kcal/mol] -4' not in result
    assert result.count('! kinbot_stereopath path2') == 2
    assert 'WellDepth[kcal/mol] 35' in result
    assert complete_mirror_channels(result).split() == result.split()


def test_whole_mc_mirror_model_is_reused_and_independent_paths_survive():
    from kinbot.mess_mirrors import _HEADER, split_model
    text = mirror_path_fixture()
    text = text.replace('R S\n ZeroEnergy[kcal/mol] -5',
        'R S\n Union\n RRHO\n Frequencies[1/cm] 1\n 500\n'
        ' ZeroEnergy[kcal/mol] -6\n End ! RRHO\n'
        ' RRHO\n Frequencies[1/cm] 1\n 600\n'
        ' ZeroEnergy[kcal/mol] -5\n End ! RRHO\n End ! Union')
    # An additional S-only route supplies the actual selected-parent reference
    # -4 through its depths. The R selected parent is -5, despite its MC minimum -6.
    text = text.replace('S R\n ZeroEnergy[kcal/mol] -5',
                        'S R\n ZeroEnergy[kcal/mol] -4')
    text = text.replace('RRHO\n ZeroEnergy[kcal/mol] 30\n WellDepth',
                        'RRHO\n ZeroEnergy[kcal/mol] 20\n End ! RRHO\n'
                        'RRHO\n ZeroEnergy[kcal/mol] 30\n WellDepth', 1)
    text = text.replace('End ! end kinetics', '''! kinbot_stereopath path3
 Barrier as3 a s
! kinbot_mirror_endpoints A A S R
 RRHO
 ZeroEnergy[kcal/mol] 40
 Tunneling Eckart
 CutoffEnergy[kcal/mol] 40
 WellDepth[kcal/mol] 40
 WellDepth[kcal/mol] 44
 End ! Eckart
 End ! RRHO
End ! end kinetics''')
    result = complete_mirror_channels(text)
    _, blocks, _ = split_model(result)
    mirrors = [b for b in blocks if _HEADER.search(b)[2].strip() in ('r', 's')]
    assert len(mirrors) == 2
    for block in mirrors:
        assert '500' in block and '600' in block
        assert 'ZeroEnergy[kcal/mol] -6' in block
        assert 'ZeroEnergy[kcal/mol] -5' in block
    for path in ('path1', 'path2', 'path3'):
        assert result.count('! kinbot_stereopath ' + path) == 2
    assert result.count('ZeroEnergy[kcal/mol] 40') == 2
    assert result.count('WellDepth[kcal/mol] 45.0') == 2
    assert 'WellDepth[kcal/mol] 46' not in result
    assert complete_mirror_channels(result).split() == result.split()


def test_reflection_reverses_torsion_direction_without_changing_samples():
    from kinbot.mess_mirrors import _reflect
    block = ' C 1.0 2.0 3.0\n Potential[kcal/mol] 4\n 0 1 4 2\n End\n'
    reflected = _reflect(block)
    assert 'C -1.000000000000 2.0 3.0' in reflected
    assert '0 2 4 1' in reflected
    assert '0 1 4 2' in _reflect(reflected)


def test_global_reflection_may_reverse_the_same_enantiomerization_channel():
    from kinbot.reaction_path import prepare_ts_context, path_geometry_allowed
    from kinbot.stereo_identity import optical_scope
    p = peroxy('C[C@H](F)Cl')
    q = copy.deepcopy(p); q.geom *= [-1., 1., 1.]
    ts = copy.deepcopy(p); ts.wellorts = 1
    ts.optical_reference = canonical_identity(p)
    prepare_ts_context(ts, p, q)
    assert path_geometry_allowed(ts, ts.geom)
    assert path_geometry_allowed(ts, ts.geom * [-1., 1., 1.])
    assert optical_scope(ts)['mirror_allowed']
    symmetry.calculate_symmetry(ts)
    assert MESS({'multi_conf_tst': 0}, p)._parent_symmetry(ts) == ts.sigma_ext / 2.
    from kinbot.conformer_counting import writer_members
    ts.conformer_index = [0, 1]
    ts.conformer_geom = [ts.geom.copy(), ts.geom * [-1., 1., 1.]]
    ts.conformer_zeroenergy = [ts.energy + ts.zpe]*2
    ts.conformer_freq = [ts.freq]*2
    records = writer_members(ts)
    assert len(records) == 2
    assert all(record.remaining_optical_weight == 1. for record in records.values())
