"""Final racemic networks: complete populations, distinct paths, no extra x2."""
import copy
import re
import unittest
from pathlib import Path

import pytest

from kinbot.mess import MESS, validate_mess_populations, finalize_mc_mess
from kinbot.mess_mirrors import complete_mirror_channels, population_keys
from kinbot.mess_racemates import _compatible_populations
from kinbot.species_routing import routing_name
from test_mess_mirror_completion import structural_transfer
import test_mess_mirror_completion as mirror_fixtures


def rrho(energy=0, symmetry=.5, frequency=500, reflected=False):
    x = -1 if reflected else 1
    return f'''    RRHO
      Geometry[angstrom] 4
 C 0 0 0
 H {x} 0 0
 F 0 1 0
 Cl 0 0 1
      Core RigidRotor
        SymmetryFactor {symmetry}
      End ! Core
      Frequencies[1/cm] 1
        {frequency}
      ElectronicLevels[1/cm] 1
        0 1
      ZeroEnergy[kcal/mol] {energy}
    End ! RRHO
'''


def well(name, own, mirror, energy=0, models=None):
    family = min(own, mirror)
    return (f'! kinbot_racemic_population {own} {family}\n Well {name}\n'
            f'! kinbot_population {own} {own}\n! kinbot_racemic_keys {own} {mirror}\n'
            + (rrho(energy=energy) if models is None else models) + ' End ! Well\n')


def barrier(name, left, right, keys, path=None, energy=30, models=None):
    a, ma, b, mb = keys
    text = (f' Barrier {name} {left} {right}\n'
            f'! kinbot_mirror_endpoints {a} {a} {b} {b}\n'
            f'! kinbot_racemic_route {a} {ma} {b} {mb}\n')
    if path:
        text += f'! kinbot_stereopath {path}\n'
    return text + (rrho(energy=energy) if models is None else models) + ' End ! Barrier\n'


def network(*blocks, reactant='r'):
    return 'Model\n Reactant ' + reactant + '\n' + ''.join(blocks) + 'End ! end kinetics\n'


def headers(text, kind):
    return re.findall(r'^\s*' + kind + r'\s+([^\n]+)', text, re.M)


def test_both_mirror_wells_and_ts_are_folded_without_an_extra_weight():
    text = network(well('r','R','S'), well('s','S','R'), well('a','A','A'),
                   barrier('tr','r','a',('R','S','A','A')),
                   barrier('ts','s','a',('S','R','A','A')), reactant='s')
    result = complete_mirror_channels(text)
    assert len(headers(result,'Well')) == 2
    assert len(headers(result,'Barrier')) == 1
    assert 'Reactant r\n' in result
    assert result.count('SymmetryFactor 0.5') == 3
    validate_mess_populations(result)
    assert complete_mirror_channels(result) == result
    # k = (2 q_TS)/(2 q_R) for a racemic entrance: no additional 2 or 1/2.
    assert (1/.5)/(1/.5) == 1


def test_achiral_entrance_to_two_mirrors_has_one_already_doubled_ts():
    result = complete_mirror_channels(network(well('a','A','A',models=rrho(symmetry=1)),
        well('r','R','S'), well('s','S','R'),
        barrier('ar','a','r',('A','A','R','S')),
        barrier('as','a','s',('A','A','S','R')), reactant='a'))
    assert len(headers(result,'Barrier')) == 1
    assert result.count('SymmetryFactor 0.5') == 2  # racemic product and TS
    assert (1/.5)/(1/1) == 2


def test_internal_enantiomerization_disappears():
    result = complete_mirror_channels(network(well('r','R','S'),well('s','S','R'),
        barrier('invert','r','s',('R','S','S','R'))))
    assert len(headers(result,'Well')) == 1
    assert not headers(result,'Barrier')


def test_two_site_paths_and_their_mirrors_give_two_routes_not_four():
    blocks=[well('r','R','S'),well('s','S','R'),well('a','A','A')]
    for p in ('h1','h2'):
        blocks += [barrier('r'+p,'r','a',('R','S','A','A'),path=p),
                   barrier('s'+p,'s','a',('S','R','A','A'),path=p)]
    result=complete_mirror_channels(network(*blocks))
    assert result.count('Union ! 2 stereochemical pathways') == 1
    assert result.count('original stereochemical path') == 2
    assert complete_mirror_channels(result) == result


def test_relative_configuration_survives_without_a_site_path_label():
    blocks=[well('r','R','S'),well('s','S','R'),well('b','B','b'),well('bb','b','B')]
    blocks += [barrier('rb','r','b',('R','S','B','b')),
               barrier('sbb','s','bb',('S','R','b','B')),
               barrier('rbb','r','bb',('R','S','b','B')),
               barrier('sb','s','b',('S','R','B','b'))]
    result=complete_mirror_channels(network(*blocks))
    assert result.count('Union ! 2 stereochemical pathways') == 1
    assert result.count('ZeroEnergy[kcal/mol] 30') == 2


def test_reverse_duplicate_retains_selected_eckart_order():
    result=complete_mirror_channels(network(well('r','R','S'),well('s','S','R'),well('a','A','A'),
        barrier('forward','r','a',('R','S','A','A'),energy=32),
        barrier('reverse','a','s',('A','A','S','R'),models=rrho(energy=30)+
                ' WellDepth[kcal/mol] 10\n WellDepth[kcal/mol] 20\n')))
    assert headers(result,'Barrier') == ['reverse a r']
    assert 'WellDepth[kcal/mol] 10\n WellDepth[kcal/mol] 20' in result


def tunneling_rrho(energy, depths):
    tunnel = ('      Tunneling Eckart\n'
              '        ImaginaryFrequency[1/cm] 1000\n'
              f'        WellDepth[kcal/mol] {depths[0]}\n'
              f'        WellDepth[kcal/mol] {depths[1]}\n'
              f'        CutoffEnergy[kcal/mol] {min(depths)}\n'
              '      End ! Tunneling\n')
    return rrho(energy=energy).replace('    End ! RRHO', tunnel + '    End ! RRHO')


@pytest.mark.parametrize('product_energy', [0., 5.])
@pytest.mark.parametrize('reverse', [False, True])
def test_folded_mirror_updates_eckart_references_without_changing_ts(product_energy, reverse):
    depths = [29., 30. - product_energy]
    endpoints, keys = ('s', 'a'), ('S', 'R', 'A', 'A')
    if reverse:
        endpoints, keys, depths = endpoints[::-1], ('A', 'A', 'S', 'R'), depths[::-1]
    model = tunneling_rrho(30., depths)
    original = network(well('r', 'R', 'S'), well('s', 'S', 'R', energy=1.),
        well('a', 'A', 'A', energy=product_energy),
        barrier('s_path', *endpoints, keys, models=model))
    result = finalize_mc_mess(complete_mirror_channels(original))
    expected = [30., 30. - product_energy]
    if reverse:
        expected.reverse()
    assert headers(result, 'Barrier') == ['s_path a r' if reverse else 's_path r a']
    assert list(map(float, re.findall(r'WellDepth\[kcal/mol\]\s+(\S+)', result))) == expected
    assert float(re.search(r'CutoffEnergy\[kcal/mol\]\s+(\S+)', result)[1]) == min(expected)
    # The complete TS calculation is retained; only tunneling references change.
    strip_depths = lambda s: re.sub(r'(?:WellDepth|CutoffEnergy)\[kcal/mol\]\s+\S+', '', s)
    rendered = re.search(r'Barrier s_path.*?(    RRHO.*?End ! RRHO)', result, re.S)[1]
    assert strip_depths(rendered).split() == strip_depths(model).split()
    assert complete_mirror_channels(result).split() == result.split()


@pytest.mark.parametrize('parent_path', [False, True])
def test_folded_mc_paths_use_selected_parent_references_not_ensemble_minima(parent_path):
    # R's selected parent is 2 despite its minimum 0; S's is 4 despite minimum 1.
    # The existing R-origin saddle supplies the selected-parent reference.
    models = lambda *parts: ' Union\n' + ''.join(parts) + ' End ! Union\n'
    blocks = [well('r', 'R', 'S', models=models(rrho(0), rrho(2))),
        well('s', 'S', 'R', models=models(rrho(1), rrho(4))), well('a', 'A', 'A', energy=5),
        barrier('r_path', 'r', 'a', ('R', 'S', 'A', 'A'), path='r_site',
                models=tunneling_rrho(40, [38, 35])),
        barrier('s_path', 's', 'a', ('S', 'R', 'A', 'A'), path='s_site',
                models=models(tunneling_rrho(30, [26, 25]), tunneling_rrho(31, [27, 26])))]
    # Without R's saddle, the saved parent is the only evidence of reference 2.
    # With it, preserve the attached Eckart convention even if the population
    # uses another level/reference, as direct inner-complex channels can do.
    saved_parent = 3 if parent_path else 2
    blocks[0] = blocks[0].replace(' Well r\n',
        f' Well r\n! kinbot_eckart_reference[kcal/mol] {saved_parent}\n')
    if not parent_path:
        blocks.pop(3)
    text = network(*blocks)
    result = finalize_mc_mess(complete_mirror_channels(text))
    assert ('Union ! 2 stereochemical pathways' in result) is parent_path
    assert list(map(float, re.findall(r'WellDepth\[kcal/mol\]\s+(\S+)', result))) == (
        ([38, 35] if parent_path else []) + [28, 25, 29, 26])
    assert list(map(float, re.findall(r'CutoffEnergy\[kcal/mol\]\s+(\S+)', result))) == (
        ([35] if parent_path else []) + [25, 26])
    assert complete_mirror_channels(result).split() == result.split()


@pytest.mark.parametrize('racemic', [False, True])
@pytest.mark.parametrize('retained_energy', [1., 2.])
def test_mirror_replacement_omits_nonpositive_tunneling_without_mc_finalization(racemic, retained_energy):
    text = network(well('r', 'R', 'S', energy=retained_energy),
        well('s', 'S', 'R'), well('a', 'A', 'A'),
        barrier('s_path', 's', 'a', ('S', 'R', 'A', 'A'), models=tunneling_rrho(1., [1., 1.])))
    if not racemic:
        text = re.sub(r'^! kinbot_racemic[^\n]*\n', '', text, flags=re.M)
        text = text.replace('kinbot_population R R', 'kinbot_population R S')
        text = text.replace('kinbot_population S S', 'kinbot_population S R')
        text = text.replace('kinbot_mirror_endpoints S S A A', 'kinbot_mirror_endpoints S R A A')
    result = complete_mirror_channels(text)  # harmonic/HIR output does not call the MC finalizer
    assert 'Tunneling Eckart' not in result
    assert 'submerged or has zero serialized tunneling depth' in result
    assert 'ZeroEnergy[kcal/mol] 1.0' in result  # TS energy is not raised


def test_folded_depth_that_rounds_to_zero_has_no_eckart_block():
    text = network(well('r', 'R', 'S', energy=.29), well('s', 'S', 'R'), well('a', 'A', 'A'),
        barrier('r_path', 'r', 'a', ('R', 'S', 'A', 'A'), path='r_site',
                models=tunneling_rrho(45.29, [45., 45.29])),
        barrier('s_path', 's', 'a', ('S', 'R', 'A', 'A'), path='s_site',
                models=tunneling_rrho(.29, [.29, .29])))
    result = complete_mirror_channels(text)
    # The inferred R reference is .28999999999999915; the raw difference is
    # positive, but the value actually written is zero.
    assert result.count('Tunneling Eckart') == 1
    assert 'submerged or has zero serialized tunneling depth' in result


def test_equivalent_mc_representations_compare_by_summed_weights():
    one=well('r','R','S',models=' Union\n'+rrho()+rrho(energy=2)+' End ! Union\n')
    two=well('s','S','R',models=' Union\n'+rrho(symmetry=1)+rrho(symmetry=1,reflected=True)+
             rrho(energy=2,symmetry=1)+rrho(energy=2,symmetry=1,reflected=True)+' End ! Union\n')
    assert _compatible_populations(one,two)
    result=complete_mirror_channels(network(one,two))
    assert result.count('End ! RRHO') == 2
    assert result.count('End ! Union') == 1


@pytest.mark.parametrize('models',[rrho(energy=.02),rrho(frequency=700),rrho()+rrho(energy=2)])
def test_incompatible_population_models_select_one_complete_model(models):
    result = complete_mirror_channels(network(well('r','R','S'),well('s','S','R',models=models)))
    assert len(headers(result, 'Well')) == 1
    assert result.count('End ! RRHO') == 1
    assert 'WARNING' in result
    assert complete_mirror_channels(result) == result


def test_conflicting_mc_ts_ensembles_are_not_silently_dropped():
    result = complete_mirror_channels(network(well('r','R','S'),well('s','S','R'),well('a','A','A'),
            barrier('tr','r','a',('R','S','A','A'),models=rrho(30)+rrho(32)),
            barrier('ts','s','a',('S','R','A','A'),models=rrho(30)+rrho(34))))
    assert 'Conflicting racemic TS conformer ensembles' in result
    assert 'ZeroEnergy[kcal/mol] 32' in result
    assert 'ZeroEnergy[kcal/mol] 34' not in result


def bimol(name, own, mirror, partner='H'):
    return (f' Bimolecular {name}\n! kinbot_population {own}_{partner} {own}_{partner}\n'
            f'! kinbot_racemic_keys {own}_{partner} {mirror}_{partner}\n'
            f'! kinbot_racemic_population {own} R\n Fragment f1\n'+rrho()+
            f' Fragment f2\n Atom\n Name {partner}\n ElectronicLevels[1/cm] 1\n 0 2\n End ! Atom\n'
            ' GroundEnergy[kcal/mol] -5\n End ! Bimolecular\n')


def test_full_fragment_composition_and_well_kind_are_part_of_grouping():
    result=complete_mirror_channels(network(well('r','R','S'),bimol('rh','R','S'),
        bimol('sh','S','R'),bimol('rf','R','S','F')))
    assert len(headers(result,'Well')) == 1
    assert len(headers(result,'Bimolecular')) == 2
    validate_mess_populations(result)


def test_specified_mode_is_unchanged():
    text=network(' Well r\n! kinbot_population R S\n'+rrho(symmetry=1),
                 ' Well a\n! kinbot_population A a\n'+rrho(symmetry=1),reactant='r')
    assert complete_mirror_channels(text) == text


def test_declared_namespace_survives_selection_of_a_mirror_geometry():
    from test_reaction_paths import peroxy
    from kinbot.stereo_identity import canonical_identity
    p=peroxy('C[C@H](F)Cl')
    p.optical_reference=canonical_identity(p)
    before=population_keys([p])
    p.geom *= [-1.,1.,1.]
    assert population_keys([p]) == before


class TestRacemicWriter(unittest.TestCase):
    setUp = mirror_fixtures.TestMirrorOutput.setUp

    def test_actual_fragment_writer_associates_racemic_models_in_both_orders(self):
        from kinbot.mess_racemates import _validate_family_models
        from kinbot.species_routing import routing_key
        from test_reaction_paths import peroxy

        for multi_conf in (0, 1):
            for pes in (0, 1):
                for atom_first in (False, True):
                    with self.subTest(multi_conf=multi_conf, pes=pes, atom_first=atom_first):
                        molecule, atom = peroxy(), peroxy('[H]')
                        molecule.name = routing_name(molecule)
                        atom.name = routing_name(atom)
                        if multi_conf:
                            molecule.conformer_index = [0]
                            molecule.conformer_geom = [molecule.geom.copy()]
                            molecule.conformer_zeroenergy = [molecule.energy + molecule.zpe]
                            molecule.conformer_freq = [molecule.freq]
                        products = [atom, molecule] if atom_first else [molecule, atom]
                        writer = MESS(dict(self.par, optical_population='racemic',
                                           multi_conf_tst=multi_conf, pes=pes), molecule)
                        writer.fragment_names = {routing_key(molecule): 'molecule',
                                                 routing_key(atom): 'hydrogen'}
                        key = '_'.join(sorted(routing_name(p) for p in products))
                        writer.bimolec_names = {key: 'products'}
                        output = writer.write_bimol(products, 0., 1., 1., 0, 0)
                        marker = re.search(r'^! kinbot_racemic_population .*$', output, re.M)
                        self.assertIsNotNone(marker)
                        following = output[marker.end():].lstrip()
                        self.assertTrue(following.startswith('Fragment '))
                        self.assertIn(molecule.name, following.splitlines()[0])
                        # Final PES assembly supplies the global energy before
                        # the racemic-network checks inspect this fragment file.
                        output = output.replace('{ground_energy}', '0.0')
                        _validate_family_models([output])
                        complete_mirror_channels(network(output))
                        # The last fragment must be checked too, not silently
                        # skipped because there is no following Fragment block.
                        changed = re.sub(r'(Frequencies\[1/cm\]\s+\d+\s+)\S+',
                                         r'\g<1>600.', output, count=1)
                        self.assertNotEqual(changed, output)
                        with self.assertRaisesRegex(ValueError, 'Conflicting intrinsic models'):
                            _validate_family_models([output, changed])

    def test_direct_and_pes_fold_the_same_real_writer_products(self):
        self.writer_roundtrip(0)

    def test_mc_direct_and_pes_keep_one_weighted_conformer_ensemble(self):
        self.writer_roundtrip(1)

    def writer_roundtrip(self, multi_conf):
        import shutil
        from kinbot.pes import create_mess_input
        from kinbot import symmetry
        p,q,ts,reaction=structural_transfer()
        qm=copy.deepcopy(q); qm.geom *= [-1.,1.,1.]
        if hasattr(qm,'optical_reference'):
            del qm.optical_reference
        qm.name=routing_name(qm)
        tm=copy.deepcopy(ts); tm.geom *= [-1.,1.,1.]; tm.name=ts.name+'_mirror'
        # A separately discovered mirror route has its own endpoint context.
        from kinbot.reaction_path import prepare_stereopath
        prepare_stereopath(tm, p, qm)
        symmetry.calculate_symmetry(qm)
        other=copy.copy(reaction); other.instance_name=tm.name; other.ts=tm
        other.products=[qm]; other.prod_opt=[type(reaction.prod_opt[0])(species=qm)]
        p.reac_obj=[reaction,other]; p.reac_ts_done=[-1,-1]
        p.reac_type=['intra_H_migration']*2; p.reac_inst=[None]*2
        p.reac_name=[ts.name,tm.name]
        Path('me').mkdir()
        par=dict(self.par,optical_population='racemic',pes=0,multi_conf_tst=multi_conf)
        if multi_conf:
            for point in (p,q,qm,ts,tm):
                point.conformer_index=[0]
                point.conformer_geom=[point.geom.copy()]
                point.conformer_zeroenergy=[point.energy+point.zpe]
                point.conformer_freq=[point.freq]
        MESS(par,p).write_input(None)
        direct=Path('me/mess_0000.inp').read_text()
        assert len(headers(direct,'Well')) == 2
        assert len(headers(direct,'Barrier')) == 1
        # Intermediate files and configured names survive final network folding.
        from kinbot.species_routing import mess_filename
        assert Path(mess_filename(q.name,0)).exists() and Path(mess_filename(qm.name,0)).exists()
        MESS(dict(par,pes=1),p).write_input(None)
        Path(p.name).mkdir()
        for file in Path('.').glob('*.mess'):
            shutil.copyfile(file,Path(p.name)/file.name)
        create_mess_input(dict(par,pes=1,smiles='[C](C)(CC)CCC',mult=2),
            [p.name,q.name,qm.name],[],
            [[p.name,ts.name,[q.name],30.],[p.name,tm.name,[qm.name],30.]],[],[],
            {p.name:0.,q.name:-5.,qm.name:-5.},{},
            {p.name:p.name,q.name:p.name,qm.name:p.name},p.mass,False)
        final=Path('me/mess_0000.inp').read_text()
        assert len(headers(final,'Well')) == 2
        assert len(headers(final,'Barrier')) == 1
        assert final.count('SymmetryFactor') == direct.count('SymmetryFactor')
        assert complete_mirror_channels(final) == final


def test_deferred_ts_offsets_are_resolved_before_ensemble_comparison():
    first=rrho(30)+rrho(30).replace('ZeroEnergy[kcal/mol] 30','ZeroEnergy[kcal/mol] 30 ! 2')
    second=rrho(30)+rrho(30).replace('ZeroEnergy[kcal/mol] 30','ZeroEnergy[kcal/mol] 30 ! 4')
    result = complete_mirror_channels(network(well('r','R','S'),well('s','S','R'),well('a','A','A'),
            barrier('tr','r','a',('R','S','A','A'),models=first),
            barrier('ts','s','a',('S','R','A','A'),models=second)))
    assert 'Conflicting racemic TS conformer ensembles' in result
    assert 'ZeroEnergy[kcal/mol] 32' in result
    assert 'ZeroEnergy[kcal/mol] 34' not in result


def test_hir_models_include_potential_symmetry_and_indexed_geometry():
    rotor=''' Rotor Hindered
 Group 2 3
 Axis 1 2
 Symmetry 2
 Potential[kcal/mol] 3
 0 1 1
 End ! Rotor
'''
    model=rrho().replace('ElectronicLevels',rotor+' ElectronicLevels')
    mirror=rrho(reflected=True).replace('ElectronicLevels',rotor+' ElectronicLevels')
    left=well('r','R','S',models=model)
    assert _compatible_populations(left,well('s','S','R',models=mirror))
    for other in (mirror.replace('0 1 1','0 2 2'),mirror.replace('Symmetry 2','Symmetry 3'),
                  mirror.replace('H -1 0 0','H -1.2 0 0')):
        assert not _compatible_populations(left,well('s','S','R',models=other))


def test_free_rotor_inner_geometry_is_orientation_independent():
    from kinbot.mess_racemates import _rrho_model
    rotor=''' Rotor Free
 Geometry[angstrom] 4
 C 0 0 0
 H 1 0 0
 F 0 1 0
 Cl 0 0 1
 ThermalPowerMax 50
 Group 2 3
 Axis 1 2
 Symmetry 2
 End ! Rotor
'''
    model=rrho().replace('ElectronicLevels',rotor+' ElectronicLevels')
    mirror=model.replace('H 1 0 0','H -1 0 0')
    assert _compatible_populations(well('r','R','S',models=model),
                                   well('s','S','R',models=mirror))


def test_reversed_mc_barrier_compares_unordered_eckart_depths_but_keeps_output():
    depths=' Tunneling Eckart\n WellDepth[kcal/mol] 10\n WellDepth[kcal/mol] 20\n End\n'
    forward=(rrho(30)+rrho(32)).replace('    End ! RRHO',depths+'    End ! RRHO')
    reverse=(rrho(30)+rrho(32)).replace('    End ! RRHO',
        depths.replace('10','TEMP').replace('20','10').replace('TEMP','20')+'    End ! RRHO')
    result=complete_mirror_channels(network(well('r','R','S'),well('s','S','R'),well('a','A','A'),
        barrier('tr','r','a',('R','S','A','A'),models=forward),
        barrier('ts','a','s',('A','A','S','R'),models=reverse)))
    assert len(headers(result,'Barrier')) == 1
    assert result.count(depths) == 2


def test_mirror_termolecular_dummy_sinks_are_combined():
    def sink(name,own,mirror):
        return (f' Bimolecular {name}\n! kinbot_population {own}_H_H {own}_H_H\n'
                f'! kinbot_racemic_keys {own}_H_H {mirror}_H_H\n Dummy\n')
    result=complete_mirror_channels(network(sink('rh2','R','S'),sink('sh2','S','R')))
    assert len(headers(result,'Bimolecular')) == 1
    assert result.count('Dummy') == 1


def test_mirror_hir_potential_can_run_in_the_opposite_direction():
    rotor=' Rotor Hindered\n Axis 1 2\n Group 2 3\n Symmetry 1\n Potential[kcal/mol] 4\n 0 1 3 2\n End ! Rotor\n'
    model=rrho().replace('ElectronicLevels',rotor+' ElectronicLevels')
    mirrored=model.replace('H 1 0 0','H -1 0 0').replace('0 1 3 2','0 2 3 1')
    assert _compatible_populations(well('r','R','S',models=model),well('s','S','R',models=mirrored))


def test_unclassified_parallel_barriers_keep_existing_guard():
    first=barrier('tr','r','a',('R','S','A','A'))
    first=re.sub(r'^! kinbot_racemic_route.*\n','',first,flags=re.M)
    with pytest.raises(ValueError,match='distinct classified'):
        complete_mirror_channels(network(well('r','R','S'),well('a','A','A'),
                                          first,first.replace('Barrier tr','Barrier tr2')))


def test_same_racemate_in_different_product_combinations_has_one_intrinsic_model():
    bad=bimol('sf','S','R','F').replace('        500','        700')
    fixed = complete_mirror_channels(network(bimol('rh','R','S'),bad))
    assert 'common-model approximation' in fixed
    assert '        700' not in fixed
    assert fixed.count('GroundEnergy[kcal/mol] -5') == 2
    assert complete_mirror_channels(fixed) == fixed
    # Global well and relative fragment energy zeros differ legitimately.
    result=complete_mirror_channels(network(well('r','R','S',energy=7),bimol('sh','S','R')))
    assert len(headers(result,'Well')) == len(headers(result,'Bimolecular')) == 1
    validate_mess_populations(result)


def test_common_hir_model_reflects_geometry_and_potential_together():
    from kinbot.mess_racemates import _harmonize_family_models
    rotor = (' Rotor Hindered\n Axis 1 2\n Group 2 3\n Symmetry 1\n'
             ' Potential[kcal/mol] 4\n 0 1 2 3\n End ! Rotor\n')
    source = bimol('rh', 'R', 'S').replace('ElectronicLevels', rotor + ' ElectronicLevels', 1)
    target = bimol('sf', 'S', 'R', 'F').replace('        500', '        700')
    output = _harmonize_family_models([source, target])[1]
    assert 'H -1.000000000000' in output
    assert '0 3 2 1' in output
    assert '        700' not in output
    assert 'GroundEnergy[kcal/mol] -5' in output
