"""Declared H2 loss uses product properties; calculated saddles take precedence."""

import copy
import itertools
import json
import logging
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from kinbot import constants, pes, postprocess
from kinbot.mess import MESS
from kinbot.parameters import Parameters
from kinbot.reaction_path import (compare_pathways, reaction_path_id,
                                  read_summary_paths, summary_path_line)
from kinbot.species_routing import routing_key
from tests import test_complex_mess as complex_output


@pytest.fixture
def output_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('me').mkdir()
    Path('input.json').write_text(json.dumps(dict(
        barrier_threshold=100., high_level=0, rotor_scan=0, multi_conf_tst=0,
        conformer_search=0, me=0, uq=0, epsilon=100., sigma=3.,
        queuing='local')))
    return Parameters('input.json', show_warnings=False).par


def declared_reaction(family):
    parent, reactions = complex_output.TestComplexMESS().reactions(True)
    reaction = reactions[0]
    reaction.instance_name = f'{routing_key(parent)}_{family}_2_3'
    parent.reac_type = [family]
    parent.reac_inst = [[1, 0, 2]]
    parent.reac_name = [reaction.instance_name]
    return parent, reaction


@pytest.mark.parametrize('family', ['hom_sci', 'h2_elim_barrierless'])
def test_declared_exit_summary_does_not_require_saddle_metadata(family):
    name = f'parent_{family}_2_3'
    reaction = SimpleNamespace(instance_name=name, ts=SimpleNamespace())
    assert reaction_path_id(reaction) is None
    assert summary_path_line(reaction) is None
    assert read_summary_paths([f'SUCCESS 10.00 {name} residue H2']) == {}
    assert compare_pathways(name, None, 10., 'actual_saddle', 'ordinary', 15.) == 'replace'
    assert compare_pathways('actual_saddle', 'ordinary', 15., name, None, 10.) == 'keep'


@pytest.mark.parametrize('family', ['hom_sci', 'h2_elim_barrierless'])
def test_pes_retains_actual_saddle_classes_instead_of_declared_exit(family):
    rows = [['R', f'parent_{family}_2_3', ['P', 'H2'], 0.],
            ['R', 'saddle_a', ['P', 'H2'], 20.],
            ['R', 'saddle_b', ['P', 'H2'], 25.]]
    paths = {'saddle_a': 'path_a', 'saddle_b': 'path_b'}
    for order in itertools.permutations(rows):
        selected = []
        for row in order:
            pes.select_summary_reaction(selected, row, paths)
        assert sorted(row[1] for row in selected) == ['saddle_a', 'saddle_b']


@pytest.mark.parametrize('family', ['hom_sci', 'h2_elim_barrierless'])
@pytest.mark.parametrize('mc', [0, 1])
def test_declared_exit_mess_uses_products_without_counting_a_placeholder_ts(
        output_directory, family, mc):
    parent, reaction = declared_reaction(family)
    writer = MESS(dict(output_directory, multi_conf_tst=mc), parent)
    original = writer._parent_symmetry

    def check_parent(point, **kwargs):
        assert point is not reaction.ts
        return original(point, **kwargs)

    with patch.object(writer, '_parent_symmetry', side_effect=check_parent), \
         patch.object(writer, 'write_barrier', side_effect=AssertionError('placeholder saddle')):
        writer.write_input(None)
    text = Path('me/mess_0000.inp').read_text()
    assert text.count('Core PhaseSpaceTheory') == 1
    assert 'Tunneling' not in text
    assert '{blessname}' not in text
    potential = {fields[0]: float(fields[1]) for line in text.splitlines()
                 if (fields := line.split()) and fields[0].startswith('Potential')}
    assert potential['PotentialPrefactor[au]'] == .6
    assert potential['PotentialPowerExponent'] == 6.


@pytest.mark.skipif(shutil.which('mess') is None, reason='Public MESS is not installed')
def test_declared_h2_exit_is_accepted_by_public_mess(output_directory):
    parent, _ = declared_reaction('h2_elim_barrierless')
    MESS(output_directory, parent).write_input(None)
    result = subprocess.run([shutil.which('mess'), 'mess_0000.inp'], cwd='me',
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert Path('me/mess_0000.out').stat().st_size > 0


@pytest.mark.parametrize('family', ['hom_sci', 'h2_elim_barrierless'])
def test_calculated_stereochemical_saddles_supersede_declared_exit_in_any_order(
        output_directory, tmp_path, monkeypatch, family):
    parent, saddles = complex_output.TestComplexMESS().reactions()
    for reaction in saddles:
        reaction.do_vdW = False
        reaction.ts.stereopath_id = reaction.instance_name
        reaction.ts.stereopath_metadata['id'] = reaction.instance_name
    declared = copy.copy(saddles[0])
    declared.instance_name = f'parent_{family}_2_3'
    declared.ts = copy.deepcopy(parent)
    for index, order in enumerate(itertools.permutations([declared, *saddles])):
        directory = tmp_path / str(index)
        directory.mkdir()
        (directory / 'me').mkdir()
        monkeypatch.chdir(directory)
        parent.reac_obj = list(order)
        parent.reac_ts_done = [-1] * len(order)
        parent.reac_type = [family if reaction is declared else 'test' for reaction in order]
        MESS(output_directory, parent).write_input(None)
        text = Path('me/mess_0000.inp').read_text()
        assert 'Union ! 2 stereochemical pathways' in text
        assert 'Core PhaseSpaceTheory' not in text


@pytest.mark.parametrize('family', ['h2_elim_barrierless'])
def test_declared_exit_l3_summary_references_the_l3_parent(output_directory, family):
    parent, reaction = declared_reaction(family)
    parent.l3energy = -.4 + parent.energy
    reaction.ts.l3energy = parent.l3energy
    reaction.products[0].l3energy = reaction.products[0].energy - .35
    reaction.products[1].l3energy = reaction.products[1].energy
    expected = (sum(product.l3energy + product.zpe for product in reaction.products)
                - parent.l3energy - parent.zpe) * constants.AUtoKCAL
    postprocess.create_summary_file(parent, SimpleNamespace(qc='fc'), output_directory)
    text = Path(f'summary_{routing_key(parent)}.out').read_text()
    success = next(line.split() for line in text.splitlines() if line.startswith('SUCCESS'))
    assert float(success[1]) == pytest.approx(expected, abs=.005)
    assert read_summary_paths(text.splitlines()) == {}


@pytest.mark.parametrize('family', ['h2_elim_barrierless'])
@pytest.mark.parametrize('l3', [0, 1])
def test_pes_saved_declared_exit_uses_product_threshold_at_the_selected_level(
        output_directory, monkeypatch, family, l3):
    base, reactant = '100', '200'
    products = ['300', '400']
    name = f'{reactant}_{family}_2_3'
    for job in (base, reactant):
        Path(job).mkdir()
        Path(f'{job}/summary_{job}.out').write_text('\n' * 5)
    Path(f'{reactant}/summary_{reactant}.out').write_text(
        '\n' * 5 + f'SUCCESS 6.28 {name} {" ".join(products)}\n')
    energies = {base: (-100., 0.), reactant: (-99.99, 0.),
                products[0]: (-99.98, 0.), products[1]: (0., 0.)}
    l3energies = {base: -110., reactant: -109.9,
                  products[0]: -109.85, products[1]: 0.}
    selected = l3energies if l3 else {key: value[0] for key, value in energies.items()}
    expected = (sum(selected[product] for product in products)
                - selected[base]) * constants.AUtoKCAL
    requested_jobs = []

    def energy_reader(jobs, job, ts, *args, **kwargs):
        requested_jobs.append((job, ts))
        assert not ts, 'Declared H2 loss must never request a saddle energy.'
        return energies[job]

    captured = {}
    monkeypatch.setattr(pes, 'logger', logging.getLogger('KinBot'), raising=False)
    monkeypatch.setattr(pes, 'ensure_current_run', lambda *args: None)
    monkeypatch.setattr(pes, 'validate_theory_restart', lambda *args: None)
    monkeypatch.setattr(pes, 'get_energy', energy_reader)
    monkeypatch.setattr(pes, 'get_l3energy', lambda job, par: (l3, l3energies[job]))
    monkeypatch.setattr(pes, 'copy_from_kinbot', lambda *args: None)
    monkeypatch.setattr(pes, 'create_rotdpy_inputs', lambda *args: None)

    def capture_output(par, wells, product_names, saddles, barrierless, *args):
        captured['barrierless'] = barrierless
        captured['saddles'] = saddles

    monkeypatch.setattr(pes, 'create_pesviewer_input', capture_output)
    par = dict(output_directory, L3_calc=l3, single_point_qc='orca', me=0)
    pes.postprocess(par, [base, reactant], 'l2threshold', [], 14.)
    assert len(captured['barrierless']) == 1
    assert captured['barrierless'][0][:3] == [reactant, name, products]
    assert captured['barrierless'][0][3] == pytest.approx(expected, abs=.01)
    assert captured['saddles'] == []
    assert all(job != name for job, _ in requested_jobs)
