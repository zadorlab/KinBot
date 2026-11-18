"""Homolytic-scission thresholds use consistent units and energy references."""
import json
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import constants, pes, postprocess
from kinbot.parameters import Parameters
from kinbot.reaction_path import read_summary_paths
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


@pytest.mark.parametrize('family', ['hom_sci'])
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


@pytest.mark.parametrize('family', ['hom_sci'])
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
        assert not ts, 'Homolytic scission must never request a saddle energy.'
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
