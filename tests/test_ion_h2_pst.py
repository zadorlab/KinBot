"""Ion--H2 PST agrees with the classical Langevin capture limit.

Geometries carry analytic positive frequencies here; these are output and
statistical-mechanics checks, not calculated cation potential-energy surfaces.
"""
import copy
from dataclasses import replace
import json
import logging
import re
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase.build import molecule
from scipy.constants import atomic_mass, physical_constants

from kinbot import constants, pes, symmetry
from kinbot.conformer_counting import representative_record
from kinbot.mess import MESS, finalize_mc_mess
from kinbot.parameters import Parameters
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint


def ion_h2_system(charge):
    """Use NH2 cations so the ionic fragment also has a symmetry divisor."""
    mult = 1 if charge == 1 else 2
    points = []
    for formula, point_charge in (('CH4', charge), ('NH2', charge), ('H2', 0)):
        atoms = molecule(formula)
        if formula == 'CH4':
            atoms[0].symbol = 'N'
        if formula == 'H2':
            atoms.positions = [[0., 0., 0.], [0., 0., .740]]
        point = StationaryPoint(formula, point_charge, mult if point_charge else 1,
            atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        point.characterize()
        point.name = routing_name(point)
        symmetry.calculate_symmetry(point)
        point.freq = point.reduced_freqs = [1000.] * max(
            0, 3 * point.natom - (5 if point.natom == 2 else 6))
        point.energy, point.zpe = -1., 0.
        points.append(point)
    parent, ion, hydrogen = points
    ion.energy, hydrogen.energy = -.49, -.49
    name = routing_name(parent) + '_h2_elim_barrierless_1_2'
    reaction = SimpleNamespace(instance_name=name, products=[ion, hydrogen],
        prod_opt=[SimpleNamespace(species=p) for p in (ion, hydrogen)],
        ts=copy.deepcopy(parent), do_vdW=False, mp2=0)
    parent.reac_obj, parent.reac_ts_done = [reaction], [-1]
    parent.reac_type, parent.reac_inst, parent.reac_name = ['h2_elim_barrierless'], [[1, 0, 2]], [name]
    return parent, reaction


def parameters(directory, parent, mc, pes_mode):
    path = directory / 'input.json'
    path.write_text(json.dumps(dict(barrier_threshold=100., high_level=0,
        rotor_scan=0, multi_conf_tst=mc, conformer_search=mc, me=0, uq=0,
        epsilon=100., sigma=3., queuing='local', pes=pes_mode,
        charge=parent.charge, mult=parent.mult,
        structure=[value for symbol, xyz in zip(parent.atom, parent.geom)
                   for value in (str(symbol), *map(float, xyz))],
        TemperatureList=[300, 500, 1000], PressureList=[760])))
    return Parameters(str(path), show_warnings=False).par


def langevin_rate(products):
    # Atomic-unit k_L = 2*pi*|q|*sqrt(alpha/mu), converted to cm^3/s.
    ion, hydrogen = products
    mass_unit = physical_constants['electron mass'][0]
    reduced_mass = ion.mass * hydrogen.mass / (ion.mass + hydrogen.mass) * atomic_mass / mass_unit
    volume = physical_constants['Bohr radius'][0] ** 3
    time = physical_constants['atomic unit of time'][0]
    return 2. * np.pi * ion.charge * np.sqrt(5.437 / reduced_mass) * volume / time * 1.e6


def capture_rates(text):
    rows = text.split('Capture/Escape Rate Coefficients:', 1)[1].splitlines()
    rows = [row.split() for row in rows if row.split() and re.fullmatch(r'[\d.eE+\-]+', row.split()[0])]
    return np.array([float(row[-1]) for row in rows[:3]])


@pytest.mark.parametrize('charge', [1, 2])
@pytest.mark.parametrize('mc', [0, 1])
@pytest.mark.parametrize('pes_mode', [0, 1])
@pytest.mark.skipif(shutil.which('mess') is None, reason='Public MESS is not installed')
def test_public_mess_ion_h2_capture_matches_langevin(
        tmp_path, monkeypatch, charge, mc, pes_mode):
    parent, reaction = ion_h2_system(charge)
    assert reaction.products[0].sigma_ext == reaction.products[1].sigma_ext == 2
    par = parameters(tmp_path, parent, mc, pes_mode)
    worker = tmp_path / routing_name(parent) if pes_mode else tmp_path
    worker.mkdir(exist_ok=True)
    (worker / 'me').mkdir()
    monkeypatch.chdir(worker)
    MESS(par, parent).write_input(None)
    if pes_mode:
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(pes, 'logger', logging.getLogger('KinBot'), raising=False)
        product_names = sorted(routing_name(p) for p in reaction.products)
        product_key = '_'.join(product_names)
        parent_key = routing_name(parent)
        threshold = (sum(p.energy + p.zpe for p in reaction.products)
                     - parent.energy - parent.zpe) * constants.AUtoKCAL
        pes.create_mess_input(par, [parent_key], [product_key], [],
            [[parent_key, reaction.instance_name, product_names, threshold]], [],
            {parent_key: 0.}, {product_key: threshold},
            {parent_key: parent_key, product_key: parent_key}, parent.mass, False)
    text = Path('me/mess_0000.inp').read_text()
    assert text.count('AtomDistanceMin[angstrom] 0.5') == 1
    assert 'PotentialPowerExponent 4.' in text
    assert re.search(r'SymmetryFactor\s+4(?:\.0)?\s', text)
    result = subprocess.run([shutil.which('mess'), 'mess_0000.inp'], cwd='me',
        capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    rates = capture_rates(Path('me/mess_0000.out').read_text())
    assert len(rates) == 3
    np.testing.assert_allclose(rates, langevin_rate(reaction.products), rtol=.005)
    np.testing.assert_allclose(rates, rates[0], rtol=.005)


def test_ion_h2_minimum_distance_finalization_preserves_explicit_setting():
    note = '! Cation + neutral H2: isotropic ion-induced-dipole capture; alpha(H2)=5.437 a.u.\n'
    ordinary = 'TemperatureList[K] 300\nModel\nEnd\n'
    assert finalize_mc_mess(ordinary) == ordinary
    assert finalize_mc_mess(ordinary + note) == ordinary + note
    geometry = 'Geometry[angstrom] 2\n H 0. 0. 0.\n H 0. 0. .740\n'
    formatted = finalize_mc_mess(ordinary + geometry)
    assert formatted.count('AtomDistanceMin[angstrom] 0.5') == 1
    assert finalize_mc_mess(formatted) == formatted
    for explicit in ('AtomDistanceMin[angstrom] 0.6\n', 'AtomDistanceMin[au] 1.0\n'):
        assert finalize_mc_mess(explicit + ordinary + geometry) == explicit + ordinary + geometry


@pytest.mark.parametrize('geometry', [
    'Geometry[angstrom] 2\n H 0. 0. 0.\n O 0. 0. .970\n',
    'Geometry[angstrom] 3\n H 0. 0. 0.\n H 0. 0. .740\n O 0. 1. 0.\n'])
def test_other_geometries_keep_the_existing_mess_distance_limit(geometry):
    assert finalize_mc_mess(geometry) == geometry


@pytest.mark.parametrize('mc', [0, 1])
@pytest.mark.parametrize('pes_mode', [0, 1])
@pytest.mark.parametrize('declared_first', [False, True])
@pytest.mark.skipif(shutil.which('mess') is None, reason='Public MESS is not installed')
def test_public_mess_accepts_h2_with_an_ordinary_saddle(
        tmp_path, monkeypatch, mc, pes_mode, declared_first):
    from kinbot.stereo_identity import canonical_identity

    parent, reaction = ion_h2_system(1)
    par = parameters(tmp_path, parent, mc, pes_mode)
    reaction.instance_name = 'analytic_h2_elim'
    parent.reac_type, parent.reac_name = ['h2_elim'], [reaction.instance_name]
    reaction.ts = copy.deepcopy(parent)
    reaction.ts.name = reaction.instance_name
    reaction.ts.energy = parent.energy + 30./constants.AUtoKCAL
    reaction.ts.wellorts = 1
    reaction.ts.optical_reference = canonical_identity(parent)
    reaction.ts.stereopath_id = 'ordinary'
    reaction.ts.stereopath_metadata = {'schema': 'kinbot.stereopath.v1', 'id': 'ordinary'}
    reaction.ts.freq = reaction.ts.reduced_freqs = [-1000.] + [1000.]*8
    declared = copy.copy(reaction)
    declared.instance_name = 'analytic_h2_elim_barrierless'
    declared.ts = copy.deepcopy(parent)
    parent.reac_obj = [declared, reaction] if declared_first else [reaction, declared]
    parent.reac_type = ['h2_elim_barrierless' if point is declared else 'h2_elim'
                        for point in parent.reac_obj]
    parent.reac_name = [point.instance_name for point in parent.reac_obj]
    parent.reac_ts_done = [-1, -1]
    worker = tmp_path/routing_name(parent) if pes_mode else tmp_path
    worker.mkdir(exist_ok=True)
    (worker/'me').mkdir()
    monkeypatch.chdir(worker)
    MESS(par, parent).write_input(None)
    if pes_mode:
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(pes, 'logger', logging.getLogger('KinBot'), raising=False)
        product_names = sorted(routing_name(p) for p in reaction.products)
        product_key, parent_key = '_'.join(product_names), routing_name(parent)
        threshold = (sum(p.energy+p.zpe for p in reaction.products)
                     - parent.energy-parent.zpe)*constants.AUtoKCAL
        pes.create_mess_input(par, [parent_key], [product_key],
            [[parent_key, reaction.instance_name, product_names, 30.]], [], [],
            {parent_key: 0.}, {product_key: threshold},
            {parent_key: parent_key, product_key: parent_key}, parent.mass, False)
    text = Path('me/mess_0000.inp').read_text()
    assert text.count('AtomDistanceMin[angstrom] 0.5') == 1
    assert '! Cation + neutral H2:' not in text
    assert 'Core PhaseSpaceTheory' not in text
    result = subprocess.run([shutil.which('mess'), 'mess_0000.inp'], cwd='me',
        capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout+result.stderr


def test_charged_h2_does_not_receive_neutral_h2_polarizability(tmp_path, monkeypatch, caplog):
    parent, reaction = ion_h2_system(1)
    products = []
    for original, charge in zip(reaction.products, (0, 1)):
        point = StationaryPoint('product', charge, 2,
            atom=original.atom, geom=original.geom.copy())
        if point.natom == 2:
            point.geom[1, 2] = 1.06
        point.characterize()
        point.name = routing_name(point)
        symmetry.calculate_symmetry(point)
        point.freq = point.reduced_freqs = list(original.freq)
        point.energy, point.zpe = original.energy, original.zpe
        products.append(point)
    reaction.products = products
    reaction.prod_opt = [SimpleNamespace(species=p) for p in products]
    par = parameters(tmp_path, parent, 0, 0)
    monkeypatch.chdir(tmp_path)
    Path('me').mkdir()
    MESS(par, parent).write_input(None)
    text = Path('me/mess_0000.inp').read_text()
    assert 'PotentialPrefactor[au] 0.6' in text
    assert 'PotentialPowerExponent 6.' in text
    assert re.search(r'Core PhaseSpaceTheory[\s\S]*?SymmetryFactor\s+1(?:\.0)?\s', text)
    assert '! Cation + neutral H2:' not in text
    message = ('Charged H2 product: retaining generic -0.6/R6 PST; '
               'ion-induced-dipole capture requires the neutral residue polarizability.')
    assert '! WARNING: ' + message in text
    assert message in caplog.text


@pytest.mark.parametrize('charge', [1, 2])
def test_ion_h2_potential_is_independent_of_product_order(tmp_path, charge):
    parent, reaction = ion_h2_system(charge)
    writer = MESS(parameters(tmp_path, parent, 0, 0), parent)
    forward = writer._phase_space_models(reaction.products, 'N1H4', 0., 0., 1., 1.)
    reverse = writer._phase_space_models(reaction.products[::-1], 'N1H4', 0., 0., 1., 1.)
    for text in (forward, reverse):
        assert f'PotentialPrefactor[au] {2.7185 * charge**2:g}' in text
        assert 'PotentialPowerExponent 4.' in text
        assert re.search(r'SymmetryFactor\s+4(?:\.0)?\s', text)
        assert '! Cation + neutral H2:' in text


@pytest.mark.parametrize('pes_mode', [0, 1])
def test_ion_h2_mc_members_keep_their_absolute_rotational_optical_divisors(
        tmp_path, monkeypatch, pes_mode):
    from kinbot.mess import writer_members

    parent, reaction = ion_h2_system(1)
    ion, _ = reaction.products
    base = representative_record(ion)
    # Analytic member weights isolate the absolute product divisor from
    # the existing ratios between members and their selected parent.
    alternative = replace(base, index=2,
        zero_energy_hartree=base.zero_energy_hartree + 1. / constants.AUtoKCAL,
        sigma_ext=1, remaining_optical_weight=2.)

    def members(point, population, **kwargs):
        return {0: base, 2: alternative} if point is ion else writer_members(
            point, population, **kwargs)

    monkeypatch.setattr('kinbot.mess.writer_members', members)
    writer = MESS(parameters(tmp_path, parent, 1, pes_mode), parent)
    energy = '{ground_energy}' if pes_mode else 10.
    text = writer._phase_space_models(reaction.products, 'N1H4', energy, 0., 1., 1.)
    divisors = [float(value) for value in re.findall(r'SymmetryFactor\s+([-+\d.eE]+)', text)]
    # H2 has divisor 2: 2*(2/1)=4 and 2*(1/2)=1 for the two ion members.
    assert divisors == [4., 1.]
    assert text.count('PotentialPowerExponent 4.') == 2


def test_pst_cannot_silently_drop_a_third_fragment():
    parent, reaction = ion_h2_system(1)
    writer = MESS(dict(multi_conf_tst=0, pes=0), parent)
    with pytest.raises(ValueError, match='exactly two fragments; got 3'):
        writer._phase_space_models(reaction.products+[reaction.products[1]],
                                   'CH4', 0., 0., 1., 1.)
