"""C--H--C connections, counterexamples and downstream chemical operations.

The reference file distinguishes published Cartesian coordinates, our saved
ethyl minimum, and constructed geometries for ring/cage/complex tests.
"""
import copy
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pytest
from ase import Atoms

from kinbot import constants, symmetry
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.reaction_finder import ReactionFinder
from kinbot.reaction_path import endpoint_snapshot
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity
from kinbot.stereo_routing import guard_well_job, StereoRoutingError, _chemical_context


REFERENCES = json.loads((Path(__file__).parent/'reference/hydrogen_bridges.json').read_text())


def point(name, order=None, **overrides):
    data = dict(REFERENCES[name], **overrides)
    atom, geom = np.array(data['atom']), np.array(data['geom'])
    if order is not None:
        atom, geom = atom[order], geom[order]
    result = StationaryPoint(name, data['charge'], data['mult'], atom=atom, geom=geom,
                             wellorts=data.get('wellorts', 0))
    result.characterize()
    return result


@pytest.mark.parametrize('name', list(REFERENCES))
def test_reference_structures_have_consistent_connectivity_and_valence(name):
    species = point(name)
    assert len(species.hydrogen_bridges) == REFERENCES[name]['expected_bridges']
    assert len(species.bonds) == len(species.rads)
    assert nx.is_connected(nx.from_numpy_array(species.bond > 0))
    assert canonical_identity(species)['status'] == 'assigned'
    for bond, rad in zip(species.bonds, species.rads):
        assert min(rad) >= 0
        assert all(np.sum(bond[i]) <= 1 for i, symbol in enumerate(species.atom) if symbol == 'H')
        assert all(np.sum(bond[i]) <= 4 for i, symbol in enumerate(species.atom) if symbol == 'C')
    for first, hydrogen, second in species.hydrogen_bridges:
        assert species.bond[first, hydrogen] == species.bond[hydrogen, second] == 1
        assert any(matrix[first, hydrogen] for matrix in species.bonds)
        assert any(matrix[second, hydrogen] for matrix in species.bonds)
        assert all(hydrogen not in ring for ring in species.cycle_chain)
    fragments, _ = species.start_multi_molecular()
    assert len(fragments) == 1
    assert fragments[0].charge == species.charge
    assert fragments[0].mult == species.mult


@pytest.mark.parametrize('name', ['ethyl', 'vinyl', 'protonated_ethane', '2butyl', 'cyclooctyl', 'cage', 'radical_complex'])
def test_bridge_identity_is_independent_of_atom_order(name):
    original = point(name)
    expected = canonical_identity(original)
    rng = np.random.default_rng(20261113)
    for _ in range(6):
        permutation = rng.permutation(original.natom)
        observed = point(name, permutation)
        assert canonical_identity(observed)['id'] == expected['id']
        assert observed.chemid == original.chemid


@pytest.mark.parametrize('name', ['ethyl', 'vinyl', 'protonated_ethane', '2butyl'])
def test_symmetric_bridged_ions_have_two_external_rotations(name):
    species = point(name)
    symmetry.calculate_symmetry(species)
    assert species.sigma_ext == 2
    assert canonical_identity(species)['is_chiral_configuration'] is False
    for torsion in species.dihed:
        assert len(set(torsion)) == 4


def test_methanium_is_outside_the_carbon_hydrogen_carbon_bridge_model():
    from kinbot.stereo_identity import require_supported_identity, UnsupportedStereochemistry

    # ATcT, Methanium, B2PLYP-D3/aug-cc-pVTZ
    # https://atct.anl.gov/Thermochemical%20Data/version%201.124/xyz/CH5+.xyz
    geom = [[-.075345, -.007652, 0.], [-.184031, 1.093287, 0.],
            [-.463537, -.393953, -.937377], [-.463537, -.393953, .937377],
            [1.049953, .375591, 0.], [.958271, -.589856, 0.]]
    species = StationaryPoint('methanium', 1, 1, atom=['C'] + ['H']*5,
                              geom=np.array(geom))
    species.characterize()
    assert species.hydrogen_bridges == []
    with pytest.raises(UnsupportedStereochemistry, match='coordination'):
        require_supported_identity(species)


def test_ethyl_bridge_prevents_false_carbon_carbon_rotation():
    species = point('ethyl')
    assert species.cycle_chain == []
    assert species.dihed == []


@pytest.mark.parametrize('substitution, expected', [
    ('none', [3, 3]), ('one_D', [1, 3]), ('three_D', [3, 3]),
    ('bridge_D', [3, 3]), ('one_13C', [3, 3])])
def test_protonated_ethane_retains_two_methyl_rotors_and_isotope_divisors(substitution, expected):
    species = point('protonated_ethane')
    first, hydrogen, second = species.hydrogen_bridges[0]
    species.isotopes = [0]*species.natom
    terminal = [i for i, symbol in enumerate(species.atom)
                if symbol == 'H' and species.bond[first, i] and i != hydrogen]
    if substitution == 'one_D':
        species.isotopes[terminal[0]] = 2
    elif substitution == 'three_D':
        for atom in terminal:
            species.isotopes[atom] = 2
    elif substitution == 'bridge_D':
        species.isotopes[hydrogen] = 2
    elif substitution == 'one_13C':
        species.isotopes[first] = 13
    symmetry.calculate_symmetry(species)
    assert len(species.dihed) == 2
    assert sorted(species.sigma_int[b][c] for a, b, c, d in species.dihed) == expected


def test_protonated_ethane_bridge_is_not_the_separated_methane_methyl_identity():
    connected = point('protonated_ethane')
    separated = copy.deepcopy(connected)
    matrix = connected.bonds[0]
    components = list(nx.connected_components(nx.from_numpy_array(matrix > 0)))
    assert len(components) == 2
    separated.geom[list(components[0])] += [10., 0., 0.]
    separated.reset_order()
    assert separated.hydrogen_bridges == []
    assert canonical_identity(separated)['id'] != canonical_identity(connected)['id']


@pytest.mark.parametrize('name', ['cyclooctyl', 'cage'])
def test_transannular_bridge_preserves_real_carbon_rings(name):
    species = point(name)
    assert species.cycle_chain
    assert all(all(species.atom[i] == 'C' for i in ring) for ring in species.cycle_chain)
    assert any(not species.bond[first, second] for first, hydrogen, second in species.hydrogen_bridges)


@pytest.mark.parametrize('charge,mult,wellorts', [(0, 2, 0), (-1, 1, 0), (1, 1, 1)])
def test_neutral_anionic_and_transition_state_contacts_keep_ordinary_treatment(charge, mult, wellorts):
    species = point('ethyl', charge=charge, mult=mult, wellorts=wellorts)
    assert species.hydrogen_bridges == []


def test_short_secondary_contact_does_not_create_a_valence_impossible_bridge():
    species = StationaryPoint('distorted ethane cation', 1, 2, atom=['C','C']+['H']*6,
        geom=np.array([[0,0,0],[1.5,0,0],[.75,1,0],[-.6,-.55,.9],[-.6,-.55,-.9],
                       [2.3,-.6,0],[1.8,.6,.85],[1.8,.6,-.85]]))
    species.characterize()
    assert species.hydrogen_bridges == []


@pytest.mark.parametrize('first, second', [(1.09, 1.8), (1.072, 1.980), (1.09, 1.35), (1.20, 1.49), (1.6, 1.6)])
def test_normal_donor_or_long_unbalanced_contacts_do_not_form_a_bridge(first, second):
    species = point('cage')
    carbon1, hydrogen, carbon2 = species.hydrogen_bridges[0]
    axis = species.geom[carbon2]-species.geom[carbon1]
    axis /= np.linalg.norm(axis)
    middle = .5*(species.geom[carbon1]+species.geom[carbon2])
    # Shift each connected carbon half of the cage while retaining its arms.
    # This is a geometric counterexample, not an optimized molecule.
    species.geom[carbon1] = middle - first*axis
    species.geom[carbon2] = middle + second*axis
    species.geom[hydrogen] = middle
    species.reset_order()
    assert species.hydrogen_bridges == []


def test_a_hydrogen_with_an_ordinary_heteroatom_bond_is_not_a_carbon_bridge():
    species = point('ethyl')
    hydrogen = species.hydrogen_bridges[0][1]
    atom = list(species.atom)+['O']
    geom = np.vstack([species.geom, species.geom[hydrogen]+[0.,0.,.98]])
    example = StationaryPoint('proton bound to oxygen', 1, 1, atom=atom, geom=geom)
    example.characterize()
    assert example.hydrogen_bridges == []


def test_snapshots_and_database_context_retain_the_three_centre_connection():
    species = point('protonated_ethane')
    snapshot = endpoint_snapshot(species)
    assert snapshot.hydrogen_bridges == species.hydrogen_bridges
    assert canonical_identity(snapshot)['id'] == canonical_identity(species)['id']
    context = _chemical_context(species)
    assert context['hydrogen_bridges'] == [list(value) for value in species.hydrogen_bridges]


@pytest.mark.parametrize('name, residue', [('ethyl','CCHHH'), ('vinyl','CCH'), ('protonated_ethane','CCHHHHH')])
def test_declared_h2_loss_from_the_shared_hydrogen_preserves_the_residue(name, residue, tmp_path):
    species = point(name)
    path = tmp_path/'input.json'
    path.write_text(json.dumps(dict(barrier_threshold=150., high_level=0,
        families=['h2_elim', 'h2_elim_barrierless'], ringrange=[3,4,5,6,7])))
    par = Parameters(str(path), show_warnings=False).par
    ReactionFinder(species, par, None).find_reactions()
    bridge_hydrogen = species.hydrogen_bridges[0][1]
    found = False
    for reaction in species.reac_obj:
        if 'barrierless' not in reaction.instance_name or bridge_hydrogen not in reaction.instance:
            continue
        found = True
        product = reaction.get_barrierless_product()
        fragments, _ = product.start_multi_molecular()
        assert len(fragments) == 2
        assert sorted(''.join(sorted(fragment.atom)) for fragment in fragments) == sorted([residue, 'HH'])
    assert found


def test_initial_ethyl_cation_must_keep_its_input_connectivity(tmp_path, monkeypatch, caplog):
    import logging
    import sys
    from kinbot import kb

    monkeypatch.chdir(tmp_path)
    observed = point('ethyl')
    donor, hydrogen, _ = observed.hydrogen_bridges[0]
    geom = observed.geom.copy()
    # Select a classical C-H bond with the other carbon beyond the bridge range.
    geom[hydrogen] = geom[donor] + [0., 0., 1.09]
    structure = [value for atom, position in zip(observed.atom, geom)
                 for value in (str(atom), *position.tolist())]
    Path('input.json').write_text(json.dumps(dict(structure=structure, charge=1, mult=1,
                                                  high_level=0, barrier_threshold=100)))
    monkeypatch.setattr(sys, 'argv', ['kinbot', 'input.json'])
    monkeypatch.setattr(kb, 'config_log', lambda *a, **k: logging.getLogger('KinBot'))
    calls = []
    monkeypatch.setattr(QuantumChemistry, 'qc_opt',
                        lambda qc, species, geom: calls.append(routing_name(species)))
    monkeypatch.setattr(QuantumChemistry, 'get_qc_geom', lambda *a, **k: (0, observed.geom))
    monkeypatch.setattr(QuantumChemistry, 'get_qc_freq', lambda *a, **k: (0, [100.] * 15))
    kb._main()
    assert len(calls) == 1
    assert calls[0] != routing_name(observed)
    assert "changed connectivity" in caplog.text
    assert calls[0] + "_well in kinbot.db" in caplog.text


def test_pes_creates_the_input_well_directory_before_qc(tmp_path, monkeypatch):
    import logging
    import sys
    from kinbot import pes

    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps(dict(smiles='[CH2+]C', charge=1, mult=1,
                                                  high_level=0, barrier_threshold=100)))
    monkeypatch.setattr(sys, 'argv', ['pes', 'input.json'])
    monkeypatch.setattr(pes, 'config_log', lambda *a, **k: logging.getLogger('KinBot'))
    monkeypatch.setattr(QuantumChemistry, 'qc_opt', lambda *a, **k: pytest.fail('PES ran QC'))
    write_input = pes.write_input

    class InputWritten(Exception):
        pass

    def check_input(path, species, *args):
        write_input(path, species, *args)
        name = routing_name(species)
        assert Path(name, name + '.json').exists()
        assert not Path('kinbot.db').exists()
        raise InputWritten

    monkeypatch.setattr(pes, 'write_input', check_input)
    with pytest.raises(InputWritten):
        pes._main()


def test_migration_of_a_bridging_h_requires_a_new_acceptor(tmp_path):
    species = point('2butyl')
    path = tmp_path / 'input.json'
    path.write_text(json.dumps(dict(barrier_threshold=100., families=['intra_H_migration'])))
    finder = ReactionFinder(species, Parameters(str(path), show_warnings=False).par, None)
    finder.find_reactions()
    first, hydrogen, second = species.hydrogen_bridges[0]
    migrations = finder.reactions['intra_H_migration']
    acceptors = {instance[0] for instance in migrations if instance[-1] == hydrogen}
    assert acceptors
    assert not acceptors.intersection((first, second))
    assert any(instance[-1] != hydrogen for instance in migrations)


def classical_attachment(species):
    requested = copy.deepcopy(species)
    requested.hydrogen_bridges = []
    requested.bond = requested.bonds[0].copy()
    requested.bonds, requested.rads = [requested.bond], [requested.rads[0]]
    requested.bond01 = (requested.bond > 0).astype(int)
    requested.calc_chemid()
    return requested


@pytest.mark.parametrize('name', ['vinyl', 'ethyl', 'protonated_ethane'])
def test_homolytic_cuts_really_separate_bridged_ions(name, tmp_path, caplog):
    species = point(name)
    path = tmp_path/'input.json'
    path.write_text(json.dumps(dict(barrier_threshold=100., families=['hom_sci'])))
    parameters = Parameters(str(path), show_warnings=False).par
    finder = ReactionFinder(species, parameters, None)
    with caplog.at_level('INFO', logger='KinBot'):
        finder.find_reactions()
    if name == 'ethyl':
        assert 'the bridged ion would remain connected' in caplog.text
    for reaction in species.reac_obj:
        product = reaction.get_barrierless_product()
        fragments, _ = product.start_multi_molecular(bond_mx=product.bond)
        assert len(fragments) == 2


def test_ordinary_and_high_level_cache_checks_do_not_accept_bridge_formation(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    path = tmp_path/'input.json'
    path.write_text(json.dumps(dict(barrier_threshold=100.)))
    parameters = Parameters(str(path), show_warnings=False).par
    qc = QuantumChemistry(parameters)
    observed = point('ethyl')
    requested = classical_attachment(observed)
    data = dict(status='normal', energy=-100./constants.EVtoHARTREE, zpe=.02,
                frequencies=[100.]*15, charge=1, multiplicity=1)
    for suffix in ('_well', '_well_high'):
        job = routing_name(requested)+suffix
        qc.db.write(Atoms(observed.atom, positions=observed.geom), name=job, data=data)
        with pytest.raises(StereoRoutingError):
            guard_well_job(qc, requested, requested.geom, job)
    with pytest.raises(StereoRoutingError):
        qc.qc_opt(requested, requested.geom, high_level=1)


@pytest.mark.parametrize('high, conf', [(0, 0), (0, 1), (1, 0)])
def test_pes_reads_bridged_force_worker_energies_without_calculator_charge(tmp_path, monkeypatch, high, conf):
    from ase.db import connect
    from types import SimpleNamespace
    import logging
    from kinbot import pes
    from kinbot.run_format import ensure_current_run

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(pes, 'logger', logging.getLogger('KinBot'), raising=False)
    species = point('ethyl')
    well = routing_name(species)
    ensure_current_run(well, create=True)
    db = connect(well+'/kinbot.db')
    qc = SimpleNamespace(db=db, par={})
    # Force-worker rows do not necessarily contain calculator parameters.
    # The initial-well input reference supplies charge for conf/_low too.
    guard_well_job(qc, species, species.geom, well+'_well')
    job = ('conf/'+well+'_low' if conf else well+'_well') + ('_high' if high else '')
    db.write(Atoms(species.atom, positions=species.geom), name=job, data=dict(
        status='normal', energy=-100./constants.EVtoHARTREE, zpe=.02, frequencies=[100.]*15))
    assert pes.get_energy([well], well, False, high, conf=conf) == pytest.approx((-100., .02))


@pytest.mark.parametrize('name', ['2butyl', 'protonated_ethane', 'cage'])
def test_only_confirmed_spectator_bridges_survive_ts_and_conformer_copies(name):
    from kinbot.hydrogen_bridges import spectator_bridges, copy_ts_bridges
    parent = point(name)
    product = copy.deepcopy(parent)
    # Reproduce TS construction after the endpoints have been identified.
    ts = StationaryPoint('remote reaction TS', parent.charge, parent.mult,
                         atom=parent.atom, geom=parent.geom.copy(), wellorts=1)
    ts.hydrogen_bridges = spectator_bridges(parent, product)
    ts.bond_mx()
    ts.bond = np.maximum(ts.bond, np.maximum(parent.bond, product.bond))
    ts.find_cycle()
    ts.find_dihedral()
    assert ts.hydrogen_bridges == parent.hydrogen_bridges
    assert ts.cycle_chain == parent.cycle_chain
    assert ts.dihed == parent.dihed
    for _, hydrogen, _ in parent.hydrogen_bridges:
        assert all(hydrogen not in ring for ring in ts.cycle_chain)
    dummy = StationaryPoint('TS conformer', parent.charge, parent.mult,
                            atom=parent.atom, geom=parent.geom.copy())
    copy_ts_bridges(ts, dummy)
    dummy.bond_mx()
    assert dummy.hydrogen_bridges == parent.hydrogen_bridges
    # A previous annotation cannot retain a contact that is no longer there.
    dummy.geom[dummy.hydrogen_bridges[0][1]] += [0.,0.,5.]
    dummy.bond_mx()
    assert dummy.hydrogen_bridges == []
    # A bridge that changes during the reaction is not a spectator.
    product.bond = product.bond.copy()
    first, hydrogen, _ = product.hydrogen_bridges[0]
    product.bond[first, hydrogen] = product.bond[hydrogen, first] = 0
    assert spectator_bridges(parent, product) == []


def test_conformer_copy_does_not_discover_a_new_bridge_on_an_h_transfer_ts():
    from kinbot.hydrogen_bridges import copy_ts_bridges
    ts = point('ethyl', wellorts=1)
    dummy = StationaryPoint('H-transfer conformer', 1, 1, atom=ts.atom, geom=ts.geom)
    copy_ts_bridges(ts, dummy)
    dummy.bond_mx()
    assert dummy.hydrogen_bridges == []


def test_confirmed_bridge_is_not_lost_when_a_different_hydrogen_is_transferred():
    from kinbot.hydrogen_bridges import spectator_bridges
    parent = point('protonated_ethane')
    product = copy.deepcopy(parent)
    first, shared, second = parent.hydrogen_bridges[0]
    terminal = next(i for i, symbol in enumerate(parent.atom)
                    if symbol == 'H' and parent.bond[first, i] and i != shared)
    product.bond[first, terminal] = product.bond[terminal, first] = 0
    assert spectator_bridges(parent, product) == parent.hydrogen_bridges


def test_confirmed_bridge_survives_excess_contacts_at_a_remote_ts_centre():
    parent = point('ethyl')
    geom = np.vstack([parent.geom, [[10.,0.,0.], [11.,0.,0.], [12.,0.,0.]]])
    ts = StationaryPoint('bridged ion and H-transfer TS', 1, 2,
                         atom=list(parent.atom)+['H']*3, geom=geom, wellorts=1)
    ts.hydrogen_bridges = list(parent.hydrogen_bridges)
    ts.bond_mx()
    ts.find_cycle()
    assert ts.hydrogen_bridges == parent.hydrogen_bridges
    assert ts.bond[-2,-3] == ts.bond[-2,-1] == 1
    assert all(parent.hydrogen_bridges[0][1] not in ring for ring in ts.cycle_chain)


@pytest.mark.parametrize('mc, pes_mode', [(0, 0), (1, 0), (0, 1), (1, 1)])
def test_public_mess_uses_bridged_ion_symmetry_in_h2_capture(tmp_path, monkeypatch, mc, pes_mode):
    import logging
    import shutil
    import subprocess
    from types import SimpleNamespace
    from kinbot import pes
    from kinbot.mess import MESS
    from tests.test_ion_h2_pst import ion_h2_system, parameters, capture_rates, langevin_rate

    executable = shutil.which('mess')
    if executable is None:
        pytest.skip('Public MESS is not installed')
    parent = point('protonated_ethane')
    ion = point('ethyl')
    _, old_reaction = ion_h2_system(1)
    hydrogen = old_reaction.products[1]
    for species in (parent, ion):
        species.name = routing_name(species)
        symmetry.calculate_symmetry(species)
        species.freq = species.reduced_freqs = [1000.]*(3*species.natom-6)
        species.energy, species.zpe = (-1. if species is parent else -.49), 0.
    assert parent.sigma_ext == ion.sigma_ext == hydrogen.sigma_ext == 2
    name = routing_name(parent)+'_h2_elim_barrierless_1_4'
    reaction = SimpleNamespace(instance_name=name, products=[ion, hydrogen],
        prod_opt=[SimpleNamespace(species=p) for p in (ion, hydrogen)],
        ts=copy.deepcopy(parent), do_vdW=False, mp2=0)
    parent.reac_obj, parent.reac_ts_done = [reaction], [-1]
    parent.reac_type, parent.reac_inst, parent.reac_name = ['h2_elim_barrierless'], [[0,1,3]], [name]
    par = parameters(tmp_path, parent, mc, pes_mode)
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
        threshold = (sum(p.energy+p.zpe for p in reaction.products)-parent.energy)*constants.AUtoKCAL
        pes.create_mess_input(par, [parent_key], [product_key], [],
            [[parent_key, name, product_names, threshold]], [],
            {parent_key: 0.}, {product_key: threshold},
            {parent_key: parent_key, product_key: parent_key}, parent.mass, False)
    run = subprocess.run([executable, 'mess_0000.inp'], cwd='me',
                         capture_output=True, text=True, timeout=30)
    assert run.returncode == 0, run.stdout+run.stderr
    rates = capture_rates(Path('me/mess_0000.out').read_text())
    np.testing.assert_allclose(rates, langevin_rate(reaction.products), rtol=.005)
