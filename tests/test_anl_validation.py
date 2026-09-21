"""The external-site validation graph is general and cannot overclaim ANL."""

from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from ase import Atoms
from ase.db import connect

from kinbot.anl.dispatch import validate_spec
from kinbot.anl.validation import (interface_validation_spec,
                                   main as validation_main,
                                   molecule_from_database)
from kinbot.anl import site
from kinbot.parameters import Parameters
from kinbot.reaction_finder import ReactionFinder
from kinbot.stationary_pt import StationaryPoint


def _molecule():
    return {'symbols': ['C', 'H', 'H', 'H', 'H'],
            'positions': [[0., 0., 0.], [0.63, 0.63, 0.63],
                          [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63],
                          [0.63, -0.63, -0.63]],
            'charge': 0, 'multiplicity': 1}


def test_non_mrcc_graph_has_geometry_barrier_and_parallel_fanout():
    spec = interface_validation_spec(
        _molecule(), max_nodes=4, partition='day-long-cpu')
    resolved = deepcopy(spec)
    for task in resolved['tasks']:
        task['resources'].update(cores=4, memory_mb=64000,
                                 partition='test')
    validate_spec(resolved)
    tasks = {task['id']: task for task in spec['tasks']}
    assert tasks['l3_geometry']['geometry_from'] == 'l2_geometry'
    for ident in ('harmonic', 'f12_tz', 'f12_qz', 'ccsdt_dz',
                  'cfour_dboc'):
        assert tasks[ident]['geometry_from'] == 'l3_geometry'
    assert tasks['gaussian_vpt2']['geometry_from'] == 'l2_geometry'
    assert tasks['gaussian_vpt2']['depends_on'] == ['l3_geometry']
    assert spec['intent']['mrcc_enabled'] is False
    assert spec['intent']['claim'] == 'interface-validation-only'
    assert all(task['resources']['cores'] == 'auto' for task in spec['tasks'])
    assert all(task['resources']['partition'] == 'day-long-cpu'
               for task in spec['tasks'])


def test_database_export_requires_one_complete_accepted_l2_record():
    with TemporaryDirectory() as temporary:
        database = Path(temporary) / 'kinbot.db'
        db = connect(database)
        atoms = Atoms('CH4', positions=_molecule()['positions'])
        db.write(atoms, name='methane_well_high', data={
            'status': 'normal', 'energy': -10., 'zpe': .1,
            'frequencies': [100., 200.]})
        molecule = molecule_from_database(
            database, 'methane_well_high', charge=0, multiplicity=1)
        assert molecule['symbols'] == ['C', 'H', 'H', 'H', 'H']
        assert molecule['charge'] == 0
        db.write(atoms, name='bad_well_high', data={'status': 'error'})
        try:
            molecule_from_database(database, 'bad_well_high', charge=0,
                                   multiplicity=1)
        except ValueError as error:
            assert 'not normal' in str(error)
        else:
            raise AssertionError('An incomplete L2 result was exported.')


def test_ethane_external_example_finds_only_requested_scission():
    root = Path(__file__).resolve().parents[1]
    parameters = Parameters(
        root / 'examples/anl/ethane_profiled_hpc/ethane.json',
        show_warnings=False).par
    species = StationaryPoint(
        'well0', parameters['charge'], parameters['mult'],
        smiles=parameters['smiles'], structure=parameters['structure'])
    species.characterize()
    assert str(species.chemid) == '301020900180000000001'
    species.name = str(species.chemid)
    ReactionFinder(species, parameters, None).find_reactions()
    assert species.reac_type == ['hom_sci']
    assert species.reac_inst == [[0, 1]]


def test_prepare_from_database_cli_stages_general_graph(monkeypatch):
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        database = root / 'kinbot.db'
        connect(database).write(
            Atoms('CH4', positions=_molecule()['positions']),
            name='methane_well_high', data={
                'status': 'normal', 'energy': -10., 'zpe': .1,
                'frequencies': [100., 200.]})
        monkeypatch.setattr(site, '_partitions', lambda: [{
            'name': 'test-cpu', 'default': True, 'cores': 32,
            'memory_mb': 128000, 'seconds': 7 * 24 * 3600}])
        run_dir = root / 'run'
        assert validation_main([
            'prepare-from-db', str(database), 'methane_well_high',
            str(run_dir), '--max-nodes', '2', '--partition', 'test-cpu']) == 0
        workflow = json.loads((run_dir / 'workflow.json').read_text())
        assert workflow['limits']['max_nodes'] == 2
        assert all(task['resources']['partition'] == 'test-cpu'
                   for task in workflow['tasks'])
        state = json.loads((run_dir / 'state.json').read_text())
        assert set(state['tasks']) == {'l2_geometry'}
