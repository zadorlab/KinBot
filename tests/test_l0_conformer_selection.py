"""Preliminary conformers must retain their own energy reference."""
import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule

from kinbot import constants
from kinbot.conformers import Conformers
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint


@pytest.fixture
def search(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('conf').mkdir()
    Path('input.json').write_text(json.dumps({'barrier_threshold': 100., 'queuing': 'local'}))
    par = Parameters('input.json').par
    atoms = molecule('H2CO')
    point = StationaryPoint('formaldehyde', 0, 1,
                            atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    point.characterize()
    qc = QuantumChemistry(par)
    search = Conformers(point, par, qc, semi_emp=1)
    search.conf, search.conf_status = 2, [0, 0]
    parent = routing_name(point) + '_well'
    qc.db.write(atoms, name=parent, data={'energy': -3115., 'zpe': .02,
                                          'frequencies': [100.] * 6, 'status': 'normal'})
    Path(parent + '.log').write_text('done\n')
    for i, energy in enumerate((-1., -2.)):
        job = search.get_job_name(i, add='semi_emp_')
        qc.db.write(atoms, name=job, data={'energy': energy, 'zpe': .02,
                                          'frequencies': [100.] * 6, 'status': 'normal'})
        Path(job + '.log').write_text('done\n')
    return search


def test_l0_selects_its_own_lowest_result_without_replacing_l1_parent(search, caplog):
    result = search.check_conformers()
    assert result[1] == '0001'
    assert result[3] == -2. * constants.EVtoHARTREE
    assert search.selected_job == search.get_job_name(1, add='semi_emp_')
    assert len(result[4]) == len(result[5]) == len(result[7]) == 2
    assert result[7] == [0, 0]
    assert not list(search.db.select(name=f'conf/{search.get_name()}_low'))
    assert 'same energy as its parent' not in caplog.text


@pytest.mark.parametrize('failure', ['failed_jobs', 'imaginary_frequencies', 'stereochemistry'])
def test_no_valid_l0_result_does_not_become_an_l1_parent(search, failure, monkeypatch):
    if failure == 'failed_jobs':
        search.conf_status = [1, 1]
    elif failure == 'imaginary_frequencies':
        for i in range(2):
            search.db.write(Atoms(search.species.atom, positions=search.species.geom),
                name=search.get_job_name(i, add='semi_emp_'),
                data={'energy': -1., 'zpe': .02, 'frequencies': [-100.] + [100.] * 5,
                      'status': 'normal'})
    else:
        monkeypatch.setattr('kinbot.conformers.configured_geometry_allowed',
                            lambda *args: False)
    result = search.check_conformers()
    assert result[0] == 1 and result[1] is None
    assert search.selected_job is None and search.selected_conf is None
    assert np.isinf(result[3])
    assert result[7] == [1, 1]
    assert len(result[4]) == len(result[5]) == len(result[7]) == 2
    assert not list(search.db.select(name=f'conf/{search.get_name()}_low'))


def test_l1_parent_recovery_is_unchanged(search):
    search.semi_emp = 0
    search.conf_status = [1, 1]
    result = search.check_conformers()
    assert search.selected_job == routing_name(search.species) + '_well'
    assert result[3] == -3115. * constants.EVtoHARTREE
    assert list(search.db.select(name=f'conf/{search.get_name()}_low'))


@pytest.mark.parametrize('second_valid', [False, True])
def test_rejected_l0_ts_cannot_supply_the_ranking_energy(search, second_valid):
    search.species.wellorts = 1
    search.species.freq = [-1000.] + [100.] * 5
    atoms = Atoms(search.species.atom, positions=search.species.geom)
    search.db.write(atoms, name=search.get_name(),
                    data={'energy': -3115., 'zpe': .02, 'status': 'normal'})
    # The lower-energy result fails the existing imaginary-mode ratio check.
    # It must neither select the parent nor block the valid higher-energy TS.
    for i, freq in enumerate((-500., -1000. if second_valid else -500.)):
        job = search.get_job_name(i, add='semi_emp_')
        search.db.write(atoms, name=job, data={'energy': -2. + i, 'zpe': .02,
            'frequencies': [freq] + [100.] * 5, 'status': 'normal'})
        Path(job + '.log').write_text('done\n')
    result = search.check_conformers()
    assert result[7] == [1, 0 if second_valid else 1]
    if second_valid:
        assert search.selected_job == search.get_job_name(1, add='semi_emp_')
        assert result[3] == -1. * constants.EVtoHARTREE
    else:
        assert search.selected_job is None and result[1] is None
        assert np.isinf(result[3])
