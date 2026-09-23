"""A second reaction may reach a saved configured product with different H indices."""
import copy
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
from ase import Atoms

from kinbot import constants
from kinbot.calculation import load_calculation_record
from kinbot.species_routing import (connect, routing_name, reusable_cached_product,
                                    same_species, prepare_pes_directory)
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity
from kinbot.stereo_routing import guard_well_job


def test_product_adopts_complete_cache_order_and_keeps_endpoint(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    saved = StationaryPoint('product', 0, 1, smiles='C[C@H](O)CC')
    saved.characterize()
    saved.optical_reference = canonical_identity(saved)
    saved.optical_population = 'specified'
    qc = SimpleNamespace(db=connect('kinbot.db'), par={})
    job = routing_name(saved) + '_well'
    guard_well_job(qc, saved, saved.geom, job)
    hessian = np.diag(np.arange(1., 3*saved.natom+1))
    qc.db.write(Atoms(saved.atom, positions=saved.geom), name=job,
        data={'status': 'normal', 'charge': 0, 'multiplicity': 1,
              'energy': -100./constants.EVtoHARTREE, 'zpe': .01,
              'frequencies': [100.]*(3*saved.natom-6), 'hess': hessian})
    # Exchange hydrogens on chemically different atoms, as independent
    # H-transfer paths do. The stereoisomer itself is unchanged.
    hydrogens = np.flatnonzero(np.asarray(saved.atom) == 'H')
    pair = next((i, j) for i in hydrogens for j in hydrogens
                if not np.array_equal(saved.bond01[i], saved.bond01[j]))
    order = np.arange(saved.natom)
    order[list(pair)] = order[list(pair)[::-1]]
    requested = StationaryPoint('product', 0, 1, atom=np.asarray(saved.atom)[order],
                                geom=saved.geom[order].copy())
    requested.characterize()
    endpoint = requested.geom.copy()
    assert same_species(saved, requested)
    assert not np.array_equal(saved.bond01, requested.bond01)
    adopted = reusable_cached_product(qc, requested)
    assert adopted is not requested
    np.testing.assert_array_equal(requested.geom, endpoint)
    np.testing.assert_array_equal(adopted.bond01, saved.bond01)
    guard_well_job(qc, adopted, adopted.geom, job)
    load_calculation_record(adopted, qc, job)
    np.testing.assert_array_equal(adopted.geom, saved.geom)
    np.testing.assert_array_equal(adopted.hess, hessian)
    assert np.isclose(adopted.energy, -100.)


def test_unsupported_graph_can_be_reused_without_inventing_stereo():
    point = StationaryPoint('substituted PAH', 0, 1, smiles='Cc1ccc2cc3ccccc3cc2c1')
    point.characterize()
    assert canonical_identity(point)['status'] == 'unsupported'
    other = copy.copy(point)
    other.geom = point.geom.copy()
    assert same_species(point, other)
    assert canonical_identity(point)['status'] == 'unsupported'
    assert 'unverified' in point.stereo_routing_status
    # A chemid collision must not license a different chemical graph.
    unrelated = StationaryPoint('other', 0, 1, smiles='CO')
    unrelated.characterize()
    unrelated.chemid = point.chemid
    assert not same_species(point, unrelated)


def test_wrong_parent_does_not_discard_a_valid_configured_conformer(tmp_path, monkeypatch):
    from kinbot.conformers import Conformers
    from kinbot.parameters import Parameters
    from kinbot.qc import QuantumChemistry
    monkeypatch.chdir(tmp_path)
    (tmp_path/'conf').mkdir()
    point = StationaryPoint('butanol', 0, 1, smiles='C[C@H](O)CC')
    point.characterize()
    point.name = routing_name(point)
    point.freq = [100.] * (3*point.natom-6)
    (tmp_path/'input.json').write_text('{"barrier_threshold": 100}')
    par = Parameters('input.json', show_warnings=False).par
    qc = QuantumChemistry(par)
    qc.qc = 'fc'
    qc.check_qc = Mock(return_value='normal')
    search = Conformers(point, par, qc)
    def record(job, geom, energy):
        qc.db.write(Atoms(point.atom, positions=geom), name=job,
            data={'status': 'normal', 'energy': energy/constants.EVtoHARTREE,
                  'zpe': .01, 'frequencies': point.freq})
    record(point.name+'_well', point.geom*[-1., 1., 1.], -100.)
    record(search.get_job_name(0), point.geom, -99.)
    search.conf, search.conf_status = 1, [0]
    result = search.check_conformers()
    assert search.selected_job == search.get_job_name(0)
    assert np.isclose(result[3], -99.)
    np.testing.assert_array_equal(result[2], point.geom)


def test_unidentified_legacy_pes_directory_does_not_block_fresh_configured_work(tmp_path):
    point = StationaryPoint('butanol', 0, 1, smiles='C[C@H](O)CC')
    point.characterize()
    legacy = tmp_path/str(point.chemid)
    legacy.mkdir()
    (legacy/'existing.log').write_text('preserved calculation')
    target = prepare_pes_directory(tmp_path, point)
    assert target.name == routing_name(point)
    assert target.is_dir() and not target.is_symlink()
    assert (legacy/'existing.log').read_text() == 'preserved calculation'
