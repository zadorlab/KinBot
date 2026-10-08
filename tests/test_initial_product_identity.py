"""The first product result defines its identity; later well reads stay strict.

The reflected result is a controlled database fixture, not a calculated
inversion of 2-butanol. No quantum-chemistry jobs are submitted.
"""
import copy
import json
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms

from kinbot import constants
from kinbot.calculation import load_calculation_record
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.reaction_generator import ReactionGenerator
from kinbot.species_routing import routing_name, configured_result_identity
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity
from kinbot.stereo_routing import guard_well_job, StereoRoutingError


def point(smiles):
    species = StationaryPoint('product fixture', 0, 1, smiles=smiles)
    species.characterize()
    return species


@pytest.mark.parametrize('family', ['hom_sci', 'intra_H_migration'])
@pytest.mark.parametrize('copies', [1, 2])
def test_first_product_optimization_identifies_final_structure(tmp_path, monkeypatch, family, copies):
    monkeypatch.chdir(tmp_path)
    (tmp_path/'input.json').write_text(json.dumps({'barrier_threshold': 100.}))
    par = Parameters('input.json', show_warnings=False).par
    par.update(high_level=0, L3_calc=0)
    qc = QuantumChemistry(par)
    def status(name):
        row = next(qc.db.select(name=name, sort='-id', limit=1), None)
        return row.data['status'] if row is not None else 0
    qc.check_qc = status
    initial = point('C[C@H](O)CC')
    initial.optical_population = 'specified'
    initial.optical_reference = canonical_identity(initial)
    source_job = routing_name(initial) + '_well'
    guard_well_job(qc, initial, initial.geom, source_job)
    final = StationaryPoint('optimized product', 0, 1, atom=initial.atom,
                            geom=initial.geom * [-1., 1., 1.])
    final.characterize()
    target_job = routing_name(final) + '_well'
    hessian = np.diag(np.arange(1., 3*final.natom+1))
    data = {'status': 'normal', 'energy': -100./constants.EVtoHARTREE,
            'zpe': .02, 'frequencies': [100.] * (3*final.natom-6),
            'hess': hessian, 'charge': 0, 'multiplicity': 1}
    source_row = qc.db.write(Atoms(final.atom, positions=final.geom), name=source_job, data=data)
    with pytest.raises(StereoRoutingError, match='different stereoisomer'):
        guard_well_job(qc, initial, initial.geom, source_job)
    submissions = []
    def submit(species, geometry, initial_product=False):
        submissions.append((routing_name(species), initial_product))
        guard_well_job(qc, species, geometry, routing_name(species)+'_well',
                       initial_product=initial_product)
    qc.qc_opt = submit
    parent = point('CCO')
    parent.energy, parent.zpe = -100., .02
    parent.start_energy, parent.start_zpe = parent.energy, parent.zpe
    endpoint = copy.copy(initial)
    endpoint.name = source_job
    endpoint.start_multi_molecular = Mock(return_value=([initial]*copies, None))
    endpoint_geom = endpoint.geom.copy()
    reaction = SimpleNamespace(instance_name=family+'_fixture', species=parent,
        instance=[0, 1], products=[], prod_opt=[], do_vdW=False, irc_prod=endpoint)
    parent.reac_obj, parent.reac_inst = [reaction], [[0, 1]]
    parent.reac_type, parent.reac_name = [family], [reaction.instance_name]
    parent.reac_step, parent.reac_ts_done = [0], [2]
    class ReadyForConformers(Exception):
        pass
    read_geometry = qc.get_qc_geom
    def read(name, *args, **kwargs):
        if parent.reac_ts_done[0] == 3:
            raise ReadyForConformers
        return read_geometry(name, *args, **kwargs)
    qc.get_qc_geom = read
    monkeypatch.setattr('kinbot.reaction_generator.time.sleep', lambda _: None)
    with pytest.raises(ReadyForConformers):
        ReactionGenerator(parent, par, qc, 'input.json').generate()
    assert len(reaction.products) == copies
    assert all(routing_name(product)+'_well' == target_job for product in reaction.products)
    assert (routing_name(initial), True) in submissions
    assert (routing_name(final), False) in submissions
    target = list(qc.db.select(name=target_job))[-1]
    assert target.data['copied_from_row_id'] == source_row
    assert configured_result_identity(qc.db, routing_name(final), target) is not None
    for product in reaction.products:
        guard_well_job(qc, product, product.geom, target_job)
        load_calculation_record(product, qc, target_job)
        np.testing.assert_array_equal(product.geom, final.geom)
        np.testing.assert_array_equal(product.hess, hessian)
        assert product.energy == pytest.approx(-100.)
        assert product.zpe == .02
    np.testing.assert_array_equal(reaction.irc_product_reference.geom, endpoint_geom)
    with pytest.raises(StereoRoutingError, match='different stereoisomer'):
        guard_well_job(qc, initial, initial.geom, source_job)
