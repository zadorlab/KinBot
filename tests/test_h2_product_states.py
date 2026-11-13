"""H2 exits compare ordinary optimized states that conserve total spin."""
import json
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms

from kinbot import constants
from kinbot.parameters import Parameters
from kinbot.product_states import (fragment_state_combinations,
    participating_fragment_states, spin_allowed)
from kinbot.qc import QuantumChemistry
from kinbot.reaction_generator import ReactionGenerator
from kinbot.reactions.reac_h2_elim_barrierless import H2ElimBarrierless
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity
from kinbot.stereo_routing import guard_well_job


def ethyl_endpoint(tmp_path, mult=3):
    filename = tmp_path / 'input.json'
    filename.write_text(json.dumps({'barrier_threshold': 100000., 'high_level': 0,
                                    'conformer_search': 0, 'L3_calc': 0}))
    par = Parameters(filename, show_warnings=False).par
    parent = StationaryPoint('ethyl cation', 1, mult, smiles='C[CH2+]')
    parent.characterize()
    reaction = H2ElimBarrierless(parent, None, par, [2, 0, 3],
        'ethyl_h2_elim_barrierless')
    return parent, reaction.get_barrierless_product(), par


@pytest.mark.parametrize('mult', [1, 3])
def test_fragment_choices_preserve_charge_atom_groups_and_parent_spin(tmp_path, mult):
    parent, endpoint, _ = ethyl_endpoint(tmp_path, mult)
    states, groups, combinations, _ = fragment_state_combinations(endpoint)
    assert len(groups) == 2
    assert sorted(len(group) for group in groups) == [2, 5]
    assert combinations
    for combination in combinations:
        fragments = [states[index] for index in combination]
        assert sum(fragment.charge for fragment in fragments) == parent.charge
        assert spin_allowed([fragment.mult for fragment in fragments], mult)
        assert sorted(atom for fragment in fragments for atom in fragment.atom) == sorted(parent.atom)
    active = participating_fragment_states(endpoint)
    charged_heavy = [point.mult for point in active if point.natom == 5 and point.charge == 1]
    assert charged_heavy == [mult]


def run_products(tmp_path, monkeypatch, *, mult=3, family='h2_elim_barrierless',
                 fail_heavy=False, fail_all=False, rearrange=False, ion_hh=None, dissociate=False):
    monkeypatch.chdir(tmp_path)
    parent, endpoint, par = ethyl_endpoint(tmp_path, mult)
    parent.energy = parent.start_energy = -12.5
    parent.zpe = parent.start_zpe = 0.
    reaction = SimpleNamespace(instance_name='ethyl_' + family, species=parent,
        instance=[2, 0, 3], products=[], prod_opt=[], do_vdW=False, irc_prod=endpoint)
    parent.reac_obj, parent.reac_inst = [reaction], [[2, 0, 3]]
    parent.reac_type, parent.reac_name = [family], [reaction.instance_name]
    parent.reac_step, parent.reac_ts_done = [0], [2]
    qc = QuantumChemistry(par)
    if family == 'h2_elim':
        qc.db.write(Atoms(endpoint.atom, positions=endpoint.geom), name=endpoint.name,
            data=dict(status='normal', energy=parent.energy/constants.EVtoHARTREE, zpe=0.))
    submitted = []
    def status(job):
        row = next(qc.db.select(name=job, sort='-id', limit=1), None)
        return row.data['status'] if row is not None else 0
    qc.check_qc = status
    rearranged = None
    if rearrange:
        rearranged = StationaryPoint('rearranged heavy cation', 1, mult, smiles='C[C+]')
        rearranged.characterize()
    def submit(point, geometry, initial_product=False):
        job = routing_name(point) + ('_product' if initial_product else '_well')
        if status(job) != 0:
            return
        submitted.append((point.natom, point.charge, point.mult))
        guard_well_job(qc, point, geometry, job, initial_product=initial_product)
        failed = fail_all or (fail_heavy and point.natom == 5 and point.charge == 1)
        if failed:
            data = dict(status='error', charge=point.charge, multiplicity=point.mult)
        else:
            energy = (-30. if point.natom == 5 and point.charge == 1 and point.mult == 1
                      else -12. if point.natom == 5 and point.charge == 1
                      else -10. if point.natom == 5
                      else -1. if point.charge == 0 else -.5)
            data = dict(status='normal', energy=energy/constants.EVtoHARTREE,
                        zpe=0., frequencies=[100.] * (3*point.natom-(5 if point.natom == 2 else 6)),
                        charge=point.charge, multiplicity=point.mult)
            if rearranged is not None and point.natom == 5 and point.charge == 1:
                geometry = rearranged.geom.copy()
            if dissociate and point.natom == 5 and point.charge == 1:
                geometry = geometry.copy()
                geometry[1] += [10., 0., 0.]
            if ion_hh is not None and point.natom == 2 and point.charge == 1:
                geometry = geometry.copy()
                axis = geometry[1]-geometry[0]
                geometry[1] = geometry[0] + ion_hh*axis/np.linalg.norm(axis)
        qc.db.write(Atoms(point.atom, positions=geometry), name=job, data=data)
    qc.qc_opt = submit
    class ReadyForConformers(Exception):
        pass
    read = qc.get_qc_geom
    def geometry(job, *args, **kwargs):
        if parent.reac_ts_done[0] == 3:
            raise ReadyForConformers
        return read(job, *args, **kwargs)
    monkeypatch.setattr(qc, 'get_qc_geom', geometry)
    monkeypatch.setattr('kinbot.reaction_generator.time.sleep', lambda _: None)
    if fail_all:
        ReactionGenerator(parent, par, qc, 'input.json').generate()
    else:
        with pytest.raises(ReadyForConformers):
            ReactionGenerator(parent, par, qc, 'input.json').generate()
    return parent, reaction, submitted, rearranged


@pytest.mark.parametrize('mult', [1, 3])
def test_controller_ranks_only_compatible_optimized_products(tmp_path, monkeypatch, mult):
    parent, reaction, submitted, _ = run_products(tmp_path, monkeypatch, mult=mult)
    assert sorted((point.natom, point.charge, point.mult) for point in reaction.products) == [
        (2, 0, 1), (5, 1, mult)]
    assert spin_allowed([point.mult for point in reaction.products], parent.mult)
    if mult == 3:
        assert (5, 1, 1) not in submitted


def test_failed_candidate_leaves_other_compatible_charge_placement(tmp_path, monkeypatch):
    parent, reaction, _, _ = run_products(tmp_path, monkeypatch, fail_heavy=True)
    assert sorted((point.natom, point.charge, point.mult) for point in reaction.products) == [
        (2, 1, 2), (5, 0, 2)]
    assert spin_allowed([point.mult for point in reaction.products], parent.mult)


def test_no_valid_optimized_combination_rejects_the_channel(tmp_path, monkeypatch, caplog):
    parent, _, _, _ = run_products(tmp_path, monkeypatch, fail_all=True)
    assert parent.reac_ts_done == [-999]
    assert 'no optimized product combination preserves' in caplog.text


def test_connected_initial_product_can_rearrange_without_losing_requested_spin(tmp_path, monkeypatch):
    parent, reaction, _, rearranged = run_products(tmp_path, monkeypatch, rearrange=True)
    heavy = next(point for point in reaction.products if point.natom == 5)
    assert canonical_identity(heavy)['id'] == canonical_identity(rearranged)['id']
    assert heavy.charge == 1 and heavy.mult == parent.mult == 3


def test_ordinary_hom_sci_product_ranking_is_unchanged(tmp_path, monkeypatch):
    _, reaction, submitted, _ = run_products(tmp_path, monkeypatch, family='hom_sci')
    # This deliberately records the older parity-only ordinary behavior; the
    # new H2-family policy must not silently alter all homolytic scissions.
    assert sorted((point.natom, point.charge, point.mult) for point in reaction.products) == [
        (2, 0, 1), (5, 1, 1)]
    assert (5, 1, 3) not in submitted


def test_bound_h2_cation_alternative_does_not_discard_the_ordinary_saddle(tmp_path, monkeypatch):
    parent, reaction, submitted, _ = run_products(tmp_path, monkeypatch,
        mult=1, family='h2_elim', ion_hh=1.1142815143)
    # A successful but unselected H2+ optimization must remain one molecule;
    # otherwise its false H+ fragment fails and deletes the valid saddle.
    assert parent.reac_ts_done == [3]
    assert (2, 1, 2) in submitted
    assert not any(natom == 1 for natom, charge, mult in submitted)
    assert sorted((point.natom, point.charge, point.mult) for point in reaction.products) == [
        (2, 0, 1), (5, 1, 1)]


def test_dissociated_residue_keeps_other_connected_charge_alternatives(tmp_path, monkeypatch, caplog):
    parent, reaction, submitted, _ = run_products(tmp_path, monkeypatch, dissociate=True)
    assert sorted((point.natom, point.charge, point.mult) for point in reaction.products) == [
        (2, 1, 2), (5, 0, 2)]
    assert not any(natom == 1 for natom, charge, mult in submitted)
    assert reaction.product_warnings
    assert 'Two-fragment H2 PST cannot assign a rate' in caplog.text
    from kinbot import postprocess
    postprocess.create_summary_file(parent, SimpleNamespace(qc='fc'), dict(high_level=0))
    summary = (tmp_path/('summary_'+routing_name(parent)+'.out')).read_text()
    assert '# WARNING: ' in summary
    assert 'calculation retained and other states continue' in summary
