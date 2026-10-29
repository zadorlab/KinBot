"""Replay native product identity changes through the reaction-owning driver.

QC completion/optimization is simulated; the input and result geometries are
exported native evidence. This verifies routing, not final channel validity.
"""
import copy
import json
from collections import Counter
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect

from kinbot import constants
from kinbot.parameters import Parameters
from kinbot.reaction_generator import ReactionGenerator
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import optical_scope
from kinbot.stereo_routing import guard_well_job

CASES = json.loads((Path(__file__).parent/'reference/product_configuration_change.json').read_text())['cases']


def point(smiles, mult=1):
    p = StationaryPoint('fixture', 0, mult, smiles=smiles)
    p.characterize()
    p.energy = p.start_energy = -100.
    p.zpe = p.start_zpe = .01
    p.freq = p.reduced_freqs = [100.] * max(0, 3*p.natom-6)
    return p


@pytest.mark.parametrize('case', CASES, ids=lambda c: c['worker'].split('-')[-1])
@pytest.mark.parametrize('delayed', [False, True])
def test_retired_product_is_not_resubmitted_or_revalidated(case, delayed, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text('{"barrier_threshold": 100.0}')
    par = Parameters('input.json', show_warnings=False).par
    par.update(pes=0, high_level=0, L3_calc=0)
    parent = point('[CH2]SS', 2)
    original = StationaryPoint('product', case['charge'], case['multiplicity'],
                              atom=case['input_symbols'], geom=np.array(case['input_positions']))
    original.characterize()
    for key, value in case['input_context'].items():
        if value is None:
            continue
        if key in ('bond', 'bond01'):
            value = np.array(value)
        elif key in ('bonds', 'rads'):
            value = [np.array(v) for v in value]
        setattr(original, key, copy.deepcopy(value))
    old_job = routing_name(original)+'_well'
    assert old_job == case['job']
    hydrogen = point('[H]', 2)
    endpoint = copy.copy(parent)
    endpoint.start_multi_molecular = Mock(return_value=([original, hydrogen], None))
    reaction = SimpleNamespace(instance_name='hom_sci_fixture', species=parent,
        instance=[1, 2], family_name='hom_sci', max_step=0, scan=False,
        products=[], prod_opt=[], do_vdW=False, irc_prod=endpoint)
    parent.reac_obj = [reaction]
    parent.reac_inst = [[1, 2]]
    parent.reac_type = ['hom_sci']
    parent.reac_name = [reaction.instance_name]
    parent.reac_step = [0]
    parent.reac_scan_energy = [[]]
    parent.reac_ts_done = [2]
    qc = SimpleNamespace(qc="gauss", job_ids={}, qc_opt=Mock(),
        get_qc_geom=Mock(), get_qc_energy=Mock(), get_qc_zpe=Mock())
    qc.db = connect('kinbot.db')
    qc.par = par
    qc._verified_well_jobs = {}
    polls, requests = Counter(), Counter()
    molecules = {routing_name(hydrogen)+'_well': hydrogen}
    original_result = []

    def submit(p, geom):
        job = routing_name(p)+'_well'
        optical_scope(p, 'specified')
        # The actual unchanged production guard rejects wrong-isomer reuse.
        guard_well_job(qc, p, geom, job)
        requests[job] += 1
        molecules[job] = p

    def read_geom(job, *args, **kwargs):
        polls[job] += 1
        p = molecules[job]
        geom = np.array(case['output_positions']) if job == old_job else p.geom.copy()
        status = 1 if delayed and job == old_job and polls[job] <= 3 else 0
        # Complete after submission, during the result-read phase, as in the
        # export. Another loop must poll the replacement without touching this
        # obsolete input contract again.
        if status == 0 and kwargs.get('reorder') and not list(qc.db.select(name=job)):
            rowid = qc.db.write(Atoms(p.atom, positions=geom), name=job,
                data={'status': 'normal', 'energy': -100./constants.EVtoHARTREE,
                      'zpe': .01, 'frequencies': [100.] * max(0, 3*p.natom-6)})
            if job == old_job:
                original_result.append(rowid)
        result = (status, geom)
        return (*result, p.atom.copy()) if kwargs.get('reorder') else result

    qc.qc_opt.side_effect = submit
    qc.get_qc_geom.side_effect = read_geom
    qc.get_qc_energy.return_value = (0, -100.)
    qc.get_qc_zpe.return_value = (0, .01)
    optimizers = []

    def optimize(p, *args):
        p.freq = p.reduced_freqs = [100.] * max(0, 3*p.natom-6)
        opt = SimpleNamespace(species=p, shigh=1, shir=1, do_optimization=Mock())
        optimizers.append(opt)
        return opt

    sleep_count = []
    def tick(*args):
        sleep_count.append(1)
        assert len(sleep_count) < 20, 'Driver did not finish the finite mock workflow'

    monkeypatch.setattr('kinbot.reaction_generator.Optimize', optimize)
    monkeypatch.setattr('kinbot.reaction_generator.time.sleep', tick)
    monkeypatch.setattr('kinbot.reaction_generator.postprocess.createPESViewerInput', Mock())
    ReactionGenerator(parent, par, qc, 'input.json').generate()

    assert parent.reac_ts_done == [-1]
    assert Counter(routing_name(p) for p in reaction.products) == Counter([
        case['observed_route'], routing_name(hydrogen)])
    assert len(optimizers) == 2
    assert requests[old_job] == (2 if delayed else 1)
    assert len(original_result) == 1
    row = qc.db.get(original_result[0])
    assert row.name == old_job and row.data['status'] == 'normal'
    np.testing.assert_array_equal(row.positions, case['output_positions'])
    assert not Path('unsupported_observations').exists()
