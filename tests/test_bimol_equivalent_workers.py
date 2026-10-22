"""Queued L0 workers receive equivalent-H choices for every calculator."""
import copy
import pickle
from pathlib import Path
from unittest.mock import Mock

import pytest

from kinbot import bimol
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from tests.test_bimol_observed_abstraction import (
    _reactants, _enumerate, _transferred_product)


@pytest.mark.parametrize('calculator, backend', [('gaussian', 'gauss'),
    ('qchem', 'qchem'), ('orca', 'orca'), ('fairchem', 'fc')])
@pytest.mark.parametrize('mode', ['search', 'validation'])
def test_generated_worker_accepts_equivalent_transfer(tmp_path, monkeypatch,
                                                      calculator, backend, mode):
    monkeypatch.chdir(tmp_path)
    species = _reactants('C', '[OH]', 1, 2)
    reaction, = _enumerate(species, tmp_path)
    data = reaction.get_bimol_search_data()
    alternative = next(choice for choice in data['equivalent_bond_changes']
                       if choice['instance'] != reaction.instance)
    product = _transferred_product(species, *alternative['instance'][:2])
    profile = (dict(calculator=calculator, model_path='test-model', device='cpu')
               if calculator == 'fairchem' else
               dict(calculator=calculator, method='hf', basis='sto-3g'))
    par = Parameters('input.json').par
    par.update(bimol=1, l0_profile=profile, queuing='local', queue_job_limit=0)
    submitted = Mock(return_value=0)
    monkeypatch.setattr(QuantumChemistry, 'submit_qc', submitted)
    qc = QuantumChemistry(par)
    if mode == 'search':
        qc.qc_bimol_search(reaction)
    else:
        qc.qc_bimol_irc(reaction, dict(validation_directory='bimol/equivalent',
            candidate={'geometry': species.geom.tolist()}))
    job = submitted.call_args.args[0]
    assert qc.for_job(job).qc == backend

    endpoints = [dict(state=state, converged=True, geometry_angstrom=point.geom.tolist())
                 for state, point in (('reactant', species), ('unexpected', product))]
    monkeypatch.setattr(bimol, '_irc_endpoints', lambda *args: copy.deepcopy(endpoints))
    calculator_factory = Mock(return_value=None)
    monkeypatch.setattr(bimol, 'create_calculator', calculator_factory)
    validate = Mock(wraps=bimol.validate_candidate)
    monkeypatch.setattr(bimol, 'validate_candidate', validate)

    def search(point, conformers_a, conformers_b, instance, formed, broken,
               encounters, calculator, workdir, settings, **kwargs):
        checked = validate(point, point.geom, formed, broken, calculator, workdir, settings)
        return 0, [point.geom], dict(candidates=[dict(validation=checked)])

    monkeypatch.setattr(bimol, 'find_ts_candidates', search)
    source = Path(job + '.py')
    exec(compile(source.read_text(), str(source), 'exec'), {})
    with Path(job + '.pkl').open('rb') as stream:
        result = pickle.load(stream)['data']
    assert result['status'] == 'normal'
    validation = (result['search']['candidates'][0]['validation'] if mode == 'search'
                  else result['validation'])
    assert validation['valid'] and validation['product_matches_expected']
    assert validate.call_args.args[6]['equivalent_bond_changes'] == data['equivalent_bond_changes']
    assert calculator_factory.call_args.args[0] == backend
    assert Path(job + '.done').read_text() == 'done\n'
