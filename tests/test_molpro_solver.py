"""Native tolerance contracts must not inherit the SCF energy as CI accuracy."""
import json
from copy import deepcopy
from pathlib import Path

import pytest

from kinbot.molpro_solver import apply_solver_options, validate_solver_options
from kinbot.molpro import _vrc_multireference_method


CONTROLS = {
    'scf': {'energy_convergence': 1.e-8, 'max_iterations': 2000},
    'casscf': {'energy_convergence': 1.e-8, 'gradient_convergence': 1.e-5, 'max_iterations': 999},
    'correlation': {'energy_convergence': 1.e-8, 'coefficient_convergence': 1.e-5,
                    'max_iterations': 2000, 'max_reference_iterations': 2000, 'max_reference_vectors': 60},
}

DECK = '''gthresh,energy=1.d-10,coeff=1.d-10
{rhf,so-sci,maxit=100;wf,25,1,1;orbital,2100.2}
{multi,so-sci;start,2140.2;occ,14;closed,10;wf,25,1,1;state,3;
accuracy,energy=1.d-8
natorb,2140.2
maxit,150;}
{rs2-f12,shift=0.3,ipea=0.25;orbital,2141.2,type=canonical;wf,25,1,1;state,1,3;option,maxit=150,maxiti=200;}
{mrci-f12,scalf12=1,cabs_singles=1;orbital,2140.2,type=natural;wf,25,1,1;state,1,3;refstate,3,1,2,3;option,maxit=150,maxiti=200;option,cluster=0;}
'''


def test_native_controls_preserve_chemistry_and_are_idempotent():
    result = apply_solver_options(DECK, CONTROLS)
    assert apply_solver_options(result, CONTROLS) == result
    assert result.count('thrvar=1.d-10') == 2
    assert result.count('thrden=1.d-8') == 2
    assert result.count('maxiti=2000') == 2
    assert result.count('maxvi=60') == 2
    assert 'gthresh,energy=1.d-8,coeff=1.d-5' in result
    assert '{rhf,so-sci,maxit=2000,energy=1.d-8;' in result
    assert 'accuracy,energy=1.d-8,gradient=1.d-5' in result
    assert 'maxit,999;' in result
    for chemistry in ('shift=0.3,ipea=0.25', 'state,1,3', 'refstate,3,1,2,3',
                      'start,2140.2;occ,14;closed,10;wf,25,1,1;state,3;',
                      'orbital,2141.2,type=canonical', 'option,cluster=0'):
        assert chemistry in result


@pytest.mark.parametrize('method', ['CASPT2(5,4)', 'CASPT2-F12(5,4)', 'MRCI+Q(5,4)', 'MRCI-F12+Q(5,4)'])
def test_production_and_recon_share_explicit_mr_tolerances(method):
    deck, _ = _vrc_multireference_method(method, 25, 1)
    assert 'thrden=1.d-7,thrvar=1.d-10' in deck
    deck, _ = _vrc_multireference_method(method, 25, 1, state_model={
        'state_average': [{'multiplicity': 2, 'roots': 1, 'weights': [1]}],
        'selected_roots': {'2': 1},
        'correlation_options': {'energy_convergence': 1.e-9, 'coefficient_convergence': 1.e-6}})
    assert 'thrden=1.d-9,thrvar=1.d-12' in deck


def test_iteration_only_override_does_not_change_explicit_tolerances_or_casci():
    deck = '{mrci-f12,thrden=1.d-9,thrvar=1.d-12;option,maxit=150;}\n{multi,cas-ci;maxit,50;}'
    result = apply_solver_options(deck, {'correlation': {'max_iterations': 2000}, 'casscf': CONTROLS['casscf']})
    assert result == deck.replace('maxit=150', 'maxit=2000')


def test_external_davidson_override_preserves_internal_space_and_cabs_method():
    options = deepcopy(CONTROLS)
    options['correlation']['max_external_vectors'] = 20
    deck = DECK.replace('cabs_singles=1', 'cabs_singles=2')
    result = apply_solver_options(deck, options)
    assert result.count('maxdav=20') == 2 and result.count('maxvi=60') == 2
    assert 'cabs_singles=2' in result and 'option,cluster=0' in result
    assert apply_solver_options(result, options) == result
    with pytest.raises(ValueError):
        validate_solver_options({'correlation': {'cabs_singles': 2}})


@pytest.mark.parametrize('cabs', [1, 2, 3])
def test_cabs_control_keeps_relaxed_davidson_and_reference_contract(cabs):
    model = {'state_average': [{'multiplicity': 2, 'roots': 3}],
             'selected_roots': {'2': 2}, 'correlation_options': {
                 'reference_policy': 'state_average', 'cabs_singles': cabs,
                 'max_external_vectors': 20, 'max_reference_vectors': 60}}
    deck, variable = _vrc_multireference_method('MRCI-F12+Q(5,4)', 25, 1, state_model=model)
    assert f'cabs_singles={cabs}' in deck
    assert 'maxdav=20' in deck and 'maxvi=60' in deck
    assert 'state,1,2;refstate,3,1,2,3;option,nstati=3;' in deck
    assert 'option,cluster=0' in deck and 'KBVRCSEL=energd(1)' in deck
    assert variable == 'KBVRCSEL'
    for method in ('CASPT2-F12(5,4)', 'MRCI+Q(5,4)'):
        with pytest.raises(ValueError, match='specific to MRCI-F12'):
            _vrc_multireference_method(method, 25, 1, state_model=model)


@pytest.mark.parametrize('controls', [{'cabs_singles': v} for v in (True, 0, 4, 2.0, '2')]
                         + [{'max_external_vectors': v} for v in (True, 0, -2, 20.0)])
def test_invalid_new_correlation_controls(controls):
    from kinbot.molpro import validate_mr_correlation_options
    with pytest.raises(ValueError):
        validate_mr_correlation_options(controls)


@pytest.mark.parametrize('options', [{'scf': {'max_iterations': True}},
    {'correlation': {'energy_convergence': 0}}, {'correlation': {'coefficient_convergence': float('nan')}},
    {'correlation': {'reference_policy': 'state_average'}}, {'casscf': {'ignore_error': True}}])
def test_invalid_numerical_overrides_are_rejected(options):
    with pytest.raises(ValueError):
        validate_solver_options(options)


@pytest.mark.parametrize('limit', [1000, 2000, 0, -1, True, 999.0])
def test_mcscf_backend_limit_is_checked_in_all_public_settings(limit, tmp_path):
    from kinbot.molpro import validate_vrc_state_model
    from kinbot.molpro_reconnaissance import _load
    with pytest.raises(ValueError):
        validate_solver_options({'casscf': {'max_iterations': limit}})
    with pytest.raises(ValueError):
        validate_vrc_state_model({'casscf_options': {'max_iterations': limit},
                                 'state_average': [{'multiplicity': 2, 'roots': 1}]})
    source = tmp_path / 'probe.json'
    source.write_text(json.dumps({'schema': 1, 'title': 'atom', 'symbols': ['H'],
        'geometry_angstrom': [[0., 0., 0.]], 'charge': 0, 'multiplicity': 2,
        'bases': [{'name': 'dz', 'basis': 'cc-pVDZ'}], 'anchors': [{'name': 'atom'}],
        'casscf_options': {'max_iterations': limit}}))
    with pytest.raises(ValueError):
        _load(source)


def test_legacy_illegal_multi_limit_requires_explicit_repair():
    deck = DECK.replace('maxit,150;', 'maxit,2000;')
    with pytest.raises(ValueError, match='999.*Molpro MULTI'):
        apply_solver_options(deck, {})
    repaired = apply_solver_options(deck, {'casscf': {'max_iterations': 999}})
    assert 'maxit,999;' in repaired
    assert 'maxiti=200' in repaired


def test_competing_native_threshold_cards_are_rejected():
    with pytest.raises(ValueError, match='Competing MR'):
        apply_solver_options('gthresh,energy=1.d-10\n{mrci-f12;thresh,coeff=5;}', {'correlation': CONTROLS['correlation']})


def test_execution_binds_requested_solver_options(tmp_path, monkeypatch):
    from kinbot import molpro_reconnaissance as recon
    from kinbot.execution import file_sha256
    from types import SimpleNamespace
    inp = tmp_path / 'probe.inp'
    inp.write_text('file,2,probe.wfu\n' + DECK)
    job = {'name': 'probe', 'directory': str(tmp_path), 'input_sha256': file_sha256(inp),
           'ppn': 1, 'retry_solver_options': deepcopy(CONTROLS)}
    monkeypatch.setattr(recon, 'qc_runtime_environment', lambda *a, **k: ({}, {'scratch_directory': str(tmp_path)}))
    monkeypatch.setattr(recon, 'cleanup_qc_runtime', lambda _: None)
    monkeypatch.setattr(recon.subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=0))
    recon._execute(job)
    execution = json.loads((tmp_path / 'execution.json').read_text())
    assert execution['retry_solver_options'] == CONTROLS
    assert recon._checkpoint_evidence(tmp_path, job, execution)[1] is None
    assert file_sha256(inp) == job['input_sha256']
    job['retry_solver_options']['correlation']['coefficient_convergence'] = 1.e-6
    assert 'differs' in recon._checkpoint_evidence(tmp_path, job, execution)[1]
