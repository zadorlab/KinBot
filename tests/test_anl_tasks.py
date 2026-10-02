"""Pinned native task routing for ANL higher-order corrections."""

import pytest

from kinbot.anl.results import validate_result_parser
from kinbot.anl.tasks import higher_order_task, mrcc_task


def _validate(task):
    validate_result_parser(
        task['result_parser'], backend=task['backend'],
        template=task['input_template'], outputs=task['required_outputs'])


def test_closed_shell_ccsdt_q_routes_to_mrcc_rhf_unrestricted_cc():
    task = higher_order_task(
        'hoe_tz_high', 'CCSDT(Q)', 'cc-pVTZ', multiplicity=1)
    assert task['backend'] == 'mrcc'
    assert task['command'] == ['dmrcc']
    assert task['result_parser']['reference'] == 'RHF'
    assert task['result_parser']['correlation'] == 'unrestricted'
    assert task['result_parser']['program'] == 'mrcc'
    assert 'scftype=RHF' in task['input_template']
    assert 'ccprog=mrcc' in task['input_template']
    _validate(task)


@pytest.mark.parametrize('multiplicity,reference', [(1, 'RHF'), (2, 'ROHF')])
def test_every_ccsdtq_p_routes_to_direct_mrcc_ucc(multiplicity, reference):
    task = higher_order_task(
        'hoe_dz_high', 'CCSDTQ(P)', 'cc-pVDZ',
        multiplicity=multiplicity)
    assert task['backend'] == 'mrcc'
    assert task['command'] == ['dmrcc']
    assert task['required_executables'] == ['scf', 'mrcc']
    assert task['result_parser']['reference'] == reference
    assert task['result_parser']['correlation'] == 'unrestricted'
    assert task['result_parser']['program'] == 'mrcc'
    assert f'scftype={reference}' in task['input_template']
    if reference == 'ROHF':
        assert 'rohftype=semicanonical' in task['input_template']
        assert 'rohfcore=semicanonical' in task['input_template']
    assert 'files_from_env' not in task
    _validate(task)


def test_open_shell_ccsdt_q_routes_to_direct_mrcc_rohf_ucc():
    task = higher_order_task(
        'hoe_tz_high', 'CCSDT(Q)', 'cc-pVTZ', multiplicity=2)
    assert task['backend'] == 'mrcc'
    assert task['result_parser']['reference'] == 'ROHF'
    assert task['result_parser']['correlation'] == 'unrestricted'
    assert task['resources']['max_cores'] == 8
    assert task['resources']['min_memory_mb_per_core'] == 4096
    _validate(task)


def test_mrcc_rejects_unknown_method_and_can_pin_published_uhf_reference():
    with pytest.raises(ValueError, match='Unsupported MRCC'):
        mrcc_task('bad', 'CCSD(T)', 'cc-pVDZ', multiplicity=1)
    task = mrcc_task(
        'paper', 'CCSDT(Q)', 'cc-pVDZ', multiplicity=2, reference='UHF')
    assert task['result_parser']['reference'] == 'UHF'
    assert task['result_parser']['correlation'] == 'unrestricted'


def test_no_same_program_or_same_reference_constraint_is_encoded():
    high = higher_order_task(
        'high', 'CCSDT(Q)', 'cc-pVTZ', multiplicity=1)
    lower = {'backend': 'molpro', 'reference': 'RHF'}
    assert (high['backend'], lower['backend']) == ('mrcc', 'molpro')
    assert high['result_parser']['method'] == 'CCSDT(Q)'
