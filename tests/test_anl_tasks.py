"""Pinned native task routing for ANL higher-order corrections."""

import pytest

from kinbot.anl.results import validate_result_parser
from kinbot.anl.tasks import cfour_energy_task, higher_order_task, mrcc_task


def _validate(task):
    validate_result_parser(
        task['result_parser'], backend=task['backend'],
        template=task['input_template'], outputs=task['required_outputs'])


def test_closed_shell_ccsdt_q_routes_to_direct_mrcc_rhf_unrestricted_cc():
    task = higher_order_task(
        'hoe_tz_high', 'CCSDT(Q)', 'cc-pVTZ', multiplicity=1)
    assert task['backend'] == 'mrcc'
    assert task['command'] == ['dmrcc']
    assert task['result_parser']['reference'] == 'RHF'
    assert task['result_parser']['correlation'] == 'unrestricted'
    assert task['result_parser']['program'] == 'mrcc'
    assert 'scftype=RHF' in task['input_template']
    assert 'ccprog=mrcc' in task['input_template']
    assert 'core=frozen' in task['input_template']
    assert 'files_from_env' not in task
    _validate(task)


def test_cfour_higher_order_rejects_open_shell_or_other_methods():
    with pytest.raises(ValueError, match='closed-shell'):
        cfour_energy_task('open', 'CCSDT(Q)', 'cc-pVDZ', multiplicity=2)
    with pytest.raises(ValueError, match='Unsupported CFOUR'):
        cfour_energy_task('wrong', 'CCSDTQ(P)', 'cc-pVDZ', multiplicity=1)


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


def test_ccsdtq_p_defaults_request_partition_max_and_eight_openmp_threads():
    task = higher_order_task(
        'hoe_dz_high', 'CCSDTQ(P)', 'cc-pVDZ', multiplicity=1)
    assert task['resources']['walltime'] == 'partition_max'
    assert task['resources']['max_cores'] == 8


def test_long_mrcc_defaults_do_not_replace_explicit_user_walltimes():
    for builder in (mrcc_task, higher_order_task):
        task = builder('high', 'CCSDT(Q)', 'cc-pVTZ', multiplicity=2)
        assert task['resources']['walltime'] == 'partition_max'
        task = builder('high', 'CCSDT(Q)', 'cc-pVTZ', multiplicity=2, walltime='02:00:00')
        assert task['resources']['walltime'] == '02:00:00'


def test_open_shell_ccsdt_q_routes_to_direct_mrcc_rohf_ucc():
    task = higher_order_task(
        'hoe_tz_high', 'CCSDT(Q)', 'cc-pVTZ', multiplicity=2)
    assert task['backend'] == 'mrcc'
    assert task['result_parser']['reference'] == 'ROHF'
    assert task['result_parser']['correlation'] == 'unrestricted'
    assert task['resources']['max_cores'] == 8
    assert task['resources']['min_memory_mb_per_core'] == 4096
    _validate(task)


def test_mrcc_rejects_uhf_determinants_and_accepts_configured_command():
    with pytest.raises(ValueError, match='Unsupported MRCC'):
        mrcc_task('bad', 'CCSD(T)', 'cc-pVDZ', multiplicity=1)
    with pytest.raises(ValueError, match='Unsupported MRCC reference'):
        mrcc_task('bad-uhf', 'CCSDT(Q)', 'cc-pVDZ', multiplicity=2,
                  reference='UHF')
    configured = mrcc_task(
        'configured', 'CCSDT(Q)', 'cc-pVDZ', multiplicity=1,
        command='/software/mrcc/dmrcc')
    assert configured['command'] == ['/software/mrcc/dmrcc']
    with pytest.raises(ValueError, match='without arguments'):
        mrcc_task('bad-command', 'CCSDT(Q)', 'cc-pVDZ', multiplicity=1,
                  command='dmrcc --flag')


def test_no_same_program_or_same_reference_constraint_is_encoded():
    high = higher_order_task(
        'high', 'CCSDT(Q)', 'cc-pVTZ', multiplicity=1)
    lower = {'backend': 'molpro', 'reference': 'RHF'}
    assert (high['backend'], lower['backend']) == ('mrcc', 'molpro')
    assert high['result_parser']['method'] == 'CCSDT(Q)'


@pytest.mark.parametrize('method,count,effective', [
    (method, count, effective)
    for method in ('CCSDT(Q)', 'CCSDTQ(P)')
    for count, effective in [(1, 'HF'), (2, 'CCSD'), (3, 'CCSDT'),
                              (4, 'CCSDTQ' if method == 'CCSDTQ(P)' else method),
                              (5, method), (6, method)]
])
def test_higher_order_rank_plan_preserves_requested_method(method, count, effective):
    task = higher_order_task('probe', method, 'cc-pVDZ',
                             multiplicity=2 if count % 2 else 1, electron_count=count)
    assert f'calc={effective}\n' in task['input_template']
    parser = task['result_parser']
    assert parser['method'] == method
    if effective != method:
        assert parser['rank_exact_electrons'] == count
        assert parser['effective_method'] == effective
    else:
        assert 'effective_method' not in parser
    _validate(task)
