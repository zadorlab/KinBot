import hashlib
import json
import math
from pathlib import Path

import pytest

from kinbot.anl.kinetics import (audit_methyl_recombination,
                                 audit_rate_tables,
                                 methyl_recombination_experimental_fit,
                                 methyl_recombination_legacy_fit,
                                 methyl_recombination_theory_fit,
                                 parse_high_pressure_rates, rate_series)


def output(rates):
    rows = '\n'.join(
        f'{temperature:g} *** {rate:.8e}'
        for temperature, rate in rates)
    return f'''preamble
Unimolecular Rate Units: 1/sec;  Bimolecular Rate Units: cm^3/sec

High Pressure Rate Coefficients (Temperature-Species Rate Tables):

Reactant = w_1
T(K) w_1 b_1
300 *** 1.0

Reactant = b_1
T(K) b_1 w_1
{rows}

Capture/Escape Rate Coefficients:
'''


def test_parse_high_pressure_rate_table():
    points = parse_high_pressure_rates(output([(300., 5.e-11), (1000., 2.e-11)]))
    series = rate_series(points, 'b_1', 'w_1')
    assert [point.temperature_k for point in series] == [300., 1000.]
    assert [point.rate for point in series] == pytest.approx([5.e-11, 2.e-11])
    with pytest.raises(ValueError, match='no high-pressure'):
        parse_high_pressure_rates('not a rate output')


def test_methyl_recombination_audit_uses_hash_and_network(tmp_path):
    directory = tmp_path / 'me'
    directory.mkdir()
    rates = [(temperature, methyl_recombination_experimental_fit(temperature))
             for temperature in (300., 1000., 1700.)]
    native = directory / 'mess_0000.out'
    native.write_text(output(rates))
    digest = hashlib.sha256(native.read_bytes()).hexdigest()
    (directory / 'mess_execution.json').write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'jobs': {'mess_0000': {'output_sha256': digest}},
    }))
    (directory / 'mess_networks.json').write_text(json.dumps([{
        'stem': 'mess_0000', 'wells': ['w_1'], 'products': ['b_1'],
    }]))
    audit = audit_methyl_recombination(tmp_path)
    assert audit['status'] == 'passed'
    assert audit['monotonically_nonincreasing_with_temperature']
    assert all(math.isclose(row['factor_difference'], 1., rel_tol=1.e-8)
               for row in audit['rates'])
    assert audit['source']['doi'] == '10.1021/acs.jpca.5b01002'
    assert all(row['similar_vrc_theory_rate'] > 0.
               and row['legacy_1996_rate'] > 0.
               for row in audit['rates'])
    # The general audit permits arbitrary networks and rate trends; no
    # implicit fit or methyl statistical factor is an acceptance criterion.
    (directory / 'mess_networks.json').write_text(json.dumps([{
        'stem': 'mess_0000', 'wells': ['w_1', 'w_2'], 'products': ['b_1', 'b_2'],
    }]))
    generic = audit_rate_tables(tmp_path)
    assert generic['status'] == 'passed'
    assert generic['accuracy_status'] == 'not_assessed' and 'source' not in generic
    assert len(generic['rates']) == 4
    native.write_text(native.read_text() + '\nchanged\n')
    with pytest.raises(ValueError, match='hash'):
        audit_methyl_recombination(tmp_path)


def test_methyl_recombination_audit_reads_legacy_bare_execution_key(tmp_path):
    directory = tmp_path / 'me'
    directory.mkdir()
    rates = [(temperature, methyl_recombination_experimental_fit(temperature))
             for temperature in (300., 1000., 1700.)]
    native = directory / 'mess_0000.out'
    native.write_text(output(rates))
    digest = hashlib.sha256(native.read_bytes()).hexdigest()
    (directory / 'mess_execution.json').write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'jobs': {'0000': {'output_sha256': digest}},
    }))
    (directory / 'mess_networks.json').write_text(json.dumps([{
        'stem': 'mess_0000', 'wells': ['w_1'], 'products': ['b_1'],
    }]))

    audit = audit_methyl_recombination(tmp_path)

    assert audit['status'] == 'passed'
    assert audit['job'] == 'mess_0000'
    assert audit['execution_job_key'] == '0000'


def test_methyl_recombination_retains_rates_outside_fit_range(tmp_path):
    directory = tmp_path / 'me'
    directory.mkdir()
    rates = [(temperature, methyl_recombination_experimental_fit(temperature))
             for temperature in (250., 300., 1000., 2000., 2100.)]
    native = directory / 'mess_0000.out'
    native.write_text(output(rates))
    digest = hashlib.sha256(native.read_bytes()).hexdigest()
    (directory / 'mess_execution.json').write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'jobs': {'mess_0000': {'output_sha256': digest}},
    }))
    (directory / 'mess_networks.json').write_text(json.dumps([{
        'stem': 'mess_0000', 'wells': ['w_1'], 'products': ['b_1'],
    }]))
    audit = audit_methyl_recombination(tmp_path)
    assert audit['status'] == 'passed'
    assert len(audit['rates']) == 5
    assert [row['compared_to_experiment'] for row in audit['rates']] == [
        False, True, True, True, False]
    assert 'reference_rate' not in audit['rates'][-1]


def test_methyl_recombination_published_fits():
    assert methyl_recombination_experimental_fit(298.) == pytest.approx(
        5.76e-11)
    assert methyl_recombination_theory_fit(300.) == pytest.approx(
        5.81792789939e-11)
    assert methyl_recombination_legacy_fit(300.) == pytest.approx(
        5.798148533973406e-11)


def test_production_style_rates_pass_current_not_superseded_fit(tmp_path):
    directory = tmp_path / 'me'
    directory.mkdir()
    rates = [
        (300., 5.55e-11), (500., 4.72e-11), (1000., 4.18e-11),
        (1500., 4.10e-11), (2000., 4.09e-11),
    ]
    native = directory / 'mess_0000.out'
    native.write_text(output(rates))
    digest = hashlib.sha256(native.read_bytes()).hexdigest()
    (directory / 'mess_execution.json').write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'jobs': {'mess_0000': {'output_sha256': digest}},
    }))
    (directory / 'mess_networks.json').write_text(json.dumps([{
        'stem': 'mess_0000', 'wells': ['w_1'], 'products': ['b_1'],
    }]))

    audit = audit_methyl_recombination(tmp_path)

    assert audit['status'] == 'passed'
    assert max(row['factor_difference'] for row in audit['rates']) < 1.5
    assert max(row['legacy_1996_factor_difference']
               for row in audit['rates']) > 5.
