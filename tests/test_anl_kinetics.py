import hashlib
import json
import math
from pathlib import Path

import pytest

from kinbot.anl.kinetics import (audit_methyl_recombination,
                                 methyl_recombination_experimental_fit,
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
    assert all(math.isclose(row['factor_difference'], 1.)
               for row in audit['rates'])
    native.write_text(native.read_text() + '\nchanged\n')
    with pytest.raises(ValueError, match='hash'):
        audit_methyl_recombination(tmp_path)
