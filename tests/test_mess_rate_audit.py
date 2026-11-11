from pathlib import Path
import hashlib
import json
import pytest

from kinbot.mess_rate_audit import audit_text, audit_output
from kinbot import mess_execution
from tests.test_mess_execution import writer


def test_captured_native_negative_rates_are_not_execution_success():
    folder = Path(__file__).parent/'reference'
    p = folder/'mess_negative_pressure.out'
    provenance = json.loads((folder/'mess_negative_pressure.json').read_text())
    report = audit_output(p)
    assert report['output_sha256'] == provenance['output_sha256']
    assert report['high_pressure']['tables'] == 2
    assert report['high_pressure']['status'] == 'passed_sign_check'
    pressure = report['pressure_dependent']
    assert pressure['tables'] == 2 and pressure['negative_count'] == 5
    assert not pressure['issues']
    assert {(v['temperature_k'], v['source'], v['target'], v['value'])
            for v in pressure['violations']} == {
                (300., 'w_2', 'w_1', -1.29e-7),
                (300., 'w_2', 'w_2', -9.18e-7),
                (300., 'w_3', 'w_2', -.00024),
                (1000., 'w_3', 'w_1', -69.5),
                (1000., 'w_3', 'w_2', -33000.)}


def native(rate='-3.30D+04', pressure='760 torr', *, escape=False):
    names = 'A B P' + (' escape' if escape else '')
    tail = ' 4.0' if escape else ''
    return f'''Species-Species Rate Tables:

Temperature = 1000 K

High Pressure Rate Coefficients:

From\\To A B P
A *** 1 2
B 3 *** 4
P 5 6 ***

Temperature = 1000 K    Pressure = {pressure}

From\\To {names}
A 10 {rate} 8.77e7{tail}
B 2 *** ***{tail}
P 1 2 ***{tail}

High Pressure Rate Coefficients (Temperature-Species Rate Tables):
ignored repeated layout
'''


def test_preserves_negative_sign_and_separates_hpl():
    r = audit_text(native())
    assert r['high_pressure']['status'] == 'passed_sign_check'
    p = r['pressure_dependent']
    assert p['status'] == 'invalid' and p['negative_count'] == 1
    v = p['violations'][0]
    assert v['source'] == 'A' and v['target'] == 'B' and v['value'] == -33000
    assert v['negative_to_positive_row_ratio'] == pytest.approx(33000/8.77e7)
    assert not r['automatic_retry'] and not r['negative_values_clipped']
    assert r['accuracy_status'] == 'not_assessed'


@pytest.mark.parametrize('pressure', ['1 bar', '1 atm', '760 torr'])
def test_native_pressure_units_and_escape_columns(pressure):
    r = audit_text(native('1D-14', pressure, escape=True))
    assert r['pressure_dependent']['status'] == 'passed_sign_check'
    assert r['pressure_dependent']['offdiagonal_entries'] == 8


@pytest.mark.parametrize('value', ['nan', '-inf', '+inf'])
def test_nonfinite_is_not_silently_dropped(value):
    r = audit_text(native(value))
    assert r['pressure_dependent']['nonfinite_count'] == 1
    json.dumps(r, allow_nan=False)


def test_tiny_negatives_remain_evidence_without_automatic_retry():
    r = audit_text(native('-1e-99'))
    assert r['pressure_dependent']['violations'][0]['value'] == -1e-99
    assert not r['automatic_retry']


def test_nonfinite_loss_diagonal_is_not_ignored():
    r = audit_text(native('1').replace('A 10 1', 'A nan 1'))
    assert r['pressure_dependent']['status'] == 'invalid'
    assert r['pressure_dependent']['violations'][0]['diagonal']


def test_negative_native_loss_diagonal_is_not_an_ordinary_generator():
    r = audit_text(native('1').replace('A 10 1', 'A -10 1'))
    assert r['pressure_dependent']['status'] == 'invalid'
    assert r['pressure_dependent']['violations'][0]['diagonal']
    assert r['pressure_dependent']['violations'][0]['value'] == -10


def test_overflowing_relative_diagnostic_still_serializes():
    r = audit_text(native('-1e300').replace('8.77e7', '1e-300'))
    json.dumps(r, allow_nan=False)
    assert r['pressure_dependent']['violations'][0]['value'] == -1e300


@pytest.mark.parametrize('change', [
    lambda s: s.replace('From\\To A B P\nA 10', 'From\\To A B B\nA 10'),
    lambda s: s.replace('B 2 *** ***\nP 1 2 ***', 'B 2 *** ***'),
    lambda s: s.replace('B 2 *** ***', 'B 2 ***'),
    lambda s: s.replace('B 2 *** ***', 'X 2 *** ***'),
])
def test_malformed_pressure_table_cannot_borrow_hpl_success(change):
    r = audit_text(change(native('1')))
    assert r['pressure_dependent']['status'] == 'invalid'
    assert r['pressure_dependent']['issues']
    assert r['high_pressure']['status'] == 'passed_sign_check'


def test_no_tables_not_a_validated_calculation():
    r = audit_text('normal exit\n')
    assert r['pressure_dependent']['status'] == 'not_assessed'


def test_native_fixed_width_fields_can_touch_without_changing_values():
    header = 'From\\To'.ljust(8) + ''.join(f'{name:>17}' for name in ('A', 'B', 'P'))
    def row(name, values):
        return f'{name:<8}' + ''.join(f'{v:>17}' for v in values)
    text = native('1')
    old = 'From\\To A B P\nA 10 1 8.77e7\nB 2 *** ***\nP 1 2 ***'
    new = '\n'.join([header, row('A', ['10', '0', '-5.666790091e-126']),
                     row('B', ['***', '-2.088097671e-101', '***']),
                     row('P', ['1', '2', '***'])])
    assert '0-5.666790091e-126' in new
    assert '***-2.088097671e-101' in new
    result = audit_text(text.replace(old, new))['pressure_dependent']
    assert not result['issues']
    assert result['negative_count'] == 2
    assert {entry['native_value'] for entry in result['violations']} == {
        '-5.666790091e-126', '-2.088097671e-101'}


def test_native_output_is_read_only_and_hash_bound(tmp_path):
    p = tmp_path/'native.out'; p.write_text(native())
    before = p.read_bytes()
    r = audit_output(p)
    assert r['output_sha256'] == hashlib.sha256(before).hexdigest()
    assert p.read_bytes() == before


def test_execution_audits_negative_tables_without_resubmitting(tmp_path, monkeypatch):
    obj = writer(tmp_path, monkeypatch, queue='local')
    calls = []
    def run(*args, **kwargs):
        from types import SimpleNamespace
        calls.append(args)
        Path('me/mess_0000.out').write_text(native())
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(mess_execution.subprocess, 'run', run)
    assert mess_execution.run_mess(obj) == 0
    p = Path('me/mess_execution.json')
    data = json.loads(p.read_text())
    assert data['status'] == 'complete'  # Execution is separate from rates.
    assert data['jobs']['mess_0000']['rate_validation']['pressure_dependent']['status'] == 'invalid'
    # Old completed records gain the audit without another solver run.
    del data['jobs']['mess_0000']['rate_validation']; p.write_text(json.dumps(data))
    assert mess_execution.run_mess(obj) == 0
    assert len(calls) == 1
    assert json.loads(p.read_text())['jobs']['mess_0000']['rate_validation']['pressure_dependent']['status'] == 'invalid'
    assert mess_execution.run_mess(obj) == 0 and len(calls) == 1
