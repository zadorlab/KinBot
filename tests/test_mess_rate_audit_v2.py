import json
from kinbot.mess_rate_audit import audit_text


def output():
    return '''Species-Species Rate Tables:
Temperature = 300 K    Pressure = 760 torr

From\\to A B
A 0 1
B 2 0
P 3 4

From\\to P
A -1e-30
B 6
P 0

Temperature-Species Rate Tables:
This repeated layout must not be counted twice.
'''


def test_v2_all_chunks_without_hpl():
    r = audit_text(output())
    assert r['high_pressure']['status'] == 'not_present'
    p = r['pressure_dependent']
    assert p['tables'] == 1 and p['offdiagonal_entries'] == 6
    assert p['negative_count'] == 1 and not p['issues']
    assert p['violations'][0]['target'] == 'P'
    assert p['violations'][0]['value'] == -1e-30
    assert not r['negative_values_clipped']
    json.dumps(r, allow_nan=False)


def test_v2_missing_last_chunk_is_invalid():
    s = output().split('From\\to P')[0]
    r = audit_text(s)
    assert r['pressure_dependent']['issues']


def test_v2_reordered_chunk_is_invalid():
    r = audit_text(output().replace('A -1e-30\nB 6', 'B 6\nA -1e-30'))
    assert r['pressure_dependent']['issues']


def test_v2_positive_is_only_sign_check():
    r = audit_text(output().replace('-1e-30', '0'))
    assert r['pressure_dependent']['status'] == 'passed_sign_check'
    assert r['accuracy_status'] == 'not_assessed'


def test_v2_no_finite_pressure_is_not_success():
    r = audit_text(output().replace('760 torr', '-1 torr'))
    assert r['pressure_dependent']['status'] == 'invalid'


def test_v2_dummy_sink_does_not_require_a_source_row():
    s = output().replace('P 3 4\n', '').replace('P 0\n', '')
    p = audit_text(s)['pressure_dependent']
    assert not p['issues']
    assert p['source_coverage_unverified'][0]['omitted_sources'] == ['P']
    assert p['negative_count'] == 1


def test_v2_negative_loss_is_not_a_valid_generator_diagonal():
    r = audit_text(output().replace('A 0 1', 'A -7 1'))
    assert r['pressure_dependent']['negative_count'] == 2
    assert any(v['diagonal'] for v in r['pressure_dependent']['violations'])


def test_v2_missing_pressure_is_reported_without_crashing():
    r = audit_text(output().replace('    Pressure = 760 torr', ''))
    assert r['pressure_dependent']['status'] == 'invalid'
    assert 'lacks a pressure' in r['pressure_dependent']['issues'][0]['reason']


def test_v2_nonfinite_last_chunk_is_preserved_and_json_safe():
    r = audit_text(output().replace('-1e-30', 'NaN'))
    assert r['pressure_dependent']['nonfinite_count'] == 1
    assert r['pressure_dependent']['violations'][0]['native_value'] == 'NaN'
    json.dumps(r, allow_nan=False)
