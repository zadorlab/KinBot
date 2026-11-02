from kinbot.mess_rate_audit import audit_text


def report(loss, a, b):
    return audit_text(rf'''Species-Species Rate Tables:
Temperature = 300 K
From\To A B P
A *** 1 2
B 3 *** 4
P 5 6 ***
Temperature = 300 K    Pressure = 760 torr
From\To A B P
A {loss} {a} {b}
B *** *** ***
P *** *** ***
''')['pressure_dependent']


def test_positive_but_inconsistent_native_rate_is_detected():
    # Native direct double-double experiment: positive isomerization vastly
    # exceeds the total loss. A sign-only test cannot detect this error class.
    r = report('164', '6.87e5', '163')
    assert r['status'] == 'passed_sign_check'
    assert r['loss_balance']['status'] == 'inconsistent'


def test_rounding_intervals_not_an_arbitrary_relative_tolerance():
    assert report('1.00', '.333', '.6674')['loss_balance']['status'] == 'consistent_with_printed_rounding'
    assert report('1.00000', '.333000', '.667400')['loss_balance']['status'] == 'inconsistent'


def test_missing_rate_remains_missing():
    r = report('10', '1', '***')['loss_balance']
    assert r['status'] == 'not_assessed'
    assert r['complete_rows'] == 0


def test_partial_positive_sum_can_already_disprove_balance():
    r = report('1.00', '2.00', '***')['loss_balance']
    assert r['status'] == 'inconsistent'
    assert r['violations'][0]['missing_targets'] == ['P']


def test_native_zero_is_not_given_half_a_unit_tolerance():
    assert report('0', '1e-90', '0')['loss_balance']['status'] == 'inconsistent'


def test_negative_rates_are_retained_even_if_the_row_balances():
    r = report('1', '-2', '3')
    assert r['negative_count'] == 1
    assert r['loss_balance']['status'] == 'consistent_with_printed_rounding'
