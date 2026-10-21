"""The benchmark must reject incomplete counting even with existing files."""
import copy

from tests.bimol_export_checks import check_thermochemistry


def atom_record():
    return dict(atoms=['H'], raw_harmonic_frequencies_cm1=[],
        **{'raw_harmonic_frequencies_cm-1': [], 'thermochemical_frequencies_cm-1': []},
        electronic_energy_hartree=-.5, zpe_hartree=0., sigma_ext=1., hir={'rotors': []},
        optical_counting=dict(status='resolved', remaining_multiplier=1.))


def check(data):
    errors, warnings = [], []
    check_thermochemistry(data, 'fixture', errors, warnings)
    return errors


def test_atomic_mc_error_cannot_pass_export_check():
    data = atom_record()
    assert not check(data)
    data['conformer_counting_error'] = 'NoneType is not iterable'
    data['optical_counting'].update(status='unresolved', remaining_multiplier=None)
    assert any('conformer counting' in error for error in check(data))
    assert any('optical counting' in error for error in check(data))


def test_measured_hir_mirror_is_not_multiplied_again():
    data = atom_record()
    data['optical_counting'].update(measured_coverage={'mirror_point': 6}, remaining_multiplier=2.)
    assert any('measured HIR mirror' in error for error in check(data))


def test_explicit_mc_mirrors_require_individual_weights_and_frequencies():
    data = atom_record()
    data['atoms'] = ['O', 'H']
    data['raw_harmonic_frequencies_cm-1'] = [3000.]
    data['thermochemical_frequencies_cm-1'] = [3000.]
    data['conformer_inventory'] = [dict(member_id=str(i), retained=True,
        frequencies_cm1=[3000.], zero_energy_hartree=-10.) for i in range(2)]
    data['optical_counting'] = dict(status='per_member', remaining_multiplier=None,
        members=[dict(member_id=str(i), remaining_optical_weight=1.,
                      mirror_partner_ids=[str(1-i)]) for i in range(2)])
    assert not check(data)
    broken = copy.deepcopy(data)
    broken['optical_counting']['members'][0]['remaining_optical_weight'] = 2.
    assert any('explicit mirror' in error for error in check(broken))
    broken = copy.deepcopy(data)
    broken['conformer_inventory'][1]['frequencies_cm1'] = []
    assert any('incomplete harmonic frequencies' in error for error in check(broken))


def test_accepted_soft_imaginary_correction_preserves_ts_reaction_coordinate():
    data = atom_record()
    data['atoms'] = ['O', 'H', 'H']
    data['raw_harmonic_frequencies_cm-1'] = [-1000., -26., 3000.]
    data['thermochemical_frequencies_cm-1'] = [-1000., 26., 3000.]
    errors = []
    check_thermochemistry(data, 'ts', errors, [], is_ts=True)
    assert not errors
    data['thermochemical_frequencies_cm-1'][0] = 1000.
    check_thermochemistry(data, 'ts', errors, [], is_ts=True)
    assert any('harmonic treatment changed frequencies' in error for error in errors)
