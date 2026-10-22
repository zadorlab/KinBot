import json
from types import SimpleNamespace

import pytest

from kinbot.execution import atomic_json
from kinbot.vrc_review_status import acknowledge_installed


def setup_review(tmp_path):
    from test_rotdpy_profiles import sulfur_correction
    from test_vrc_multistate import model
    from kinbot import run_policy
    par, _, _ = sulfur_correction(tmp_path)
    par.update(vrc_tst_correction_model='direct', vrc_tst_spin_splitting=None,
               vrc_tst_surface_models={}, vrc_tst_surface_model=model(),
               vrc_tst_surface_model_file='', vrc_tst_surface_model_scope='123')
    species = SimpleNamespace(name='123', wellorts=1, atom=['C','S','H','H','H'],
                              charge=0, mult=2, smiles='C[S]')
    directory = tmp_path/'123'
    request = directory/'vrctst/orbital_review/123_hom_sci_1_2/request.json'
    record = dict(status='needs_review', parent='123', reaction='123_hom_sci_1_2',
                  symbols=species.atom, charge=0, multiplicity=2, bond=[1,2],
                  products=[dict(smiles='[CH3]', charge=0, multiplicity=2),
                            dict(smiles='[S]', charge=0, multiplicity=3)])
    atomic_json(request, record)
    run_policy.update(tmp_path, max_jobs=23, orbital_review_reserve=3)
    return par, species, directory, request


def test_installed_validated_channel_releases_capacity_before_scan(tmp_path):
    from kinbot import run_policy
    par, species, directory, request = setup_review(tmp_path)
    assert run_policy.read(tmp_path)['max_jobs'] == 20
    scan = directory/'vrctst/scan_review/channel/scan_review.json'
    atomic_json(scan, dict(status='needs_review', error='connectivity'))
    before = scan.read_bytes()
    reports = acknowledge_installed(par, species, directory)
    assert reports == [dict(reaction=request.parent.name, status='model_accepted')]
    assert run_policy.read(tmp_path)['max_jobs'] == 23
    assert scan.read_bytes() == before
    accepted = request.read_bytes()
    assert acknowledge_installed(par, species, directory) == []
    assert request.read_bytes() == accepted


@pytest.mark.parametrize('invalid', ['parent','symbols','bond','unreviewed','checkpoint','spin_review'])
def test_unvalidated_models_do_not_release_capacity(tmp_path, invalid):
    from kinbot import run_policy
    par, species, directory, request = setup_review(tmp_path)
    record = json.loads(request.read_text())
    if invalid == 'parent': record['parent'] = '456'
    if invalid == 'symbols': record['symbols'] = list(reversed(species.atom))
    if invalid == 'bond': record['bond'] = [1,9]
    if invalid == 'unreviewed': par['vrc_tst_surface_model']['reviewed'] = False
    if invalid == 'checkpoint': (tmp_path/'avtz-f12.wfu').write_bytes(b'changed')
    if invalid == 'spin_review': par['vrc_tst_correction_model'] = 'spin_splitting'
    atomic_json(request, record)
    before = request.read_bytes()
    reports = acknowledge_installed(par, species, directory)
    assert reports[0]['status'] == 'needs_review'
    assert request.read_bytes() == before
    assert run_policy.read(tmp_path)['max_jobs'] == 20
