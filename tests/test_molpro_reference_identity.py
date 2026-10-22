import json

import numpy as np
import pytest

from kinbot.molpro_reference_identity import compare, matrix, report
from kinbot.execution import atomic_json, file_sha256


@pytest.mark.parametrize('multiplicity', [1, 2, 3, 4, 5])
def test_native_spin_conventions_and_legacy_repair(multiplicity):
    from kinbot.molpro_reference_identity import natural_cards, export_deck, repair_legacy_selectors
    model = {'state_average': [{'multiplicity':multiplicity, 'roots':2}]}
    ms2 = multiplicity-1
    cards = natural_cards(model, multiplicity)
    deck = cards+export_deck('probe', model, multiplicity)
    assert f'natorb,3001.2,state=1.1,spin={ms2}' in deck
    assert f'load,dr1,den,3001.2,state=1.1,ms2={ms2}' in deck
    assert f'load,dr2,den,3002.2,state=2.1,ms2={ms2}' in deck
    legacy = deck.replace(',ms2=', ',spin=')
    # A general/user MATROP selector really means S and must not be rewritten.
    unrelated = '\nload,user,den,2140.2,state=1.1,spin=1\n'
    fixed, repairs = repair_legacy_selectors(legacy+unrelated)
    assert fixed == deck+unrelated
    assert [r['ms2'] for r in repairs] == [ms2, ms2]
    assert repair_legacy_selectors(fixed) == (fixed, [])


def test_native_density_matching_detects_wrong_root_and_is_rotation_invariant(tmp_path):
    # Four occupied spatial orbitals, five electrons, quartet. Reference has
    # its three SOMOs inside the CAS; root 2, not labelled root 1, matches it.
    rng = np.random.default_rng(19)
    q, _ = np.linalg.qr(rng.normal(size=(7, 7)))
    s = np.diag(np.linspace(.7, 1.6, 7))
    hf = np.diag(1/np.sqrt(np.diag(s)))@q
    occ = np.array([2.,1,1,1,0,0,0])
    alternate = np.array([2.,1,1,0,1,0,0])
    a = {'s':s, 'hf':hf, 'cas':hf, 'dhf':hf@np.diag(occ)@hf.T,
         'dr1':hf@np.diag(alternate)@hf.T, 'dr2':hf@np.diag(occ)@hf.T}
    result = compare(a,electrons=5,active_electrons=3,active_orbitals=4,multiplicity=4)
    assert result['closest_density_root'] == 2
    assert result['reference_states'][0]['density_distance'] < 1.e-12
    assert result['rohf_somo_active_projection_eigenvalues'] == pytest.approx([1,1,1])
    # State-average rotations and signs must not look like orbital flips.
    u, _ = np.linalg.qr(rng.normal(size=(4,4)))
    a['cas'] = hf.copy(); a['cas'][:,1:5] = hf[:,1:5]@u
    rotated = compare(a,electrons=5,active_electrons=3,active_orbitals=4,multiplicity=4)
    assert rotated['rohf_somo_active_projection_eigenvalues'] == pytest.approx([1,1,1])
    paths = {}
    for key, value in a.items():
        p = tmp_path/f'job.identity_{key}.out'
        p.write_text('# MATROP\nBEGIN_DATA\n'+'\n'.join(' '.join(f'{v:.16E}' for v in row) for row in value)+'\nEND_DATA\n')
        assert matrix(p) == pytest.approx(value)
        paths[p.name] = file_sha256(p)
    atomic_json(tmp_path/'execution.json',dict(reference_identity_matrices=paths))
    assert report(tmp_path,'job',dict(active_electrons=3,active_orbitals=4),5,4)['closest_density_root'] == 2
    p.write_text('changed')
    with pytest.raises(ValueError,match='changed'): report(tmp_path,'job',dict(active_electrons=3,active_orbitals=4),5,4)
    a['hf'] = a['hf'].T
    with pytest.raises(ValueError,match='orthonormality'): compare(a,electrons=5,active_electrons=3,active_orbitals=4,multiplicity=4)
