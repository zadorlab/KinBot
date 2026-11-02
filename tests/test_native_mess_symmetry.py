"""Native R8 regressions are evidence, never production molecule rules."""
import copy
import hashlib
import json
from pathlib import Path
import re
import numpy as np
import pytest
from kinbot.stationary_pt import StationaryPoint
from kinbot import symmetry
from kinbot.optical import evaluate_optical
from kinbot.gaussian_modes import read_normal_modes
from kinbot.transition_symmetry import mirror_mode_parity, annotate_surface_symmetry
from kinbot.mess_mirrors import complete_mirror_channels
from kinbot.mess_networks import split_model, _HEADER

DATA = json.loads((Path(__file__).parent/'reference/mess_symmetry_review.json').read_text())
LOG = Path(__file__).parent/'reference/mirror_even_transition.log'


def point(name, *, smiles=''):
    entry = DATA['states'][name]
    p = StationaryPoint(name, 0, 1 if name.startswith('pr') else 2, smiles=smiles,
        atom=entry['atoms'], geom=np.array(entry['geometry']), wellorts=int(name.startswith('rxn')))
    p.characterize()
    symmetry.calculate_symmetry(p)
    return p


def test_bound_assembly_does_not_inherit_free_fragment_rotations():
    p = point('w_10')
    assert p.rigid_rotation_evidence['previous_graph_sigma'] == 2
    assert p.sigma_ext == 1
    assert evaluate_optical(p)['remaining_multiplier'] == 2


def test_asymmetric_saddle_does_not_inherit_endpoint_arm_exchange():
    p = point('rxn_29')
    assert p.rigid_rotation_evidence['previous_graph_sigma'] == 2
    assert p.sigma_ext == 1
    # A TS has no independent stable-species stereochemical population.
    # Its geometric mirror is distinct; endpoint scope is set by the writer.
    from kinbot.molecular_symmetry import geometric_mirror_states
    assert geometric_mirror_states(p) == 2


@pytest.mark.parametrize('name', ['w_14', 'w_15'])
def test_twofold_chiral_minima_keep_sigma_two_and_mirror_weight_two(name):
    from kinbot.molecular_symmetry import rigid_rotation_maps
    p = point(name)
    assert len(rigid_rotation_maps(p)) == p.sigma_ext == 2
    assert evaluate_optical(p)['remaining_multiplier'] == 2


def test_saved_full_rotor_covers_mirror_without_extra_factor():
    p = point('w_6')
    rows = sorted(DATA['w6_scan'], key=lambda r: r['name'])
    assert len(rows) == 12 and {r['status'] for r in rows} == {'normal'}
    rotor = dict(index=0, axis=[0, 1], dihedral=p.dihed[0], sigma_int=1,
        represented_domain_degrees=[0, 360], usable=True,
        points=[dict(index=i, angle_offset_degrees=30*i, status='successful',
            geometry_angstrom=r['geometry'], electronic_energy_hartree=r['energy_hartree'])
            for i, r in enumerate(rows)])
    answer = evaluate_optical(p, rotors=[rotor])
    assert answer['status'] == 'resolved'
    assert answer['states_covered_by_hir'] == 2
    assert answer['remaining_multiplier'] == 1


def test_native_cyclic_product_cannot_keep_open_chain_display_smiles():
    p = point('pr_6', smiles=DATA['states']['pr_11']['input_header'].split('!')[-1].split()[-1])
    q = point('pr_11')
    assert '1' in p.smiles and '1' not in q.smiles
    assert p.connectivity_smiles_history[-1]['current'] == p.smiles
    from kinbot.stereo_identity import canonical_identity
    assert canonical_identity(p)['mirror_id'] == canonical_identity(point('pr_18'))['id']
    assert canonical_identity(q)['id'] not in (canonical_identity(p)['id'], canonical_identity(p)['mirror_id'])


def test_native_mode_fit_and_mirror_parity_are_rotation_invariant():
    p = point('rxn_20')
    assert hashlib.sha256(LOG.read_bytes()).hexdigest() == DATA['rxn20_modes']['sha256']
    native = read_normal_modes(LOG, p.atom, p.geom)
    i = np.argmin(native['frequencies_cm1'])
    assert native['frequencies_cm1'][i] == -850.9602
    answer = mirror_mode_parity(p, native['modes'][i])
    assert answer['status'] == 'self' and answer['parity'] == 'even'
    assert min(m['reactive_mode_overlap'] for m in answer['mappings']) > .9999
    rot = np.linalg.qr(np.array([[1., 2., 3.], [4., 3., 1.], [2., 5., 3.]]))[0]
    p.geom = p.geom @ rot + 12.
    rotated = read_normal_modes(LOG, p.atom, p.geom)
    assert mirror_mode_parity(p, rotated['modes'][i])['parity'] == 'even'
    p.geom[0, 0] += .1
    with pytest.raises(ValueError, match='geometry differs'):
        read_normal_modes(LOG, p.atom, p.geom)


@pytest.mark.parametrize('name,expected', [('w_10', .5), ('rxn_29', .5), ('w_14', 1.)])
def test_legacy_rigid_fragment_upgrade_keeps_optical_division_and_is_idempotent(name, expected):
    from kinbot.mess_symmetry import upgrade_cached_rigid_rotor, mark_current
    p = point(name)
    original = ('RRHO\nGeometry[angstrom] '+str(p.natom)+'\n'+''.join(
        str(a)+' '+' '.join(map(str, g))+'\n' for a, g in zip(p.atom, p.geom))+
        'Core RigidRotor\nSymmetryFactor '+str(DATA['states'][name]['mess_symmetry_factor'])+
        '\nEnd\nEnd ! RRHO\n')
    def upgrade(text):
        return upgrade_cached_rigid_rotor(text, 0, 2, saddle=name.startswith('rxn'))
    result = upgrade(original)
    assert float(re.search(r'SymmetryFactor\s+(\S+)', result)[1]) == expected
    assert upgrade(result) == result
    assert upgrade(mark_current(original)) == mark_current(original)
    assert upgrade(original.replace('End ! RRHO', 'Rotor Hindered\nEnd\nEnd ! RRHO')) == original.replace(
        'End ! RRHO', 'Rotor Hindered\nEnd\nEnd ! RRHO')


def test_self_mirror_harmonic_surface_shares_flux_and_replays_once(tmp_path):
    p = point('rxn_20')
    native = tmp_path/'selected.log'; native.write_bytes(LOG.read_bytes())
    p.source_job = str(native.with_suffix(''))
    populations = ('Model\n Well a\n! kinbot_population A A\n End\n'
                   ' Well r\n! kinbot_population R S\n End\n'
                   ' Well s\n! kinbot_population S R\n End\n')
    barrier = ('! kinbot_stereopath ordinary\n Barrier {name} a r\n'
               '! kinbot_mirror_endpoints A A R S\n RRHO\n'
               ' Core RigidRotor\n SymmetryFactor 1.0\n End\n'
               ' ZeroEnergy[kcal/mol] 10\n End ! RRHO\n End\n')
    annotated = annotate_surface_symmetry(barrier, p).format(name='rxn')
    model = complete_mirror_channels(populations+annotated+'End ! end kinetics\n')
    blocks = [b for b in split_model(model)[1] if _HEADER.search(b)[1] == 'Barrier']
    assert len(blocks) == 2
    assert sum(1/float(re.search(r'SymmetryFactor\s+(\S+)', b)[1]) for b in blocks) == 1.
    assert all('kinbot_shared_mirror_surface_weight 0.5' in b for b in blocks)
    assert complete_mirror_channels(model) == model


def test_native_modes_reject_incomplete_frequency_output(tmp_path):
    p = point('rxn_20'); target = tmp_path/'incomplete.log'
    target.write_text(LOG.read_text().replace('Normal termination of Gaussian', 'Error termination of Gaussian'))
    with pytest.raises(ValueError, match='normally completed'):
        read_normal_modes(target, p.atom, p.geom)


def test_real_writer_and_cached_template_use_selected_native_modes(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot.mess import MESS
    from kinbot.parameters import Parameters
    from kinbot.species_routing import routing_key
    from kinbot.transition_symmetry import annotate_cached_surface, SURFACE_MARKER
    monkeypatch.chdir(tmp_path)
    native = tmp_path/'selected_high.log'; native.write_bytes(LOG.read_bytes())
    p, q, ts = point('w_8'), point('w_3'), point('rxn_20')
    ts.source_job = str(native.with_suffix(''))
    ts.freq = ts.reduced_freqs = DATA['rxn20_modes']['frequencies_cm1']
    ts.stereopath_id = 'ordinary'
    ts.stereopath_metadata = dict(id='ordinary', schema='kinbot.stereopath.v1')
    for x in (p, q, ts):
        x.energy = -836.; x.zpe = .03
    p.reac_type = ['R_Addition_MultipleBond']
    (tmp_path/'input.json').write_text(json.dumps(dict(barrier_threshold=150., me=0)))
    par = dict(Parameters('input.json', show_warnings=False).par, pes=0, multi_conf_tst=0,
               rotor_scan=0, high_level=1, charge=0, mult=2)
    writer = MESS(par, p)
    writer.well_names = {routing_key(p): 'r', routing_key(q): 'a'}
    writer.ts_names = {'selected': 'rxn'}
    reaction = SimpleNamespace(ts=ts, products=[q], instance_name='selected', do_vdW=False)
    rendered = writer.write_barrier(reaction, 0, 30., 20., 0., 1., 1., 0)[0]
    assert SURFACE_MARKER.search(rendered)
    legacy = SURFACE_MARKER.sub('', rendered)
    assert annotate_cached_surface(legacy, tmp_path, 'selected', dict(par, high_level=0)) == legacy
    upgraded = annotate_cached_surface(legacy, tmp_path, 'selected', par)
    assert SURFACE_MARKER.search(upgraded)
    assert annotate_cached_surface(upgraded, tmp_path, 'selected', par) == upgraded
    assert native.read_bytes() == LOG.read_bytes()


def test_unknown_native_mode_cannot_license_flux_sharing():
    from kinbot.transition_symmetry import selected_surface_symmetry
    from kinbot.calculation import array_fingerprint
    p = point('rxn_20')
    p.hess = np.eye(p.geom.size)
    p.hess[0, 0] = -1.
    p.optical_hessian_reference = dict(geometry_sha256=array_fingerprint(p.geom),
        hessian_sha256=array_fingerprint(p.hess), hessian_massweighted=False)
    # Geometry/hash coincidence alone is not a selected calculation reference.
    assert selected_surface_symmetry(p)['status'] == 'unavailable'


def test_each_mc_surface_is_weighted_individually_and_only_with_a_mirror(tmp_path):
    from kinbot.mess_mirrors import _share_self_mirror_surfaces
    from kinbot.transition_symmetry import SURFACE_MARKER
    native = tmp_path/'native.log'; native.write_bytes(LOG.read_bytes())
    p = point('rxn_20'); p.source_job = str(native.with_suffix(''))
    rrho = ' RRHO\n Core RigidRotor\n SymmetryFactor 1\n End\n ZeroEnergy[kcal/mol] 10\n End ! RRHO\n'
    shared = annotate_surface_symmetry(rrho, p)
    first = ' Barrier first a r\n! kinbot_mirror_endpoints A A R S\n Union\n'+shared+rrho+' End\n'
    second = first.replace('first a r', 'second a s').replace('A A R S', 'A A S R')
    assert _share_self_mirror_surfaces([first]) == [first]
    changed = _share_self_mirror_surfaces([first, second])
    assert all(re.findall(r'SymmetryFactor\s+(\S+)', b) == ['2', '1'] for b in changed)
    assert _share_self_mirror_surfaces(changed) == changed
    # A mirror-odd coordinate reverses direction. It does not license this split.
    import base64
    evidence = json.loads(base64.urlsafe_b64decode(SURFACE_MARKER.search(shared)[1]))
    evidence['parity'] = 'odd'
    odd = SURFACE_MARKER.sub('! kinbot_transition_symmetry '+base64.urlsafe_b64encode(
        json.dumps(evidence).encode()).decode(), shared)
    original = [b.replace(shared, odd) for b in (first, second)]
    assert _share_self_mirror_surfaces(original) == original
