"""Read-only checks of complete bimolecular calculation and counting exports."""
import math
import sqlite3
import json

import numpy as np
from ase.db import connect
from ase.io import read

from kinbot import constants, frequencies
from kinbot.calculation import array_fingerprint


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def check_thermochemistry(data, label, errors, warnings, database=None, *, is_ts=False,
                         imagfreq_threshold=50.):
    def fail(message):
        errors.append(f'{label}: {message}')
    atoms = data.get('atoms', [])
    raw = data.get('raw_harmonic_frequencies_cm-1', [])
    projected = data.get('thermochemical_frequencies_cm-1', [])
    if not atoms:
        fail('missing atoms')
    expected = {0} if len(atoms) == 1 else {1} if len(atoms) == 2 else {3 * len(atoms) - 6, 3 * len(atoms) - 5}
    if len(raw) not in expected or not all(finite(v) for v in raw):
        fail('raw frequency count or values are invalid')
    if not all(finite(v) for v in projected):
        fail('projected frequencies contain invalid values')
    for field in ('electronic_energy_hartree', 'zpe_hartree', 'sigma_ext'):
        if not finite(data.get(field)):
            fail(f'{field} is missing or nonfinite')
    if data.get('conformer_counting_error'):
        fail('conformer counting error: ' + str(data['conformer_counting_error']))
    optical = data.get('optical_counting', {})
    if optical.get('status') == 'per_member':
        members = optical.get('members', [])
        if not members or optical.get('remaining_multiplier') is not None:
            fail('MC counting lacks individual weights or adds an ensemble multiplier')
        inventory = {item['member_id']: item for item in data.get('conformer_inventory', [])}
        retained = {key for key, item in inventory.items() if item.get('retained')}
        if {item['member_id'] for item in members} != retained:
            fail('MC weights do not correspond to the retained conformers')
        for member in members:
            key = member['member_id']
            if member.get('remaining_optical_weight') not in (1., 2.):
                fail(f'{key}: missing optical weight')
            if member.get('mirror_partner_ids') and member.get('remaining_optical_weight') != 1.:
                fail(f'{key}: explicit mirror receives an additional optical multiplier')
            record = inventory.get(key, {})
            if not finite(record.get('zero_energy_hartree')):
                fail(f'{key}: missing total energy')
            if len(record.get('frequencies_cm1', [])) not in expected:
                fail(f'{key}: incomplete harmonic frequencies')
    elif (optical.get('status') not in ('resolved', 'assumed')
          or optical.get('remaining_multiplier') not in (1., 2.)):
        fail('optical counting is unresolved or has no usable weight')
    if optical.get('status') == 'assumed':
        warnings.append(f'{label}: optical weight uses {optical.get("heuristic", "an explicit approximation")}')
    if optical.get('measured_coverage') and optical.get('remaining_multiplier') != 1.:
        fail('measured HIR mirror coverage receives an additional optical multiplier')
    if len(atoms) == 1 and optical.get('remaining_multiplier') != 1.:
        fail('an atom must have one optical state')
    if optical.get('site_degeneracy_included'):
        fail('site degeneracy has been included in the optical weight')
    rotors = data.get('hir', {}).get('rotors', [])
    usable = [rotor for rotor in rotors if rotor.get('usable')]
    projection = data.get('frequency_projection') or {}
    if usable:
        if optical.get('status') == 'per_member':
            fail('MC-RRHO and HIR are both active')
        if not projection.get('reference'):
            fail('usable rotors lack a documented Hessian/projection reference')
        if len(raw) - len(projected) != projection.get('internal_rank'):
            fail('rotor projection and exported frequency counts disagree')
        projected_indices = {item['rotor_index'] for item in projection.get('rotors', []) if item.get('projected')}
        if {rotor['index'] for rotor in usable} != projected_indices:
            fail('projected motions and usable rotor contributions disagree')
        for rotor in usable:
            points = rotor.get('points', [])
            success = [point for point in points if point.get('status') == 'successful']
            if not success or not all(finite(point.get('relative_energy_kcal_mol')) for point in success):
                fail(f'rotor {rotor["index"]}: missing measured energies')
            if len(success) != len(points):
                warnings.append(f'{label}: rotor {rotor["index"]} uses a partial scan ({len(success)}/{len(points)})')
    elif len(raw) != len(projected) or (not projection and not np.allclose(
            frequencies.thermochemical_frequencies(raw, int(is_ts), imagfreq_threshold),
            projected, rtol=0., atol=1.e-8)):
        fail('harmonic treatment changed frequencies without a usable rotor')
    if projection:
        reference = projection.get('reference') or {}
        if (reference.get('geometry_sha256') != array_fingerprint(data.get('geometry_angstrom', []))
                or reference.get('source_job') != data.get('source_job')
                or reference.get('source_row_id') != data.get('source_row_id')):
            fail('projection reference does not match the selected geometry/job')
    if database is not None:
        try:
            row = database.get(id=data['source_row_id'])
            if row.name != data['source_job'] or row.data.get('status') != 'normal':
                fail('selected source job/status does not match the exported record')
            comparisons = (
                (row.positions, data['geometry_angstrom'], 'geometry'),
                (row.data['energy'] * constants.EVtoHARTREE, data['electronic_energy_hartree'], 'energy'),
                (row.data['zpe'], data['zpe_hartree'], 'ZPE'),
                (row.data['frequencies'], raw, 'raw frequencies'))
            for original, exported, name in comparisons:
                if not np.allclose(original, exported, rtol=0., atol=1.e-8):
                    fail(f'{name} differs from the selected database observation')
        except (KeyError, ValueError, TypeError) as error:
            fail(f'cannot check selected database observation: {error}')


def validate_network_export(directory, network):
    errors, warnings = [], []
    input_path = directory / 'input.json'
    parameters = json.loads(input_path.read_text()) if input_path.is_file() else {}
    threshold = parameters.get('imagfreq_threshold', 50.)
    db_path = directory / 'kinbot.db'
    database = None
    if db_path.is_file():
        database = connect(str(db_path), use_lock_file=False)
        # Enforce read-only access even when inspecting a live benchmark.
        database._connect = lambda: sqlite3.connect(db_path.resolve().as_uri() + '?mode=ro', uri=True)
    else:
        errors.append('kinbot.db is absent; selected records cannot be verified')
    if network.get('rate_model_ready') is not False:
        errors.append('export incorrectly claims a completed rate model')
    for node in network.get('nodes', []):
        species = ([node['stationary_point']] if 'stationary_point' in node else node.get('fragment_data', []))
        for index, point in enumerate(species):
            check_thermochemistry(point.get('thermochemistry', {}), f'{node["id"]}/fragment{index}',
                                 errors, warnings, database, imagfreq_threshold=threshold)
    for ts in network.get('transitions', []):
        label = ts['id']
        data = ts.get('stationary_point', {}).get('thermochemistry', {})
        check_thermochemistry(data, label, errors, warnings, database,
                             is_ts=True, imagfreq_threshold=threshold)
        if not ts.get('product_matches_expected'):
            errors.append(f'{label}: final IRC does not give the intended abstraction products')
        site = ts.get('site_counting', {})
        if site.get('unlabelled_model_conversion', {}).get('additional_site_multiplier') != 1.:
            errors.append(f'{label}: residual site counting convention is missing')
        for direction in ('forward', 'reverse'):
            path = directory / 'bimol/unique_ts' / label / f'irc_{direction}.xyz'
            try:
                frames = read(path, index=':')
                if not frames or any(frame.get_chemical_symbols() != data['atoms']
                                     or not np.all(np.isfinite(frame.positions)) for frame in frames):
                    errors.append(f'{label}: {direction} IRC trajectory is invalid')
            except (OSError, ValueError, IndexError):
                errors.append(f'{label}: {direction} IRC trajectory is missing or unreadable')
    return errors, warnings
