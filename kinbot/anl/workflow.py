"""Read method-aware components from verified dispatcher task records."""

from __future__ import annotations

import json
import hashlib
import math
import re
from dataclasses import replace

from ase.data import atomic_numbers

from kinbot.anl.dispatch import _load, _verify_execution, _verify_stage_files
from kinbot.anl.extrapolation import (core_valence_correction, two_point_cbs)
from kinbot.anl.model import (ComponentRequirement, ComponentResult,
                              IncompleteRecipeError)
from kinbot.anl.results import parse_result, parsed_result_extends


def _verified_task_result(run_dir, task_id):
    """Reparse one complete task after verifying every recorded artifact."""
    run_dir, spec, state = _load(run_dir)
    tasks = {task['id']: task for task in spec['tasks']}
    if task_id not in tasks:
        raise KeyError(task_id)
    task = tasks[task_id]
    entry = state['tasks'].get(task_id)
    if entry is None or entry['status'] != 'complete':
        raise IncompleteRecipeError(f'{task_id}: task is not complete.')
    request = task.get('result_parser')
    if request is None:
        raise IncompleteRecipeError(f'{task_id}: no method-aware parser was declared.')
    _verify_stage_files(run_dir, task, entry)
    directory = run_dir / 'tasks' / task_id
    execution = json.loads((directory / 'execution.json').read_text())
    _verify_execution(run_dir, task, entry, execution)
    if execution['status'] != 'executed':
        raise IncompleteRecipeError(f'{task_id}: QC execution failed.')
    native = directory / request['file']
    parsed = parse_result(native.read_text(errors='replace'), request)
    if not parsed_result_extends(
            execution.get('details', {}).get('parsed_result'), parsed):
        raise ValueError(f'{task_id}: saved parsed result differs from native output.')
    return run_dir, spec, task, execution, native, parsed


def task_component(run_dir, task_id, *, key, state_id):
    """Return one validated native QC component from a completed task.

    The source is reparsed after checking all staged and output artifact hashes.
    Derived core-valence, relativistic, and higher-order providers are
    still separate work; a dispatch success alone cannot create those terms.
    """
    run_dir, spec, task, execution, native, parsed = _verified_task_result(
        run_dir, task_id)
    request = task['result_parser']
    kind = parsed['kind']
    settings = {}
    review_required = parsed.get('review_required', False)
    if kind in ('molpro_energy', 'mrcc_energy', 'cfour_energy'):
        quantity = 'electronic'
        value = parsed['energy_hartree']
        method = parsed['method']
        basis = parsed['basis']
        if 'reference' in parsed:
            settings['reference'] = parsed['reference']
        if 'program_variant' in parsed:
            settings['program_variant'] = parsed['program_variant']
        if 'effective_program_variant' in parsed:
            settings['effective_program_variant'] = \
                parsed['effective_program_variant']
        if 'rank_exact' in parsed:
            settings['rank_exact'] = parsed['rank_exact']
        if 'correlation' in parsed:
            settings['correlation'] = parsed['correlation']
        if kind == 'molpro_energy':
            if 'f12_auxiliary' in parsed:
                settings['f12_auxiliary'] = parsed['f12_auxiliary']
            if 'core' in parsed:
                settings['core'] = parsed['core']
            if 'relativistic' in parsed:
                settings['relativistic'] = parsed['relativistic']
        if kind in ('mrcc_energy', 'cfour_energy'):
            settings.update(reference=parsed['reference'], core=parsed['core'],
                            driver=parsed['driver'], program=parsed['program'])
        if method == 'CCSD(T)-F12b':
            settings['scale_trip'] = 1
        if (settings.get('reference') in ('RHF', 'ROHF')
                and settings.get('correlation') == 'unrestricted'):
            settings['reference_convention'] = 'RHF/ROHF-UCC'
    elif kind in ('molpro_harmonic', 'gaussian_harmonic'):
        quantity = 'zpe'
        value = parsed['zpe']['hartree']
        method = parsed['method']
        basis = parsed['basis']
        if kind == 'gaussian_harmonic':
            settings['dispersion'] = parsed['dispersion']
            if parsed['dispersion'] == 'GD3BJ':
                method += '-D3BJ'
        if method == 'CCSD(T)-F12b':
            settings['scale_trip'] = 1
        if 'reference' in parsed:
            settings['reference'] = parsed['reference']
        if 'program_variant' in parsed:
            settings['program_variant'] = parsed['program_variant']
        if 'effective_program_variant' in parsed:
            settings['effective_program_variant'] = \
                parsed['effective_program_variant']
        if 'rank_exact' in parsed:
            settings['rank_exact'] = parsed['rank_exact']
        if 'correlation' in parsed:
            settings['correlation'] = parsed['correlation']
        if (settings.get('reference') in ('RHF', 'ROHF')
                and settings.get('correlation') == 'unrestricted'):
            settings['reference_convention'] = 'RHF/ROHF-UCC'
    elif kind == 'gaussian_vpt2':
        if parsed['optimized_in_job']:
            raise IncompleteRecipeError(
                f'{task_id}: VPT2 optimized in the frequency job; '
                'use a separate accepted L2 geometry.')
        quantity = 'correction'
        value = parsed['anharmonic_correction_hartree']
        method = (parsed['method'] + '-D3BJ'
                  if parsed['dispersion'] == 'GD3BJ' else parsed['method'])
        basis = parsed['basis']
        settings['dispersion'] = parsed['dispersion']
        if parsed.get('framework_group_cap'):
            settings['framework_group_cap'] = parsed['framework_group_cap']
    elif kind == 'cfour_dboc':
        if 'basis' not in request:
            raise IncompleteRecipeError(
                f'{task_id}: CFOUR DBOC task lacks a basis-specific parser.')
        quantity = 'correction'
        value = parsed['selected']['hartree']
        method = parsed['selected_level']
        basis = request['basis']
        if 'reference' in parsed:
            settings['reference'] = parsed['reference']
    else:
        raise ValueError(f'{task_id}: unsupported component parser {kind!r}.')
    return ComponentResult(
        key=key, value_hartree=value, quantity=quantity, method=method,
        basis=basis, backend=task['backend'].lower(), state_id=state_id,
        charge=spec['molecule'].get('charge', 0),
        multiplicity=spec['molecule'].get('multiplicity', 1),
        geometry_sha256=execution['geometry_sha256'],
        source_sha256=execution['artifacts'][request['file']],
        source=str(native), settings=settings, review_required=review_required,
    )


def rank_exact_higher_order_component(
        molecule, low: ComponentResult, requirement: ComponentRequirement,
        *, key: str) -> ComponentResult:
    """Represent a rigorously zero higher-excitation correction.

    A state containing at most two electrons cannot form triple or quadruple
    excitations.  Its CCSD(T), CCSDT(Q), and CCSDTQ(P) total energies are
    therefore identical within a fixed one-particle basis.  This provider is
    deliberately restricted to that provable case; an arbitrary failed QC
    calculation can never be converted into a zero correction.
    """
    if not isinstance(molecule, dict):
        raise TypeError('molecule must be an object.')
    symbols = molecule.get('symbols')
    charge = molecule.get('charge', 0)
    if (not isinstance(symbols, list) or not symbols
            or isinstance(charge, bool) or not isinstance(charge, int)):
        raise ValueError('Rank-exact component needs symbols and integer charge.')
    try:
        electrons = sum(atomic_numbers[symbol] for symbol in symbols) - charge
    except KeyError as exc:
        raise ValueError(f'Unknown element symbol {exc.args[0]!r}.') from exc
    if electrons < 1 or electrons > 2:
        raise ValueError('Higher-order rank exactness requires at most two electrons.')
    expected_low = {
        'CCSDT(Q)': 'CCSD(T)',
        'CCSDTQ(P)': 'CCSDT(Q)',
    }.get(requirement.method)
    if (low.quantity != 'electronic' or low.method != expected_low
            or low.basis != requirement.basis
            or low.charge != charge
            or low.multiplicity != molecule.get('multiplicity', 1)
            or requirement.method not in ('CCSDT(Q)', 'CCSDTQ(P)')):
        raise ValueError('Low component is incompatible with rank exactness.')
    settings = dict(requirement.settings)
    settings['rank_exact'] = {
        'electron_count': electrons,
        'reason': ('triple and higher excitations are absent for a state '
                   'containing at most two electrons'),
        'reused_low_method': low.method,
        'reused_low_source_sha256': low.source_sha256,
    }
    provenance = {
        'schema': 1, 'provider': 'known_zero_excitation_rank',
        'electron_count': electrons, 'high_method': requirement.method,
        'low_method': low.method, 'basis': requirement.basis,
        'low_source_sha256': low.source_sha256,
        'geometry_sha256': low.geometry_sha256,
    }
    digest = hashlib.sha256(json.dumps(
        provenance, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return ComponentResult(
        key=key, value_hartree=low.value_hartree,
        quantity=requirement.quantity, method=requirement.method,
        basis=requirement.basis, backend='known_zero',
        state_id=low.state_id, charge=low.charge,
        multiplicity=low.multiplicity,
        geometry_sha256=low.geometry_sha256, source_sha256=digest,
        source=('rank-exact higher-order identity derived from '
                f'{low.source}'), settings=settings)


def rank_exact_higher_order_correction(
        molecule, requirement: ComponentRequirement, *, state_id: str,
        geometry_sha256: str) -> ComponentResult:
    """Return the exact-zero post-CCSD increment without scheduling QC.

    One- and two-electron states cannot contain triple or higher excitation
    amplitudes.  The ANL higher-order difference is therefore zero by rank,
    and neither side of the otherwise cancelling subtraction needs to run.
    """
    if not isinstance(molecule, dict):
        raise TypeError('molecule must be an object.')
    symbols = molecule.get('symbols')
    charge = molecule.get('charge', 0)
    if (not isinstance(symbols, list) or not symbols
            or isinstance(charge, bool) or not isinstance(charge, int)):
        raise ValueError('Rank-exact correction needs symbols and integer charge.')
    try:
        electrons = sum(atomic_numbers[symbol] for symbol in symbols) - charge
    except KeyError as exc:
        raise ValueError(f'Unknown element symbol {exc.args[0]!r}.') from exc
    if electrons < 1 or electrons > 2:
        raise ValueError('Higher-order rank exactness requires at most two electrons.')
    if (requirement.key != 'higher_order_correction'
            or requirement.quantity != 'correction'
            or requirement.backends != ('known_zero',)
            or requirement.settings.get('rank_exact_electron_maximum') != 2):
        raise ValueError('Invalid rank-exact higher-order requirement.')
    if (not isinstance(geometry_sha256, str)
            or not re.fullmatch(r'[0-9a-f]{64}', geometry_sha256)):
        raise ValueError('Rank-exact correction needs a geometry hash.')
    provenance = {
        'schema': 1, 'provider': 'known_zero_excitation_rank',
        'symbols': symbols, 'charge': charge,
        'multiplicity': molecule.get('multiplicity', 1),
        'electron_count': electrons, 'component': requirement.key,
        'method': requirement.method, 'basis': requirement.basis,
        'geometry_sha256': geometry_sha256,
    }
    digest = hashlib.sha256(json.dumps(
        provenance, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    settings = dict(requirement.settings)
    settings['rank_exact'] = {
        'electron_count': electrons,
        'reason': ('triple and higher excitations are absent for a state '
                   'containing at most two electrons'),
    }
    return ComponentResult(
        key=requirement.key, value_hartree=0.,
        quantity=requirement.quantity, method=requirement.method,
        basis=requirement.basis, backend='known_zero', state_id=state_id,
        charge=charge, multiplicity=molecule.get('multiplicity', 1),
        geometry_sha256=geometry_sha256, source_sha256=digest,
        source='exact zero: triple and higher excitation ranks are absent',
        settings=settings)


def rank_exact_core_valence_component(
        molecule, requirement: ComponentRequirement, *, state_id: str,
        geometry_sha256: str) -> ComponentResult:
    """Represent the exact-zero core-valence term for H/He-only species."""
    if not isinstance(molecule, dict):
        raise TypeError('molecule must be an object.')
    symbols = molecule.get('symbols')
    if not isinstance(symbols, list) or not symbols:
        raise ValueError('Core-valence exactness needs molecular symbols.')
    try:
        numbers = [atomic_numbers[symbol] for symbol in symbols]
    except KeyError as exc:
        raise ValueError(f'Unknown element symbol {exc.args[0]!r}.') from exc
    if any(number > 2 for number in numbers):
        raise ValueError('Core-valence exactness requires only H/He atoms.')
    if (requirement.key != 'core_valence_cbs'
            or requirement.quantity != 'correction'):
        raise ValueError('Invalid core-valence recipe requirement.')
    if (not isinstance(geometry_sha256, str)
            or not re.fullmatch(r'[0-9a-f]{64}', geometry_sha256)):
        raise ValueError('Core-valence exactness needs a geometry hash.')
    provenance = {
        'schema': 1, 'provider': 'known_zero_no_frozen_core',
        'symbols': symbols, 'charge': molecule.get('charge', 0),
        'multiplicity': molecule.get('multiplicity', 1),
        'component': requirement.key, 'method': requirement.method,
        'basis': requirement.basis, 'geometry_sha256': geometry_sha256,
    }
    digest = hashlib.sha256(json.dumps(
        provenance, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    settings = dict(requirement.settings)
    settings['rank_exact'] = {
        'reason': ('H/He-only species have no frozen inner-shell orbitals; '
                   'all-electron and frozen-core CBS energies are identical'),
        'atomic_numbers': numbers,
    }
    settings['basis_uncertainty_series'] = {
        requirement.key: {
            'lower_hartree': 0., 'upper_hartree': 0.,
            'lower_basis': requirement.settings.get('native_basis_pair', ('cc-pCVTZ', 'cc-pCVQZ'))[0],
            'upper_basis': requirement.settings.get('native_basis_pair', ('cc-pCVTZ', 'cc-pCVQZ'))[1],
        },
    }
    return ComponentResult(
        key=requirement.key, value_hartree=0.,
        quantity=requirement.quantity, method=requirement.method,
        basis=requirement.basis, backend='known_zero', state_id=state_id,
        charge=molecule.get('charge', 0),
        multiplicity=molecule.get('multiplicity', 1),
        geometry_sha256=geometry_sha256, source_sha256=digest,
        source='exact zero: no frozen inner-shell orbitals',
        settings=settings)


def atomic_zero_vibrational_component(
        molecule, requirement: ComponentRequirement, *, state_id: str,
        geometry_sha256: str) -> ComponentResult:
    """Provide a rigorously zero harmonic or anharmonic atomic term."""
    if (not isinstance(molecule, dict)
            or len(molecule.get('symbols', ())) != 1):
        raise ValueError('Zero vibrational components require one atom.')
    if requirement.key not in ('harmonic_zpe', 'vpt2_correction'):
        raise ValueError('Only atomic vibrational terms are identically zero.')
    if (not isinstance(geometry_sha256, str)
            or not re.fullmatch(r'[0-9a-f]{64}', geometry_sha256)):
        raise ValueError('Atomic vibrational term needs a geometry hash.')
    provenance = {
        'schema': 1, 'provider': 'known_zero_atomic_vibration',
        'symbol': molecule['symbols'][0], 'charge': molecule.get('charge', 0),
        'multiplicity': molecule.get('multiplicity', 1),
        'component': requirement.key, 'method': requirement.method,
        'basis': requirement.basis, 'geometry_sha256': geometry_sha256,
    }
    digest = hashlib.sha256(json.dumps(
        provenance, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    settings = dict(requirement.settings)
    if ('upper_cardinal' in settings
            and 'extrapolation_power' in settings):
        settings['basis_uncertainty_series'] = {
            requirement.key: {
                'lower_hartree': 0., 'upper_hartree': 0.,
                'lower_basis': 'cc-pVTZ', 'upper_basis': 'cc-pVQZ',
            },
        }
    return ComponentResult(
        key=requirement.key, value_hartree=0.,
        quantity=requirement.quantity, method=requirement.method,
        basis=requirement.basis, backend='known_zero', state_id=state_id,
        charge=molecule.get('charge', 0),
        multiplicity=molecule.get('multiplicity', 1),
        geometry_sha256=geometry_sha256, source_sha256=digest,
        source=('exact zero: a monatomic species has no vibrational '
                f'{requirement.key} contribution'),
        settings=settings)


def attach_task_vpt2_frequencies(species, run_dir, task_id, *,
                                 review_file=None, **match_options):
    """Attach verified mode-resolved VPT2 corrections to a KinBot species."""
    from kinbot.energy import attach_anharmonic_frequencies
    from kinbot.anl.validation import _apply_quality_review

    _, spec, task, execution, native, parsed = _verified_task_result(run_dir, task_id)
    if parsed.get('kind') != 'gaussian_vpt2' or task['backend'].lower() != 'gaussian':
        raise ValueError(f'{task_id}: expected a Gaussian VPT2 task.')
    component = task_component(
        run_dir, task_id, key='vpt2_correction',
        state_id='thermochemistry-frequency-handoff')
    component = _apply_quality_review(component, task_id, review_file)
    if component.review_required:
        raise ValueError(
            f'{task_id}: VPT2 fundamentals require explicit quality review.')
    parsed = dict(parsed)
    parsed['review_required'] = False
    from kinbot.anl import anharmonic
    from ase.io import read
    _, _, state = _load(run_dir)
    paired = task_id == 'gaussian_vpt2' and anharmonic.available(spec, state)
    if paired and not (getattr(species, 'anl_thermal_model', None) or {}).get('harmonic'):
        from kinbot.anl.thermal import attach_l2_harmonic_reference
        attach_l2_harmonic_reference(species, parsed, native,
                                    read(native.parent/'geometry.xyz'), run=run_dir)
    if (getattr(species, 'anl_thermal_model', None) or {}).get('harmonic'):
        from kinbot.anl.thermal import attach_vpt2
        frequencies = attach_vpt2(species, parsed, native,
            read(native.parent/'geometry.xyz'), **match_options)
        if paired:
            # The scalar verifier also proves both basis-matched optimized
            # geometries and rejects either unreviewed native VPT2 result.
            anharmonic.component(run_dir, state_id='thermochemistry-frequency-handoff')
            upper_id = spec['intent']['anharmonic_basis_extension']['upper_task']
            _, _, upper_task, upper_execution, upper_native, upper_parsed = _verified_task_result(run_dir, upper_id)
            from kinbot.anl.thermal import attach_vpt2_upper
            frequencies = attach_vpt2_upper(species, upper_parsed, upper_native,
                read(upper_native.parent/'geometry.xyz'),
                power=spec['intent']['anharmonic_basis_extension']['extrapolation_power'],
                evidence=dict(task_id=upper_id, native_output=str(upper_native),
                    source_sha256=upper_execution['artifacts'][upper_task['result_parser']['file']]))
    else:
        frequencies = attach_anharmonic_frequencies(
            species, parsed, **match_options)
    species.anl_thermochemistry_frequency_source.update({
        'task_id': task_id,
        'native_output': str(native),
        'source_sha256': execution['artifacts'][task['result_parser']['file']],
        'quality_review': component.settings.get('quality_review'),
    })
    thermal = getattr(species, 'anl_thermal_model', None)
    if thermal and thermal.get('vpt2'):
        thermal['vpt2']['evidence'] = {
            key: species.anl_thermochemistry_frequency_source[key]
            for key in ('task_id', 'native_output', 'source_sha256', 'quality_review')}
    return frequencies


def core_valence_task_component(
        run_dir, *, all_electron_lower, all_electron_upper,
        frozen_core_lower, frozen_core_upper,
        requirement: ComponentRequirement, state_id: str) -> ComponentResult:
    """Build the all-electron minus frozen-core CCSD(T)/CBS correction."""
    task_ids = (all_electron_lower, all_electron_upper,
                frozen_core_lower, frozen_core_upper)
    parts = [task_component(run_dir, task_id, key=task_id, state_id=state_id)
             for task_id in task_ids]
    ae_lower, ae_upper, fc_lower, fc_upper = parts
    if (requirement.key != 'core_valence_cbs'
            or requirement.quantity != 'correction'
            or 'composite' not in requirement.backends):
        raise ValueError('Invalid core-valence recipe requirement.')
    if len({part.geometry_sha256 for part in parts}) != 1:
        raise ValueError('Core-valence tasks use different geometries.')
    if len({(part.state_id, part.charge, part.multiplicity)
            for part in parts}) != 1:
        raise ValueError('Core-valence tasks use different electronic states.')
    method = requirement.settings.get('component_method', 'CCSD(T)')
    bases = tuple(requirement.settings.get('native_basis_pair', ('cc-pCVTZ', 'cc-pCVQZ')))
    if any(part.method != method or part.quantity != 'electronic'
           or part.backend != 'molpro' for part in parts):
        raise ValueError('Core-valence inputs must be Molpro CCSD(T) energies.')
    if ((ae_lower.basis, ae_upper.basis, fc_lower.basis, fc_upper.basis)
            != bases + bases):
        raise ValueError('Core-valence inputs differ from the recipe basis pair.')
    if [part.settings.get('core') for part in parts] != [
            'all-electron', 'all-electron', 'frozen', 'frozen']:
        raise ValueError('Core-valence inputs have incorrect core treatments.')
    auxiliary = requirement.settings.get('f12_auxiliary_by_basis')
    if auxiliary is not None:
        for part in parts:
            rank = part.settings.get('rank_exact', {})
            # One-electron ions use the native-validated HF reduction. No
            # F12 fitting calculation exists in that deck. Two-electron CCSD
            # still requires the explicitly selected fitting sets.
            hf_exact = (rank.get('electron_count') == 1 and rank.get('effective_method') == 'HF'
                        and part.settings.get('effective_program_variant') == 'RHF')
            expected = None if hf_exact else auxiliary.get(part.basis)
            if part.settings.get('f12_auxiliary') != expected:
                raise ValueError('Core-valence F12 auxiliaries differ from the recipe.')
    comparable = [{key: value for key, value in part.settings.items()
                   if key not in ('core', 'f12_auxiliary')} for part in parts]
    if any(item != comparable[0] for item in comparable[1:]):
        raise ValueError('Core-valence input settings differ.')
    cardinal = requirement.settings.get('upper_cardinal')
    power = requirement.settings.get('extrapolation_power')
    value = core_valence_correction(
        all_electron_lower=ae_lower.value_hartree,
        all_electron_upper=ae_upper.value_hartree,
        frozen_core_lower=fc_lower.value_hartree,
        frozen_core_upper=fc_upper.value_hartree,
        upper_cardinal=cardinal, power=power)
    provenance = {
        'formula': 'CBS(all-electron)-CBS(frozen-core)',
        'task_ids': task_ids,
        'source_sha256': [part.source_sha256 for part in parts],
        'upper_cardinal': cardinal, 'power': power,
    }
    digest = hashlib.sha256(json.dumps(
        provenance, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    settings = dict(requirement.settings)
    settings['basis_uncertainty_series'] = {
        requirement.key: {
            'lower_hartree': float(
                ae_lower.value_hartree - fc_lower.value_hartree),
            'upper_hartree': float(
                ae_upper.value_hartree - fc_upper.value_hartree),
            'lower_basis': bases[0],
            'upper_basis': bases[1],
        },
    }
    first = parts[0]
    return ComponentResult(
        key=requirement.key, value_hartree=value,
        quantity=requirement.quantity, method=requirement.method,
        basis=requirement.basis, backend='composite', state_id=state_id,
        charge=first.charge, multiplicity=first.multiplicity,
        geometry_sha256=first.geometry_sha256, source_sha256=digest,
        source='core-valence(' + ','.join(part.source for part in parts) + ')',
        settings=settings)


def scalar_relativistic_task_component(
        run_dir, relativistic_task, nonrelativistic_task, *,
        requirement: ComponentRequirement, state_id: str) -> ComponentResult:
    """Build the all-electron DKH2 minus nonrelativistic CCSD(T) correction."""
    bases = requirement.settings.get('native_basis_pair')
    if bases is not None:
        if (len(bases) != 2 or not isinstance(relativistic_task, (tuple, list))
                or not isinstance(nonrelativistic_task, (tuple, list))
                or len(relativistic_task) != 2 or len(nonrelativistic_task) != 2
                or len(set((*relativistic_task, *nonrelativistic_task))) != 4):
            raise ValueError('Relativistic CBS needs two distinct matched DKH/NR pairs.')
        settings = {k: v for k, v in requirement.settings.items()
                    if k not in ('native_basis_pair', 'upper_cardinal', 'extrapolation_power')}
        parts = [scalar_relativistic_task_component(
            run_dir, dkh, nr, requirement=replace(requirement, basis=basis,
                                                 settings=settings), state_id=state_id)
            for dkh, nr, basis in zip(relativistic_task, nonrelativistic_task, bases)]
        if len({(p.geometry_sha256, p.state_id, p.charge, p.multiplicity) for p in parts}) != 1:
            raise ValueError('Relativistic CBS inputs use different geometries or states.')
        value = two_point_cbs(parts[0].value_hartree, parts[1].value_hartree,
                             upper_cardinal=requirement.settings['upper_cardinal'],
                             power=requirement.settings['extrapolation_power'])
        provenance = dict(formula='CBS[E(DKH2)-E(nonrelativistic)]',
                          source_sha256=[p.source_sha256 for p in parts],
                          settings=dict(requirement.settings))
        digest = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
        settings = dict(requirement.settings)
        settings['basis_uncertainty_series'] = {requirement.key: dict(
            lower_hartree=parts[0].value_hartree, upper_hartree=parts[1].value_hartree,
            lower_basis=bases[0], upper_basis=bases[1])}
        return replace(parts[1], value_hartree=value, basis=requirement.basis,
                       source_sha256=digest, source='CBS('+','.join(p.source for p in parts)+')',
                       settings=settings)
    relativistic = task_component(
        run_dir, relativistic_task, key=relativistic_task, state_id=state_id)
    nonrelativistic = task_component(
        run_dir, nonrelativistic_task, key=nonrelativistic_task,
        state_id=state_id)
    if (requirement.key != 'scalar_relativistic'
            or requirement.quantity != 'correction'
            or 'composite' not in requirement.backends):
        raise ValueError('Invalid scalar-relativistic recipe requirement.')
    pair = (relativistic, nonrelativistic)
    if (relativistic.geometry_sha256 != nonrelativistic.geometry_sha256
            or (relativistic.state_id, relativistic.charge,
                relativistic.multiplicity) !=
               (nonrelativistic.state_id, nonrelativistic.charge,
                nonrelativistic.multiplicity)):
        raise ValueError('Relativistic tasks use different geometry or state.')
    if any(part.method != 'CCSD(T)' or part.quantity != 'electronic'
           or part.backend != 'molpro' or part.basis != requirement.basis
           or part.settings.get('core') != 'all-electron' for part in pair):
        raise ValueError('Relativistic inputs need all-electron Molpro CCSD(T).')
    if (relativistic.settings.get('relativistic'),
            nonrelativistic.settings.get('relativistic')) != ('DKH2', 'none'):
        raise ValueError('Relativistic inputs do not form DKH2/nonrel pair.')
    comparable = [{key: value for key, value in part.settings.items()
                   if key != 'relativistic'} for part in pair]
    if comparable[0] != comparable[1]:
        raise ValueError('Relativistic input settings differ.')
    value = math.fsum((relativistic.value_hartree,
                       -nonrelativistic.value_hartree))
    provenance = {
        'formula': 'E(DKH2)-E(nonrelativistic)',
        'relativistic_sha256': relativistic.source_sha256,
        'nonrelativistic_sha256': nonrelativistic.source_sha256,
    }
    digest = hashlib.sha256(json.dumps(
        provenance, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return ComponentResult(
        key=requirement.key, value_hartree=value,
        quantity=requirement.quantity, method=requirement.method,
        basis=requirement.basis, backend='composite', state_id=state_id,
        charge=relativistic.charge, multiplicity=relativistic.multiplicity,
        geometry_sha256=relativistic.geometry_sha256, source_sha256=digest,
        source=f'DKH2({relativistic.source})-nonrel({nonrelativistic.source})',
        settings=dict(requirement.settings))


def state_correction_component(*, requirement: ComponentRequirement,
                               value_hartree: float, state_id: str,
                               charge: int, multiplicity: int, source: str,
                               backend='manual') -> ComponentResult:
    """Create an explicit, provenance-bearing state-only correction."""
    if (requirement.geometry_role != 'state'
            or backend not in requirement.backends
            or isinstance(value_hartree, bool)
            or not isinstance(value_hartree, (int, float))
            or not math.isfinite(value_hartree) or not source.strip()):
        raise ValueError('Invalid state-only correction.')
    provenance = {'key': requirement.key, 'value_hartree': value_hartree,
                  'state_id': state_id, 'charge': charge,
                  'multiplicity': multiplicity, 'source': source,
                  'backend': backend}
    digest = hashlib.sha256(json.dumps(
        provenance, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return ComponentResult(
        key=requirement.key, value_hartree=float(value_hartree),
        quantity=requirement.quantity, method=requirement.method,
        basis=requirement.basis, backend=backend, state_id=state_id,
        charge=charge, multiplicity=multiplicity, geometry_sha256=None,
        source_sha256=digest, source=source,
        settings=dict(requirement.settings))


def cbs_task_component(run_dir, lower_task_id, upper_task_id, *,
                       requirement: ComponentRequirement, state_id: str,
                       lower_basis: str, upper_basis: str) -> ComponentResult:
    """Extrapolate two verified electronic or zero-point task results.

    Both native outputs are reparsed and checked against their execution
    records by ``task_component``. Electronic single points must share one
    geometry. Harmonic and VPT2 pairs must each trace to a completed,
    basis-matched optimization, with the higher-basis geometry retained as
    the composite result's geometry. The source digest records both inputs,
    geometry sources, and exact extrapolation parameters.
    """
    if lower_task_id == upper_task_id:
        raise ValueError('A CBS pair needs two distinct tasks.')
    if not isinstance(requirement, ComponentRequirement):
        raise TypeError('CBS requirement must be a ComponentRequirement.')
    lower = task_component(run_dir, lower_task_id, key=lower_task_id,
                           state_id=state_id)
    upper = task_component(run_dir, upper_task_id, key=upper_task_id,
                           state_id=state_id)
    run_dir, spec, state = _load(run_dir)
    tasks = {task['id']: task for task in spec['tasks']}
    if requirement.quantity == 'electronic':
        settings = requirement.settings
        if (not settings.get('geometry_method')
                or not settings.get('geometry_basis')):
            raise ValueError('Electronic CBS recipe lacks a geometry level.')
        if (tasks[lower_task_id].get('geometry_from') !=
                tasks[upper_task_id].get('geometry_from')):
            raise ValueError('Electronic CBS inputs have different geometry sources.')
        geometry_sources = tuple(
            _verified_optimized_geometry(
                run_dir, tasks, state, tasks[task_id], component,
                expected_method=settings['geometry_method'],
                expected_basis=settings['geometry_basis'])
            for task_id, component in ((lower_task_id, lower),
                                       (upper_task_id, upper)))
    else:
        geometry_sources = (
            _verified_optimized_geometry(run_dir, tasks, state,
                                         tasks[lower_task_id], lower),
            _verified_optimized_geometry(run_dir, tasks, state,
                                         tasks[upper_task_id], upper),
        )
    try:
        electron_count = (sum(atomic_numbers[symbol]
                              for symbol in spec['molecule']['symbols'])
                          - spec['molecule'].get('charge', 0))
    except (KeyError, TypeError) as exc:
        raise ValueError('CBS pair lacks a valid molecular electron count.') from exc
    return _cbs_components(
        lower, upper, requirement=requirement,
        lower_basis=lower_basis, upper_basis=upper_basis,
        geometry_sources=geometry_sources, electron_count=electron_count)


def _verified_optimized_geometry(run_dir, tasks, state, task, component, *,
                                 expected_method=None, expected_basis=None):
    source_id = state['tasks'].get(task['id'], {}).get(
        'geometry_from', task.get('geometry_from', 'initial'))
    if source_id == 'initial':
        # Automatic post-L2 graphs bind their initial geometry to the exact
        # accepted KinBot record and its resolved theory.  This is the only
        # valid optimized-geometry substitute for an unavailable L3 point.
        try:
            _, specification, _ = _load(run_dir)
            spec_source = specification['molecule']['source']
            theory = spec_source['l2_theory']
        except (KeyError, TypeError):
            spec_source = None
            theory = None
        method = expected_method or component.method.removesuffix('-D3BJ')
        basis = expected_basis or component.basis
        if (not isinstance(spec_source, dict)
                or not isinstance(theory, dict)
                or not re.fullmatch(
                    r'[0-9a-f]{64}', str(spec_source.get('record_sha256')))
                or theory.get('calculator') != 'gaussian'
                or theory.get('method', '').casefold() != method.casefold()
                or theory.get('basis', '').casefold() != basis.casefold()
                or spec_source.get('geometry_sha256') !=
                   component.geometry_sha256):
            raise ValueError(
                f"{task['id']}: accepted L2 geometry provenance differs.")
        return {
            'task_id': 'initial',
            'geometry_sha256': component.geometry_sha256,
            'record_sha256': spec_source.get('record_sha256'),
        }
    source_task = tasks.get(source_id)
    if source_task is None or source_task.get('kind') != 'ase_optimize':
        raise ValueError(f"{task['id']}: CBS input needs an optimized geometry.")
    profile = source_task.get('profile', {})
    backend = ('gaussian' if component.quantity == 'correction' else 'molpro')
    method = expected_method or (
        component.method.removesuffix('-D3BJ')
        if component.quantity == 'correction' else component.method)
    basis = expected_basis or component.basis
    if (not isinstance(profile, dict)
            or profile.get('calculator', '').lower() not in (
                ('gaussian', 'gauss') if backend == 'gaussian' else ('molpro',))
            or profile.get('method', '').casefold() != method.casefold()
            or profile.get('basis', '').casefold() != basis.casefold()):
        raise ValueError(f"{task['id']}: optimized geometry level differs from CBS input.")
    if component.quantity == 'correction':
        keywords = profile.get('calculator_kwargs', {})
        if (not isinstance(keywords, dict) or
                keywords.get('EmpiricalDispersion', '').upper()
                != component.settings.get('dispersion', '')):
            raise ValueError(f"{task['id']}: optimized geometry dispersion differs.")
    source_entry = state['tasks'].get(source_id)
    if source_entry is None or source_entry.get('status') != 'complete':
        raise IncompleteRecipeError(f'{source_id}: optimized geometry is not complete.')
    _verify_stage_files(run_dir, source_task, source_entry)
    execution = json.loads((run_dir / 'tasks' / source_id /
                            'execution.json').read_text())
    geometry_hash = _verify_execution(run_dir, source_task, source_entry,
                                      execution)
    if (execution['status'] != 'executed' or geometry_hash is None
            or geometry_hash != source_entry.get('final_geometry_sha256')
            or geometry_hash != component.geometry_sha256):
        raise ValueError(f"{task['id']}: accepted optimized geometry differs.")
    return {'task_id': source_id, 'geometry_sha256': geometry_hash,
            'final_xyz_sha256': execution['artifacts']['final.xyz']}


def _one_electron_f12_settings_match(lower, upper, electron_count):
    """Ignore only legacy rank-recovery metadata for one-electron F12 pairs."""
    if (electron_count != 1 or lower.method != 'CCSD(T)-F12b'
            or upper.method != 'CCSD(T)-F12b'):
        return False

    def normalized(component):
        settings = dict(component.settings)
        rank = settings.pop('rank_exact', None)
        effective = settings.pop('effective_program_variant', None)
        if (rank is not None
                and (not isinstance(rank, dict)
                     or rank.get('electron_count') != 1)):
            return None
        # Older normal outputs did not record an effective variant.  The
        # recoverable NaN path records the physically identical UCCSD-F12b
        # one-electron limit.  Bare-HF rank-exact tasks are intentionally not
        # mixed with those archived F12 calculations.
        if effective not in (None, 'RHF-UCCSD-F12b'):
            return None
        return settings

    lower_settings = normalized(lower)
    upper_settings = normalized(upper)
    return (lower_settings is not None
            and lower_settings == upper_settings)


def _cbs_components(lower, upper, *, requirement, lower_basis, upper_basis,
                    geometry_sources=None, electron_count=None):
    allowed = (requirement.quantity == 'electronic'
               or requirement.quantity == 'zpe'
               or (requirement.quantity == 'correction'
                   and requirement.key == 'vpt2_correction'))
    if not allowed or 'composite' not in requirement.backends:
        raise ValueError('Recipe requirement is not a supported CBS term.')
    basis_pair = f'CBS({lower_basis},{upper_basis})'
    native_pair = requirement.settings.get('native_basis_pair')
    if native_pair is not None:
        if (not isinstance(native_pair, (tuple, list))
                or tuple(native_pair) != (lower_basis, upper_basis)):
            raise ValueError('CBS native basis pair disagrees with the recipe.')
    if (not lower_basis or not upper_basis or lower_basis == upper_basis
            or (native_pair is None
                and not (requirement.basis == basis_pair
                         or requirement.basis.startswith(basis_pair + '//')))):
        raise ValueError('CBS basis pair disagrees with the recipe.')
    if (lower.basis, upper.basis) != (lower_basis, upper_basis):
        raise ValueError('CBS native basis order differs from the recipe.')
    if (lower.method != requirement.method or upper.method != requirement.method
            or lower.quantity != requirement.quantity
            or upper.quantity != requirement.quantity):
        raise ValueError('CBS method or quantity differs between inputs.')
    backend = 'gaussian' if requirement.quantity == 'correction' else 'molpro'
    if lower.backend != backend or upper.backend != backend:
        raise ValueError(f'CBS inputs must come from {backend} native parsers.')
    if ((lower.state_id, lower.charge, lower.multiplicity) !=
            (upper.state_id, upper.charge, upper.multiplicity)):
        raise ValueError('CBS inputs have different electronic states.')
    if lower.geometry_sha256 is None or upper.geometry_sha256 is None:
        raise ValueError('CBS inputs lack geometry hashes.')
    if requirement.quantity == 'electronic':
        if lower.geometry_sha256 != upper.geometry_sha256:
            raise ValueError('Electronic CBS inputs have different geometries.')
        if (not isinstance(geometry_sources, tuple)
                or len(geometry_sources) != 2
                or geometry_sources[0]['task_id'] != geometry_sources[1]['task_id']
                or any(source['geometry_sha256'] != lower.geometry_sha256
                       for source in geometry_sources)):
            raise ValueError('Electronic CBS lacks one verified highest-level geometry.')
    else:
        if (requirement.settings.get('geometry_mode') != 'basis_optimized'
                or not isinstance(geometry_sources, tuple)
                or len(geometry_sources) != 2
                or geometry_sources[0]['geometry_sha256'] != lower.geometry_sha256
                or geometry_sources[1]['geometry_sha256'] != upper.geometry_sha256
                or geometry_sources[0]['task_id'] == geometry_sources[1]['task_id']):
            raise ValueError('CBS zero-point inputs lack basis-matched optimized geometries.')
    if lower.review_required or upper.review_required:
        raise IncompleteRecipeError('CBS input requires native quality review.')
    settings_match = lower.settings == upper.settings
    legacy_rank_equivalent = (
        not settings_match
        and _one_electron_f12_settings_match(
            lower, upper, electron_count))
    if not settings_match and not legacy_rank_equivalent:
        raise ValueError('CBS input calculation settings differ.')
    settings = dict(requirement.settings)
    if 'upper_cardinal' not in settings or 'extrapolation_power' not in settings:
        raise ValueError('CBS recipe lacks cardinal number or exponent.')
    for key, value in settings.items():
        if key not in ('extrapolation_power', 'upper_cardinal', 'geometry_mode',
                       'geometry_method', 'geometry_basis',
                       'native_basis_pair') \
                and lower.settings.get(key) != value:
            raise ValueError(f'CBS input {key} setting differs from the recipe.')
    for item in (lower, upper):
        if (not item.source or not isinstance(item.source_sha256, str)
                or re.fullmatch(r'[0-9a-f]{64}', item.source_sha256) is None
                or not math.isfinite(item.value_hartree)):
            raise ValueError('CBS input lacks valid native provenance.')
    cardinal = settings['upper_cardinal']
    power = settings['extrapolation_power']
    value = two_point_cbs(lower.value_hartree, upper.value_hartree,
                          upper_cardinal=cardinal, power=power)
    settings['basis_uncertainty_series'] = {
        requirement.key: {
            'lower_hartree': float(lower.value_hartree),
            'upper_hartree': float(upper.value_hartree),
            'lower_basis': lower_basis,
            'upper_basis': upper_basis,
        },
    }
    provenance = {
        'formula': 'E_n + alpha_n*(E_n-E_(n-1))',
        'lower_sha256': lower.source_sha256,
        'upper_sha256': upper.source_sha256,
        'lower_basis': lower_basis, 'upper_basis': upper_basis,
        'upper_cardinal': cardinal, 'power': power,
        'geometry_sources': geometry_sources,
        'settings_compatibility': (
            'one-electron F12 legacy rank metadata'
            if legacy_rank_equivalent else 'exact'),
    }
    digest = hashlib.sha256(json.dumps(provenance, sort_keys=True,
                                       separators=(',', ':')).encode()).hexdigest()
    return ComponentResult(
        key=requirement.key, value_hartree=value,
        quantity=requirement.quantity, method=requirement.method,
        basis=requirement.basis, backend='composite',
        state_id=lower.state_id, charge=lower.charge,
        multiplicity=lower.multiplicity,
        geometry_sha256=upper.geometry_sha256,
        source_sha256=digest,
        source=f'CBS({lower.source},{upper.source})', settings=settings,
    )
