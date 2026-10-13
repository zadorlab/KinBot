"""Automatic KinBot network to CBH/ANL production planning.

The ordinary KinBot search remains the authority for accepted wells,
products, transition states, frequencies, and hindered rotors.  This module
freezes that live network, expands every requested CBH tree, deduplicates the
required electronic states, and selects the most expensive calculation tier
allowed by one explicit size policy.  Execution and assembly consume this
record; users never enumerate CBH reference species or balance reactions.
"""

from __future__ import annotations

from dataclasses import asdict
from collections import Counter
import hashlib
from itertools import combinations
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np
from ase import Atoms

from kinbot.anl.atct import _canonical_smiles, load_atct_references
from kinbot.anl.cbh import (cbh_network_requirements, cbh_network_states,
                            cbh_size_metrics, FormationEnthalpy,
                            generate_cbh_reaction, select_cbh_ladder,
                            zero_k_from_composite_record)
from kinbot.anl.dispatch import (_file_hash, _geometry_hash, _load, advance,
                                 preflight, prepare, _verify_execution,
                                 _verify_stage_files)
from kinbot.anl.model import ComponentResult, IncompleteRecipeError
from kinbot.anl.recipes import ANL_REDUCIBLE_FAMILIES
from kinbot.anl.validation import (assemble_profiled_anl0,
                                   assemble_profiled_anl0_f12,
                                   automatic_composite_spec)
from kinbot.anl.workflow import _verified_task_result, task_component
from kinbot.species_routing import routing_name
from kinbot.theory import (DEFAULT_ANL_SIZE_POLICY,
                           automatic_cbh_anl_requested, resolve_profiles,
                           validate_anl_size_policy)


# Public compatibility alias used by validation helpers and existing scripts.
DEFAULT_SIZE_POLICY = DEFAULT_ANL_SIZE_POLICY

_TIERS = ('ANL0-F12', 'ANL0', 'L3', 'L3-SP@L2', 'L2')


def _resolve_user_executable(command):
    """Resolve a command without requiring a login-shell PATH mutation.

    User-managed scientific programs are commonly installed into isolated
    prefixes below ``~/.local``.  Search only their conventional ``bin``
    directories and accept the result only when it is unambiguous.
    """
    executable = shutil.which(command)
    if executable is not None:
        return executable
    if os.path.basename(command) != command:
        return None
    local = Path.home() / '.local'
    candidates = [local / 'bin' / command]
    if local.is_dir():
        candidates.extend(sorted(local.glob(f'*/bin/{command}')))
    matches = []
    for candidate in candidates:
        if (candidate.is_file() and os.access(candidate, os.X_OK)
                and candidate not in matches):
            matches.append(candidate)
    if len(matches) > 1:
        raise FileNotFoundError(
            f'Executable {command!r} has multiple user installations: '
            + ', '.join(str(path) for path in matches)
            + '. Select one with its command parameter or environment variable.')
    return str(matches[0].resolve()) if matches else None


def automatic_requested(parameters):
    """Whether this ordinary KinBot input requests a new automatic solve."""
    return automatic_cbh_anl_requested(parameters)


def validate_automatic_environment(parameters):
    """Fail before reaction/QC work when a requested final solver is absent."""
    if not automatic_requested(parameters):
        return {}
    checked = {}

    def require_command(key, needed):
        if not needed:
            return
        environment_key = {
            'mess_command': 'KINBOT_MESS_COMMAND',
            'messpf_command': 'KINBOT_MESSPF_COMMAND',
        }.get(key)
        value = (os.environ.get(environment_key, '').strip()
                 if environment_key else '')
        if not value:
            value = parameters.get(key, key.removesuffix('_command'))
        try:
            command = shlex.split(str(value))
        except ValueError as error:
            raise ValueError(f'Invalid {key}: {error}') from error
        if not command:
            raise ValueError(f'{key} must name an executable.')
        executable = _resolve_user_executable(command[0])
        if executable is None:
            raise FileNotFoundError(
                f'{key} executable {command[0]!r} is unavailable; configure '
                'it before starting automatic CBH/ANL work.')
        command[0] = executable
        parameters[key] = shlex.join(command)
        checked[key] = executable

    require_command('messpf_command', bool(parameters.get('run_messpf')))
    require_command(
        'mess_command', parameters.get('me') == 1
        and bool(parameters.get('run_me')))
    if parameters.get('rotdpy_run'):
        from kinbot.rotdpy import ensure_available
        checked['rotdpy'] = ensure_available()
    if parameters.get('queuing') == 'slurm':
        for command in ('sbatch', 'squeue'):
            executable = shutil.which(command)
            if executable is None:
                raise FileNotFoundError(
                    f'{command} is unavailable for the requested Slurm run.')
            checked[command] = executable
    try:
        import atct.api  # noqa: F401
    except ImportError as error:
        raise ImportError(
            'Automatic CBH/ANL requires the ATcT Python client.') from error
    checked['atct_client'] = 'available'
    return checked


def _sha256_json(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _finite(value, label):
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f'{label} is unavailable.') from error
    if not math.isfinite(result):
        raise ValueError(f'{label} is unavailable.')
    return result


def _resolved_l2_theory(parameters):
    """Return the L2 identity used by both KinBot and automatic assembly."""
    values = dict(parameters)
    for key, default in {
            'qc': 'gauss', 'method': 'b3lyp', 'basis': '6-31G',
            'high_level_method': 'b3lyp', 'high_level_basis': '6-31G',
            'qc_command': 'g16', 'profiled_theory': 1,
            'theory_preset': '', 'l1': '', 'l2': '',
            'l1_profile': {}, 'l2_profile': {}, 'l3_overrides': {},
            'l3_resource_overrides': {}, 'fc_model_path': '__planning__',
            'fc_task_name': 'omol', 'fc_device': 'cpu'}.items():
        values.setdefault(key, default)
    profiles = resolve_profiles(values)
    try:
        profile = profiles['l2']
    except KeyError as error:
        raise ValueError('Automatic CBH/ANL requires a resolved L2 profile.') \
            from error
    method = profile.method.upper()
    basis = profile.basis
    if (profile.calculator != 'gaussian' or method not in ('B3LYP', 'B2PLYP')
            or basis.casefold() != 'cc-pvtz'):
        raise ValueError(
            'Automatic ANL0/ANL0-F12 requires Gaussian B3LYP/cc-pVTZ '
            'or B2PLYP-D3(BJ)/cc-pVTZ at L2.')
    kwargs = {str(key).casefold(): str(value).upper()
              for key, value in profile.calculator_kwargs.items()}
    dispersion = kwargs.get('empiricaldispersion', '')
    if method == 'B2PLYP' and dispersion != 'GD3BJ':
        raise ValueError(
            'Automatic B2PLYP L2 requires EmpiricalDispersion=GD3BJ.')
    if method == 'B3LYP' and dispersion:
        raise ValueError(
            'The automatic B3LYP ANL surface must not add dispersion.')
    composite_method = ('B2PLYP-D3BJ' if method == 'B2PLYP' else 'B3LYP')
    return {
        'calculator': 'gaussian', 'method': method, 'basis': 'cc-pVTZ',
        'dispersion': dispersion, 'composite_method': composite_method,
        'label': profile.label or (
            'B2PLYP-D3(BJ)/cc-pVTZ' if method == 'B2PLYP'
            else 'B3LYP/cc-pVTZ'),
    }


def _l2_species_record(species, l2_theory, *, identity_smiles=None,
                       routing_identity=None):
    """Freeze one accepted L2 point with an explicit energy identity.

    Stable points use their own canonical SMILES. A transition state has no
    unique equilibrium connectivity, so its identity is anchored to the
    reactant SMILES and its routing identity is the reaction name. Keeping
    that choice explicit avoids inventing a TS connectivity while still
    binding every component in the barrier difference to the same state.
    """
    symbols = list(map(str, species.atom))
    positions = [[float(value) for value in row] for row in species.geom]
    identity_smiles = _canonical_smiles(
        species.smiles if identity_smiles is None else identity_smiles)
    if not identity_smiles:
        raise ValueError('Accepted L2 point needs an explicit SMILES identity.')
    molecule = {
        'symbols': symbols, 'positions': positions,
        'charge': int(species.charge), 'multiplicity': int(species.mult),
    }
    geometry_sha256 = _geometry_hash(Atoms(symbols=symbols,
                                           positions=positions))
    payload = {
        'schema': 1, 'provider': 'kinbot-accepted-l2',
        'routing_name': (routing_name(species) if routing_identity is None
                         else str(routing_identity)),
        'smiles': identity_smiles,
        'molecule': molecule,
        'energy_hartree': _finite(species.energy, 'accepted L2 energy'),
        'harmonic_zpe_hartree': _finite(species.zpe,
                                        'accepted L2 harmonic ZPE'),
        'frequencies_cm_inverse': [
            _finite(value, 'accepted L2 frequency')
            for value in getattr(species, 'freq', ())],
        'geometry_sha256': geometry_sha256,
        'source_job': getattr(species, 'source_job', None),
        'source_row_id': getattr(species, 'source_row_id', None),
        'l2_theory': dict(l2_theory),
    }
    payload['record_sha256'] = _sha256_json(payload)
    molecule['smiles'] = payload['smiles']
    molecule['source'] = {
        'provider': payload['provider'],
        'geometry_sha256': geometry_sha256,
        'record_sha256': payload['record_sha256'],
        'energy_hartree': payload['energy_hartree'],
        'harmonic_zpe_hartree': payload['harmonic_zpe_hartree'],
        'source_job': payload['source_job'],
        'source_row_id': payload['source_row_id'],
        'l2_theory': dict(l2_theory),
    }
    payload['molecule'] = molecule
    return payload


def _validate_policy(policy):
    """Compatibility wrapper for the lightweight configuration validator."""
    return validate_anl_size_policy(policy)


def highest_tier(smiles, charge, multiplicity, policy,
                 requested_method='ANL0-F12'):
    """Select a deterministic upper calculation tier from CBH connectivity."""
    if requested_method not in ('ANL0', 'ANL0-F12'):
        raise ValueError('Automatic highest tier must be ANL0 or ANL0-F12.')
    policy = _validate_policy(policy)
    metrics = cbh_size_metrics(
        smiles, charge=charge, multiplicity=multiplicity)
    for policy_name, tier in (
            ('full_anl0_f12', requested_method),
            ('l3_geometry', 'L3'),
            ('l3_single_point', 'L3-SP@L2')):
        limits = policy[policy_name]
        if all(metrics[key] <= limits['max_' + key]
               for key in ('atoms', 'heavy_atoms', 'heavy_bonds')):
            return tier, metrics
    return 'L2', metrics


def _highest_reaction(smiles, rung, charge, multiplicity):
    for candidate in range(rung, -1, -1):
        reaction = generate_cbh_reaction(
            smiles, candidate, charge=charge,
            multiplicity=multiplicity)
        if reaction is not None:
            return reaction
    # Small saturated species such as methane are legitimate CBH-0 terminal
    # references: their only connectivity reaction would contain the target
    # itself. Invalid or state-inconsistent connectivity raises above and is
    # not converted into a terminal.
    return None


def build_network_request(root, parameters):
    """Freeze one completed KinBot reaction network and expand its CBH trees."""
    requested = str(parameters.get('composite_method', '')).upper()
    if requested not in ('ANL0', 'ANL0-F12'):
        raise ValueError(
            'Automatic production supports composite_method ANL0 or ANL0-F12.')
    rung_policy = parameters.get('cbh_rung', 'auto')
    if rung_policy == 'auto':
        max_rung = 3
    elif (not isinstance(rung_policy, bool)
          and isinstance(rung_policy, int) and rung_policy in range(4)):
        max_rung = rung_policy
    else:
        raise ValueError(
            'Automatic CBH rung policy must be auto, 0, 1, 2, or 3.')
    policy = _validate_policy(parameters.get(
        'anl_size_policy', DEFAULT_SIZE_POLICY))
    l2_theory = _resolved_l2_theory(parameters)

    stable = {}
    route_to_smiles = {}

    def add_stable(species, role):
        record = _l2_species_record(species, l2_theory)
        smiles = record['smiles']
        state = (record['molecule']['charge'],
                 record['molecule']['multiplicity'])
        key = (smiles, *state)
        existing = stable.get(key)
        if existing is not None:
            roles = sorted(set(existing['roles']) | {role})
            # Several reaction paths can reach the same product.  Use the
            # lowest accepted L2 E+ZPE record exactly once, matching KinBot's
            # conformer/product selection rule, and retain every network role.
            if (record['energy_hartree'] + record['harmonic_zpe_hartree'] <
                    existing['energy_hartree'] +
                    existing['harmonic_zpe_hartree']):
                record['roles'] = roles
                stable[key] = record
            else:
                existing['roles'] = roles
        else:
            record['roles'] = [role]
            stable[key] = record
        route = routing_name(species)
        previous = route_to_smiles.setdefault(route, smiles)
        if previous != smiles:
            raise ValueError(f'{route}: routing identity maps to two SMILES.')

    add_stable(root, 'initial_well')
    reactions = []
    transition_states = {}
    for index, reaction in enumerate(getattr(root, 'reac_obj', ())):
        statuses = getattr(root, 'reac_ts_done', ())
        if index >= len(statuses) or statuses[index] != -1:
            continue
        products = []
        for product in reaction.products:
            add_stable(product, 'reaction_product')
            products.append(routing_name(product))
        family = getattr(reaction, 'family_name', None) or getattr(
            root, 'reac_type', ['unknown'])[index]
        barrierless = family == 'hom_sci'
        ts_name = None
        if not barrierless:
            ts = reaction.ts
            ts_name = str(reaction.instance_name)
            if ts_name in transition_states:
                raise ValueError(f'Duplicate transition-state identity {ts_name}.')
            transition_states[ts_name] = _l2_species_record(
                ts, l2_theory, identity_smiles=root.smiles,
                routing_identity=ts_name)
            transition_states[ts_name]['stationary_point'] = 'transition_state'
            transition_states[ts_name]['identity_kind'] = \
                'reactant_connectivity_anchor'
        reactions.append({
            'name': str(reaction.instance_name), 'family': str(family),
            'reactant': routing_name(root), 'products': sorted(products),
            'transition_state': ts_name, 'barrierless': barrierless,
        })

    cbh_targets = {}
    calculated_states = {}
    terminals = {}
    for record in stable.values():
        molecule = record['molecule']
        reaction = _highest_reaction(
            record['smiles'], max_rung, molecule['charge'],
            molecule['multiplicity'])
        if reaction is None:
            calculated_states.setdefault(
                record['smiles'], [molecule['charge'],
                                   molecule['multiplicity']])
            formula = dict(Counter(record['molecule']['symbols']))
            terminals.setdefault(record['smiles'], {
                'formula': formula,
                'state': [molecule['charge'], molecule['multiplicity']]})
            cbh_targets[routing_name_from_record(record)] = {
                'target_smiles': record['smiles'],
                'requested_rung': 0, 'selected_structural_rung': None,
                'rung_policy': rung_policy,
                'terminal_atct': True,
            }
            continue
        # Plan every structurally valid rung through the selected maximum.
        # Runtime failures can therefore lower the rung without discovering
        # an unplanned reference species after licensed work has begun.
        available_rungs = []
        for candidate in range(reaction.rung, -1, -1):
            candidate_reaction = generate_cbh_reaction(
                record['smiles'], candidate, charge=molecule['charge'],
                multiplicity=molecule['multiplicity'])
            if candidate_reaction is None:
                continue
            available_rungs.append(candidate)
            calculated, leaves = cbh_network_requirements(candidate_reaction)
            states = cbh_network_states(candidate_reaction)
            if set(calculated) != set(states):
                raise RuntimeError(
                    'CBH calculated-species/state expansion disagrees.')
            for smiles, state in states.items():
                previous = calculated_states.setdefault(smiles, list(state))
                if previous != list(state):
                    raise ValueError(
                        f'{smiles}: combined CBH trees require conflicting states.')
            for smiles, (formula, state) in leaves.items():
                value = {'formula': formula, 'state': list(state)}
                previous = terminals.setdefault(smiles, value)
                if previous != value:
                    raise ValueError(
                        f'{smiles}: conflicting terminal CBH identity.')
        cbh_targets[routing_name_from_record(record)] = {
            'target_smiles': reaction.target_smiles,
            'requested_rung': reaction.rung,
            'selected_structural_rung': reaction.rung,
            'available_structural_rungs': available_rungs,
            'rung_policy': rung_policy,
            'reaction': asdict(reaction),
        }

    states = []
    live_by_identity = {
        (record['smiles'], record['molecule']['charge'],
         record['molecule']['multiplicity']): record
        for record in stable.values()}
    for smiles, (charge, multiplicity) in sorted(calculated_states.items()):
        live = live_by_identity.get((smiles, charge, multiplicity))
        tier, metrics = highest_tier(
            smiles, charge, multiplicity, policy, requested)
        states.append({
            'smiles': smiles, 'charge': charge,
            'multiplicity': multiplicity, 'highest_tier': tier,
            'fallback_tiers': list(_TIERS[_TIERS.index(tier):]),
            'size_metrics': metrics,
            'source': ('network_l2' if live is not None
                       else 'automatic_reference_l2'),
            'network_record': live,
        })
    for name, record in transition_states.items():
        molecule = record['molecule']
        tier, metrics = highest_tier(
            record['smiles'], molecule['charge'], molecule['multiplicity'],
            policy, requested)
        record['highest_tier'] = tier
        record['fallback_tiers'] = list(_TIERS[_TIERS.index(tier):])
        record['size_metrics'] = metrics

    thermochemistry_states = {}
    for entry in states:
        key = _state_key(entry['smiles'], entry['charge'],
                         entry['multiplicity'])
        thermochemistry_states[key] = {
            'smiles': entry['smiles'], 'charge': entry['charge'],
            'multiplicity': entry['multiplicity'],
            'source': entry['source'],
        }
    for smiles, terminal in terminals.items():
        charge, multiplicity = terminal['state']
        key = _state_key(smiles, charge, multiplicity)
        live = live_by_identity.get((smiles, charge, multiplicity))
        thermochemistry_states.setdefault(key, {
            'smiles': smiles, 'charge': charge,
            'multiplicity': multiplicity,
            'source': ('network_l2' if live is not None
                       else 'automatic_reference_l2'),
        })

    max_jobs = parameters.get('queue_job_limit', -1)
    if isinstance(max_jobs, bool) or not isinstance(max_jobs, int):
        raise ValueError('queue_job_limit must be an integer.')
    if max_jobs < 1:
        max_jobs = 1
    return {
        'schema': 1, 'status': 'planned',
        'requested_method': requested,
        'requested_cbh_rung': rung_policy,
        'l2_theory': l2_theory,
        'atct_version': str(parameters.get('atct_version', '1.220')),
        'resources': {
            'partition': parameters.get('queue_name'),
            'cores_per_job': int(parameters.get('ppn', 1)),
            'global_max_jobs': max_jobs,
        },
        'size_policy': policy,
        'network': {
            'initial': routing_name(root),
            'route_to_smiles': dict(sorted(route_to_smiles.items())),
            'stable_species': sorted(stable.values(),
                                     key=lambda item: item['routing_name']),
            'transition_states': transition_states,
            'reactions': sorted(reactions, key=lambda item: item['name']),
        },
        'cbh': {
            'targets': cbh_targets, 'calculated_states': states,
            'terminal_atct_species': terminals,
            'thermochemistry_states': [
                thermochemistry_states[key]
                for key in sorted(thermochemistry_states)],
        },
    }


def routing_name_from_record(record):
    value = record.get('routing_name')
    if not isinstance(value, str) or not value:
        raise ValueError('Network L2 record has no routing identity.')
    return value


def write_network_request(root, parameters, path='anl_network_request.json'):
    """Write the immutable automatic plan beside the KinBot calculation."""
    result = build_network_request(root, parameters)
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    return destination


def _component_dict(component):
    value = asdict(component)
    value['settings'] = dict(value.get('settings', {}))
    return value


def _ase_energy_component(run_dir, task_id, *, key, method, basis, state_id):
    run_dir, spec, state = _load(run_dir)
    task = next((item for item in spec['tasks'] if item['id'] == task_id), None)
    entry = state['tasks'].get(task_id)
    if (task is None or task.get('kind') != 'ase_optimize'
            or entry is None or entry.get('status') != 'complete'):
        raise IncompleteRecipeError(f'{task_id}: completed optimizer is unavailable.')
    _verify_stage_files(run_dir, task, entry)
    execution_path = run_dir / 'tasks' / task_id / 'execution.json'
    execution = json.loads(execution_path.read_text())
    final_hash = _verify_execution(run_dir, task, entry, execution)
    value = execution.get('details', {}).get('energy_hartree')
    if (not isinstance(value, (int, float)) or isinstance(value, bool)
            or not math.isfinite(value) or final_hash is None):
        raise ValueError(f'{task_id}: optimizer energy record is invalid.')
    return ComponentResult(
        key=key, value_hartree=float(value), quantity='electronic',
        method=method, basis=basis,
        backend=task['profile']['calculator'].lower(), state_id=state_id,
        charge=spec['molecule'].get('charge', 0),
        multiplicity=spec['molecule'].get('multiplicity', 1),
        geometry_sha256=final_hash,
        source_sha256=_file_hash(execution_path),
        source=str(execution_path),
        settings={'optimizer': execution['details'].get('optimizer')},
        review_required=False)


def _l2_source_component(source_record, *, key, quantity, value, method,
                         basis, state_id):
    molecule = source_record['molecule']
    source = molecule['source']
    return ComponentResult(
        key=key, value_hartree=float(value), quantity=quantity,
        method=method, basis=basis, backend='kinbot', state_id=state_id,
        charge=molecule['charge'], multiplicity=molecule['multiplicity'],
        geometry_sha256=source['geometry_sha256'],
        source_sha256=source['record_sha256'],
        source=(f"KinBot accepted L2 record {source.get('source_job')}"
                f"#{source.get('source_row_id')}"),
        settings={'provider': source.get('provider')}, review_required=False)


def _vpt2_zpe_components(run_dir, *, state_id):
    run_dir, spec, task, execution, native, parsed = _verified_task_result(
        run_dir, 'gaussian_vpt2')
    if parsed.get('review_required'):
        return None
    try:
        harmonic_cm = float(parsed['zpe_cm_inverse']['harmonic'])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError('gaussian_vpt2: harmonic ZPE is missing.') from error
    from ase.units import Hartree, invcm
    harmonic = harmonic_cm * invcm / Hartree
    common = dict(
        basis=parsed['basis'], backend='gaussian', state_id=state_id,
        charge=spec['molecule'].get('charge', 0),
        multiplicity=spec['molecule'].get('multiplicity', 1),
        geometry_sha256=execution['geometry_sha256'],
        source_sha256=execution['artifacts'][task['result_parser']['file']],
        source=str(native), review_required=False)
    harmonic_component = ComponentResult(
        key='harmonic_zpe', value_hartree=harmonic, quantity='zpe',
        method=(parsed['method'] + ('-D3BJ' if
                parsed.get('dispersion') == 'GD3BJ' else '')),
        settings={'dispersion': parsed.get('dispersion')}, **common)
    correction = task_component(
        run_dir, 'gaussian_vpt2', key='vpt2_correction',
        state_id=state_id)
    return harmonic_component, correction


def _simple_composite(recipe_name, molecule, components, electronic_keys,
                      zero_point_keys):
    if not components or not electronic_keys or not zero_point_keys:
        raise ValueError('Reduced composite expression is incomplete.')
    if set(electronic_keys) | set(zero_point_keys) != set(components):
        raise ValueError('Reduced composite expression does not use each component.')
    electronic = math.fsum(components[key].value_hartree
                           for key in electronic_keys)
    zero_point = math.fsum(components[key].value_hartree
                           for key in zero_point_keys)
    return {
        'schema': 1, 'status': 'complete', 'recipe': recipe_name,
        'recipe_version': 1,
        'charge': molecule.get('charge', 0),
        'multiplicity': molecule.get('multiplicity', 1),
        'electronic_hartree': electronic,
        'zero_point_hartree': zero_point,
        'zero_k_hartree': math.fsum((electronic, zero_point)),
        'components': {key: _component_dict(value)
                       for key, value in components.items()},
        'electronic_terms': [
            {'component': key, 'coefficient': 1}
            for key in electronic_keys],
        'zero_point_terms': [
            {'component': key, 'coefficient': 1}
            for key in zero_point_keys],
    }


def _component_family_omissions(availability):
    """Return every common-component fallback, best expression first."""
    if (not isinstance(availability, dict)
            or set(availability) != set(ANL_REDUCIBLE_FAMILIES)
            or any(not isinstance(value, bool)
                   for value in availability.values())):
        raise ValueError('ANL component-family availability is invalid.')
    missing = {
        family for family in ANL_REDUCIBLE_FAMILIES
        if not availability[family]}
    available = tuple(
        family for family in ANL_REDUCIBLE_FAMILIES
        if availability[family])
    result = []
    for count in range(len(available) + 1):
        for extra in combinations(available, count):
            omitted = missing | set(extra)
            result.append(tuple(
                family for family in ANL_REDUCIBLE_FAMILIES
                if family in omitted))
    return tuple(result)


def assemble_available_profiles(run_dir, source_record, destination, *,
                                state_id, stationary_point='minimum'):
    """Assemble every complete tier without converting bad data to fallback."""
    run_dir, spec, state = _load(run_dir)
    if spec.get('name') != 'anl-automatic-post-l2':
        raise ValueError('Automatic profile assembly needs an automatic graph.')
    if stationary_point not in ('minimum', 'transition_state'):
        raise ValueError('Invalid stationary point.')
    molecule = spec['molecule']
    theory = source_record.get('l2_theory')
    if (not isinstance(theory, dict)
            or theory.get('composite_method') not in ('B3LYP',
                                                       'B2PLYP-D3BJ')
            or theory.get('basis') != 'cc-pVTZ'):
        raise ValueError('Automatic L2 theory provenance is incomplete.')
    l2_method = theory['composite_method']
    l2_label = theory['label']
    statuses = {task['id']: state['tasks'].get(task['id'], {}).get(
        'status', 'waiting') for task in spec['tasks']}
    unfinished = {key: value for key, value in statuses.items()
                  if value not in ('complete', 'unavailable', 'skipped')}
    if unfinished:
        raise RuntimeError(f'Automatic graph is unfinished: {unfinished}')
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    profiles = {'thermochemistry': {}, 'barrier': {}}

    def save(role, label, record, filename):
        path = destination / filename
        path.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
        profiles[role][label] = str(path)

    declared = set(statuses)
    electrons = cbh_size_metrics(
        source_record['smiles'], charge=molecule.get('charge', 0),
        multiplicity=molecule.get('multiplicity', 1))['electrons']
    vpt2 = None
    if statuses.get('gaussian_vpt2') == 'complete':
        vpt2 = _vpt2_zpe_components(run_dir, state_id=state_id)
    source_harmonic = _l2_source_component(
        source_record, key='harmonic_zpe', quantity='zpe',
        value=source_record['harmonic_zpe_hartree'], method=l2_method,
        basis=theory['basis'], state_id=state_id)
    source_electronic = _l2_source_component(
        source_record, key='electronic', quantity='electronic',
        value=source_record['energy_hartree'], method=l2_method,
        basis=theory['basis'], state_id=state_id)

    family_tasks = {
        'higher_order': (() if electrons <= 2
                         else ('ccsdt_dz', 'ccsdtq_dz')),
        'dboc': ('cfour_dboc',),
        'scalar_relativistic': ('rel_dkh', 'rel_nonrel'),
        'core_valence': (() if all(
            symbol in ('H', 'He') for symbol in molecule['symbols']) else
            ('cv_ae_tz', 'cv_ae_qz', 'cv_fc_tz', 'cv_fc_qz')),
    }
    family_available = {
        family: all(task in declared and statuses[task] == 'complete'
                    for task in tasks)
        for family, tasks in family_tasks.items()
    }
    l3_status = statuses.get('l3_geometry')
    if l3_status in ('complete', 'unavailable', 'skipped'):
        common = dict(
            state_id=state_id, spin_orbit_hartree=0.,
            spin_orbit_source=(
                'balanced-CBH gauge zero; composition-balanced reaction '
                'energies are invariant to this per-species convention'),
            spin_orbit_backend='known_zero', base_run=run_dir,
            vpt2_method=l2_method)
        for family, assembler, reference_ids, stem in (
                ('ANL0-F12', assemble_profiled_anl0_f12,
                 ('f12_tz', 'f12_qz'), 'anl0_f12'),
                ('ANL0', assemble_profiled_anl0,
                 ('anl0_qz', 'anl0_5z'), 'anl0')):
            if not (set(reference_ids) <= declared and all(
                    statuses[key] == 'complete' for key in reference_ids)):
                continue
            # Emit every scientifically valid subset of complete additive
            # families.  The CBH solver can then choose the best *common*
            # expression across a reaction without rerunning a species that
            # completed more corrections than another member.  Missing
            # paired corrections are always omitted atomically.
            omission_sets = _component_family_omissions(family_available)
            for index, omitted in enumerate(omission_sets):
                omission_statuses = {
                    item: {task: statuses.get(task, 'not_declared')
                           for task in family_tasks[item]}
                    for item in omitted
                }

                def annotate(record):
                    if omitted:
                        record['reduced_recipe'].update({
                            'reason': (
                                'automatic common-component fallback; no '
                                'partial correction differences are used'),
                            'omitted_family_task_statuses': omission_statuses,
                        })
                    return record

                suffix = '' if not omitted else f'_reduced_{index:02d}'
                # Build the barrier expression without the VPT2 term. The
                # separate label prevents thermochemical and barrier records
                # from ever being confused.
                harmonic = annotate(assembler(
                    run_dir, run_dir, run_dir, harmonic_only=True,
                    omitted_families=omitted,
                    l2_geometry_fallback=source_harmonic, **common))
                save('barrier', harmonic['recipe'], harmonic,
                     f'{stem}{suffix}_harmonic.json')
                if stationary_point == 'minimum' and vpt2 is not None:
                    thermal = annotate(assembler(
                        run_dir, run_dir, run_dir, harmonic_only=False,
                        omitted_families=omitted,
                        l2_geometry_fallback=source_harmonic, **common))
                    save('thermochemistry', thermal['recipe'], thermal,
                         f'{stem}{suffix}.json')

    if statuses.get('l3_geometry') == 'complete' and \
            statuses.get('harmonic') == 'complete':
        electronic = _ase_energy_component(
            run_dir, 'l3_geometry', key='electronic', method='CCSD(T)',
            basis='cc-pVTZ', state_id=state_id)
        harmonic = task_component(
            run_dir, 'harmonic', key='harmonic_zpe', state_id=state_id)
        harmonic_components = {
            'electronic': electronic, 'harmonic_zpe': harmonic}
        harmonic_label = (
            'L3:CCSD(T)/cc-pVTZ//CCSD(T)/cc-pVTZ:harmonic-ZPE')
        harmonic_record = _simple_composite(
            harmonic_label, molecule, harmonic_components, ('electronic',),
            ('harmonic_zpe',))
        save('barrier', harmonic_label, harmonic_record, 'l3_harmonic.json')
        components = dict(harmonic_components)
        suffix = 'harmonic-ZPE'
        if stationary_point == 'minimum' and vpt2 is not None:
            components['vpt2_correction'] = vpt2[1]
            suffix = l2_method + '-anharmonic-correction'
        label = f'L3:CCSD(T)/cc-pVTZ//CCSD(T)/cc-pVTZ:{suffix}'
        save('thermochemistry', label, _simple_composite(
            label, molecule, components, ('electronic',),
            tuple(key for key in components if key != 'electronic')),
             'l3.json')

    if statuses.get('l3_sp_l2') == 'complete':
        electronic = task_component(
            run_dir, 'l3_sp_l2', key='electronic', state_id=state_id)
        harmonic_components = {
            'electronic': electronic, 'harmonic_zpe': source_harmonic}
        harmonic_label = (
            f'L3:CCSD(T)/cc-pVTZ//{l2_label}:harmonic-ZPE')
        save('barrier', harmonic_label, _simple_composite(
            harmonic_label, molecule, harmonic_components, ('electronic',),
            ('harmonic_zpe',)), 'l3_sp_l2_harmonic.json')
        components = {'electronic': electronic}
        if stationary_point == 'minimum' and vpt2 is not None:
            components.update(harmonic_zpe=vpt2[0],
                              vpt2_correction=vpt2[1])
            suffix = l2_method + '-anharmonic-ZPE'
        else:
            components['harmonic_zpe'] = source_harmonic
            suffix = l2_method + '-harmonic-ZPE'
        label = f'L3:CCSD(T)/cc-pVTZ//{l2_label}:{suffix}'
        save('thermochemistry', label, _simple_composite(
            label, molecule, components, ('electronic',),
            tuple(key for key in components if key != 'electronic')),
             'l3_sp_l2.json')

    harmonic_components = {
        'electronic': source_electronic, 'harmonic_zpe': source_harmonic}
    harmonic_label = f'L2:{l2_label}:harmonic-ZPE'
    save('barrier', harmonic_label, _simple_composite(
        harmonic_label, molecule, harmonic_components, ('electronic',),
        ('harmonic_zpe',)), 'l2_harmonic.json')
    components = {'electronic': source_electronic}
    if stationary_point == 'minimum' and vpt2 is not None:
        components.update(harmonic_zpe=vpt2[0], vpt2_correction=vpt2[1])
        suffix = 'anharmonic-ZPE'
    else:
        components['harmonic_zpe'] = source_harmonic
        suffix = 'harmonic-ZPE'
    label = f'L2:{l2_label}:{suffix}'
    save('thermochemistry', label, _simple_composite(
        label, molecule, components, ('electronic',),
        tuple(key for key in components if key != 'electronic')),
         'l2.json')
    # A missing or rejected anharmonic correction must still leave a complete
    # same-method harmonic profile available to every member of a CBH
    # reaction.  Exact labels below prevent ambiguity between the two.
    for label, path in profiles['barrier'].items():
        profiles['thermochemistry'].setdefault(label, path)
    (destination / 'profiles.json').write_text(json.dumps({
        'schema': 1, 'status': 'complete', 'state_id': state_id,
        'stationary_point': stationary_point, 'task_statuses': statuses,
        'profiles': profiles,
    }, indent=2, sort_keys=True) + '\n')
    return profiles


def _read_json(path, label):
    try:
        value = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f'{label} is not valid JSON: {path}') from error
    if not isinstance(value, dict):
        raise ValueError(f'{label} must be a JSON object: {path}')
    return value


def _state_key(smiles, charge, multiplicity):
    return f'{_canonical_smiles(smiles)}|{int(charge)}|{int(multiplicity)}'


def _safe_state_directory(smiles, charge, multiplicity):
    digest = hashlib.sha256(_state_key(
        smiles, charge, multiplicity).encode()).hexdigest()[:16]
    return f'state_{digest}'


def _reference_input(source_input, entry, directory, model_path):
    """Create one ordinary, restartable KinBot L2/HIR reference input."""
    source = _read_json(source_input, 'parent KinBot input')
    source.update({
        'title': 'automatic_cbh_reference_' + directory.name,
        'smiles': entry['smiles'], 'structure': [], 'bimol': 0,
        'charge': entry['charge'], 'mult': entry['multiplicity'],
        'reaction_search': 0, 'selected_reactions': {},
        'pes': 0, 'specific_reaction': 0, 'skip_chemids': [],
        'stop_product_reaction_search': 1,
        'conformer_search': 1, 'high_level': 1, 'rotor_scan': 1,
        'vrc_tst_scan': 0, 'vrc_tst_noscan': 0, 'rotdpy_run': 0,
        'me': 0, 'run_me': 0, 'run_messpf': 0,
        'cbh_rung': 0, 'anl_handoff_manifest': '',
        'cbh_thermochemistry_models': '',
        'thermochemistry_model_output': str(model_path),
        # Several references may run concurrently, but each reference owns
        # at most one live L1/L2 scheduler slot.
        'queue_job_limit': 1, 'do_clean': False,
    })
    for key in ('fc_model_path', 'q_temp', 'q_temp_hi', 'q_temp_mp2',
                'q_temp_l3'):
        value = source.get(key)
        if isinstance(value, str) and value and not os.path.isabs(value):
            candidate = (Path(source_input).resolve().parent / value).resolve()
            if candidate.exists():
                source[key] = str(candidate)
    path = directory / 'input.json'
    path.write_text(json.dumps(source, indent=2, sort_keys=True) + '\n')
    return path


def _model_record(model_path, entry, l2_theory):
    library = _read_json(model_path, 'automatic L2 thermochemistry model')
    models = library.get('species')
    if library.get('schema') != 1 or not isinstance(models, dict):
        raise ValueError(f'{model_path}: invalid thermochemistry model library.')
    canonical = _canonical_smiles(entry['smiles'])
    matches = [value for key, value in models.items()
               if _canonical_smiles(key) == canonical]
    if len(matches) != 1:
        raise ValueError(f'{model_path}: expected one model for {canonical}.')
    model = matches[0]
    if (model.get('charge') != entry['charge']
            or model.get('multiplicity') != entry['multiplicity']):
        raise ValueError(f'{canonical}: automatic L2 model state differs.')
    symbols = list(model.get('symbols', ()))
    positions = model.get('positions')
    atoms = Atoms(symbols=symbols, positions=positions)
    source = model.get('source', {})
    energy = _finite(source.get('electronic_energy_hartree'),
                     f'{canonical} L2 electronic energy')
    zpe = _finite(source.get('harmonic_zpe_hartree'),
                  f'{canonical} L2 harmonic ZPE')
    digest = _file_hash(model_path)
    geometry_hash = _geometry_hash(atoms)
    source_geometry_hash = source.get('geometry_sha256')
    legacy_geometry_hash = hashlib.sha256(np.asarray(
        atoms.positions, dtype='<f8').tobytes()).hexdigest()
    if source_geometry_hash not in (geometry_hash, legacy_geometry_hash):
        raise ValueError(
            f'{canonical}: automatic L2 model geometry hash differs.')
    molecule = {
        'symbols': symbols, 'positions': atoms.positions.tolist(),
        'charge': entry['charge'], 'multiplicity': entry['multiplicity'],
        'smiles': canonical,
        'source': {
            'provider': 'automatic-kinbot-reference-l2',
            'geometry_sha256': geometry_hash, 'record_sha256': digest,
            'exported_geometry_sha256': source_geometry_hash,
            'geometry_hash_format': (
                'symbols-and-positions-json'
                if source_geometry_hash == geometry_hash
                else 'legacy-position-bytes'),
            'energy_hartree': energy, 'harmonic_zpe_hartree': zpe,
            'model_file': str(Path(model_path).resolve()),
            'l2_theory': dict(l2_theory),
        },
    }
    return {
        'schema': 1, 'provider': 'automatic-kinbot-reference-l2',
        'routing_name': model.get('name', _safe_state_directory(
            canonical, entry['charge'], entry['multiplicity'])),
        'smiles': canonical, 'molecule': molecule,
        'energy_hartree': energy, 'harmonic_zpe_hartree': zpe,
        'frequencies_cm_inverse': list(model.get(
            'frequencies_cm_inverse', ())),
        'l2_theory': dict(l2_theory),
        'geometry_sha256': geometry_hash, 'record_sha256': digest,
        'thermochemistry_model': model,
        'thermochemistry_model_file': str(Path(model_path).resolve()),
    }


def _pid_matches_reference(pid, input_path, directory):
    """Reject stale or reused PIDs before attaching to a child driver."""
    proc = Path('/proc') / str(pid)
    try:
        command = (proc / 'cmdline').read_bytes().replace(b'\0', b' ').decode()
        cwd = (proc / 'cwd').resolve()
    except (OSError, UnicodeDecodeError):
        return False
    return (cwd == directory.resolve() and 'kinbot.kb' in command
            and str(input_path.resolve()) in command)


def _run_reference_l2(entries, source_input, root, max_jobs, l2_theory):
    """Run external CBH references through ordinary KinBot, with one slot each."""
    root.mkdir(parents=True, exist_ok=True)
    pending = []
    records = {}
    for entry in entries:
        key = _state_key(entry['smiles'], entry['charge'],
                         entry['multiplicity'])
        directory = root / _safe_state_directory(
            entry['smiles'], entry['charge'], entry['multiplicity'])
        model = directory / 'thermochemistry_model.json'
        if model.is_file():
            records[key] = _model_record(model, entry, l2_theory)
            continue
        directory.mkdir(parents=True, exist_ok=True)
        input_path = _reference_input(source_input, entry, directory, model)
        pid_path = directory / 'driver.pid'
        if pid_path.is_file():
            try:
                pid = int(pid_path.read_text().strip())
                os.kill(pid, 0)
                if not _pid_matches_reference(pid, input_path, directory):
                    raise ProcessLookupError
            except (OSError, ValueError):
                pass
            else:
                pending.append((key, entry, directory, input_path, model, pid))
                continue
        pending.append((key, entry, directory, input_path, model, None))

    active = []
    while pending or active:
        while pending and len(active) < max_jobs:
            key, entry, directory, input_path, model, existing_pid = \
                pending.pop(0)
            if existing_pid is None:
                stdout = (directory / 'driver.stdout').open('a')
                stderr = (directory / 'driver.stderr').open('a')
                process = subprocess.Popen(
                    [sys.executable, '-m', 'kinbot.kb', str(input_path)],
                    cwd=directory, stdout=stdout, stderr=stderr)
                (directory / 'driver.pid').write_text(f'{process.pid}\n')
                pid = process.pid
            else:
                process, pid, stdout, stderr = None, existing_pid, None, None
            active.append((process, pid, key, entry, directory, model,
                           stdout, stderr))
        if not active:
            break
        time.sleep(2.)
        remaining = []
        for process, pid, key, entry, directory, model, stdout, stderr in active:
            if process is not None:
                returncode = process.poll()
                running = returncode is None
            else:
                try:
                    os.kill(pid, 0)
                except OSError:
                    running, returncode = False, None
                else:
                    running, returncode = True, None
            if running:
                remaining.append((process, pid, key, entry, directory, model,
                                  stdout, stderr))
                continue
            if stdout is not None:
                stdout.close()
            if stderr is not None:
                stderr.close()
            if (returncode not in (None, 0)) or not model.is_file():
                diagnostic = directory / 'driver.stderr'
                tail = (diagnostic.read_text(errors='replace').splitlines()[-60:]
                        if diagnostic.is_file() else [])
                raise RuntimeError(
                    f'{entry["smiles"]}: automatic L2 reference run failed '
                    f'with status {returncode}.\n' + '\n'.join(tail))
            records[key] = _model_record(model, entry, l2_theory)
        active = remaining
    return records


def _prepare_graphs(plan, sources, root):
    root.mkdir(parents=True, exist_ok=True)
    graphs = {}
    for entry in plan['cbh']['calculated_states']:
        key = _state_key(entry['smiles'], entry['charge'],
                         entry['multiplicity'])
        record = sources[key]
        directory = root / _safe_state_directory(
            entry['smiles'], entry['charge'], entry['multiplicity'])
        graph = directory / 'composite'
        if not (graph / 'workflow.json').is_file():
            directory.mkdir(parents=True, exist_ok=True)
            specification = automatic_composite_spec(
                record['molecule'], highest_tier=entry['highest_tier'],
                max_nodes=plan['resources']['global_max_jobs'],
                partition=plan['resources']['partition'],
                vpt2_method=plan['l2_theory']['method'],
                vpt2_basis=plan['l2_theory']['basis'],
                vpt2_dispersion=plan['l2_theory']['dispersion'])
            prepare(specification, graph)
        preflight(graph, skip_optional=True)
        graphs['stable:' + key] = graph
    for name, entry in plan['network']['transition_states'].items():
        molecule = entry['molecule']
        directory = root / ('ts_' + hashlib.sha256(
            name.encode()).hexdigest()[:16])
        graph = directory / 'composite'
        if not (graph / 'workflow.json').is_file():
            directory.mkdir(parents=True, exist_ok=True)
            specification = automatic_composite_spec(
                molecule, highest_tier=entry['highest_tier'],
                max_nodes=plan['resources']['global_max_jobs'],
                partition=plan['resources']['partition'],
                stationary_point='transition_state',
                vpt2_method=plan['l2_theory']['method'],
                vpt2_basis=plan['l2_theory']['basis'],
                vpt2_dispersion=plan['l2_theory']['dispersion'])
            prepare(specification, graph)
        preflight(graph, skip_optional=True)
        graphs['ts:' + name] = graph
    return graphs


def _drive_graphs(graphs, max_jobs, interval=20.):
    """Advance many dispatcher DAGs under one global live-job cap."""
    names = sorted(graphs)
    cursor = 0
    while True:
        states = {name: advance(graphs[name], submit=False) for name in names}
        hard = []
        active = 0
        complete = True
        for name, state in states.items():
            for task, entry in state['tasks'].items():
                status = entry['status']
                active += status == 'submitted'
                complete &= status in ('complete', 'unavailable', 'skipped')
                if status in ('failed', 'blocked', 'submitting'):
                    hard.append(f'{name}/{task}={status}')
        if hard:
            raise RuntimeError('Automatic composite graph failed: ' +
                               ', '.join(hard))
        if complete:
            return states
        slots = max_jobs - active
        for offset in range(len(names)):
            if slots <= 0:
                break
            name = names[(cursor + offset) % len(names)]
            state = states[name]
            staged = next((task for task, entry in state['tasks'].items()
                           if entry['status'] == 'staged'), None)
            if staged is None:
                continue
            states[name] = advance(
                graphs[name], submit=True, submit_only=[staged])
            slots -= 1
        cursor = (cursor + 1) % max(len(names), 1)
        time.sleep(float(interval))


def _terminal_formation(smiles, state, formula, table, energy):
    record = table.gas_by_smiles(smiles, formula=formula)
    if (record.formation_0k_kj_mol is None
            or record.formation_298k_kj_mol is None):
        raise ValueError(f'{smiles}: terminal ATcT record lacks 0/298 K values.')
    return FormationEnthalpy(
        target_smiles=smiles, rung=0, method=energy.method,
        reaction_energy_0k_kj_mol=0.,
        formation_0k_kj_mol=record.formation_0k_kj_mol,
        references={smiles: record}, atct_version=table.version,
        atct_source_sha256=table.source_sha256,
        energy_sources={smiles: energy.source},
        stoichiometry={}, formulas={smiles: dict(formula)},
        states={smiles: tuple(state)})


def _write_formation(path, selection, *, terminal=False):
    if terminal:
        formation, reaction, skipped = selection, None, ()
    else:
        formation, reaction, skipped = (
            selection.formation, selection.reaction, selection.skipped)
    payload = {
        'schema': 1, 'status': 'complete',
        'reaction': None if reaction is None else asdict(reaction),
        'formation': asdict(formation), 'skipped_ladder_entries': list(skipped),
    }
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    return formation


def _reference_thermochemistry_models(reference_entries, sources, graphs):
    """Return external CBH models, enriching only states with an ANL graph.

    Terminal CBH-0 leaves need an L2 partition-function model so that the
    same CBH reaction can be evaluated at finite temperature.  They do not
    need an electronic composite graph because ATcT supplies their formation
    enthalpy.  Calculated intermediate references may have both records and
    can use a completed VPT2 fundamental set when it is dimensionally valid.
    """
    models = {}
    for entry in reference_entries:
        key = _state_key(entry['smiles'], entry['charge'],
                         entry['multiplicity'])
        model = dict(sources[key]['thermochemistry_model'])
        graph = graphs.get('stable:' + key)
        if graph is not None:
            state = _read_json(graph / 'state.json', 'dispatcher state')
            if state['tasks'].get('gaussian_vpt2', {}).get(
                    'status') == 'complete':
                _, _, task, execution, _, parsed = _verified_task_result(
                    graph, 'gaussian_vpt2')
                fundamentals = parsed.get(
                    'anharmonic_fundamentals_cm_inverse')
                if (not parsed.get('review_required')
                        and isinstance(fundamentals, list)
                        and len(fundamentals) == len(
                            model.get('frequencies_cm_inverse', ()))
                        and all(float(value) > 0. for value in fundamentals)):
                    model['frequencies_cm_inverse'] = list(
                        map(float, fundamentals))
                    model['frequency_source'] = {
                        'kind': 'gaussian_vpt2', 'run': str(graph),
                        'native_output_sha256': execution['artifacts'][
                            task['result_parser']['file']]}
        models[entry['smiles']] = model
    return models


def _artifact(path, label):
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f'{label} is missing: {path}')
    return {'path': str(path), 'sha256': _file_hash(path),
            'bytes': path.stat().st_size}


def _audit_automatic_graph(path):
    """Bind every terminal dispatcher state to its staged and native files."""
    run_dir, spec, state = _load(path)
    tasks = {}
    for task in spec['tasks']:
        ident = task['id']
        entry = state['tasks'].get(ident)
        if not isinstance(entry, dict):
            raise ValueError(f'{run_dir}: task {ident} has no final state.')
        status = entry.get('status')
        if status not in ('complete', 'unavailable', 'skipped'):
            raise ValueError(
                f'{run_dir}: task {ident} is not terminal ({status!r}).')
        record = {'status': status}
        if status == 'complete':
            execution = run_dir / 'tasks' / ident / 'execution.json'
            record['execution'] = _artifact(
                execution, f'{ident} execution')
        else:
            for key in ('unavailable_reason', 'skip_reason', 'skipped_by'):
                if key in entry:
                    record[key] = entry[key]
        tasks[ident] = record
    return {
        'workflow': _artifact(run_dir / 'workflow.json', 'workflow'),
        'state': _artifact(run_dir / 'state.json', 'dispatcher state'),
        'tasks': tasks,
    }


def finalize_automatic_run(automatic, parameters, *, run_dir='.'):
    """Write the self-contained audit record for a one-input automatic run.

    This is independent of the historical hand-authored end-to-end validation
    manifest. It audits the products requested by the ordinary KinBot input
    and refuses to call a partial run complete.
    """
    run_dir = Path(run_dir).expanduser().resolve()
    if not isinstance(automatic, dict) or automatic.get('status') != 'complete':
        raise ValueError('Automatic CBH/ANL execution is incomplete.')
    execution_path = Path(automatic.get(
        'execution_record', Path(automatic['plan']).with_name(
            'execution.json'))).resolve()
    execution = _read_json(execution_path, 'automatic execution')
    if execution.get('status') != 'complete':
        raise ValueError('Automatic execution record is incomplete.')
    plan_path = Path(execution['plan']).resolve()
    plan = _read_json(plan_path, 'automatic network request')

    formations = {}
    for route, value in sorted(execution.get('formations', {}).items()):
        artifact = _artifact(value, f'{route} formation')
        payload = _read_json(value, f'{route} formation')
        if (payload.get('schema') != 1 or payload.get('status') != 'complete'
                or not isinstance(payload.get('formation'), dict)):
            raise ValueError(f'{route}: automatic formation is incomplete.')
        formations[route] = {
            **artifact, 'reaction': payload.get('reaction'),
            'formation': payload['formation'],
            'skipped_ladder_entries': payload.get(
                'skipped_ladder_entries', []),
        }

    graphs = {
        name: _audit_automatic_graph(path)
        for name, path in sorted(execution.get('graphs', {}).items())}
    handoff = _artifact(execution['handoff_manifest'], 'handoff manifest')
    applied = _read_json(
        run_dir / 'me' / 'anl_handoff.json', 'applied ANL handoff')
    if (applied.get('status') != 'complete'
            or applied.get('manifest_sha256') != handoff['sha256']):
        raise ValueError('Applied ANL handoff does not match the plan.')

    messpf = None
    if parameters.get('run_messpf'):
        from kinbot.anl.thermochemistry import audit_messpf_run
        messpf = audit_messpf_run(run_dir)

    mess = {'status': 'not_requested'}
    if parameters.get('me') == 1:
        inputs = {
            path.name: _artifact(path, 'MESS input')
            for path in sorted((run_dir / 'me').glob('mess_*.inp'))}
        mess = {'status': ('requested' if parameters.get('run_me')
                           else 'prepared'), 'inputs': inputs}
        if parameters.get('run_me'):
            record_path = run_dir / 'me' / 'mess_execution.json'
            record = _read_json(record_path, 'MESS execution')
            if record.get('status') != 'complete' or not record.get('jobs'):
                raise ValueError('MESS execution record is incomplete.')
            for name, job in record['jobs'].items():
                output = Path(job.get('output', '')).expanduser().resolve()
                if (not output.is_file()
                        or _file_hash(output) != job.get('output_sha256')):
                    raise ValueError(f'{name}: MESS native output hash differs.')
            mess.update(status='complete', execution=_artifact(
                record_path, 'MESS execution'), jobs=record['jobs'])

    rotdpy = {}
    if parameters.get('rotdpy_run'):
        from kinbot.rotdpy import read_result, result_path
        required = parameters.get('rotdpy_validation_profile', 'interface')
        directory = run_dir / parameters.get('rotdpy_directory', 'rotdPy')
        for reaction in plan['network']['reactions']:
            if not reaction['barrierless']:
                continue
            input_path = directory / f'{reaction["name"]}.py'
            result = read_result(input_path, required_profile=required)
            correction = run_dir / 'vrctst' / f'corr_{reaction["name"]}.json'
            rotdpy[reaction['name']] = {
                'input': _artifact(input_path, 'ROTD_py input'),
                'result': _artifact(result_path(input_path), 'ROTD_py result'),
                'correction': _artifact(correction, 'VRC correction'),
                'surface_count': result['surface_count'],
                'validation': result.get('validation'),
            }

    report = {
        'schema': 1, 'status': 'complete',
        'requested_method': plan['requested_method'],
        'requested_cbh_rung': plan['requested_cbh_rung'],
        'network': {
            'initial': plan['network']['initial'],
            'reactions': plan['network']['reactions'],
            'stable_species': [entry['routing_name'] for entry in
                               plan['network']['stable_species']],
            'transition_states': sorted(
                plan['network']['transition_states']),
        },
        'plan': _artifact(plan_path, 'automatic network request'),
        'execution': _artifact(execution_path, 'automatic execution'),
        'handoff': {**handoff, 'applied': applied},
        'graphs': graphs, 'formations': formations,
        'messpf': messpf, 'rotdpy': rotdpy, 'mess': mess,
        'scope': {
            'implemented_composites': ['ANL0', 'ANL0-F12'],
            'cbh_rungs': [0, 1, 2, 3],
            'graph_state': [
                'integral Kekule bond order', 'formal charge',
                'atom-local radical electrons'],
        },
    }
    destination = execution_path.with_name('final_report.json')
    serialized = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if destination.is_file() and destination.read_text() != serialized:
        raise ValueError(
            'Automatic final report differs from the completed artifacts; '
            'use a fresh run directory or remove only the stale report.')
    destination.write_text(serialized)
    return str(destination)


def _choose_common_barrier(profile_sets):
    common = set.intersection(*(set(item) for item in profile_sets))
    ordered = _ordered_methods({'barrier': {label: None
                                            for label in common}})
    if ordered:
        return ordered[0]
    raise ValueError('No common harmonic energy is available for tight barriers.')


def _ordered_methods(energies):
    """Return exact automatic labels from most to least complete."""
    labels = set()
    for profiles in energies.values():
        labels.update(profiles)

    def rank(label):
        harmonic = 'harmonic-ZPE' in label
        if label.startswith('profiled:ANL0-F12'):
            tier = 0
        elif label.startswith('profiled:ANL0') or label == 'ANL0':
            tier = 1
        elif label.startswith('L3:') and '//CCSD(T)/cc-pVTZ' in label:
            tier = 2
        elif label.startswith('L3:'):
            tier = 3
        elif label.startswith('L2:'):
            tier = 4
        else:
            tier = 5
        marker = 'reduced[without='
        if marker in label:
            body = label.split(marker, 1)[1].split(']', 1)[0]
            omitted = tuple(body.split('+'))
        else:
            omitted = ()
        family_order = {name: index for index, name in enumerate(
            ANL_REDUCIBLE_FAMILIES)}
        omission_order = tuple(family_order[item] for item in omitted)
        # Preserve the requested electronic tier first. Within it, retain as
        # many additive families as possible; for equal counts, the declared
        # family order drops higher-order increments before larger baseline
        # corrections. Harmonic ZPE follows the same electronic expression's
        # accepted anharmonic profile.
        return tier, len(omitted), omission_order, harmonic, label

    return tuple(sorted(labels, key=rank))


def run_automatic_cbh_anl(root_species, parameters, input_file, *,
                          directory='anl_automatic', interval=20.):
    """Build, execute, solve, and publish a complete automatic CBH/ANL DAG."""
    base = Path(directory).expanduser().resolve()
    base.mkdir(parents=True, exist_ok=True)
    plan = build_network_request(root_species, parameters)
    plan_path = base / 'network_request.json'
    if plan_path.is_file():
        previous = _read_json(plan_path, 'automatic network request')
        if _sha256_json(previous) != _sha256_json(plan):
            raise ValueError(
                'Accepted KinBot network differs from the existing automatic plan; '
                'use a fresh calculation directory.')
    else:
        plan_path.write_text(json.dumps(
            plan, indent=2, sort_keys=True) + '\n')

    # Resolve the pinned terminal thermochemistry before launching any
    # licensed electronic-structure calculation. A proxy, cache, identity,
    # or ATcT-version problem must fail cheaply at the planning boundary.
    terminals = plan['cbh']['terminal_atct_species']
    table = load_atct_references(
        terminals, plan['atct_version'],
        Path(parameters.get('atct_cache',
                            '~/.cache/kinbot/atct')).expanduser())

    sources = {}
    for entry in plan['cbh']['calculated_states']:
        record = entry.get('network_record')
        if record is not None:
            sources[_state_key(entry['smiles'], entry['charge'],
                               entry['multiplicity'])] = record
    reference_entries = [
        entry for entry in plan['cbh']['thermochemistry_states']
        if entry['source'] == 'automatic_reference_l2']
    sources.update(_run_reference_l2(
        reference_entries, input_file, base / 'reference_l2',
        plan['resources']['global_max_jobs'], plan['l2_theory']))
    missing = [entry['smiles'] for entry in plan['cbh']['calculated_states']
               if _state_key(entry['smiles'], entry['charge'],
                             entry['multiplicity']) not in sources]
    if missing:
        raise RuntimeError('Automatic L2 source records are missing: ' +
                           ', '.join(missing))

    graphs = _prepare_graphs(plan, sources, base / 'states')
    _drive_graphs(graphs, plan['resources']['global_max_jobs'], interval)

    profiles = {}
    profile_root = base / 'profiles'
    for entry in plan['cbh']['calculated_states']:
        key = _state_key(entry['smiles'], entry['charge'],
                         entry['multiplicity'])
        profiles[key] = assemble_available_profiles(
            graphs['stable:' + key], sources[key],
            profile_root / _safe_state_directory(
                entry['smiles'], entry['charge'], entry['multiplicity']),
            state_id=key)
    for name, record in plan['network']['transition_states'].items():
        profiles['ts:' + name] = assemble_available_profiles(
            graphs['ts:' + name], record,
            profile_root / ('ts_' + hashlib.sha256(
                name.encode()).hexdigest()[:16]),
            state_id='ts:' + name, stationary_point='transition_state')

    energies = {}
    for entry in plan['cbh']['calculated_states']:
        key = _state_key(entry['smiles'], entry['charge'],
                         entry['multiplicity'])
        energies[entry['smiles']] = {
            label: zero_k_from_composite_record(entry['smiles'], path)
            for label, path in profiles[key]['thermochemistry'].items()}
    method_ladder = _ordered_methods(energies)
    formations = {}
    formation_paths = {}
    for record in plan['network']['stable_species']:
        route = record['routing_name']
        smiles = record['smiles']
        target = plan['cbh']['targets'][route]
        path = base / 'formations' / (hashlib.sha256(
            route.encode()).hexdigest()[:16] + '.json')
        path.parent.mkdir(parents=True, exist_ok=True)
        if target.get('terminal_atct'):
            key = _state_key(smiles, record['molecule']['charge'],
                             record['molecule']['multiplicity'])
            try:
                method = next(label for label in method_ladder
                              if label in energies[smiles])
            except StopIteration as error:
                raise ValueError(
                    f'{smiles}: no accepted energy exists for an ATcT '
                    'network state.') from error
            energy = energies[smiles][method]
            formation = _terminal_formation(
                smiles, (record['molecule']['charge'],
                         record['molecule']['multiplicity']),
                Counter(record['molecule']['symbols']), table, energy)
            _write_formation(path, formation, terminal=True)
        else:
            point = SimpleNamespace(
                smiles=smiles, charge=record['molecule']['charge'],
                mult=record['molecule']['multiplicity'], wellorts=0,
                atom=record['molecule']['symbols'])
            selection = select_cbh_ladder(
                point, energies, table,
                max_rung=target['requested_rung'],
                method_ladder=method_ladder)
            formation = _write_formation(path, selection)
        formations[route] = formation
        formation_paths[route] = path

    root_route = plan['network']['initial']
    root_smiles = plan['network']['route_to_smiles'][root_route]
    root_state = next(record for record in plan['network']['stable_species']
                      if record['routing_name'] == root_route)['molecule']
    root_key = _state_key(root_smiles, root_state['charge'],
                          root_state['multiplicity'])
    tight = [item for item in plan['network']['reactions']
             if item['transition_state'] is not None]
    barrier_label = None
    if tight:
        barrier_label = _choose_common_barrier([
            profiles[root_key]['barrier'],
            *(profiles['ts:' + item['transition_state']]['barrier']
              for item in tight)])

    handoff = {
        'schema': 1, 'mode': 'cbh-anl', 'strict': True,
        'species': {}, 'transition_states': {},
    }
    for record in plan['network']['stable_species']:
        route = record['routing_name']
        molecule = record['molecule']
        key = _state_key(record['smiles'], molecule['charge'],
                         molecule['multiplicity'])
        method = formations[route].method
        composite = profiles[key]['thermochemistry'].get(method)
        if composite is None:
            raise ValueError(f'{route}: selected CBH method has no composite file.')
        entry = {
            'composite': composite, 'formation': str(formation_paths[route]),
            'geometry': molecule,
        }
        if route == root_route and barrier_label is not None:
            entry['barrier_composite'] = profiles[key]['barrier'][barrier_label]
        graph = graphs['stable:' + key]
        if ('harmonic-ZPE' not in method
                and _read_json(
                    graph / 'state.json', 'dispatcher state')['tasks']
                .get('gaussian_vpt2', {}).get('status') == 'complete'):
            entry['vpt2'] = {'run': str(graph),
                             'task_id': 'gaussian_vpt2'}
        handoff['species'][route] = entry
    for reaction in tight:
        name = reaction['transition_state']
        handoff['transition_states'][name] = {
            'composite': profiles['ts:' + name]['barrier'][barrier_label],
            'identity_smiles': plan['network']['transition_states'][name][
                'smiles'],
        }

    handoff_path = base / 'handoff.json'
    handoff_path.write_text(json.dumps(handoff, indent=2, sort_keys=True) + '\n')

    models = _reference_thermochemistry_models(
        reference_entries, sources, graphs)
    model_path = base / 'cbh_thermochemistry_models.json'
    model_path.write_text(json.dumps({
        'schema': 1, 'species': models}, indent=2, sort_keys=True) + '\n')

    result = {
        'schema': 1, 'status': 'complete', 'plan': str(plan_path),
        'handoff_manifest': str(handoff_path),
        'thermochemistry_models': str(model_path),
        'graphs': {key: str(value) for key, value in graphs.items()},
        'formations': {key: str(value)
                       for key, value in formation_paths.items()},
        'selected_barrier_method': barrier_label,
    }
    result_path = base / 'execution.json'
    result['execution_record'] = str(result_path)
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    return result
