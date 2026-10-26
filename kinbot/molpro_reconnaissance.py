"""Prepare and audit basis-specific Molpro orbital reconnaissance jobs.

The utility deliberately stops before choosing an active space.  It produces
ROHF/RHF orbital listings, Molden files, and hashable WFU records at every
requested geometry.  A reviewed Molpro record can then be converted to the
``vrc_tst_state_model`` contract with the ``select`` command.
"""

from argparse import ArgumentParser
from contextlib import contextmanager, ExitStack
from copy import deepcopy
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import socket
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
from ase.data import atomic_numbers

from kinbot.anl.site import render_site_setup
from kinbot.molpro import (validate_vrc_state_model, automatic_vrc_state_model,
                           _render_vrc_state_model, full_valence_active_space,
                           validate_input_blocks, compact_molpro_titles,
                           molpro_wfu_name, _vrc_multireference_method,
                           vrc_multireference_energy_contract, validate_mr_correlation_options)
from kinbot.theory import fortran_threshold
from kinbot.molpro_solver import (apply_solver_options, validate_solver_options,
                                 validate_casscf_iterations)
from kinbot.execution import atomic_json, file_sha256, slurm_job_record, SLURM_TERMINAL_STATES
from kinbot.anl.runtime import (qc_runtime_environment, cleanup_qc_runtime,
                               MolproScratchError, molpro_io_failure)


_SAFE = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')


class ProbeFailure(RuntimeError):
    """A submitted probe ended without a validated scientific result."""


@contextmanager
def _manifest_lock(root, *, wait=0.):
    if not np.isfinite(wait) or wait < 0:
        raise ValueError('Owner wait must be finite and nonnegative.')
    root = Path(root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    with (root / '.reconnaissance.lock').open('a+') as lock:
        deadline = time.monotonic() + wait
        announced = False
        while True:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError as error:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise RuntimeError('Another reconnaissance controller is updating this directory; '
                                       'read reconnaissance.json for live progress.') from error
                if not announced:
                    print('Waiting for the existing review owner to finish; its jobs continue.', flush=True)
                    announced = True
                time.sleep(min(5., remaining))
        yield


def _coordinates(data, source):
    """Expand independently defined dissociations from molecular connectivity."""
    import networkx as nx
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.electronic_states import apply_electronic_state_model, state_average_roots
    from kinbot.vrc_correction import coupled_multiplicities
    parent = StationaryPoint('reconnaissance', data['charge'], data['multiplicity'],
                             atom=np.asarray(data['symbols']),
                             geom=np.asarray(data['geometry_angstrom'], dtype=float))
    parent.characterize()
    graph = nx.from_numpy_array(parent.bond)
    coordinates = data['coordinates']
    if coordinates == 'auto':
        coordinates = [{'name': f'bond_{left + 1}_{right + 1}',
                        'bond': [left + 1, right + 1]}
                       for left, right in nx.bridges(graph)]
    if not isinstance(coordinates, list) or not coordinates:
        raise ValueError('coordinates must define at least one separable bond.')
    parameters = {}
    if data.get('kinbot_input'):
        path = Path(data['kinbot_input']).expanduser()
        if not path.is_absolute():
            path = source.parent / path
        parameters = json.loads(path.read_text())
    policy = {
        'minimum_distance_scale': 0.9, 'near_distance_scale': 2.,
        'near_step_angstrom': 0.2, 'far_step_factor': 1.4,
        'asymptote_angstrom': 30., 'asymptote_tolerance_hartree': 1.e-6,
        'maximum_energy_step_hartree': 0.02,
        'maximum_occupation_change': 0.2,
        'occupation_lower_bound': 0.02, 'occupation_upper_bound': 1.98,
    }
    overrides = data.get('scan_policy', {})
    if not isinstance(overrides, dict) or set(overrides) - set(policy):
        raise ValueError('Unknown reconnaissance scan_policy controls.')
    policy.update(overrides)
    if (any(isinstance(value, bool) or not isinstance(value, (int, float))
            or not np.isfinite(value) or value <= 0. for value in policy.values())
            or policy['far_step_factor'] <= 1.
            or policy['minimum_distance_scale'] >= policy['near_distance_scale']
            or not 0. < policy['occupation_lower_bound'] < policy['occupation_upper_bound'] < 2.):
        raise ValueError('Reconnaissance scan_policy ranges must be finite and ordered.')
    data['scan_policy'] = policy
    anchors = []
    orientation_controls = data.get('orientation_probes')
    if orientation_controls is not None:
        if (not isinstance(orientation_controls, dict)
                or set(orientation_controls) - {'distance_scales', 'angles_degrees', 'axes', 'include_asymptote'}):
            raise ValueError('Unknown orientation_probes controls.')
        for key in ('distance_scales', 'angles_degrees'):
            values = orientation_controls.get(key)
            if (not isinstance(values, list) or not values
                    or any(isinstance(value, bool) or not isinstance(value, (int, float))
                           or not np.isfinite(value) or value <= 0. for value in values)):
                raise ValueError(f'orientation_probes.{key} needs positive finite values.')
        if (not isinstance(orientation_controls.get('axes'), list)
                or not orientation_controls['axes']
                or set(orientation_controls['axes']) - {'bond', 'perpendicular_1', 'perpendicular_2'}
                or not isinstance(orientation_controls.get('include_asymptote', True), bool)):
            raise ValueError('Invalid orientation probe axes or asymptote control.')
    for coordinate in coordinates:
        if (not isinstance(coordinate, dict)
                or set(coordinate) - {'name', 'bond', 'distances_angstrom', 'fragment_states'}
                or not _SAFE.fullmatch(str(coordinate.get('name', '')))):
            raise ValueError('Each coordinate needs a safe name and bond.')
        bond = coordinate.get('bond')
        if (not isinstance(bond, list) or len(bond) != 2
                or any(isinstance(index, bool) or not isinstance(index, int)
                       or not 1 <= index <= len(data['symbols']) for index in bond)):
            raise ValueError('Coordinate bond indices must be valid and one based.')
        left, right = (index - 1 for index in bond)
        split = graph.copy()
        if not split.has_edge(left, right):
            raise ValueError('Reconnaissance coordinate must be a molecular bond.')
        split.remove_edge(left, right)
        components = list(nx.connected_components(split))
        if len(components) != 2:
            raise ValueError('A single reconnaissance bond must separate two '
                             'fragments; ring opening needs an explicit multi-bond model.')
        declarations = coordinate.get('fragment_states')
        if parent.charge != 0 and declarations is None:
            raise ValueError('Charged dissociation requires explicit fragment '
                             'charge/state evidence; automatic neutral splitting is invalid.')
        parent.bond[left, right] = parent.bond[right, left] = 0
        products, maps = parent.start_multi_molecular()
        parent.bond[left, right] = parent.bond[right, left] = graph[left][right]['weight']
        if declarations is not None:
            if not isinstance(declarations, list) or len(declarations) != 2:
                raise ValueError('fragment_states must describe both separated fragments.')
            from rdkit import Chem
            for product, mapping in zip(products, maps):
                matches = [record for record in declarations if isinstance(record, dict)
                           and record.get('atom') in [int(index) + 1 for index in mapping]]
                if len(matches) != 1:
                    raise ValueError('Each fragment state needs one distinguishing parent atom.')
                record = matches[0]
                if (set(record) != {'atom', 'charge', 'multiplicity', 'smiles'}
                        or any(isinstance(record[key], bool) or not isinstance(record[key], int)
                               for key in ('atom', 'charge', 'multiplicity'))):
                    raise ValueError('A fragment state requires atom, charge, multiplicity, and smiles.')
                molecule = Chem.MolFromSmiles(str(record['smiles']))
                if molecule is None:
                    raise ValueError('Invalid fragment-state SMILES.')
                molecule = Chem.AddHs(molecule)
                evidence = nx.Graph()
                evidence.add_nodes_from((atom.GetIdx(), {'element': atom.GetSymbol()})
                                        for atom in molecule.GetAtoms())
                evidence.add_edges_from((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                                        for bond in molecule.GetBonds())
                actual = nx.from_numpy_array(product.bond01)
                nx.set_node_attributes(actual, dict(enumerate(map(str, product.atom))), 'element')
                electrons = sum(atomic_numbers[str(atom)] for atom in product.atom) - record['charge']
                spin = record['multiplicity'] - 1
                if (Chem.GetFormalCharge(molecule) != record['charge']
                        or spin < 0 or electrons < 0 or spin > electrons or (electrons - spin) % 2
                        or not nx.is_isomorphic(actual, evidence,
                            node_match=lambda first, second: first['element'] == second['element'])):
                    raise ValueError('Fragment-state charge, spin, or connectivity is inconsistent.')
                product.charge = record['charge']
                product.mult = record['multiplicity']
                product.smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule))
            if sum(product.charge for product in products) != parent.charge:
                raise ValueError('Fragment charges do not sum to the parent charge.')
        for product in products:
            declared_spin = product.mult
            apply_electronic_state_model(product, parameters)
            if declarations is not None and product.mult != declared_spin:
                raise ValueError('Fragment multiplicity conflicts with the electronic-state model.')
        multiplicities = coupled_multiplicities(*(int(product.mult) for product in products))
        if data['multiplicity'] not in multiplicities:
            raise ValueError('Dissociation fragment states cannot couple to the parent spin.')
        roots = state_average_roots(products, multiplicities, parameters)
        moving = sorted(next(part for part in components if right in part))
        initial = float(np.linalg.norm(parent.geom[right] - parent.geom[left]))
        distances = coordinate.get('distances_angstrom')
        if distances is None:
            near = initial * policy['near_distance_scale']
            asymptote = policy['asymptote_angstrom']
            if asymptote <= near:
                raise ValueError('The reconnaissance asymptote must exceed the near region.')
            distances = list(np.arange(initial * policy['minimum_distance_scale'],
                                       near, policy['near_step_angstrom']))
            distance = near
            while distance < asymptote:
                distances.append(distance)
                distance *= policy['far_step_factor']
            distances.extend([asymptote, initial])
        minimum_points = 4 if data.get('probe_active_space', False) else 2
        if (not isinstance(distances, list) or len(distances) < minimum_points
                or any(isinstance(value, bool) or not isinstance(value, (int, float))
                       or not np.isfinite(value) or value <= 0. for value in distances)):
            raise ValueError(f'Each coordinate needs at least {minimum_points} finite positive distances.')
        distances = sorted({float(f'{value:.10g}') for value in distances})
        if len(distances) < minimum_points:
            raise ValueError(f'Each coordinate needs at least {minimum_points} distinct distances.')
        context = {'multiplicities': list(multiplicities), 'roots': roots,
                   'fragment_state_source': 'declared' if declarations is not None else 'connectivity and electronic_state_models',
                   'fragment_states': [{'smiles': product.smiles,
                                        'multiplicity': int(product.mult),
                                        'charge': int(product.charge)} for product in products]}
        for distance in distances:
            token = f'{distance:.10g}'.replace('.', 'p').replace('+', '')
            anchors.append({'name': f'{coordinate["name"]}_r_{token}', 'bond': bond,
                            'moving_atoms': [index + 1 for index in moving],
                            'distance_angstrom': distance, '_channel': coordinate['name'],
                            '_context': context})
        if orientation_controls:
            direction = (parent.geom[right] - parent.geom[left]) / initial
            transverse = np.cross(direction, np.eye(3)[np.argmin(np.abs(direction))])
            transverse /= np.linalg.norm(transverse)
            axes = {'bond': direction, 'perpendicular_1': transverse,
                    'perpendicular_2': np.cross(direction, transverse)}
            orient_distances = [initial * scale for scale in orientation_controls['distance_scales']]
            if orientation_controls.get('include_asymptote', True):
                orient_distances.append(policy['asymptote_angstrom'])
            for component in components:
                if len(component) == 1:
                    continue
                pivot = left if left in component else right
                seen_geometries = set()
                for axis_name in orientation_controls['axes']:
                    for angle in orientation_controls['angles_degrees']:
                        for distance in sorted(set(orient_distances)):
                            suffix = f'f{pivot + 1}_{axis_name}_a{angle:g}_r{distance:.10g}'.replace('.', 'p')
                            anchor = {'name': f'{coordinate["name"]}_orient_{suffix}', 'bond': bond,
                                'moving_atoms': [index + 1 for index in moving],
                                'distance_angstrom': distance, '_context': context,
                                '_rotation': {'atoms': sorted(component), 'pivot': pivot,
                                              'axis': axes[axis_name].tolist(), 'angle_degrees': angle},
                                '_orientation': {'channel': coordinate['name'], 'pivot_atom': pivot + 1,
                                                 'axis': axis_name, 'angle_degrees': angle}}
                            oriented = _anchor_geometry(data, anchor)
                            baseline = _anchor_geometry(data, {key: value for key, value in anchor.items()
                                                               if key != '_rotation'})
                            # Rotation about a linear fragment's own axis, or
                            # about an atomic fragment, creates no new geometry.
                            digest = tuple(np.round(oriented.ravel(), 8))
                            if np.allclose(oriented, baseline, atol=1.e-8) or digest in seen_geometries:
                                continue
                            seen_geometries.add(digest)
                            anchors.append(anchor)
    return anchors


def require_explicit_probe_space(data):
    """Do not turn a request for orbital evidence into an implicit valence CI.

    Reading legacy manifests remains supported. New preparation/planning needs
    counted chemical candidates or an explicit full-valence diagnostic choice.
    """
    controls = data.get('probe_space', {})
    if (data.get('probe_active_space') and not data.get('probe_definition')
            and (not isinstance(controls, dict) or controls.get('mode') != 'full_valence')):
        raise ValueError('An active-space probe needs an explicit probe_definition; '
                         'a deliberate valence diagnostic requires '
                         'probe_space.mode=full_valence. No active space is inferred for new jobs.')


def _load(path):
    path = Path(path).expanduser().resolve()
    data = json.loads(path.read_text())
    required = {'title', 'symbols', 'geometry_angstrom', 'charge',
                'multiplicity', 'bases'}
    if not isinstance(data, dict) or data.get('schema') != 1 \
            or not required <= set(data):
        raise ValueError('Reconnaissance input must use schema 1 and define '
                         + ', '.join(sorted(required)) + '.')
    symbols = list(map(str, data['symbols']))
    if not isinstance(data['title'], str) or not _SAFE.fullmatch(data['title']):
        raise ValueError('Reconnaissance title must be a safe single-line name.')
    geometry = np.asarray(data['geometry_angstrom'], dtype=float)
    if (not symbols or geometry.shape != (len(symbols), 3)
            or not np.isfinite(geometry).all()
            or any(symbol not in atomic_numbers for symbol in symbols)):
        raise ValueError('Invalid reconnaissance atoms or geometry.')
    for key in ('charge', 'multiplicity'):
        value = data[key]
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(f'{key} must be an integer.')
    if data['multiplicity'] < 1:
        raise ValueError('multiplicity must be positive.')
    electrons = sum(atomic_numbers[symbol] for symbol in symbols) \
        - data['charge']
    if (electrons < 1 or data['multiplicity'] - 1 > electrons
            or (electrons - data['multiplicity'] + 1) % 2):
        raise ValueError('Charge and multiplicity have inconsistent parity.')
    bases = data['bases']
    if (not isinstance(bases, list) or not bases
            or any(not isinstance(item, dict)
                   or set(item) != {'name', 'basis'}
                   or not _SAFE.fullmatch(str(item['name']))
                   or not isinstance(item['basis'], str)
                   or not re.fullmatch(r'[A-Za-z0-9_+./*()-]+', item['basis']) for item in bases)):
        raise ValueError('bases must contain safe name/basis objects.')
    if len({item['name'] for item in bases}) != len(bases):
        raise ValueError('Basis job names must be unique.')
    if data.get('basis_source', 'library') != 'library':
        raise ValueError('Reconnaissance uses the native Molpro basis library.')
    ppn = data.get('ppn', 8)
    memory = data.get('memory_mw', 1000)
    if any(isinstance(value, bool) or not isinstance(value, int) or value < 1
           for value in (ppn, memory)):
        raise ValueError('ppn and memory_mw must be positive integers.')
    for key in ('memory_mb', 'queue_job_limit', 'scratch_min_mb'):
        if key in data and (isinstance(data[key], bool) or not isinstance(data[key], int)
                            or data[key] < 1):
            raise ValueError(f'{key} must be a positive integer.')
    if not isinstance(data.get('probe_active_space', False), bool):
        raise ValueError('probe_active_space must be a boolean.')
    character_only = data.get('orbital_character_only', False)
    if not isinstance(character_only, bool) or character_only and data.get('probe_active_space'):
        raise ValueError('orbital_character_only requires SCF-only diagnostics.')
    virtuals = data.get('orbital_print_virtuals', 0)
    if isinstance(virtuals, bool) or not isinstance(virtuals, int) or virtuals < 0:
        raise ValueError('orbital_print_virtuals must be a nonnegative integer.')
    threshold = data.get('orbital_print_threshold', .25)
    if (isinstance(threshold, bool) or not isinstance(threshold, (int, float))
            or not np.isfinite(threshold) or threshold < 0):
        raise ValueError('orbital_print_threshold must be finite and nonnegative.')
    for key, choices, default in (
            ('probe_method', ('casci', 'casscf'), 'casci'),
            ('avas_reference', ('atomic', 'rhf'), 'atomic'),
            ('orbital_initialization', ('propagated', 'independent'), 'propagated')):
        if data.get(key, default) not in choices:
            raise ValueError(f'{key} must be one of {choices}.')
    for key, allowed in (
            ('scf_options', {'solver', 'max_iterations'}),
            ('casscf_options', {'solver', 'exact_diagonalization', 'quasi_newton', 'max_iterations',
                               'energy_convergence', 'gradient_convergence', 'verbosity'})):
        options = data.get(key, {})
        if not isinstance(options, dict) or set(options) - allowed:
            raise ValueError(f'Unknown {key} controls.')
        if 'max_iterations' in options and (isinstance(options['max_iterations'], bool)
                or not isinstance(options['max_iterations'], int) or options['max_iterations'] < 1):
            raise ValueError(f'{key}.max_iterations must be a positive integer.')
        if key == 'casscf_options' and 'max_iterations' in options:
            validate_casscf_iterations(options['max_iterations'])
        for control in ('exact_diagonalization', 'quasi_newton'):
            if control in options and not isinstance(options[control], bool):
                raise ValueError(f'{key}.{control} must be boolean.')
        for control in ('energy_convergence', 'gradient_convergence'):
            if control in options:
                fortran_threshold(options[control], f'{key}.{control}')
        if 'verbosity' in options and (type(options['verbosity']) is not int
                                      or not 0 <= options['verbosity'] <= 4):
            raise ValueError(f'{key}.verbosity must be an integer from 0 to 4.')
    if data.get('scf_options', {}).get('solver', 'diis') not in ('diis', 'second_order', 'so_sci'):
        raise ValueError('scf_options.solver must be diis, second_order or so_sci.')
    if data.get('casscf_options', {}).get('solver', 'second_order') not in ('second_order', 'so_sci'):
        raise ValueError('casscf_options.solver must be second_order or so_sci.')
    if data.get('partition') is not None and not _SAFE.fullmatch(str(data['partition'])):
        raise ValueError('partition must be a safe Slurm partition name.')
    if 'coordinates' in data:
        if 'anchors' in data:
            raise ValueError('Configure coordinates or legacy anchors, not both.')
        data['anchors'] = _coordinates(data, path)
    anchors = data.get('anchors')
    if not isinstance(anchors, list) or not anchors:
        raise ValueError('anchors must be a nonempty list.')
    if data.get('scan_policy') and any(isinstance(a, dict) and
            ('geometry_angstrom' in a or 'orbital_chain' in a) for a in anchors):
        raise ValueError('Explicit orbital chains use a fixed geometry sequence, not adaptive scan_policy.')
    names = set()
    for anchor in anchors:
        allowed = {'name', 'bond', 'moving_atoms', 'distance_angstrom', '_channel', '_context',
                   '_rotation', '_orientation', 'geometry_angstrom', 'orbital_chain'}
        if not isinstance(anchor, dict) or set(anchor) - allowed \
                or not _SAFE.fullmatch(str(anchor.get('name', ''))):
            raise ValueError('Each anchor requires a safe name.')
        if anchor['name'] in names:
            raise ValueError('Anchor names must be unique.')
        names.add(anchor['name'])
        if 'orbital_chain' in anchor and (not isinstance(anchor['orbital_chain'], str)
                                         or not _SAFE.fullmatch(anchor['orbital_chain'])):
            raise ValueError('orbital_chain must be a safe nonempty name.')
        if 'geometry_angstrom' in anchor:
            explicit = np.asarray(anchor['geometry_angstrom'], dtype=float)
            if (explicit.shape != geometry.shape or not np.isfinite(explicit).all()
                    or set(anchor) & {'bond', 'moving_atoms', 'distance_angstrom', '_rotation'}):
                raise ValueError('Explicit anchor geometry must be finite N-by-3 coordinates '
                                 'without a competing displacement or rotation.')
            continue
        if set(anchor) <= {'name', 'orbital_chain', '_channel', '_context', '_orientation'}:
            continue
        if not {'bond', 'moving_atoms', 'distance_angstrom'} <= set(anchor):
            raise ValueError('A displaced anchor requires bond, moving_atoms, '
                             'and distance_angstrom.')
        bond = anchor['bond']
        moving = anchor['moving_atoms']
        distance = anchor['distance_angstrom']
        if (not isinstance(bond, list) or len(bond) != 2
                or not isinstance(moving, list) or not moving
                or any(isinstance(index, bool) or not isinstance(index, int)
                       or index < 1 or index > len(symbols)
                       for index in [*bond, *moving])
                or len(moving) != len(set(moving))
                or bond[0] == bond[1] or (bond[0] in moving) == (bond[1] in moving)
                or isinstance(distance, bool)
                or not isinstance(distance, (int, float))
                or not np.isfinite(distance) or distance <= 0.):
            raise ValueError('Invalid displaced anchor definition.')
    if data.get('probe_active_space'):
        controls = data.get('probe_space', {})
        if (not isinstance(controls, dict)
                or set(controls) - {'mode', 'include_light_atoms', 'additional_centers'}
                or ('mode' in controls and controls['mode'] != 'full_valence')
                or ('include_light_atoms' in controls and not isinstance(controls['include_light_atoms'], bool))):
            raise ValueError('probe_space accepts mode=full_valence, include_light_atoms '
                             'and additional_centers.')
        # Legacy parsing retains the former valence definition for evidence
        # readers. New jobs must explicitly opt in; this is not a compact CAS.
        hydrogen_coordinate = any(any(symbols[index - 1] in ('H', 'He')
                                     for index in anchor.get('bond', [])) for anchor in anchors)
        definition = data.get('probe_definition')
        if definition is None:
            data['_probe_space'] = full_valence_active_space(
                symbols, data['charge'], additional_centers=controls.get('additional_centers'),
                include_light_atoms=controls.get('include_light_atoms', hydrogen_coordinate))
        if definition is not None:
            required_definition = {
                'active_electrons', 'active_orbitals', 'state_model', 'chemical_note'}
            if (not isinstance(definition, dict) or not required_definition <= set(definition)
                    or set(definition) - required_definition - {'orbital_partition', 'orbital_source'}
                    or not isinstance(definition['chemical_note'], str)
                    or not definition['chemical_note'].strip()
                    or any(isinstance(definition[key], bool)
                           or not isinstance(definition[key], int) or definition[key] < 1
                           for key in ('active_electrons', 'active_orbitals'))):
                raise ValueError('probe_definition needs counted orbitals, state_model and chemical_note.')
            model = validate_vrc_state_model(definition['state_model'],
                active_electrons=definition['active_electrons'],
                active_orbitals=definition['active_orbitals'], nelectron=electrons,
                atom_count=len(symbols))
            partition = definition.get('orbital_partition')
            scf_seed = definition.get('orbital_source') == 'scf'
            if 'orbital_source' in definition and not scf_seed:
                raise ValueError('orbital_source must be scf when specified.')
            if scf_seed and (partition is not None or model.get('orbital_selection') is not None
                             or data.get('probe_method') != 'casscf'):
                raise ValueError('SCF seeds require CASSCF and no competing orbital selector/partition.')
            if partition is not None:
                keys = {'source_active_electrons', 'source_active_orbitals',
                        'close_lowest', 'discard_highest'}
                # Explicit model refinement: promote the highest occupied
                # inactive natural orbitals, preserving their checkpoint order.
                # This is not an automatic expansion or a canonical selector.
                opened = partition.get('open_highest_closed', 0) if isinstance(partition, dict) else 0
                if (not isinstance(partition, dict) or not keys <= set(partition)
                        or set(partition) - keys - {'orbital_swaps', 'open_highest_closed'}
                        or any(type(partition[key]) is not int or partition[key] < 0 for key in keys)
                        or type(opened) is not int or opened < 0
                        or (opened and partition['close_lowest'])
                        or partition['source_active_electrons'] < 1
                        or partition['source_active_orbitals'] < 1
                        or partition['source_active_electrons'] > 2 * partition['source_active_orbitals']
                        or partition['source_active_electrons'] - 2 * partition['close_lowest'] + 2 * opened
                            != definition['active_electrons']
                        or partition['source_active_orbitals'] - partition['close_lowest']
                            - partition['discard_highest'] + opened != definition['active_orbitals']
                        or model.get('orbital_selection') is not None
                        or data.get('probe_method') != 'casscf'):
                    raise ValueError('An orbital partition requires consistent source/target counts, '
                                     'optimized CASSCF and no competing orbital selector.')
                swaps = partition.get('orbital_swaps', [])
                if (not isinstance(swaps, list) or any(
                        not isinstance(pair, list) or len(pair) != 2
                        or any(type(i) is not int or not 1 <= i <= partition['source_active_orbitals']
                               for i in pair) or pair[0] == pair[1] for pair in swaps)):
                    raise ValueError('orbital_swaps must contain distinct one-based source-active index pairs.')
                if (electrons - partition['source_active_electrons'] < 0
                        or (electrons - partition['source_active_electrons']) % 2
                        or opened > (electrons - partition['source_active_electrons']) // 2):
                    raise ValueError('Source active electron count is incompatible with the molecule.')
            elif not scf_seed and model.get('orbital_selection', {}).get('mode') != 'avas':
                raise ValueError('A reconnaissance candidate requires explicit AVAS targets.')
            if data['multiplicity'] not in {entry['multiplicity'] for entry in model['state_average']}:
                raise ValueError('Probe multiplicity must occur in its requested state average.')
            data['_probe_model'] = model
            data['_probe_space'] = {key: definition[key] for key in
                                    ('active_electrons', 'active_orbitals', 'chemical_note')}
            if scf_seed:
                from kinbot.recon_scf_seed import validate
                data['_probe_space']['orbital_source'] = 'scf'
                seeds = data.get('initial_records', {})
                for basis in data['bases']:
                    for anchor in anchors:
                        seed = seeds.get(basis['name']+'__'+anchor['name'], {})
                        if seed.get('basis') != basis['basis'] or seed.get('project_orbitals'):
                            raise ValueError('SCF selections require a separate native seed in each matching basis.')
                        validate(seed, electrons, definition['active_electrons'], definition['active_orbitals'])
            elif partition is None:
                data['_probe_space']['orbital_selection'] = model['orbital_selection']
            else:
                data['_probe_space']['orbital_partition'] = partition
                seeds = data.get('initial_records', {})
                if not isinstance(seeds, dict):
                    raise ValueError('Partition initial_records must be a mapping.')
                for basis in data['bases']:
                    seen_chains = set()
                    for anchor in anchors:
                        chain = anchor.get('orbital_chain', anchor.get('_channel'))
                        propagated = (chain and chain in seen_chains
                                      and data.get('orbital_initialization', 'propagated') == 'propagated')
                        if not seeds.get(basis['name'] + '__' + anchor['name']) and not propagated:
                            raise ValueError('Every independent partition or new orbital chain needs '
                                             'its own verified natural-orbital checkpoint in each basis.')
                        if chain:
                            seen_chains.add(chain)
    elif data.get('probe_definition') or data.get('correlation_checks'):
        raise ValueError('Candidate and correlation checks require probe_active_space=true.')
    checks = data.get('correlation_checks', [])
    if not isinstance(checks, list):
        raise ValueError('correlation_checks must be a list.')
    if type(data.get('reference_identity_check', False)) is not bool:
        raise ValueError('reference_identity_check must be boolean.')
    if data.get('reference_identity_check'):
        refs = [c.get('multiplicity') for c in checks
                if isinstance(c, dict) and c.get('method', '').lower() == 'uccsd(t)-f12b']
        if (len(refs) != 1 or refs[0] not in
                {e['multiplicity'] for e in data.get('_probe_model', {}).get('state_average', [])}):
            raise ValueError('Reference identity check requires one CC reference in the MR state average.')
    if checks and data.get('probe_method') != 'casscf':
        raise ValueError('Correlated validation requires optimized CASSCF orbitals.')
    single_reference_seen = False
    for check in checks:
        if not isinstance(check, dict) or not isinstance(check.get('method'), str):
            raise ValueError('Each correlation check requires a method string.')
        if check['method'].lower() == 'uccsd(t)-f12b':
            mult = check.get('multiplicity')
            if (set(check) != {'method', 'multiplicity'} or type(mult) is not int or mult < 1
                    or mult - 1 > electrons or (electrons - mult + 1) % 2):
                raise ValueError('UCCSD(T)-F12b check needs an explicit compatible multiplicity and no MR controls.')
            single_reference_seen = True
            continue
        if single_reference_seen:
            raise ValueError('Single-reference checks must follow all multireference checks.')
        if (not isinstance(check, dict) or set(check) - {'method', 'shift', 'ipea', 'reference_policy'}
                or not vrc_multireference_energy_contract(check.get('method', ''))):
            raise ValueError('correlation_checks requires supported counted CASPT2 or MRCI+Q methods.')
        if check.get('reference_policy', 'selected_root') not in ('selected_root', 'state_average'):
            raise ValueError('Invalid correlated reference_policy.')
        contract = vrc_multireference_energy_contract(check['method'])
        if any(contract[key] != data['_probe_space'][key]
               for key in ('active_electrons', 'active_orbitals')):
            raise ValueError('Correlation method counts must match the probe space.')
        for key in ('shift', 'ipea'):
            if key in check and (isinstance(check[key], bool)
                    or not isinstance(check[key], (float, int)) or not np.isfinite(check[key])):
                raise ValueError('Correlation shifts must be finite numbers.')
    data = deepcopy(data)
    data['_source'] = str(path)
    data['_electrons'] = electrons
    return data


def _anchor_geometry(data, anchor):
    if 'geometry_angstrom' in anchor:
        return np.asarray(anchor['geometry_angstrom'], dtype=float).copy()
    geometry = np.asarray(data['geometry_angstrom'], dtype=float).copy()
    if 'bond' not in anchor:
        return geometry
    first, second = (value - 1 for value in anchor['bond'])
    vector = geometry[second] - geometry[first]
    length = np.linalg.norm(vector)
    if not np.isfinite(length) or length < 1.e-8:
        raise ValueError(f'{anchor["name"]}: reference bond has zero length.')
    shift = vector / length * (float(anchor['distance_angstrom']) - length)
    moving = np.asarray([value - 1 for value in anchor['moving_atoms']], dtype=int)
    geometry[moving] += shift if second in moving else -shift
    rotation = anchor.get('_rotation')
    if rotation:
        axis = np.asarray(rotation['axis'], dtype=float)
        angle = np.deg2rad(rotation['angle_degrees'])
        points = np.asarray(rotation['atoms'], dtype=int)
        origin = geometry[rotation['pivot']].copy()
        vectors = geometry[points] - origin
        geometry[points] = (vectors * np.cos(angle)
            + np.cross(axis, vectors) * np.sin(angle)
            + np.outer(vectors @ axis, axis) * (1. - np.cos(angle)) + origin)
    return geometry


def _probe_model(data, anchor):
    """Keep native energy-array order and reported electronic labels together."""
    if data.get('_probe_model'):
        return deepcopy(data['_probe_model'])
    context = anchor.get('_context', {})
    multiplicities = context.get('multiplicities', [data['multiplicity']])
    multiplicities = [data['multiplicity'], *[mult for mult in multiplicities
                                             if mult != data['multiplicity']]]
    return automatic_vrc_state_model(
        data['symbols'], data['charge'], multiplicities,
        roots=context.get('roots'), active_space=data['_probe_space'])


def _correlation_labels(checks, states, state_model=None):
    labels = [{**check, **state} for check in checks for state in
            ([{'multiplicity': check['multiplicity'], 'root': 1}]
             if check['method'].lower() == 'uccsd(t)-f12b' else states)]
    if (state_model or {}).get('correlation_options', {}).get('state_policy') == 'simultaneous':
        for label in labels:
            if label['method'].lower().startswith('mrci'):
                label['root_convention'] = 'correlated_variational_energy_order'
    return labels


def _validate_cc_check_energies(output, labels, values):
    """Match each F12b marker within its native CC program's A/B summaries.

    UCCSD(T)-F12 prints both approximations, including a triples contribution
    for each. A whole-file triples count therefore is not a count of checks.
    Keep the native variant labels and ordered totals together; never accept
    an F12a marker or an inconsistent total merely because the job exited zero.
    """
    indices = [index for index, label in enumerate(labels, 1)
               if label['method'].lower() == 'uccsd(t)-f12b']
    if not indices:
        return
    from kinbot.anl.results import validate_molpro_triples
    number = r'([-+0-9.]+(?:[EDed][-+]?\d+)?)'
    headers = list(re.finditer(r'(?im)^\s*PROGRAM\s*\*\s*([^\r\n]+)', output))
    native = []
    for pos, header in enumerate(headers):
        if not re.match(r'CCSD\b', header[1], re.I):
            continue
        block = output[header.end():headers[pos+1].start() if pos+1 < len(headers) else len(output)]
        bases = re.findall(r'(?im)^\s*!?RHF-UCCSD-F12([ab]) energy\s+' + number + r'\s*$', block)
        if not bases:
            continue
        triples = validate_molpro_triples(block, scaled=False)
        totals = re.findall(r'(?im)^\s*!RHF-UCCSD\(T\)-F12(?: total)? energy\s+' + number + r'\s*$', block)
        if ([variant.lower() for variant, _ in bases] not in (['b'], ['a', 'b'])
                or len(bases) != len(triples) or len(bases) != len(totals)):
            raise ValueError('Single-reference checks lack uniquely matched native F12a/F12b/(T) summaries.')
        for (variant, base), trip, total in zip(bases, triples, totals):
            base, total = [float(v.replace('D', 'E').replace('d', 'e')) for v in (base, total)]
            if not np.isfinite([base, trip, total]).all() or abs(total - base - trip) > 2e-9:
                raise ValueError('Single-reference native F12 total differs from its CCSD plus (T) energy.')
            if variant.lower() == 'b':
                native.append(total)
    if len(native) != len(indices):
        raise ValueError('Single-reference checks lack uniquely matched native F12b/(T) energies.')
    if any(not np.isfinite(values[i]) or abs(values[i] - total) > 1e-7
           for i, total in zip(indices, native)):
        raise ValueError('Single-reference marker differs from the native F12b plus (T) energy.')


def _recon_rhf(command, electrons, spin, record):
    """Specify the C1 determinant independently of inherited AVAS/MULTI state.

    Restart files retain global SPIN as well as orbital occupations. Supply
    both global and local spin, and the matching determinant partition.
    Local WF cards in subsequent MULTI/CI blocks retain precedence.
    """
    if spin < 0 or electrons < spin or (electrons - spin) % 2:
        raise ValueError('Incompatible RHF electron count and spin.')
    return (f'set,spin={spin}\n{{{command};wf,{electrons},1,{spin};'
            f'occ,{(electrons + spin) // 2};closed,{(electrons - spin) // 2};'
            f'orbital,{record}}}\n')


def _deck(data, basis, anchor, stem, previous=None):
    geometry = _anchor_geometry(data, anchor)
    rows = '\n'.join(
        f'{symbol} {x:.12f} {y:.12f} {z:.12f}'
        for symbol, (x, y, z) in zip(data['symbols'], geometry))
    scf = fortran_threshold(data.get('scf_convergence', 1.e-10),
                            'scf_convergence')
    coeff = fortran_threshold(data.get('correlation_convergence', 1.e-5), 'correlation_convergence')
    thresholds = (f'gthresh,energy={scf},coeff={coeff},civec={coeff}\n'
                  if 'correlation_convergence' in data else f'gthresh,energy={scf},coeff=1.d-5\n')
    scf_options = data.get('scf_options', {})
    scf_command = 'rhf' + {'diis': '', 'second_order': ',so', 'so_sci': ',so-sci'}[
        scf_options.get('solver', 'diis')]
    if 'max_iterations' in scf_options:
        scf_command += f',maxit={scf_options["max_iterations"]}'
    basis_input = f'basis={basis}\n'
    printing = 'gprint,orbitals,orben,civector\n'
    if 'orbital_print_virtuals' in data:
        printing = f'gprint,orbitals={data["orbital_print_virtuals"]},orben,civector\n'
    if 'orbital_print_threshold' in data:
        thresholds += f'gthresh,thrprint={data["orbital_print_threshold"]}\n'
    deck = (
        f'***,{data["title"]} orbital reconnaissance: {stem}\n'
        f'memory,{int(data.get("memory_mw", 1000))},m\n'
        f'file,2,{stem}.wfu\n'
        'angstrom\norient,noorient\nsymmetry,nosym\n'
        f'geometry={{\n{len(data["symbols"])}\n{stem}\n{rows}\n}}\n'
        + basis_input + thresholds + printing
        + _recon_rhf(scf_command, data['_electrons'], data['multiplicity'] - 1, '2100.2')
        + printing +
        f'put,molden,{stem}.molden\n'
        'reconenergy=energy\nshow,reconenergy\n')
    if data.get('probe_active_space', False):
        space = data['_probe_space']
        model = _probe_model(data, anchor)
        # The recon renderer adds its natural exports and local controls below.
        # Inherit reviewed production controls without rendering them twice.
        inherited_options = model.pop('casscf_options', {})
        if previous is not None:
            model.pop('orbital_selection', None)
        avas, multi, _ = _render_vrc_state_model(
            model, nelectron=data['_electrons'],
            active_electrons=space['active_electrons'],
            active_orbitals=space['active_orbitals'], spin=data['multiplicity'] - 1)
        method = data.get('probe_method', 'casci')
        opening = '{multi,cas-ci;' if method == 'casci' else '{multi,'
        multi = multi.replace('{multi,', opening)
        # Always emit natural orbitals, including when canonical_orbitals=false.
        multi = multi.replace('canonical,2141.2\n', '')
        if data.get('correlation_checks'):
            # A reused WFU can contain canonical orbitals for an older/larger
            # space. Save this optimization separately and select it explicitly
            # for CASPT2 (especially important with the non-invariant IPEA shift).
            multi = multi.replace('natorb,2140.2\n',
                                  'canonical,2141.2\nnatorb,2140.2\n')
        options = {**inherited_options, **data.get('casscf_options', {})}
        if method == 'casscf':
            # SO-SCI derives a 100x tighter energy target from GTHRESH than
            # second-order MULTI. Set the intended local target explicitly;
            # changing the optimizer must not silently change convergence.
            accuracy = 'energy=' + fortran_threshold(
                options.get('energy_convergence', data.get('scf_convergence', 1.e-10)),
                'casscf_options.energy_convergence')
            if 'gradient_convergence' in options:
                accuracy += ',gradient=' + fortran_threshold(
                    options['gradient_convergence'], 'casscf_options.gradient_convergence')
            multi = multi.replace('natorb,2140.2\n', f'accuracy,{accuracy}\nnatorb,2140.2\n')
            controls = []
            if 'solver' in options:
                controls.append({'second_order': 'so', 'so_sci': 'so-sci'}[options['solver']])
            for key, option in (('exact_diagonalization', 'exact_diag'), ('quasi_newton', 'qn')):
                if key in options:
                    controls.append(f'{option}={int(options[key])}')
            if 'verbosity' in options:
                controls.append(f'verbosity={options["verbosity"]}')
            if 'max_iterations' in options:
                multi = multi.replace('maxit,50;', f'maxit,{options["max_iterations"]};')
        if avas and data.get('avas_reference') == 'rhf':
            avas = avas.replace('}\n', 'start,2100.2\n}\n', 1)
        if previous is not None:
            start = ('start,2140.2,type=natural' if method == 'casscf'
                     else 'start,2140.2')
            if isinstance(previous, dict) and previous.get('project_orbitals'):
                # An explicitly declared cross-basis seed is a projected guess,
                # not a reviewed wavefunction in the destination basis.
                start += ',ovl'
            # A source checkpoint may need a one-time partition/reordering.
            # A propagated predecessor already contains the target CAS.
            initial_partition = bool(space.get('orbital_partition') and isinstance(previous, dict))
            diab = '' if initial_partition else '\ndiab,2140.2'
            partition = space.get('orbital_partition', {})
            if space.get('orbital_source') == 'scf':
                from kinbot.recon_scf_seed import start as scf_start
                if not isinstance(previous, dict):
                    raise ValueError('SCF source indices require an explicit initial record.')
                start = scf_start(previous, data['_electrons'], space['active_electrons'])
                diab = ''
            if initial_partition and partition.get('orbital_swaps'):
                source_closed = (data['_electrons'] - partition['source_active_electrons']) // 2
                # Reorder the verified source before closing/discarding its
                # ends. Occupation rank alone does not identify sigma/pi MOs.
                # MULTI documents angle=0 as its explicit exchange operation.
                start += ''.join(f'\nrotate,{source_closed + a}.1,{source_closed + b}.1,0'
                                 for a, b in partition['orbital_swaps'])
            multi = multi.replace(opening, opening + f'\n{start}{diab}', 1)
            deck = '! Initial orbital checkpoint SHA-256: {{PREVIOUS_WFU_SHA256}}\n' + deck
        if method == 'casscf' and controls:
            multi = multi.replace('{multi,', '{multi,' + ','.join(controls) + ';', 1)
        if data.get('reference_identity_check'):
            from kinbot.molpro_reference_identity import natural_cards
            reference = next(c['multiplicity'] for c in data['correlation_checks']
                             if c['method'].lower() == 'uccsd(t)-f12b')
            # Diagnostic saves only: keep the original SA natural/canonical
            # records last, and do not change optimization roots or weights.
            multi = multi.replace('canonical,2141.2\n',
                natural_cards(model, reference)+'canonical,2141.2\n')
        deck += (avas + multi + f'{{put,molden,{stem}.natural.molden;orbital,2140.2,type=natural}}\n'
                 'reconenergy=energy(1)\nshow,reconenergy\n')
        for index in range(1, sum(entry['roots'] for entry in model['state_average']) + 1):
            deck += f'reconr{index}=energy({index})\nshow,reconr{index}\n'
        if data.get('reference_identity_check'):
            from kinbot.molpro_reference_identity import export_deck
            deck += export_deck(stem, model, reference, part='mr')
        # Reuse the production selected-state energy contract, including
        # Davidson-corrected rather than variational MRCI energies.
        counter = 0
        for check in data.get('correlation_checks', []):
            if check['method'].lower() == 'uccsd(t)-f12b':
                counter += 1
                # Fresh restricted reference for this explicitly requested spin.
                # Preserve the reviewed MULTI orbital records for later selection.
                deck += _recon_rhf(scf_command, data['_electrons'], check['multiplicity'] - 1, '2150.2')
                if data.get('reference_identity_check'):
                    from kinbot.molpro_reference_identity import export_deck
                    deck += (f'put,molden,{stem}.reference_c{counter}_m{check["multiplicity"]}.molden\n'
                             + export_deck(stem, model, check['multiplicity'], part='hf'))
                if data.get('reference_identity_check'):
                    from kinbot.theory import molpro_cc_thresholds
                    den, var = molpro_cc_thresholds(data.get('correlation_convergence', 1.e-7))
                    deck += (f'{{uccsd(t)-f12b,thrden={den},thrvar={var}}}\n'
                             f'reconc{counter}=energy\nshow,reconc{counter}\n')
                else:
                    deck += ('{uccsd(t)-f12}\n'
                             f'reconc{counter}=energy(2)\nshow,reconc{counter}\n')
                continue
            for entry in model['state_average']:
                simultaneous = (model.get('correlation_options', {}).get('state_policy')
                                == 'simultaneous')
                # One common eigensolve per spin protects the requested roots
                # from collapse into separately optimized lower solutions.
                roots = [1] if simultaneous else range(1, entry['roots'] + 1)
                for root in roots:
                    selected = deepcopy(model)
                    # Reconnaissance records each labelled root. A production
                    # minimum selector must not relabel a minimum as root 1.
                    selected.pop('energy_selection', None)
                    if 'reference_policy' in check:
                        selected.setdefault('correlation_options', {})['reference_policy'] = check['reference_policy']
                    selected['selected_roots'][entry['multiplicity']] = root
                    avas_prefix, multi_prefix, _ = _render_vrc_state_model(selected,
                        nelectron=data['_electrons'], active_electrons=space['active_electrons'],
                        active_orbitals=space['active_orbitals'], spin=entry['multiplicity'] - 1)
                    body, variable = _vrc_multireference_method(check['method'], data['_electrons'],
                        entry['multiplicity'] - 1, state_model=selected,
                        shift=check.get('shift'), ipea=check.get('ipea'))
                    prefix = avas_prefix + multi_prefix
                    if not body.startswith(prefix):
                        raise ValueError('Unexpected production multireference renderer prefix.')
                    counter += 1
                    correlation = body[len(prefix):]
                    orbital = ('2141.2,type=canonical' if check['method'].lower().startswith('caspt2')
                               else '2140.2,type=natural')
                    correlation = re.sub(r';orbital,[^;]+;', ';', correlation, count=1)
                    correlation = correlation.replace(';wf,', f';orbital,{orbital};wf,', 1)
                    deck += correlation + f'reconc{counter}={variable}\nshow,reconc{counter}\n'
                    if simultaneous:
                        for slot in range(2, entry['roots'] + 1):
                            counter += 1
                            deck += (f'! Simultaneous MRCI energy-ordered state {slot}; relaxed Davidson\n'
                                     f'KBVRCVAR=energy({slot})\nKBVRCSEL=energd({slot})\n'
                                     f'reconc{counter}=KBVRCSEL\nshow,reconc{counter}\n')
    return deck + 'recondone=1\nshow,recondone\n'


def _repair_legacy_energy_selectors(text):
    """Correct only the old generated MR contract in an immutable input copy.

    CCSD-F12's A/B indexing never applied to RS2-F12/MRCI-F12. Keep prepared
    inputs and their hashes intact; execution records bind this native copy.
    """
    pattern = (r'(?im)^! KinBot energy contract: (CASPT2|MRCI)-F12B(\+Q)? = '
               r'(ENERGY|ENERGD)\(2\)(; Davidson relaxed reference)?\n'
               r'(?:KBVRCVAR=energy\(2\)\n)?KBVRCSEL=(?:energy|energd)\(2\)\n')

    def replace(match):
        mrci = match[1].upper() == 'MRCI'
        if (mrci != bool(match[2]) or mrci != (match[3].upper() == 'ENERGD')
                or mrci != bool(match[4])):
            raise ValueError('Inconsistent legacy multireference energy contract.')
        variable = 'energd' if mrci else 'energy'
        method = 'MRCI-F12+Q' if mrci else 'CASPT2-F12'
        return (f'! KinBot energy contract: {method} = {variable.upper()}(1)'
                + ('; Davidson relaxed reference' if mrci else '') + '\n'
                + ('KBVRCVAR=energy(1)\n' if mrci else '')
                + f'KBVRCSEL={variable}(1)\n')

    return re.subn(pattern, replace, text)


def _recover_legacy_caspt2_energies(output, native_input, labels, electron_count):
    """Recover printed CASPT2 states only for the exact historical selector bug.

    Callers still require successful execution, immutable input hashes,
    convergence, exports and checkpoint evidence. Never modify native output
    or use an unlabelled final energy as a substitute for a requested root.
    """
    from kinbot.molpro import vrc_multireference_energy_contract
    if not labels or any(not label['method'].lower().startswith('caspt2-f12') for label in labels):
        raise ValueError('Legacy recovery supports only labelled RS2-F12 results.')
    _, repairs = _repair_legacy_energy_selectors(native_input)
    if repairs != len(labels):
        raise ValueError('Not the known generated CASPT2 energy-selector defect.')
    undefined = re.findall(r'(?im)^\s*VARIABLE\s+(\w+)\s+UNDEFINED,\s*ASSUMING\s+0', output)
    if undefined != ['ENERGY'] * len(labels):
        raise ValueError('Unexpected undefined variables in legacy output.')
    commands = re.findall(r'(?is)\{(rs2-f12\b[^{}]*)\}', native_input)
    sections = re.split(r'(?im)^\s*PROGRAM\s*\*\s*', output)[1:]
    blocks = [section for section in sections if re.match(r'RS2-F12\b', section, re.I)]
    if len(commands) != len(labels) or len(blocks) != len(labels):
        raise ValueError('Native CASPT2 invocation count differs from the requested roots.')
    spins = ('Singlet', 'Doublet', 'Triplet', 'Quartet', 'Quintet', 'Sextet',
             'Septet', 'Octet', 'Nonet', 'Decet')
    values = {}
    number = r'([-+0-9.]+(?:[EDed][-+]?\d+)?)'
    for index, (command, block, label) in enumerate(zip(commands, blocks, labels), 1):
        contract = vrc_multireference_energy_contract(label['method'])
        mult, root = label['multiplicity'], label['root']
        wf = re.search(r'(?i)\bwf,(\d+),1,(\d+)\s*;', command)
        state = re.search(r'(?i)\bstate,1,(\d+)\s*;', command)
        energy = re.findall(r'(?m)^\s*!RSPT2 STATE\s+(\d+)\.1\s+Energy\s+' + number + r'\s*$', block)
        roots = re.search(r'Number of optimized states:\s*1\s+Roots:\s*(\d+)\s*\n', block)
        electrons = re.search(r'(?m)^\s*Number of electrons:\s*(\d+)\s*$', block)
        active = re.search(r'Number of active\s+orbitals:\s*(\d+)', block)
        core = re.search(r'Number of core orbitals:\s*(\d+)', block)
        closed = re.search(r'Number of closed-shell orbitals:\s*(\d+)', block)
        if (not wf or tuple(map(int, wf.groups())) != (electron_count, mult - 1)
                or not state or int(state[1]) != root or not roots or int(roots[1]) != root
                or not electrons or int(electrons[1]) != electron_count
                or not active or int(active[1]) != contract['active_orbitals']
                # Molpro omits a zero-sized core/closed block entirely.
                or electron_count - 2 * ((int(core[1]) if core else 0)
                                        + (int(closed[1]) if closed else 0)) != contract['active_electrons']
                or not 1 <= mult <= len(spins)
                or not re.search(r'Reference symmetry:\s*1\s+' + spins[mult - 1] + r'\b', block)
                or len(energy) != 1 or int(energy[0][0]) != root):
            raise ValueError(f'Native CASPT2 state identity mismatch at check {index}.')
        for option, heading, default in (('shift', 'Level shift', .3), ('ipea', 'IPEA shift', 0.)):
            requested = float(label.get(option, default))
            specified = re.search(r'(?i)\b' + option + r'=' + number, command)
            printed = re.search(heading + r'\s*=\s*' + number, block)
            for found in (specified, printed):
                actual = float(found[1].replace('D', 'E').replace('d', 'e')) if found else 0.
                if not np.isclose(actual, requested, rtol=0., atol=1e-10):
                    raise ValueError(f'Native CASPT2 {option} mismatch at check {index}.')
        value = float(energy[0][1].replace('D', 'E').replace('d', 'e'))
        if not np.isfinite(value):
            raise ValueError('Nonfinite native CASPT2 energy.')
        values[index] = value
    return values


def prepare(source, output):
    require_explicit_probe_space(json.loads(Path(source).expanduser().read_text()))
    with _manifest_lock(output):
        return _prepare(source, output)


def _prepare(source, output):
    data = _load(source)
    require_explicit_probe_space(data)
    root = Path(output).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    setup = root / 'site_setup.sh'
    if not setup.is_file():
        setup.write_text(render_site_setup({'molpro': {'molpro'}}))
    ppn = int(data.get('ppn', 8))
    partition = data.get('partition')
    # Molpro's MEMORY is a stack per rank, not the process's total RSS.
    # Its parallel manual estimates another 300 MW/rank for the program;
    # reserve 25% beyond that for auxiliary allocations and I/O buffering.
    # This is a starting budget, not an upper bound on every MR calculation.
    memory_mb = int(data.get('memory_mb', ppn * (data.get('memory_mw', 1000) + 300) * 8 * 1.25))
    walltime = str(data.get('walltime', 'partition_max'))
    if memory_mb < ppn * data.get('memory_mw', 1000) * 8:
        raise ValueError('Slurm memory_mb must cover every Molpro rank memory stack.')
    if walltime != 'partition_max' and not re.fullmatch(r'(?:\d+-)?\d+:\d{2}:\d{2}', walltime):
        raise ValueError('Reconnaissance walltime must be partition_max or a Slurm time string.')
    # Resolve at submission: offline preparation must not guess site limits.
    time_directive = ('# KINBOT_WALLTIME=partition_max\n' if walltime == 'partition_max'
                      else f'#SBATCH --time={walltime}\n')
    records = []
    manifest_file = root / 'reconnaissance.json'
    previous = (json.loads(manifest_file.read_text()) if manifest_file.is_file() else {})
    # A stopped campaign may explicitly exclude costly basis families while
    # retaining its immutable source configuration and all native evidence.
    excluded_bases = previous.get('excluded_bases', [])
    data['bases'] = [basis for basis in data['bases'] if basis['basis'] not in excluded_bases]
    if not data['bases']:
        raise ValueError('Basis selection would leave no reconnaissance jobs.')
    previous_jobs = {job['name']: job for job in previous.get('jobs', [])}
    predecessors = {}
    for basis in data['bases']:
        for anchor in data['anchors']:
            stem = f'{basis["name"]}__{anchor["name"]}'
            directory = root / stem
            directory.mkdir(exist_ok=True)
            inp = directory / f'{stem}.inp'
            channel = anchor.get('orbital_chain', anchor.get('_channel'))
            old_job = previous_jobs.get(stem, {})
            if old_job and (not inp.is_file() or file_sha256(inp)
                            != old_job['input_sha256']):
                raise ValueError(f'{stem}: existing probe input changed after preparation.')
            # Refinement inserts points into the distance ordering. Existing
            # probes must retain the checkpoint they actually consumed.
            predecessor = (old_job.get('orbital_predecessor') if old_job else
                           predecessors.get((basis['name'], channel))
                           if channel and data.get('probe_active_space')
                           and data.get('orbital_initialization', 'propagated') == 'propagated' else None)
            seed = data.get('initial_records', {}).get(stem)
            if data.get('_probe_space', {}).get('orbital_partition') and not (seed or predecessor):
                raise ValueError('A partitioned candidate needs a verified seed or a propagated predecessor.')
            if seed:
                scf_fields = {'active_indices', 'source_output', 'source_output_sha256',
                              'source_input', 'source_input_sha256'}
                if (isinstance(seed, dict) and set(seed) & scf_fields
                        and (data.get('_probe_space') or {}).get('orbital_source') != 'scf'):
                    raise ValueError('Explicit SCF seed evidence requires orbital_source=scf.')
                if (not isinstance(seed, dict) or not {'path', 'sha256'} <= set(seed)
                        or set(seed) - {'path', 'sha256', 'basis', 'project_orbitals',
                                       'active_indices', 'source_output', 'source_output_sha256',
                                       'source_input', 'source_input_sha256'}
                        or ('project_orbitals' in seed and type(seed['project_orbitals']) is not bool)
                        or ('basis' in seed and (not isinstance(seed['basis'], str)
                                                or not seed['basis'].strip()))
                        or (seed.get('project_orbitals') and not seed.get('basis'))
                        or (seed.get('basis') and seed['basis'].lower() != basis['basis'].lower()
                            and not seed.get('project_orbitals'))
                        or data.get('probe_method') != 'casscf'
                        or not Path(seed['path']).is_file() or file_sha256(seed['path']) != seed['sha256']):
                    raise ValueError('Initial CASSCF record requires a verified checkpoint; '
                                     'a different declared basis requires project_orbitals=true '
                                     'and a declared source basis.')
                predecessor = None
            deck = _deck(data, basis['basis'], anchor, stem, previous=seed or predecessor)
            if seed:
                deck = deck.replace('{{PREVIOUS_WFU_SHA256}}', seed['sha256'])
            validate_input_blocks(deck, stem)
            digest = hashlib.sha256(deck.encode()).hexdigest()
            if old_job and old_job.get('template_sha256', old_job['input_sha256']) != digest:
                raise ValueError(f'{stem}: probe definition changed; choose a '
                                 'fresh output directory to preserve prior work.')
            if not old_job:
                inp.write_text(deck)
            sbatch = directory / f'{stem}.sbatch'
            sbatch.write_text(
                '#!/bin/bash\n'
                f'#SBATCH -N 1\n#SBATCH --exclusive\n#SBATCH -n {ppn}\n#SBATCH -c 1\n'
                + (f'#SBATCH -p {partition}\n' if partition else '')
                + f'#SBATCH --mem={memory_mb}M\n' + time_directive +
                f'#SBATCH -J kb-orb-{stem[:40]}\n'
                f'#SBATCH -o {stem}.stdout\n#SBATCH -e {stem}.stderr\n'
                'set -euo pipefail\n'
                'export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\n'
                'export KINBOT_BACKEND=molpro\n'
                'set +u\n'
                f'source {shlex.quote(str((root / "site_setup.sh").resolve()))}\n'
                'set -u\n'
                'command -v molpro >/dev/null || { echo "Molpro unavailable after site setup" >&2; exit 127; }\n'
                f'exec {shlex.quote(sys.executable)} -m kinbot.molpro_reconnaissance execute '
                f'{shlex.quote(str(manifest_file))} {shlex.quote(stem)}\n')
            sbatch.chmod(0o700)
            from kinbot.slurm_retry import wrap_script
            sbatch.write_text(wrap_script(sbatch.read_text(), identity=stem,
                                         directories=[str(directory)]))
            records.append({
                **old_job,
                'name': stem, 'basis': basis['basis'], 'anchor': anchor['name'],
                'directory': str(directory),
                'input_sha256': file_sha256(inp),
                'template_sha256': digest,
                'orbital_predecessor': predecessor,
                'initial_record': seed,
                'channel': channel,
                'orientation': anchor.get('_orientation'),
                'distance_angstrom': anchor.get('distance_angstrom'),
                'spin_context': anchor.get('_context'),
                'root_states': ([{'multiplicity': entry['multiplicity'], 'root': root}
                    for entry in _probe_model(data, anchor)['state_average']
                    for root in range(1, entry['roots'] + 1)]
                    if data.get('probe_active_space') else []),
                'natural_orbital_probe': bool(data.get('probe_active_space', False)),
                'electron_count': data['_electrons'],
                'probe_active_space': data.get('_probe_space'),
                'probe_method': data.get('probe_method', 'casci'),
                'avas_reference': data.get('avas_reference', 'atomic'),
                'geometry_angstrom': _anchor_geometry(data, anchor).tolist(),
                'state_model': (_probe_model(data, anchor) if data.get('probe_active_space') else None),
                'correlation_checks': data.get('correlation_checks', []),
                'basis_source': 'library',
                'symbols': data['symbols'],
                'ppn': ppn,
                'scratch_min_mb': data.get('scratch_min_mb'),
                **{k:data[k] for k in ('reference_review_contract_sha256',
                                      'reference_review_distance_angstrom') if k in data},
            })
            if channel:
                predecessors[(basis['name'], channel)] = stem
    manifest = {
        **{key: previous[key] for key in ('excluded_bases', 'excluded_jobs', 'slurm_excluded_nodes',
                                                         'slurm_inherited_exclusion_allow', 'startup_retries',
                                                         'reference_solver_options', 'solver_options', 'shared_policy_file') if key in previous},
        'schema': 1, 'source': data['_source'], 'title': data['title'],
        'electron_count': data['_electrons'], 'jobs': records,
        'scan_policy': data.get('scan_policy'),
        **({'orbital_character_only': True} if data.get('orbital_character_only') else {}),
        'status': 'prepared',
    }
    path = root / 'reconnaissance.json'
    atomic_json(path, manifest)
    (root / 'submit.sh').write_text(
        '#!/bin/bash\nset -euo pipefail\n'
        f'exec {shlex.quote(sys.executable)} -m kinbot.molpro_reconnaissance run '
        f'{shlex.quote(str(path))} --max-jobs {int(data.get("queue_job_limit", 8))}\n')
    (root / 'submit.sh').chmod(0o700)
    return path


def execute(manifest_path, name):
    """Run one allocated probe with the same isolated runtime as production."""
    manifest = json.loads(Path(manifest_path).read_text())
    job = next((job for job in manifest['jobs'] if job['name'] == name), None)
    if job is None:
        raise ValueError(f'Unknown reconnaissance job: {name}.')
    # Publish before staging seeds or validating inputs: failures in those
    # steps previously left no execution record and looked like missing CAS
    # occupations/energies. The native worker fills in the remaining metadata.
    record = Path(job['directory']) / 'execution.json'
    execution = {'schema': 1, 'status': 'running', 'host': socket.gethostname(),
                 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
                 'input_sha256': job['input_sha256']}
    atomic_json(record, execution)
    started = time.monotonic()
    try:
        _execute(job)
    except BaseException as error:
        execution = json.loads(record.read_text())
        execution.update(status='failed', error=str(error), error_type=type(error).__name__,
                         elapsed_seconds=time.monotonic() - started)
        atomic_json(record, execution)
        raise


def _reference_solver_override(text, vectors=None, iterations=None):
    """Change only explicit MR reference-CI limits in hash-bound native inputs."""
    controls = {key: value for key, value in (
        ('max_reference_vectors', vectors), ('max_reference_iterations', iterations))
        if value is not None}
    if not controls:
        return text
    validate_mr_correlation_options(controls)
    def change(match):
        body = match[0][:-1].rstrip()
        for value, option, card in ((vectors, 'maxvi', 'MAXDAV'),
                                     (iterations, 'maxiti', 'MAXITER')):
            if value is None:
                continue
            # Generated decks use OPTION; refuse competing foreign cards.
            if re.search(r'(?:[;\n])\s*' + card + r'\s*,', body, re.I):
                raise ValueError(f'Cannot override an explicit {card} card; prepare a reviewed input.')
            if re.search(r'\b' + option + r'\s*=', body, re.I):
                body = re.sub(r'\b' + option + r'\s*=\s*[^,;\s}]+',
                              f'{option}={value}', body, flags=re.I)
            else:
                body += ('' if body.endswith(';') else ';') + f'option,{option}={value};'
        return body + '}'
    result, count = re.subn(r'\{\s*(?:rs2(?:c|-f12)?|mrci(?:-f12)?)(?=[,;\s}])[^{}]*\}',
                            change, text, flags=re.I)
    if not count:
        raise ValueError('No supported MR correlation block for reference solver override.')
    return result


def _separate_cross_basis_seed(text, seed_name):
    """Project a foreign-basis seed from FILE,3 into a fresh FILE,2.

    Source natural orbitals stay at 2140.3; the established destination
    records (2140.2/2141.2) remain valid for exports and subsequent jobs.
    Apply to a native copy, preserving old prepared inputs and fingerprints.
    """
    if not _SAFE.fullmatch(seed_name) or len(seed_name) > 31:
        raise ValueError('Unsafe cross-basis seed filename.')
    if re.search(r'(?im)^\s*file\s*,\s*3\s*,', text):
        raise ValueError('Cross-basis transformation requires unused FILE,3.')
    text, count = re.subn(r'(?im)^(file,2,[A-Za-z0-9_.-]+)\s*$',
                          lambda m: m[1] + ',new\nfile,3,' + seed_name, text)
    if count != 1:
        raise ValueError('Cross-basis transformation requires one destination FILE,2.')
    text, count = re.subn(r'(?im)^(\s*start,)2140\.2(,type=natural,ovl\s*)$',
                          r'\g<1>2140.3\2', text)
    if count != 1:
        raise ValueError('Cross-basis transformation requires one projected natural START.')
    # With a fresh FILE,2 there is no implicit SCF restart. Give the SCF
    # guess the same explicit seed, without changing its final HF equations.
    text, count = re.subn(r'(?i)(\{rhf(?:,[^;{}]*)?;)(wf,)',
                          r'\1start,2140.3,type=natural,ovl;\2', text, count=1)
    if count != 1:
        raise ValueError('Cross-basis transformation requires the generated initial RHF block.')
    return text


def _execute(job):
    name = job['name']
    directory = Path(job['directory'])
    seed = job.get('initial_record') or {}
    cross_basis = bool(seed.get('project_orbitals') and
                       seed.get('basis', '').lower() != job['basis'].lower())
    scf_seed = (job.get('probe_active_space') or {}).get('orbital_source') == 'scf'
    separate_seed = cross_basis or scf_seed
    seed_name = molpro_wfu_name(name + '_source') if separate_seed else None
    if scf_seed:
        from kinbot.recon_scf_seed import validate
        space = job['probe_active_space']
        if seed.get('basis', '').lower() != job['basis'].lower():
            raise ValueError('SCF seed coefficient restart requires the same basis.')
        validate(seed, job['electron_count'], space['active_electrons'], space['active_orbitals'])
    if job.get('initial_record'):
        import shutil
        seed = job['initial_record']
        if file_sha256(seed['path']) != seed['sha256']:
            raise ValueError('Initial CASSCF checkpoint changed before execution.')
        # Keep foreign-basis records off FILE,2. MULTI cannot overwrite a
        # dump record with another basis dimension, even after START,OVL.
        shutil.copyfile(seed['path'], directory / (seed_name or f'{name}.wfu'))
    inp = directory / f'{name}.inp'
    if file_sha256(inp) != job['input_sha256']:
        raise ValueError(f'{name}: input changed before native execution.')
    validate_input_blocks(inp.read_text(), name)
    # Molden cannot encode h and higher functions in Molpro 2024. Use the
    # native XML exporter for every basis, without guessing a library's lmax.
    # Keep the prepared input immutable so running/finished campaigns and
    # refinement checkpoints retain their original fingerprints.
    native = directory / f'{name}.native.inp'
    native_text = re.sub(r'(?i)put,molden,([^;\s}]+)\.molden',
                         r'put,xml,\1.orbitals.xml', inp.read_text())
    native_text = compact_molpro_titles(native_text)
    native_text, energy_repairs = _repair_legacy_energy_selectors(native_text)
    from kinbot.molpro_reference_identity import repair_legacy_selectors
    native_text, identity_repairs = repair_legacy_selectors(native_text)
    vectors = job.get('retry_max_reference_vectors')
    iterations = job.get('retry_max_reference_iterations')
    native_text = _reference_solver_override(native_text, vectors, iterations)
    native_text = apply_solver_options(native_text, job.get('retry_solver_options', {}))
    checkpoint = directory / molpro_wfu_name(name)
    native_text, changed = re.subn(r'(?im)^file,2,[^\n]+$',
                                  f'file,2,{checkpoint.name}', native_text)
    if changed != 1:
        raise ValueError('Expected exactly one reconnaissance FILE,2 directive.')
    if cross_basis:
        native_text = _separate_cross_basis_seed(native_text, seed_name)
    elif scf_seed:
        from kinbot.recon_scf_seed import separate
        native_text = separate(native_text, seed_name)
    recovery = job.get('recovered_scf')
    if recovery:
        source = Path(recovery['path'])
        staged = directory / f'{name}.wfu'
        if (recovery['input_sha256'] != job['input_sha256']
                or not source.is_file() or not staged.is_file()
                or file_sha256(source) != recovery['sha256']
                or file_sha256(staged) != recovery['sha256']):
            raise ValueError('Recovered SCF record no longer matches its archived source.')
        native_text, changed = re.subn(r'(?i)\{rhf;[^{}]*\}',
            '! Reusing converged SCF record from the archived export-only failure',
            native_text, count=1)
        if changed != 1:
            raise ValueError('Cannot identify the original SCF block for recovery.')
        native_text = native_text.replace('reconenergy=energy\n',
            f'reconenergy={float(recovery["energy_hartree"])!r}\n', 1)
        # Explicitly select the existing canonical record for the SCF export.
        native_text = native_text.replace(f'put,xml,{name}.orbitals.xml',
            f'{{put,xml,{name}.orbitals.xml;orbital,{recovery["record"]},type=canonical}}')
    validate_input_blocks(native_text, name)
    native.write_text(native_text)
    # Submission stages predecessor/export-recovery records under the old
    # canonical path. Copy to the short native name before launching Molpro.
    staged = directory / f'{name}.wfu'
    if not separate_seed and staged != checkpoint and staged.is_file():
        import shutil
        shutil.copyfile(staged, checkpoint)
    exports = {'scf': str(directory / f'{name}.orbitals.xml')}
    if job.get('natural_orbital_probe'):
        exports['natural'] = str(directory / f'{name}.natural.orbitals.xml')
    runtime = {}
    execution = {'schema': 1, 'status': 'running',
                 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
                 'host': socket.gethostname(),
                 'input_sha256': job['input_sha256'],
                 'native_input': str(native),
                 'native_input_sha256': file_sha256(native),
                 'input_transform': 'native_xml_mr_state_energy_v4', 'orbital_exports': exports,
                 'legacy_mr_energy_selectors_repaired': energy_repairs}
    if identity_repairs:
        execution['reference_density_selectors_repaired'] = identity_repairs
    for key in ('retry_max_reference_vectors', 'retry_max_reference_iterations'):
        if job.get(key) is not None:
            execution[key] = job[key]
    if job.get('retry_solver_options'):
        execution['retry_solver_options'] = deepcopy(job['retry_solver_options'])
        execution['input_transform'] = 'native_xml_solver_controls_v5'
    if cross_basis:
        execution['input_transform'] = 'native_xml_cross_basis_file_v6'
        execution['projected_seed'] = {
            'path': str(directory / seed_name), 'sha256': seed['sha256'],
            'source_basis': seed['basis'], 'destination_basis': job['basis'],
            'source_record': '2140.3', 'natural_destination_record': '2140.2',
            'canonical_destination_record': '2141.2'}
    if scf_seed:
        execution['input_transform'] = 'native_xml_scf_seed_same_basis_v8'
        execution['scf_seed'] = dict(seed, source_record='2100.3',
                                     natural_destination_record='2140.2',
                                     transport='same_basis_coefficients')
    if recovery:
        execution['recovered_scf'] = recovery
    record = directory / 'execution.json'
    atomic_json(record, execution)
    started = time.monotonic()
    try:
        env, runtime = qc_runtime_environment('molpro', 'molpro', work_directory=directory,
                                              scratch_min_mb=job.get('scratch_min_mb'),
                                              avoid_scratch_roots=job.get('io_failed_scratch_roots', []))
        command = ['molpro', '-d', runtime['scratch_directory'], '-W', str(directory), '-t', '1',
                   '-n', str(job['ppn']), '-o', f'{name}.out', native.name]
        launcher_stdout = directory / f'{name}.launcher.stdout'
        launcher_stderr = directory / f'{name}.launcher.stderr'
        execution.update(command=command, runtime=runtime,
                         launcher_stdout=str(launcher_stdout), launcher_stderr=str(launcher_stderr))
        atomic_json(record, execution)
        from kinbot.qc_scratch import workspace, molpro_command
        with workspace(directory, runtime, backend='molpro', inputs=[native.name], label=name) as local:
            command = molpro_command(command, runtime, local)
            execution.update(command=command, runtime=runtime)
            atomic_json(record, execution)
            with (local / launcher_stdout.name).open('w') as stdout, \
                    (local / launcher_stderr.name).open('w') as stderr:
                result = subprocess.run(command, cwd=local, env=env, check=False,
                                        stdout=stdout, stderr=stderr)
        execution['returncode'] = result.returncode
        if job.get('reference_review_contract_sha256'):
            execution['reference_review_contract_sha256'] = job['reference_review_contract_sha256']
        identity = [directory/p for p in re.findall(
            r'(?im)^write,\w+,([^,\n]+\.identity_\w+\.out),new,scientific\s*$', native_text)]
        if identity:
            execution['reference_identity_matrices'] = {p.name: file_sha256(p) for p in identity if p.is_file()}
            if all(p.is_file() for p in identity):
                from kinbot.molpro_reference_identity import validate_density_saves
                execution['reference_density_selection'] = validate_density_saves(
                    native_text, (directory/(name+'.out')).read_text(errors='replace'))
        from kinbot.molpro_diagnostics import write_native_report
        report = write_native_report(directory / (name + '.out'), native)
        if report is not None:
            execution['orbital_audit'] = {'path': str(report), 'sha256': file_sha256(report)}
        if result.returncode:
            output = directory / (name + '.out')
            detail = ('\n' + '\n'.join(output.read_text(errors='replace').splitlines()[-20:])
                      if output.is_file() else f' No native output found at {output}.')
            for label, path in (('Launcher stderr', launcher_stderr), ('Launcher stdout', launcher_stdout)):
                tail = '\n'.join(path.read_text(errors='replace').splitlines()[-20:])
                if tail.strip():
                    detail += f'\n{label}:\n{tail}'
            execution['retryable_io'] = molpro_io_failure(detail)
            raise RuntimeError(f'{name}: Molpro exited with status {result.returncode}.{detail}')
        if any(not p.is_file() for p in identity):
            raise RuntimeError('Native calculation did not export every requested reference-identity matrix.')
        execution['status'] = 'complete'
    except BaseException as error:
        execution.update(status='failed', error=str(error))
        if isinstance(error, MolproScratchError):
            execution['retryable_io'] = error.retryable_io
        raise
    finally:
        execution['elapsed_seconds'] = time.monotonic() - started
        atomic_json(record, execution)
        cleanup_qc_runtime(runtime)


def audit(manifest_path):
    with _manifest_lock(Path(manifest_path).expanduser().resolve().parent):
        return _audit(manifest_path)


def _native_occupations(text, root_states=()):
    """Read a natural table for the declared state average or single state.

    A state-specific table cannot stand in for a multi-state average. Stop
    at its own dump boundary rather than reading later numeric output rows.
    """
    headers = list(re.finditer(
        r'(?im)^[ \t]*NATURAL ORBITALS(?:[ \t]*\(([^\n)]*)\))?[ \t]*$', text))
    if not headers:
        return []
    header = headers[-1]
    label = (header[1] or '').upper().strip()
    if label != 'STATE AVERAGED':
        if len(root_states) != 1:
            return []
        state = root_states[0]
        expected = f'STATE {state["root"]}.{state.get("symmetry", 1)}'
        if label not in ('', expected):
            return []
    table = re.split(r'(?im)^\s*(?:Natural orbital dump|PROGRAM\s*\*|'
                     r'(?:NATURAL|CANONICAL|DIABATIC) ORBITALS)',
                     text[header.end():], maxsplit=1)[0]
    return [float(value.replace('D', 'E').replace('d', 'e')) for value in re.findall(
        r'(?m)^\s*\d+\.\d+\s+([-+0-9.]+(?:[EDed][-+]?\d+)?)\s+', table)]


def _xml_evidence(directory, job):
    record = directory / 'execution.json'
    if not record.is_file():
        return None
    execution = json.loads(record.read_text())
    native_name = execution.get('native_input')
    if not isinstance(native_name, str) or not native_name:
        return None
    native = Path(native_name)
    exports = execution.get('orbital_exports', {})
    if (execution.get('status') != 'complete' or not native.is_file()
            or file_sha256(native) != execution.get('native_input_sha256')
            or execution.get('input_sha256') != job['input_sha256']
            or 'scf' not in exports
            or job.get('natural_orbital_probe') and 'natural' not in exports):
        return None
    evidence = {}
    for kind, name in exports.items():
        path = Path(name)
        if path.parent.resolve() != directory.resolve() or not path.is_file() or not path.stat().st_size:
            return None
        try:
            ET.parse(path)
        except ET.ParseError:
            return None
        evidence[kind] = {'path': str(path), 'sha256': file_sha256(path)}
    return evidence


def _checkpoint_evidence(directory, job, execution):
    """Locate the exact native FILE,2 result, including old truncated names.

    Never guess from a glob: a leftover staged predecessor is not evidence of
    this calculation's final orbitals. Legacy recovery requires a completed,
    hash-bound executed input and an unambiguous 31/32-character FILE basename.
    Pre-XML workers executed the prepared input itself; later workers recorded
    a separate transformed native input. Never downgrade partial modern records.
    """
    default = directory / (job['name'] + '.wfu')
    overrides = {key: job.get(key) for key in
                 ('retry_max_reference_vectors', 'retry_max_reference_iterations', 'retry_solver_options')}
    numerical_retry = any(value is not None for value in overrides.values())
    if execution is None:
        if numerical_retry:
            return default, 'Numerical retry requires an execution record.'
        return default, None
    if any(execution.get(key) != value for key, value in overrides.items()):
        return default, 'Execution reference-solver override differs from the requested retry.'
    if execution.get('input_sha256') != job['input_sha256']:
        return default, 'Execution prepared-input SHA-256 is missing or differs from the manifest.'
    modern_fields = {'native_input', 'native_input_sha256', 'input_transform', 'orbital_exports'}
    if modern_fields.intersection(execution):
        name = execution.get('native_input')
        if not isinstance(name, str) or not name:
            return default, 'Transformed execution record is missing its native_input path.'
        native = Path(name)
        if not native.is_file():
            return default, f'Recorded native input is missing: {native}.'
        if not execution.get('native_input_sha256'):
            return default, 'Transformed execution record is missing its native_input_sha256.'
        if file_sha256(native) != execution['native_input_sha256']:
            return default, f'Native input SHA-256 mismatch: {native}.'
        if numerical_retry:
            try:
                if _reference_solver_override(native.read_text(), job.get('retry_max_reference_vectors'),
                                              job.get('retry_max_reference_iterations')) != native.read_text():
                    return default, 'Native input does not contain the requested reference-solver override.'
                if apply_solver_options(native.read_text(), job.get('retry_solver_options', {})) != native.read_text():
                    return default, 'Native input does not contain the requested solver tolerances.'
            except ValueError as error:
                return default, str(error)
    else:
        # Workers before the XML-export change did not create .native.inp.
        # Their command and input_sha256 bind execution to the original .inp.
        if numerical_retry:
            return default, 'Numerical retry requires a hash-bound transformed native input.'
        native = directory / (job['name'] + '.inp')
        command = execution.get('command')
        if (execution.get('status') != 'complete' or execution.get('returncode') != 0
                or not isinstance(command, list) or not command
                or command[-1] not in (native.name, str(native.resolve()))):
            return default, 'Legacy execution does not confirm successful use of the prepared input.'
        if not native.is_file() or file_sha256(native) != job['input_sha256']:
            return default, f'Legacy executed input is missing or its SHA-256 changed: {native}.'
    names = re.findall(r'(?im)^\s*file\s*,\s*2\s*,\s*([A-Za-z0-9_.-]+)\s*(?:,\s*new\s*)?$',
                       native.read_text())
    if len(names) != 1:
        return default, 'Native input does not identify one safe FILE,2 checkpoint.'
    requested = names[0].lower()
    # Old builds truncate at 32; the manual documents a 31-character limit.
    if len(requested) > 31:
        candidates = [directory / name for name in dict.fromkeys(
            [requested, requested[:31], requested[:32]])]
        present = [path for path in candidates if path.is_file() and path.stat().st_size]
        if job.get('initial_wfu_sha256'):
            present = [path for path in present if path.name != requested or
                       file_sha256(path) != job['initial_wfu_sha256']]
        if len(present) > 1:
            return default, 'Ambiguous legacy FILE,2 checkpoints; preserve and reconcile all copies.'
        if present:
            return present[0], None
        return default, 'Missing final native checkpoint; a staged predecessor is not a completed FILE,2 result.'
    return directory / requested, None


def convergence_diagnostics(output_text, job):
    """Check solver completion independently of Molpro's zero return code.

    A failed RHF starting guess does not invalidate subsequently converged
    CASSCF energies. It does invalidate an RHF-only probe or a frozen CASCI
    comparison advertised as starting from converged RHF orbitals.
    """
    errors, warnings = [], []
    casscf = job.get('probe_method') == 'casscf'
    blocks = list(re.finditer(r'PROGRAM\s*\*\s*MULTI\b', output_text, re.I))
    before_multi = output_text[:blocks[0].start()] if blocks else output_text
    multi = output_text[blocks[-1].end():] if blocks else output_text
    multi = re.split(r'PROGRAM\s*\*', multi, maxsplit=1, flags=re.I)[0]
    failed = re.search(r'NO CONVERGENCE\s*!|No convergence in (?:MULTI|MCSCF)|'
                      r'MCSCF[^\n]*not converged', multi, re.I)
    converged = bool(re.search(r'(?im)^\s*CONVERGENCE REACHED!', multi))
    if failed or (casscf and (not converged or re.search(
            r'No orbitals (?:are )?optimized|No orbital optimization', multi, re.I))):
        errors.append('Requested CASSCF orbital optimization did not converge.' if casscf
                      else 'Native MULTI calculation did not converge.')
        thresholds = re.search(r'(?im)^\s*Convergence thresholds[^\r\n]*', multi)
        if thresholds:
            errors[-1] += ' Native ' + thresholds[0].strip() + '.'
    if ((job.get('probe_active_space') or {}).get('orbital_partition')
            or (job.get('probe_active_space') or {}).get('orbital_source') == 'scf'):
        space = job['probe_active_space']
        orbitals = re.search(r'Number of active\s+orbitals:\s+(\d+)', multi)
        electrons = re.findall(r'Number of active electrons:\s+(\d+)', multi)
        if (not orbitals or int(orbitals[1]) != space['active_orbitals'] or not electrons
                or any(int(value) != space['active_electrons'] for value in electrons)):
            errors.append('Native CASSCF orbital partition does not match the requested candidate.')
    if job.get('correlation_checks') and blocks:
        following = output_text[blocks[-1].end():]
        if re.search(r'NO CONVERGENCE\s*!|No convergence in (?:REFERENCE\s+CI|CI|MRCI|RS2|CASPT2)\b|'
                     r'(?:MCSCF|CI|MRCI|RS2|CASPT2)[^\n]*not converged', following, re.I):
            errors.append('Correlated validation contains an unconverged native calculation.')
        if any(check['method'].lower() == 'uccsd(t)-f12b' for check in job['correlation_checks']):
            from kinbot.molpro import validate_molpro_convergence
            from kinbot.anl.results import validate_molpro_triples
            try:
                validate_molpro_convergence(following)
                validate_molpro_triples(following)
            except ValueError as error:
                errors.append('Single-reference validation: ' + str(error))
    # RHF,SO reports the generic MULTI-style iteration-limit message in the
    # HF section, without the conventional "No convergence in rhfpr" line.
    if (re.search(r'No convergence in (?:rhfpr|[ru]?hf)|[RU]?HF[^\n]*not converged', output_text, re.I)
            or re.search(r'NO CONVERGENCE\s*!', before_multi, re.I)):
        message = 'RHF starting calculation did not converge.'
        if casscf and converged and not failed:
            warnings.append(message + ' Subsequent CASSCF optimization converged; RHF energy is not validated.')
        elif not job.get('natural_orbital_probe') or job.get('avas_reference') == 'rhf':
            errors.append(message)
        else:
            warnings.append(message)
    return {'errors': errors, 'warnings': warnings}


def _launch_diagnostics(directory, job):
    """Surface batch/guard failures before diagnosing absent chemistry.

    Read only this allocation's latest restart journal. Older archived
    attempts must not be presented as the reason a new execution failed.
    """
    paths = [directory / (job['name'] + '.stderr')]
    messages = []
    identifier = str(job.get('job_id', ''))
    if re.fullmatch(r'[0-9]+', identifier):
        journal = directory / '.kinbot_slurm_io_retries'
        for allocation in sorted(journal.glob('*/' + identifier)):
            folders = [allocation, *[p for p in allocation.glob('restart-*')
                                    if p.is_dir() and re.fullmatch(r'restart-[0-9]+', p.name)]]
            latest = max(folders, key=lambda p: int(p.name.split('-')[1]) if p != allocation else 0)
            record = latest / 'execution.json'
            if record.is_file():
                try:
                    detail = json.loads(record.read_text())
                    attempts = detail.get('attempts', [])
                    code = attempts[-1].get('returncode') if attempts else None
                    messages.append(f'Slurm retry guard status={detail.get("status")!r}, '
                                    f'returncode={code!r}; journal: {record}.')
                except (OSError, ValueError):
                    messages.append(f'Unreadable Slurm retry journal: {record}.')
            attempts = [p for p in latest.glob('attempt-*')
                        if p.is_dir() and re.fullmatch(r'attempt-[0-9]+', p.name)]
            if attempts:
                paths.append(max(attempts, key=lambda p: int(p.name.split('-')[1])) / 'stderr.txt')
    seen = set()
    for path in paths:
        if not path.is_file():
            continue
        with path.open('rb') as stream:
            stream.seek(0, 2)
            stream.seek(max(0, stream.tell() - 2000))
            tail = stream.read().decode(errors='replace').strip()
        if tail and tail not in seen:
            messages.append(f'Batch diagnostic ({path}):\n{tail}')
            seen.add(tail)
    return messages


def _audit(manifest_path, *, live_job_ids=()):
    path = Path(manifest_path).expanduser().resolve()
    manifest = json.loads(path.read_text())
    complete = True
    for job in manifest['jobs']:
        for key in ('energy_hartree', 'root_energies_hartree', 'state_energies_hartree',
                    'natural_occupations', 'natural_molden_sha256',
                    'orbital_error', 'orbital_xml', 'occupation_source', 'validation_errors',
                    'validation_warnings', 'terminal_failure',
                    'correlated_energies_hartree', 'correlated_energy_source'):
            job.pop(key, None)
        directory = Path(job['directory'])
        execution_path = directory / 'execution.json'
        # The allocation owns every attempt until the shared retry guard
        # exits. Do not parse/hash a WFU that it can still replace, or allow
        # legacy Molden files to override an explicitly unfinished attempt.
        live = str(job.get('job_id')) in live_job_ids
        execution = (json.loads(execution_path.read_text())
                     if not live and execution_path.is_file() else None)
        if live or execution is not None and execution.get('status') != 'complete':
            job['status'] = 'incomplete'
            job['validation_errors'] = [
                'Slurm allocation is still live.' if live else
                f'Execution status is {execution.get("status")!r}; a complete record is required.']
            if not live:
                job['validation_errors'].extend(_launch_diagnostics(directory, job))
            for key in ('output_sha256', 'wfu', 'wfu_sha256', 'molden', 'molden_sha256'):
                job.pop(key, None)
            continue
        stem = job['name']
        output = directory / f'{stem}.out'
        wfu, checkpoint_error = _checkpoint_evidence(directory, job, execution)
        molden = directory / f'{stem}.molden'
        natural = directory / f'{stem}.natural.molden'
        inp = directory / f'{stem}.inp'
        output_text = output.read_text(errors='replace') if output.is_file() else ''
        if not output_text.strip():
            job['status'] = 'incomplete'
            job['validation_errors'] = [f'Native output is missing or empty: {output}.',
                                        *_launch_diagnostics(directory, job)]
            for key in ('output_sha256', 'wfu', 'wfu_sha256', 'molden', 'molden_sha256'):
                job.pop(key, None)
            continue
        xml_evidence = _xml_evidence(directory, job)
        errors = []
        if not inp.is_file() or file_sha256(inp) != job['input_sha256']:
            errors.append('Prepared input is missing or its SHA-256 changed.')
        if 'Molpro calculation terminated' not in output_text:
            errors.append('Native output lacks normal termination.')
        if not re.search(r'SETTING\s+RECONDONE\s*=\s*1(?:\.0*)?\s*(?:\n|$)', output_text):
            errors.append('Native output lacks RECONDONE=1.')
        if 'ERROR EXIT' in output_text:
            errors.append('Native output contains ERROR EXIT.')
        convergence = convergence_diagnostics(output_text, job)
        errors.extend(convergence['errors'])
        job['validation_warnings'] = convergence['warnings']
        if checkpoint_error:
            errors.append(checkpoint_error)
        if not wfu.is_file() or not wfu.stat().st_size:
            errors.append(f'Missing or empty FILE,2 checkpoint: {wfu}.')
        if not (xml_evidence or molden.is_file() and molden.stat().st_size):
            errors.append('No valid hash-bound XML exports or nonempty legacy Molden file.')
        normal = (inp.is_file()
                  and file_sha256(inp) == job['input_sha256']
                  and 'Molpro calculation terminated' in output_text
                  and re.search(r'SETTING\s+RECONDONE\s*=\s*1(?:\.0*)?\s*(?:\n|$)', output_text)
                  and 'ERROR EXIT' not in output_text and not convergence['errors'])
        job['status'] = ('complete' if normal and not checkpoint_error and wfu.is_file()
                         and wfu.stat().st_size and (xml_evidence or
                         molden.is_file() and molden.stat().st_size) else 'incomplete')
        energies = re.findall(r'SETTING\s+RECONENERGY\s*=\s*'
                              r'([-+0-9.]+(?:[EDed][-+]?\d+)?)\s+(?:AU|HARTREE)\b',
                              output_text, re.IGNORECASE)
        if energies:
            value = float(energies[-1].replace('D', 'E').replace('d', 'e'))
            if np.isfinite(value):
                job['energy_hartree'] = value
        if job.get('channel') and 'energy_hartree' not in job:
            job['status'] = 'incomplete'
            errors.append('Missing finite RECONENERGY for the distance profile.')
        if job.get('natural_orbital_probe'):
            roots = {int(index): float(value.replace('D', 'E').replace('d', 'e'))
                     for index, value in re.findall(r'SETTING\s+RECONR(\d+)\s*=\s*'
                         r'([-+0-9.]+(?:[EDed][-+]?\d+)?)\s+(?:AU|HARTREE)\b',
                         output_text, re.IGNORECASE)}
            expected = (len(job.get('root_states', [])) or
                        sum((job.get('spin_context') or {}).get('roots', {}).values()))
            if (not roots or any(not np.isfinite(value) for value in roots.values())
                    or (expected and sorted(roots) != list(range(1, expected + 1)))):
                job['status'] = 'incomplete'
                errors.append(f'Missing/nonfinite RECONR roots; expected {expected}, found {sorted(roots)}.')
            else:
                job['root_energies_hartree'] = [roots[index] for index in sorted(roots)]
                if job.get('root_states'):
                    if len(job['root_states']) != len(roots):
                        job['status'] = 'incomplete'
                        errors.append('Root-state labels do not match the number of energies.')
                    else:
                        job['state_energies_hartree'] = [
                            {**state, 'energy_hartree': value}
                            for state, value in zip(job['root_states'], job['root_energies_hartree'])]
            occupations = ([] if xml_evidence or not natural.is_file() else
                           [float(value.replace('D', 'E').replace('d', 'e'))
                            for value in re.findall(r'(?im)^\s*Occup\s*=\s*'
                                r'([-+0-9.]+(?:[EDed][-+]?\d+)?)', natural.read_text())])
            if not occupations and xml_evidence:
                occupations = _native_occupations(output_text, job.get('root_states', []))
            if (not occupations or any(not np.isfinite(value) or not 0. <= value <= 2.
                                       for value in occupations)
                    or abs(sum(occupations) - manifest['electron_count']) > .02):
                job['status'] = 'incomplete'
                errors.append('Invalid natural occupations: '
                              f'count={len(occupations)}, sum={sum(occupations):.8g}, '
                              f'expected electrons={manifest["electron_count"]}.')
            else:
                job['natural_occupations'] = occupations
                if not xml_evidence and natural.is_file():
                    job['natural_molden_sha256'] = file_sha256(natural)
                job['occupation_source'] = ('native_natural_table' if xml_evidence else 'molden')
            checks = job.get('correlation_checks', [])
            if checks:
                values = {int(index): float(value.replace('D', 'E').replace('d', 'e'))
                    for index, value in re.findall(r'SETTING\s+RECONC(\d+)\s*=\s*'
                        r'([-+0-9.]+(?:[EDed][-+]?\d+)?)\s+(?:AU|HARTREE)\b', output_text, re.I)}
                labels = _correlation_labels(checks, job['root_states'], job.get('state_model'))
                if (not values and execution and execution.get('returncode') == 0
                        and execution.get('native_input')):
                    try:
                        values = _recover_legacy_caspt2_energies(output_text,
                            Path(execution['native_input']).read_text(), labels, manifest['electron_count'])
                    except (ValueError, OSError) as error:
                        errors.append('Native energy recovery unavailable: ' + str(error))
                    else:
                        job['correlated_energy_source'] = {
                            'kind': 'native_rspt2_state_lines_legacy_selector_recovery',
                            'output_sha256': file_sha256(output),
                            'native_input_sha256': file_sha256(execution['native_input'])}
                        job['validation_warnings'].append(
                            'Recovered printed CASPT2-F12 state energies after the legacy ENERGY(2) bug; '
                            'the historical F12b spelling is not a distinct Molpro MR approximation.')
                if (sorted(values) != list(range(1, len(labels) + 1))
                        or not all(np.isfinite(value) for value in values.values())):
                    job['status'] = 'incomplete'
                    errors.append('Missing/nonfinite correlated validation energies.')
                else:
                    from kinbot.molpro import vrc_multireference_energy_contract
                    try:
                        _validate_cc_check_energies(output_text, labels, values)
                    except ValueError as error:
                        job['status'] = 'incomplete'
                        errors.append(str(error))
                    job['correlated_energies_hartree'] = [
                        {**label, 'energy_hartree': values[index],
                         'native_method': (vrc_multireference_energy_contract(label['method']) or
                                           {'method': 'UCCSD(T)-F12b'})['method']}
                        for index, label in enumerate(labels, 1)]
                if ('UNDEFINED, ASSUMING' in output_text and
                        job.get('correlated_energy_source', {}).get('kind') !=
                        'native_rspt2_state_lines_legacy_selector_recovery'):
                    job['status'] = 'incomplete'
                    errors.append('Native output used an undefined energy variable.')
        complete &= job['status'] == 'complete'
        if job['status'] == 'complete':
            job.pop('missing_since', None)
            job['output_sha256'] = file_sha256(output)
            job['wfu'] = str(wfu.resolve())
            job['wfu_sha256'] = file_sha256(wfu)
            if molden.is_file() and molden.stat().st_size:
                job['molden'] = str(molden.resolve())
                job['molden_sha256'] = file_sha256(molden)
            if xml_evidence:
                job['orbital_xml'] = xml_evidence
        else:
            job['validation_errors'] = errors
            for key in ('output_sha256', 'wfu', 'wfu_sha256', 'molden',
                        'molden_sha256'):
                job.pop(key, None)
    indexed = {job['name']: job for job in manifest['jobs']}
    for job in manifest['jobs']:
        predecessor = job.get('orbital_predecessor')
        if job['status'] == 'complete' and predecessor:
            prior = indexed[predecessor]
            if (prior['status'] != 'complete'
                    or job.get('initial_wfu_sha256') != prior.get('wfu_sha256')):
                job['status'] = 'incomplete'
                job['orbital_error'] = 'Predecessor WFU no longer matches the staged orbital source.'
                job['validation_errors'] = [job['orbital_error']]
                for key in ('output_sha256', 'wfu', 'wfu_sha256', 'molden', 'molden_sha256'):
                    job.pop(key, None)
    complete = all(job['status'] == 'complete' for job in manifest['jobs'])
    manifest['coverage'] = _coverage(manifest)
    manifest['orientation_profiles'] = [{
        'job': job['name'], 'basis': job['basis'],
        'distance_angstrom': job['distance_angstrom'], **job['orientation'],
        'status': job['status'], 'states': job.get('state_energies_hartree', []),
        'natural_occupations': job.get('natural_occupations', []),
    } for job in manifest['jobs'] if job.get('orientation')]
    sufficient = all(item['status'] == 'covered' for item in manifest['coverage'].values())
    manifest['status'] = ('complete_review_required' if complete and sufficient else
                          'coverage_refinement_required' if complete else 'incomplete')
    manifest['review'] = ('Inspect orbital characters and occupations in every '
                          'native XML/WFU or legacy Molden/output pair; no active orbital identity is '
                          'selected automatically.')
    atomic_json(path, manifest)
    return manifest


def _coverage(manifest):
    """Check asymptotic plateaus and flag intervals needing denser probes."""
    if manifest.get('orbital_character_only'):
        # This panel identifies candidate orbitals. Its SCF energies are not
        # a reviewed interaction potential and must not drive radial refinement.
        return {}
    policy = manifest.get('scan_policy')
    if not policy:
        return {}
    groups = {}
    for job in manifest['jobs']:
        if job.get('channel'):
            groups.setdefault(f'{job["basis"]}:{job["channel"]}', []).append(job)
    result = {}
    for name, jobs in groups.items():
        jobs.sort(key=lambda job: job['distance_angstrom'])
        missing = any(job['status'] != 'complete' or 'energy_hartree' not in job for job in jobs)
        refinements = []
        if not missing:
            for left, right in zip(jobs, jobs[1:]):
                energy_jump = abs(right['energy_hartree'] - left['energy_hartree'])
                if left.get('root_energies_hartree') and right.get('root_energies_hartree'):
                    if len(left['root_energies_hartree']) != len(right['root_energies_hartree']):
                        raise ValueError('Reconnaissance root count changed within one channel.')
                    energy_jump = max(energy_jump, max(abs(first - second) for first, second in
                        zip(left['root_energies_hartree'], right['root_energies_hartree'])))
                if energy_jump > policy['maximum_energy_step_hartree']:
                    refinements.append((left['distance_angstrom'] + right['distance_angstrom']) / 2.)
                if left.get('natural_occupations') and right.get('natural_occupations'):
                    # Sorted populations detect a changing subspace dimension;
                    # orbital identity itself still needs overlap/character review.
                    lo = sorted(left['natural_occupations'])
                    hi = sorted(right['natural_occupations'])
                    if len(lo) != len(hi) or np.max(np.abs(np.array(lo) - hi)) > policy['maximum_occupation_change']:
                        refinements.append((left['distance_angstrom'] + right['distance_angstrom']) / 2.)
        def plateau_energies(job):
            return [job['energy_hartree'], *job.get('root_energies_hartree', [])]

        # Every sampled electronic component must have reached separation.
        # A flat ground-state curve alone can hide an unsettled excited root
        # and therefore give a false green light to the outer-state model.
        plateau = (not missing and len(jobs) >= 3
                   and jobs[-1]['distance_angstrom'] >= policy['asymptote_angstrom']
                   and all(len(plateau_energies(left)) == len(plateau_energies(right))
                           and all(abs(second - first) <= policy['asymptote_tolerance_hartree']
                                   for first, second in zip(plateau_energies(left), plateau_energies(right)))
                           for left, right in zip(jobs[-3:], jobs[-2:])))
        result[name] = {
            'status': 'covered' if plateau and not refinements else 'refine',
            'asymptote_plateau': plateau,
            'state_energy_profiles': [{
                'job': job['name'], 'distance_angstrom': job['distance_angstrom'],
                'states': job.get('state_energies_hartree', []),
            } for job in jobs if job.get('state_energies_hartree')],
            'suggested_midpoints_angstrom': sorted(set(refinements)),
            'suggested_far_distance_angstrom': (None if plateau else
                jobs[-1]['distance_angstrom'] * policy['far_step_factor']),
            'occupation_candidate_orbitals': max((sum(policy['occupation_lower_bound'] < value < policy['occupation_upper_bound']
                for value in job.get('natural_occupations', [])) for job in jobs), default=0),
            'occupation_candidates_by_geometry': [{
                'job': job['name'], 'distance_angstrom': job['distance_angstrom'],
                'orbitals_one_based': [index for index, value in enumerate(
                    job.get('natural_occupations', []), 1)
                    if policy['occupation_lower_bound'] < value < policy['occupation_upper_bound']],
                'electron_population': sum(value for value in job.get('natural_occupations', [])
                    if policy['occupation_lower_bound'] < value < policy['occupation_upper_bound']),
            } for job in jobs if job.get('natural_occupations')],
            'limitation': 'Energy/occupation coverage is a diagnostic; it does '
                          'not certify orbital overlaps, root identity, or a minimal CAS.',
        }
    return result


def refine(manifest_path):
    """Add only coverage-requested midpoints and longer-range probes."""
    with _manifest_lock(Path(manifest_path).expanduser().resolve().parent):
        return _refine(manifest_path)


def _refine(manifest_path):
    path = Path(manifest_path).expanduser().resolve()
    manifest = _audit(path)
    if manifest['status'] != 'coverage_refinement_required':
        raise ValueError('Refinement requires completed probes with unresolved coverage.')
    data = _load(manifest['source'])
    require_explicit_probe_space(data)
    channels = {}
    for anchor in data['anchors']:
        if anchor.get('_channel'):
            original = next((item for item in data['coordinates']
                             if isinstance(item, dict) and item['name'] == anchor['_channel']), {})
            entry = channels.setdefault(anchor['_channel'], {
                **deepcopy(original), 'name': anchor['_channel'], 'bond': anchor['bond'],
                'distances_angstrom': []})
            entry['distances_angstrom'].append(anchor['distance_angstrom'])
    for name, coverage in manifest['coverage'].items():
        channel = name.rsplit(':', 1)[1]
        values = channels[channel]['distances_angstrom']
        values.extend(coverage['suggested_midpoints_angstrom'])
        far = coverage['suggested_far_distance_angstrom']
        if far is not None:
            values.append(far)
    source = json.loads(Path(manifest['source']).read_text())
    source.pop('anchors', None)
    source['coordinates'] = list(channels.values())
    for entry in source['coordinates']:
        entry['distances_angstrom'] = sorted(set(entry['distances_angstrom']))
    if source.get('kinbot_input'):
        source['kinbot_input'] = str((Path(manifest['source']).parent /
                                     source['kinbot_input']).resolve())
    request = path.parent / 'refinement_request.json'
    request.write_text(json.dumps(source, indent=2) + '\n')
    return _prepare(request, path.parent)


def submit(manifest_path, *, max_jobs=8):
    """Submit new probes once, under a manifest-owned Slurm job cap."""
    with _manifest_lock(Path(manifest_path).expanduser().resolve().parent):
        return _submit(manifest_path, max_jobs=max_jobs)


def retry_failed(manifest_path, *, include_interrupted=False, job_ids=None, max_reference_vectors=None,
                 max_reference_iterations=None, solver_options=None):
    """Archive selected stopped attempts while preserving inputs and live peers."""
    path = Path(manifest_path).expanduser().resolve()
    with _manifest_lock(path.parent):
        return _retry_failed(path, include_interrupted=include_interrupted, job_ids=job_ids,
                             max_reference_vectors=max_reference_vectors,
                             max_reference_iterations=max_reference_iterations, solver_options=solver_options)


def _retry_failed(path, *, include_interrupted=False, job_ids=None, live_job_ids=None,
                  max_reference_vectors=None, max_reference_iterations=None, solver_options=None,
                  recheck_scf_seed=False):
    import shutil
    solver_options = validate_solver_options(solver_options or {})
    controls = {key: value for key, value in (
        ('max_reference_vectors', max_reference_vectors),
        ('max_reference_iterations', max_reference_iterations)) if value is not None}
    if controls or solver_options:
        validate_mr_correlation_options(controls)
        if not job_ids:
            raise ValueError('A numerical retry override requires explicit --job-id selection.')
    if live_job_ids is None:
        queued = subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                                text=True, check=True, timeout=30)
        live_job_ids = queued.stdout.split()
    live = set(live_job_ids)
    manifest = _audit(path, live_job_ids=live)
    selected = [j for j in manifest['jobs'] if job_ids is None or str(j.get('job_id')) in job_ids]
    rechecks = {}
    if recheck_scf_seed:
        from kinbot.recon_scf_seed import projection_failure
        if not job_ids or not selected or controls or solver_options or include_interrupted:
            raise ValueError('Seed recheck requires explicit completed job IDs without other retry changes.')
        for job in selected:
            proof = projection_failure(job, recheck_completed=True)
            if not proof:
                raise ValueError('Seed recheck needs an intact completed legacy forced-OVL SCF probe.')
            if any(j.get('orbital_predecessor') == job['name'] for j in manifest['jobs']):
                raise ValueError('Seed recheck cannot invalidate an orbital dependency chain.')
            rechecks[job['name']] = proof
    if any(str(job.get('job_id')) in live for job in selected):
        raise RuntimeError('Wait for all manifest allocations to leave before retrying.')
    if any(job.get('submission_status') == 'submitting' for job in manifest['jobs']):
        raise RuntimeError('Reconcile uncertain Slurm submissions before retrying.')
    retried = []
    for job in selected:
        recheck = rechecks.get(job['name'])
        if job['status'] == 'complete' and not recheck:
            continue
        directory = Path(job['directory'])
        record = directory / 'execution.json'
        execution = json.loads(record.read_text()) if record.is_file() else {}
        interrupted = (include_interrupted and
                       execution.get('status') in (None, 'running'))
        invalid_final = bool(recheck or job_ids and solver_options and execution.get('status') == 'complete')
        if (not job.get('job_id') or
                execution.get('status') != 'failed' and not interrupted and not invalid_final):
            continue
        scheduler = slurm_job_record(str(job['job_id']))
        if scheduler.get('JobId') == str(job['job_id']):
            if (scheduler.get('WorkDir') != str(directory)
                    or scheduler.get('JobState') not in SLURM_TERMINAL_STATES):
                raise RuntimeError('Slurm still owns the retry target, or its identity differs.')
            job['scheduler_record'] = scheduler
        inp = directory / (job['name'] + '.inp')
        if file_sha256(inp) != job['input_sha256']:
            raise ValueError('Failed probe input changed; do not retry under the old manifest.')
        if controls:
            _reference_solver_override(inp.read_text(), max_reference_vectors, max_reference_iterations)
        if solver_options:
            apply_solver_options(inp.read_text(), solver_options)
        attempt = int(job.get('attempt', 1))
        archive = directory / 'attempts' / str(attempt)
        # A failed visualization export does not invalidate the preceding
        # converged RHF record. Reuse only an explicitly evidenced record,
        # bound to this exact input; never reuse a failed SCF calculation.
        output = directory / (job['name'] + '.out')
        output_text = output.read_text(errors='replace') if output.is_file() else ''
        scf = re.search(r'(?i)\{rhf;[^{}]*orbital,(\d+\.\d+)\}', inp.read_text())
        energies = re.findall(r'(?im)^\s*!RHF\s+STATE\s+\d+\.\d+\s+Energy\s+'
                              r'([-+0-9.]+(?:[EDed][-+]?\d+)?)', output_text)
        wfu = directory / (job['name'] + '.wfu')
        output_lower = output_text.lower()
        xml_title_overflow = all(marker in output_lower for marker in (
            'dump information in style xml', 'string overflow inserting',
            'strings::string_insert'))
        export_only = (('basis set contains h or higher' in output_lower
                        or xml_title_overflow)
                       and execution.get('input_sha256') == job['input_sha256']
                       and scf and energies and wfu.is_file() and wfu.stat().st_size
                       and re.search(r'PROGRAMS\s+\*[^\n]*RHF-SCF', output_text))
        if energies:
            scf_energy = float(energies[-1].replace('D', 'E').replace('d', 'e'))
            export_only = export_only and np.isfinite(scf_energy)
        recovery = None
        if export_only:
            recovery = {'path': str(archive / wfu.name),
                        'sha256': file_sha256(wfu),
                        'record': scf[1], 'input_sha256': job['input_sha256'],
                        'energy_hartree': scf_energy}
        archive.mkdir(parents=True, exist_ok=False)
        keep = {job['name'] + '.inp', job['name'] + '.sbatch'}
        for artifact in directory.iterdir():
            if artifact.is_file():
                if artifact.name in keep:
                    shutil.copy2(artifact, archive / artifact.name)
                else:
                    artifact.rename(archive / artifact.name)
        job.setdefault('attempt_history', []).append({
            'attempt': attempt, 'job_id': job['job_id'],
            'archive': str(archive),
            'scheduler_record': job.get('scheduler_record'),
            'error': execution.get('error') or ('Explicitly retried interrupted allocation'
                                               if interrupted else None),
            'input_sha256': job['input_sha256'],
            'input_transform': execution.get('input_transform'),
            **({'recovery': recheck} if recheck else {}),
            **{key: job.get(key) for key in
               ('retry_max_reference_vectors', 'retry_max_reference_iterations', 'retry_solver_options')}})
        for key in ('job_id', 'missing_since', 'terminal_failure', 'scheduler_record'):
            job.pop(key, None)
        job.update(attempt=attempt + 1, submission_status='retry_ready', status='incomplete')
        for key, value in controls.items():
            job['retry_' + key] = value
        if solver_options:
            job['retry_solver_options'] = deepcopy(solver_options)
            for key in ('max_reference_vectors', 'max_reference_iterations'):
                if 'retry_' + key in job and key in solver_options.get('correlation', {}):
                    job['retry_' + key] = solver_options['correlation'][key]
        if recovery:
            job['recovered_scf'] = recovery
        retried.append(job['name'])
        atomic_json(path, manifest)
    return {'manifest': str(path), 'retry_ready': retried,
            'completed_preserved': sum(job['status'] == 'complete' for job in manifest['jobs'])}


def _retry_seed_projection_failures(path, manifest, live):
    """One migration retry for evidenced forced-OVL failures; keep all peers."""
    from kinbot.recon_scf_seed import projection_failure
    for job in manifest['jobs']:
        if str(job.get('job_id')) in live:
            continue
        proof = projection_failure(job)
        if proof is None:
            continue
        identifier = str(job['job_id'])
        result = _retry_failed(path, job_ids={identifier}, live_job_ids=live)
        if result['retry_ready']:
            manifest = json.loads(path.read_text())
            current = next(j for j in manifest['jobs'] if j['name'] == job['name'])
            current['attempt_history'][-1]['recovery'] = proof
            atomic_json(path, manifest)
            print(f'{job["name"]}: archived job {identifier}; retrying same-basis seed transport once.', flush=True)
    return manifest


def _recover_slurm_startup_failures(path, manifest, live):
    """Recover positively identified pre-shell I/O failures, never QC failures.

    Slurm's ESLURMD_IO_ERROR (4021) can appear as ExitCode=0:53. The
    batch-script retry guard cannot run when slurmstepd cannot open stdio.
    Retain scontrol evidence before it expires; absence from squeue alone is
    never authorization to resubmit. Caller holds the manifest lock.
    """
    for job in manifest['jobs']:
        job.pop('startup_failure', None)
        identifier = str(job.get('job_id', ''))
        if not identifier or identifier in live or job.get('status') == 'complete':
            continue
        directory = Path(job['directory'])
        # Capture termination evidence for started jobs too. Accounting may
        # be disabled and scontrol records can expire before manual review.
        record = slurm_job_record(identifier)
        if record.get('JobId') == identifier:
            job['scheduler_record'] = record
        else:
            record = job.get('scheduler_record', {})
            # A cached RUNNING/PENDING observation cannot keep a vanished
            # allocation alive forever after the scheduler purges it.
            if record.get('JobState') not in SLURM_TERMINAL_STATES:
                record = {}
        if (record.get('JobId') != identifier or
                record.get('WorkDir') != str(directory)):
            continue
        state = record.get('JobState')
        if state and state not in SLURM_TERMINAL_STATES:
            # scontrol can already know about a requeue absent from the
            # preceding squeue snapshot. Keep its allocation reserved.
            live.add(identifier)
            job.pop('missing_since', None)
            continue
        # Runtime evidence rules out this narrow automatic startup recovery,
        # but must not prevent recording an OOM/cancellation/requeue above.
        if any((directory / name).exists() for name in
               ('execution.json', job['name'] + '.native.inp', job['name'] + '.out')):
            continue
        node = record.get('BatchHost', '')
        if (state != 'FAILED' or record.get('ExitCode') != '0:53'
                or not _SAFE.fullmatch(node) or node == 'None'):
            continue
        # Empty stdio files are possible when only one stream opened. A
        # nonempty stream or retry-guard journal means the shell may have run.
        if (any(p.is_file() and p.stat().st_size for suffix in
                ('.stdout', '.stderr', '.launcher.stdout', '.launcher.stderr')
                for p in [directory / (job['name'] + suffix)])
                or any((directory / '.kinbot_slurm_io_retries').glob(
                    '*/' + identifier))):
            continue
        manifest['slurm_excluded_nodes'] = sorted(
            set(manifest.get('slurm_excluded_nodes', [])) | {node})
        job['startup_failure'] = (
            f'Slurm batch I/O initialization failed on {node} '
            f'(job {identifier}, ExitCode=0:53); no native execution was started.')
        # Allow filesystem publication before deciding the batch never ran.
        if time.time() - job.setdefault('missing_since', time.time()) <= 60.:
            continue
        history = job.setdefault('startup_retry_history', [])
        limit = job.get('startup_retry_limit', manifest.get('startup_retries', 2))
        if type(limit) is not int or limit < 0:
            raise ValueError('Startup retry limit must be a nonnegative integer.')
        if len(history) >= limit:
            job['startup_failure'] += f' Automatic startup retry limit ({limit}) exhausted.'
            continue
        inp = directory / (job['name'] + '.inp')
        if file_sha256(inp) != job['input_sha256']:
            raise ValueError('Probe input changed; refusing automatic Slurm startup retry.')
        history.append({'job_id': identifier, 'scheduler_record': deepcopy(record),
                        'input_sha256': job['input_sha256'], 'time': time.time(),
                        'reason': job['startup_failure']})
        for key in ('job_id', 'missing_since', 'terminal_failure', 'startup_failure',
                    'scheduler_record'):
            job.pop(key, None)
        job.update(status='incomplete', submission_status='retry_ready')
        print(f'{job["name"]}: retrying Slurm startup failure from job {identifier}; '
              f'excluding {node} (retry {len(history)}/{limit}).', flush=True)
    # Evidence and retry readiness are committed together, before sbatch.
    atomic_json(Path(path), manifest)


def _slurm_submission_options(batch, excluded, *, allow_inherited=()):
    """Pass exclusions through Slurm's CLI, preserving script/site choices.

    SBATCH_EXCLUDE is not a Slurm input variable. Read it only for compatibility
    with older KinBot instructions and translate it into the real --exclude.
    """
    if any(not isinstance(node, str) or not _SAFE.fullmatch(node) for node in excluded):
        raise ValueError('Invalid recorded Slurm exclusion node.')
    if any(not isinstance(node, str) or not _SAFE.fullmatch(node) for node in allow_inherited):
        raise ValueError('Invalid inherited exclusion override node.')
    values = [os.environ.get('SBATCH_EXCLUDE', ''), *([] if allow_inherited else sorted(excluded))]
    auto_time = False
    partition = ''
    for line in batch.splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            break  # Slurm stops reading directives at the first shell command.
        if line.strip() == '# KINBOT_WALLTIME=partition_max':
            auto_time = True
        if not line.lstrip().startswith('#SBATCH '):
            continue
        tokens = shlex.split(line.lstrip()[len('#SBATCH '):])
        for i, token in enumerate(tokens):
            if token in ('--partition', '-p') and i + 1 < len(tokens):
                partition = tokens[i + 1]
            elif token.startswith('--partition='):
                partition = token.split('=', 1)[1]
            elif token.startswith('-p') and len(token) > 2:
                partition = token[2:]
            if token in ('--exclude', '-x') and i + 1 < len(tokens):
                values.append(tokens[i + 1])
            elif token.startswith('--exclude='):
                values.append(token.split('=', 1)[1])
            elif token.startswith('-x') and len(token) > 2:
                values.append(token[2:])
    if allow_inherited:
        inherited = set()
        for value in dict.fromkeys(v for v in values if v):
            expanded = subprocess.run(['scontrol', 'show', 'hostnames', value],
                capture_output=True, text=True, check=True, timeout=20)
            inherited.update(expanded.stdout.split())
        # The explicit override applies to inherited script/environment values.
        # Newly learned allocation failures in the manifest still win.
        value = ','.join(sorted((inherited - set(allow_inherited)) | set(excluded)))
    else:
        value = ','.join(dict.fromkeys(v for v in values if v))
    options = ['--exclude=' + value] if value or allow_inherited else []
    if auto_time and not os.environ.get('SBATCH_TIMELIMIT'):
        from kinbot.run_walltime import partition_walltime
        options.append('--time=' + partition_walltime(
            os.environ.get('SBATCH_PARTITION', partition)))
    return options


def _retry_allocation_failures(path, manifest, live, limit):
    """Opt-in, bounded recovery after a worker/guard allocation is lost.

    A finalized I/O failure or an unfinished known worker can be retried.
    Cancellation, resource limits, explicit QC errors and uncertain sbatch
    outcomes are never treated as transient I/O. Caller holds manifest lock.
    """
    from kinbot.slurm_retry import _text, io_failure
    for name in [j['name'] for j in manifest['jobs']]:
        job = next(j for j in manifest['jobs'] if j['name'] == name)
        identifier = str(job.get('job_id', ''))
        if not identifier or identifier in live or job.get('status') == 'complete':
            continue
        history = job.get('allocation_retry_history', [])
        if len(history) >= limit:
            continue
        directory = Path(job['directory'])
        execution_path = directory / 'execution.json'
        if not execution_path.is_file():
            continue  # Pre-shell failures have their own evidence-based path.
        execution = json.loads(execution_path.read_text())
        scheduler = job.get('scheduler_record', {})
        if scheduler.get('JobState') in ('CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'DEADLINE'):
            continue
        detail = '\n'.join(_text(directory / (job['name'] + suffix))
                           for suffix in ('.out', '.stderr', '.launcher.stderr', '.launcher.stdout'))
        explicit_io = execution.get('retryable_io') is True or io_failure(detail)
        interrupted = execution.get('status') == 'running'
        # A Slurm requeue can fail before opening its stdio, leaving all native
        # files and the running guard journal from the previous host. Bind this
        # exception to the exact allocation/restart; old QC text is not evidence
        # that the replacement shell started. Never infer it from a missing job.
        restarts = str(scheduler.get('Restarts', ''))
        node = scheduler.get('BatchHost', '')
        requeue_io = bool(interrupted and restarts.isdigit() and int(restarts) > 0
            and scheduler.get('JobId') == identifier
            and scheduler.get('WorkDir') == str(directory)
            and scheduler.get('JobState') == 'FAILED' and scheduler.get('ExitCode') == '0:53'
            and _SAFE.fullmatch(node) and node != 'None'
            and execution.get('host') and node != execution['host']
            and scheduler.get('RunTime') in ('00:00:00', '00:00:01')
            and not any((directory / '.kinbot_slurm_io_retries').glob(
                '*/' + identifier + '/restart-' + restarts)))
        if requeue_io:
            explicit_io = True
            manifest['slurm_excluded_nodes'] = sorted(set(manifest.get('slurm_excluded_nodes', [])) | {node})
        if not (explicit_io and execution.get('status') == 'failed' or interrupted):
            continue
        # An unfinished worker needs a publication grace period. A finalized
        # I/O failure can be retried immediately after scheduler exit; waiting
        # here would let _finished terminate a single-job campaign first.
        if interrupted and time.time() - job.setdefault('missing_since', time.time()) <= 60.:
            continue
        if not explicit_io and re.search(
                r'no convergence|not converged|error exit|\?\s*error|global error', detail, re.I):
            continue
        history = [*history, {'job_id': identifier, 'time': time.time(),
            'reason': ('requeued_startup_io' if requeue_io else
                       'native_io' if explicit_io else 'interrupted_allocation'),
            'scheduler_record': deepcopy(scheduler)}]
        job['allocation_retry_history'] = history
        atomic_json(Path(path), manifest)
        # Uses exactly the same archival/seed restoration as an explicit retry.
        result = _retry_failed(Path(path), include_interrupted=True,
                               job_ids={identifier}, live_job_ids=live)
        manifest = json.loads(Path(path).read_text())
        if job['name'] in result['retry_ready']:
            print(f'{job["name"]}: retrying {history[-1]["reason"]} from job {identifier} '
                  f'(allocation retry {len(history)}/{limit}).', flush=True)
        else:
            # A late valid result appeared during the second audit. It did
            # not consume a retry and must not be announced as resubmitted.
            current = next(j for j in manifest['jobs'] if j['name'] == name)
            current['allocation_retry_history'] = history[:-1]
            atomic_json(Path(path), manifest)
    return manifest


def _submit(manifest_path, *, max_jobs=8, max_new_jobs=None,
            live_job_ids=None, startup_checked=False, audited_manifest=None):
    if isinstance(max_jobs, bool) or not isinstance(max_jobs, int) or max_jobs < 1:
        raise ValueError('max_jobs must be positive.')
    path = Path(manifest_path).expanduser().resolve()
    pending = json.loads(path.read_text())
    from kinbot import run_policy
    if run_policy.find(path.parent) is not None and not pending.get('shared_policy_file'):
        raise ValueError('Reconnaissance inside a resource pool requires explicit shared-pool staging.')
    # An unknown allocation consumes unknown capacity. Check the whole
    # manifest before the startup/dependency/cap filters can skip its entry.
    if any(job.get('submission_status') == 'submitting' for job in pending['jobs']):
        raise RuntimeError('Interrupted Slurm submission requires scheduler '
                           'reconciliation before any new probes are submitted.')
    if max_new_jobs is not None and (isinstance(max_new_jobs, bool)
            or not isinstance(max_new_jobs, int) or max_new_jobs < 0):
        raise ValueError('max_new_jobs must be nonnegative.')
    if live_job_ids is None:
        queued = subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                                text=True, check=True, timeout=30)
        live_job_ids = queued.stdout.split()
    live = set(live_job_ids)
    # The shared controller already audited under the same held locks and
    # scheduler snapshot; avoid rehashing every large checkpoint twice.
    manifest = audited_manifest if audited_manifest is not None else _audit(path, live_job_ids=live)
    if audited_manifest is None:
        _recover_slurm_startup_failures(path, manifest, live)
    active = sum(str(job.get('job_id')) in live for job in manifest['jobs'])
    # Validate one native probe before widening a fresh manifest's launch.
    # This catches shared input/runtime failures without duplicating them
    # across every basis. It adds no jobs and retains the configured cap.
    startup = None
    if not startup_checked and not any(job['status'] == 'complete' for job in manifest['jobs']):
        startup = manifest['jobs'][0]['name']
        max_jobs = min(max_jobs, 1)
    for job in manifest['jobs']:
        record = Path(job['directory']) / 'execution.json'
        if job.get('job_id') and str(job['job_id']) not in live and record.is_file():
            execution = json.loads(record.read_text())
            if execution.get('status') == 'failed':
                raise RuntimeError(execution.get('error') or f'{job["name"]}: native execution failed.')
    submitted = 0
    for job in manifest['jobs']:
        if startup is not None and job['name'] != startup:
            continue
        if (job['status'] == 'complete' or job.get('job_id') or active >= max_jobs
                or max_new_jobs is not None and submitted >= max_new_jobs):
            continue
        if job.get('submission_status') == 'submitting':
            raise RuntimeError(f'{job["name"]}: interrupted Slurm submission '
                               'requires scheduler reconciliation before resubmission.')
        controls = manifest.get('reference_solver_options', {})
        if controls and any(vrc_multireference_energy_contract(c['method'])
                            for c in job.get('correlation_checks', [])):
            validate_mr_correlation_options(controls)
            if set(controls) - {'max_reference_vectors', 'max_reference_iterations'}:
                raise ValueError('Only reference solver limits may override a prepared campaign.')
            for key, value in controls.items():
                # Explicit job retry controls take precedence over campaign defaults.
                job.setdefault('retry_' + key, value)
            inp = Path(job['directory']) / (job['name'] + '.inp')
            _reference_solver_override(inp.read_text(), job.get('retry_max_reference_vectors'),
                                        job.get('retry_max_reference_iterations'))
        if manifest.get('solver_options'):
            job.setdefault('retry_solver_options', deepcopy(manifest['solver_options']))
            for key in ('max_reference_vectors', 'max_reference_iterations'):
                if 'retry_' + key in job and key in job['retry_solver_options'].get('correlation', {}):
                    job['retry_' + key] = job['retry_solver_options']['correlation'][key]
            apply_solver_options((Path(job['directory']) / (job['name'] + '.inp')).read_text(),
                                 job['retry_solver_options'])
        predecessor = job.get('orbital_predecessor')
        if job.get('recovered_scf'):
            import shutil
            recovery = job['recovered_scf']
            source = Path(recovery['path'])
            if file_sha256(source) != recovery['sha256']:
                raise RuntimeError('Archived SCF checkpoint changed before recovery.')
            shutil.copyfile(source, Path(job['directory']) / (job['name'] + '.wfu'))
        elif predecessor:
            prior = next(item for item in manifest['jobs'] if item['name'] == predecessor)
            if prior['status'] != 'complete':
                continue
            source = Path(prior['wfu'])
            if file_sha256(source) != prior['wfu_sha256']:
                raise RuntimeError('Reconnaissance predecessor WFU changed after audit.')
            import shutil
            checkpoint = Path(job['directory']) / (job['name'] + '.wfu')
            shutil.copyfile(source, checkpoint)
            inp = Path(job['directory']) / (job['name'] + '.inp')
            inp.write_text(inp.read_text().replace('{{PREVIOUS_WFU_SHA256}}', prior['wfu_sha256']))
            job['initial_wfu_sha256'] = prior['wfu_sha256']
            job['input_sha256'] = file_sha256(inp)
        # Existing snapshotted campaigns gain the shared batch retry guard
        # without changing the chemically meaningful prepared input/hash.
        from kinbot.slurm_retry import wrap_script
        batch = Path(job['directory']) / (job['name'] + '.sbatch')
        batch.write_text(wrap_script(batch.read_text(), identity=job['name'],
                                     directories=[job['directory']]))
        options = _slurm_submission_options(batch.read_text(), manifest.get('slurm_excluded_nodes', []),
            **({'allow_inherited': manifest['slurm_inherited_exclusion_allow']}
               if manifest.get('slurm_inherited_exclusion_allow') else {}))
        shared = manifest.get('shared_policy_file')
        command = ['sbatch', '--export=ALL', '--parsable',
                   *([] if shared else options), job['name'] + '.sbatch']
        def journal():
            job['submission_status'] = 'submitting'
            job['submission_command'] = command
            job['submission_overrides'] = {key: os.environ[key] for key in
                ('SBATCH_MEM_PER_NODE', 'SBATCH_MEM_PER_CPU', 'SBATCH_TIMELIMIT',
                 'SBATCH_PARTITION') if key in os.environ}
            atomic_json(path, manifest)
        try:
            if shared:
                if run_policy.find(job['directory']) != Path(shared):
                    raise ValueError('Shared reconnaissance policy/path changed; no submission.')
                # Current pool exclusions supersede inherited script settings;
                # native allocation failures add exclusions without editing policy.
                ceiling = next((o.split('=', 1)[1] for o in options
                                if o.startswith('--time=')), None)
                result = run_policy.submit(command, cwd=job['directory'], before_submit=journal,
                    extra_excluded_nodes=manifest.get('slurm_excluded_nodes', []),
                    walltime_ceiling=ceiling, capture_output=True, text=True, timeout=30)
                if result is None:
                    break  # No admission: leave the probe prepared, not submitting.
            else:
                journal()
                result = subprocess.run(
                    command, cwd=job['directory'], capture_output=True, text=True, timeout=30)
        except OSError as error:
            if getattr(error, 'kinbot_submission_not_started', False):
                job['submission_status'] = 'retry_ready'
                atomic_json(path, manifest)
            raise
        except subprocess.TimeoutExpired as error:
            raise RuntimeError(f'{job["name"]}: sbatch timed out; submission outcome is '
                               'uncertain. Retain the submitting journal for scheduler '
                               'reconciliation.') from error
        if result.returncode:
            job['submission_status'] = 'failed'
            atomic_json(path, manifest)
            raise RuntimeError(f'{job["name"]}: Slurm submission failed: {result.stderr.strip()}')
        identifier = result.stdout.strip().split(';')[0]
        if not identifier.isdigit():
            raise RuntimeError(f'{job["name"]}: unrecognized Slurm job identifier; '
                               'reconcile the submitting journal before retrying.')
        job['job_id'] = identifier
        job['submission_status'] = 'submitted'
        active += 1
        submitted += 1
        atomic_json(path, manifest)
    return manifest


def run(manifest_path, *, max_jobs=8, interval=20.):
    with _manifest_lock(Path(manifest_path).expanduser().resolve().parent):
        return _run(manifest_path, max_jobs=max_jobs, interval=interval)


def _run(manifest_path, *, max_jobs=8, interval=20.):
    if not np.isfinite(interval) or not 0. < interval <= 60.:
        raise ValueError('Polling interval must be in (0, 60] seconds.')
    while True:
        manifest = _submit(manifest_path, max_jobs=max_jobs)
        queued = subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                                text=True, check=True, timeout=30)
        live = set(queued.stdout.split())
        if _finished(manifest_path, manifest, live):
            return manifest
        time.sleep(interval)


def _finished(manifest_path, manifest, live):
    if (manifest['status'] in ('complete_review_required', 'coverage_refinement_required')
            and not any(str(job.get('job_id')) in live for job in manifest['jobs'])):
        return True
    failures = []
    for job in manifest['jobs']:
        job.pop('terminal_failure', None)
        if job.get('job_id') and str(job['job_id']) not in live and job['status'] != 'complete':
            execution_path = Path(job['directory']) / 'execution.json'
            if execution_path.is_file():
                execution = json.loads(execution_path.read_text())
                if execution.get('status') == 'failed':
                    reason = execution.get('error') or 'Native execution failed.'
                    job['terminal_failure'] = (f'{job["name"]}: {reason} '
                                               f'Inspect {execution_path}.')
            # Allow shared-filesystem publication after queue exit. Until it
            # validates, a missing allocation still reserves campaign capacity.
            deadline = job.setdefault('missing_since', time.time())
            if not job.get('terminal_failure') and time.time() - deadline > 60.:
                reason = '; '.join(job.get('validation_errors', [])) or 'No validated result was published.'
                if job.get('startup_failure'):
                    reason = job['startup_failure']
                elif job.get('scheduler_record', {}).get('JobId') == str(job['job_id']):
                    record = job['scheduler_record']
                    reason += (f' Slurm: {record.get("JobState")}, '
                               f'{record.get("Reason")}, ExitCode={record.get("ExitCode")}, '
                               f'BatchHost={record.get("BatchHost")}.')
                job['terminal_failure'] = (f'{job["name"]}: queued probe ended without '
                    f'valid output: {reason} Inspect {execution_path} '
                    'and native output before retrying.')
            if job.get('terminal_failure'):
                failures.append(job['terminal_failure'])
        else:
            job.pop('missing_since', None)
    if failures:
        manifest['status'] = 'failed'
    atomic_json(Path(manifest_path), manifest)
    if failures:
        raise ProbeFailure('\n'.join(failures))
    return False


def run_many(manifest_paths, *, max_jobs=8, interval=20., startup_checked=False,
             capacity_manifests=(), keep_going=False, allocation_retries=0):
    """Poll independent manifests under one journaled campaign-wide job cap.

    A campaign may waive the per-manifest startup check only after its native
    pilots pass. Existing jobs from another campaign phase also consume slots.
    On probe failure, audit already submitted work before raising all observed
    failures. keep_going also submits other independent manifests. Opt-in
    allocation retries recover interrupted workers and explicit I/O failures.
    Scheduler/journal errors still stop at once.
    No worker threads or competing manifest writers are created.
    """
    if isinstance(max_jobs, bool) or not isinstance(max_jobs, int) or max_jobs < 1:
        raise ValueError('max_jobs must be positive.')
    if not np.isfinite(interval) or not 0. < interval <= 60.:
        raise ValueError('Polling interval must be in (0, 60] seconds.')
    if type(allocation_retries) is not int or not 0 <= allocation_retries <= 10:
        raise ValueError('allocation_retries must be an integer from 0 to 10.')
    paths = [Path(p).expanduser().resolve() for p in manifest_paths]
    if not paths or len(set(paths)) != len(paths):
        raise ValueError('Expected distinct reconnaissance manifests.')
    all_paths = list(dict.fromkeys(paths + [Path(p).expanduser().resolve() for p in capacity_manifests]))
    announced = set()
    with ExitStack() as locks:
        for parent in sorted({p.parent for p in all_paths}):
            locks.enter_context(_manifest_lock(parent))
        while True:
            snapshots = {p: json.loads(p.read_text()) for p in all_paths}
            if any(job.get('submission_status') == 'submitting'
                   for m in snapshots.values() for job in m['jobs'] + m.get('excluded_jobs', [])):
                raise RuntimeError('Interrupted Slurm submission requires scheduler reconciliation '
                                   'before any new campaign probes are submitted.')
            queued = subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                                    text=True, check=True, timeout=30)
            live = set(queued.stdout.split())
            # Save ephemeral scheduler records before hashing large completed
            # WFU/XML files. On sites without accounting they are the only
            # evidence of an allocation-level termination.
            for path, manifest in snapshots.items():
                _recover_slurm_startup_failures(path, manifest, live)
            # Re-audit every phase before allocating slots. An apparently
            # finished output remains untrusted while its allocation is live.
            snapshots = {p: _audit(p, live_job_ids=live) for p in all_paths}
            snapshots = {p: _retry_seed_projection_failures(p, m, live)
                         for p, m in snapshots.items()}
            if allocation_retries:
                snapshots = {p: _retry_allocation_failures(p, m, live, allocation_retries)
                             for p, m in snapshots.items()}
            # Persist the union in all phases, so later phases and resumed
            # controllers also avoid identified failing nodes. No global
            # environment or cluster configuration is changed.
            excluded = sorted({node for m in snapshots.values()
                               for node in m.get('slurm_excluded_nodes', [])})
            if excluded:
                for path, manifest in snapshots.items():
                    manifest['slurm_excluded_nodes'] = excluded
                    atomic_json(path, manifest)
            finished, failures = {}, []
            for path, manifest in snapshots.items():
                try:
                    finished[path] = _finished(path, manifest, live)
                except ProbeFailure as error:
                    finished[path] = False
                    failures.append(str(error))
            if failures and not keep_going:
                for failure in failures:
                    if failure not in announced:
                        print('Probe failed; no new jobs will be submitted. '
                              'Monitoring remaining allocations.\n' + failure, flush=True)
                        announced.add(failure)
                # Include other phases and excluded live allocations, but do
                # not wait for unsubmitted probes or failed predecessors.
                waiting = any(j.get('job_id') and (str(j['job_id']) in live or
                              j['status'] != 'complete' and not j.get('terminal_failure'))
                              for m in snapshots.values() for j in m['jobs'])
                waiting |= any(str(j.get('job_id')) in live
                               for m in snapshots.values() for j in m.get('excluded_jobs', []))
                if not waiting:
                    raise ProbeFailure('\n'.join(failures))
                time.sleep(interval)
                continue
            if failures and keep_going:
                # A failed manifest is isolated; never run its dependent
                # points. Other independent manifests retain the shared cap.
                for failure in failures:
                    if failure not in announced:
                        print('Probe failed; continuing independent manifests.\n' + failure, flush=True)
                        announced.add(failure)
                unresolved = {p for p, m in snapshots.items()
                              if any(j.get('terminal_failure') for j in m['jobs'])}
                busy = any(str(j.get('job_id')) in live or
                           j.get('job_id') and j['status'] != 'complete' and not j.get('terminal_failure')
                           for m in snapshots.values() for j in m['jobs'] + m.get('excluded_jobs', []))
                if not busy and all(finished[p] or p in unresolved for p in paths):
                    raise ProbeFailure('\n'.join(failures))
            if not failures and all(finished[p] for p in paths):
                return [snapshots[p] for p in paths]
            reserved = {str(j['job_id']) for m in snapshots.values() for j in m['jobs']
                        if j.get('job_id') and (str(j['job_id']) in live or
                            j['status'] != 'complete' and not j.get('terminal_failure'))}
            reserved.update(str(j['job_id']) for m in snapshots.values() for j in m.get('excluded_jobs', [])
                            if str(j.get('job_id')) in live)
            for path in paths:
                if finished[path] or any(j.get('terminal_failure') for j in snapshots[path]['jobs']):
                    continue
                before = {str(j['job_id']) for j in snapshots[path]['jobs'] if j.get('job_id')}
                report = _submit(path, max_jobs=max_jobs, max_new_jobs=max(0, max_jobs - len(reserved)),
                                 live_job_ids=live, startup_checked=startup_checked,
                                 audited_manifest=snapshots[path])
                new_ids = {str(j['job_id']) for j in report['jobs'] if j.get('job_id')} - before
                # Reserve immediately, even before squeue exposes a new job.
                reserved.update(new_ids)
                live.update(new_ids)
            time.sleep(interval)


def select_record(wfu, output, *, record, active_electrons, active_orbitals,
                  state_average, selected_roots, review_note, basis=None,
                  state_average_policy='reviewed', casscf_options=None, correlation_options=None,
                  active_indices=None, orbital_type='auto', align_active_orbitals=False):
    path = Path(wfu).expanduser().resolve()
    if not path.is_file() or not path.stat().st_size:
        raise FileNotFoundError(f'Reviewed Molpro WFU is unavailable: {path}.')
    if not str(review_note).strip():
        raise ValueError('A nonempty orbital review note is required.')
    if state_average_policy not in ('reviewed', 'channel_derived'):
        raise ValueError('Unknown state_average_policy.')
    model = validate_vrc_state_model({
        'orbital_selection': {
            'mode': 'molpro_record', 'path': str(path),
            'record': record,
            'sha256': file_sha256(path),
            'orbital_type': orbital_type,
            'align_active_orbitals': align_active_orbitals,
            **({'active_indices': active_indices} if active_indices is not None else {}),
        },
        'state_average': state_average,
        'selected_roots': selected_roots,
        'canonical_orbitals': True,
        **({'casscf_options': casscf_options} if casscf_options is not None else {}),
        **({'correlation_options': correlation_options} if correlation_options is not None else {}),
    }, active_electrons=active_electrons,
       active_orbitals=active_orbitals)
    active_space = {
            'mode': 'selected', 'active_electrons': active_electrons,
            'active_orbitals': active_orbitals,
            'source': str(review_note).strip(),
    }
    destination = Path(output).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if basis is None:
        payload = {'schema': 1, 'active_space': active_space,
                   'state_model': model}
    else:
        basis = str(basis).strip()
        if not basis:
            raise ValueError('basis cannot be empty.')
        if destination.is_file():
            payload = json.loads(destination.read_text())
            if (payload.get('schema') != 2
                    or payload.get('active_space', {}).get('active_electrons')
                    != active_electrons
                    or payload.get('active_space', {}).get('active_orbitals')
                    != active_orbitals):
                raise ValueError('Existing basis-specific state-model file '
                                 'has an incompatible schema or active space.')
        else:
            payload = {'schema': 2, 'active_space': active_space,
                       'basis_models': {}}
        payload['basis_models'][basis] = model
    if payload.get('state_average_policy', state_average_policy) != state_average_policy:
        raise ValueError('Existing state-model policy differs from this selection.')
    payload['state_average_policy'] = state_average_policy
    atomic_json(destination, payload)
    return destination


def main(argv=None):
    parser = ArgumentParser(prog='python -m kinbot.molpro_reconnaissance')
    commands = parser.add_subparsers(dest='command', required=True)
    prep = commands.add_parser('prepare')
    prep.add_argument('configuration')
    prep.add_argument('output')
    check = commands.add_parser('audit')
    check.add_argument('manifest')
    native = commands.add_parser('execute')
    native.add_argument('manifest')
    native.add_argument('name')
    enqueue = commands.add_parser('submit')
    enqueue.add_argument('manifest')
    enqueue.add_argument('--max-jobs', type=int, default=8)
    drive = commands.add_parser('run')
    drive.add_argument('manifest')
    drive.add_argument('--max-jobs', type=int, default=8)
    drive.add_argument('--interval', type=float, default=20.)
    drive.add_argument('--max-refinements', type=int, default=0)
    densify = commands.add_parser('refine')
    densify.add_argument('manifest')
    select = commands.add_parser('select')
    select.add_argument('wfu')
    select.add_argument('output')
    select.add_argument('--record', default='2140.2')
    select.add_argument('--basis')
    select.add_argument('--active-electrons', type=int, required=True)
    select.add_argument('--active-orbitals', type=int, required=True)
    select.add_argument('--active-indices', type=int, nargs='+',
                        help='One-based C1 indices in this exact WFU/record, in desired initial order.')
    select.add_argument('--orbital-type', choices=['auto', 'natural', 'canonical'], default='auto')
    select.add_argument('--align-active-orbitals', action='store_true',
                        help='Use DIAB only with a previously packed, reviewed active block.')
    select.add_argument('--state-average-json', required=True)
    select.add_argument('--selected-roots-json', default='{}')
    select.add_argument('--casscf-options-json', help='Reviewed CASSCF solver/tolerances to preserve in production.')
    select.add_argument('--correlation-options-json', help='Reviewed MR correlation iteration limits to preserve in production.')
    select.add_argument('--review-note', required=True)
    select.add_argument('--state-average-policy', choices=['reviewed', 'channel_derived'], default='reviewed')
    args = parser.parse_args(argv)
    if args.command == 'prepare':
        print(prepare(args.configuration, args.output))
    elif args.command == 'audit':
        print(json.dumps(audit(args.manifest), indent=2))
    elif args.command == 'execute':
        execute(args.manifest, args.name)
    elif args.command == 'submit':
        print(json.dumps(submit(args.manifest, max_jobs=args.max_jobs), indent=2))
    elif args.command == 'refine':
        print(refine(args.manifest))
    elif args.command == 'run':
        if args.max_refinements < 0:
            raise ValueError('max-refinements must be nonnegative.')
        for attempt in range(args.max_refinements + 1):
            result = run(args.manifest, max_jobs=args.max_jobs, interval=args.interval)
            if result['status'] != 'coverage_refinement_required' or attempt == args.max_refinements:
                break
            refine(args.manifest)
        print(json.dumps(result, indent=2))
    else:
        print(select_record(
            args.wfu, args.output, record=args.record,
            active_electrons=args.active_electrons,
            active_orbitals=args.active_orbitals,
            active_indices=args.active_indices, orbital_type=args.orbital_type,
            align_active_orbitals=args.align_active_orbitals,
            state_average=json.loads(args.state_average_json),
            selected_roots=json.loads(args.selected_roots_json),
            review_note=args.review_note, basis=args.basis,
            casscf_options=json.loads(args.casscf_options_json) if args.casscf_options_json else None,
            correlation_options=json.loads(args.correlation_options_json) if args.correlation_options_json else None,
            state_average_policy=args.state_average_policy))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
