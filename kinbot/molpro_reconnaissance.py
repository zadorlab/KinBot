"""Prepare and audit basis-specific Molpro orbital reconnaissance jobs.

The utility deliberately stops before choosing an active space.  It produces
ROHF/RHF orbital listings, Molden files, and hashable WFU records at every
requested geometry.  A reviewed Molpro record can then be converted to the
``vrc_tst_state_model`` contract with the ``select`` command.
"""

from argparse import ArgumentParser
from contextlib import contextmanager
from copy import deepcopy
import fcntl
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time

import numpy as np
from ase.data import atomic_numbers

from kinbot.anl.site import render_site_setup
from kinbot.molpro import (validate_vrc_state_model, automatic_vrc_state_model,
                           _render_vrc_state_model, full_valence_active_space,
                           validate_input_blocks)
from kinbot.theory import fortran_threshold
from kinbot.execution import atomic_json
from kinbot.anl.runtime import qc_runtime_environment, cleanup_qc_runtime


_SAFE = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')


@contextmanager
def _manifest_lock(root):
    root = Path(root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    with (root / '.reconnaissance.lock').open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError('Another reconnaissance controller is updating this directory; '
                               'read reconnaissance.json for live progress.') from error
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
        if (not isinstance(distances, list) or len(distances) < 4
                or any(isinstance(value, bool) or not isinstance(value, (int, float))
                       or not np.isfinite(value) or value <= 0. for value in distances)):
            raise ValueError('Each coordinate needs at least four finite positive distances.')
        distances = sorted({float(f'{value:.10g}') for value in distances})
        if len(distances) < 4:
            raise ValueError('Each coordinate needs at least four distinct distances.')
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
    for key in ('memory_mb', 'queue_job_limit'):
        if key in data and (isinstance(data[key], bool) or not isinstance(data[key], int)
                            or data[key] < 1):
            raise ValueError(f'{key} must be a positive integer.')
    if not isinstance(data.get('probe_active_space', False), bool):
        raise ValueError('probe_active_space must be a boolean.')
    if data.get('partition') is not None and not _SAFE.fullmatch(str(data['partition'])):
        raise ValueError('partition must be a safe Slurm partition name.')
    if 'coordinates' in data:
        if 'anchors' in data:
            raise ValueError('Configure coordinates or legacy anchors, not both.')
        data['anchors'] = _coordinates(data, path)
    anchors = data.get('anchors')
    if not isinstance(anchors, list) or not anchors:
        raise ValueError('anchors must be a nonempty list.')
    names = set()
    for anchor in anchors:
        allowed = {'name', 'bond', 'moving_atoms', 'distance_angstrom', '_channel', '_context',
                   '_rotation', '_orientation'}
        if not isinstance(anchor, dict) or set(anchor) - allowed \
                or not _SAFE.fullmatch(str(anchor.get('name', ''))):
            raise ValueError('Each anchor requires a safe name.')
        if anchor['name'] in names:
            raise ValueError('Anchor names must be unique.')
        names.add(anchor['name'])
        if set(anchor) == {'name'}:
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
                or set(controls) - {'include_light_atoms', 'additional_centers'}
                or ('include_light_atoms' in controls and not isinstance(controls['include_light_atoms'], bool))):
            raise ValueError('probe_space accepts include_light_atoms and additional_centers.')
        # Hydrogen cleavage needs its independent 1s projection, including
        # the corresponding antibonding space. Heavy-bond probes retain the
        # compact default, with a user-controlled all-atom extension.
        hydrogen_coordinate = any(any(symbols[index - 1] in ('H', 'He')
                                     for index in anchor.get('bond', [])) for anchor in anchors)
        data['_probe_space'] = full_valence_active_space(
            symbols, data['charge'], additional_centers=controls.get('additional_centers'),
            include_light_atoms=controls.get('include_light_atoms', hydrogen_coordinate))
    data = deepcopy(data)
    data['_source'] = str(path)
    data['_electrons'] = electrons
    return data


def _anchor_geometry(data, anchor):
    geometry = np.asarray(data['geometry_angstrom'], dtype=float).copy()
    if set(anchor) == {'name'}:
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
    context = anchor.get('_context', {})
    multiplicities = context.get('multiplicities', [data['multiplicity']])
    multiplicities = [data['multiplicity'], *[mult for mult in multiplicities
                                             if mult != data['multiplicity']]]
    return automatic_vrc_state_model(
        data['symbols'], data['charge'], multiplicities,
        roots=context.get('roots'), active_space=data['_probe_space'])


def _deck(data, basis, anchor, stem, previous=None):
    geometry = _anchor_geometry(data, anchor)
    rows = '\n'.join(
        f'{symbol} {x:.12f} {y:.12f} {z:.12f}'
        for symbol, (x, y, z) in zip(data['symbols'], geometry))
    scf = fortran_threshold(data.get('scf_convergence', 1.e-10),
                            'scf_convergence')
    basis_input = f'basis={basis}\n'
    deck = (
        f'***,{data["title"]} orbital reconnaissance: {stem}\n'
        f'memory,{int(data.get("memory_mw", 1000))},m\n'
        f'file,2,{stem}.wfu\n'
        'angstrom\norient,noorient\nsymmetry,nosym\n'
        f'geometry={{\n{len(data["symbols"])}\n{stem}\n{rows}\n}}\n'
        + basis_input +
        f'gthresh,energy={scf},coeff={scf}\n'
        'gprint,orbitals,orben,civector\n'
        f'{{rhf;wf,{data["_electrons"]},1,{data["multiplicity"] - 1};'
        'orbital,2100.2}\n'
        'gprint,orbitals,orben,civector\n'
        f'put,molden,{stem}.molden\n'
        'reconenergy=energy\nshow,reconenergy\n')
    if data.get('probe_active_space', False):
        space = data['_probe_space']
        model = _probe_model(data, anchor)
        if previous is not None:
            model.pop('orbital_selection', None)
        avas, multi, _ = _render_vrc_state_model(
            model, nelectron=data['_electrons'],
            active_electrons=space['active_electrons'],
            active_orbitals=space['active_orbitals'], spin=data['multiplicity'] - 1)
        multi = multi.replace('{multi,', '{multi,cas-ci;').replace(
            'canonical,2140.2', 'natorb,2140.2')
        if previous is not None:
            multi = multi.replace('{multi,cas-ci;',
                                  '{multi,cas-ci;\nstart,2140.2\ndiab,2140.2')
            deck = '! Initial orbital checkpoint SHA-256: {{PREVIOUS_WFU_SHA256}}\n' + deck
        deck += (avas + multi + f'{{put,molden,{stem}.natural.molden;orbital,2140.2,type=natural}}\n'
                 'reconenergy=energy(1)\nshow,reconenergy\n')
        for index in range(1, sum(entry['roots'] for entry in model['state_average']) + 1):
            deck += f'reconr{index}=energy({index})\nshow,reconr{index}\n'
    return deck + 'recondone=1\nshow,recondone\n'


def prepare(source, output):
    with _manifest_lock(output):
        return _prepare(source, output)


def _prepare(source, output):
    data = _load(source)
    root = Path(output).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    (root / 'site_setup.sh').write_text(
        render_site_setup({'molpro': {'molpro'}}))
    ppn = int(data.get('ppn', 8))
    partition = data.get('partition')
    memory_mb = int(data.get('memory_mb', ppn * data.get('memory_mw', 1000) * 8 * 1.25))
    walltime = str(data.get('walltime', '24:00:00'))
    if memory_mb < ppn * data.get('memory_mw', 1000) * 8:
        raise ValueError('Slurm memory_mb must cover every Molpro rank memory stack.')
    if not re.fullmatch(r'(?:\d+-)?\d+:\d{2}:\d{2}', walltime):
        raise ValueError('Reconnaissance walltime must be a Slurm time string.')
    records = []
    manifest_file = root / 'reconnaissance.json'
    previous = (json.loads(manifest_file.read_text()) if manifest_file.is_file() else {})
    previous_jobs = {job['name']: job for job in previous.get('jobs', [])}
    predecessors = {}
    for basis in data['bases']:
        for anchor in data['anchors']:
            stem = f'{basis["name"]}__{anchor["name"]}'
            directory = root / stem
            directory.mkdir(exist_ok=True)
            inp = directory / f'{stem}.inp'
            channel = anchor.get('_channel')
            old_job = previous_jobs.get(stem, {})
            if old_job and (not inp.is_file() or hashlib.sha256(inp.read_bytes()).hexdigest()
                            != old_job['input_sha256']):
                raise ValueError(f'{stem}: existing probe input changed after preparation.')
            # Refinement inserts points into the distance ordering. Existing
            # probes must retain the checkpoint they actually consumed.
            predecessor = (old_job.get('orbital_predecessor') if old_job else
                           predecessors.get((basis['name'], channel))
                           if channel and data.get('probe_active_space') else None)
            deck = _deck(data, basis['basis'], anchor, stem, previous=predecessor)
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
                f'#SBATCH -N 1\n#SBATCH -n {ppn}\n#SBATCH -c 1\n'
                + (f'#SBATCH -p {partition}\n' if partition else '')
                + f'#SBATCH --mem={memory_mb}M\n#SBATCH --time={walltime}\n'
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
            records.append({
                **old_job,
                'name': stem, 'basis': basis['basis'], 'anchor': anchor['name'],
                'directory': str(directory),
                'input_sha256': hashlib.sha256(inp.read_bytes()).hexdigest(),
                'template_sha256': digest,
                'orbital_predecessor': predecessor,
                'channel': anchor.get('_channel'),
                'orientation': anchor.get('_orientation'),
                'distance_angstrom': anchor.get('distance_angstrom'),
                'spin_context': anchor.get('_context'),
                'root_states': ([{'multiplicity': entry['multiplicity'], 'root': root}
                    for entry in _probe_model(data, anchor)['state_average']
                    for root in range(1, entry['roots'] + 1)]
                    if data.get('probe_active_space') else []),
                'natural_orbital_probe': bool(data.get('probe_active_space', False)),
                'probe_active_space': data.get('_probe_space'),
                'basis_source': 'library',
                'symbols': data['symbols'],
                'ppn': ppn,
            })
            if channel:
                predecessors[(basis['name'], channel)] = stem
    manifest = {
        'schema': 1, 'source': data['_source'], 'title': data['title'],
        'electron_count': data['_electrons'], 'jobs': records,
        'scan_policy': data.get('scan_policy'),
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
    directory = Path(job['directory'])
    inp = directory / f'{name}.inp'
    if hashlib.sha256(inp.read_bytes()).hexdigest() != job['input_sha256']:
        raise ValueError(f'{name}: input changed before native execution.')
    validate_input_blocks(inp.read_text(), name)
    env, runtime = qc_runtime_environment('molpro', 'molpro', work_directory=directory)
    command = ['molpro', '-d', runtime['scratch_directory'], '-W', str(directory), '-t', '1',
               '-n', str(job['ppn']), inp.name]
    execution = {'schema': 1, 'status': 'running', 'command': command,
                 'input_sha256': job['input_sha256'], 'runtime': runtime}
    record = directory / 'execution.json'
    atomic_json(record, execution)
    try:
        result = subprocess.run(command, cwd=directory, env=env, check=False)
        execution['returncode'] = result.returncode
        if result.returncode:
            output = directory / (name + '.out')
            detail = ('\n' + '\n'.join(output.read_text(errors='replace').splitlines()[-20:])
                      if output.is_file() else f' No native output found at {output}.')
            raise RuntimeError(f'{name}: Molpro exited with status {result.returncode}.{detail}')
        execution['status'] = 'complete'
    except BaseException as error:
        execution.update(status='failed', error=str(error))
        raise
    finally:
        atomic_json(record, execution)
        cleanup_qc_runtime(runtime)


def audit(manifest_path):
    with _manifest_lock(Path(manifest_path).expanduser().resolve().parent):
        return _audit(manifest_path)


def _audit(manifest_path):
    path = Path(manifest_path).expanduser().resolve()
    manifest = json.loads(path.read_text())
    complete = True
    for job in manifest['jobs']:
        for key in ('energy_hartree', 'root_energies_hartree', 'state_energies_hartree',
                    'natural_occupations', 'natural_molden_sha256',
                    'orbital_error'):
            job.pop(key, None)
        directory = Path(job['directory'])
        stem = job['name']
        output = directory / f'{stem}.out'
        wfu = directory / f'{stem}.wfu'
        molden = directory / f'{stem}.molden'
        natural = directory / f'{stem}.natural.molden'
        inp = directory / f'{stem}.inp'
        output_text = output.read_text(errors='replace') if output.is_file() else ''
        normal = (inp.is_file()
                  and hashlib.sha256(inp.read_bytes()).hexdigest() == job['input_sha256']
                  and 'Molpro calculation terminated' in output_text
                  and re.search(r'SETTING\s+RECONDONE\s*=\s*1(?:\.0*)?\s*(?:\n|$)', output_text)
                  and 'ERROR EXIT' not in output_text)
        job['status'] = ('complete' if normal and wfu.is_file()
                         and wfu.stat().st_size and molden.is_file()
                         and molden.stat().st_size else 'incomplete')
        energies = re.findall(r'SETTING\s+RECONENERGY\s*=\s*'
                              r'([-+0-9.]+(?:[EDed][-+]?\d+)?)\s+(?:AU|HARTREE)\b',
                              output_text, re.IGNORECASE)
        if energies:
            value = float(energies[-1].replace('D', 'E').replace('d', 'e'))
            if np.isfinite(value):
                job['energy_hartree'] = value
        if job.get('channel') and 'energy_hartree' not in job:
            job['status'] = 'incomplete'
        if job.get('natural_orbital_probe'):
            roots = {int(index): float(value.replace('D', 'E').replace('d', 'e'))
                     for index, value in re.findall(r'SETTING\s+RECONR(\d+)\s*=\s*'
                         r'([-+0-9.]+(?:[EDed][-+]?\d+)?)\s+(?:AU|HARTREE)\b',
                         output_text, re.IGNORECASE)}
            expected = sum((job.get('spin_context') or {}).get('roots', {}).values())
            if (not roots or any(not np.isfinite(value) for value in roots.values())
                    or (expected and sorted(roots) != list(range(1, expected + 1)))):
                job['status'] = 'incomplete'
            else:
                job['root_energies_hartree'] = [roots[index] for index in sorted(roots)]
                if job.get('root_states'):
                    if len(job['root_states']) != len(roots):
                        job['status'] = 'incomplete'
                    else:
                        job['state_energies_hartree'] = [
                            {**state, 'energy_hartree': value}
                            for state, value in zip(job['root_states'], job['root_energies_hartree'])]
            occupations = ([] if not natural.is_file() else
                           [float(value.replace('D', 'E').replace('d', 'e'))
                            for value in re.findall(r'(?im)^\s*Occup\s*=\s*'
                                r'([-+0-9.]+(?:[EDed][-+]?\d+)?)', natural.read_text())])
            if (not occupations or any(not np.isfinite(value) or not 0. <= value <= 2.
                                       for value in occupations)
                    or abs(sum(occupations) - manifest['electron_count']) > .02):
                job['status'] = 'incomplete'
            else:
                job['natural_occupations'] = occupations
                job['natural_molden_sha256'] = hashlib.sha256(natural.read_bytes()).hexdigest()
        complete &= job['status'] == 'complete'
        if job['status'] == 'complete':
            job['output_sha256'] = hashlib.sha256(output.read_bytes()).hexdigest()
            job['wfu'] = str(wfu.resolve())
            job['wfu_sha256'] = hashlib.sha256(wfu.read_bytes()).hexdigest()
            job['molden'] = str(molden.resolve())
            job['molden_sha256'] = hashlib.sha256(molden.read_bytes()).hexdigest()
        else:
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
                          'Molden/output pair; no active orbital identity is '
                          'selected automatically.')
    atomic_json(path, manifest)
    return manifest


def _coverage(manifest):
    """Check asymptotic plateaus and flag intervals needing denser probes."""
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


def retry_failed(manifest_path):
    """Archive known failed native attempts after every allocation has left.

    Completed calculations and immutable inputs/checkpoint chains stay intact.
    Unknown submissions and parser-only failures require reconciliation rather
    than an automatic duplicate calculation.
    """
    import shutil
    path = Path(manifest_path).expanduser().resolve()
    with _manifest_lock(path.parent):
        manifest = _audit(path)
        queued = subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                                text=True, check=True, timeout=30)
        live = set(queued.stdout.split())
        if any(str(job.get('job_id')) in live for job in manifest['jobs']):
            raise RuntimeError('Wait for all manifest allocations to leave before retrying.')
        if any(job.get('submission_status') == 'submitting' for job in manifest['jobs']):
            raise RuntimeError('Reconcile uncertain Slurm submissions before retrying.')
        retried = []
        for job in manifest['jobs']:
            if job['status'] == 'complete':
                continue
            directory = Path(job['directory'])
            record = directory / 'execution.json'
            if not record.is_file():
                continue
            execution = json.loads(record.read_text())
            if execution.get('status') != 'failed' or not job.get('job_id'):
                continue
            inp = directory / (job['name'] + '.inp')
            if hashlib.sha256(inp.read_bytes()).hexdigest() != job['input_sha256']:
                raise ValueError('Failed probe input changed; do not retry under the old manifest.')
            attempt = int(job.get('attempt', 1))
            archive = directory / 'attempts' / str(attempt)
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
                'archive': str(archive), 'error': execution.get('error'),
                'input_sha256': job['input_sha256']})
            for key in ('job_id', 'missing_since'):
                job.pop(key, None)
            job.update(attempt=attempt + 1, submission_status='retry_ready', status='incomplete')
            retried.append(job['name'])
            atomic_json(path, manifest)
        return {'manifest': str(path), 'retry_ready': retried,
                'completed_preserved': sum(job['status'] == 'complete' for job in manifest['jobs'])}


def _submit(manifest_path, *, max_jobs=8):
    if isinstance(max_jobs, bool) or not isinstance(max_jobs, int) or max_jobs < 1:
        raise ValueError('max_jobs must be positive.')
    path = Path(manifest_path).expanduser().resolve()
    manifest = _audit(path)
    queued = subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                            text=True, check=True)
    live = set(queued.stdout.split())
    active = sum(str(job.get('job_id')) in live for job in manifest['jobs'])
    # Validate one native probe before widening a fresh manifest's launch.
    # This catches shared input/runtime failures without duplicating them
    # across every basis. It adds no jobs and retains the configured cap.
    startup = None
    if not any(job['status'] == 'complete' for job in manifest['jobs']):
        startup = manifest['jobs'][0]['name']
        max_jobs = min(max_jobs, 1)
    for job in manifest['jobs']:
        record = Path(job['directory']) / 'execution.json'
        if job.get('job_id') and str(job['job_id']) not in live and record.is_file():
            execution = json.loads(record.read_text())
            if execution.get('status') == 'failed':
                raise RuntimeError(execution.get('error') or f'{job["name"]}: native execution failed.')
    for job in manifest['jobs']:
        if startup is not None and job['name'] != startup:
            continue
        if job['status'] == 'complete' or job.get('job_id') or active >= max_jobs:
            continue
        if job.get('submission_status') == 'submitting':
            raise RuntimeError(f'{job["name"]}: interrupted Slurm submission '
                               'requires scheduler reconciliation before resubmission.')
        predecessor = job.get('orbital_predecessor')
        if predecessor:
            prior = next(item for item in manifest['jobs'] if item['name'] == predecessor)
            if prior['status'] != 'complete':
                continue
            source = Path(prior['wfu'])
            if hashlib.sha256(source.read_bytes()).hexdigest() != prior['wfu_sha256']:
                raise RuntimeError('Reconnaissance predecessor WFU changed after audit.')
            import shutil
            checkpoint = Path(job['directory']) / (job['name'] + '.wfu')
            shutil.copyfile(source, checkpoint)
            inp = Path(job['directory']) / (job['name'] + '.inp')
            inp.write_text(inp.read_text().replace('{{PREVIOUS_WFU_SHA256}}', prior['wfu_sha256']))
            job['initial_wfu_sha256'] = prior['wfu_sha256']
            job['input_sha256'] = hashlib.sha256(inp.read_bytes()).hexdigest()
        job['submission_status'] = 'submitting'
        atomic_json(path, manifest)
        result = subprocess.run(['sbatch', '--export=ALL', '--parsable', job['name'] + '.sbatch'],
                                cwd=job['directory'], capture_output=True, text=True)
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
                                text=True, check=True)
        live = set(queued.stdout.split())
        if (manifest['status'] in ('complete_review_required', 'coverage_refinement_required')
                and not any(str(job.get('job_id')) in live for job in manifest['jobs'])):
            return manifest
        for job in manifest['jobs']:
            if job.get('job_id') and str(job['job_id']) not in live and job['status'] != 'complete':
                # Account for shared-filesystem publication after queue exit.
                deadline = job.setdefault('missing_since', time.time())
                if time.time() - deadline > 60.:
                    raise RuntimeError(f'{job["name"]}: queued probe ended without '
                                       'valid output; inspect its stdout/stderr before retrying.')
            else:
                job.pop('missing_since', None)
        atomic_json(Path(manifest_path), manifest)
        time.sleep(interval)


def select_record(wfu, output, *, record, active_electrons, active_orbitals,
                  state_average, selected_roots, review_note, basis=None):
    path = Path(wfu).expanduser().resolve()
    if not path.is_file() or not path.stat().st_size:
        raise FileNotFoundError(f'Reviewed Molpro WFU is unavailable: {path}.')
    if not str(review_note).strip():
        raise ValueError('A nonempty orbital review note is required.')
    model = validate_vrc_state_model({
        'orbital_selection': {
            'mode': 'molpro_record', 'path': str(path),
            'record': record,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        },
        'state_average': state_average,
        'selected_roots': selected_roots,
        'canonical_orbitals': True,
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
    drive.add_argument('--max-refinements', type=int, default=3)
    densify = commands.add_parser('refine')
    densify.add_argument('manifest')
    select = commands.add_parser('select')
    select.add_argument('wfu')
    select.add_argument('output')
    select.add_argument('--record', default='2140.2')
    select.add_argument('--basis')
    select.add_argument('--active-electrons', type=int, required=True)
    select.add_argument('--active-orbitals', type=int, required=True)
    select.add_argument('--state-average-json', required=True)
    select.add_argument('--selected-roots-json', default='{}')
    select.add_argument('--review-note', required=True)
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
            state_average=json.loads(args.state_average_json),
            selected_roots=json.loads(args.selected_roots_json),
            review_note=args.review_note, basis=args.basis))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
