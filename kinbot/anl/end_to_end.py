"""Restartable CBH/ANL -> MESSPF -> MESS verification orchestration.

The individual QC, composite, CBH, handoff, MESSPF, ROTD_py, and MESS
layers remain independently usable.  This module joins already completed,
hash-verified records through a declarative manifest; it never substitutes a
published value for a calculated one and never launches electronic-structure
work itself.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

import numpy as np

from kinbot.anl.atct import _canonical_smiles
from kinbot.anl.cbh import (HARTREE_TO_KJ_MOL, formation_from_dict,
                            formation_species, solve_composite_records,
                            zero_k_from_composite_record)
from kinbot.anl.kinetics import audit_rate_tables
from kinbot.anl.literature import compare_formation, compare_run
from kinbot.anl.thermochemistry import audit_messpf_run
from kinbot.anl.validation import (assemble_profiled_anl0_f12,
                                   molecule_from_completed_run)
from kinbot.anl.workflow import (_verified_task_result,
                                 attach_task_vpt2_frequencies)
from kinbot.execution import (atomic_json as _write_json,
                              file_sha256 as _sha256,
                              resolve_path as _path)
from kinbot.rotdpy import number_of_states_file, read_result
from kinbot.stationary_pt import StationaryPoint
from kinbot.symmetry import calculate_symmetry


HARTREE_TO_KCAL_MOL = HARTREE_TO_KJ_MOL / 4.184


def _read_json(path, label):
    try:
        value = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f'{label} is not valid JSON: {path}.') from error
    if not isinstance(value, dict):
        raise ValueError(f'{label} must contain a JSON object.')
    return value


def _optional_path(base, value, label):
    return None if value is None else _path(base, value, label)


def _manifest(source):
    path = Path(source).expanduser().resolve()
    data = _read_json(path, 'End-to-end verification manifest')
    if data.get('schema') != 1:
        raise ValueError('End-to-end verification manifest needs schema 1.')
    for key in ('run_dir', 'kinbot_input', 'rotdpy', 'composites', 'cbh',
                'handoff', 'mess'):
        if key not in data:
            raise ValueError(f'End-to-end verification manifest lacks {key}.')
    return path, data


def _composite_value(record, key):
    try:
        value = float(record['components'][key]['value_hartree'])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f'Composite record lacks component {key}.') from error
    if not math.isfinite(value):
        raise ValueError(f'Composite component {key} is nonfinite.')
    return value


def component_breakdown(reaction, energy_files):
    """Return signed contributions from each record's actual expression.

    Composite records at one CBH node share an accepted method label, but
    they need not have identical component inventories.  Rank-exact zeros and
    lower fallback recipes can therefore coexist without fabricating missing
    calculations.  The CBH solver, rather than this diagnostic report,
    enforces the one-common-level rule for the reaction.
    """
    records = {smiles: _read_json(path, f'{smiles} composite')
               for smiles, path in energy_files.items()}

    def category(component):
        if component == 'vpt2_correction':
            return 'anharmonic_zpe_correction'
        if component == 'higher_order_correction' or component.startswith('hoe_'):
            return 'higher_order'
        return component

    def terms(record):
        declared = record.get('electronic_terms')
        zero_point = record.get('zero_point_terms')
        if isinstance(declared, list) and isinstance(zero_point, list):
            return declared + zero_point
        # Read old completed records that predate explicit expression lists.
        components = record.get('components', {})
        legacy = []
        for key in ('reference_cbs', 'core_valence_cbs',
                    'scalar_relativistic', 'dboc', 'spin_orbit',
                    'harmonic_zpe', 'vpt2_correction',
                    'higher_order_correction'):
            if key in components:
                legacy.append({'component': key, 'coefficient': 1})
        for key, coefficient in (
                ('hoe_high', 1), ('hoe_low', -1),
                ('hoe_tz_high', 1), ('hoe_tz_low', -1),
                ('hoe_dz_high', 1), ('hoe_dz_low', -1)):
            if key in components:
                legacy.append({'component': key, 'coefficient': coefficient})
        return legacy

    values = {label: [] for label in (
        'reference_cbs', 'higher_order', 'core_valence_cbs',
        'scalar_relativistic', 'dboc', 'spin_orbit', 'harmonic_zpe',
        'anharmonic_zpe_correction')}
    for smiles, stoichiometric_coefficient in reaction.stoichiometry.items():
        record = records[smiles]
        components = record.get('components', {})
        used = set()
        for term in terms(record):
            if (not isinstance(term, dict)
                    or set(term) != {'component', 'coefficient'}
                    or term['component'] in used
                    or term['component'] not in components
                    or term['coefficient'] not in (-1, 1)):
                raise ValueError(
                    f'{smiles}: invalid composite expression in breakdown.')
            used.add(term['component'])
            label = category(term['component'])
            values.setdefault(label, []).append(
                stoichiometric_coefficient * term['coefficient']
                * _composite_value(record, term['component']))
        if not used:
            raise ValueError(f'{smiles}: composite expression is empty.')

    result = {}
    for label, contributions in values.items():
        hartree = math.fsum(contributions)
        result[label] = {
            'hartree': hartree,
            'kcal_mol': hartree * HARTREE_TO_KCAL_MOL,
        }
    total = math.fsum(item['hartree'] for item in result.values())
    result['total_zero_k_reaction'] = {
        'hartree': total, 'kcal_mol': total * HARTREE_TO_KCAL_MOL}
    return result


def _assemble_composites(base, specification):
    if not isinstance(specification, dict) or not specification:
        raise ValueError('composites must be a nonempty object.')
    records = {}
    for name, entry in specification.items():
        if not isinstance(entry, dict):
            raise ValueError(f'composites.{name} must be an object.')
        smiles = entry.get('smiles')
        output = _path(base, entry.get('path'), f'composites.{name}.path')
        assembly = entry.get('assembly')
        if assembly is not None:
            if not isinstance(assembly, dict):
                raise ValueError(f'composites.{name}.assembly must be an object.')
            try:
                result = assemble_profiled_anl0_f12(
                    _path(base, assembly.get('interface_run'),
                          f'{name}.interface_run'),
                    _path(base, assembly.get('higher_order_run'),
                          f'{name}.higher_order_run'),
                    _path(base, assembly.get('corrections_run'),
                          f'{name}.corrections_run'),
                    state_id=assembly.get('state_id'),
                    spin_orbit_hartree=assembly.get('spin_orbit_hartree'),
                    spin_orbit_source=assembly.get('spin_orbit_source'),
                    spin_orbit_backend=assembly.get(
                        'spin_orbit_backend', 'manual'),
                    vpt2_review=_optional_path(
                        base, assembly.get('vpt2_review'),
                        f'{name}.vpt2_review'),
                    base_run=_optional_path(base, assembly.get('base_run'),
                                            f'{name}.base_run'))
            except (ValueError, RuntimeError) as exc:
                raise type(exc)(f'composites.{name}: {exc}') from exc
            _write_json(output, result)
        final = zero_k_from_composite_record(smiles, output)
        l2_geometry = None
        if assembly is not None:
            l2_geometry = molecule_from_completed_run(
                _path(base, assembly.get('interface_run'),
                      f'{name}.interface_run'), 'l2_geometry')
        records[name] = {
            'smiles': smiles, 'path': output, 'method': final.method,
            'charge': final.charge, 'multiplicity': final.multiplicity,
            'zero_k_hartree': final.hartree,
            'l2_geometry': l2_geometry,
        }
        if assembly is not None:
            from kinbot.anl.thermal import handoff_options
            thermal = handoff_options(_path(base,
                assembly.get('base_run') or assembly.get('interface_run'), f'{name}.thermal_run'))
            if thermal is not None:
                records[name]['thermal'] = thermal
    return records


def _solve_cbh(base, specification, composites, atct):
    if not isinstance(specification, dict) or not specification:
        raise ValueError('cbh must be a nonempty object.')
    solved = {}
    for name, entry in specification.items():
        if not isinstance(entry, dict):
            raise ValueError(f'cbh.{name} must be an object.')
        mapping = entry.get('energies')
        if not isinstance(mapping, dict) or not mapping:
            raise ValueError(f'cbh.{name}.energies must be a nonempty object.')
        try:
            energy_files = {}
            for smiles, selected in mapping.items():
                if isinstance(selected, str):
                    energy_files[smiles] = composites[selected]['path']
                elif isinstance(selected, dict) and selected:
                    energy_files[smiles] = {
                        method: composites[record]['path']
                        for method, record in selected.items()}
                else:
                    raise ValueError(
                        f'cbh.{name}.energies.{smiles} must name one composite '
                        'or a method-to-composite mapping.')
        except KeyError as error:
            raise ValueError(
                f'cbh.{name} names unknown composite {error.args[0]!r}.') from error
        reaction, formation = solve_composite_records(
            entry.get('target_smiles'), entry.get('rung', 0), energy_files,
            atct_version=atct.get('version'),
            atct_cache=_path(base, atct.get('cache'), 'atct.cache'),
            reference_ids=entry.get('reference_ids', {}),
            refresh_atct=bool(atct.get('refresh', False)))
        output = _path(base, entry.get('output'), f'cbh.{name}.output')
        payload = {
            'schema': 1, 'status': 'complete',
            'reaction': {
                'rung': reaction.rung,
                'target_smiles': reaction.target_smiles,
                'stoichiometry': dict(reaction.stoichiometry),
                'formulas': dict(reaction.formulas),
                'states': {key: list(value)
                           for key, value in reaction.states.items()},
            },
            'formation': {
                **asdict(formation),
                'formation_0k_kcal_mol':
                    formation.formation_0k_kj_mol / 4.184,
                'reaction_energy_0k_kcal_mol':
                    formation.reaction_energy_0k_kj_mol / 4.184,
            },
            'component_breakdown': component_breakdown(reaction, {
                smiles: (paths[formation.method]
                         if isinstance(paths, dict) else paths)
                for smiles, paths in energy_files.items()
                if smiles in reaction.stoichiometry}),
        }
        calculated = payload['component_breakdown'][
            'total_zero_k_reaction']['kcal_mol']
        if not math.isclose(
                calculated, payload['formation'][
                    'reaction_energy_0k_kcal_mol'],
                rel_tol=0., abs_tol=2e-8):
            raise ValueError(
                f'cbh.{name}: component reaction sum is inconsistent.')
        _write_json(output, payload)
        solved[name] = {'path': output, **payload}
    return solved


def _vpt2_entry(base, value, label):
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ValueError(f'{label} must be an object.')
    result = dict(value)
    result['run'] = str(_path(base, value.get('run'), f'{label}.run'))
    if value.get('review') is not None:
        result['review'] = str(_path(
            base, value['review'], f'{label}.review'))
    return result


def _write_handoff(base, specification, composites, solved):
    if not isinstance(specification, dict):
        raise ValueError('handoff must be an object.')
    species = {}
    for routing, entry in (specification.get('species') or {}).items():
        if not isinstance(entry, dict):
            raise ValueError(f'handoff species {routing} must be an object.')
        try:
            record = {
                'composite': str(composites[entry['composite']]['path']),
                'formation': str(solved[entry['formation']]['path']),
            }
        except KeyError as error:
            raise ValueError(
                f'handoff species {routing} has unknown record '
                f'{error.args[0]!r}.') from error
        if entry.get('barrier_composite') is not None:
            try:
                record['barrier_composite'] = str(
                    composites[entry['barrier_composite']]['path'])
            except KeyError as error:
                raise ValueError(
                    f'handoff species {routing} has unknown barrier '
                    f'composite {error.args[0]!r}.') from error
        vpt2 = _vpt2_entry(base, entry.get('vpt2'),
                           f'handoff.species.{routing}.vpt2')
        if vpt2 is not None:
            record['vpt2'] = vpt2
        if composites[entry['composite']].get('thermal'):
            record['thermal'] = composites[entry['composite']]['thermal']
        if entry.get('require_l2_geometry_match'):
            geometry = composites[entry['composite']].get('l2_geometry')
            if geometry is None:
                raise ValueError(
                    f'handoff species {routing} has no accepted L2 geometry.')
            record['geometry'] = geometry
        species[str(routing)] = record
    transition_states = {}
    for name, entry in (specification.get('transition_states') or {}).items():
        if not isinstance(entry, dict) or entry.get('composite') not in composites:
            raise ValueError(f'handoff transition state {name} is invalid.')
        record = {'composite': str(composites[entry['composite']]['path'])}
        vpt2 = _vpt2_entry(base, entry.get('vpt2'),
                           f'handoff.transition_states.{name}.vpt2')
        if vpt2 is not None:
            record['vpt2'] = vpt2
        if composites[entry['composite']].get('thermal'):
            record['thermal'] = composites[entry['composite']]['thermal']
        transition_states[name] = record
    output = _path(base, specification.get('output'), 'handoff.output')
    payload = {
        'schema': 1, 'mode': 'cbh-anl',
        'strict': specification.get('strict', True),
        'species': species, 'transition_states': transition_states,
    }
    _write_json(output, payload)
    return output, payload


def _write_cbh_thermochemistry_models(base, specification, composites,
                                      solved, handoff_specification):
    """Build hash-traceable RRHO models for CBH references outside the PES.

    Species present in the live KinBot network retain their normal L2
    geometry, accepted VPT2 frequencies, and L2 hindered rotors.  This file
    supplies only the other species required by the recursive CBH reaction
    trees.  A rotatable external reference fails closed because its HIR model
    cannot be inferred from a geometry/frequency calculation alone.
    """
    required = set()
    for entry in solved.values():
        # Current formation records contain the full recursive reaction tree.
        # Decoding it here also validates every 0 K formation identity before
        # any finite-temperature reference model is staged.
        required.update(formation_species(
            formation_from_dict(entry['formation'])))
    network = set()
    for entry in (handoff_specification.get('species') or {}).values():
        if isinstance(entry, dict) and entry.get('composite') in composites:
            network.add(_canonical_smiles(
                composites[entry['composite']]['smiles']))
    external = {_canonical_smiles(smiles) for smiles in required} - network
    by_smiles = {}
    for name, record in composites.items():
        canonical = _canonical_smiles(record['smiles'])
        if canonical in by_smiles:
            raise ValueError(
                f'Duplicate composite thermochemistry identity {canonical!r}.')
        by_smiles[canonical] = (name, record)
    missing = sorted(smiles for smiles in external if smiles not in by_smiles)
    if missing:
        raise ValueError(
            'CBH reference thermochemistry lacks composite models for: '
            + ', '.join(missing))

    models = {}
    for canonical in sorted(external):
        name, record = by_smiles[canonical]
        geometry = record.get('l2_geometry')
        if not isinstance(geometry, dict):
            raise ValueError(
                f'{name}: an external CBH reference needs an accepted L2 geometry.')
        symbols = list(geometry.get('symbols', ()))
        positions = geometry.get('positions')
        point = StationaryPoint(
            name=name, charge=record['charge'], mult=record['multiplicity'],
            smiles=record['smiles'], natom=len(symbols), atom=symbols,
            geom=np.asarray(positions, dtype=float))
        point.characterize()
        calculate_symmetry(point)
        if point.dihed:
            raise ValueError(
                f'{name}: external CBH reference has rotors; provide it as a '
                'KinBot network species so its L2 hindered rotors are retained.')
        assembly = specification[name].get('assembly')
        if not isinstance(assembly, dict):
            raise ValueError(
                f'{name}: external CBH reference lacks assembly provenance.')
        frequency_source = {
            'kind': 'atomic_exact_zero',
            'geometry': geometry['source'],
        }
        frequencies = []
        if len(symbols) > 1:
            interface = _path(
                base, assembly.get('interface_run'), f'{name}.interface_run')
            _, _, task, execution, native, parsed = _verified_task_result(
                interface, 'gaussian_vpt2')
            harmonic = parsed.get('harmonic_fundamentals_cm_inverse')
            if (not isinstance(harmonic, list) or not harmonic
                    or any(not isinstance(value, (int, float))
                           or not math.isfinite(value) or value <= 0.
                           for value in harmonic)):
                raise ValueError(
                    f'{name}: accepted VPT2 task lacks positive harmonic modes.')
            point.reduced_freqs = list(harmonic)
            point.freq = list(harmonic)
            if record.get('thermal'):
                from kinbot.anl.thermal import attach_model
                attach_model(point, record['thermal'])
            review = _optional_path(
                base, assembly.get('vpt2_review'), f'{name}.vpt2_review')
            frequencies = list(attach_task_vpt2_frequencies(
                point, interface, 'gaussian_vpt2', review_file=review))
            frequency_source = dict(
                point.anl_thermochemistry_frequency_source)
            frequency_source['native_output'] = str(native)
            frequency_source['native_output_sha256'] = execution[
                'artifacts'][task['result_parser']['file']]
        thermal_geometry = getattr(point, 'anl_thermal_model', {}).get('geometry')
        if thermal_geometry:
            positions = thermal_geometry['positions']
        models[record['smiles']] = {
            'name': name, 'smiles': record['smiles'], 'symbols': symbols,
            'positions': positions, 'charge': record['charge'],
            'multiplicity': record['multiplicity'],
            'symmetry_factor': float(thermal_geometry['sigma_ext'] if thermal_geometry else point.sigma_ext),
            'frequencies_cm_inverse': frequencies,
            'hindered_rotor_source': 'none-no-rotatable-bonds',
            'source': {
                'composite': str(record['path']),
                'composite_sha256': _sha256(record['path']),
                'geometry': geometry['source'],
                'frequencies': frequency_source,
            },
        }
    output = base / 'cbh_thermochemistry_models.json'
    return _write_json(output, {
        'schema': 1, 'temperature_model': 'RRHO-with-accepted-VPT2-modes',
        'species': models})


def _write_mess_input(base, manifest, handoff, rotdpy_input,
                      thermochemistry_models=None):
    run_dir = _path(base, manifest.get('run_dir'), 'run_dir')
    source = _path(base, manifest.get('kinbot_input'), 'kinbot_input')
    parameters = _read_json(source, 'KinBot MESS source input')
    options = manifest.get('mess')
    if not isinstance(options, dict):
        raise ValueError('mess must be an object.')
    parameters.update({
        'me': 1,
        'run_me': 1,
        'me_code': 'mess',
        'run_messpf': 1,
        'anl_handoff_manifest': str(handoff),
        'rotdpy_directory': str(rotdpy_input.parent),
        'rotdpy_mess_mode': 'production',
        'rotdpy_run': 0,
        'vrc_tst_scan': 0,
        'vrc_tst_noscan': 0,
        # ReactionFinder must recreate the in-memory reaction objects that
        # MESS consumes. The profiled router makes this a read-only replay:
        # every L1 calculation must already exist, while new L2/HIR jobs are
        # still permitted.
        'reaction_search': 1,
        'profiled_l1_reuse_only': 1,
        'profiled_route_file': '.kinbot_theory_jobs.full_verification.json',
        # Reconstruct the stable species on their matching L2 surface and
        # generate the HIR data MESS needs. These inexpensive jobs are the
        # only electronic-structure work in this continuation.
        'high_level': 1,
        'allow_l2_without_conf': 1,
        'rotor_scan': 1,
        'do_clean': False,
    })
    if thermochemistry_models is not None:
        parameters['cbh_thermochemistry_models'] = str(
            thermochemistry_models)
    for key in ('mess_command', 'messpf_command', 'queue_name', 'ppn',
                'queue_job_limit', 'uq_max_runs', 'collider', 'epsilon',
                'epsilon_unit', 'sigma'):
        if key in options:
            parameters[key] = options[key]
    output = _path(base, options.get(
        'input', str(run_dir / 'full_cbh_anl_mess.json')), 'mess.input')
    if output.parent != run_dir:
        raise ValueError('The generated KinBot MESS input must be in run_dir.')
    _write_json(output, parameters)
    return run_dir, output


def prepare(source):
    manifest_path, data = _manifest(source)
    base = manifest_path.parent
    rotdpy = data['rotdpy']
    if not isinstance(rotdpy, dict):
        raise ValueError('rotdpy must be an object.')
    rotdpy_input = _path(base, rotdpy.get('input'), 'rotdpy.input')
    result = read_result(rotdpy_input, required_profile='production')
    composites = _assemble_composites(base, data['composites'])
    atct = data.get('atct')
    if not isinstance(atct, dict):
        raise ValueError('atct must be an object.')
    solved = _solve_cbh(base, data['cbh'], composites, atct)
    handoff, handoff_payload = _write_handoff(
        base, data['handoff'], composites, solved)
    thermochemistry_models = _write_cbh_thermochemistry_models(
        base, data['composites'], composites, solved, data['handoff'])
    run_dir, kinbot_input = _write_mess_input(
        base, data, handoff, rotdpy_input, thermochemistry_models)
    literature = _literature_checks(base, data.get('literature'))
    failed_literature = [
        result.get('benchmark', '<unknown>') for result in literature
        if result.get('status') != 'passed']
    if failed_literature:
        raise RuntimeError(
            'ANL/CBH literature preflight failed for: '
            + ', '.join(failed_literature) + '.')
    payload = {
        'schema': 1, 'status': 'prepared',
        'manifest': str(manifest_path),
        'run_dir': str(run_dir),
        'kinbot_input': str(kinbot_input),
        'rotdpy': {
            'input': str(rotdpy_input),
            'surface_count': result['surface_count'],
            'validation': result['validation'],
        },
        'composites': {
            key: {**value, 'path': str(value['path'])}
            for key, value in composites.items()},
        'cbh': {
            key: {
                'path': str(value['path']),
                'formation': value['formation'],
                'component_breakdown': value['component_breakdown'],
            } for key, value in solved.items()},
        'handoff': str(handoff),
        'cbh_thermochemistry_models': str(thermochemistry_models),
        'handoff_species': sorted(handoff_payload['species']),
        'literature': literature,
    }
    preparation = _path(base, data.get(
        'preparation_report', 'full_verification_preparation.json'),
        'preparation_report')
    _write_json(preparation, payload)
    return payload


def _literature_checks(base, entries):
    checks = []
    for index, entry in enumerate(entries or []):
        if not isinstance(entry, dict):
            raise ValueError(f'literature[{index}] must be an object.')
        reference = _path(base, entry.get('reference'),
                          f'literature[{index}].reference')
        kind = entry.get('kind')
        if kind == 'run':
            result = compare_run(
                entry.get('benchmark'),
                _path(base, entry.get('source'), f'literature[{index}].source'),
                reference=reference)
        elif kind == 'formation':
            source = _path(
                base, entry.get('source'), f'literature[{index}].source')
            method = entry.get('method', 'auto')
            if method == 'auto':
                method = _read_json(source, 'CBH formation result')[
                    'formation']['method']
            result = compare_formation(
                entry.get('benchmark'), source, reference=reference, method=method,
                absolute_tolerance_kcal_mol=float(entry.get(
                    'absolute_tolerance_kcal_mol', 0.1)))
        else:
            raise ValueError(f'literature[{index}] has unsupported kind {kind!r}.')
        required = entry.get('required_checks', ())
        if required:
            if (not isinstance(required, list)
                    or not all(isinstance(item, str) and item
                               for item in required)):
                raise ValueError(
                    f'literature[{index}].required_checks must be a list '
                    'of names.')
            missing = sorted(set(required) - set(result.get('checks', {})))
            if missing:
                raise ValueError(
                    f'literature[{index}] lacks required comparisons: '
                    f'{missing}.')
            result['required_checks'] = list(required)
        checks.append(result)
    return checks


def _audit_handoff(base, data, run_dir):
    """Prove which composite and CBH files the reconstructed network used."""
    manifest_path = _path(
        base, data['handoff'].get('output'), 'handoff.output')
    applied_path = run_dir / 'me' / 'anl_handoff.json'
    manifest = _read_json(manifest_path, 'ANL handoff manifest')
    applied = _read_json(applied_path, 'Applied ANL handoff')
    if (manifest.get('schema') != 1 or applied.get('schema') != 1
            or applied.get('status') != 'complete'
            or not applied.get('strict')
            or applied.get('manifest_sha256') != _sha256(manifest_path)):
        raise ValueError('Applied ANL handoff does not match its strict manifest.')
    records = {}
    for group, required_files in (
            ('species', ('composite', 'formation')),
            ('transition_states', ('composite',))):
        expected = set((manifest.get(group) or {}).keys())
        observed = set((applied.get(group) or {}).keys())
        if expected != observed:
            raise ValueError(
                f'Applied ANL handoff has different {group} coverage.')
        group_records = {}
        for ident, record in sorted((applied.get(group) or {}).items()):
            if not isinstance(record, dict):
                raise ValueError(
                    f'{ident}: applied ANL handoff record is invalid.')
            files = {}
            for key in required_files:
                path = Path(record.get(key, '')).expanduser().resolve()
                expected_hash = record.get(f'{key}_sha256')
                if not path.is_file() or expected_hash != _sha256(path):
                    raise ValueError(
                        f'{ident}: applied {key} file failed hash verification.')
                files[key] = {
                    'path': str(path), 'sha256': expected_hash,
                }
            if record.get('barrier_composite') is not None:
                path = Path(record['barrier_composite']).expanduser().resolve()
                expected_hash = record.get('barrier_composite_sha256')
                if not path.is_file() or expected_hash != _sha256(path):
                    raise ValueError(
                        f'{ident}: applied barrier composite failed hash '
                        'verification.')
                files['barrier_composite'] = {
                    'path': str(path), 'sha256': expected_hash}
            group_records[ident] = {
                'method': record.get('method'),
                'zero_k_hartree': record.get('zero_k_hartree'),
                'formation_0k_kj_mol': record.get(
                    'formation_0k_kj_mol'),
                'files': files,
                'vpt2': record.get('vpt2'),
            }
        records[group] = group_records
    return {
        'schema': 1, 'status': 'passed', 'strict': True,
        'manifest': str(manifest_path),
        'manifest_sha256': applied['manifest_sha256'],
        'applied_record': str(applied_path), **records,
    }


def _chemistry_report(base, data):
    """Expose calculated ANL components and CBH values in the final report."""
    report = {}
    for name, entry in sorted(data['cbh'].items()):
        path = _path(base, entry.get('output'), f'cbh.{name}.output')
        record = _read_json(path, f'cbh.{name} output')
        if record.get('schema') != 1 or record.get('status') != 'complete':
            raise ValueError(f'cbh.{name}: result is incomplete.')
        report[name] = {
            'path': str(path), 'sha256': _sha256(path),
            'reaction': record.get('reaction'),
            'formation': record.get('formation'),
            'component_breakdown': record.get('component_breakdown'),
        }
    return report


def _mess_rotd_statistical_factors(path, staged_name):
    """Read the statistical factors for one rendered scalar ROTD barrier."""
    path = Path(path)
    text = path.read_text(errors='replace')
    number = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?'
    top = r'(?=^[ \t]*(?:Well|Bimolecular|Barrier)\s|\Z)'
    barriers = re.findall(r'(?ms)^[ \t]*Barrier\s+.*?' + top, text)
    barriers = [block for block in barriers
                if re.search(rf'(?m)^\s*File\s+{re.escape(staged_name)}\s*$',
                             block)]
    if len(barriers) != 1:
        raise ValueError(
            'Expected one MESS ROTD barrier using the staged surface.')
    barrier = barriers[0]
    endpoints = barrier.splitlines()[0].split()[2:]
    core = re.search(
        rf'(?ms)Core\s+Rotd\b.*?SymmetryFactor\s+({number}).*?'
        r'End(?:\s*!\s*Core)?', barrier)
    barrier_electronic = re.search(
        rf'(?ms)ElectronicLevels\[1/cm\]\s+\d+\s*\n\s*'
        rf'{number}\s+(\d+)', barrier)
    if core is None or barrier_electronic is None:
        raise ValueError('MESS ROTD barrier statistical factors are incomplete.')

    bimolecular = re.findall(
        r'(?ms)^[ \t]*Bimolecular\s+.*?' + top, text)
    candidates = []
    for block in bimolecular:
        if block.splitlines()[0].split()[1] not in endpoints:
            continue
        fragments = re.findall(
            r'(?ms)^[ \t]*Fragment\s+.*?'
            r'(?=^[ \t]*(?:Fragment\s|GroundEnergy\b))', block)
        if len(fragments) == 2:
            candidates.append(fragments)
    if len(candidates) != 1:
        raise ValueError(
            'Expected one two-fragment MESS model connected to the selected ROTD barrier.')
    symmetries = []
    multiplicities = []
    for fragment in candidates[0]:
        symmetry = re.search(
            rf'(?ms)Core\s+RigidRotor\b.*?SymmetryFactor\s+({number})',
            fragment)
        fragment_electronic = re.search(
            rf'(?ms)ElectronicLevels\[1/cm\]\s+\d+\s*\n\s*'
            rf'{number}\s+(\d+)', fragment)
        if symmetry is None or fragment_electronic is None:
            raise ValueError(
                'MESS bimolecular fragment statistical factors are incomplete.')
        symmetries.append(float(symmetry.group(1)))
        multiplicities.append(int(fragment_electronic.group(1)))
    return {
        'mess_input': str(path.resolve()),
        'mess_input_sha256': _sha256(path),
        'fragment_rotational_symmetries': symmetries,
        'fragment_multiplicities': multiplicities,
        'reactive_surface_multiplicity': int(barrier_electronic.group(1)),
        'mess_symmetry_factor': float(core.group(1)),
    }


def audit(source):
    manifest_path, data = _manifest(source)
    base = manifest_path.parent
    run_dir = _path(base, data['run_dir'], 'run_dir')
    rotdpy_input = _path(base, data['rotdpy'].get('input'), 'rotdpy.input')
    rotdpy_result = read_result(rotdpy_input, required_profile='production')
    selected, selected_provenance = number_of_states_file(
        rotdpy_input, -1, required_profile='production')
    messpf = audit_messpf_run(run_dir)
    for ident in data['mess'].get('required_hindered_rotor_species', ()):
        try:
            count = messpf['species'][str(ident)]['hindered_rotor_blocks']
        except KeyError as error:
            raise ValueError(
                f'Required MESSPF rotor species {ident!r} is absent.') from error
        if count < 1:
            raise ValueError(
                f'{ident}: MESSPF has no accepted hindered-rotor block.')
    handoff = _audit_handoff(base, data, run_dir)
    chemistry = _chemistry_report(base, data)
    kinetics = audit_rate_tables(run_dir, job=data['mess'].get('job', 'mess_0000'))
    sources_path = run_dir / 'me' / 'rotdpy_sources.json'
    sources = _read_json(sources_path, 'MESS ROTD_py provenance')
    reaction = data['rotdpy'].get('reaction')
    try:
        consumed = sources['reactions'][reaction]
    except (KeyError, TypeError) as error:
        raise ValueError('MESS did not record the requested ROTD_py reaction.') from error
    if (consumed.get('source_sha256') != selected_provenance['source_sha256']
            or Path(consumed.get('source', '')).resolve() != selected.resolve()):
        raise ValueError('MESS did not consume the selected production ROTD_py surface.')
    staged_name = consumed.get('mess_file')
    staged = run_dir / 'me' / str(staged_name)
    preprocessing = consumed.get('mess_preprocessing')
    boundary = preprocessing.get('row') if isinstance(preprocessing, dict) else None
    if (not isinstance(staged_name, str)
            or Path(staged_name).name != staged_name
            or not staged.is_file()
            or consumed.get('mess_file_sha256') != _sha256(staged)
            or not isinstance(preprocessing, dict)
            or preprocessing.get('name') != 'explicit_zero_flux_origin'
            or not isinstance(boundary, list) or len(boundary) != 2
            or any(not isinstance(value, (int, float)) or not math.isfinite(value) for value in boundary)
            or boundary[1] != 0.):
        raise ValueError(
            'MESS-facing ROTD_py surface provenance is incomplete.')
    first_row = next((line.split() for line in staged.read_text().splitlines()
                      if line.split()), None)
    try:
        explicit_origin = ([float(value) for value in first_row] == boundary)
    except (TypeError, ValueError):
        explicit_origin = False
    if not explicit_origin:
        raise ValueError(
            'MESS-facing ROTD_py surface differs from its declared zero-flux boundary.')
    job = data['mess'].get('job', 'mess_0000')
    mess_input = run_dir / 'me' / f'{job}.inp'
    rendered_factors = _mess_rotd_statistical_factors(
        mess_input, staged.name)
    selected_correction = selected_provenance.get('dynamical_correction')
    if selected_correction is None:
        raise ValueError(
            'Selected production ROTD_py surface lacks dynamical correction '
            'provenance.')
    statistical_factors = {
        **rendered_factors,
        'dynamical_correction': selected_correction,
        'dynamical_correction_application': selected_provenance.get(
            'dynamical_correction_application'),
        'rotational_symmetry_product': math.prod(
            rendered_factors['fragment_rotational_symmetries']),
        'indistinguishability_factor': consumed.get(
            'indistinguishability_factor'),
        'uncertainty_factor': consumed.get('uncertainty_factor'),
    }
    statistical_factors['raw_flux_equivalent_mess_symmetry_factor'] = (
        statistical_factors['mess_symmetry_factor'] / selected_correction)
    statistical_factors['electronic_statistical_ratio'] = (
        statistical_factors['reactive_surface_multiplicity'] /
        math.prod(statistical_factors['fragment_multiplicities']))
    for key in ('rotational_symmetry_product', 'mess_symmetry_factor'):
        if key in consumed and not math.isclose(
                float(consumed[key]), float(statistical_factors[key]),
                rel_tol=1.e-12, abs_tol=1.e-12):
            raise ValueError(
                f'MESS ROTD_py provenance disagrees with rendered {key}.')
    for key in ('fragment_multiplicities',
                'reactive_surface_multiplicity'):
        if key in consumed and consumed[key] != statistical_factors[key]:
            raise ValueError(
                f'MESS ROTD_py provenance disagrees with rendered {key}.')
    statistical_expectations = {
        'dynamical_correction': 'dynamical_correction',
        'expected_fragment_rotational_symmetries': (
            'fragment_rotational_symmetries'),
        'expected_fragment_multiplicities': 'fragment_multiplicities',
        'expected_reactive_surface_multiplicity': (
            'reactive_surface_multiplicity'),
        'expected_electronic_statistical_ratio': (
            'electronic_statistical_ratio'),
        'expected_rotational_symmetry_product': (
            'rotational_symmetry_product'),
        'expected_indistinguishability_factor': (
            'indistinguishability_factor'),
        'expected_uncertainty_factor': 'uncertainty_factor',
        'expected_mess_symmetry_factor': 'mess_symmetry_factor',
        'expected_raw_flux_equivalent_mess_symmetry_factor': (
            'raw_flux_equivalent_mess_symmetry_factor'),
    }
    for manifest_key, provenance_key in statistical_expectations.items():
        if manifest_key not in data['mess']:
            continue
        expected = data['mess'][manifest_key]
        observed = statistical_factors.get(provenance_key)
        if isinstance(expected, (int, float)) and not isinstance(expected, bool):
            try:
                matches = math.isclose(
                    float(observed), float(expected), rel_tol=1.e-12,
                    abs_tol=1.e-12)
            except (TypeError, ValueError):
                matches = False
        else:
            matches = observed == expected
        if not matches:
            raise ValueError(
                f'MESS ROTD_py {provenance_key} is {observed!r}; expected '
                f'{expected!r}.')
    literature = _literature_checks(base, data.get('literature'))
    statuses = [kinetics['status'], *(item['status'] for item in literature)]
    status = 'passed' if all(value == 'passed' for value in statuses) else 'failed'
    payload = {
        'schema': 1, 'status': status,
        'manifest': str(manifest_path), 'run_dir': str(run_dir),
        'rotdpy': {
            'status': rotdpy_result['status'],
            'surface_count': rotdpy_result['surface_count'],
            'selected_number_of_states': str(selected),
            'selected_sha256': selected_provenance['source_sha256'],
            'mess_number_of_states': str(staged),
            'mess_sha256': consumed['mess_file_sha256'],
            'mess_preprocessing': preprocessing,
            'mess_source_record': str(sources_path),
            'statistical_factors': statistical_factors,
        },
        'handoff': handoff, 'chemistry': chemistry,
        'messpf': messpf, 'kinetics': kinetics,
        'mess_model': {
            key: data['mess'].get(key) for key in (
                'collider', 'epsilon', 'epsilon_unit', 'sigma',
                'collision_source')
        },
        'literature': literature,
    }
    report = _path(base, data.get(
        'final_report', 'full_verification_report.json'), 'final_report')
    _write_json(report, payload)
    payload['report'] = str(report)
    return payload


def _command_available(value, label):
    try:
        words = shlex.split(str(value))
    except ValueError as error:
        raise ValueError(f'Invalid {label}.') from error
    if not words or shutil.which(words[0]) is None:
        raise FileNotFoundError(f'{label} executable is unavailable: {value!r}.')


def run(source):
    prepared = prepare(source)
    parameters = _read_json(prepared['kinbot_input'], 'generated KinBot input')
    _command_available(
        os.environ.get('KINBOT_MESS_COMMAND', '').strip()
        or parameters['mess_command'], 'mess_command')
    _command_available(
        os.environ.get('KINBOT_MESSPF_COMMAND', '').strip()
        or parameters['messpf_command'], 'messpf_command')
    run_dir = Path(prepared['run_dir'])
    result = subprocess.run(
        [sys.executable, '-m', 'kinbot.kb', Path(
            prepared['kinbot_input']).name], cwd=run_dir, check=False)
    if result.returncode:
        log = run_dir / 'kinbot.log'
        detail = ''
        if log.is_file():
            lines = log.read_text(errors='replace').splitlines()
            detail = '\nLatest KinBot log entries:\n' + '\n'.join(lines[-40:])
        raise RuntimeError(
            f'KinBot/MESS continuation exited with status '
            f'{result.returncode}.{detail}')
    return audit(source)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Run the hash-verified CBH/ANL, MESSPF, and MESS gate')
    commands = parser.add_subparsers(dest='action', required=True)
    for action in ('prepare', 'run', 'audit'):
        command = commands.add_parser(action)
        command.add_argument('manifest', type=Path)
    args = parser.parse_args(argv)
    function = {'prepare': prepare, 'run': run, 'audit': audit}[args.action]
    payload = function(args.manifest)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return int(payload['status'] == 'failed')


if __name__ == '__main__':
    raise SystemExit(main())
