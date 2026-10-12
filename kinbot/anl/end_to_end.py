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
import hashlib
import json
import math
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from kinbot.anl.cbh import (HARTREE_TO_KJ_MOL, solve_composite_records,
                            zero_k_from_composite_record)
from kinbot.anl.kinetics import audit_methyl_recombination
from kinbot.anl.literature import compare_formation, compare_run
from kinbot.anl.thermochemistry import audit_messpf_run
from kinbot.anl.validation import (assemble_profiled_anl0_f12,
                                   molecule_from_completed_run)
from kinbot.rotdpy import number_of_states_file, read_result


HARTREE_TO_KCAL_MOL = HARTREE_TO_KJ_MOL / 4.184


def _read_json(path, label):
    try:
        value = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f'{label} is not valid JSON: {path}.') from error
    if not isinstance(value, dict):
        raise ValueError(f'{label} must contain a JSON object.')
    return value


def _write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    return path


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _path(base, value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{label} must be a nonempty path.')
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


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
    """Return signed component contributions for one generated CBH reaction."""
    records = {smiles: _read_json(path, f'{smiles} composite')
               for smiles, path in energy_files.items()}

    def higher_order(record):
        if 'higher_order_correction' in record.get('components', {}):
            return _composite_value(record, 'higher_order_correction')
        return math.fsum((_composite_value(record, 'hoe_high'),
                          -_composite_value(record, 'hoe_low')))

    categories = {
        'reference_cbs': lambda record: _composite_value(
            record, 'reference_cbs'),
        'higher_order': higher_order,
        'core_valence_cbs': lambda record: _composite_value(
            record, 'core_valence_cbs'),
        'scalar_relativistic': lambda record: _composite_value(
            record, 'scalar_relativistic'),
        'dboc': lambda record: _composite_value(record, 'dboc'),
        'spin_orbit': lambda record: _composite_value(record, 'spin_orbit'),
        'harmonic_zpe': lambda record: _composite_value(
            record, 'harmonic_zpe'),
        'anharmonic_zpe_correction': lambda record: _composite_value(
            record, 'vpt2_correction'),
    }
    result = {}
    for label, getter in categories.items():
        hartree = math.fsum(
            coefficient * getter(records[smiles])
            for smiles, coefficient in reaction.stoichiometry.items())
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
            energy_files = {
                smiles: composites[record]['path']
                for smiles, record in mapping.items()}
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
            'component_breakdown': component_breakdown(
                reaction, energy_files),
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
        vpt2 = _vpt2_entry(base, entry.get('vpt2'),
                           f'handoff.species.{routing}.vpt2')
        if vpt2 is not None:
            record['vpt2'] = vpt2
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
        transition_states[name] = record
    output = _path(base, specification.get('output'), 'handoff.output')
    payload = {
        'schema': 1, 'mode': 'cbh-anl',
        'strict': specification.get('strict', True),
        'species': species, 'transition_states': transition_states,
    }
    _write_json(output, payload)
    return output, payload


def _write_mess_input(base, manifest, handoff, rotdpy_input):
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
    run_dir, kinbot_input = _write_mess_input(
        base, data, handoff, rotdpy_input)
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
        kind = entry.get('kind')
        if kind == 'run':
            result = compare_run(
                entry.get('benchmark'),
                _path(base, entry.get('source'), f'literature[{index}].source'))
        elif kind == 'formation':
            source = _path(
                base, entry.get('source'), f'literature[{index}].source')
            method = entry.get('method', 'auto')
            if method == 'auto':
                method = _read_json(source, 'CBH formation result')[
                    'formation']['method']
            result = compare_formation(
                entry.get('benchmark'), source, method=method,
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
    kinetics = audit_methyl_recombination(
        run_dir, job=data['mess'].get('job', 'mess_0000'),
        max_factor=float(data['mess'].get('maximum_rate_factor', 3.)))
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
            'mess_source_record': str(sources_path),
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
    _command_available(parameters['mess_command'], 'mess_command')
    _command_available(parameters['messpf_command'], 'messpf_command')
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
