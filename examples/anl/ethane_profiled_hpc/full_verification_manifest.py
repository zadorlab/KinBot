#!/usr/bin/env python3
"""Write the completed ethane/methyl CBH/ANL + MESS verification manifest."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


REACTION = '301020900180000000001_hom_sci_1_2'


def _require_workflow(path, label, allowed):
    try:
        workflow = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise SystemExit(f'{label} is not a readable workflow: {path}') from error
    name = workflow.get('name') if isinstance(workflow, dict) else None
    if (not isinstance(workflow, dict) or workflow.get('schema') != 1
            or name not in allowed):
        expected = ', '.join(sorted(allowed))
        raise SystemExit(
            f'{label} has workflow type {name!r}; expected one of: {expected}.')


def _review(path):
    return str(path) if path.is_file() else None


def assembly(interface, higher, corrections, state, spin_source, *,
             base=None, review=None):
    value = {
        'interface_run': str(interface),
        'higher_order_run': str(higher),
        'corrections_run': str(corrections),
        'state_id': state,
        'spin_orbit_hartree': 0.0,
        'spin_orbit_backend': 'known_zero',
        'spin_orbit_source': spin_source,
    }
    if base is not None:
        value['base_run'] = str(base)
    if review is not None:
        value['vpt2_review'] = str(review)
    return value


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('base', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--ethane-base', type=Path)
    parser.add_argument('--ethane-post', type=Path)
    parser.add_argument('--partition', default='day-long-cpu')
    parser.add_argument('--ppn', type=int, default=8)
    parser.add_argument('--max-jobs', type=int, default=8)
    args = parser.parse_args(argv)
    base = args.base.expanduser().resolve()
    refs = base / 'cbh0_reference_graphs'
    run = base / 'vrc_production_run'
    ethane_interface = base / 'anl_interface'
    ethane_base = (args.ethane_base.expanduser().resolve()
                   if args.ethane_base else
                   base / 'anl_current_base_ethane_449cfbd')
    ethane_post = (args.ethane_post.expanduser().resolve()
                   if args.ethane_post else
                   base / 'anl_post_geometry_ethane_449cfbd')
    methane_interface, methane_post = (
        refs / 'methane_interface', refs / 'methane_post')
    hydrogen_interface, hydrogen_post = (
        refs / 'hydrogen_interface', refs / 'hydrogen_post')
    methyl = refs / 'methyl_composite'
    atomic_h = refs / 'atomic_h_composite'
    required = {
        'production KinBot input': run / 'ethane_vrc_production.json',
        'production ROTD_py input': run / 'rotdPy' / f'{REACTION}.py',
        'ethane interface': ethane_interface / 'workflow.json',
        'ethane current base': ethane_base / 'workflow.json',
        'ethane post geometry': ethane_post / 'workflow.json',
        'methane interface': methane_interface / 'workflow.json',
        'methane post geometry': methane_post / 'workflow.json',
        'hydrogen interface': hydrogen_interface / 'workflow.json',
        'hydrogen post geometry': hydrogen_post / 'workflow.json',
        'methyl composite': methyl / 'workflow.json',
        'atomic H composite': atomic_h / 'workflow.json',
        'methane source-geometry comparison':
            refs / 'methane_source_qz' / 'workflow.json',
        'methyl source-geometry comparison':
            refs / 'methyl_source_qz' / 'workflow.json',
    }
    missing = [f'{label}: {path}' for label, path in required.items()
               if not path.is_file()]
    if missing:
        raise SystemExit('Missing completed-run prerequisite(s):\n  '
                         + '\n  '.join(missing))
    interface_types = {
        'anl1-f12-non-mrcc-interface-validation',
        'anl-composite-validation',
    }
    base_types = {'anl-current-base-validation',
                  'anl-composite-validation'}
    post_types = {'anl-post-geometry-validation',
                  'anl-composite-validation'}
    _require_workflow(ethane_interface / 'workflow.json',
                      'ethane interface', interface_types)
    _require_workflow(ethane_base / 'workflow.json',
                      'ethane current base', base_types)
    _require_workflow(ethane_post / 'workflow.json',
                      'ethane post geometry', post_types)
    for label, interface, post in (
            ('methane', methane_interface, methane_post),
            ('hydrogen', hydrogen_interface, hydrogen_post)):
        _require_workflow(interface / 'workflow.json',
                          f'{label} interface', interface_types)
        _require_workflow(post / 'workflow.json',
                          f'{label} post geometry', post_types)
    for label, composite in (('methyl', methyl), ('atomic H', atomic_h)):
        _require_workflow(composite / 'workflow.json', f'{label} composite',
                          {'anl-composite-validation'})
    output = (args.output.expanduser().resolve() if args.output else
              base / 'full_verification_manifest.json')
    literature = [
        {'kind': 'run', 'benchmark': 'methane-qz-2017',
         'source': str(refs / 'methane_source_qz'),
         'required_checks': [
             'ccsdt_dz', 'ccsdtq_dz', 'ccsdtqp_dz',
             'delta_q_dz', 'delta_p_dz']},
        {'kind': 'run', 'benchmark': 'methyl-qz-2017',
         'source': str(refs / 'methyl_source_qz'),
         'required_checks': [
             'ccsdt_dz', 'ccsdtq_dz', 'ccsdtqp_dz',
             'delta_q_dz', 'delta_p_dz']},
        # This same-geometry difference remains meaningful even when the
        # production L3 geometry is not byte-for-byte the source geometry.
        {'kind': 'run', 'benchmark': 'ethane-tz-2017',
         'source': str(ethane_post),
         'required_checks': ['delta_q_dz']},
        {'kind': 'formation', 'benchmark': 'ethane-tz-2017',
         'source': str(base / 'ethane_cbh0_anl0_f12.json'),
         'method': 'auto', 'absolute_tolerance_kcal_mol': 0.1},
        {'kind': 'formation', 'benchmark': 'methyl-qz-2017',
         'source': str(refs / 'methyl_cbh0_anl0_f12.json'),
         'method': 'auto', 'absolute_tolerance_kcal_mol': 0.1},
    ]
    manifest = {
        'schema': 1,
        'run_dir': str(run),
        'kinbot_input': str(run / 'ethane_vrc_production.json'),
        'preparation_report': str(base / 'full_verification_preparation.json'),
        'final_report': str(base / 'full_verification_report.json'),
        'rotdpy': {
            'input': str(run / 'rotdPy' / f'{REACTION}.py'),
            'reaction': REACTION,
        },
        'atct': {
            'version': '1.220',
            'cache': str(Path.home() / '.cache/kinbot/atct'),
            'refresh': False,
        },
        'composites': {
            'ethane': {
                'smiles': 'CC',
                'path': str(base / 'ethane_profiled_anl0_f12.json'),
                'assembly': assembly(
                    ethane_interface, ethane_post, ethane_post,
                    'ethane-X1Ag',
                    'zero for nondegenerate closed-shell ethane ground state',
                    base=ethane_base,
                    review=_review(base / 'ethane_vpt2_review.json')),
            },
            'methane': {
                'smiles': 'C',
                'path': str(refs / 'methane_profiled_anl0_f12.json'),
                'assembly': assembly(
                    methane_interface, methane_post, methane_post,
                    'methane-X1A1',
                    'zero for nondegenerate closed-shell methane ground state',
                    review=_review(refs / 'methane_vpt2_review.json')),
            },
            'hydrogen': {
                'smiles': '[H][H]',
                'path': str(refs / 'hydrogen_profiled_anl0_f12.json'),
                'assembly': assembly(
                    hydrogen_interface, hydrogen_post, hydrogen_post,
                    'hydrogen-X1Sigma-g-plus',
                    'zero for the H2 X1Sigma-g-plus ground state',
                    review=_review(refs / 'hydrogen_vpt2_review.json')),
            },
            'methyl': {
                'smiles': '[CH3]',
                'path': str(refs / 'methyl_profiled_anl0_f12.json'),
                'assembly': assembly(
                    methyl, methyl, methyl, 'methyl-X2A2-double-prime',
                    'zero first-order spin-orbit term for orbital-singlet CH3',
                    review=_review(refs / 'methyl_vpt2_review.json')),
            },
            'atomic_h': {
                'smiles': '[H]',
                'path': str(refs / 'atomic_h_profiled_anl0_f12.json'),
                'assembly': assembly(
                    atomic_h, atomic_h, atomic_h, 'hydrogen-2S-one-half',
                    'zero spin-orbit splitting for atomic H 2S ground state'),
            },
        },
        'cbh': {
            'ethane': {
                'target_smiles': 'CC', 'rung': 0,
                'energies': {'CC': 'ethane', 'C': 'methane',
                             '[H][H]': 'hydrogen'},
                'output': str(base / 'ethane_cbh0_anl0_f12.json'),
            },
            'methyl': {
                'target_smiles': '[CH3]', 'rung': 0,
                'energies': {'[CH3]': 'methyl', 'C': 'methane',
                             '[H][H]': 'hydrogen', '[H]': 'atomic_h'},
                'reference_ids': {'[H]': '12385-13-6*0'},
                'output': str(refs / 'methyl_cbh0_anl0_f12.json'),
            },
        },
        'handoff': {
            'output': str(base / 'anl_handoff.json'), 'strict': True,
            'species': {
                '301020900180000000001': {
                    'composite': 'ethane', 'formation': 'ethane',
                    'require_l2_geometry_match': True,
                    'vpt2': {
                        'run': str(ethane_interface),
                        **({'review': str(base / 'ethane_vpt2_review.json')}
                           if (base / 'ethane_vpt2_review.json').is_file()
                           else {}),
                    },
                },
                '150390060000000000002': {
                    'composite': 'methyl', 'formation': 'methyl',
                    'require_l2_geometry_match': True,
                    'vpt2': {
                        'run': str(methyl),
                        **({'review': str(refs / 'methyl_vpt2_review.json')}
                           if (refs / 'methyl_vpt2_review.json').is_file()
                           else {}),
                    },
                },
            },
            'transition_states': {},
        },
        'mess': {
            'input': str(run / 'full_cbh_anl_mess.json'),
            'mess_command': os.environ.get('KINBOT_MESS_COMMAND', 'mess'),
            'messpf_command': os.environ.get(
                'KINBOT_MESSPF_COMMAND', 'messpf'),
            'queue_name': args.partition, 'ppn': args.ppn,
            'queue_job_limit': args.max_jobs,
            # Standard CHEMKIN/GRI-Mech C2H6 transport parameters. The
            # literature gate below compares the high-pressure limit, while
            # these declared values make the pressure-dependent MESS model
            # complete and reproducible.
            'collider': 'He', 'epsilon': 252.3,
            'epsilon_unit': 'K', 'sigma': 4.302,
            'collision_source': {
                'name': 'GRI-Mech/CHEMKIN transport data for C2H6',
                'url': ('https://www.cantera.org/dev/examples/input/'
                        'gri30.html'),
            },
            'uq_max_runs': 1, 'job': 'mess_0000',
            'required_hindered_rotor_species': [
                '301020900180000000001'],
            'maximum_rate_factor': 3.0,
        },
        'literature': literature,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    print(output)


if __name__ == '__main__':
    main()
