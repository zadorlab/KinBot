"""Build and audit portable ANL interface validation graphs.

This module deliberately calls the result an *interface validation*.  It runs
the available Gaussian, Molpro, CFOUR, and optional MRCC calculation types and
verifies their native outputs without labeling a partial graph as a complete
ANL1 or ANL1-F12 composite energy.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path

from ase.db import connect
from ase.io import read

from kinbot.anl.dispatch import (_load, _verify_execution, _verify_stage_files,
                                 prepare)
from kinbot.anl.extrapolation import two_point_cbs
from kinbot.anl.model import ComponentResult
from kinbot.anl.recipes import recipe
from kinbot.anl.results import legacy_molpro_parser
from kinbot.anl.tasks import higher_order_task, molpro_task
from kinbot.anl.workflow import (_verified_task_result,
                                 cbs_task_component,
                                 core_valence_task_component,
                                 scalar_relativistic_task_component,
                                 state_correction_component,
                                 task_component)


MOLPRO_HEADER = """***,KinBot ANL interface validation
symmetry,nosym
orient,noorient
geomtyp=xyz
geometry={
{{XYZ}}
}
set,charge={{CHARGE}}
set,spin={{SPIN}}
"""


def _resources(walltime, *, max_cores=8, min_stack_mw=None, partition=None):
    result = {'cores': 'auto', 'memory_mb': 'node', 'walltime': walltime,
              'max_cores': max_cores}
    if min_stack_mw is not None:
        result['min_stack_mw'] = min_stack_mw
    if partition is not None:
        result['partition'] = partition
    return result


def _molpro_task(ident, body, *, geometry_from='l3_geometry',
                  walltime='04:00:00', max_cores=8, parser=None,
                  partition=None):
    task = {
        'id': ident, 'kind': 'external', 'backend': 'molpro',
        'geometry_from': geometry_from,
        'resources': _resources(walltime, max_cores=max_cores,
                                min_stack_mw=1024, partition=partition),
        'input_name': f'{ident}.inp',
        'input_template': MOLPRO_HEADER + body,
        'command': ['molpro', '-n', '{cores}', '-m',
                    '{molpro_stack_mw}', '{input}'],
        'stdout': 'launcher.stdout', 'stderr': 'launcher.stderr',
        'required_outputs': [f'{ident}.out'],
        'success_marker': {'file': f'{ident}.out',
                           'contains': 'Molpro calculation terminated'},
    }
    if parser:
        task['result_parser'] = {'file': f'{ident}.out', **parser}
    return task


def interface_validation_spec(molecule, *, max_nodes=3, partition=None):
    """Return a molecule-independent Gaussian/Molpro/CFOUR task graph."""
    if not isinstance(molecule, dict):
        raise TypeError('molecule must be an object.')
    l2_profile = {
        'calculator': 'gaussian', 'method': 'B2PLYP',
        'basis': 'cc-pVTZ', 'command': 'g16',
        'calculator_kwargs': {
            'EmpiricalDispersion': 'GD3BJ', 'Symm': 'None',
            'scf': 'xqc', 'integral': 'UltraFine'},
        'optimizer': 'sella',
    }
    reference = 'RHF' if molecule.get('multiplicity', 1) == 1 else 'ROHF'
    conventional = 'uccsd(t),uhf_uccsd=1'
    f12 = 'uccsd(t)-f12b'
    tasks = [
        {
            'id': 'l2_geometry', 'kind': 'ase_optimize',
            'geometry_from': 'initial', 'geometry_output': 'final.xyz',
            'resources': _resources('06:00:00', max_cores=8,
                                    partition=partition),
            'profile': l2_profile,
            'optimizer': {'fmax': 0.0005, 'steps': 160,
                          'sella_kwargs': {'internal': True}},
        },
        {
            'id': 'l3_geometry', 'kind': 'ase_optimize',
            'geometry_from': 'l2_geometry', 'geometry_output': 'final.xyz',
            'resources': _resources('24:00:00', max_cores=12,
                                    min_stack_mw=1024, partition=partition),
            'profile': {
                'calculator': 'molpro', 'method': 'CCSD(T)',
                'basis': 'cc-pVTZ', 'command': 'molpro',
                'optimizer': 'sella'},
            'optimizer': {'fmax': 0.03, 'steps': 100,
                          'sella_kwargs': {'internal': True}},
        },
        _molpro_task(
            'harmonic',
            f'basis=cc-pVTZ\nrhf\n{conventional}\nfrequencies,numerical\n',
            walltime='24:00:00', max_cores=8,
            parser={'kind': 'molpro_harmonic', 'basis': 'cc-pVTZ',
                    'reference': reference},
            partition=partition),
        _molpro_task(
            'f12_tz',
            f'basis=cc-pVTZ-F12\nrhf\n{f12},scale_trip=1\n'
            'kb_f12b=energy\n',
            walltime='12:00:00', max_cores=12,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)-F12b',
                    'basis': 'cc-pVTZ-F12', 'reference': reference},
            partition=partition),
        _molpro_task(
            'f12_qz',
            f'basis=cc-pVQZ-F12\nrhf\n{f12},scale_trip=1\n'
            'kb_f12b=energy\n',
            walltime='24:00:00', max_cores=12,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)-F12b',
                    'basis': 'cc-pVQZ-F12', 'reference': reference},
            partition=partition),
        _molpro_task(
            'ccsdt_dz',
            f'basis=cc-pVDZ\nrhf\n{conventional}\nkb_dz_energy=energy\n',
            walltime='08:00:00', max_cores=8,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': 'cc-pVDZ', 'reference': reference},
            partition=partition),
        {
            'id': 'cfour_dboc', 'kind': 'external', 'backend': 'cfour',
            'geometry_from': 'l3_geometry',
            'resources': _resources('08:00:00', max_cores=4,
                                    partition=partition),
            'input_name': 'ZMAT',
            'input_template': (
                'KinBot DBOC interface validation\n{{CARTESIAN}}\n\n'
                '*CFOUR(CALC=SCF\nBASIS=cc-pVTZ\nDBOC=ON\n'
                'COORDINATES=CARTESIAN\nUNITS=ANGSTROM\nCHARGE={{CHARGE}}\n'
                'MULTIPLICITY={{MULT}}\nMEM_UNIT=MB\n'
                'MEMORY_SIZE={{WORK_MEMORY_MB}})\n'),
            'command': ['xcfour'], 'stdout': 'cfour.out',
            'stderr': 'cfour.err', 'required_outputs': ['cfour.out'],
            'files_from_env': {'GENBAS': 'CFOUR_GENBAS'},
            'success_marker': {
                'file': 'cfour.out',
                'contains': 'The total diagonal Born-Oppenheimer correction (DBOC) is:'},
            'result_parser': {'kind': 'cfour_dboc', 'file': 'cfour.out',
                              'level': 'HF', 'basis': 'cc-pVTZ'},
        },
        {
            'id': 'gaussian_vpt2', 'kind': 'external',
            'backend': 'gaussian', 'geometry_from': 'l2_geometry',
            'depends_on': ['l3_geometry'],
            'resources': _resources('24:00:00', max_cores=8,
                                    partition=partition),
            'input_name': 'vpt2.com',
            'input_template': (
                '%nprocshared={{CORES}}\n%mem={{WORK_MEMORY_MB}}MB\n'
                '#p B2PLYP/cc-pVTZ Freq=Anharmonic NoSymm SCF=XQC '
                'EmpiricalDispersion=GD3BJ Integral=UltraFine\n\n'
                'KinBot frequency-only VPT2 interface validation\n\n'
                '{{CHARGE}} {{MULT}}\n{{CARTESIAN}}\n\n'),
            'command': ['g16'], 'stdin': 'vpt2.com',
            'stdout': 'vpt2.log', 'stderr': 'vpt2.err',
            'required_outputs': ['vpt2.log'],
            'success_marker': {'file': 'vpt2.log',
                               'contains': 'Normal termination of Gaussian'},
            'result_parser': {'kind': 'gaussian_vpt2', 'file': 'vpt2.log',
                              'method': 'B2PLYP', 'basis': 'cc-pVTZ',
                              'dispersion': 'GD3BJ'},
        },
    ]
    return {
        'schema': 1, 'name': 'anl1-f12-non-mrcc-interface-validation',
        'molecule': molecule,
        'limits': {'max_nodes': max_nodes},
        'tasks': tasks,
        'intent': {
            'requested_ladder_head': 'ANL1-F12',
            'mrcc_enabled': False,
            'claim': 'interface-validation-only',
            'unavailable_components': [
                'ANL1-F12 pinned electronic equation',
                'CCSDTQ(P)/cc-pVDZ (MRCC)',
                'core-valence CBS provider',
                'scalar-relativistic provider',
                'state-specific spin-orbit provider',
                'complete recipe assembly and CBH/ATcT solve'],
        },
    }


def current_base_validation_spec(molecule, *, max_nodes=3, partition=None):
    """Build current unrestricted L3 geometry and base ANL0-F12 jobs.

    ``molecule`` is normally a verified L2 geometry imported from an older
    completed interface run.  This focused continuation avoids repeating the
    L2 optimization and Gaussian VPT2 calculation while ensuring no legacy
    restricted Molpro result enters the current recipe.
    """
    full = interface_validation_spec(
        molecule, max_nodes=max_nodes, partition=partition)
    keep = {'l3_geometry', 'harmonic', 'f12_tz', 'f12_qz', 'cfour_dboc'}
    tasks = [deepcopy(task) for task in full['tasks'] if task['id'] in keep]
    for task in tasks:
        if task['id'] == 'l3_geometry':
            task['geometry_from'] = 'initial'
    return {
        'schema': 1, 'name': 'anl-current-base-validation',
        'molecule': molecule, 'limits': {'max_nodes': max_nodes},
        'tasks': tasks,
        'intent': {
            'claim': 'current-unrestricted-base-interface-validation',
            'source_geometry_role': 'l2',
            'provides': ['current unrestricted L3 geometry',
                         'CCSD(T)-F12b TZ/QZ reference pair',
                         'CCSD(T)/cc-pVTZ harmonic ZPE',
                         'HF/cc-pVTZ DBOC'],
        },
    }


_HIGHER_ORDER_TASK_IDS = (
    'ccsdt_tz', 'ccsdt_dz', 'ccsdtq_tz', 'ccsdtq_dz', 'ccsdtqp_dz')


def higher_order_validation_spec(molecule, *, max_nodes=3, partition=None,
                                 mrcc_command='dmrcc', task_ids=None):
    """Build the five-job ANL1 higher-order interface probe.

    This graph checks native input generation, execution, parsing, and the
    cross-program correction arithmetic on one supplied geometry.  It is a
    focused backend test and does not claim a complete composite energy.
    """
    if not isinstance(molecule, dict):
        raise TypeError('molecule must be an object.')
    multiplicity = molecule.get('multiplicity', 1)
    restricted_reference = 'RHF' if multiplicity == 1 else 'ROHF'
    conventional = 'uccsd(t),uhf_uccsd=1'
    common = dict(geometry_from='initial', partition=partition)
    tasks = [
        molpro_task(
            'ccsdt_tz',
            f'basis=cc-pVTZ\nrhf\n{conventional}\nkb_ccsdt=energy\n',
            walltime='08:00:00', max_cores=8,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': 'cc-pVTZ',
                    'reference': restricted_reference}, **common),
        molpro_task(
            'ccsdt_dz',
            f'basis=cc-pVDZ\nrhf\n{conventional}\nkb_ccsdt=energy\n',
            walltime='04:00:00', max_cores=8,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': 'cc-pVDZ',
                    'reference': restricted_reference}, **common),
        higher_order_task(
            'ccsdtq_tz', 'CCSDT(Q)', 'cc-pVTZ',
            multiplicity=multiplicity, walltime='24:00:00',
            max_cores=8, command=mrcc_command, **common),
        higher_order_task(
            'ccsdtq_dz', 'CCSDT(Q)', 'cc-pVDZ',
            multiplicity=multiplicity, walltime='12:00:00',
            max_cores=8, command=mrcc_command, **common),
        higher_order_task(
            'ccsdtqp_dz', 'CCSDTQ(P)', 'cc-pVDZ',
            multiplicity=multiplicity, walltime='7-00:00:00',
            max_cores=8, command=mrcc_command, **common),
    ]
    if task_ids is not None:
        requested = tuple(dict.fromkeys(task_ids))
        unknown = sorted(set(requested) - set(_HIGHER_ORDER_TASK_IDS))
        if not requested or unknown:
            detail = f': {unknown}' if unknown else ''
            raise ValueError(f'Invalid higher-order task selection{detail}.')
        tasks = [task for task in tasks if task['id'] in requested]
    else:
        requested = _HIGHER_ORDER_TASK_IDS
    return {
        'schema': 1, 'name': 'anl1-higher-order-interface-validation',
        'molecule': molecule, 'limits': {'max_nodes': max_nodes},
        'tasks': tasks,
        'intent': {
            'claim': 'higher-order-interface-validation-only',
            'equation': ('CCSDT(Q)/TZ - CCSD(T)/TZ + CCSDTQ(P)/DZ '
                         '- CCSDT(Q)/DZ'),
            'backend_policy': {
                'closed_shell_ccsdt_q': 'cfour-uhf-vcc-unrestricted-cc',
                'open_shell_ccsdt_q': 'direct-mrcc-semicanonical-rohf',
                'closed_shell_ccsdtq_p': 'direct-mrcc-rhf-unrestricted-cc',
                'open_shell_ccsdtq_p': 'direct-mrcc-semicanonical-rohf',
            },
            'selected_tasks': list(requested),
        },
    }


def common_corrections_validation_spec(molecule, *, max_nodes=3,
                                       partition=None):
    """Build core-valence and scalar-relativistic ANL correction jobs."""
    if not isinstance(molecule, dict):
        raise TypeError('molecule must be an object.')
    reference = 'RHF' if molecule.get('multiplicity', 1) == 1 else 'ROHF'
    common = dict(geometry_from='initial', partition=partition)

    def energy_task(ident, basis, *, core, relativistic='none', walltime):
        dkh = 'set,dkho=2\n' if relativistic == 'DKH2' else ''
        command = ('{uccsd(t),uhf_uccsd=1;core}'
                   if core == 'all-electron'
                   else 'uccsd(t),uhf_uccsd=1')
        body = (f'basis={basis}\n{dkh}rhf\n{command}\n'
                f'kb_{ident}=energy\n')
        return molpro_task(
            ident, body, walltime=walltime, max_cores=8,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': basis, 'reference': reference, 'core': core,
                    'relativistic': relativistic}, **common)

    tasks = [
        energy_task('cv_ae_tz', 'cc-pCVTZ', core='all-electron',
                    walltime='12:00:00'),
        energy_task('cv_ae_qz', 'cc-pCVQZ', core='all-electron',
                    walltime='24:00:00'),
        energy_task('cv_fc_tz', 'cc-pCVTZ', core='frozen',
                    walltime='12:00:00'),
        energy_task('cv_fc_qz', 'cc-pCVQZ', core='frozen',
                    walltime='24:00:00'),
        energy_task('rel_dkh', 'aug-cc-pCVTZ-DK', core='all-electron',
                    relativistic='DKH2', walltime='12:00:00'),
        energy_task('rel_nonrel', 'aug-cc-pCVTZ-DK', core='all-electron',
                    relativistic='none', walltime='12:00:00'),
    ]
    return {
        'schema': 1, 'name': 'anl-common-corrections-validation',
        'molecule': molecule, 'limits': {'max_nodes': max_nodes},
        'tasks': tasks,
        'intent': {
            'claim': 'common-corrections-interface-validation-only',
            'core_valence_equation': (
                'CBS[CCSD(T,all-electron),TZ/QZ] - '
                'CBS[CCSD(T,frozen-core),TZ/QZ]'),
            'scalar_relativistic_equation': (
                'CCSD(T,all-electron,DKH2)/aug-cc-pCVTZ-DK - '
                'CCSD(T,all-electron,nonrel)/aug-cc-pCVTZ-DK'),
        },
    }


def post_geometry_validation_spec(molecule, *, max_nodes=3, partition=None,
                                  mrcc_command='dmrcc'):
    """Build one globally throttled fan-out after the accepted L3 geometry."""
    higher = higher_order_validation_spec(
        molecule, max_nodes=max_nodes, partition=partition,
        mrcc_command=mrcc_command)
    corrections = common_corrections_validation_spec(
        molecule, max_nodes=max_nodes, partition=partition)
    tasks = higher['tasks'] + corrections['tasks']
    identifiers = [task['id'] for task in tasks]
    if len(identifiers) != len(set(identifiers)):
        raise RuntimeError('Post-geometry validation task identifiers overlap.')
    return {
        'schema': 1, 'name': 'anl-post-geometry-validation',
        'molecule': molecule, 'limits': {'max_nodes': max_nodes},
        'tasks': tasks,
        'intent': {
            'claim': 'post-geometry-interface-validation-only',
            'higher_order_equation': higher['intent']['equation'],
            'backend_policy': higher['intent']['backend_policy'],
            'core_valence_equation':
                corrections['intent']['core_valence_equation'],
            'scalar_relativistic_equation':
                corrections['intent']['scalar_relativistic_equation'],
            'scheduling': ('one shared max_nodes limit for every independent '
                           'post-geometry task'),
        },
    }


def audit_higher_order_run(run_dir):
    """Reparse a complete full or targeted higher-order interface probe."""
    _, spec, state = _load(run_dir)
    if spec.get('name') not in ('anl1-higher-order-interface-validation',
                                 'anl-post-geometry-validation'):
        raise ValueError('Run is not an ANL1 higher-order validation graph.')
    declared = {task['id'] for task in spec['tasks']}
    identifiers = tuple(ident for ident in _HIGHER_ORDER_TASK_IDS
                        if ident in declared)
    if not identifiers:
        raise ValueError('Higher-order run declares no recognized tasks.')
    statuses = {ident: state['tasks'].get(ident, {}).get('status', 'waiting')
                for ident in identifiers}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'Higher-order run is incomplete: {statuses}')
    state_id = 'higher-order-validation-state'
    components = {
        ident: task_component(run_dir, ident, key=ident, state_id=state_id)
        for ident in identifiers
    }
    corrections = {}
    pairs = {
        'delta_q_dz': ('ccsdtq_dz', 'ccsdt_dz'),
        'delta_q_tz': ('ccsdtq_tz', 'ccsdt_tz'),
        'delta_p_dz': ('ccsdtqp_dz', 'ccsdtq_dz'),
    }
    for name, (high, low) in pairs.items():
        if high in components and low in components:
            corrections[name] = math.fsum((
                components[high].value_hartree,
                -components[low].value_hartree))
    if {'delta_q_tz', 'delta_p_dz'} <= corrections.keys():
        corrections['anl1_higher_order'] = math.fsum((
            corrections['delta_q_tz'], corrections['delta_p_dz']))
    result = {
        'status': ('higher_order_interface_complete'
                   if set(identifiers) == set(_HIGHER_ORDER_TASK_IDS)
                   else 'targeted_higher_order_interface_complete'),
        'task_statuses': statuses,
        'corrections_hartree': corrections,
        'components': {
            key: {
                'energy_hartree': value.value_hartree,
                'method': value.method, 'basis': value.basis,
                'backend': value.backend,
                'reference': value.settings.get('reference'),
                'correlation': value.settings.get('correlation'),
                'driver': value.settings.get('driver'),
                'program': value.settings.get('program'),
                'source_sha256': value.source_sha256,
            } for key, value in components.items()
        },
        'claim': 'interface-validation-only',
    }
    if set(identifiers) == set(_HIGHER_ORDER_TASK_IDS):
        result['correction_hartree'] = corrections['anl1_higher_order']
    return result


def audit_anl0_post_geometry_run(run_dir):
    """Audit only the terms required by an ANL0 or ANL0-F12 assembly.

    A combined validation graph may also contain ANL1-only CCSDT(Q)/TZ and
    CCSDTQ(P)/DZ probes.  Their failure or deliberate cancellation must not
    prevent an independently complete ANL0-F12 result from being audited.
    """
    run_dir, spec, state = _load(run_dir)
    if spec.get('name') != 'anl-post-geometry-validation':
        raise ValueError('Run is not a combined post-geometry validation graph.')
    required = ('ccsdt_dz', 'ccsdtq_dz', 'cv_ae_tz', 'cv_ae_qz',
                'cv_fc_tz', 'cv_fc_qz', 'rel_dkh', 'rel_nonrel')
    declared = {task['id'] for task in spec['tasks']}
    missing = sorted(set(required) - declared)
    if missing:
        raise RuntimeError(f'ANL0 post-geometry tasks are absent: {missing}')
    statuses = {ident: state['tasks'].get(ident, {}).get('status', 'waiting')
                for ident in required}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'ANL0 post-geometry run is incomplete: {statuses}')
    state_id = 'anl0-post-geometry-validation-state'
    high = task_component(
        run_dir, 'ccsdtq_dz', key='hoe_high', state_id=state_id)
    low = task_component(
        run_dir, 'ccsdt_dz', key='hoe_low', state_id=state_id)
    corrections = audit_common_corrections_run(run_dir)
    all_statuses = {
        task['id']: state['tasks'].get(task['id'], {}).get('status', 'waiting')
        for task in spec['tasks']}
    return {
        'status': 'anl0_post_geometry_interface_complete',
        'task_statuses': statuses,
        'other_task_statuses': {
            key: value for key, value in all_statuses.items()
            if key not in statuses},
        'higher_order_dz_hartree': math.fsum((
            high.value_hartree, -low.value_hartree)),
        'higher_order_components': {
            'ccsdtq_dz': asdict(high), 'ccsdt_dz': asdict(low)},
        'common_corrections': corrections,
        'geometry_sha256': spec.get('molecule', {}).get(
            'source', {}).get('geometry_sha256'),
    }


def audit_common_corrections_run(run_dir):
    """Reparse and combine completed core-valence and DKH correction jobs."""
    _, spec, state = _load(run_dir)
    if spec.get('name') not in ('anl-common-corrections-validation',
                                 'anl-post-geometry-validation'):
        raise ValueError('Run is not an ANL common-corrections graph.')
    identifiers = ('cv_ae_tz', 'cv_ae_qz', 'cv_fc_tz', 'cv_fc_qz',
                   'rel_dkh', 'rel_nonrel')
    statuses = {ident: state['tasks'].get(ident, {}).get('status', 'waiting')
                for ident in identifiers}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'Common-corrections run is incomplete: {statuses}')
    equation = recipe(
        'ANL0-F12', vpt2_method='B2PLYP-D3BJ',
        multiplicity=spec['molecule'].get('multiplicity', 1))
    requirements = {item.key: item for item in equation.requirements}
    state_id = 'common-corrections-validation-state'
    core_valence = core_valence_task_component(
        run_dir, all_electron_lower='cv_ae_tz',
        all_electron_upper='cv_ae_qz', frozen_core_lower='cv_fc_tz',
        frozen_core_upper='cv_fc_qz',
        requirement=requirements['core_valence_cbs'], state_id=state_id)
    relativistic = scalar_relativistic_task_component(
        run_dir, 'rel_dkh', 'rel_nonrel',
        requirement=requirements['scalar_relativistic'], state_id=state_id)
    return {
        'status': 'common_corrections_interface_complete',
        'task_statuses': statuses,
        'core_valence_hartree': core_valence.value_hartree,
        'scalar_relativistic_hartree': relativistic.value_hartree,
        'components': {
            'core_valence_cbs': {
                'source_sha256': core_valence.source_sha256,
                'method': core_valence.method, 'basis': core_valence.basis},
            'scalar_relativistic': {
                'source_sha256': relativistic.source_sha256,
                'method': relativistic.method, 'basis': relativistic.basis},
        },
        'claim': 'interface-validation-only',
    }


def audit_post_geometry_run(run_dir):
    """Audit the shared post-geometry fan-out and both correction groups."""
    run_dir, spec, state = _load(run_dir)
    if spec.get('name') != 'anl-post-geometry-validation':
        raise ValueError('Run is not a combined post-geometry validation graph.')
    statuses = {task['id']: state['tasks'].get(task['id'], {}).get(
        'status', 'waiting') for task in spec['tasks']}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'Post-geometry run is incomplete: {statuses}')
    source = spec.get('molecule', {}).get('source', {})
    higher = audit_higher_order_run(run_dir)
    corrections = audit_common_corrections_run(run_dir)
    return {
        'status': 'post_geometry_interface_complete',
        'task_statuses': statuses,
        'geometry_sha256': source.get('geometry_sha256'),
        'higher_order': higher,
        'common_corrections': corrections,
    }


def molecule_from_database(database, job, *, charge, multiplicity):
    rows = list(connect(str(database)).select(name=job))
    if not rows:
        raise ValueError(f'No KinBot database row named {job!r}.')
    row = rows[-1]
    if row.data.get('status') != 'normal':
        raise ValueError(f'{job}: latest KinBot record is not normal.')
    missing = [key for key in ('energy', 'frequencies', 'zpe')
               if row.data.get(key) is None]
    if missing:
        raise ValueError(f'{job}: incomplete accepted L2 record: {missing}.')
    return {'symbols': list(row.symbols), 'positions': row.positions.tolist(),
            'charge': charge, 'multiplicity': multiplicity}


def molecule_from_completed_run(run_dir, geometry_task='l3_geometry'):
    """Import one hash-verified optimized geometry from a dispatcher run."""
    run_dir, spec, state = _load(run_dir)
    tasks = {task['id']: task for task in spec['tasks']}
    task = tasks.get(geometry_task)
    entry = state['tasks'].get(geometry_task)
    if (task is None or task.get('kind') != 'ase_optimize'
            or entry is None or entry.get('status') != 'complete'):
        raise ValueError(f'{geometry_task}: completed optimization is unavailable.')
    _verify_stage_files(run_dir, task, entry)
    execution = json.loads(
        (run_dir / 'tasks' / geometry_task / 'execution.json').read_text())
    geometry_hash = _verify_execution(run_dir, task, entry, execution)
    if (execution.get('status') != 'executed' or geometry_hash is None
            or geometry_hash != entry.get('final_geometry_sha256')):
        raise ValueError(f'{geometry_task}: optimized geometry provenance failed.')
    path = run_dir / 'tasks' / geometry_task / task['geometry_output']
    atoms = read(path)
    profile = task.get('profile', {})
    source = {
        'run_dir': str(run_dir), 'task_id': geometry_task,
        'geometry_sha256': geometry_hash,
        'artifact_sha256': execution['artifacts'][task['geometry_output']],
    }
    if isinstance(profile, dict):
        source['profile'] = {
            key: profile[key] for key in ('calculator', 'method', 'basis')
            if key in profile}
    return {
        'symbols': atoms.get_chemical_symbols(),
        'positions': atoms.get_positions().tolist(),
        'charge': spec['molecule'].get('charge', 0),
        'multiplicity': spec['molecule'].get('multiplicity', 1),
        'source': source,
    }


def _same_imported_geometry(source_run, run_dir, geometry_hash, role,
                            source_task_id):
    """Require a derived run to identify one verified source geometry."""
    run_dir, spec, _ = _load(run_dir)
    source = spec.get('molecule', {}).get('source')
    if (not isinstance(source, dict)
            or source.get('geometry_sha256') != geometry_hash
            or source.get('task_id') != source_task_id):
        raise ValueError(
            f'{run_dir}: molecule is not imported from the accepted {role} geometry.')
    source_run = Path(source_run).resolve()
    recorded = source.get('run_dir')
    if (not isinstance(recorded, str)
            or Path(recorded).resolve() != source_run):
        raise ValueError(
            f'{run_dir}: molecule source does not identify the source run.')


def _same_imported_l3_geometry(interface_run, run_dir, l3_hash):
    """Backward-compatible L3 provenance helper."""
    _same_imported_geometry(
        interface_run, run_dir, l3_hash, 'L3', 'l3_geometry')


def _apply_quality_review(component, task_id, review_file):
    """Accept one flagged native component through a hash-bound review file."""
    if not component.review_required:
        if review_file is not None:
            raise ValueError(f'{task_id}: a review was supplied for an '
                             'unflagged component.')
        return component
    if review_file is None:
        return component
    path = Path(review_file).resolve()
    try:
        review = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f'{task_id}: invalid quality-review file.') from exc
    if (set(review) != {'schema', 'task_id', 'native_output_sha256',
                        'decision', 'reviewer', 'rationale'}
            or review.get('schema') != 1
            or review.get('task_id') != task_id
            or review.get('native_output_sha256') != component.source_sha256
            or review.get('decision') != 'accept'
            or not isinstance(review.get('reviewer'), str)
            or not review['reviewer'].strip()
            or not isinstance(review.get('rationale'), str)
            or not review['rationale'].strip()):
        raise ValueError(f'{task_id}: quality review does not accept this '
                         'exact native output.')
    canonical = json.dumps(review, sort_keys=True,
                           separators=(',', ':')).encode()
    review_sha256 = hashlib.sha256(canonical).hexdigest()
    provenance = hashlib.sha256(json.dumps({
        'native_output_sha256': component.source_sha256,
        'review_sha256': review_sha256,
    }, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return ComponentResult(
        key=component.key, value_hartree=component.value_hartree,
        quantity=component.quantity, method=component.method,
        basis=component.basis, backend=component.backend,
        state_id=component.state_id, charge=component.charge,
        multiplicity=component.multiplicity,
        geometry_sha256=component.geometry_sha256,
        source_sha256=provenance,
        source=f'{component.source}; reviewed by {path}',
        settings={**component.settings,
                  'quality_review': {
                      'reviewer': review['reviewer'].strip(),
                      'rationale': review['rationale'].strip(),
                      'review_file': str(path),
                      'review_sha256': review_sha256,
                      'native_output_sha256': component.source_sha256,
                  }},
        review_required=False)


def assemble_profiled_anl0_f12(
        interface_run, higher_order_run, corrections_run, *, state_id,
        spin_orbit_hartree, spin_orbit_source, spin_orbit_backend='manual',
        vpt2_review=None, base_run=None):
    """Assemble one provenance-checked profiled ANL0-F12 0 K energy.

    Every geometry-dependent term must trace to the accepted L2 geometry in
    ``interface_run`` and the accepted current L3 geometry in ``base_run``.
    For a fully current interface graph, omitting ``base_run`` uses that graph
    for both roles. Higher-order and correction graphs must identify the same
    completed L3 result.
    """
    interface_run, interface_spec, _ = _load(interface_run)
    if interface_spec.get('name') != 'anl1-f12-non-mrcc-interface-validation':
        raise ValueError('Interface run has the wrong workflow type.')
    if not isinstance(state_id, str) or not state_id.strip():
        raise ValueError('state_id must be a nonempty string.')
    molecule = interface_spec['molecule']
    charge = molecule.get('charge', 0)
    multiplicity = molecule.get('multiplicity', 1)
    l2 = molecule_from_completed_run(interface_run, 'l2_geometry')
    l2_hash = l2['source']['geometry_sha256']
    if base_run is None:
        base_run = interface_run
    else:
        base_run, base_spec, _ = _load(base_run)
        if base_spec.get('name') != 'anl-current-base-validation':
            raise ValueError('Base run has the wrong workflow type.')
        _same_imported_geometry(
            interface_run, base_run, l2_hash, 'L2', 'l2_geometry')
    l3 = molecule_from_completed_run(base_run, 'l3_geometry')
    l3_hash = l3['source']['geometry_sha256']
    _same_imported_l3_geometry(base_run, higher_order_run, l3_hash)
    _same_imported_l3_geometry(base_run, corrections_run, l3_hash)

    equation = recipe('ANL0-F12', vpt2_method='B2PLYP-D3BJ',
                      multiplicity=multiplicity)
    required = {item.key: item for item in equation.requirements}
    components = {
        'reference_cbs': cbs_task_component(
            base_run, 'f12_tz', 'f12_qz',
            requirement=required['reference_cbs'], state_id=state_id,
            lower_basis='cc-pVTZ-F12', upper_basis='cc-pVQZ-F12'),
        'harmonic_zpe': task_component(
            base_run, 'harmonic', key='harmonic_zpe', state_id=state_id),
        'vpt2_correction': task_component(
            interface_run, 'gaussian_vpt2', key='vpt2_correction',
            state_id=state_id),
        'dboc': task_component(
            base_run, 'cfour_dboc', key='dboc', state_id=state_id),
        'hoe_high': task_component(
            higher_order_run, 'ccsdtq_dz', key='hoe_high', state_id=state_id),
        'hoe_low': task_component(
            higher_order_run, 'ccsdt_dz', key='hoe_low', state_id=state_id),
    }
    components['vpt2_correction'] = _apply_quality_review(
        components['vpt2_correction'], 'gaussian_vpt2', vpt2_review)
    components['core_valence_cbs'] = core_valence_task_component(
        corrections_run, all_electron_lower='cv_ae_tz',
        all_electron_upper='cv_ae_qz', frozen_core_lower='cv_fc_tz',
        frozen_core_upper='cv_fc_qz',
        requirement=required['core_valence_cbs'], state_id=state_id)
    components['scalar_relativistic'] = scalar_relativistic_task_component(
        corrections_run, 'rel_dkh', 'rel_nonrel',
        requirement=required['scalar_relativistic'], state_id=state_id)
    components['spin_orbit'] = state_correction_component(
        requirement=required['spin_orbit'], value_hartree=spin_orbit_hartree,
        state_id=state_id, charge=charge, multiplicity=multiplicity,
        source=spin_orbit_source, backend=spin_orbit_backend)
    result = equation.evaluate(
        components, state_id=state_id, charge=charge,
        multiplicity=multiplicity,
        geometry_hashes={'l2': l2_hash, 'l3': l3_hash})
    return {
        'schema': 1, 'status': 'complete', 'state_id': state_id,
        'charge': charge, 'multiplicity': multiplicity,
        'recipe': result.recipe, 'recipe_version': result.recipe_version,
        'electronic_hartree': result.electronic_hartree,
        'zero_point_hartree': result.zero_point_hartree,
        'zero_k_hartree': result.zero_k_hartree,
        'geometry_sha256': {'l2': l2_hash, 'l3': l3_hash},
        'components': {key: asdict(value)
                       for key, value in result.components.items()},
        'electronic_terms': [asdict(term) for term in result.electronic_terms],
        'zero_point_terms': [asdict(term) for term in result.zero_point_terms],
    }


def audit_interface_run(run_dir):
    """Reparse every declared native result and report the honest boundary."""
    _, spec, state = _load(run_dir)
    statuses = {task['id']: state['tasks'].get(task['id'], {}).get(
        'status', 'waiting') for task in spec['tasks']}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'Interface run is incomplete: {statuses}')
    parsed = {}
    for task in spec['tasks']:
        if task.get('result_parser'):
            *_, result = _verified_task_result(run_dir, task['id'])
            parsed[task['id']] = result
    equation = recipe(
        'ANL0-F12', vpt2_method='B2PLYP-D3BJ',
        multiplicity=spec['molecule'].get('multiplicity', 1))
    requirement = next(item for item in equation.requirements
                       if item.key == 'reference_cbs')
    legacy = sorted(
        task['id'] for task in spec['tasks']
        if task.get('result_parser', {}).get('kind') in
        ('molpro_energy', 'molpro_harmonic')
        and legacy_molpro_parser(task['result_parser'],
                                 task.get('input_template', '')))
    if legacy:
        lower = parsed['f12_tz']['energy_hartree']
        upper = parsed['f12_qz']['energy_hartree']
        diagnostic = two_point_cbs(
            lower, upper,
            upper_cardinal=requirement.settings['upper_cardinal'],
            power=requirement.settings['extrapolation_power'])
        unavailable = list(spec['intent']['unavailable_components'])
        unavailable.append(
            'current unrestricted Molpro geometry, reference, and harmonic components')
        return {
            'status': 'interface_complete_legacy_recipe_incompatible',
            'requested_ladder_head': 'ANL1-F12',
            'mrcc_enabled': False,
            'task_statuses': statuses,
            'parsed_kinds': {key: value['kind']
                             for key, value in parsed.items()},
            'legacy_molpro_tasks': legacy,
            'legacy_f12_cbs_hartree_diagnostic': diagnostic,
            'unavailable_components': unavailable,
        }
    reference = cbs_task_component(
        run_dir, 'f12_tz', 'f12_qz', requirement=requirement,
        state_id='validation-state', lower_basis='cc-pVTZ-F12',
        upper_basis='cc-pVQZ-F12')
    return {
        'status': 'interface_complete_recipe_incomplete',
        'requested_ladder_head': 'ANL1-F12',
        'mrcc_enabled': False,
        'task_statuses': statuses,
        'parsed_kinds': {key: value['kind'] for key, value in parsed.items()},
        'verified_f12_cbs_hartree': reference.value_hartree,
        'unavailable_components': spec['intent']['unavailable_components'],
    }


def audit_current_base_run(run_dir):
    """Audit the current unrestricted L3/base continuation graph."""
    _, spec, state = _load(run_dir)
    if spec.get('name') != 'anl-current-base-validation':
        raise ValueError('Run is not a current ANL base validation graph.')
    statuses = {task['id']: state['tasks'].get(task['id'], {}).get(
        'status', 'waiting') for task in spec['tasks']}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'Current base run is incomplete: {statuses}')
    equation = recipe(
        'ANL0-F12', vpt2_method='B2PLYP-D3BJ',
        multiplicity=spec['molecule'].get('multiplicity', 1))
    required = {item.key: item for item in equation.requirements}
    reference = cbs_task_component(
        run_dir, 'f12_tz', 'f12_qz',
        requirement=required['reference_cbs'], state_id='base-validation-state',
        lower_basis='cc-pVTZ-F12', upper_basis='cc-pVQZ-F12')
    harmonic = task_component(
        run_dir, 'harmonic', key='harmonic_zpe',
        state_id='base-validation-state')
    dboc = task_component(
        run_dir, 'cfour_dboc', key='dboc',
        state_id='base-validation-state')
    return {
        'status': 'current_base_interface_complete',
        'task_statuses': statuses,
        'reference_cbs_hartree': reference.value_hartree,
        'harmonic_zpe_hartree': harmonic.value_hartree,
        'dboc_hartree': dboc.value_hartree,
        'geometry_sha256': reference.geometry_sha256,
    }


def audit_kinbot_run(run_dir, reaction, *, parent=None, hir_points=0,
                     require_rotdpy=False,
                     require_rotdpy_execution=False):
    """Gate downstream ANL work on an accepted KinBot reaction result."""
    run_dir = Path(run_dir).resolve()
    monitor = run_dir / 'kinbot_monitor.out'
    if not monitor.is_file():
        raise RuntimeError('KinBot did not write kinbot_monitor.out.')
    matches = []
    for line in monitor.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[2] == reaction:
            matches.append(fields)
    if len(matches) != 1 or matches[0][0] != '-1':
        raise RuntimeError(f'{reaction}: expected one accepted channel in '
                           f'kinbot_monitor.out, found {matches!r}.')
    products = matches[0][3:]
    if len(products) != 2:
        raise RuntimeError(f'{reaction}: expected two product entries, found '
                           f'{products!r}.')

    normal_hir = []
    if hir_points:
        database = run_dir / 'kinbot.db'
        if not database.is_file():
            raise RuntimeError('KinBot database is missing.')
        prefix = f'hir/{parent}_hir_' if parent else 'hir/'
        for row in connect(str(database)).select():
            if (getattr(row, 'name', '').startswith(prefix)
                    and row.data.get('status') == 'normal'):
                normal_hir.append(row.name)
        if len(set(normal_hir)) < hir_points:
            raise RuntimeError(f'Expected at least {hir_points} accepted '
                               f'hindered-rotor points for {parent}, found '
                               f'{len(set(normal_hir))}.')
        log = run_dir / 'kinbot.log'
        if (log.is_file()
                and 'will be treated as harmonic oscillators' in log.read_text()):
            raise RuntimeError('KinBot demoted a requested hindered rotor to '
                               'a harmonic oscillator.')

    rotdpy = None
    correction = None
    if require_rotdpy_execution:
        require_rotdpy = True
    if require_rotdpy:
        correction = run_dir / 'vrctst' / f'corr_{reaction}.json'
        rotdpy = run_dir / 'rotdPy' / f'{reaction}.py'
        if not correction.is_file() or not rotdpy.is_file():
            raise RuntimeError(f'{reaction}: VRC correction or rotdPy input '
                               'is missing.')
        payload = json.loads(correction.read_text())
        required = {'dist', 'e_samp', 'e_high', 'scan_ref', 'ra',
                    'e_inf_samp', 'e_inf_high', 'frags_atom', 'frags_geom',
                    'frags_mult'}
        if required - payload.keys():
            raise RuntimeError(f'{reaction}: incomplete VRC correction record.')
        if (len(payload['dist']) != len(payload['e_samp'])
                or len(payload['dist']) != len(payload['e_high'])):
            raise RuntimeError(f'{reaction}: inconsistent VRC correction arrays.')
        if not rotdpy.read_text().strip():
            raise RuntimeError(f'{reaction}: rotdPy input is empty.')
        try:
            compile(rotdpy.read_text(), str(rotdpy), 'exec')
        except SyntaxError as error:
            raise RuntimeError(f'{reaction}: rotdPy input is invalid: '
                               f'{error}') from error
        if require_rotdpy_execution:
            execution_file = rotdpy.with_name(f'{reaction}.execution.json')
            if not execution_file.is_file():
                raise RuntimeError(f'{reaction}: rotdPy was not executed.')
            execution = json.loads(execution_file.read_text())
            input_hash = hashlib.sha256(rotdpy.read_bytes()).hexdigest()
            if (execution.get('status') != 'complete'
                    or execution.get('returncode') != 0
                    or execution.get('input_sha256') != input_hash):
                raise RuntimeError(f'{reaction}: rotdPy execution is incomplete '
                                   'or does not match its input.')
            from kinbot.rotdpy import read_result
            rotdpy_result = read_result(rotdpy)
        else:
            execution_file = None
            rotdpy_result = None
    else:
        execution_file = None
        rotdpy_result = None

    return {
        'status': 'kinbot_reaction_complete',
        'reaction': reaction,
        'products': products,
        'normal_hir_points': len(set(normal_hir)),
        'vrc_correction': str(correction) if correction else None,
        'rotdpy_input': str(rotdpy) if rotdpy else None,
        'rotdpy_execution': (str(execution_file)
                             if execution_file else None),
        'rotdpy_surfaces': (rotdpy_result['surface_count']
                            if rotdpy_result else 0),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Prepare or audit ANL native-interface validations')
    commands = parser.add_subparsers(dest='action', required=True)
    build = commands.add_parser('from-db')
    build.add_argument('database', type=Path)
    build.add_argument('job')
    build.add_argument('spec', type=Path)
    build.add_argument('--charge', type=int, default=0)
    build.add_argument('--multiplicity', type=int, default=1)
    build.add_argument('--max-nodes', type=int, default=3)
    build.add_argument('--partition')
    stage = commands.add_parser('prepare-from-db')
    stage.add_argument('database', type=Path)
    stage.add_argument('job')
    stage.add_argument('run_dir', type=Path)
    stage.add_argument('--charge', type=int, default=0)
    stage.add_argument('--multiplicity', type=int, default=1)
    stage.add_argument('--max-nodes', type=int, default=3)
    stage.add_argument('--partition')
    audit = commands.add_parser('audit')
    audit.add_argument('run_dir', type=Path)
    prepare_base = commands.add_parser('prepare-current-base-from-run')
    prepare_base.add_argument('source_run', type=Path)
    prepare_base.add_argument('run_dir', type=Path)
    prepare_base.add_argument('--geometry-task', default='l2_geometry')
    prepare_base.add_argument('--max-nodes', type=int, default=3)
    prepare_base.add_argument('--partition')
    audit_base = commands.add_parser('audit-current-base')
    audit_base.add_argument('run_dir', type=Path)
    higher = commands.add_parser('higher-order-from-db')
    higher.add_argument('database', type=Path)
    higher.add_argument('job')
    higher.add_argument('spec', type=Path)
    higher.add_argument('--charge', type=int, default=0)
    higher.add_argument('--multiplicity', type=int, default=1)
    higher.add_argument('--max-nodes', type=int, default=3)
    higher.add_argument('--partition')
    higher.add_argument('--mrcc-command', default='dmrcc')
    higher.add_argument('--task', dest='task_ids', action='append',
                        choices=_HIGHER_ORDER_TASK_IDS)
    prepare_higher = commands.add_parser('prepare-higher-order-from-db')
    prepare_higher.add_argument('database', type=Path)
    prepare_higher.add_argument('job')
    prepare_higher.add_argument('run_dir', type=Path)
    prepare_higher.add_argument('--charge', type=int, default=0)
    prepare_higher.add_argument('--multiplicity', type=int, default=1)
    prepare_higher.add_argument('--max-nodes', type=int, default=3)
    prepare_higher.add_argument('--partition')
    prepare_higher.add_argument('--mrcc-command', default='dmrcc')
    prepare_higher.add_argument('--task', dest='task_ids', action='append',
                                choices=_HIGHER_ORDER_TASK_IDS)
    prepare_higher_run = commands.add_parser('prepare-higher-order-from-run')
    prepare_higher_run.add_argument('source_run', type=Path)
    prepare_higher_run.add_argument('run_dir', type=Path)
    prepare_higher_run.add_argument('--geometry-task', default='l3_geometry')
    prepare_higher_run.add_argument('--max-nodes', type=int, default=3)
    prepare_higher_run.add_argument('--partition')
    prepare_higher_run.add_argument('--mrcc-command', default='dmrcc')
    prepare_higher_run.add_argument('--task', dest='task_ids', action='append',
                                    choices=_HIGHER_ORDER_TASK_IDS)
    audit_higher = commands.add_parser('audit-higher-order')
    audit_higher.add_argument('run_dir', type=Path)
    prepare_corrections = commands.add_parser(
        'prepare-common-corrections-from-run')
    prepare_corrections.add_argument('source_run', type=Path)
    prepare_corrections.add_argument('run_dir', type=Path)
    prepare_corrections.add_argument('--geometry-task', default='l3_geometry')
    prepare_corrections.add_argument('--max-nodes', type=int, default=3)
    prepare_corrections.add_argument('--partition')
    audit_corrections = commands.add_parser('audit-common-corrections')
    audit_corrections.add_argument('run_dir', type=Path)
    prepare_post = commands.add_parser('prepare-post-geometry-from-run')
    prepare_post.add_argument('source_run', type=Path)
    prepare_post.add_argument('run_dir', type=Path)
    prepare_post.add_argument('--geometry-task', default='l3_geometry')
    prepare_post.add_argument('--max-nodes', type=int, default=3)
    prepare_post.add_argument('--partition')
    prepare_post.add_argument('--mrcc-command', default='dmrcc')
    audit_post = commands.add_parser('audit-post-geometry')
    audit_post.add_argument('run_dir', type=Path)
    audit_anl0_post = commands.add_parser('audit-anl0-post-geometry')
    audit_anl0_post.add_argument('run_dir', type=Path)
    assemble = commands.add_parser('assemble-anl0-f12')
    assemble.add_argument('interface_run', type=Path)
    assemble.add_argument('higher_order_run', type=Path)
    assemble.add_argument('corrections_run', type=Path)
    assemble.add_argument('output', type=Path)
    assemble.add_argument('--state-id', required=True)
    assemble.add_argument('--spin-orbit-hartree', required=True, type=float)
    assemble.add_argument('--spin-orbit-source', required=True)
    assemble.add_argument('--spin-orbit-backend', default='manual',
                          choices=('manual', 'known_zero', 'table',
                                   'calculated'))
    assemble.add_argument('--vpt2-review', type=Path)
    assemble.add_argument('--base-run', type=Path)
    gate = commands.add_parser('gate-kinbot')
    gate.add_argument('run_dir', type=Path)
    gate.add_argument('reaction')
    gate.add_argument('--parent')
    gate.add_argument('--hir-points', type=int, default=0)
    gate.add_argument('--require-rotdpy', action='store_true')
    gate.add_argument('--require-rotdpy-execution', action='store_true')
    args = parser.parse_args(argv)
    if args.action == 'audit':
        print(json.dumps(audit_interface_run(args.run_dir), indent=2,
                         sort_keys=True))
        return 0
    if args.action == 'audit-current-base':
        print(json.dumps(audit_current_base_run(args.run_dir), indent=2,
                         sort_keys=True))
        return 0
    if args.action == 'audit-higher-order':
        print(json.dumps(audit_higher_order_run(args.run_dir), indent=2,
                         sort_keys=True))
        return 0
    if args.action == 'audit-common-corrections':
        print(json.dumps(audit_common_corrections_run(args.run_dir), indent=2,
                         sort_keys=True))
        return 0
    if args.action == 'audit-post-geometry':
        print(json.dumps(audit_post_geometry_run(args.run_dir), indent=2,
                         sort_keys=True))
        return 0
    if args.action == 'audit-anl0-post-geometry':
        print(json.dumps(audit_anl0_post_geometry_run(args.run_dir), indent=2,
                         sort_keys=True))
        return 0
    if args.action == 'assemble-anl0-f12':
        payload = assemble_profiled_anl0_f12(
            args.interface_run, args.higher_order_run, args.corrections_run,
            state_id=args.state_id,
            spin_orbit_hartree=args.spin_orbit_hartree,
            spin_orbit_source=args.spin_orbit_source,
            spin_orbit_backend=args.spin_orbit_backend,
            vpt2_review=args.vpt2_review, base_run=args.base_run)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)
                               + '\n')
        print(args.output.resolve())
        return 0
    if args.action == 'gate-kinbot':
        print(json.dumps(audit_kinbot_run(
            args.run_dir, args.reaction, parent=args.parent,
            hir_points=args.hir_points,
            require_rotdpy=args.require_rotdpy,
            require_rotdpy_execution=args.require_rotdpy_execution),
            indent=2, sort_keys=True))
        return 0
    from_run = args.action in ('prepare-current-base-from-run',
                               'prepare-higher-order-from-run',
                               'prepare-common-corrections-from-run',
                               'prepare-post-geometry-from-run')
    molecule = (molecule_from_completed_run(
        args.source_run, geometry_task=args.geometry_task) if from_run
        else molecule_from_database(
            args.database, args.job, charge=args.charge,
            multiplicity=args.multiplicity))
    is_higher = args.action in ('higher-order-from-db',
                                'prepare-higher-order-from-db',
                                'prepare-higher-order-from-run')
    is_corrections = args.action == 'prepare-common-corrections-from-run'
    is_post = args.action == 'prepare-post-geometry-from-run'
    is_base = args.action == 'prepare-current-base-from-run'
    if is_base:
        spec = current_base_validation_spec(
            molecule, max_nodes=args.max_nodes, partition=args.partition)
    elif is_higher:
        spec = higher_order_validation_spec(
            molecule, max_nodes=args.max_nodes, partition=args.partition,
            mrcc_command=args.mrcc_command, task_ids=args.task_ids)
    elif is_corrections:
        spec = common_corrections_validation_spec(
            molecule, max_nodes=args.max_nodes, partition=args.partition)
    elif is_post:
        spec = post_geometry_validation_spec(
            molecule, max_nodes=args.max_nodes, partition=args.partition,
            mrcc_command=args.mrcc_command)
    else:
        spec = interface_validation_spec(
            molecule, max_nodes=args.max_nodes, partition=args.partition)
    if args.action in ('from-db', 'higher-order-from-db'):
        args.spec.write_text(json.dumps(spec, indent=2) + '\n')
        print(args.spec.resolve())
    else:
        print(prepare(spec, args.run_dir))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
