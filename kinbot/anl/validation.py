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
import re

from ase.db import connect
from ase.data import atomic_numbers
from ase.io import read

from kinbot.anl.dispatch import (_load, _verify_execution, _verify_stage_files,
                                 prepare)
from kinbot.anl.extrapolation import two_point_cbs
from kinbot.anl.model import (ComponentRequirement, ComponentResult,
                              IncompleteRecipeError)
from kinbot.anl.recipes import (ANL0_NATIVE_REFERENCE_BASES,
                                ANL1_NATIVE_REFERENCE_BASES,
                                f12_reference_bases, recipe)
from kinbot.anl.results import legacy_molpro_parser
from kinbot.anl.tasks import (higher_order_task, molpro_ccsdt_command,
                              molpro_task as native_molpro_task)
from kinbot.theory import (convergence_exponent, fortran_threshold, molpro_cc_thresholds)
from kinbot.anl.workflow import (_verified_task_result,
                                 atomic_zero_vibrational_component,
                                 cbs_task_component,
                                 core_valence_task_component,
                                 rank_exact_core_valence_component,
                                 rank_exact_higher_order_correction,
                                 rank_exact_higher_order_component,
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
                  partition=None, scf_convergence=1.e-10,
                  cc_convergence=1.e-8):
    return native_molpro_task(
        ident, body, geometry_from=geometry_from, walltime=walltime,
        max_cores=max_cores, parser=parser, partition=partition,
        scf_convergence=scf_convergence,
        cc_convergence=cc_convergence)


def interface_validation_spec(molecule, *, max_nodes=3, partition=None,
                              stationary_point='minimum'):
    """Return a molecule-independent Gaussian/Molpro/CFOUR task graph."""
    if not isinstance(molecule, dict):
        raise TypeError('molecule must be an object.')
    if stationary_point not in ('minimum', 'transition_state'):
        raise ValueError(
            'stationary_point must be minimum or transition_state.')
    saddle_order = 1 if stationary_point == 'transition_state' else 0
    l2_profile = {
        'calculator': 'gaussian', 'method': 'B2PLYP',
        'basis': 'cc-pVTZ', 'command': 'g16',
        'calculator_kwargs': {
            'EmpiricalDispersion': 'GD3BJ', 'Symm': 'None',
            'scf': 'xqc', 'integral': 'UltraFine'},
        'optimizer': 'sella',
    }
    reference = 'RHF' if molecule.get('multiplicity', 1) == 1 else 'ROHF'
    internal_coordinates = len(molecule.get('symbols', ())) > 2
    conventional = molpro_ccsdt_command(molecule.get('multiplicity', 1))
    electrons = _electron_count(molecule)
    effective = ('HF' if electrons == 1 else
                 'CCSD' if electrons == 2 else None)
    conventional_command = ({'HF': '', 'CCSD': 'uccsd'}[effective]
                            if effective else conventional)
    f12 = ({'HF': '', 'CCSD': 'uccsd-f12b'}[effective]
           if effective else 'uccsd(t)-f12b,scale_trip=1')
    rank_parser = ({'rank_exact_electrons': electrons,
                    'effective_method': effective}
                   if effective else {})
    f12_rank_parser = ({**rank_parser,
                        'effective_method': ('HF' if electrons == 1
                                             else 'CCSD-F12b')}
                       if effective else {})
    tasks = [
        {
            'id': 'l2_geometry', 'kind': 'ase_optimize',
            'geometry_from': 'initial', 'geometry_output': 'final.xyz',
            'resources': _resources('06:00:00', max_cores=8,
                                    partition=partition),
            'profile': l2_profile,
            'optimizer': {'fmax': 0.0005, 'steps': 160,
                          'order': saddle_order,
                          'sella_kwargs': {
                              'internal': internal_coordinates}},
        },
        {
            'id': 'l3_geometry', 'kind': 'ase_optimize',
            'geometry_from': 'l2_geometry', 'geometry_output': 'final.xyz',
            'resources': _resources('24:00:00', max_cores=12,
                                    min_stack_mw=1024, partition=partition),
            'profile': {
                'calculator': 'molpro', 'method': 'CCSD(T)',
                'basis': 'cc-pVTZ', 'command': 'molpro',
                'calculator_kwargs': ({'effective_method': effective}
                                      if effective else {}),
                'optimizer': 'sella'},
            'optimizer': {'fmax': 0.03, 'steps': 100,
                          'order': saddle_order,
                          'sella_kwargs': {
                              'internal': internal_coordinates}},
        },
        _molpro_task(
            'harmonic',
            f'basis=cc-pVTZ\nrhf\n{conventional_command}\nfrequencies,numerical\n',
            walltime='24:00:00', max_cores=8,
            parser={'kind': 'molpro_harmonic', 'basis': 'cc-pVTZ',
                    'reference': reference,
                    'stationary_point': stationary_point, **rank_parser},
            partition=partition),
        _molpro_task(
            'f12_tz',
            f'basis=cc-pVTZ-F12\nrhf\n{f12}\n'
            'kb_f12b=energy\n',
            walltime='12:00:00', max_cores=12,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)-F12b',
                    'basis': 'cc-pVTZ-F12', 'reference': reference,
                    **f12_rank_parser},
            partition=partition),
        _molpro_task(
            'f12_qz',
            f'basis=cc-pVQZ-F12\nrhf\n{f12}\n'
            'kb_f12b=energy\n',
            walltime='24:00:00', max_cores=12,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)-F12b',
                    'basis': 'cc-pVQZ-F12', 'reference': reference,
                    **f12_rank_parser},
            partition=partition),
        _molpro_task(
            'ccsdt_dz',
            f'basis=cc-pVDZ\nrhf\n{conventional_command}\nkb_dz_energy=energy\n',
            walltime='08:00:00', max_cores=8,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': 'cc-pVDZ', 'reference': reference},
            partition=partition),
        _molpro_task(
            'l3_sp_l2',
            f'basis=cc-pVTZ\nrhf\n{conventional_command}\n'
            'kb_l3_sp_l2=energy\n',
            walltime='12:00:00', max_cores=12,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': 'cc-pVTZ', 'reference': reference,
                    **rank_parser},
            geometry_from='l2_geometry',
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
                'MULTIPLICITY={{MULT}}\nSCF_CONV=10\nMEM_UNIT=MB\n'
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
            'resources': _resources('24:00:00', max_cores=8,
                                    partition=partition),
            'input_name': 'vpt2.com',
            'input_template': (
                '%nprocshared={{CORES}}\n%mem={{WORK_MEMORY_MB}}MB\n'
                '#p B2PLYP/cc-pVTZ Freq=Anharmonic Symmetry=None SCF=XQC '
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
    if effective:
        # CCSD(T)-and-higher increments are identically zero for one- and
        # two-electron states.  Reference, geometry, frequency, DBOC, and
        # relativistic work remain independent physical contributions.
        tasks = [task for task in tasks if task['id'] != 'ccsdt_dz']
    if len(molecule.get('symbols', ())) == 1:
        # An atom has no vibrational normal modes.  Keep both electronic
        # geometry evaluations and every electronic correction, while the
        # assembler supplies hash-bound exact zeros for ZPE and VPT2.
        tasks = [task for task in tasks
                 if task['id'] not in ('harmonic', 'gaussian_vpt2')]
    if stationary_point == 'transition_state':
        # VPT2 is a minimum thermochemistry correction.  A tight saddle uses
        # its one-imaginary-mode harmonic ZPE, paired with a harmonic-ZPE
        # reactant record from the same electronic recipe.
        tasks = [task for task in tasks if task['id'] != 'gaussian_vpt2']
    # Every electronic correction is a single point and can use the accepted
    # L2 geometry when the optional L3 optimization is unavailable. Only the
    # L3 harmonic calculation intrinsically belongs to the failed L3
    # structure; the accepted L2 frequencies remain available separately.
    for task in tasks:
        if (task.get('geometry_from') == 'l3_geometry'
                and task['id'] != 'harmonic'):
            task['geometry_fallback'] = 'l2_geometry'
    by_id = {task['id']: task for task in tasks}
    if 'l3_geometry' in by_id:
        _component_fallback(
            by_id['l3_geometry'], 'l3_optimized_geometry',
            failure_skips=tuple(
                task['id'] for task in tasks
                if (task.get('geometry_from') == 'l3_geometry'
                    and 'geometry_fallback' not in task)))
    for family, identifiers in {
            'reference_cbs_f12': ('f12_tz', 'f12_qz'),
            'harmonic_zpe': ('harmonic',),
            'higher_order_anl0': ('ccsdt_dz',),
            'dboc': ('cfour_dboc',),
            'vpt2_correction': ('gaussian_vpt2',),
            'l3_single_point_l2_geometry': ('l3_sp_l2',),
    }.items():
        present = tuple(ident for ident in identifiers if ident in by_id)
        for ident in present:
            _component_fallback(
                by_id[ident], family,
                failure_skips=tuple(other for other in present
                                    if other != ident))
    return {
        'schema': 1, 'name': 'anl1-f12-non-mrcc-interface-validation',
        'molecule': molecule,
        'limits': {'max_nodes': max_nodes},
        'tasks': tasks,
        'intent': {
            'requested_ladder_head': 'ANL1-F12',
            'mrcc_enabled': False,
            'claim': 'interface-validation-only',
            'stationary_point': stationary_point,
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
    if len(molecule.get('symbols', ())) == 1:
        keep.remove('harmonic')
    tasks = [deepcopy(task) for task in full['tasks'] if task['id'] in keep]
    for task in tasks:
        if task['id'] == 'l3_geometry':
            task['geometry_from'] = 'initial'
        if task.get('geometry_fallback') == 'l2_geometry':
            task['geometry_fallback'] = 'initial'
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


def _component_fallback(task, *families, failure_skips=()):
    """Mark optional component work and its still-useful failure boundary."""
    task['failure_policy'] = 'skip_component'
    task['component_families'] = list(dict.fromkeys(
        [*task.get('component_families', ()), *families]))
    task['failure_skips'] = list(dict.fromkeys(
        [*task.get('failure_skips', ()), *failure_skips]))
    return task


def _electron_count(molecule):
    """Return a validated molecular electron count from dispatcher symbols."""
    symbols = molecule.get('symbols') if isinstance(molecule, dict) else None
    charge = molecule.get('charge', 0) if isinstance(molecule, dict) else None
    if (not isinstance(symbols, list) or not symbols
            or isinstance(charge, bool) or not isinstance(charge, int)):
        raise ValueError('Molecule needs symbols and an integer charge.')
    try:
        electrons = sum(atomic_numbers[symbol] for symbol in symbols) - charge
    except KeyError as exc:
        raise ValueError(f'Unknown element symbol {exc.args[0]!r}.') from exc
    if electrons < 1:
        raise ValueError('Molecule must contain at least one electron.')
    return electrons


def _rank_exact_higher_order(molecule):
    """Whether triple and higher excitations are rigorously absent."""
    return _electron_count(molecule) <= 2


def _rank_exact_core_valence(molecule):
    """Whether the ANL frozen-core space is empty."""
    return all(atomic_numbers[symbol] <= 2 for symbol in molecule['symbols'])


def _anl0_reference_tasks(molecule, *, partition=None,
                          geometry_from='l3_geometry'):
    """Build the fully augmented QZ/5Z reference pair for ANL0."""
    reference = 'RHF' if molecule.get('multiplicity', 1) == 1 else 'ROHF'
    electrons = _electron_count(molecule)
    effective = 'HF' if electrons == 1 else 'CCSD' if electrons == 2 else None
    command = ({'HF': '', 'CCSD': 'uccsd'}[effective]
               if effective else molpro_ccsdt_command(
                   molecule.get('multiplicity', 1)))
    rank_parser = ({'rank_exact_electrons': electrons,
                    'effective_method': effective}
                   if effective else {})
    tasks = []
    for ident, basis, walltime, cores in (
            ('anl0_qz', ANL0_NATIVE_REFERENCE_BASES[0], '24:00:00', 12),
            ('anl0_5z', ANL0_NATIVE_REFERENCE_BASES[1], '24:00:00', 12)):
        task = _molpro_task(
            ident, f'basis={basis}\nrhf\n{command}\nkb_{ident}=energy\n',
            geometry_from=geometry_from, walltime=walltime, max_cores=cores,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': basis, 'reference': reference, **rank_parser},
            partition=partition)
        tasks.append(task)
    for task in tasks:
        _component_fallback(
            task, 'reference_cbs_anl0',
            failure_skips=tuple(other['id'] for other in tasks
                                if other is not task))
    return tasks


def _anl1_reference_tasks(molecule, *, partition=None,
                          geometry_from='l3_geometry'):
    """Build the fully augmented 5Z/6Z reference pair."""
    reference = 'RHF' if molecule.get('multiplicity', 1) == 1 else 'ROHF'
    electrons = _electron_count(molecule)
    effective = 'HF' if electrons == 1 else 'CCSD' if electrons == 2 else None
    command = ({'HF': '', 'CCSD': 'uccsd'}[effective]
               if effective else molpro_ccsdt_command(
                   molecule.get('multiplicity', 1)))
    rank_parser = ({'rank_exact_electrons': electrons,
                    'effective_method': effective}
                   if effective else {})
    tasks = []
    for ident, basis in zip(
            ('anl1_5z', 'anl1_6z'), ANL1_NATIVE_REFERENCE_BASES):
        task = _molpro_task(
            ident, f'basis={basis}\nrhf\n{command}\nkb_{ident}=energy\n',
            geometry_from=geometry_from, walltime='7-00:00:00',
            max_cores=12,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': basis, 'reference': reference, **rank_parser},
            partition=partition)
        tasks.append(task)
    for task in tasks:
        _component_fallback(
            task, 'reference_cbs_anl1',
            failure_skips=tuple(other['id'] for other in tasks
                                if other is not task))
    return tasks


def _anl1_f12_reference_tasks(molecule, *, partition=None,
                              geometry_from='l3_geometry',
                              f12_augmented=False):
    """Build the 5Z-F12 task completing the Ram-style QZ/5Z pair.

    The ordinary interface graph already contains the physically identical
    QZ-F12 task.  Reusing that task avoids a duplicate calculation.
    """
    reference = 'RHF' if molecule.get('multiplicity', 1) == 1 else 'ROHF'
    electrons = _electron_count(molecule)
    effective = 'HF' if electrons == 1 else 'CCSD-F12b' \
        if electrons == 2 else None
    command = ({'HF': '', 'CCSD-F12b': 'uccsd-f12b'}[effective]
               if effective else 'uccsd(t)-f12b,scale_trip=1')
    rank_parser = ({'rank_exact_electrons': electrons,
                    'effective_method': effective}
                   if effective else {})
    ident = 'f12_5z'
    basis = f12_reference_bases('ANL1-F12', f12_augmented)[1]
    task = _molpro_task(
        ident, f'basis={basis}\nrhf\n{command}\nkb_{ident}=energy\n',
        geometry_from=geometry_from, walltime='7-00:00:00', max_cores=12,
        parser={'kind': 'molpro_energy', 'method': 'CCSD(T)-F12b',
                'basis': basis, 'reference': reference, **rank_parser},
        partition=partition)
    return [_component_fallback(task, 'reference_cbs_dual_f12')]


def _anl1_harmonic_tasks(molecule, *, partition=None,
                         tz_geometry='l3_geometry_tz',
                         qz_geometry='l3_geometry',
                         stationary_point='minimum'):
    """Build the basis-matched TZ/QZ harmonic pair required by ANL1."""
    if len(molecule.get('symbols', ())) == 1:
        return []
    reference = 'RHF' if molecule.get('multiplicity', 1) == 1 else 'ROHF'
    electrons = _electron_count(molecule)
    effective = 'HF' if electrons == 1 else 'CCSD' if electrons == 2 else None
    command = ({'HF': '', 'CCSD': 'uccsd'}[effective]
               if effective else molpro_ccsdt_command(
                   molecule.get('multiplicity', 1)))
    rank_parser = ({'rank_exact_electrons': electrons,
                    'effective_method': effective}
                   if effective else {})
    tasks = []
    for ident, basis, geometry in (
            ('harmonic_tz', 'cc-pVTZ', tz_geometry),
            ('harmonic_qz', 'cc-pVQZ', qz_geometry)):
        task = _molpro_task(
            ident,
            f'basis={basis}\nrhf\n{command}\nfrequencies,numerical\n',
            geometry_from=geometry, walltime='7-00:00:00', max_cores=8,
            parser={'kind': 'molpro_harmonic', 'basis': basis,
                    'reference': reference,
                    'stationary_point': stationary_point, **rank_parser},
            partition=partition)
        tasks.append(task)
    for task in tasks:
        _component_fallback(
            task, 'harmonic_zpe_anl1',
            failure_skips=tuple(other['id'] for other in tasks
                                if other is not task))
    return tasks


def higher_order_validation_spec(molecule, *, max_nodes=3, partition=None,
                                 mrcc_command='dmrcc', task_ids=None):
    """Build the ANL1 higher-order interface probe when ranks exist.

    This graph checks native input generation, execution, parsing, and the
    cross-program correction arithmetic on one supplied geometry.  It is a
    focused backend test and does not claim a complete composite energy.
    One- and two-electron states return a valid calculation-free graph because
    every requested triple-and-higher increment is exactly zero.
    """
    if not isinstance(molecule, dict):
        raise TypeError('molecule must be an object.')
    multiplicity = molecule.get('multiplicity', 1)
    restricted_reference = 'RHF' if multiplicity == 1 else 'ROHF'
    conventional = molpro_ccsdt_command(multiplicity)
    common = dict(geometry_from='initial', partition=partition)
    tasks = [
        _molpro_task(
            'ccsdt_tz',
            f'basis=cc-pVTZ\nrhf\n{conventional}\nkb_ccsdt=energy\n',
            walltime='08:00:00', max_cores=8,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': 'cc-pVTZ',
                    'reference': restricted_reference, 'core': 'frozen'},
            **common),
        _molpro_task(
            'ccsdt_dz',
            f'basis=cc-pVDZ\nrhf\n{conventional}\nkb_ccsdt=energy\n',
            walltime='04:00:00', max_cores=8,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': 'cc-pVDZ',
                    'reference': restricted_reference, 'core': 'frozen'},
            **common),
        higher_order_task(
            'ccsdtq_tz', 'CCSDT(Q)', 'cc-pVTZ',
            multiplicity=multiplicity,
            walltime='24:00:00',
            max_cores=8, command=mrcc_command, **common),
        higher_order_task(
            'ccsdtq_dz', 'CCSDT(Q)', 'cc-pVDZ',
            multiplicity=multiplicity,
            walltime='12:00:00',
            max_cores=8, command=mrcc_command, **common),
        higher_order_task(
            'ccsdtqp_dz', 'CCSDTQ(P)', 'cc-pVDZ',
            multiplicity=multiplicity,
            walltime='7-00:00:00',
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
    by_id = {task['id']: task for task in tasks}
    if 'ccsdt_tz' in by_id:
        _component_fallback(
            by_id['ccsdt_tz'], 'higher_order_anl1_tz',
            failure_skips=tuple(
                ident for ident in ('ccsdtq_tz',) if ident in by_id))
    if 'ccsdtq_tz' in by_id:
        _component_fallback(by_id['ccsdtq_tz'], 'higher_order_anl1_tz')
    if 'ccsdt_dz' in by_id:
        _component_fallback(
            by_id['ccsdt_dz'], 'higher_order_anl0',
            failure_skips=tuple(
                ident for ident in ('ccsdtq_dz', 'ccsdtqp_dz')
                if ident in by_id))
    if 'ccsdtq_dz' in by_id:
        _component_fallback(
            by_id['ccsdtq_dz'], 'higher_order_anl0',
            'higher_order_anl1_dz',
            failure_skips=tuple(
                ident for ident in ('ccsdtqp_dz',) if ident in by_id))
    if 'ccsdtqp_dz' in by_id:
        _component_fallback(
            by_id['ccsdtqp_dz'], 'higher_order_anl1_dz')
    rank_exact = _rank_exact_higher_order(molecule)
    if rank_exact:
        tasks = []
    return {
        'schema': 1, 'name': 'anl1-higher-order-interface-validation',
        'molecule': molecule, 'limits': {'max_nodes': max_nodes},
        'tasks': tasks,
        'intent': {
            'claim': 'higher-order-interface-validation-only',
            'equation': ('CCSDT(Q)/TZ - CCSD(T)/TZ + CCSDTQ(P)/DZ '
                         '- CCSDT(Q)/DZ'),
            'backend_policy': {
                'closed_shell_ccsdt_q': 'direct-mrcc-rhf-unrestricted-cc',
                'open_shell_ccsdt_q': 'direct-mrcc-semicanonical-rohf',
                'closed_shell_ccsdtq_p': 'direct-mrcc-rhf-unrestricted-cc',
                'open_shell_ccsdtq_p': 'direct-mrcc-semicanonical-rohf',
            },
            'selected_tasks': list(requested),
            'scheduled_tasks': [task['id'] for task in tasks],
            'rank_exact_higher_order': ({
                'electron_count': _electron_count(molecule),
                'correction_hartree': 0.0,
                'omitted_tasks': sorted(set(requested) - {
                    task['id'] for task in tasks}),
            } if rank_exact else None),
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
        electrons = _electron_count(molecule)
        effective = ('HF' if electrons == 1 else
                     'CCSD' if electrons == 2 else None)
        command = (('' if effective == 'HF' else
                    '{uccsd;core}' if core == 'all-electron' else 'uccsd')
                   if effective else molpro_ccsdt_command(
                       molecule.get('multiplicity', 1),
                       all_electron=(core == 'all-electron')))
        body = (f'basis={basis}\n{dkh}rhf\n{command}\n'
                f'kb_{ident}=energy\n')
        return _molpro_task(
            ident, body, walltime=walltime, max_cores=8,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': basis, 'reference': reference, 'core': core,
                    'relativistic': relativistic,
                    **({'rank_exact_electrons': electrons,
                        'effective_method': effective}
                       if effective else {})}, **common)

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
    core_valence_exact = _rank_exact_core_valence(molecule)
    if core_valence_exact:
        tasks = [task for task in tasks if not task['id'].startswith('cv_')]
    by_id = {task['id']: task for task in tasks}
    cv_ids = tuple(ident for ident in (
        'cv_ae_tz', 'cv_ae_qz', 'cv_fc_tz', 'cv_fc_qz') if ident in by_id)
    for ident in cv_ids:
        _component_fallback(
            by_id[ident], 'core_valence_cbs',
            failure_skips=tuple(other for other in cv_ids if other != ident))
    rel_ids = tuple(ident for ident in ('rel_dkh', 'rel_nonrel')
                    if ident in by_id)
    for ident in rel_ids:
        _component_fallback(
            by_id[ident], 'scalar_relativistic',
            failure_skips=tuple(other for other in rel_ids if other != ident))
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
            'rank_exact_core_valence': ({
                'correction_hartree': 0.0,
                'reason': 'H/He-only species have no frozen inner-shell orbitals',
            } if core_valence_exact else None),
        },
    }


def post_geometry_validation_spec(molecule, *, max_nodes=3, partition=None,
                                  mrcc_command='dmrcc', anl0_only=False):
    """Build one globally throttled fan-out after the accepted L3 geometry.

    ``anl0_only`` omits the ANL1-only TZ quadruples and DZ pentuples probes.
    It leaves the complete ANL0/ANL0-F12 higher-order and common-correction
    graph, so an affordable chemistry validation cannot accidentally submit
    the week-long CCSDTQ(P) task.
    """
    rank_exact = _rank_exact_higher_order(molecule)
    higher_task_ids = (('ccsdt_dz', 'ccsdtq_dz')
                       if anl0_only else None)
    higher = higher_order_validation_spec(
        molecule, max_nodes=max_nodes, partition=partition,
        mrcc_command=mrcc_command, task_ids=higher_task_ids)
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
            'selected_higher_order_tasks':
                higher['intent']['selected_tasks'],
            'rank_exact_higher_order': ({
                'electron_count': _electron_count(molecule),
                'correction_hartree': 0.0,
                'reason': ('triple and higher excitations are absent for a '
                           'state containing at most two electrons'),
            } if rank_exact else None),
            'core_valence_equation':
                corrections['intent']['core_valence_equation'],
            'scalar_relativistic_equation':
                corrections['intent']['scalar_relativistic_equation'],
            'scheduling': ('one shared max_nodes limit for every independent '
                           'post-geometry task'),
        },
    }


def composite_validation_spec(molecule, *, max_nodes=3, partition=None,
                              mrcc_command='dmrcc', anl0_only=False,
                              stationary_point='minimum'):
    """Build one L2 -> L3 -> parallel ANL validation graph.

    Every electronic, harmonic, anharmonic, DBOC, core-valence, and
    scalar-relativistic task belongs to this graph. Gaussian VPT2 and the
    L3 single point on L2 geometry branch directly from the accepted L2
    geometry. They therefore remain usable if L3 optimization is unavailable.

    The older interface/base/post builders remain available for reading and
    continuing immutable validation runs.  New end-to-end runs should use
    this graph so no external driver has to export L3 and prepare a second
    dispatcher directory.
    """
    interface = interface_validation_spec(
        molecule, max_nodes=max_nodes, partition=partition,
        stationary_point=stationary_point)
    post = post_geometry_validation_spec(
        molecule, max_nodes=max_nodes, partition=partition,
        mrcc_command=mrcc_command, anl0_only=anl0_only)

    # The interface probe carries a conventional CCSD(T)/DZ task.  The
    # higher-order group supplies the same physical component with the task
    # identity used by recipe assembly, so retain exactly one copy.
    tasks = [deepcopy(task) for task in interface['tasks']
             if task['id'] != 'ccsdt_dz']
    for task in post['tasks']:
        task = deepcopy(task)
        task['geometry_from'] = 'l3_geometry'
        task['geometry_fallback'] = 'l2_geometry'
        tasks.append(task)
    identifiers = [task['id'] for task in tasks]
    if len(identifiers) != len(set(identifiers)):
        raise RuntimeError('Composite validation task identifiers overlap.')
    by_id = {task['id']: task for task in tasks}
    if 'l3_geometry' in by_id:
        by_id['l3_geometry']['failure_skips'] = [
            task['id'] for task in tasks
            if (task.get('geometry_from') == 'l3_geometry'
                and 'geometry_fallback' not in task)]

    return {
        'schema': 1, 'name': 'anl-composite-validation',
        'molecule': molecule, 'limits': {'max_nodes': max_nodes},
        'tasks': tasks,
        'intent': {
            'claim': 'unified-anl-interface-validation-only',
            'requested_ladder_head': ('ANL0-F12' if anl0_only
                                      else 'ANL1-F12'),
            'anl0_only': anl0_only,
            'stationary_point': stationary_point,
            'higher_order_equation': post['intent']['higher_order_equation'],
            'backend_policy': post['intent']['backend_policy'],
            'selected_higher_order_tasks':
                post['intent']['selected_higher_order_tasks'],
            'rank_exact_higher_order':
                post['intent']['rank_exact_higher_order'],
            'core_valence_equation': post['intent']['core_valence_equation'],
            'scalar_relativistic_equation':
                post['intent']['scalar_relativistic_equation'],
            'scheduling': (
                'L2 geometry, then independent VPT2, L3-on-L2 single-point, '
                'and L3-optimization branches; L3 properties fan out only '
                'after the optimized geometry is accepted'),
        },
    }


def configure_task_convergence(task, *, scf_convergence, cc_convergence):
    """Update numerical convergence only; keep method, geometry and resources."""
    scf_exponent = convergence_exponent(scf_convergence, 'SCF convergence')
    cc_exponent = convergence_exponent(cc_convergence, 'CC convergence')
    scf_threshold = fortran_threshold(scf_convergence, 'SCF convergence')
    cc_threshold, cc_variance = molpro_cc_thresholds(cc_convergence)
    profile = task.get('profile')
    if isinstance(profile, dict):
        calculator = str(profile.get('calculator', '')).casefold()
        kwargs = dict(profile.get('calculator_kwargs', {}))
        if calculator == 'gaussian':
            kwargs['scf'] = f'xqc,conver={scf_exponent}'
        elif calculator == 'molpro':
            kwargs.update(
                scf_convergence=scf_convergence,
                correlation_convergence=cc_convergence)
        profile['calculator_kwargs'] = kwargs
    template = task.get('input_template')
    if not isinstance(template, str):
        return
    if task.get('backend') == 'gaussian':
        template = re.sub(
            r'\bSCF\s*=\s*(?:XQC|\([^\n)]*XQC[^\n)]*\))',
            f'SCF=(XQC,Conver={scf_exponent})', template,
            flags=re.IGNORECASE)
    elif task.get('backend') == 'molpro':
        template = re.sub(
            r'^\s*gthresh\s*,[^\n]*$',
            (f'gthresh,energy={scf_threshold},coeff={cc_threshold},'
             f'civec={cc_threshold}'), template,
            flags=re.IGNORECASE | re.MULTILINE)
        def cc_controls(match):
            body = match.group(1)
            # Migrate old generated cards too; preserve CORE, scale_trip and
            # all method options. Never apply CC semantics to RS2/MRCI.
            body = re.sub(r';\s*thresh\s*,\s*energy\s*=\s*1\.d-\d+\s*,\s*'
                          r'coeff\s*=\s*1\.d-\d+', '', body, flags=re.I)
            command, separator, cards = body.partition(';')
            command = re.sub(r',\s*(?:thrden|thrvar)\s*=\s*1\.d-\d+',
                             '', command, flags=re.I)
            return (f'{{{command},thrden={cc_threshold},thrvar={cc_variance}'
                    f'{separator}{cards}}}')
        template = re.sub(r'\{\s*([ru]?ccsd[^{}]*)\}', cc_controls,
                          template, flags=re.IGNORECASE)
    elif task.get('backend') == 'mrcc':
        template = re.sub(r'^sctol=\d+$', f'sctol={scf_exponent}',
                          template, flags=re.MULTILINE)
        template = re.sub(r'^cctol=\d+$', f'cctol={cc_exponent}',
                          template, flags=re.MULTILINE)
    elif task.get('backend') == 'cfour':
        template = re.sub(r'^SCF_CONV=\d+$',
                          f'SCF_CONV={scf_exponent}', template,
                          flags=re.MULTILINE)
        template = re.sub(r'^(?:CC_CONV|LINEQ_CONV)=\d+$',
                          lambda match: (match.group(0).split('=', 1)[0]
                                         + f'={cc_exponent}'), template,
                          flags=re.MULTILINE)
    task['input_template'] = template


def require_primary_geometries(spec, *, preserve_tasks=()):
    """Require requested L3 geometries; never silently stage an L2 fallback."""
    spec.setdefault('intent', {})['require_l3_geometry'] = True
    for task in spec['tasks']:
        if task['id'] in preserve_tasks:
            continue
        task.pop('geometry_fallback', None)
        if task['kind'] == 'ase_optimize':
            for key in ('failure_policy', 'failure_skips', 'component_families'):
                task.pop(key, None)


def automatic_composite_spec(molecule, *, highest_tier='ANL0-F12',
                             max_nodes=3, partition=None,
                             mrcc_command='dmrcc',
                             stationary_point='minimum',
                             vpt2_method='B2PLYP',
                             vpt2_basis='cc-pVTZ',
                             vpt2_dispersion='GD3BJ',
                             scf_convergence=1.e-8,
                             cc_convergence=1.e-7,
                             geometry_fmax=5.e-4,
                             minimum_recovery=None,
                             f12_augmented=False, require_l3_geometry=True,
                             l3_max_retries=1):
    """Build the automatic post-L2 graph for one stable point or tight TS.

    The caller supplies the already accepted L2 geometry and its provenance
    in ``molecule.source``. No L1 or L2 optimization is repeated. The graph
    contains only work that can contribute to the requested highest tier or
    one of its declared fallbacks.
    """
    from kinbot.minimum_recovery import validate_anl_recovery
    recovery = validate_anl_recovery(
        {} if minimum_recovery is None else minimum_recovery)
    ladders = {
        'ANL1-F12': ('ANL1-F12', 'ANL1', 'L3', 'L3-SP@L2', 'L2'),
        'ANL1': ('ANL1', 'L3', 'L3-SP@L2', 'L2'),
        'ANL0-F12': ('ANL0-F12', 'ANL0', 'L3', 'L3-SP@L2', 'L2'),
        'ANL0': ('ANL0', 'L3', 'L3-SP@L2', 'L2'),
        'L3': ('L3', 'L3-SP@L2', 'L2'),
        'L3-SP@L2': ('L3-SP@L2', 'L2'),
        'L2': ('L2',),
    }
    allowed = tuple(ladders)
    if highest_tier not in allowed:
        raise ValueError(f'highest_tier must be one of {allowed}.')
    if not isinstance(f12_augmented, bool):
        raise ValueError('f12_augmented must be boolean.')
    if not isinstance(require_l3_geometry, bool):
        raise ValueError('require_l3_geometry must be boolean.')
    if type(l3_max_retries) is not int or l3_max_retries < 0:
        raise ValueError('l3_max_retries must be a nonnegative integer.')
    scf_exponent = convergence_exponent(
        scf_convergence, 'SCF convergence')
    cc_exponent = convergence_exponent(
        cc_convergence, 'CC convergence')
    scf_threshold = fortran_threshold(scf_convergence, 'SCF convergence')
    cc_threshold, cc_coefficient = molpro_cc_thresholds(cc_convergence)
    if isinstance(geometry_fmax, bool) or not isinstance(geometry_fmax, (int, float)) or not math.isfinite(geometry_fmax) or geometry_fmax <= 0:
        raise ValueError('geometry_fmax must be a positive finite force threshold in eV/angstrom.')
    source = molecule.get('source') if isinstance(molecule, dict) else None
    if (not isinstance(source, dict)
            or not re.fullmatch(r'[0-9a-f]{64}', str(
                source.get('geometry_sha256', '')))
            or not re.fullmatch(r'[0-9a-f]{64}', str(
                source.get('record_sha256', '')))
            or not isinstance(source.get('energy_hartree'), (int, float))
            or isinstance(source.get('energy_hartree'), bool)
            or not math.isfinite(source['energy_hartree'])):
        raise ValueError(
            'Automatic composite input needs a hash-bound accepted L2 source.')
    anl1_requested = highest_tier in ('ANL1', 'ANL1-F12')
    full = composite_validation_spec(
        molecule, max_nodes=max_nodes, partition=partition,
        mrcc_command=mrcc_command, anl0_only=not anl1_requested,
        stationary_point=stationary_point)
    if (vpt2_method not in ('B3LYP', 'B2PLYP')
            or vpt2_basis.casefold() != 'cc-pvtz'
            or vpt2_dispersion not in ('', 'GD3BJ')
            or (vpt2_method == 'B2PLYP') != (vpt2_dispersion == 'GD3BJ')):
        raise ValueError('Unsupported automatic L2/VPT2 profile.')
    tasks = []
    for original in full['tasks']:
        if original['id'] == 'l2_geometry':
            continue
        task = deepcopy(original)
        if anl1_requested and task['id'] == 'l3_geometry':
            task['id'] = 'l3_geometry_tz'
        elif anl1_requested and task['id'] == 'harmonic':
            # ANL1 uses a separately generated basis-matched TZ/QZ pair.
            continue
        if task['id'] == 'gaussian_vpt2':
            dispersion = (f' EmpiricalDispersion={vpt2_dispersion}'
                          if vpt2_dispersion else '')
            task['input_template'] = (
                '%nprocshared={{CORES}}\n%mem={{WORK_MEMORY_MB}}MB\n'
                f'#p {vpt2_method}/{vpt2_basis} Freq=Anharmonic '
                f'Symmetry=None SCF=XQC{dispersion} Integral=UltraFine\n\n'
                'KinBot frequency-only VPT2 automatic calculation\n\n'
                '{{CHARGE}} {{MULT}}\n{{CARTESIAN}}\n\n')
            task['result_parser'] = {
                'kind': 'gaussian_vpt2', 'file': 'vpt2.log',
                'method': vpt2_method, 'basis': vpt2_basis,
                'dispersion': vpt2_dispersion}
        if task.get('geometry_from') == 'l2_geometry':
            task['geometry_from'] = 'initial'
        if task.get('geometry_fallback') == 'l2_geometry':
            task['geometry_fallback'] = 'initial'
        task['depends_on'] = [
            dependency for dependency in task.get('depends_on', ())
            if dependency != 'l2_geometry']
        if not task['depends_on']:
            task.pop('depends_on', None)
        tasks.append(task)
    if f12_augmented:
        for task in tasks:
            if task['id'] not in ('f12_tz', 'f12_qz'):
                continue
            old_basis = task['result_parser']['basis']
            new_basis = 'aug-' + old_basis
            task['input_template'] = task['input_template'].replace(
                f'basis={old_basis}', f'basis={new_basis}')
            task['result_parser']['basis'] = new_basis
    if anl1_requested:
        tz_geometry = next(task for task in tasks
                           if task['id'] == 'l3_geometry_tz')
        qz_geometry = deepcopy(tz_geometry)
        qz_geometry['id'] = 'l3_geometry'
        qz_geometry['geometry_from'] = 'l3_geometry_tz'
        qz_geometry['profile']['basis'] = 'cc-pVQZ'
        qz_geometry['profile']['label'] = 'CCSD(T)/cc-pVQZ'
        qz_geometry['resources'] = _resources(
            '7-00:00:00', max_cores=12, min_stack_mw=1024,
            partition=partition)
        tasks.append(qz_geometry)
        tasks.extend(_anl1_harmonic_tasks(
            molecule, partition=partition,
            stationary_point=stationary_point))
        tasks.extend(_anl1_reference_tasks(
            molecule, partition=partition,
            geometry_from='l3_geometry'))
        if highest_tier == 'ANL1-F12':
            tasks.extend(_anl1_f12_reference_tasks(
                molecule, partition=partition,
                geometry_from='l3_geometry',
                f12_augmented=f12_augmented))
            by_id = {task['id']: task for task in tasks}
            _component_fallback(
                by_id['f12_qz'], 'reference_cbs_dual_f12',
                failure_skips=('f12_5z',))
    if not anl1_requested:
        tasks.extend(_anl0_reference_tasks(
            molecule, partition=partition, geometry_from='l3_geometry'))
        if highest_tier == 'ANL0-F12':
            by_id = {task['id']: task for task in tasks}
            _component_fallback(
                by_id['f12_tz'], 'reference_cbs_dual_f12',
                failure_skips=('f12_qz',))
    for task in tasks:
        if (task.get('geometry_from') == 'l3_geometry'
                and task['id'] != 'harmonic'):
            task['geometry_fallback'] = 'initial'

    keep = {
        'ANL1-F12': ({task['id'] for task in tasks}
                     - {'f12_tz', 'anl0_qz', 'anl0_5z'}),
        'ANL1': ({task['id'] for task in tasks}
                 - {'f12_tz', 'f12_qz', 'f12_5z',
                    'anl0_qz', 'anl0_5z'}),
        'ANL0-F12': ({task['id'] for task in tasks}
                     - {'f12_5z', 'anl1_5z', 'anl1_6z'}),
        'ANL0': ({task['id'] for task in tasks}
                 - {'f12_tz', 'f12_qz'}),
        'L3': {'l3_geometry', 'harmonic', 'gaussian_vpt2', 'l3_sp_l2'},
        'L3-SP@L2': {'gaussian_vpt2', 'l3_sp_l2'},
        'L2': {'gaussian_vpt2'},
    }[highest_tier]
    if stationary_point == 'transition_state' and highest_tier == 'L2':
        # The accepted KinBot L2 TS record already contains its harmonic ZPE.
        keep = set()
    tasks = [task for task in tasks if task['id'] in keep]
    identifiers = {task['id'] for task in tasks}
    if 'l3_sp_l2' in identifiers and 'l3_geometry' in identifiers:
        by_id = {task['id']: task for task in tasks}
        by_id['l3_sp_l2']['fallback_on'] = ['l3_geometry']
    for task in tasks:
        task['depends_on'] = [dependency for dependency in
                              task.get('depends_on', ())
                              if dependency in identifiers]
        if not task['depends_on']:
            task.pop('depends_on', None)
        task['failure_skips'] = [target for target in
                                 task.get('failure_skips', ())
                                 if target in identifiers]
    by_id = {task['id']: task for task in tasks}
    if 'l3_geometry' in by_id:
        by_id['l3_geometry']['failure_skips'] = [
            task['id'] for task in tasks
            if (task.get('geometry_from') == 'l3_geometry'
                and 'geometry_fallback' not in task)]
    for task in tasks:
        if task.get('kind') == 'ase_optimize':
            task['optimizer']['fmax'] = geometry_fmax
            if (not task['optimizer'].get('order', 0)
                    and not task['optimizer'].get('constraints')):
                task['optimizer']['minimum_recovery'] = dict(recovery)
        configure_task_convergence(
            task, scf_convergence=scf_convergence,
            cc_convergence=cc_convergence)
    spec = {
        'schema': 1, 'name': 'anl-automatic-post-l2',
        'molecule': molecule, 'limits': {'max_nodes': max_nodes},
        'tasks': tasks,
        'intent': {
            'convergence': {
                'scf': scf_convergence,
                'coupled_cluster_or_ci': cc_convergence,
                'geometry_fmax_ev_angstrom': geometry_fmax,
            },
            'f12_augmented': f12_augmented,
            'l3_max_retries': l3_max_retries,
            'claim': 'automatic-cbh-anl-production',
            # A user may deliberately choose a short partition. Native jobs
            # then request that partition's limit, so scheduler timeouts are
            # observable execution failures and optional recipe components
            # can follow the normal recorded fallback ladder.
            'cap_walltime_to_explicit_partition': True,
            'highest_tier': highest_tier,
            'fallback_ladder': list(ladders[highest_tier]),
            'stationary_point': stationary_point,
            'source_geometry_role': 'accepted-l2',
            'anl0_only': not anl1_requested,
        },
    }

    if require_l3_geometry:
        require_primary_geometries(spec)
    return spec


def audit_higher_order_run(run_dir):
    """Reparse a complete full or targeted higher-order interface probe."""
    _, spec, state = _load(run_dir)
    if spec.get('name') not in ('anl1-higher-order-interface-validation',
                                 'anl-post-geometry-validation',
                                 'anl-composite-validation'):
        raise ValueError('Run is not an ANL1 higher-order validation graph.')
    declared = {task['id'] for task in spec['tasks']}
    identifiers = tuple(ident for ident in _HIGHER_ORDER_TASK_IDS
                        if ident in declared)
    rank_exact = _rank_exact_higher_order(spec['molecule'])
    selected = tuple(spec.get('intent', {}).get(
        'selected_tasks', identifiers))
    if not identifiers and rank_exact:
        corrections = {}
        if set(selected) & {'ccsdt_dz', 'ccsdtq_dz'}:
            corrections['delta_q_dz'] = 0.
        if set(selected) & {'ccsdt_tz', 'ccsdtq_tz'}:
            corrections['delta_q_tz'] = 0.
        if set(selected) & {'ccsdtq_dz', 'ccsdtqp_dz'}:
            corrections['delta_p_dz'] = 0.
        if {'delta_q_tz', 'delta_p_dz'} <= corrections.keys():
            corrections['anl1_higher_order'] = 0.
        return {
            'status': 'rank_exact_higher_order_not_scheduled',
            'task_statuses': {}, 'corrections_hartree': corrections,
            'components': {}, 'claim': 'interface-validation-only',
            'rank_exact_higher_order': {
                'electron_count': _electron_count(spec['molecule']),
                'correction_hartree': 0.,
            },
        }
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
    if rank_exact:
        anl0_requirements = {item.key: item for item in recipe(
            'ANL0-F12', vpt2_method='B2PLYP-D3BJ',
            multiplicity=spec['molecule'].get('multiplicity', 1)
        ).requirements}
        anl1_requirements = {item.key: item for item in recipe(
            'ANL1', vpt2_method='B2PLYP-D3BJ',
            multiplicity=spec['molecule'].get('multiplicity', 1)
        ).requirements}
        if 'ccsdtq_tz' in selected:
            components['ccsdtq_tz'] = rank_exact_higher_order_component(
                spec['molecule'], components['ccsdt_tz'],
                anl1_requirements['hoe_tz_high'], key='ccsdtq_tz')
        if set(selected) & {'ccsdtq_dz', 'ccsdtqp_dz'}:
            components['ccsdtq_dz'] = rank_exact_higher_order_component(
                spec['molecule'], components['ccsdt_dz'],
                anl0_requirements['hoe_high'], key='ccsdtq_dz')
        if 'ccsdtqp_dz' in selected:
            components['ccsdtqp_dz'] = rank_exact_higher_order_component(
                spec['molecule'], components['ccsdtq_dz'],
                anl1_requirements['hoe_dz_high'], key='ccsdtqp_dz')
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
                   if set(selected) == set(_HIGHER_ORDER_TASK_IDS)
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
    if set(selected) == set(_HIGHER_ORDER_TASK_IDS):
        result['correction_hartree'] = corrections['anl1_higher_order']
    return result


def audit_anl0_post_geometry_run(run_dir):
    """Audit only the terms required by an ANL0 or ANL0-F12 assembly.

    A combined validation graph may also contain ANL1-only CCSDT(Q)/TZ and
    CCSDTQ(P)/DZ probes.  Their failure or deliberate cancellation must not
    prevent an independently complete ANL0-F12 result from being audited.
    """
    run_dir, spec, state = _load(run_dir)
    if spec.get('name') not in ('anl-post-geometry-validation',
                                 'anl-composite-validation'):
        raise ValueError('Run is not a combined post-geometry validation graph.')
    rank_exact = _rank_exact_higher_order(spec['molecule'])
    core_valence_exact = _rank_exact_core_valence(spec['molecule'])
    required = (() if rank_exact else ('ccsdt_dz', 'ccsdtq_dz')) + (
                    (() if core_valence_exact else
                     ('cv_ae_tz', 'cv_ae_qz', 'cv_fc_tz', 'cv_fc_qz'))
                    + ('rel_dkh', 'rel_nonrel'))
    declared = {task['id'] for task in spec['tasks']}
    missing = sorted(set(required) - declared)
    if missing:
        raise RuntimeError(f'ANL0 post-geometry tasks are absent: {missing}')
    statuses = {ident: state['tasks'].get(ident, {}).get('status', 'waiting')
                for ident in required}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'ANL0 post-geometry run is incomplete: {statuses}')
    state_id = 'anl0-post-geometry-validation-state'
    if rank_exact:
        requirement = next(
            item for item in recipe(
                'ANL0-F12', vpt2_method='B2PLYP-D3BJ',
                multiplicity=spec['molecule'].get('multiplicity', 1),
                rank_exact_higher_order=True
            ).requirements if item.key == 'higher_order_correction')
        geometry_hash = (
            spec.get('molecule', {}).get('source', {}).get('geometry_sha256')
            or state['tasks'].get('l3_geometry', {}).get(
                'final_geometry_sha256'))
        correction = rank_exact_higher_order_correction(
            spec['molecule'], requirement, state_id=state_id,
            geometry_sha256=geometry_hash)
        higher_components = {'higher_order_correction': asdict(correction)}
        higher_value = correction.value_hartree
    else:
        low = task_component(
            run_dir, 'ccsdt_dz', key='hoe_low', state_id=state_id)
        high = task_component(
            run_dir, 'ccsdtq_dz', key='hoe_high', state_id=state_id)
        higher_components = {
            'ccsdtq_dz': asdict(high), 'ccsdt_dz': asdict(low)}
        higher_value = math.fsum((high.value_hartree, -low.value_hartree))
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
        'higher_order_dz_hartree': higher_value,
        'higher_order_components': higher_components,
        'rank_exact_higher_order': rank_exact,
        'common_corrections': corrections,
        'geometry_sha256': (
            spec.get('molecule', {}).get('source', {}).get('geometry_sha256')
            or state['tasks'].get('l3_geometry', {}).get(
                'final_geometry_sha256')),
    }


def audit_common_corrections_run(run_dir):
    """Reparse and combine completed core-valence and DKH correction jobs."""
    _, spec, state = _load(run_dir)
    if spec.get('name') not in ('anl-common-corrections-validation',
                                 'anl-post-geometry-validation',
                                 'anl-composite-validation'):
        raise ValueError('Run is not an ANL common-corrections graph.')
    molecule = spec['molecule']
    core_valence_exact = _rank_exact_core_valence(molecule)
    identifiers = ((() if core_valence_exact else
                    ('cv_ae_tz', 'cv_ae_qz', 'cv_fc_tz', 'cv_fc_qz'))
                   + ('rel_dkh', 'rel_nonrel'))
    statuses = {ident: state['tasks'].get(ident, {}).get('status', 'waiting')
                for ident in identifiers}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'Common-corrections run is incomplete: {statuses}')
    equation = recipe(
        'ANL0-F12', vpt2_method='B2PLYP-D3BJ',
        multiplicity=spec['molecule'].get('multiplicity', 1))
    requirements = {item.key: item for item in equation.requirements}
    state_id = 'common-corrections-validation-state'
    relativistic = scalar_relativistic_task_component(
        run_dir, 'rel_dkh', 'rel_nonrel',
        requirement=requirements['scalar_relativistic'], state_id=state_id)
    core_valence = (rank_exact_core_valence_component(
        molecule, requirements['core_valence_cbs'], state_id=state_id,
        geometry_sha256=relativistic.geometry_sha256)
        if core_valence_exact else core_valence_task_component(
            run_dir, all_electron_lower='cv_ae_tz',
            all_electron_upper='cv_ae_qz', frozen_core_lower='cv_fc_tz',
            frozen_core_upper='cv_fc_qz',
            requirement=requirements['core_valence_cbs'], state_id=state_id))
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
    if spec.get('name') not in ('anl-post-geometry-validation',
                                 'anl-composite-validation'):
        raise ValueError('Run is not a combined post-geometry validation graph.')
    post_ids = set(_HIGHER_ORDER_TASK_IDS) | {
        'cv_ae_tz', 'cv_ae_qz', 'cv_fc_tz', 'cv_fc_qz',
        'rel_dkh', 'rel_nonrel'}
    statuses = {
        task['id']: state['tasks'].get(task['id'], {}).get('status', 'waiting')
        for task in spec['tasks'] if task['id'] in post_ids}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'Post-geometry run is incomplete: {statuses}')
    source = spec.get('molecule', {}).get('source', {})
    higher = audit_higher_order_run(run_dir)
    corrections = audit_common_corrections_run(run_dir)
    return {
        'status': 'post_geometry_interface_complete',
        'task_statuses': statuses,
        'geometry_sha256': (
            source.get('geometry_sha256')
            or state['tasks'].get('l3_geometry', {}).get(
                'final_geometry_sha256')),
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


def molecule_from_smiles(smiles, *, charge, multiplicity, random_seed=312):
    """Build a deterministic initial Cartesian structure for an ANL graph.

    This is an initial geometry only.  The prepared graph still performs the
    requested L2 and L3 optimizations before using any geometry-dependent
    energy or frequency result.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError('SMILES must be a nonempty string.')
    if (isinstance(charge, bool) or not isinstance(charge, int)
            or isinstance(multiplicity, bool)
            or not isinstance(multiplicity, int) or multiplicity < 1):
        raise ValueError('SMILES charge and multiplicity must be valid integers.')
    if (isinstance(random_seed, bool) or not isinstance(random_seed, int)
            or random_seed < 0 or random_seed > 0x7fffffff):
        raise ValueError('RDKit random seed must be a nonnegative 32-bit integer.')
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f'RDKit could not parse SMILES {smiles!r}.')
    molecule = Chem.AddHs(molecule)
    if Chem.GetFormalCharge(molecule) != charge:
        raise ValueError('SMILES formal charge disagrees with the requested charge.')
    radical_electrons = sum(atom.GetNumRadicalElectrons()
                            for atom in molecule.GetAtoms())
    unpaired = multiplicity - 1
    if radical_electrons < unpaired or (radical_electrons - unpaired) % 2:
        raise ValueError('SMILES radical electrons disagree with multiplicity.')

    if molecule.GetNumAtoms() == 1:
        conformer = Chem.Conformer(1)
        conformer.SetAtomPosition(0, (0., 0., 0.))
        molecule.AddConformer(conformer, assignId=True)
    else:
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = random_seed
        status = AllChem.EmbedMolecule(molecule, parameters)
        if status != 0:
            molecule.RemoveAllConformers()
            parameters.useRandomCoords = True
            status = AllChem.EmbedMolecule(molecule, parameters)
        if status != 0:
            raise ValueError(f'RDKit could not embed SMILES {smiles!r}.')
        if molecule.GetNumAtoms() > 2:
            if AllChem.MMFFHasAllMoleculeParams(molecule):
                AllChem.MMFFOptimizeMolecule(molecule, maxIters=2000)
            elif AllChem.UFFHasAllMoleculeParams(molecule):
                AllChem.UFFOptimizeMolecule(molecule, maxIters=2000)

    conformer = molecule.GetConformer()
    positions = []
    for index in range(molecule.GetNumAtoms()):
        point = conformer.GetAtomPosition(index)
        positions.append([float(point.x), float(point.y), float(point.z)])
    result = {
        'symbols': [atom.GetSymbol() for atom in molecule.GetAtoms()],
        'positions': positions, 'charge': charge,
        'multiplicity': multiplicity, 'smiles': smiles,
        'initial_geometry': {
            'provider': ('RDKit-ETKDGv3' if molecule.GetNumAtoms() > 1
                         else 'RDKit-origin'),
            'random_seed': random_seed,
        },
    }
    # Reuse the dispatcher's electron-count and Cartesian validation before a
    # workflow can be written or submitted.
    from kinbot.anl.dispatch import _atoms
    _atoms(result)
    return result


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


def _dual_f12_reference_component(run_dir, requirement, *, state_id):
    """Assemble and average the F12 and fully augmented CBS limits."""
    settings = dict(requirement.settings)
    geometry = {
        key: settings[key] for key in ('geometry_method', 'geometry_basis')
        if key in settings
    }
    f12_pair = tuple(settings['f12_native_basis_pair'])
    conventional_pair = tuple(settings['conventional_native_basis_pair'])
    if settings['f12_upper_cardinal'] == 4:
        f12_tasks = ('f12_tz', 'f12_qz')
    elif settings['f12_upper_cardinal'] == 5:
        f12_tasks = ('f12_qz', 'f12_5z')
    else:
        raise ValueError('Unsupported dual-CBS F12 cardinal pair.')
    if settings['conventional_upper_cardinal'] == 5:
        conventional_tasks = ('anl0_qz', 'anl0_5z')
    elif settings['conventional_upper_cardinal'] == 6:
        conventional_tasks = ('anl1_5z', 'anl1_6z')
    else:
        raise ValueError('Unsupported dual-CBS conventional cardinal pair.')
    f12_requirement = ComponentRequirement(
        'reference_cbs_f12', 'electronic', 'CCSD(T)-F12b',
        f'CBS({f12_pair[0]},{f12_pair[1]})', 'l3', ('composite',),
        {'scale_trip': settings['scale_trip'],
         'correlation': settings['correlation'],
         'reference': settings['reference'],
         'reference_convention': settings['reference_convention'],
         'upper_cardinal': settings['f12_upper_cardinal'],
         'extrapolation_power': settings['f12_extrapolation_power'],
         'native_basis_pair': f12_pair,
         **geometry})
    conventional_requirement = ComponentRequirement(
        'reference_cbs_conventional', 'electronic', 'CCSD(T)',
        f'CBS({conventional_pair[0]},{conventional_pair[1]})',
        'l3', ('composite',),
        {'correlation': settings['correlation'],
         'reference': settings['reference'],
         'reference_convention': settings['reference_convention'],
         'upper_cardinal': settings['conventional_upper_cardinal'],
         'extrapolation_power':
             settings['conventional_extrapolation_power'],
         'native_basis_pair':
             settings['conventional_native_basis_pair'],
         **geometry})
    f12 = cbs_task_component(
        run_dir, *f12_tasks, requirement=f12_requirement,
        state_id=state_id,
        lower_basis=f12_pair[0], upper_basis=f12_pair[1])
    conventional = cbs_task_component(
        run_dir, *conventional_tasks,
        requirement=conventional_requirement, state_id=state_id,
        lower_basis=conventional_pair[0], upper_basis=conventional_pair[1])
    if (f12.state_id, f12.charge, f12.multiplicity,
            f12.geometry_sha256) != (
            conventional.state_id, conventional.charge,
            conventional.multiplicity, conventional.geometry_sha256):
        raise ValueError('Dual-CBS F12 reference estimates describe different states.')
    weights = settings.get('component_weights')
    if (not isinstance(weights, (tuple, list)) or len(weights) != 2
            or any(isinstance(value, bool)
                   or not isinstance(value, (int, float))
                   or not math.isfinite(value) or value < 0.
                   for value in weights)
            or abs(math.fsum(weights) - 1.) > 1.e-12):
        raise ValueError('Dual-CBS reference weights must be two nonnegative values summing to one.')
    f12_weight, conventional_weight = map(float, weights)
    uncertainty_series = {}
    for prefix, component, weight in (
            ('f12', f12, f12_weight),
            ('conventional', conventional, conventional_weight)):
        series = component.settings.get('basis_uncertainty_series', {})
        if len(series) != 1:
            raise ValueError(
                f'Dual-CBS {prefix} reference lacks native endpoint data.')
        endpoint = next(iter(series.values()))
        uncertainty_series[f'{requirement.key}:{prefix}'] = {
            'lower_hartree': weight * float(endpoint['lower_hartree']),
            'upper_hartree': weight * float(endpoint['upper_hartree']),
            'lower_basis': endpoint['lower_basis'],
            'upper_basis': endpoint['upper_basis'],
            'weight': weight,
        }
    result_settings = dict(settings)
    result_settings['basis_uncertainty_series'] = uncertainty_series
    provenance = {
        'formula': (
            f'{f12_weight:g}*CBS(F12,{f12_pair[0]}/{f12_pair[1]})+'
            f'{conventional_weight:g}*CBS(conventional,'
            f'{conventional_pair[0]}/{conventional_pair[1]})'),
        'f12_sha256': f12.source_sha256,
        'conventional_sha256': conventional.source_sha256,
        'settings': settings,
    }
    digest = hashlib.sha256(json.dumps(
        provenance, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return ComponentResult(
        key=requirement.key,
        value_hartree=math.fsum(
            (f12_weight * f12.value_hartree,
             conventional_weight * conventional.value_hartree)),
        quantity=requirement.quantity, method=requirement.method,
        basis=requirement.basis, backend='composite',
        state_id=f12.state_id, charge=f12.charge,
        multiplicity=f12.multiplicity,
        geometry_sha256=f12.geometry_sha256,
        source_sha256=digest,
        source=f'mean({f12.source},{conventional.source})',
        settings=result_settings)


def _assemble_profiled_anl(
        recipe_name, interface_run, higher_order_run, corrections_run, *, state_id,
        spin_orbit_hartree, spin_orbit_source, spin_orbit_backend='manual',
        vpt2_review=None, base_run=None, harmonic_only=False,
        vpt2_method='B2PLYP-D3BJ', l2_geometry_fallback=None,
        omitted_families=(), f12_augmented=None):
    """Assemble one provenance-checked supported ANL 0 K energy.

    Every geometry-dependent term must trace to the accepted L2 geometry in
    ``interface_run`` and the accepted current L3 geometry in ``base_run``.
    For a fully current interface graph, omitting ``base_run`` uses that graph
    for both roles. Higher-order and correction graphs must identify the same
    completed L3 result.
    """
    if recipe_name not in ('ANL0', 'ANL0-F12', 'ANL1', 'ANL1-F12'):
        raise ValueError('Unsupported profiled ANL assembly recipe.')
    anl1 = recipe_name in ('ANL1', 'ANL1-F12')
    interface_run, interface_spec, interface_state = _load(interface_run)
    if interface_spec.get('name') not in (
            'anl1-f12-non-mrcc-interface-validation',
            'anl-composite-validation', 'anl-automatic-post-l2'):
        raise ValueError('Interface run has the wrong workflow type.')
    if not isinstance(state_id, str) or not state_id.strip():
        raise ValueError('state_id must be a nonempty string.')
    molecule = interface_spec['molecule']
    charge = molecule.get('charge', 0)
    multiplicity = molecule.get('multiplicity', 1)
    if interface_spec.get('name') == 'anl-automatic-post-l2':
        l2 = molecule
        l2_hash = molecule.get('source', {}).get('geometry_sha256')
        if not re.fullmatch(r'[0-9a-f]{64}', str(l2_hash)):
            raise ValueError('Automatic run lacks its accepted L2 geometry hash.')
    else:
        l2 = molecule_from_completed_run(interface_run, 'l2_geometry')
        l2_hash = l2['source']['geometry_sha256']
    requested_base = (None if base_run is None else Path(base_run).resolve())
    if requested_base is None or requested_base == interface_run:
        # Supplying the interface explicitly is equivalent to omitting the
        # optional base. This matters for completed full-interface graphs,
        # whose workflow type is valid for both roles but is intentionally
        # distinct from a split current-base graph.
        base_run = interface_run
        base_spec, base_state = interface_spec, interface_state
    else:
        base_run, base_spec, base_state = _load(requested_base)
        if base_spec.get('name') not in ('anl-current-base-validation',
                                         'anl-composite-validation',
                                         'anl-automatic-post-l2'):
            raise ValueError('Base run has the wrong workflow type.')
        if interface_spec.get('name') == 'anl-automatic-post-l2':
            if (base_spec.get('molecule', {}).get('source', {}).get(
                    'geometry_sha256') != l2_hash):
                raise ValueError('Automatic base run uses a different L2 geometry.')
        else:
            _same_imported_geometry(
                interface_run, base_run, l2_hash, 'L2', 'l2_geometry')
    recorded_f12_augmented = bool(
        base_spec.get('intent', {}).get('f12_augmented', False))
    if f12_augmented is None:
        f12_augmented = recorded_f12_augmented
    elif (not isinstance(f12_augmented, bool)
          or f12_augmented != recorded_f12_augmented):
        raise ValueError(
            'Requested F12 augmentation differs from the workflow intent.')
    l3_status = base_state['tasks'].get('l3_geometry', {}).get('status')
    geometry_fallback = l3_status in ('unavailable', 'skipped')
    harmonic_ids = (('harmonic_tz', 'harmonic_qz') if anl1
                    else ('harmonic',))
    harmonic_statuses = tuple(
        base_state['tasks'].get(task_id, {}).get('status')
        for task_id in harmonic_ids)
    harmonic_fallback = (len(molecule.get('symbols', ())) > 1
                         and any(status in ('unavailable', 'skipped')
                                 for status in harmonic_statuses))
    if geometry_fallback:
        harmonic_fallback = True
    needs_l2_fallback = geometry_fallback or harmonic_fallback
    if (needs_l2_fallback
            and (interface_spec.get('name') != 'anl-automatic-post-l2'
                 or not isinstance(l2_geometry_fallback, ComponentResult)
                 or l2_geometry_fallback.geometry_sha256 != l2_hash)):
        raise IncompleteRecipeError(
            'An accepted-L2 geometry/frequency fallback is required.')
    if geometry_fallback:
        if (Path(higher_order_run).resolve() != base_run
                or Path(corrections_run).resolve() != base_run):
            raise IncompleteRecipeError(
                'L3 geometry is unavailable and no complete accepted-L2 '
                'geometry/frequency fallback was supplied.')
        l3_hash = l2_hash
        l2_theory = molecule.get('source', {}).get('l2_theory', {})
        fallback_geometry_method = l2_theory.get('method')
        fallback_basis = l2_theory.get('basis')
        if not fallback_geometry_method or not fallback_basis:
            raise ValueError('Accepted L2 fallback provenance is incomplete.')
    else:
        l3 = molecule_from_completed_run(base_run, 'l3_geometry')
        l3_hash = l3['source']['geometry_sha256']
        fallback_geometry_method = fallback_basis = None
        if Path(higher_order_run).resolve() != base_run:
            _same_imported_l3_geometry(base_run, higher_order_run, l3_hash)
        if Path(corrections_run).resolve() != base_run:
            _same_imported_l3_geometry(base_run, corrections_run, l3_hash)

    atomic = len(molecule.get('symbols', ())) == 1
    rank_exact = _rank_exact_higher_order(molecule)
    equation = recipe(
        recipe_name, vpt2_method=vpt2_method,
        multiplicity=multiplicity, atomic=atomic,
        f12_augmented=f12_augmented,
        rank_exact_higher_order=rank_exact,
        harmonic_only=harmonic_only,
        omitted_families=tuple(omitted_families),
        fallback_geometry_method=fallback_geometry_method,
        fallback_harmonic_method=(l2_geometry_fallback.method
                                  if harmonic_fallback else None),
        fallback_basis=((l2_geometry_fallback.basis
                         if harmonic_fallback else fallback_basis)
                        or 'cc-pVTZ'))
    required = {item.key: item for item in equation.requirements}
    if recipe_name in ('ANL0-F12', 'ANL1-F12'):
        reference_tasks = reference_bases = None
    elif recipe_name == 'ANL1':
        reference_tasks = ('anl1_5z', 'anl1_6z')
        reference_bases = ANL1_NATIVE_REFERENCE_BASES
    else:
        reference_tasks = ('anl0_qz', 'anl0_5z')
        reference_bases = ANL0_NATIVE_REFERENCE_BASES
    if recipe_name in ('ANL0-F12', 'ANL1-F12'):
        reference = _dual_f12_reference_component(
            base_run, required['reference_cbs'], state_id=state_id)
    else:
        reference = cbs_task_component(
            base_run, *reference_tasks,
            requirement=required['reference_cbs'], state_id=state_id,
            lower_basis=reference_bases[0], upper_basis=reference_bases[1])
    components = {'reference_cbs': reference}
    if 'dboc' in required:
        components['dboc'] = task_component(
            base_run, 'cfour_dboc', key='dboc', state_id=state_id)
    if atomic:
        components['harmonic_zpe'] = atomic_zero_vibrational_component(
            molecule, required['harmonic_zpe'], state_id=state_id,
            geometry_sha256=l3_hash)
        if not harmonic_only:
            components['vpt2_correction'] = atomic_zero_vibrational_component(
                molecule, required['vpt2_correction'], state_id=state_id,
                geometry_sha256=l2_hash)
    elif harmonic_fallback:
        components['harmonic_zpe'] = l2_geometry_fallback
        if not harmonic_only:
            components['vpt2_correction'] = task_component(
                interface_run, 'gaussian_vpt2', key='vpt2_correction',
                state_id=state_id)
    elif anl1:
        components['harmonic_zpe'] = cbs_task_component(
            base_run, 'harmonic_tz', 'harmonic_qz',
            requirement=required['harmonic_zpe'], state_id=state_id,
            lower_basis='cc-pVTZ', upper_basis='cc-pVQZ')
        if not harmonic_only:
            components['vpt2_correction'] = task_component(
                interface_run, 'gaussian_vpt2', key='vpt2_correction',
                state_id=state_id)
    else:
        components['harmonic_zpe'] = task_component(
            base_run, 'harmonic', key='harmonic_zpe', state_id=state_id)
        if not harmonic_only:
            components['vpt2_correction'] = task_component(
                interface_run, 'gaussian_vpt2', key='vpt2_correction',
                state_id=state_id)
    if rank_exact and 'higher_order_correction' in required:
        components['higher_order_correction'] = \
            rank_exact_higher_order_correction(
                molecule, required['higher_order_correction'],
                state_id=state_id, geometry_sha256=l3_hash)
    elif 'hoe_low' in required:
        components['hoe_low'] = task_component(
            higher_order_run, 'ccsdt_dz', key='hoe_low', state_id=state_id)
        components['hoe_high'] = task_component(
            higher_order_run, 'ccsdtq_dz', key='hoe_high',
            state_id=state_id)
    elif 'hoe_tz_low' in required:
        for key, task_id in (
                ('hoe_tz_low', 'ccsdt_tz'),
                ('hoe_tz_high', 'ccsdtq_tz'),
                ('hoe_dz_low', 'ccsdtq_dz'),
                ('hoe_dz_high', 'ccsdtqp_dz')):
            components[key] = task_component(
                higher_order_run, task_id, key=key, state_id=state_id)
    if harmonic_only:
        if vpt2_review is not None:
            raise ValueError(
                'A harmonic-only energy must not receive a VPT2 review.')
    elif atomic:
        if vpt2_review is not None:
            raise ValueError('A monatomic species must not receive a VPT2 review.')
    else:
        components['vpt2_correction'] = _apply_quality_review(
            components['vpt2_correction'], 'gaussian_vpt2', vpt2_review)
    if 'core_valence_cbs' in required:
        components['core_valence_cbs'] = (
            rank_exact_core_valence_component(
                molecule, required['core_valence_cbs'], state_id=state_id,
                geometry_sha256=l3_hash)
            if _rank_exact_core_valence(molecule) else
            core_valence_task_component(
                corrections_run, all_electron_lower='cv_ae_tz',
                all_electron_upper='cv_ae_qz', frozen_core_lower='cv_fc_tz',
                frozen_core_upper='cv_fc_qz',
                requirement=required['core_valence_cbs'], state_id=state_id))
    if 'scalar_relativistic' in required:
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
    record = {
        'schema': 1, 'status': 'complete', 'state_id': state_id,
        'charge': charge, 'multiplicity': multiplicity,
        'recipe': result.recipe, 'recipe_version': result.recipe_version,
        'electronic_hartree': result.electronic_hartree,
        'zero_point_hartree': result.zero_point_hartree,
        'zero_k_hartree': result.zero_k_hartree,
        'geometry_sha256': {'l2': l2_hash, 'l3': l3_hash},
        'geometry_provenance': {
            'electronic_single_points': (
                'accepted_l2_fallback' if geometry_fallback else
                'accepted_l3_optimization'),
            'harmonic_zpe': (
                'accepted_l2_fallback' if harmonic_fallback else
                ('l3_harmonic_cbs' if anl1 else 'l3_harmonic')),
        },
        'components': {key: asdict(value)
                       for key, value in result.components.items()},
        'electronic_terms': [asdict(term) for term in result.electronic_terms],
        'zero_point_terms': [asdict(term) for term in result.zero_point_terms],
    }
    if omitted_families:
        record['reduced_recipe'] = {
            'parent_recipe': recipe(
                recipe_name, vpt2_method=vpt2_method,
                multiplicity=multiplicity, atomic=atomic,
                f12_augmented=f12_augmented,
                rank_exact_higher_order=rank_exact,
                harmonic_only=harmonic_only,
                fallback_geometry_method=fallback_geometry_method,
                fallback_harmonic_method=(
                    l2_geometry_fallback.method if harmonic_fallback else None),
                fallback_basis=((l2_geometry_fallback.basis
                                 if harmonic_fallback else fallback_basis)
                                or 'cc-pVTZ')).name,
            'omitted_families': list(omitted_families),
        }
    return record


def assemble_profiled_anl0_f12(
        interface_run, higher_order_run, corrections_run, **kwargs):
    """Assemble the explicitly correlated profiled ANL0 variant."""
    return _assemble_profiled_anl(
        'ANL0-F12', interface_run, higher_order_run, corrections_run,
        **kwargs)


def assemble_profiled_anl0(
        interface_run, higher_order_run, corrections_run, **kwargs):
    """Assemble profiled ANL0 with the fully augmented QZ/5Z reference."""
    return _assemble_profiled_anl(
        'ANL0', interface_run, higher_order_run, corrections_run, **kwargs)


def assemble_profiled_anl1_f12(
        interface_run, higher_order_run, corrections_run, **kwargs):
    """Assemble the KinBot ANL1-F12 profiled extension."""
    return _assemble_profiled_anl(
        'ANL1-F12', interface_run, higher_order_run, corrections_run,
        **kwargs)


def assemble_profiled_anl1(
        interface_run, higher_order_run, corrections_run, **kwargs):
    """Assemble profiled ANL1 with the fully augmented 5Z/6Z reference."""
    return _assemble_profiled_anl(
        'ANL1', interface_run, higher_order_run, corrections_run, **kwargs)


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
    diagnostic_requirement = ComponentRequirement(
        'reference_cbs', 'electronic', 'CCSD(T)-F12b',
        'CBS(cc-pVTZ-F12,cc-pVQZ-F12)', 'l3', ('composite',),
        {'scale_trip': 1, 'correlation': 'unrestricted',
         'reference': requirement.settings['reference'],
         'reference_convention': requirement.settings['reference_convention'],
         'upper_cardinal': 4,
         'extrapolation_power':
             requirement.settings['f12_extrapolation_power'],
         'geometry_method': requirement.settings['geometry_method'],
         'geometry_basis': requirement.settings['geometry_basis']})
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
            upper_cardinal=4,
            power=requirement.settings['f12_extrapolation_power'])
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
        run_dir, 'f12_tz', 'f12_qz', requirement=diagnostic_requirement,
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
    if spec.get('name') not in ('anl-current-base-validation',
                                 'anl-composite-validation'):
        raise ValueError('Run is not a current ANL base validation graph.')
    molecule = spec['molecule']
    atomic = len(molecule.get('symbols', ())) == 1
    required_ids = {'l3_geometry', 'f12_tz', 'f12_qz', 'cfour_dboc'}
    if not atomic:
        required_ids.add('harmonic')
    declared = {task['id'] for task in spec['tasks']}
    missing = sorted(required_ids - declared)
    if missing:
        raise RuntimeError(f'Current base tasks are absent: {missing}')
    statuses = {
        ident: state['tasks'].get(ident, {}).get('status', 'waiting')
        for ident in sorted(required_ids)}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'Current base run is incomplete: {statuses}')
    equation = recipe(
        'ANL0-F12', vpt2_method='B2PLYP-D3BJ',
        multiplicity=molecule.get('multiplicity', 1), atomic=atomic)
    required = {item.key: item for item in equation.requirements}
    reference = cbs_task_component(
        run_dir, 'f12_tz', 'f12_qz',
        requirement=required['reference_cbs'], state_id='base-validation-state',
        lower_basis='cc-pVTZ-F12', upper_basis='cc-pVQZ-F12')
    harmonic = (atomic_zero_vibrational_component(
        molecule, required['harmonic_zpe'], state_id='base-validation-state',
        geometry_sha256=reference.geometry_sha256) if atomic else
        task_component(run_dir, 'harmonic', key='harmonic_zpe',
                       state_id='base-validation-state'))
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


def audit_composite_run(run_dir):
    """Audit one unified geometry and ANL property graph."""
    run_dir, spec, state = _load(run_dir)
    if spec.get('name') != 'anl-composite-validation':
        raise ValueError('Run is not a unified ANL validation graph.')
    statuses = {
        task['id']: state['tasks'].get(task['id'], {}).get('status', 'waiting')
        for task in spec['tasks']}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'Composite ANL run is incomplete: {statuses}')
    atomic = len(spec['molecule'].get('symbols', ())) == 1
    stationary_point = spec.get('intent', {}).get(
        'stationary_point', 'minimum')
    if stationary_point == 'transition_state':
        vpt2_summary = {
            'review_required': False, 'warnings': [],
            'provider': 'not-applicable-transition-state-harmonic-zpe',
        }
    elif atomic:
        vpt2_summary = {
            'review_required': False, 'warnings': [],
            'anharmonic_correction_hartree': 0.,
            'provider': 'known_zero_atomic_vibration',
        }
    else:
        *_, vpt2 = _verified_task_result(run_dir, 'gaussian_vpt2')
        vpt2_summary = {
            'review_required': vpt2.get('review_required', False),
            'warnings': vpt2.get('warnings', []),
            'anharmonic_correction_hartree':
                vpt2['anharmonic_correction_hartree'],
        }
    post = (audit_anl0_post_geometry_run(run_dir)
            if spec.get('intent', {}).get('anl0_only')
            else audit_post_geometry_run(run_dir))
    return {
        'status': 'composite_anl_interface_complete',
        'requested_ladder_head':
            spec.get('intent', {}).get('requested_ladder_head'),
        'task_statuses': statuses,
        'base': audit_current_base_run(run_dir),
        'post_geometry': post,
        'vpt2': vpt2_summary,
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
    smiles_build = commands.add_parser('from-smiles')
    smiles_build.add_argument('smiles')
    smiles_build.add_argument('spec', type=Path)
    smiles_build.add_argument('--charge', type=int, default=0)
    smiles_build.add_argument('--multiplicity', type=int, default=1)
    smiles_build.add_argument('--random-seed', type=int, default=312)
    smiles_build.add_argument('--max-nodes', type=int, default=3)
    smiles_build.add_argument('--partition')
    smiles_stage = commands.add_parser('prepare-from-smiles')
    smiles_stage.add_argument('smiles')
    smiles_stage.add_argument('run_dir', type=Path)
    smiles_stage.add_argument('--charge', type=int, default=0)
    smiles_stage.add_argument('--multiplicity', type=int, default=1)
    smiles_stage.add_argument('--random-seed', type=int, default=312)
    smiles_stage.add_argument('--max-nodes', type=int, default=3)
    smiles_stage.add_argument('--partition')
    composite_db = commands.add_parser('prepare-composite-from-db')
    composite_db.add_argument('database', type=Path)
    composite_db.add_argument('job')
    composite_db.add_argument('run_dir', type=Path)
    composite_db.add_argument('--charge', type=int, default=0)
    composite_db.add_argument('--multiplicity', type=int, default=1)
    composite_db.add_argument('--max-nodes', type=int, default=3)
    composite_db.add_argument('--partition')
    composite_db.add_argument('--mrcc-command', default='dmrcc')
    composite_db.add_argument('--anl0-only', action='store_true')
    composite_smiles = commands.add_parser('prepare-composite-from-smiles')
    composite_smiles.add_argument('smiles')
    composite_smiles.add_argument('run_dir', type=Path)
    composite_smiles.add_argument('--charge', type=int, default=0)
    composite_smiles.add_argument('--multiplicity', type=int, default=1)
    composite_smiles.add_argument('--random-seed', type=int, default=312)
    composite_smiles.add_argument('--max-nodes', type=int, default=3)
    composite_smiles.add_argument('--partition')
    composite_smiles.add_argument('--mrcc-command', default='dmrcc')
    composite_smiles.add_argument('--anl0-only', action='store_true')
    audit = commands.add_parser('audit')
    audit.add_argument('run_dir', type=Path)
    audit_composite = commands.add_parser('audit-composite')
    audit_composite.add_argument('run_dir', type=Path)
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
    prepare_post.add_argument('--anl0-only', action='store_true')
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
    if args.action == 'audit-composite':
        print(json.dumps(audit_composite_run(args.run_dir), indent=2,
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
    from_smiles = args.action in (
        'from-smiles', 'prepare-from-smiles',
        'prepare-composite-from-smiles')
    if from_run:
        molecule = molecule_from_completed_run(
            args.source_run, geometry_task=args.geometry_task)
    elif from_smiles:
        molecule = molecule_from_smiles(
            args.smiles, charge=args.charge, multiplicity=args.multiplicity,
            random_seed=args.random_seed)
    else:
        molecule = molecule_from_database(
            args.database, args.job, charge=args.charge,
            multiplicity=args.multiplicity)
    is_higher = args.action in ('higher-order-from-db',
                                'prepare-higher-order-from-db',
                                'prepare-higher-order-from-run')
    is_corrections = args.action == 'prepare-common-corrections-from-run'
    is_post = args.action == 'prepare-post-geometry-from-run'
    is_base = args.action == 'prepare-current-base-from-run'
    is_composite = args.action in ('prepare-composite-from-db',
                                   'prepare-composite-from-smiles')
    if is_composite:
        spec = composite_validation_spec(
            molecule, max_nodes=args.max_nodes, partition=args.partition,
            mrcc_command=args.mrcc_command, anl0_only=args.anl0_only)
    elif is_base:
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
            mrcc_command=args.mrcc_command, anl0_only=args.anl0_only)
    else:
        spec = interface_validation_spec(
            molecule, max_nodes=args.max_nodes, partition=args.partition)
    if args.action in ('from-db', 'from-smiles', 'higher-order-from-db'):
        args.spec.write_text(json.dumps(spec, indent=2) + '\n')
        print(args.spec.resolve())
    else:
        print(prepare(spec, args.run_dir))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
