"""Reusable, molecule-independent native QC task builders for ANL workflows."""

from __future__ import annotations


MOLPRO_HEADER = """***,KinBot ANL calculation
symmetry,nosym
orient,noorient
geomtyp=xyz
geometry={
{{XYZ}}
}
set,charge={{CHARGE}}
set,spin={{SPIN}}
"""


def resources(walltime, *, max_cores=8, min_stack_mw=None,
              min_memory_mb_per_core=None, partition=None):
    """Request one exclusive node with portable, memory-aware core sizing."""
    result = {'cores': 'auto', 'memory_mb': 'node', 'walltime': walltime,
              'max_cores': max_cores}
    if min_stack_mw is not None:
        result['min_stack_mw'] = min_stack_mw
    if min_memory_mb_per_core is not None:
        result['min_memory_mb_per_core'] = min_memory_mb_per_core
    if partition is not None:
        result['partition'] = partition
    return result


def molpro_task(ident, body, *, geometry_from='l3_geometry',
                walltime='04:00:00', max_cores=8, parser=None,
                partition=None):
    """Build one conventional Molpro task with native `.out` provenance."""
    task = {
        'id': ident, 'kind': 'external', 'backend': 'molpro',
        'geometry_from': geometry_from,
        'resources': resources(walltime, max_cores=max_cores,
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


def mrcc_task(ident, method, basis, *, multiplicity, reference=None,
              geometry_from='l3_geometry', walltime='24:00:00',
              max_cores=8, partition=None):
    """Build a direct MRCC higher-order task.

    Direct ``dmrcc`` uses RHF for a closed shell and semicanonical ROHF for an
    open shell by default.  A recipe that reproduces the older Ram/Elliott
    UHF higher-order convention can request ``reference='UHF'`` explicitly.
    One process uses OpenMP threads selected from node memory and the explicit
    performance cap; MRCC's replicated-memory MPI mode is not enabled.
    """
    if method not in ('CCSDT(Q)', 'CCSDTQ(P)'):
        raise ValueError(f'Unsupported MRCC method {method!r}.')
    if isinstance(multiplicity, bool) or not isinstance(multiplicity, int) \
            or multiplicity < 1:
        raise ValueError('MRCC multiplicity must be a positive integer.')
    if reference is None:
        reference = 'RHF' if multiplicity == 1 else 'ROHF'
    if reference not in ('RHF', 'ROHF', 'UHF'):
        raise ValueError(f'Unsupported MRCC reference {reference!r}.')
    if ((reference == 'RHF' and multiplicity != 1)
            or (reference == 'ROHF' and multiplicity == 1)):
        raise ValueError('MRCC reference conflicts with multiplicity.')
    rohf = ('rohftype=semicanonical\nrohfcore=semicanonical\n'
            if reference == 'ROHF' else '')
    output = f'{ident}.out'
    template = f"""basis={basis}
calc={method}
ccprog=mrcc
scftype={reference}
{rohf}core=frozen
charge={{{{CHARGE}}}}
mult={{{{MULT}}}}
cctol=10
mem={{{{WORK_MEMORY_MB}}}}MB
molden=off
symm=off
unit=angs
geom=xyz
{{{{MRCC_XYZ}}}}
"""
    return {
        'id': ident, 'kind': 'external', 'backend': 'mrcc',
        'geometry_from': geometry_from,
        'resources': resources(
            walltime, max_cores=max_cores,
            min_memory_mb_per_core=4096, partition=partition),
        'input_name': 'MINP', 'input_template': template,
        'command': ['dmrcc'],
        'required_executables': ['scf', 'mrcc'],
        'stdout': output, 'stderr': f'{ident}.err',
        'required_outputs': [output],
        'success_marker': {'file': output,
                           'contains': 'Normal termination of mrcc.'},
        'failure_markers': [
            {'file': output, 'contains': 'Error at the termination of mrcc.'},
            {'file': output, 'contains': 'Fatal error'},
        ],
        'result_parser': {
            'kind': 'mrcc_energy', 'file': output,
            'method': method, 'basis': basis,
            'reference': reference, 'correlation': 'unrestricted',
            'core': 'frozen', 'program': 'mrcc',
        },
    }


def higher_order_task(ident, method, basis, *, multiplicity, reference=None,
                      geometry_from='l3_geometry', walltime='24:00:00',
                      max_cores=8, partition=None):
    """Select the pinned ANL higher-order backend for one electronic state.

    Every higher-order calculation uses direct MRCC's general spin-orbital CC
    implementation. RHF is the default determinant for a closed shell and
    semicanonical ROHF for an open shell; an older-paper reproduction can pin
    UHF. No partially spin-adapted or restricted CC ansatz is generated.
    """
    options = dict(
        geometry_from=geometry_from, walltime=walltime,
        max_cores=max_cores, partition=partition)
    if method in ('CCSDT(Q)', 'CCSDTQ(P)'):
        return mrcc_task(
            ident, method, basis, multiplicity=multiplicity,
            reference=reference, **options)
    raise ValueError(f'Unsupported higher-order method {method!r}.')
