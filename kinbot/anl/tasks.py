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


def molpro_ccsdt_command(multiplicity, *, all_electron=False):
    """Return the Molpro CCSD(T) command for the ANL reference convention.

    Molpro's ``RHF`` command supplies an RHF determinant for a closed shell and
    an ROHF determinant for an open shell.  Bare ``UCCSD(T)`` then selects the
    RHF/ROHF-UCCSD(T) formalism used by the ANL work.  Do not add
    ``UHF_UCCSD=1``: that forces Molpro's separate UHF-UCC engine and Molpro
    2024.1 can then print an RHF-UCCSD(T) label with a zero (T) contribution
    for an ordinary closed-shell molecule.
    """
    if isinstance(multiplicity, bool) or not isinstance(multiplicity, int) \
            or multiplicity < 1:
        raise ValueError('Multiplicity must be a positive integer.')
    method = 'uccsd(t)'
    return f'{{{method};core}}' if all_electron else method


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


def cfour_energy_task(ident, method, basis, *, multiplicity,
                      geometry_from='l3_geometry', walltime='24:00:00',
                      max_cores=8, partition=None):
    """Build a legacy CFOUR VCC probe for diagnosing site capabilities.

    This builder is retained so old run directories remain readable.  It is
    not selected by :func:`higher_order_task`: CFOUR 2.1 routes CCSDT(Q) to
    ``xncc`` even when ``CC_PROG=VCC`` is present, and that closed-shell solver
    does not satisfy the unrestricted-CC ANL convention.  Its strict result
    parser therefore rejects such output instead of silently consuming it.
    """
    if method != 'CCSDT(Q)':
        raise ValueError(f'Unsupported CFOUR higher-order method {method!r}.')
    if multiplicity != 1:
        raise ValueError('Native CFOUR higher-order routing is closed-shell only.')
    output = f'{ident}.out'
    template = f"""KinBot ANL {method} calculation
{{{{CARTESIAN}}}}

*CFOUR(CALC={method}
BASIS={basis}
REFERENCE=RHF
CC_PROG=VCC
FROZEN_CORE=ON
COORDINATES=CARTESIAN
UNITS=ANGSTROM
SYMMETRY=OFF
CHARGE={{{{CHARGE}}}}
MULTIPLICITY={{{{MULT}}}}
SCF_CONV=10
LINEQ_CONV=10
CC_MAXCYC=300
MEM_UNIT=MB
MEMORY_SIZE={{{{WORK_MEMORY_MB}}}})
"""
    return {
        'id': ident, 'kind': 'external', 'backend': 'cfour',
        'geometry_from': geometry_from,
        'resources': resources(
            walltime, max_cores=max_cores,
            min_memory_mb_per_core=4096, partition=partition),
        'input_name': 'ZMAT', 'input_template': template,
        'command': ['xcfour'], 'stdout': output, 'stderr': f'{ident}.err',
        'required_outputs': [output],
        'files_from_env': {'GENBAS': 'CFOUR_GENBAS'},
        'success_marker': {'file': output,
                           'contains': 'The final electronic energy is'},
        'failure_markers': [
            {'file': output, 'contains': 'ERROR ERROR ERROR'},
            {'file': output, 'contains': 'Job has terminated with error flag'},
        ],
        'result_parser': {
            'kind': 'cfour_energy', 'file': output,
            'method': method, 'basis': basis, 'reference': 'RHF',
            'correlation': 'unrestricted', 'core': 'frozen',
            'program': 'cfour', 'driver': 'VCC',
        },
    }


def mrcc_task(ident, method, basis, *, multiplicity, reference=None,
              geometry_from='l3_geometry', walltime='24:00:00',
              max_cores=8, partition=None, command='dmrcc'):
    """Build a direct MRCC higher-order task.

    Direct ``dmrcc`` uses RHF for a closed shell and semicanonical ROHF for an
    open shell by default.  An explicit UHF determinant is supported for
    reproducing source data whose numerical convention requires it.  The
    correlation treatment remains unrestricted.  One process uses OpenMP
    threads selected from node memory and the explicit performance cap;
    MRCC's replicated-memory MPI mode is not enabled.
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
            or (reference == 'ROHF' and multiplicity == 1)
            or (reference == 'UHF' and multiplicity == 1)):
        raise ValueError('MRCC reference conflicts with multiplicity.')
    if (not isinstance(command, str) or not command.strip()
            or any(character.isspace() for character in command)):
        raise ValueError('MRCC command must name one executable without arguments.')
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
        'command': [command],
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
                      max_cores=8, partition=None, command='dmrcc'):
    """Select the pinned ANL higher-order backend for one electronic state.

    Every CCSDT(Q) and CCSDTQ(P) calculation uses direct MRCC.  Closed shells
    use an RHF determinant and open shells use semicanonical ROHF by default;
    the correlation treatment remains unrestricted in both cases.  CFOUR 2.1
    forcibly selects its closed-shell ``xncc`` implementation for CCSDT(Q), so
    it cannot supply the requested RHF-UCCSDT(Q) component.
    """
    options = dict(
        geometry_from=geometry_from, walltime=walltime,
        max_cores=max_cores, partition=partition, command=command)
    if method in ('CCSDT(Q)', 'CCSDTQ(P)'):
        return mrcc_task(
            ident, method, basis, multiplicity=multiplicity,
            reference=reference, **options)
    raise ValueError(f'Unsupported higher-order method {method!r}.')
