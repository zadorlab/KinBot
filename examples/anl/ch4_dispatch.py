"""CH4 input for the general ANL dispatch smoke test.

This tests execution and artifact collection, not the full ANL0-F12 energy
expression. See docs/composite_qc_validation.md for the external acceptance
checks and the vendor documentation behind each input.
"""

import json
from pathlib import Path
import sys


MOLPRO_HEADER = """***,KinBot CH4 dispatch
symmetry,nosym
orient,noorient
geomtyp=xyz
geometry={
{{XYZ}}
}
"""


def resources(cores=4, memory_mb=16000, walltime='04:00:00'):
    return {'cores': cores, 'memory_mb': memory_mb, 'walltime': walltime}


def molpro_task(ident, body, *, geometry_from='l3_geometry',
                geometry_output=None, cores=4, memory_mb=16000,
                walltime='04:00:00'):
    task = {
        'id': ident, 'kind': 'external', 'backend': 'molpro',
        'geometry_from': geometry_from,
        'resources': resources(cores, memory_mb, walltime),
        'input_name': f'{ident}.inp',
        'input_template': MOLPRO_HEADER + body,
        'command': ['molpro', '-n', '{cores}', '-M',
                    '{molpro_total_mw}', '{input}'],
        'stdout': 'launcher.stdout', 'stderr': 'launcher.stderr',
        'required_outputs': [f'{ident}.out'],
        'success_marker': {'file': f'{ident}.out',
                           'contains': 'Molpro calculation terminated'},
    }
    if geometry_output:
        task['geometry_output'] = geometry_output
        task['required_outputs'] += [f'{ident}.log', geometry_output]
    return task


def ch4_spec():
    # A near-tetrahedral, slightly imperfect starting geometry gives the two
    # optimizers something to do while keeping this a single-species test.
    return {
        'schema': 1,
        'name': 'ch4-qc-dispatch',
        'molecule': {
            'symbols': ['C', 'H', 'H', 'H', 'H'],
            'positions': [
                [0.0, 0.0, 0.0],
                [0.638, 0.627, 0.634],
                [-0.630, -0.630, 0.630],
                [-0.630, 0.630, -0.630],
                [0.630, -0.630, -0.630],
            ],
            'charge': 0, 'multiplicity': 1,
        },
        'limits': {
            'max_nodes': 3,
            'max_cores_per_node': 16,
            'max_memory_mb_per_node': 64000,
        },
        'tasks': [
            {
                'id': 'l2_geometry', 'kind': 'ase_optimize',
                'geometry_from': 'initial', 'geometry_output': 'final.xyz',
                'resources': resources(4, 16000, '02:00:00'),
                'profile': {
                    'calculator': 'gaussian', 'method': 'B2PLYP',
                    'basis': 'cc-pVTZ', 'command': 'g16',
                    'calculator_kwargs': {
                        'EmpiricalDispersion': 'GD3BJ', 'Symm': 'None',
                        'scf': 'xqc',
                    },
                    'optimizer': 'sella',
                },
                'optimizer': {'fmax': 0.03, 'steps': 80,
                              'sella_kwargs': {'internal': True}},
            },
            molpro_task(
                'l3_geometry',
                'basis=cc-pVTZ\nhf\nccsd(t)\n'
                'optg,numerical,savexyz=l3_geometry.xyz\n',
                geometry_from='l2_geometry',
                geometry_output='l3_geometry.xyz',
                cores=8, memory_mb=32000, walltime='08:00:00',
            ),
            molpro_task(
                'harmonic',
                'basis=cc-pVTZ\nhf\nccsd(t)\nfrequencies,numerical\n',
                cores=8, memory_mb=32000, walltime='12:00:00',
            ),
            molpro_task(
                'f12_tz',
                'basis=cc-pVTZ-F12\nhf\nccsd(t)-f12,scale_trip=1\n'
                'kb_f12b=energy(2)\n',
                cores=8, memory_mb=32000, walltime='06:00:00',
            ),
            molpro_task(
                'f12_qz',
                'basis=cc-pVQZ-F12\nhf\nccsd(t)-f12,scale_trip=1\n'
                'kb_f12b=energy(2)\n',
                cores=8, memory_mb=48000, walltime='12:00:00',
            ),
            molpro_task(
                'ccsdt_dz',
                'basis=cc-pVDZ\nhf\nccsd(t)\nkb_ccsdt=energy\n',
            ),
            {
                'id': 'cfour_dboc', 'kind': 'external', 'backend': 'cfour',
                'geometry_from': 'l3_geometry',
                'resources': resources(4, 16000, '04:00:00'),
                'input_name': 'ZMAT',
                'input_template': (
                    'KinBot CH4 DBOC\n{{CARTESIAN}}\n\n'
                    '*CFOUR(CALC=SCF,BASIS=cc-pVTZ,DBOC=ON,'
                    'COORD=CARTESIAN,UNITS=ANGSTROM,'
                    'CHARGE={{CHARGE}},MULTIPLICITY={{MULT}},'
                    'MEM_UNIT=MB,MEMORY_SIZE={{WORK_MEMORY_MB}})\n'),
                'command': ['xcfour'], 'stdout': 'cfour.out',
                'stderr': 'cfour.err', 'required_outputs': ['cfour.out'],
                'files_from_env': {'GENBAS': 'CFOUR_GENBAS'},
                'success_marker': {
                    'file': 'cfour.out',
                    'contains': 'The total diagonal Born-Oppenheimer correction (DBOC) is:'},
            },
            {
                'id': 'mrcc_ccsdtq', 'kind': 'external', 'backend': 'mrcc',
                'geometry_from': 'l3_geometry',
                'resources': resources(8, 32000, '12:00:00'),
                'input_name': 'MINP',
                'input_template': (
                    'calc=CCSDT(Q)\nbasis=cc-pVDZ\nscftype=RHF\n'
                    'core=frozen\ngauss=spher\n'
                    'mem={{WORK_MEMORY_MB}}MB\n'
                    'charge={{CHARGE}}\nmult={{MULT}}\nunit=angs\n'
                    'geom=xyz\n{{MRCC_XYZ}}\n'),
                'command': ['dmrcc'], 'stdout': 'mrcc.out',
                'stderr': 'mrcc.err', 'required_outputs': ['mrcc.out'],
                'setup': ['export MKL_NUM_THREADS="$OMP_NUM_THREADS"'],
                'success_marker': {'file': 'mrcc.out',
                                   'contains': 'Normal termination of mrcc.'},
                'failure_markers': [
                    {'file': 'mrcc.out',
                     'contains': 'Error at the termination of mrcc.'},
                    {'file': 'mrcc.out', 'contains': 'Fatal error'},
                ],
            },
            {
                'id': 'gaussian_vpt2', 'kind': 'external',
                'backend': 'gaussian', 'geometry_from': 'l3_geometry',
                'resources': resources(4, 16000, '08:00:00'),
                'input_name': 'vpt2.com',
                'input_template': (
                    '%nprocshared={{CORES}}\n%mem={{WORK_MEMORY_MB}}MB\n'
                    '#p B3LYP/cc-pVTZ Opt=(Tight,CalcFC) '
                    'Freq=Anharmonic NoSymm SCF=XQC\n\n'
                    'KinBot CH4 VPT2\n\n{{CHARGE}} {{MULT}}\n'
                    '{{CARTESIAN}}\n\n'),
                'command': ['g16'], 'stdin': 'vpt2.com',
                'stdout': 'vpt2.log', 'stderr': 'vpt2.err',
                'required_outputs': ['vpt2.log'],
                'success_marker': {'file': 'vpt2.log',
                                   'contains': 'Normal termination of Gaussian'},
            },
        ],
    }


if __name__ == '__main__':
    target = Path(sys.argv[1] if len(sys.argv) > 1 else 'ch4_dispatch.json')
    target.write_text(json.dumps(ch4_spec(), indent=2) + '\n')
    print(target)
