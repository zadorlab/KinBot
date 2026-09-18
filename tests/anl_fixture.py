"""Small synthetic task graph for dispatcher unit tests; never submitted to QC."""


MOLPRO_HEADER = """***,KinBot dispatch fixture
symmetry,nosym
orient,noorient
geomtyp=xyz
geometry={
{{XYZ}}
}
set,charge={{CHARGE}}
set,spin={{SPIN}}
"""


def resources(cores=4, memory_mb=16000, walltime='04:00:00'):
    return {'cores': cores, 'memory_mb': memory_mb, 'walltime': walltime}


def molpro_task(ident, body, *, geometry_from='l3_geometry',
                geometry_output=None, cores=4, memory_mb=16000,
                walltime='04:00:00', result_parser=None):
    task = {
        'id': ident, 'kind': 'external', 'backend': 'molpro',
        'geometry_from': geometry_from,
        'resources': resources(cores, memory_mb, walltime),
        'input_name': f'{ident}.inp',
        'input_template': MOLPRO_HEADER + body,
        'command': ['molpro', '-n', '{cores}', '-m',
                    '{molpro_stack_mw}', '{input}'],
        'stdout': 'launcher.stdout', 'stderr': 'launcher.stderr',
        'required_outputs': [f'{ident}.out'],
        'success_marker': {'file': f'{ident}.out',
                           'contains': 'Molpro calculation terminated'},
    }
    if geometry_output:
        task['geometry_output'] = geometry_output
        task['required_outputs'] += [f'{ident}.log', geometry_output]
    if result_parser:
        task['result_parser'] = {'file': f'{ident}.out', **result_parser}
    return task


def dispatch_spec(*, auto_resources=False):
    # Slightly imperfect methane coordinates exercise geometry propagation.
    spec = {
        'schema': 1,
        'name': 'synthetic-qc-dispatch',
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
                        'scf': 'xqc', 'integral': 'UltraFine',
                    },
                    'optimizer': 'sella',
                },
                'optimizer': {'fmax': 0.0005, 'steps': 120,
                              'sella_kwargs': {'internal': True}},
            },
            {
                'id': 'l3_geometry', 'kind': 'ase_optimize',
                'geometry_from': 'l2_geometry',
                'geometry_output': 'final.xyz',
                'resources': resources(8, 32000, '24:00:00'),
                'profile': {
                    'calculator': 'molpro', 'method': 'CCSD(T)',
                    'basis': 'cc-pVTZ', 'command': 'molpro',
                    'optimizer': 'sella',
                },
                'optimizer': {'fmax': 0.03, 'steps': 80,
                              'sella_kwargs': {'internal': True}},
            },
            molpro_task(
                'harmonic',
                'basis=cc-pVTZ\nhf\nccsd(t)\nfrequencies,numerical\n',
                cores=8, memory_mb=32000, walltime='12:00:00',
                result_parser={'kind': 'molpro_harmonic', 'basis': 'cc-pVTZ'},
            ),
            molpro_task(
                'f12_tz',
                'basis=cc-pVTZ-F12\nhf\nccsd(t)-f12,scale_trip=1\n'
                'kb_f12b=energy(2)\n',
                cores=8, memory_mb=32000, walltime='06:00:00',
                result_parser={'kind': 'molpro_energy',
                               'method': 'CCSD(T)-F12b', 'basis': 'cc-pVTZ-F12'},
            ),
            molpro_task(
                'f12_qz',
                'basis=cc-pVQZ-F12\nhf\nccsd(t)-f12,scale_trip=1\n'
                'kb_f12b=energy(2)\n',
                cores=8, memory_mb=48000, walltime='12:00:00',
                result_parser={'kind': 'molpro_energy',
                               'method': 'CCSD(T)-F12b', 'basis': 'cc-pVQZ-F12'},
            ),
            molpro_task(
                'molpro_dz_sp',
                'basis=cc-pVDZ\nhf\nccsd(t)\nkb_dz_energy=energy\n',
                result_parser={'kind': 'molpro_energy',
                               'method': 'CCSD(T)', 'basis': 'cc-pVDZ'},
            ),
            {
                'id': 'cfour_dboc', 'kind': 'external', 'backend': 'cfour',
                'geometry_from': 'l3_geometry',
                'resources': resources(4, 16000, '04:00:00'),
                'input_name': 'ZMAT',
                'input_template': (
                    'KinBot test DBOC\n{{CARTESIAN}}\n\n'
                    '*CFOUR(CALC=SCF\n'
                    'BASIS=cc-pVTZ\n'
                    'DBOC=ON\n'
                    'COORD=CARTESIAN\n'
                    'UNITS=ANGSTROM\n'
                    'CHARGE={{CHARGE}}\n'
                    'MULTIPLICITY={{MULT}}\n'
                    'MEM_UNIT=MB\n'
                    'MEMORY_SIZE={{WORK_MEMORY_MB}})\n'),
                'command': ['xcfour'], 'stdout': 'cfour.out',
                'stderr': 'cfour.err', 'required_outputs': ['cfour.out'],
                'files_from_env': {'GENBAS': 'CFOUR_GENBAS'},
                'success_marker': {
                    'file': 'cfour.out',
                    'contains': 'The total diagonal Born-Oppenheimer correction (DBOC) is:'},
                'result_parser': {
                    'kind': 'cfour_dboc', 'file': 'cfour.out',
                    'level': 'HF', 'basis': 'cc-pVTZ'},
            },
            # Frequency-only VPT2 uses the accepted, tier-matched L2 geometry.
            # The L3 dependency keeps it in the post-geometry fan-out.
            {
                'id': 'gaussian_vpt2', 'kind': 'external',
                'backend': 'gaussian', 'geometry_from': 'l2_geometry',
                'depends_on': ['l3_geometry'],
                'resources': resources(4, 16000, '24:00:00'),
                'input_name': 'vpt2.com',
                'input_template': (
                    '%nprocshared={{CORES}}\n%mem={{WORK_MEMORY_MB}}MB\n'
                    '#p B2PLYP/cc-pVTZ Freq=Anharmonic NoSymm SCF=XQC '
                    'EmpiricalDispersion=GD3BJ Integral=UltraFine\n\n'
                    'KinBot test VPT2\n\n{{CHARGE}} {{MULT}}\n'
                    '{{CARTESIAN}}\n\n'),
                'command': ['g16'], 'stdin': 'vpt2.com',
                'stdout': 'vpt2.log', 'stderr': 'vpt2.err',
                'required_outputs': ['vpt2.log'],
                'success_marker': {'file': 'vpt2.log',
                                   'contains': 'Normal termination of Gaussian'},
                'result_parser': {'kind': 'gaussian_vpt2', 'file': 'vpt2.log',
                                  'method': 'B2PLYP', 'basis': 'cc-pVTZ',
                                  'dispersion': 'GD3BJ'},
            },
        ],
    }
    if auto_resources:
        spec['limits'] = {'max_nodes': spec['limits']['max_nodes']}
        for task in spec['tasks']:
            task['resources'].pop('cores')
            task['resources'].pop('memory_mb')
    return spec
