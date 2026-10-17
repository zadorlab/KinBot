"""Exercise the generated worker with native-format logs and a model surface.

Only Gaussian and optimization are substituted. Real frequency/energy/mode
readers, geometry alignment, checkpoint branching and result publication run.
This is a workflow test, not native Gaussian validation of a chemical PES.
"""
import json
import pickle
import ast
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write
from ase.units import Hartree

from kinbot.minimum_recovery import aligned_mode, worker_sella_kwargs, pop_worker_controls
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry


def test_align_mode_rotates_without_reflecting_or_reordering_atoms():
    source = np.array([[0., 0., 0.], [1., .2, .3], [.1, 1., .5], [.2, .3, 1.]])
    rotation = np.array([[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]])
    mode = np.arange(12.).reshape(4, 3)
    target = source @ rotation + 7.
    np.testing.assert_allclose(aligned_mode(source, target, mode),
                               mode @ rotation / np.linalg.norm(mode), atol=1e-14)
    with pytest.raises(ValueError, match='differs'):
        aligned_mode(source, -target, mode)
    with pytest.raises(ValueError, match='differs'):
        aligned_mode(source, target[[0, 2, 1, 3]], mode)
    with pytest.raises(ValueError, match='zero'):
        aligned_mode(source, target, mode * 0)
    with pytest.raises(ValueError, match='Invalid'):
        aligned_mode(source, target, mode * np.nan)


def test_controls_are_removed_before_sella_and_old_formatters_need_no_new_fields():
    from string import Formatter
    options = {'gamma': 0., 'internal': True}
    old_options = dict(options)
    assert pop_worker_controls(old_options)['minimum_recovery_attempts'] == 2
    assert old_options == options
    wrapped = worker_sella_kwargs({'sella_kwargs': options,
                                  'minimum_recovery_attempts': 0,
                                  'minimum_recovery_displacement': .03,
                                  'imagfreq_threshold': 75.})
    controls = pop_worker_controls(wrapped)
    assert controls == {'minimum_recovery_attempts': 0,
                        'minimum_recovery_displacement': .03, 'imagfreq_threshold': 75.}
    assert wrapped == options == {'gamma': 0., 'internal': True}
    for name in ('ase_sella_opt_well.tpl.py', 'ase_sella_ts_end.tpl.py'):
        text = (Path(__file__).parents[1]/'kinbot/tpl'/name).read_text()
        fields = {field for _, field, _, _ in Formatter().parse(text) if field}
        assert not fields.intersection(controls)


def test_private_controls_only_reach_templates_that_consume_them(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps({'theory_preset': 'uma-b2plyp-anl', 'barrier_threshold': 1000}))
    qc = QuantumChemistry(Parameters('input.json').par).l1
    monkeypatch.setattr(qc, 'submit_qc', lambda *a, **k: 1)
    species = SimpleNamespace(chemid=123, name='probe', mult=1, charge=0, nel=10,
        natom=3, wellorts=0, atom=['O', 'H', 'H'],
        geom=np.array([[0., 0., 0.], [.9, 0., 0.], [-.2, .8, 0.]]))
    Path('conf').mkdir()
    qc.qc_opt(species, species.geom)
    qc.qc_conf(species, species.geom, 0)
    qc.qc_opt_ts(species, species.geom)
    for path in (Path('123_well.py'), Path('conf/123_0000.py'), Path('probe.py')):
        tree = ast.parse(path.read_text())
        assignments = [node for node in ast.walk(tree) if isinstance(node, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id == 'sella_kwargs' for t in node.targets)]
        assert len(assignments) == 1
        options = ast.literal_eval(assignments[0].value)
        assert options == qc.par['sella_kwargs']


def native_log(atoms, energy, frequencies, direction):
    """Different native orientation, same atom order and physical geometry."""
    rotation = np.array([[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]])
    xyz = atoms.positions @ rotation + [2., 3., 4.]
    mode = direction @ rotation
    lines = [' Standard orientation:', ' -----', ' Center Atomic Atomic Coordinates',
             ' Number Number Type X Y Z', ' -----']
    for i, (number, position) in enumerate(zip(atoms.numbers, xyz), 1):
        lines.append(f' {i} {number} 0 ' + ' '.join(f'{x:.10f}' for x in position))
    lines += [' -----', f' SCF Done: E(UB2PLYP) = {energy / Hartree + .1:.14f} A.U.',
              f' E2(B2PLYP) = -0.1 E(B2PLYP) = {energy / Hartree:.14f}',
              ' Frequencies -- ' + ' '.join(map(str, frequencies)),
              '  Atom  AN      X      Y      Z']
    for i, (number, vector) in enumerate(zip(atoms.numbers, mode), 1):
        lines.append(f' {i} {number} ' + ' '.join(f'{x:.10f}' for x in vector))
    return '\n'.join(lines) + ('\n Zero-point correction= 0.01 (Hartree/Particle)'
                              '\n Normal termination of Gaussian 16\n')


@pytest.mark.parametrize('conformer', [False, True])
@pytest.mark.parametrize('legacy_formatter', [False, True])
@pytest.mark.parametrize('outcome', ['minimum', 'still_saddle', 'uphill', 'probe_failure'])
def test_generated_minimum_recovers_or_rejects_with_fresh_hessian(
        tmp_path, monkeypatch, conformer, legacy_formatter, outcome):
    monkeypatch.chdir(tmp_path)
    from kinbot.ase_modules.calculators import gaussian
    import sella
    base = np.array([[0., 0., 0.], [.9, 0., 0.], [-.2, .8, 0.]])
    direction = np.array([[0., 0., 0.], [0., 0., 1.], [0., 0., -1.]]) / np.sqrt(2.)
    calls, optimization_starts, probe_checkpoints = [], [], []

    class ModelGaussian(Calculator):
        implemented_properties = ['energy', 'forces']

        def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            q = np.sum((atoms.positions-base)*direction)
            probe = self.prefix == 'probe'
            checkpoint = Path(self.directory) / (str(self.parameters['chk']).removesuffix('.chk') + '.chk')
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            if probe:
                probe_checkpoints.append(checkpoint.read_text())
                if outcome == 'probe_failure':
                    raise RuntimeError('probe QC failed')
            energy = -10. - q*q + 10.*q**4
            if probe and outcome == 'uphill':
                energy = -9.
            self.results = {'energy': energy, 'forces': (2.*q-40.*q**3)*direction}
            frequencies = ([100., 200., 300.] if abs(q) > .2
                           else [-587.334, 200., 300.])
            Path(self.label+'.log').write_text(native_log(atoms, energy, frequencies, direction))
            Path(self.label+'.com').write_text(str(dict(self.parameters)))
            checkpoint.write_text('hessian checkpoint' if 'freq' in self.parameters
                                  else 'probe checkpoint' if probe else 'force checkpoint')
            calls.append((self.label, dict(self.parameters)))

    class ModelOptimizer:
        def __init__(self, atoms, trajectory, logfile, **kwargs):
            self.atoms, self.trajectory, self.logfile = atoms, trajectory, logfile

        def run(self, **kwargs):
            optimization_starts.append(self.atoms.positions.copy())
            assert 'freq' not in self.atoms.calc.parameters
            q = np.sum((self.atoms.positions-base)*direction)
            if q and outcome == 'minimum':
                self.atoms.positions = base + np.sign(q)*np.sqrt(.05)*direction
            elif outcome == 'still_saddle':
                self.atoms.positions = base.copy()
            self.atoms.get_forces()
            write(self.trajectory, self.atoms)
            with Path(self.logfile).open('a') as stream:
                stream.write('Sella   0 19:11:56 -10.0 0.001 0.0\n')
            return True

    monkeypatch.setattr(gaussian, 'Gaussian', ModelGaussian)
    monkeypatch.setattr(sella, 'Sella', ModelOptimizer)
    if legacy_formatter:
        monkeypatch.setattr('kinbot.qc.worker_sella_kwargs', lambda par: dict(par['sella_kwargs']))
    Path('input.json').write_text(json.dumps({'theory_preset': 'uma-b2plyp-anl',
                                            'barrier_threshold': 1000}))
    qc = QuantumChemistry(Parameters('input.json').par).l2
    monkeypatch.setattr(qc, 'submit_qc', lambda *a, **k: 1)
    species = SimpleNamespace(chemid=123, name='probe', mult=1, charge=0, nel=10,
                             natom=3, wellorts=0, atom=['O', 'H', 'H'], geom=base)
    if conformer:
        Path('conf').mkdir()
        qc.qc_conf(species, base, 0)
        label = 'conf/123_0000'
    else:
        qc.qc_opt(species, base, high_level=1)
        label = '123_well_high'
    script = Path(label+'.py')
    tree = ast.parse(script.read_text())
    if legacy_formatter and conformer:
        # Older qc_conf omitted chk but subsequently requested Guess=Read.
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'kwargs'
                                                   for t in node.targets):
                kwargs = ast.literal_eval(node.value)
                kwargs.pop('chk')
                node.value = ast.parse(repr(kwargs), mode='eval').body
        ast.fix_missing_locations(tree)
    exec(compile(tree, str(script), 'exec'), {'__name__': '__main__'})
    with Path(label+'.pkl').open('rb') as stream:
        payload = pickle.load(stream)
    data = payload['data']
    assert json.loads(Path(label+'.result.json').read_text())['status'] == data['status']
    reports = [json.loads(p.read_text()) for p in Path(label+'_minimum_recovery').glob('*/recovery.json')]
    hessians = [p for _, p in calls if 'freq' in p]
    assert all(p['method'] == 'B2PLYP' and p['basis'] == 'cc-pVTZ' for _, p in calls)
    assert all(p['EmpiricalDispersion'] == 'GD3BJ' for _, p in calls)
    # Both branches start from the SAME rejected reference, independently.
    assert set(probe_checkpoints) == {'hessian checkpoint'}
    np.testing.assert_allclose(optimization_starts[0], base)
    for report in reports:
        assert Path(report['archive'], Path(label).name+'.log').is_file()
        assert Path(report['archive'], Path(label).name+'.chk').read_text() == 'hessian checkpoint'
    if outcome == 'minimum':
        assert data['status'] == 'normal' and data['minimum_recovery_attempts'] == 1
        assert len(hessians) == 2 and len(optimization_starts) == 2
        assert all(f > 0 for f in data['frequencies'])
        assert np.linalg.norm(optimization_starts[1] - base) == pytest.approx(.02)
        assert data['energy'] < -10.
    elif outcome == 'still_saddle':
        assert data['status'] == 'error' and data['minimum_recovery_attempts'] == 2
        assert len(hessians) == 3 and len(reports) == 2
        assert 'Frequency validation failed' in data['error']
        assert sorted(r['displacement_angstrom'] for r in reports) == [.02, .04]
    else:
        assert data['status'] == 'error' and len(hessians) == 1
        assert data['minimum_recovery_attempts'] == 0
        assert reports[0]['status'] == 'failed'
        assert 'Neither' in data['error'] if outcome == 'uphill' else 'probe QC failed' in data['error']
        np.testing.assert_allclose(payload['pos'], base)
        assert Path(label+'.chk').read_text() == 'hessian checkpoint'


@pytest.mark.parametrize('key,value', [('minimum_recovery_attempts', -1),
    ('minimum_recovery_attempts', True), ('minimum_recovery_attempts', 1.5),
    ('minimum_recovery_displacement', 0), ('minimum_recovery_displacement', float('nan')),
    ('minimum_recovery_displacement', True)])
def test_invalid_recovery_settings_fail_before_submission(tmp_path, key, value):
    path = tmp_path/'input.json'
    path.write_text(json.dumps({key: value}))
    with pytest.raises(ValueError, match=key):
        Parameters(str(path))
