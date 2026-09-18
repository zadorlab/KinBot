"""Molpro ASE/Sella contract tests; the Molpro executable is simulated.

The numerical-gradient table matches the format seen in the real Molpro 2024
CH4 first step. The revised PUT,XYZ deck remains an offsite acceptance gate.
"""

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

from ase import Atoms
from ase.units import Bohr, Hartree
import numpy as np
import pytest

from examples.anl.ch4_dispatch import ch4_spec
from kinbot.anl.dispatch import advance, prepare, run_task
from kinbot.ase_modules.calculators.molpro import (
    Molpro, check_process_count, parse_numerical_gradient, parse_output, render_input,
)


FAKE_MOLPRO = r'''#!{python}
import math
from pathlib import Path
import sys

args = sys.argv[1:]
assert args[:5] == ['-g', '-n', '8', '-m', '225'], args
inp = Path(args[5])
text = inp.read_text()
assert 'forces,numerical,variable=kb_geom_energy,startcmd=rhf' in text.lower()
assert 'put,xyz,' in text.lower() and 'xyzgrad' not in text.lower()
assert ('kb_geom_energy=energy\n'
        'forces,numerical,variable=kb_geom_energy,startcmd=rhf') in text.lower()
assert 'optg' not in text.lower()
assert 'memory,' not in text.lower()
assert 'basis=cc-pVTZ' in text
assert 'orient,noorient' in text
block = text.split('geometry={', 1)[1].split('}', 1)[0].strip().splitlines()
n = int(block[0])
rows = [(line.split()[0], [float(v) for v in line.split()[1:4]])
        for line in block[2:2+n]]
assert len(rows) == n == 5
target = [[0., 0., 0.], [.700, .627, .634], [-.630, -.630, .630],
          [-.630, .630, -.630], [.630, -.630, -.630]]
k = 0.1
energy = -40.0 + k * sum((value - goal)**2
                          for (_, xyz), goals in zip(rows, target)
                          for value, goal in zip(xyz, goals))
grad = [[2*k*1.8897261254578281*(value - goal)
         for value, goal in zip(xyz, goals)]
        for (_, xyz), goals in zip(rows, target)]
stem = inp.stem
Path(stem + '.out').write_text(
    f' SETTING KB_GEOM_ENERGY = {energy:.14f} AU\n'
    f' SETTING KB_GEOM_ENERGY = {energy + 0.1:.14f} AU\n'
    ' Numerical gradient for KB_GEOM_ENERGY\n'
    ' Atom          dE/dx               dE/dy               dE/dz'
    '                  d2E/dx2             d2E/dy2             d2E/dz2\n' +
    ''.join(f' {index:3d} {force[0]: .12f} {force[1]: .12f} '
            f'{force[2]: .12f} 0.0 0.0 0.0\n'
            for index, force in enumerate(grad, 1)) +
    ' Molpro calculation terminated\n')
Path(stem + '.log').write_text('Molpro gradient log\n')
Path(stem + '.xyz').write_text(
    str(n) + '\nMolpro PUT,XYZ geometry\n' +
    ''.join(f'{symbol} {x:.12f} {y:.12f} {z:.12f}\n'
            for symbol, (x, y, z) in rows))
'''


def test_documented_molpro_force_input_and_parser_contract():
    atoms = Atoms('CH2', positions=[[0, 0, 0], [0, 0, 1], [0, 1, 0]])
    deck = render_input(atoms, basis='cc-pVTZ', geometry_name='step.xyz')
    assert 'set,charge=0\nset,spin=0\ngthresh,energy=1.d-9\n' in deck
    assert ('rhf\nccsd(t)\nkb_geom_energy=energy\n'
            'forces,numerical,variable=kb_geom_energy,startcmd=rhf\n'
            'put,xyz,step.xyz') in deck
    assert 'optg' not in deck.lower()
    assert '3\nKinBot ASE geometry (angstrom)\nC' in deck
    with TemporaryDirectory() as temporary:
        output = Path(temporary) / 'step.out'
        output.write_text(' SETTING KB_GEOM_ENERGY        =       -40.123456789  AU\n'
                          ' SETTING KB_GEOM_ENERGY = -39.0 AU\n'
                          ' Molpro calculation terminated\n')
        assert parse_output(output) == pytest.approx(-40.123456789 * Hartree)
        geometry = Path(temporary) / 'step.xyz'
        geometry.write_text(
            '3\nMolpro XYZ\n'
            'H 0 1 0\n'
            'C 0 0 0\n'
            'H 0 0 1\n')
        output.write_text(output.read_text() +
                          ' Numerical gradient for KB_GEOM_ENERGY\n'
                          ' Atom dE/dx dE/dy dE/dz d2E/dx2 d2E/dy2 d2E/dz2\n'
                          ' 1 0.1 0.2 0.3 0 0 0\n'
                          ' 2 -1 -2 -3 0 0 0\n'
                          ' 3 0.4 0.5 0.6 0 0 0\n')
        np.testing.assert_allclose(parse_numerical_gradient(output, geometry, atoms),
                                   np.array([[1, 2, 3], [-.4, -.5, -.6],
                                             [-.1, -.2, -.3]]) * Hartree / Bohr)
        output.write_text(output.read_text().replace('Molpro calculation terminated',
                                                    'ERROR EXIT'))
        with pytest.raises(ValueError, match='terminate normally'):
            parse_output(output)
        # Molpro 2024.1 CH4 OPTG example uses a different F12 energy label.
        # It must never be accepted as conventional CCSD(T) geometry energy.
        output.write_text(' Version 2024.1\n'
                          ' !RHF-UCCSD(T)-F12 energy             -40.448936939114\n'
                          ' Molpro calculation terminated\n')
        with pytest.raises(ValueError, match='KB_GEOM_ENERGY missing'):
            parse_output(output)


def test_molpro_launcher_process_count_must_fit_allocation():
    with TemporaryDirectory() as temporary:
        output = Path(temporary) / 'step.out'
        output.write_text(' Distribution of processes:   nprocs(total)=   12   '
                          'nprocs(compute)=   11   nprocs(helper)=    1\n')
        with pytest.raises(RuntimeError, match='12 MPI processes.*8 ranks'):
            check_process_count(output, 8)
        assert check_process_count(output, 12) == 12


def test_molpro_numerical_gradient_rejects_missing_or_mismatched_data():
    atoms = Atoms('CH', positions=[[0, 0, 0], [0, 0, 1]])
    with TemporaryDirectory() as temporary:
        output = Path(temporary) / 'step.out'
        geometry = Path(temporary) / 'step.xyz'
        geometry.write_text('2\nMolpro XYZ\nC 0 0 0\nH 0 0 1\n')
        output.write_text('Molpro calculation terminated\n')
        with pytest.raises(ValueError, match='exactly one'):
            parse_numerical_gradient(output, geometry, atoms)
        output.write_text('Numerical gradient for KB_GEOM_ENERGY\n'
                          'Atom dE/dx dE/dy dE/dz d2E/dx2 d2E/dy2 d2E/dz2\n'
                          '1 0.1 0.2 0.3 0 0 0\n'
                          'Molpro calculation terminated\n')
        with pytest.raises(ValueError, match='incomplete'):
            parse_numerical_gradient(output, geometry, atoms)
        output.write_text(output.read_text().replace(
            'Molpro calculation terminated\n',
            '2 0.1 0.2 0.3 0 0 0\nMolpro calculation terminated\n'))
        geometry.write_text('2\nMolpro XYZ\nC 0 0 0\nH 0 0 2\n')
        with pytest.raises(ValueError, match='cannot be mapped'):
            parse_numerical_gradient(output, geometry, atoms)


def test_molpro_2024_ch4_gradient_rows_from_first_live_step():
    """Use the five actual numerical-gradient rows from Blodgett's .out."""
    atoms = Atoms('CH4', positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0],
                                   [0, 0, 1], [-1, 0, 0]])
    rows = np.array([
        [-0.000388912, 0.000324638, -0.000065297],
        [-0.000809204, -0.001219598, -0.000976978],
        [0.001072353, 0.001001028, -0.000987025],
        [0.000982618, -0.001151056, 0.000947788],
        [-0.000856854, 0.001044988, 0.001081512],
    ])
    with TemporaryDirectory() as temporary:
        output = Path(temporary) / 'step.out'
        geometry = Path(temporary) / 'step.xyz'
        output.write_text(
            'SETTING KB_GEOM_ENERGY = -40.43808112 AU\n'
            'Numerical gradient for KB_GEOM_ENERGY\n'
            'Total Energy           -40.43808112\n'
            'Atom dE/dx dE/dy dE/dz d2E/dx2 d2E/dy2 d2E/dz2\n' +
            ''.join(f'{i} {x} {y} {z} 0 0 0\n'
                    for i, (x, y, z) in enumerate(rows, 1)) +
            'Molpro calculation terminated\n')
        geometry.write_text('5\nMolpro geometry\nC 0 0 0\nH 0 0 1\n'
                            'H 1 0 0\nH 0 1 0\nH -1 0 0\n')
        forces = parse_numerical_gradient(output, geometry, atoms)
        np.testing.assert_allclose(forces, -rows[[0, 2, 3, 1, 4]] * Hartree / Bohr)
        assert forces[0, 0] > 0
        output.write_text(output.read_text().replace('Molpro calculation terminated',
                                                    'ERROR EXIT'))
        with pytest.raises(ValueError, match='terminate normally'):
            parse_numerical_gradient(output, geometry, atoms)


def test_ch4_l3_sella_calls_molpro_at_each_geometry_and_hashes_native_files():
    spec = ch4_spec()
    spec['tasks'] = [dict(spec['tasks'][1], geometry_from='initial')]
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        fake = root / 'molpro'
        fake.write_text(FAKE_MOLPRO.replace('{python}', sys.executable))
        fake.chmod(0o755)
        spec['tasks'][0]['profile']['command'] = str(fake)
        spec_file = root / 'spec.json'
        spec_file.write_text(json.dumps(spec))
        run_dir = prepare(spec_file, root / 'run')
        result = run_task(run_dir / 'tasks' / 'l3_geometry' / 'task.json')
        assert result['status'] == 'executed'
        assert result['details']['optimizer'] == 'sella'
        artifacts = result['artifacts']
        steps = sorted(name for name in artifacts if name.endswith('.inp'))
        assert len(steps) >= 2
        for inp in steps:
            stem = inp[:-4]
            assert {stem + suffix for suffix in ('.inp', '.out', '.log', '.xyz',
                                                  '.stdout', '.stderr')} <= set(artifacts)
        assert 'final.xyz' in artifacts
        assert 'optimization.traj' in artifacts
        assert advance(run_dir)['tasks']['l3_geometry']['status'] == 'complete'


def test_molpro_calculator_rejects_unsupported_spin():
    with TemporaryDirectory() as temporary:
        with pytest.raises(ValueError, match='neutral singlet'):
            Molpro(directory=temporary, label='ch3', method='CCSD(T)',
                   basis='cc-pVTZ', charge=0, mult=2)
