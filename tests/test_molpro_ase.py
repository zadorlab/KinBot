"""Molpro ASE/Sella contract tests; the Molpro executable is simulated.

The XYZGRAD rows are synthetic in the format documented by Molpro. The real
program/version remains an offsite acceptance gate.
"""

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

from ase import Atoms
from ase.units import Hartree
import numpy as np
import pytest

from examples.anl.ch4_dispatch import ch4_spec
from kinbot.anl.dispatch import advance, prepare, run_task
from kinbot.ase_modules.calculators.molpro import (
    Molpro, parse_output, parse_xyzgrad, render_input,
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
grad = [[-2*k*27.211386245988*(value - goal)
         for value, goal in zip(xyz, goals)]
        for (_, xyz), goals in zip(rows, target)]
stem = inp.stem
Path(stem + '.out').write_text(
    f' SETTING KB_GEOM_ENERGY = {energy:.14f} AU\n'
    f' SETTING KB_GEOM_ENERGY = {energy + 0.1:.14f} AU\n'
    ' Molpro calculation terminated\n')
Path(stem + '.log').write_text('Molpro gradient log\n')
Path(stem + '.xyz').write_text(
    str(n) + '\nMolpro XYZGRAD forces (-eV/Angstrom)\n' +
    ''.join(f'{symbol} {x:.12f} {y:.12f} {z:.12f} '
            f'{6 if symbol == "C" else 1} '
            f'{force[0]:.12f} {force[1]:.12f} {force[2]:.12f}\n'
            for (symbol, (x, y, z)), force in zip(rows, grad)))
'''


def test_documented_molpro_force_input_and_parser_contract():
    atoms = Atoms('CH2', positions=[[0, 0, 0], [0, 0, 1], [0, 1, 0]])
    deck = render_input(atoms, basis='cc-pVTZ', gradient_name='step.xyz')
    assert 'set,charge=0\nset,spin=0\ngthresh,energy=1.d-9\n' in deck
    assert ('rhf\nccsd(t)\nkb_geom_energy=energy\n'
            'forces,numerical,variable=kb_geom_energy,startcmd=rhf\n'
            'put,xyzgrad,step.xyz') in deck
    assert 'optg' not in deck.lower()
    assert '3\nKinBot ASE geometry (angstrom)\nC' in deck
    with TemporaryDirectory() as temporary:
        output = Path(temporary) / 'step.out'
        output.write_text(' SETTING KB_GEOM_ENERGY        =       -40.123456789  AU\n'
                          ' SETTING KB_GEOM_ENERGY = -39.0 AU\n'
                          ' Molpro calculation terminated\n')
        assert parse_output(output) == pytest.approx(-40.123456789 * Hartree)
        gradient = Path(temporary) / 'step.xyz'
        gradient.write_text(
            '3\nMolpro XYZGRAD\n'
            'H 0 1 0 1 0.1 0.2 0.3\n'
            'C 0 0 0 6 -1 -2 -3\n'
            'H 0 0 1 1 0.4 0.5 0.6\n')
        np.testing.assert_allclose(parse_xyzgrad(gradient, atoms),
                                   [[-1, -2, -3], [.4, .5, .6], [.1, .2, .3]])
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
