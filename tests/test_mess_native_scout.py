"""Native MESS regressions from saved scouting evidence, without QC retries.

This checks serialization and native solver execution, not the accuracy of the
historical kinetic model or the numerical reliability of every pressure table.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

import pytest

from kinbot.mess import MESS
from tests.test_mess_pes_artifacts import at, write_worker
from tests.test_mess_population_identity import stage_collision_network, assemble_collision_network


FIXTURE = Path(__file__).parent / 'reference/mess_scout_h2_geometry.inp'


def test_direct_and_pes_headers_preserve_short_bond_geometry_guard(tmp_path):
    cases = stage_collision_network(tmp_path)
    for case in cases:
        write_worker(case)
    combined = assemble_collision_network(cases)
    case = cases[0]
    with at(case.worker):
        writer = MESS(dict(case.par, pes=0), case.parent)
        writer.write_input(None)
        direct = Path('me/mess_0000.inp').read_text()
    for text in (direct, *combined.values()):
        header = text.split('\nModel\n', 1)[0]
        values = re.findall(r'^AtomDistanceMin\[angstrom\]\s+(\S+)', header, re.M)
        assert len(values) == 1 and float(values[0]) == .32


def test_native_complete_scout_input_accepts_unchanged_h2_geometry(tmp_path):
    executable = os.environ.get('KINBOT_TEST_MESS')
    if not executable:
        pytest.skip('Set KINBOT_TEST_MESS to check the saved complete scout network')
    from kinbot import kb_path
    original = FIXTURE.read_text()
    metadata = json.loads(FIXTURE.with_suffix('.json').read_text())
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == metadata['sha256']
    template = Path(kb_path, 'tpl/mess_header.tpl').read_text()
    # Use the real common writer setting. Keep the entire native model,
    # temperature/pressure grids and numerical method from the captured input.
    guard = re.search(r'^AtomDistanceMin\[angstrom\].+$', template, re.M)
    assert guard is not None
    text = guard[0] + '\n' + original
    assert text.split('\n', 1)[1] == original
    path = tmp_path / 'scout.inp'
    path.write_text(text)
    result = subprocess.run([executable, path.name], cwd=tmp_path,
        env=dict(os.environ, OMP_NUM_THREADS='1'), text=True,
        capture_output=True, timeout=120)
    (tmp_path / 'console.log').write_text(result.stdout + '\n' + result.stderr)
    assert result.returncode == 0, result.stderr[-4000:]
    log = (tmp_path / 'scout.log').read_text()
    output = (tmp_path / 'scout.out').read_text()
    assert 'rate calculation done' in log
    assert 'High Pressure Rate Coefficients' in output
    micro = re.search(r'^MicroRateOutput\s+(\S+)', text, re.M)[1]
    assert (tmp_path / micro).stat().st_size > 0
    # Native execution is necessary but not sufficient for rate validation.
    # Low-temperature pressure-dependent matrices still need a precision and
    # reduction check before they can be used to prune a production network.
