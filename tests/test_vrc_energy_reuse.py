from copy import deepcopy
import json
from types import SimpleNamespace

import pytest
from ase.io import read, write
from ase.units import Hartree

from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from kinbot.vrc_energy_reuse import accepted_energy, _same_subspace
from tests.anl_fixture import dispatch_spec
from tests.test_anl_dispatch import _mock_execution


@pytest.fixture
def point(tmp_path):
    spec = dispatch_spec()
    task = deepcopy(spec['tasks'][1])
    task.update(id='relaxed_scan', geometry_from='initial')
    task['profile'].update(method='caspt2-f12(6,4)', basis='avdz-f12',
        calculator_kwargs={'caspt2_shift':.1, 'caspt2_ipea':.25})
    spec['tasks'] = [task]
    graph = dispatch.prepare(spec, tmp_path/'point')
    _mock_execution(graph, 'relaxed_scan', geometry='final.xyz')
    d = graph/'tasks/relaxed_scan'
    a = read(d/'final.xyz')
    write(d/'central.xyz', a, format='xyz')
    (d/'relaxed_scan_step_0001.inp').write_text('geometry={\n'+(d/'central.xyz').read_text()+'}\n')
    (d/'relaxed_scan_step_0001.out').write_text('KB_GEOM_ENERGY = -40.123456789 AU\nMolpro calculation terminated\n')
    p = d/'execution.json'; e = json.loads(p.read_text())
    e['details'].update(energy_ev=-40.123456789*Hartree)
    for name in ('relaxed_scan_step_0001.inp', 'relaxed_scan_step_0001.out'):
        e['artifacts'][name] = file_sha256(d/name)
    atomic_json(p, e)
    dispatch.advance(graph)
    row = SimpleNamespace(data=dict(vrc_scan_role='relaxed', vrc_scan_dispatch=str(graph), energy=123.))
    component = dict(method='caspt2-f12(6,4)', basis='avdz-f12', multiplicity=None, state_model=None)
    species = SimpleNamespace(atom=a.get_chemical_symbols(), geom=a.positions, charge=0, mult=1)
    return graph, row, component, species


def test_verified_reuse_uses_native_energy_not_legacy_row_energy(point):
    graph, row, component, species = point
    proof = accepted_energy(row, component, species)
    assert proof['energy_hartree'] == pytest.approx(-40.123456789)
    assert proof['source_kind'] == 'accepted_relaxed_scan'
    assert proof['execution_sha256'] == file_sha256(graph/'tasks/relaxed_scan/execution.json')
    assert proof['geometry_sha256'] == dispatch._load(graph)[2]['tasks']['relaxed_scan']['final_geometry_sha256']


@pytest.mark.parametrize('change', ['geometry', 'basis', 'method', 'shift', 'ipea', 'multiplicity', 'charge', 'model'])
def test_different_physical_component_requires_its_own_energy(point, change):
    _, row, component, species = point
    if change == 'geometry': species.geom[1,0] += .01
    elif change == 'charge': species.charge = 1
    else: component[change] = {'basis':'avtz-f12', 'method':'caspt2(6,4)', 'shift':.2,
        'ipea':0., 'multiplicity':3, 'model':{'selected_roots':{'1':2}}}[change]
    if change == 'model': component['state_model'] = component.pop('model')
    assert accepted_energy(row, component, species) is None


def test_frozen_row_reuses_only_when_geometry_is_actually_identical(point):
    _, row, component, species = point
    row.data['vrc_scan_role'] = 'frozen'
    assert accepted_energy(row, component, species)
    species.geom[1,0] += .02
    assert accepted_energy(row, component, species) is None


def test_changed_native_output_cannot_supply_cached_energy(point):
    graph, row, component, species = point
    (graph/'tasks/relaxed_scan/relaxed_scan_step_0001.out').write_text('changed')
    with pytest.raises(ValueError, match='source changed'):
        accepted_energy(row, component, species)


def test_orbital_selection_cannot_be_dropped_from_model_comparison():
    old = {'orbital_selection': {'mode':'avas', 'centers':[{'atom':1, 'orbitals':['2p']}]}}
    new = deepcopy(old); new['orbital_selection']['centers'][0]['atom'] = 2
    assert _same_subspace(old, old)
    assert not _same_subspace(old, new)
