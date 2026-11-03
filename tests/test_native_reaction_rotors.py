"""Saved TS/IRC evidence exercises the production rotor writer, without QC."""
import ast
import json
from pathlib import Path
from unittest.mock import Mock

import numpy as np

from kinbot import frequencies
from kinbot.hindered_rotors import HIR
from kinbot.optimize import Optimize
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


DATA = json.loads((Path(__file__).parent/'reference/rotor_reaction_endpoints.json').read_text())


def native_ts(order=None):
    order = np.arange(6) if order is None else np.asarray(order)
    points = {}
    for label, entry in DATA['points'].items():
        p = StationaryPoint(label, 0, 2,
            atom=np.asarray(entry['atoms'])[order],
            geom=np.asarray(entry['geometry'])[order],
            wellorts=int(label.startswith('ts')))
        p.characterize()
        points[label] = p
    ts, parent, endpoint = [points[k] for k in ('ts_l1', 'reactant', 'product')]
    # Same TS graph construction as ReactionGenerator stage 3. In particular,
    # keep both H-transfer bonds; the near-saddle distance graph misses H6.
    ts.bond = np.maximum(ts.bond, np.maximum(parent.bond, endpoint.bonds[0]))
    ts.reac_bond = (endpoint.bond > 0).astype(int) - parent.bond01
    from kinbot.stereo_identity import canonical_identity
    from kinbot.reaction_path import prepare_stereopath
    ts.optical_reference = canonical_identity(parent)
    ts.find_cycle()
    ts.find_conf_dihedral()
    prepare_stereopath(ts, parent, endpoint)
    ts.geom = points['ts_l2'].geom.copy()
    return ts


def test_native_spectator_sulfur_axis_selected_but_reaction_loop_excluded():
    rng = np.random.default_rng(31)
    for order in [np.arange(6), np.arange(6)[::-1]] + [rng.permutation(6) for _ in range(12)]:
        ts = native_ts(order)
        assert len(ts.dihed) == 1
        assert set(order[ts.dihed[0][1:3]]) == {1, 2}  # S-S, not C-S
        assert {0, 1, 5} in [set(order[cycle]) for cycle in ts.cycle_chain]
        left, right = frequencies.partition(ts, ts.dihed[0], ts.natom)
        assert {frozenset(order[left]), frozenset(order[right])} == {
            frozenset([0, 1, 4, 5]), frozenset([2, 3])}


def test_native_axis_reaches_real_hir_writer_and_live_points_are_adopted(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('hir').mkdir()
    Path('input.json').write_text(json.dumps(dict(barrier_threshold=150.,
        queuing='local', rotor_scan=1, conformer_search=1, high_level=1,
        use_sella=True, nrotation=12)))
    par = Parameters('input.json', show_warnings=False).par
    qc = QuantumChemistry(par)
    running, submissions = set(), []
    qc.check_qc = Mock(side_effect=lambda job: 'running' if job in running else 0)
    def submitted(job, nproc):
        running.add(job)
        submissions.append(job)
    qc.submit_qc = Mock(side_effect=submitted)
    ts = native_ts()
    opt = Optimize(ts, par, qc)  # real re-characterization must retain the axis
    assert opt.species.dihed == [[5, 1, 2, 3]]
    hir = HIR(ts, qc, par)
    original = ts.geom.copy()
    hir.generate_hir_geoms(original.copy(), False)
    assert len(submissions) == 12
    scripts = {job: Path(job+'.py').read_bytes() for job in submissions}
    coordinates = []
    for job in submissions:
        source = scripts[job].decode()
        tree = ast.parse(source)
        atoms = next(node.value for node in tree.body if isinstance(node, ast.Assign)
                     and isinstance(node.value, ast.Call)
                     and getattr(node.value.func, 'id', None) == 'Atoms')
        xyz = np.asarray(next(ast.literal_eval(k.value) for k in atoms.keywords
                              if k.arg == 'positions'))
        coordinates.append(xyz)
        assert 'order=1' in source  # constrained saddle, not minimum optimization
        assert '[6, 2, 3, 4]' in source
        for group in frequencies.partition(ts, ts.dihed[0], ts.natom):
            old, new = original[group], xyz[group]
            np.testing.assert_allclose(np.linalg.norm(old[:, None]-old[None, :], axis=-1),
                                       np.linalg.norm(new[:, None]-new[None, :], axis=-1), atol=1.e-7)
    assert not np.allclose(coordinates[0], coordinates[6])
    # A repeated controller construction adopts the exact saved/live inputs.
    # It cannot overwrite them or submit the same point again.
    HIR(ts, qc, par).generate_hir_geoms(original.copy(), False)
    assert len(submissions) == 12
    assert scripts == {job: Path(job+'.py').read_bytes() for job in submissions}
    np.testing.assert_array_equal(ts.geom, original)


def test_hir_selection_does_not_establish_optical_coverage():
    from kinbot.optical import evaluate_optical
    ts = native_ts()
    assert ts.dihed
    result = evaluate_optical(ts)
    assert result.get('states_covered_by_hir', 1) == 1
