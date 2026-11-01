"""Worker artifacts must cover the route actually selected by PES assembly."""
from contextlib import contextmanager
from copy import deepcopy
import json
import os
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from ase import Atoms
from ase.db import connect
import pytest

from kinbot import constants, pes, postprocess
from kinbot.execution import file_sha256
from kinbot.mess import MESS
from kinbot.parameters import Parameters
from kinbot.run_format import ensure_current_run
from kinbot.species_routing import routing_name
from kinbot.stereo_routing import guard_well_job


@contextmanager
def at(directory):
    old = Path.cwd()
    os.chdir(directory)
    try:
        yield
    finally:
        os.chdir(old)


def stage_artifact_handoff(root, *, reverse=False, distinct=False, mc=False, tied=True):
    # Reuse structural fixtures; assigned energies test serialization/selection,
    # not the chemistry of these structures.
    from tests.test_complex_mess import TestComplexMESS
    parent, reactions = TestComplexMESS().reactions()
    for reaction in reactions:
        reaction.instance_name = routing_name(parent)+'_test_'+reaction.instance_name
        reaction.ts.name = reaction.instance_name
        reaction.irc_prod.name = reaction.instance_name+'_IRC_F_prod'
    if tied:
        reactions[1].ts.energy = reactions[0].ts.energy
    if distinct:
        for reaction in reactions:
            reaction.ts.stereopath_id = reaction.instance_name
            reaction.ts.stereopath_metadata['id'] = reaction.instance_name
    if reverse:
        reactions.reverse()
    parent.reac_obj = reactions
    parent.reac_inst = [None] * len(reactions)
    parent.reac_name = [r.instance_name for r in reactions]
    settings = dict(barrier_threshold=100., high_level=1, rotor_scan=0,
        multi_conf_tst=int(mc), conformer_search=1, me=2, uq=0,
        epsilon=100., sigma=3., queuing='local', smiles='C', pes=1)
    root.mkdir(parents=True, exist_ok=True)
    source = root/'artifact_test_input.json'
    source.write_text(json.dumps(settings))
    par = Parameters(str(source), show_warnings=False).par
    worker = root/routing_name(parent)
    worker.mkdir()
    ensure_current_run(worker, create=True)
    (worker/'me').mkdir()
    db = connect(str(worker/'kinbot.db'))
    points = [(parent, routing_name(parent)+'_well')]
    points.extend((p, routing_name(p)+'_well') for p in reactions[0].products)
    for reaction in reactions:
        points.extend([(reaction.ts, reaction.instance_name),
                       (reaction.irc_prod_opt.species, reaction.irc_prod.name)])
    for point, name in points:
        name += '_high'
        if 'IRC' not in name and not point.wellorts:
            guard_well_job(SimpleNamespace(db=db), point, point.geom, name)
        db.write(Atoms(point.atom, positions=point.geom), name=name,
                 data=dict(status='normal', energy=point.energy/constants.EVtoHARTREE,
                           zpe=point.zpe, frequencies=point.freq))
    return SimpleNamespace(root=root, worker=worker, par=par, parent=parent,
                           reactions=reactions, db_hash=file_sha256(worker/'kinbot.db'))


def write_worker(t):
    with at(t.worker):
        MESS(t.par, t.parent).write_input(None)
        postprocess.create_summary_file(t.parent, SimpleNamespace(qc='fc'), t.par)


def assemble(t, *, skip=False):
    # Keep the actual summary parser, DB lookup, route selection, complex
    # selection, MESS fragments, finalization and network splitting. Only
    # unrelated plotting/L3/RotdPy exports are suppressed.
    with at(t.root), patch('kinbot.pes.copy_from_kinbot'), \
            patch('kinbot.pes.create_pesviewer_input'), \
            patch('kinbot.pes.create_rotdpy_inputs'), \
            patch('kinbot.pes.t1_analysis'), \
            patch('kinbot.pes.get_l3energy', return_value=(0, -1)):
        pes.postprocess(dict(t.par, me_skip_vdW=int(skip)), [t.worker.name],
                        'all', [], t.parent.mass)
    return (t.root/'me/mess_0000.inp').read_text()


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('distinct,tied', [(False, True), (False, False), (True, True)])
@pytest.mark.parametrize('mc', [False, True])
@pytest.mark.parametrize('mode', ['off', 'auto'])
def test_worker_summary_to_pes_assembly_exports_every_accepted_complex(
        tmp_path, reverse, distinct, tied, mc, mode):
    t = stage_artifact_handoff(tmp_path, reverse=reverse, distinct=distinct, mc=mc, tied=tied)
    t.par['rotdpy_mess_mode'] = mode
    original = deepcopy([(r.ts.energy, r.ts.zpe, r.ts.geom.tolist(),
                          r.irc_prod_opt.species.energy, r.irc_prod_opt.species.freq)
                         for r in t.reactions])
    previous = None
    for _ in range(2):
        write_worker(t)
        # Exercise assembly before inspecting artifacts: the old code must
        # reproduce the observed FileNotFoundError for the equal-energy tie.
        text = assemble(t)
        assert text.count('Core PhaseSpaceTheory') == 1
        assert text.count('ImaginaryFrequency') == (2 if distinct else 1)
        assert sum(line.strip().startswith('Well ') for line in text.splitlines()) == 2
        assert '{name}' not in text and '{wellname}' not in text
        for reaction in t.reactions:
            assert (t.worker/(reaction.irc_prod.name+'_0000.mess')).is_file()
        artifacts = {p.name: p.read_bytes() for p in t.worker.glob('*.mess')}
        omitted = assemble(t, skip=True)
        assert 'Core PhaseSpaceTheory' not in omitted
        assert omitted.count('ImaginaryFrequency') == (2 if distinct else 1)
        assert sum(line.strip().startswith('Well ') for line in omitted.splitlines()) == 1
        assert artifacts == {p.name: p.read_bytes() for p in t.worker.glob('*.mess')}
        assert t.parent.reac_ts_done == [-1, -1]
        assert file_sha256(t.worker/'kinbot.db') == t.db_hash
        assert original == [(r.ts.energy, r.ts.zpe, r.ts.geom.tolist(),
                             r.irc_prod_opt.species.energy, r.irc_prod_opt.species.freq)
                            for r in t.reactions]
        if previous is not None:
            assert (text, omitted, artifacts) == previous
        previous = text, omitted, artifacts


def test_legacy_single_exit_fallback_and_wrong_reaction_context(tmp_path):
    from kinbot.species_routing import mess_filename
    t = stage_artifact_handoff(tmp_path, tied=False, reverse=True)
    write_worker(t)
    products = '_'.join(sorted(routing_name(p) for p in t.reactions[0].products))
    # The first route has the lower TS, but the shared legacy product file
    # belongs to the last route. Do not silently reinterpret its ROTD context.
    source = t.worker/mess_filename(products, 0, reaction=t.reactions[0].instance_name)
    good = source.read_bytes()
    source.unlink()
    with pytest.raises(ValueError, match='context does not match'):
        assemble(t)
    (t.worker/mess_filename(products, 0)).write_bytes(good)
    assert assemble(t).count('Core PhaseSpaceTheory') == 1
