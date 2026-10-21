"""A final IRC cannot reuse another level's path at unchanged TS coordinates."""
import copy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms

from kinbot.bimol import final_irc_job_name
from kinbot.irc import IRC
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from tests.test_bimol_sequential_poses import reactants


@pytest.mark.parametrize('initial_high_level', [0, 1])
def test_changing_selected_level_preserves_old_irc_and_requests_new_path(
        tmp_path, monkeypatch, initial_high_level):
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps(dict(
        barrier_threshold=100., high_level=initial_high_level,
        conformer_search=1, queuing='local',
        l1_profile=dict(calculator='gaussian', method='HF', basis='STO-3G',
                        optimizer='sella'),
        l2_profile=dict(calculator='gaussian', method='B3LYP', basis='6-31G(d)',
                        optimizer='sella'))))
    par = Parameters('input.json', show_warnings=False).par
    par['bimol'] = 1
    qc = QuantumChemistry(par)
    species = reactants()
    ts = copy.copy(species)
    level = 'l2' if initial_high_level else 'l1'
    ts.source_job = 'abstraction_high' if initial_high_level else 'abstraction'
    qc._record(ts.source_job, level)
    reaction = SimpleNamespace(qc=qc, par=par, species=species, ts=ts,
                               instance_name='abstraction')
    reaction.final_irc_name = final_irc_job_name(reaction)
    old_name = reaction.final_irc_name
    for direction in ('F', 'R'):
        job = f'{old_name}_IRC_{direction}'
        qc._record(job, level)
        Path(job + '_sella.log').write_text('done\n')
        qc.db.write(Atoms(species.atom, positions=species.geom), name=job,
                    data=dict(status='normal', energy=-10.))
    assert IRC(reaction, par).check_irc() == ['normal', 'normal']

    new_par = dict(par, high_level=1 - initial_high_level)
    restarted = QuantumChemistry(new_par)
    new_level = 'l2' if new_par['high_level'] else 'l1'
    ts.source_job = 'abstraction_high' if new_par['high_level'] else 'abstraction'
    restarted._record(ts.source_job, new_level)
    reaction.qc, reaction.par = restarted, new_par
    reaction.final_irc_name = final_irc_job_name(reaction)
    assert reaction.final_irc_name != old_name
    assert IRC(reaction, new_par).qc.theory_level == new_level
    assert IRC(reaction, new_par).check_irc() == [0, 0]
    assert restarted.db.count(name=f'{old_name}_IRC_F') == 1
    assert Path(f'{old_name}_IRC_F_sella.log').read_text() == 'done\n'

    # The selected calculation's recorded level takes precedence over the
    # global switch, just as it does when IRC chooses its calculator.
    reaction.par = par
    assert f'_final_{new_level}_' in final_irc_job_name(reaction)


def test_old_geometry_only_name_is_not_reused():
    ts = SimpleNamespace(geom=np.zeros((3, 3)))
    reaction = SimpleNamespace(ts=ts, par={'high_level': 0}, instance_name='abstraction')
    old_hash = hashlib.sha1(np.round(ts.geom, 6).tobytes()).hexdigest()[:10]
    assert final_irc_job_name(reaction) != f'abstraction_final_{old_hash}'
    assert final_irc_job_name(reaction) == f'abstraction_final_l1_{old_hash}'
