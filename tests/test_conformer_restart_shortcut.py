"""A finished conformer search is skipped on restart for every QC code.

The restart shortcut looked for conf/<name>_low.log and <last job>.log, the
Gaussian names, so Q-Chem (.out) and the Sella-driven codes (_sella.log)
repeated the whole conformer sweep on every restart, and the file removal
that forces a new search after a lower rotor point never applied to them.
It also derived the number of conformers from grid ** dihedrals * ring
conformers, which undercounts a search started from the original geometry
plus its ring conformers, and from the first random batch only; conformers
beyond that count, possibly the lowest one, were dropped on restart. The
count now comes from the finished jobs in kinbot.db.
"""
import json
from pathlib import Path

import pytest
from ase.build import molecule

from kinbot.conformers import Conformers
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


EXTENSIONS = {'gauss': '.log', 'qchem': '.out', 'fc': '_sella.log', 'orca': '_sella.log'}


@pytest.fixture
def search(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('conf').mkdir()
    Path('input.json').write_text(json.dumps({'barrier_threshold': 100., 'queuing': 'local'}))
    par = Parameters('input.json').par
    atoms = molecule('H2O2')  # one torsion: conf_grid ** 1 = 3 grid conformers
    point = StationaryPoint('hydrogen_peroxide', 0, 1,
                            atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    point.characterize()
    search = Conformers(point, par, QuantumChemistry(par))
    assert len(point.conf_dihed) == 1
    search.atoms = atoms

    def finish(count):
        """Record `count` finished conformer jobs, as a completed search leaves them."""
        for i in range(count):
            search.db.write(atoms, name=search.get_job_name(i),
                            data={'energy': -1. - i, 'zpe': .02, 'frequencies': [100.] * 6, 'status': 'normal'})
        search.last = search.get_job_name(count - 1)
    search.finish = finish
    finish(search.grid)
    return search


@pytest.mark.parametrize('code', sorted(EXTENSIONS))
def test_finished_search_is_skipped_for_every_code(search, code):
    search.qc.qc = code
    ext = EXTENSIONS[code]
    assert search.low_log() == f'conf/{search.get_name()}_low{ext}'
    Path(search.last + ext).write_text('done\n')
    Path(search.low_log()).write_text('done\n')
    assert search.generate_conformers(0, search.species.geom) == 1
    assert search.conf == search.grid


def test_count_comes_from_the_finished_jobs_not_from_the_grid(search):
    # A ring species starts from its original geometry and every ring
    # conformer, and random sampling adds a batch per geometry: more jobs
    # than grid ** dihedrals. All of them must survive a restart.
    search.finish(12)
    Path(search.last + '.log').write_text('done\n')
    Path(search.low_log()).write_text('done\n')
    assert search.generate_conformers(0, search.species.geom) == 1
    assert search.conf == 12


def test_unfinished_search_is_not_skipped(search):
    Path(search.last + '.log').write_text('done\n')
    assert search.completed_conformers() == 0  # no conf/<name>_low file


def test_saddle_imaginary_mode_follows_the_low_copy(search, monkeypatch):
    search.qc.qc = 'fc'
    monkeypatch.setattr(search.qc, 'qc_conf', lambda *args, **kwargs: None)
    Path(search.last + '_sella.log').write_text('done\n')
    Path(search.last + '_vib.xyz').write_text('4\nimaginary mode\n')
    search.copy_low_log(search.last)
    assert Path(search.low_vib()).read_text().startswith('4')
    search.add_new_conf_from_hir(search.species.geom)
    assert not Path(search.low_vib()).exists()


@pytest.mark.parametrize('code', sorted(EXTENSIONS))
def test_copy_and_removal_use_the_code_extension(search, code, monkeypatch):
    search.qc.qc = code
    monkeypatch.setattr(search.qc, 'qc_conf', lambda *args, **kwargs: None)  # no submission
    job = search.last
    Path(job + EXTENSIONS[code]).write_text('done\n')
    search.copy_low_log(job)
    assert Path(search.low_log()).read_text() == 'done\n'
    search.add_new_conf_from_hir(search.species.geom)
    assert not Path(search.low_log()).exists()


def test_missing_output_warns_instead_of_failing(search, caplog):
    search.qc.qc = 'fc'
    search.copy_low_log(search.last)
    assert not Path(search.low_log()).exists()
    assert 'restart will check every conformer again' in caplog.text
