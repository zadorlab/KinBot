"""Exercise queued worker paths and completed endpoint records without QC jobs."""
import json
import os
from pathlib import Path
import subprocess
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.db import connect

from kinbot import kb_path, reader_qchem, utils
from kinbot.ase_modules.calculators import qchem
from kinbot.bimol import load_selected_species_properties
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.reaction_generator import ReactionGenerator
from kinbot.stationary_pt import StationaryPoint


@pytest.mark.parametrize('queue', ['slurm', 'pbs'])
@pytest.mark.parametrize('job', ['bimol/ch4_oh/search', 'product_spin/ch3_state'])
def test_orca_nested_worker_starts_from_calculation_directory(
        tmp_path, monkeypatch, queue, job):
    root = tmp_path / 'calculation directory'
    root.mkdir()
    monkeypatch.chdir(root)
    scratch = tmp_path / 'scratch'
    scratch.mkdir()
    Path('header').write_text(f'#!/bin/sh\ncd {scratch}\n')
    Path('input.json').write_text(json.dumps(dict(barrier_threshold=50.,
        qc='orca', qc_command='orca', queuing=queue, queue_template='header',
        queue_job_limit=0)))
    par = Parameters('input.json', show_warnings=False).par
    par['bimol'] = 1
    qc = QuantumChemistry(par).l1
    monkeypatch.setattr(qc, 'check_qc', lambda name: 0)
    Path(job).parent.mkdir(parents=True)
    Path(job + '.py').write_text(
        "from pathlib import Path\nPath('worker-cwd').write_text(str(Path.cwd()))\n")
    # Only the scheduler submission is simulated. Execute the generated shell
    # script with a harmless Python worker to check its actual path and cwd.
    with pytest.MonkeyPatch.context() as patch:
        process = Mock()
        process.communicate.return_value = (
            b'Submitted batch job 123\n' if queue == 'slurm' else b'123.server\n', b'')
        patch.setattr('kinbot.qc.subprocess.Popen', Mock(return_value=process))
        qc.submit_qc(job, 1)
    script = Path(job + ('.sbatch' if queue == 'slurm' else '.pbs'))
    source = script.read_text()
    assert f'python {job}.py' in source
    assert 'rm -rf' not in source
    subprocess.run(['sh', str(script.resolve())], cwd=scratch,
                   env=dict(os.environ, PBS_O_WORKDIR=str(root)), check=True)
    assert Path('worker-cwd').read_text() == str(root)


@pytest.mark.parametrize('converged', [None, False, True])
def test_l1_fairchem_complex_classification_checks_saved_convergence(
        tmp_path, monkeypatch, converged):
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps(dict(barrier_threshold=50.,
        high_level=0, rotor_scan=0, conformer_search=0, L3_calc=0,
        l1_profile=dict(calculator='fairchem', model_path='test-model'))))
    par = Parameters('input.json', show_warnings=False).par
    par['bimol'] = 1
    qc = QuantumChemistry(par)
    monkeypatch.setattr(QuantumChemistry, 'submit_qc',
                        Mock(side_effect=AssertionError('No additional jobs expected')))
    endpoints = []
    for direction in ('F', 'R'):
        point = StationaryPoint('ts_IRC_' + direction + '_prod', 0, 1,
            atom=['H', 'H'], geom=np.array([[0., 0., 0.], [.74, 0., 0.]]))
        point.characterize()
        data = dict(status='normal', energy=-28., zpe=.01, frequencies=[4400.])
        if converged is not None:
            data['endpoint_converged'] = converged
        qc.db.write(Atoms(point.atom, positions=point.geom), name=point.name, data=data)
        Path(point.name + '_sella.log').write_text('done\n')
        endpoints.append(point)
    fragment = SimpleNamespace(energy=-.5, zpe=0.)
    reaction = SimpleNamespace(instance_name='abstraction',
        irc_reactant=endpoints[0], irc_prod=endpoints[1],
        prod_opt=[SimpleNamespace(species=fragment)] * 2)
    generator = ReactionGenerator(SimpleNamespace(fragA=fragment, fragB=fragment),
                                  par, qc, 'input.json')
    generator._classify_bimol_complexes(reaction)
    assert reaction.do_vdW is (converged is True)
    assert reaction.do_reac_vdW is (converged is True)
    if converged:
        for opt in (reaction.irc_prod_opt, reaction.irc_reactant_opt):
            load_selected_species_properties(opt, qc)
    generator._finalize_bimol_complexes(reaction)
    assert reaction.do_vdW is (converged is True)
    assert reaction.do_reac_vdW is (converged is True)


@pytest.mark.parametrize('bimol', [False, True])
@pytest.mark.parametrize('failure', [None, 'opt', 'freq', 'missing_zpe', 'missing_freq', 'partial_irc'])
def test_qchem_endpoint_worker_publishes_harmonic_properties(
        tmp_path, monkeypatch, bimol, failure):
    monkeypatch.chdir(tmp_path)
    calls = []
    optimized = np.array([[0., 0., 0.], [.75, 0., 0.]])

    class FakeQChem(Calculator):
        implemented_properties = ['energy']

        def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            stage = self.parameters['jobtype']
            calls.append((stage, self.label, atoms.positions.copy()))
            output = 'Frequency: 4400.0\nZero point vibrational energy: 6.29 kcal/mol\n'
            if failure == 'missing_zpe':
                output = 'Frequency: 4400.0\n'
            if failure == 'missing_freq':
                output = 'Zero point vibrational energy: 6.29 kcal/mol\n'
            Path(self.label + '.out').write_text(output)
            if ((stage == failure and self.label.startswith('irc_F_prod'))
                    or (failure == 'partial_irc' and stage == 'rpath')):
                raise RuntimeError('Simulated stopped calculation')
            self.results['energy'] = -28.

    monkeypatch.setattr(qchem, 'QChem', FakeQChem)
    monkeypatch.setattr(utils, 'iowait', lambda *args: None)
    monkeypatch.setattr(reader_qchem, 'read_geom', lambda *args, **kwargs: optimized.copy())
    script = Path(kb_path, 'tpl/ase_qchem_irc.tpl.py').read_text().format(
        working_dir=str(tmp_path), label='irc_F', atom=['H', 'H'],
        geom=[[0., 0., 0.], [1., 0., 0.]],
        kwargs=dict(label='irc_F', jobtype='freq'),
        prod_kwargs=dict(label='irc_F_prod', jobtype='opt', method='HF', basis='sto-3g'),
        qc_command='unused', ppn=1, bimol=bimol)
    exec(compile(script, 'irc_F.py', 'exec'), {})
    db = connect('kinbot.db')
    rows = list(db.select(name='irc_F_prod'))
    assert len(rows) == 1
    row = rows[0]
    endpoint_frequency_calls = [call for call in calls if call[1] == 'irc_F_prod_freq']
    if bimol and failure in ('opt', 'freq', 'missing_zpe', 'missing_freq'):
        assert row.data.status == 'error'
    elif bimol:
        assert row.data.status == 'normal'
        assert row.data.energy == -28.
        assert list(row.data.frequencies) == [4400.]
        assert row.data.zpe > 0.
        assert len(endpoint_frequency_calls) == 1
        np.testing.assert_allclose(endpoint_frequency_calls[0][2], optimized)
        assert Path('irc_F_prod_freq.out').read_text().endswith('done\n')
    else:
        assert not endpoint_frequency_calls
        assert 'frequencies' not in row.data and 'zpe' not in row.data
    if failure == 'opt':
        assert not endpoint_frequency_calls
    assert Path('irc_F_prod.out').read_text().endswith('done\n')
