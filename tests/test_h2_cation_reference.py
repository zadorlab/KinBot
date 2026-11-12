"""Bound H2+ stays molecular, and its UMA correction keeps one energy reference."""
import json
import logging
import pickle
import runpy
import sys
from types import SimpleNamespace, ModuleType

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from kinbot import constants, symmetry
from kinbot.fairchem_utils import h2_cation_reference, uma_model
from kinbot.frequencies import get_frequencies
from kinbot.stationary_pt import StationaryPoint


class NeutralHydrogen(Calculator):
    """Analytic neutral-H2 potential; raises if asked to calculate an ion."""
    implemented_properties = ['energy', 'forces']
    task_name = 'omol'

    def __init__(self, offset=-30., **kwargs):
        super().__init__(**kwargs)
        self.offset = offset

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        assert atoms.info['charge'] == 0
        assert atoms.info['spin'] == 1
        displacement = atoms.positions[1] - atoms.positions[0]
        distance = np.linalg.norm(displacement)
        force = 10. * (distance - .74) * displacement / distance
        self.results = {'energy': self.offset + 5. * (distance - .74)**2,
                        'forces': np.array([force, -force])}


def hydrogen(distance, charge=1, mult=2):
    point = StationaryPoint('hydrogen', charge, mult, atom=['H', 'H'],
        geom=np.array([[1., 2., 3.], [1., 2., 3.+distance]]))
    point.characterize()
    return point


@pytest.mark.parametrize('distance,charge,mult,bonded', [
    (1.1142815143, 1, 2, True), (1.29, 1, 2, True), (1.31, 1, 2, False),
    (1.1142815143, 0, 1, False), (.74, 0, 1, True),
    (1.1142815143, -1, 2, False)])
def test_only_isolated_doublet_h2_cation_has_the_longer_cutoff(distance, charge, mult, bonded):
    point = hydrogen(distance, charge, mult)
    assert bool(point.bond[0, 1]) == bonded
    fragments, _ = point.start_multi_molecular(vary_charge=True)
    if bonded:
        assert len(fragments) == 1
        assert fragments[0].charge == charge
        assert fragments[0].mult == mult


def test_h2_cation_cutoff_does_not_join_hydrogens_inside_larger_ions():
    point = StationaryPoint('three hydrogens', 1, 2, atom=['H']*3,
        geom=np.array([[0., 0., 0.], [1.114, 0., 0.], [5., 0., 0.]]))
    point.characterize()
    assert point.bond[0, 1] == 0


@pytest.mark.parametrize('offset', [-30., -31.])
def test_uma_h2_cation_energy_geometry_and_hessian_agree(tmp_path, monkeypatch, offset):
    monkeypatch.chdir(tmp_path)
    point = hydrogen(.544)
    mol = Atoms('H2', positions=point.geom)
    mol.info.update(charge=1, spin=2)
    mol.calc = NeutralHydrogen(offset)
    center = mol.positions.mean(axis=0)
    data = h2_cation_reference(mol, 'ion', 'uma-m-1p1', .005, 50)
    reference = data['h2_cation_reference']
    assert data['status'] == 'normal'
    np.testing.assert_allclose(mol.positions.mean(axis=0), center)
    assert mol.get_distance(0, 1) == pytest.approx(1.058)
    assert data['frequencies'] == pytest.approx([2321.7])
    assert data['zpe'] == pytest.approx(.5*2321.7/constants.AUtoCM)
    delta = (data['energy'] + data['zpe']/constants.EVtoHARTREE
             - reference['neutral_h2_energy_ev']
             - reference['neutral_h2_zpe_hartree']/constants.EVtoHARTREE)
    assert delta == pytest.approx(15.42593)
    assert reference['neutral_h2_energy_ev'] == offset
    point.geom = mol.positions.copy()
    point.characterize()
    frequencies, _ = get_frequencies(point, data['hess'], point.geom)
    assert frequencies == pytest.approx([2321.7], rel=1.e-10)
    assert np.linalg.matrix_rank(data['hess']) == 1
    symmetry.calculate_symmetry(point)
    assert point.sigma_ext == 2
    assert not (tmp_path/'ion_neutral_h2_vib').exists()


@pytest.mark.parametrize('charge,mult,symbols,task,model', [
    (0, 1, 'H2', 'omol', 'uma-m-1p1'),
    (1, 2, 'H3', 'omol', 'uma-m-1p1'),
    (1, 2, 'H2', 'omat', 'uma-m-1p1'),
    (1, 2, 'H2', 'omol', 'another-model')])
def test_other_species_and_calculators_do_not_receive_the_reference(
        charge, mult, symbols, task, model):
    mol = Atoms(symbols)
    mol.info.update(charge=charge, spin=mult)
    calc = NeutralHydrogen()
    calc.task_name = task
    mol.calc = calc
    assert h2_cation_reference(mol, 'unused', model, .005, 10) is None
    assert mol.calc is calc


def test_local_uma_checkpoint_is_identified_by_its_backbone():
    backbone = type('Backbone', (), {'__module__': 'fairchem.core.models.uma.escn_md'})()
    calc = SimpleNamespace(predictor=SimpleNamespace(
        model=SimpleNamespace(module=SimpleNamespace(backbone=backbone))))
    assert uma_model('/models/checkpoint.pt', calc)
    assert not uma_model('/models/checkpoint.pt')


def test_failed_neutral_reference_is_not_reported_as_corrected_energy(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    mol = Atoms('H2', positions=[[0., 0., 0.], [.74, 0., 0.]])
    mol.info.update(charge=1, spin=2)
    mol.calc = NeutralHydrogen()
    monkeypatch.setattr('ase.optimize.BFGS.run', lambda *args, **kwargs: False)
    with pytest.raises(RuntimeError, match='did not converge'):
        h2_cation_reference(mol, 'ion', 'uma-m-1p1', .005, 0)


@pytest.mark.parametrize('worker', ['product', 'conformer', 'high'])
def test_generated_uma_workers_publish_the_complete_corrected_record(tmp_path, monkeypatch, worker):
    from kinbot.parameters import Parameters
    from kinbot.qc import QuantumChemistry
    from kinbot.species_routing import routing_name

    monkeypatch.chdir(tmp_path)
    (tmp_path/'conf').mkdir()
    (tmp_path/'input.json').write_text(json.dumps(dict(
        qc='fc', fc_model_path='uma-m-1p1', queuing='local', high_level=0,
        barrier_threshold=100.)))
    par = Parameters('input.json', show_warnings=False).par
    qc = QuantumChemistry(par)
    monkeypatch.setattr(qc, 'submit_qc', lambda *args, **kwargs: None)
    point = hydrogen(.74)
    name = routing_name(point)
    if worker == 'conformer':
        qc.qc_conf(point, point.geom, 0)
        label = 'conf/' + name + '_0000'
    else:
        qc.qc_opt(point, point.geom, initial_product=worker == 'product',
                  high_level=worker == 'high')
        label = name + ('_product' if worker == 'product' else '_well_high')
    core = ModuleType('fairchem.core')
    core.FAIRChemCalculator = lambda *args, **kwargs: NeutralHydrogen()
    monkeypatch.setitem(sys.modules, 'fairchem.core', core)
    monkeypatch.setattr('kinbot.fairchem_utils.load_predictor', lambda *args, **kwargs: None)
    with pytest.raises(SystemExit) as stopped:
        runpy.run_path(label+'.py')
    assert stopped.value.code == 0
    with open(label+'.pkl', 'rb') as handle:
        data = pickle.load(handle)['data']
    assert data['charge'] == 1 and data['multiplicity'] == 2
    assert data['h2_cation_reference']['ionization_energy_ev'] == 15.42593
    assert data['hess'].shape == (6, 6)
    assert qc.check_qc(label) == 'normal'
    assert not (tmp_path/(label+'.pkl')).exists()
    assert qc.get_qc_freq(label, 2)[1] == pytest.approx([2321.7])


@pytest.mark.parametrize('source', ['row', 'input_reference', 'worker', 'unknown'])
def test_old_uma_h2_cation_results_cannot_be_silently_reused(tmp_path, monkeypatch, source):
    from ase.db import connect
    from kinbot.fairchem_utils import check_h2_cation_cache

    monkeypatch.chdir(tmp_path)
    db = connect('kinbot.db')
    mol = Atoms('H2', positions=[[0., 0., 0.], [.544, 0., 0.]])
    data = {'status': 'normal', 'energy': -30.}
    if source == 'row':
        data.update(charge=1, multiplicity=2)
    elif source == 'input_reference':
        db.write(mol, name='stereochemistry/ion',
                 data={'input_charge': 1, 'input_multiplicity': 2})
    elif source == 'worker':
        (tmp_path/'ion.py').write_text("kwargs = {'charge': 1, 'mult': 2}\n")
    db.write(mol, name='ion', data=data)
    qc = SimpleNamespace(db=db, par={'fc_task_name': 'omol', 'fc_model_path': 'uma-m-1p1'})
    with pytest.raises(ValueError, match='fresh directory'):
        check_h2_cation_cache(qc, 'ion')
    assert len(list(db.select(name='ion'))) == 1


@pytest.mark.parametrize('uma_marker', [None, True])
def test_removing_local_checkpoint_does_not_allow_uncorrected_h2_cation(tmp_path, uma_marker):
    from ase.db import connect
    from kinbot.fairchem_utils import check_h2_cation_cache

    checkpoint = tmp_path/'local_weights.pt'
    checkpoint.write_text('Only file availability is tested; no model is loaded.\n')
    db = connect(str(tmp_path/'kinbot.db'))
    data = {'status': 'normal', 'charge': 1, 'multiplicity': 2, 'energy': -30.}
    if uma_marker is not None:
        data['uma_model'] = uma_marker
    db.write(Atoms('H2', positions=[[0., 0., 0.], [.544, 0., 0.]]), name='ion', data=data)
    qc = SimpleNamespace(db=db, par={'fc_task_name': 'omol', 'fc_model_path': str(checkpoint)})
    with pytest.raises(ValueError, match='fresh directory'):
        check_h2_cation_cache(qc, 'ion')
    checkpoint.unlink()
    with pytest.raises(ValueError, match='fresh directory'):
        check_h2_cation_cache(qc, 'ion')
    assert len(list(db.select(name='ion'))) == 1


@pytest.mark.parametrize('model', ['uma-m-1p1', 'missing_weights.pt', 'another-model'])
@pytest.mark.parametrize('case', ['neutral', 'non_uma', 'corrected', 'other_molecule'])
def test_corrected_or_unaffected_results_remain_reusable(tmp_path, case, model):
    from ase.db import connect
    from kinbot.fairchem_utils import check_h2_cation_cache

    db = connect(str(tmp_path/'kinbot.db'))
    mol = Atoms('HO' if case == 'other_molecule' else 'H2')
    data = {'status': 'normal', 'charge': 1, 'multiplicity': 2}
    if case == 'neutral':
        data.update(charge=0, multiplicity=1)
    elif case == 'non_uma':
        data['uma_model'] = False
    elif case == 'corrected':
        data['h2_cation_reference'] = {'ionization_energy_ev': 15.42593}
    db.write(mol, name='job', data=data)
    qc = SimpleNamespace(db=db, par={'fc_task_name': 'omol', 'fc_model_path': model})
    check_h2_cation_cache(qc, 'job')


@pytest.mark.parametrize('high_level', [0, 1])
@pytest.mark.parametrize('conformers', [0, 1])
@pytest.mark.parametrize('model', ['uma-m-1p1', 'missing_weights.pt'])
def test_pes_checks_old_h2_cation_reference_at_the_selected_profile(
        tmp_path, monkeypatch, high_level, conformers, model):
    from ase.db import connect
    from kinbot import pes
    from kinbot.parameters import Parameters
    from kinbot.species_routing import routing_name

    monkeypatch.chdir(tmp_path)
    (tmp_path/'input.json').write_text(json.dumps(dict(
        barrier_threshold=100., high_level=high_level,
        l1='gauss' if high_level else 'fc',
        l2='fc' if high_level else 'gauss',
        fc_model_path=model, conformer_search=conformers or high_level)))
    par = Parameters('input.json', show_warnings=False).par
    point = hydrogen(1.058)
    name = routing_name(point)
    job = (name+'_well_high' if high_level else
           'conf/'+name+'_low' if conformers else name+'_well')
    worker = tmp_path/'parent'
    worker.mkdir()
    db = connect(str(worker/'kinbot.db'))
    mol = Atoms('H2', positions=point.geom)
    data = {'status': 'normal', 'charge': 1, 'multiplicity': 2, 'energy': -15., 'zpe': .005}
    db.write(mol, name=job, data=data)
    monkeypatch.setattr(pes, 'ensure_current_run', lambda *args: None)
    monkeypatch.setattr(pes, 'logger', logging.getLogger('KinBot'), raising=False)
    with pytest.raises(ValueError, match='corrected UMA H2'):
        pes.get_energy(['parent'], name, 0, high_level, conf=conformers, par=par)
    data['h2_cation_reference'] = {'ionization_energy_ev': 15.42593}
    db.write(mol, name=job, data=data)
    energy, zpe = pes.get_energy(['parent'], name, 0, high_level, conf=conformers, par=par)
    assert energy == pytest.approx(-15.*constants.EVtoHARTREE)
    assert zpe == .005


@pytest.mark.parametrize('published', [False, True])
def test_pes_cache_check_reads_old_neutral_worker_from_its_directory(tmp_path, monkeypatch, published):
    from ase.db import connect
    from kinbot.fairchem_utils import check_h2_cation_cache

    worker = tmp_path/'parent'
    worker.mkdir()
    db = connect(str(worker/'kinbot.db'))
    data = {'status': 'normal'}
    source = 'conformer' if published else 'neutral'
    if published:
        data['copied_from_job'] = source
    db.write(Atoms('H2'), name='neutral', data=data)
    (worker/(source+'.py')).write_text("kwargs = {'charge': 0, 'mult': 1}\n")
    monkeypatch.chdir(tmp_path)
    qc = SimpleNamespace(db=db, par={'fc_task_name': 'omol', 'fc_model_path': 'uma-m-1p1'})
    check_h2_cation_cache(qc, 'neutral', str(worker))


@pytest.mark.parametrize('source', ['row', 'input_reference', 'calculator'])
def test_pes_keeps_the_charge_of_a_saved_long_bonded_h2_cation(tmp_path, monkeypatch, source):
    from ase.db import connect
    from ase.calculators.singlepoint import SinglePointCalculator
    from kinbot import pes
    from kinbot.species_routing import routing_name

    monkeypatch.chdir(tmp_path)
    worker = tmp_path/'parent'
    worker.mkdir()
    db = connect(str(worker/'kinbot.db'))
    point = hydrogen(1.1142815143)
    name, job = routing_name(point), routing_name(point)+'_well'
    mol = Atoms('H2', positions=point.geom)
    data = {'status': 'normal', 'energy': -15., 'zpe': .005}
    if source == 'row':
        data.update(charge=1, multiplicity=2)
    elif source == 'input_reference':
        db.write(mol, name='stereochemistry/'+job,
                 data={'input_charge': 1, 'input_multiplicity': 2})
    else:
        mol.calc = SinglePointCalculator(mol, energy=-15.)
        mol.calc.parameters.update(charge=1, mult=2)
    db.write(mol, name=job, data=data)
    monkeypatch.setattr(pes, 'ensure_current_run', lambda *args: None)
    monkeypatch.setattr(pes, 'logger', logging.getLogger('KinBot'), raising=False)
    energy, zpe = pes.get_energy(['parent'], name, 0, 0)
    assert energy == pytest.approx(-15.*constants.EVtoHARTREE)
    assert zpe == .005


@pytest.mark.parametrize('profiled', [False, True])
@pytest.mark.parametrize('charge', [0, 1])
@pytest.mark.parametrize('selection', ['sample', 'failed_samples', 'lower_parent'])
def test_actual_conformer_selection_preserves_h2_state_and_reference(
        tmp_path, monkeypatch, profiled, charge, selection):
    from pathlib import Path
    from kinbot.conformers import Conformers
    from kinbot.parameters import Parameters
    from kinbot.qc import QuantumChemistry
    from kinbot.species_routing import routing_name
    from kinbot.stereo_routing import guard_well_job

    monkeypatch.chdir(tmp_path)
    Path('conf').mkdir()
    settings = dict(qc='fc', fc_model_path='uma-s-1p1', high_level=0,
                    queuing='local', barrier_threshold=100.)
    if profiled:
        settings['l1_profile'] = dict(calculator='fairchem', model_path='uma-s-1p1',
                                     device='cpu', task_name='omol')
    Path('input.json').write_text(json.dumps(settings))
    par = Parameters('input.json', show_warnings=False).par
    qc = QuantumChemistry(par)
    monkeypatch.setattr(qc, 'submit_qc', lambda *args, **kwargs: pytest.fail('submitted calculation'))
    point = hydrogen(1.058 if charge else .74, charge=charge, mult=charge+1)
    point.freq = [2321.7 if charge else 4400.]
    name = routing_name(point)
    parent = name+'_well'
    guard_well_job(qc, point, point.geom, parent)
    data = dict(status='normal', charge=charge, multiplicity=charge+1, uma_model=True,
                energy=-15., zpe=.005, frequencies=point.freq, hess=np.eye(6))
    if charge:
        data['h2_cation_reference'] = dict(ionization_energy_ev=15.42593,
                                          neutral_h2_energy_ev=-30., neutral_h2_zpe_hartree=.01)
    qc.db.write(Atoms('H2', positions=point.geom), name=parent, data=data)
    Path(parent+'_sella.log').write_text('done\n')
    search = Conformers(point, par, qc)
    sample = search.get_job_name(0)
    sample_data = dict(data, energy=-14. if selection=='lower_parent' else -15.)
    qc.db.write(Atoms('H2', positions=point.geom), name=sample, data=sample_data)
    Path(sample+'_sella.log').write_text('done\n')
    search.conf = 1
    search.conf_status = [1 if selection=='failed_samples' else 0]
    result = search.check_conformers()
    assert result[0] == 1
    expected = parent if selection!='sample' else sample
    assert search.selected_job == expected
    selected = next(qc.db.select(name=expected, sort='-id', limit=1))
    target = 'conf/'+name+'_low'
    saved = next(qc.db.select(name=target, sort='-id', limit=1))
    for field, value in selected.data.items():
        if field=='hess':
            np.testing.assert_array_equal(saved.data[field], value)
        else:
            assert saved.data[field] == value
    assert saved.data['copied_from_job'] == expected
    assert saved.data['copied_from_row_id'] == selected.id
    assert saved.data['nconf'] == 1
    assert qc.check_qc(target) == 'normal'
    assert qc.get_qc_energy(target)[1] == pytest.approx(data['energy']*constants.EVtoHARTREE)
    assert qc.get_qc_freq(target, 2)[1] == point.freq
    assert search.completed_conformers() == 1


@pytest.mark.parametrize('fails', [False, True])
def test_atomic_fairchem_failure_is_a_completed_error_record(tmp_path, monkeypatch, fails):
    from pathlib import Path
    from kinbot.parameters import Parameters
    from kinbot.qc import QuantumChemistry
    from kinbot.species_routing import routing_name

    class AtomicReference(Calculator):
        implemented_properties = ['energy', 'forces']

        def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            if fails:
                raise ValueError('This model has not stored this element with this charge.')
            self.results = dict(energy=-10., forces=np.zeros((1,3)))

    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps(dict(qc='fc', fc_model_path='uma-s-1p1',
        high_level=0, queuing='local', barrier_threshold=100.)))
    qc = QuantumChemistry(Parameters('input.json', show_warnings=False).par)
    monkeypatch.setattr(qc, 'submit_qc', lambda *args, **kwargs: None)
    carbon = StationaryPoint('carbon cation', 1, 2, atom=['C'], geom=np.zeros((1,3)))
    carbon.characterize()
    qc.qc_opt(carbon, carbon.geom, initial_product=True)
    job = routing_name(carbon)+'_product'
    core = ModuleType('fairchem.core')
    core.FAIRChemCalculator = lambda *args, **kwargs: AtomicReference()
    monkeypatch.setitem(sys.modules, 'fairchem.core', core)
    monkeypatch.setattr('kinbot.fairchem_utils.load_predictor', lambda *args, **kwargs: None)
    with pytest.raises(SystemExit) as stopped:
        runpy.run_path(job+'.py')
    assert stopped.value.code == 0
    assert qc.check_qc(job) == ('error' if fails else 'normal')
    row = next(qc.db.select(name=job, sort='-id', limit=1))
    if fails:
        assert 'not stored this element' in row.data['reason']
        assert 'energy' not in row.data
    else:
        assert row.data['energy'] == -10.
        assert row.data['zpe'] == 0.
        assert len(row.data['frequencies']) == 0
    assert Path(job+'_sella.log').read_text().endswith('done\n')
