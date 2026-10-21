"""Fragment spins stay coupled to the IRC, including after product relaxation."""
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule
from ase.calculators.calculator import Calculator, all_changes
from ase.db import connect

from kinbot import constants
from kinbot.bimol_products import (ProductSpinSelection, check_product_spins,
                                   multiplicities, spin_allowed, spin_copy)
from kinbot.calculation import load_calculation_record
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.reaction_generator import ReactionGenerator
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint


def point(formula, mult=1, charge=0):
    atoms = molecule(formula)
    result = StationaryPoint(formula, charge, mult,
        atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    result.characterize()
    result.calc_chemid()
    return result


def endpoint(*fragments, mult, charge=0):
    atoms = np.concatenate([p.atom for p in fragments])
    geom = np.concatenate([p.geom + [8. * i, 0., 0.]
                           for i, p in enumerate(fragments)])
    result = StationaryPoint('endpoint', charge, mult, atom=atoms, geom=geom)
    result.characterize()
    return result


@pytest.mark.parametrize('spins,total,allowed', [
    ([2, 2], 1, True), ([2, 2], 3, True), ([2, 2], 2, False),
    ([1, 1], 3, False), ([1, 3], 3, True), ([2, 3], 2, True),
    ([2, 3], 4, True), ([2, 2, 2], 2, True), ([2, 2, 2], 4, True),
    ([2, 2, 2], 6, False), ([3, 3, 3], 1, True)])
def test_angular_momentum_coupling(spins, total, allowed):
    assert spin_allowed(spins, total) == allowed


@pytest.mark.parametrize('formula,mult,charge,expected', [
    ('C2H4', 1, 0, [1]), ('C2H6', 1, 0, [1]), ('CH3', 2, 0, [2]),
    ('O2', 3, 0, [1, 3]), ('O', 3, 0, [1, 3]), ('O', 2, 1, [2, 4]),
    ('N', 2, 0, [2, 4]), ('NH', 1, 0, [1, 3]), ('NH3', 1, 0, [1]),
    ('S', 3, 0, [1, 3])])
def test_candidate_electron_parity_and_valence(formula, mult, charge, expected):
    assert multiplicities(point(formula, mult, charge)) == expected


@pytest.mark.parametrize('total,needed,expected', [
    (1, False, [1, 1]), (3, True, [1, 3])])
def test_ethyl_disproportionation_preserves_total_spin(total, needed, expected):
    products = endpoint(point('C2H6'), point('C2H4'), mult=total)
    original = products.geom.copy(), products.bond.copy(), products.mult
    search = ProductSpinSelection(products)
    assert search.needed == needed
    assert [[search.states[i].mult for i in comb]
            for comb in search.combinations] == [expected]
    assert search.evidence['multiple_bond_excitation_added'] == needed
    np.testing.assert_array_equal(products.geom, original[0])
    np.testing.assert_array_equal(products.bond, original[1])
    assert products.mult == original[2]


def test_ordinary_doublet_and_connected_products_need_no_comparison():
    assert not ProductSpinSelection(endpoint(point('CH3', 2), point('H2'), mult=2)).needed
    assert not ProductSpinSelection(point('C2H4', 3)).needed


def test_charged_variants_keep_one_state_per_physical_fragment():
    search = ProductSpinSelection(endpoint(point('O', 3), point('H', 2), mult=1, charge=1))
    assert len(search.maps) == 2
    for comb in search.combinations:
        chosen = [search.states[i] for i in comb]
        assert sum(p.charge for p in chosen) == 1
        assert sorted(a for p in chosen for a in p.atom) == ['H', 'O']
        assert spin_allowed([p.mult for p in chosen], 1)


def test_unique_allowed_charged_combination_cannot_use_ordinary_charge_selection():
    search = ProductSpinSelection(endpoint(point('H', 2), point('CH4'), mult=3, charge=1))
    assert search.needed
    assert len(search.combinations) == 1
    selected = [search.states[i] for i in search.combinations[0]]
    assert [(p.charge, p.mult) for p in selected] == [(0, 2), (1, 2)]


def test_charge_variants_with_same_ordinary_job_name_refuse_reuse():
    with pytest.raises(ValueError, match='same ordinary job name'):
        ProductSpinSelection(endpoint(point('CO'), point('H2'), mult=1, charge=2))


def test_final_validation_rejects_incompatible_or_changed_products():
    reference = endpoint(point('C2H6'), point('C2H4'), mult=3)
    with pytest.raises(ValueError, match='multiplicity 3'):
        check_product_spins([point('C2H6'), point('C2H4')], reference)
    check_product_spins([point('C2H6'), point('C2H4', 3)], reference)
    broken = point('H2')
    broken.geom[1] += [8., 0., 0.]
    with pytest.raises(ValueError, match='connectivity'):
        ProductSpinSelection._check_product(point('H2'), broken)


class CompletedQC:
    """Completed energy records; never launch an electronic-structure program."""
    qc = 'gauss'

    def __init__(self, path, high_level=0):
        self.db = connect(path / 'kinbot.db')
        self.par = dict(bimol=1, high_level=high_level, calc_kwargs={},
                        multi_conf_tst=0, rotor_scan=0, conformer_search=0)
        self.sp_calls, self.opt_calls = [], []
        self.failed = set()
        self.pending = set()

    def for_level(self, level):
        return self

    def write(self, species, job, energy, zpe=0.):
        self.db.write(Atoms(species.atom, positions=species.geom), name=job,
            data=dict(status='normal', energy=energy / constants.EVtoHARTREE,
                      zpe=zpe, frequencies=[]))

    def qc_product_spin(self, species, geom):
        from kinbot.calculation import array_fingerprint
        job = f'product_spin/{routing_name(species)}_{array_fingerprint(geom)}_l1'
        self.sp_calls.append(job)
        if not list(self.db.select(name=job)):
            self.write(species, job, -20. if species.mult == 1 else -10.)
        return job

    def qc_opt(self, species, geom):
        job = routing_name(species) + '_well'
        if not list(self.db.select(name=job)):
            self.opt_calls.append(job)
            self.write(species, job, -10. - species.mult)

    def check_qc(self, job):
        return 'error' if job in self.failed else 'running' if job in self.pending else 'normal'

    def get_qc_energy(self, job):
        return 0, list(self.db.select(name=job))[-1].data.energy * constants.EVtoHARTREE


class CompletedOptimization:
    def __init__(self, species, par, qc):
        self.species, self.par, self.qc = species, par, qc
        self.shigh = self.shir = 1
        self.just_high = self.defer_hir = False
        self.selected_job = routing_name(species) + ('_high' if par['high_level'] else '_well')
        # L2 deliberately reverses L1: the decision must use the requested final level.
        energy = (-10. + species.mult) if par['high_level'] else (-10. - species.mult)
        qc.write(species, self.selected_job, energy)
        load_calculation_record(species, qc, self.selected_job)
        self.calls = 0

    def do_optimization(self):
        self.calls += 1


@pytest.mark.parametrize('high_level,expected', [(0, 3), (1, 1)])
@pytest.mark.parametrize('backend', ['gauss', 'fc'])
def test_relaxed_energies_select_states_without_single_point_cutoff(
        tmp_path, monkeypatch, high_level, expected, backend):
    monkeypatch.setattr('kinbot.bimol_products.Optimize', CompletedOptimization)
    qc = CompletedQC(tmp_path, high_level)
    qc.qc = backend
    search = ProductSpinSelection(endpoint(point('OH', 2), point('O', 3), mult=2))
    assert search.advance(qc, qc.par)
    products = [opt.species for opt in search.selected_optimizations()]
    assert [p.mult for p in products] == [2, expected]
    assert len(qc.sp_calls) == len(qc.opt_calls) == 3
    assert search.evidence['level'] == ('l2' if high_level else 'l1')
    assert len(search.evidence['allowed_combinations']) == 2
    calls = list(qc.sp_calls), list(qc.opt_calls)
    assert search.advance(qc, qc.par)
    assert (qc.sp_calls, qc.opt_calls) == calls
    json.dumps(search.evidence, allow_nan=False)


def test_duplicate_fragments_keep_both_slots_and_reuse_optimizations(tmp_path, monkeypatch):
    monkeypatch.setattr('kinbot.bimol_products.Optimize', CompletedOptimization)
    qc = CompletedQC(tmp_path)
    search = ProductSpinSelection(endpoint(point('O', 3), point('O', 3), mult=1))
    assert search.advance(qc, qc.par)
    selected = search.selected_optimizations()
    assert len(selected) == 2 and selected[0] is selected[1]
    assert len(qc.opt_calls) == 2
    assert len(search.combinations) == 2


@pytest.mark.parametrize('failure', ['single_point', 'optimization', 'refinement'])
def test_failed_candidate_does_not_select_from_incomplete_comparison(tmp_path, monkeypatch, failure):
    monkeypatch.setattr('kinbot.bimol_products.Optimize', CompletedOptimization)
    qc = CompletedQC(tmp_path)
    search = ProductSpinSelection(endpoint(point('OH', 2), point('O', 3), mult=2))
    target = next(i for i in search.active if search.states[i].mult == 3)
    if failure == 'single_point':
        job = qc.qc_product_spin(search.states[target], search.states[target].geom)
        qc.failed.add(job)
    elif failure == 'optimization':
        qc.failed.add(routing_name(search.states[target]) + '_well')
    else:
        def failed(*args):
            opt = CompletedOptimization(*args)
            if opt.species.mult == 3:
                opt.shigh = -999
            return opt
        monkeypatch.setattr('kinbot.bimol_products.Optimize', failed)
    with pytest.raises(ValueError, match='unresolved'):
        search.advance(qc, qc.par)
    assert search.selected is None


def test_pending_jobs_resume_and_reuse_completed_product_models(tmp_path, monkeypatch):
    monkeypatch.setattr('kinbot.bimol_products.Optimize', CompletedOptimization)
    qc = CompletedQC(tmp_path)
    products = endpoint(point('OH', 2), point('O', 3), mult=2)
    search = ProductSpinSelection(products)
    job = qc.qc_product_spin(search.states[0], search.states[0].geom)
    qc.pending.add(job)
    assert not search.advance(qc, qc.par)
    assert not qc.opt_calls
    qc.pending.clear()
    assert search.advance(qc, qc.par)
    completed_calls = [opt.calls for opt in search.optimizations.values()]
    repeat = ProductSpinSelection(products)
    assert repeat.advance(qc, qc.par, search.optimizations.values())
    assert repeat.selected_optimizations() == search.selected_optimizations()
    assert len(qc.opt_calls) == 3
    assert [opt.calls for opt in search.optimizations.values()] == completed_calls


@pytest.mark.parametrize('backend', ['gauss', 'fc'])
@pytest.mark.parametrize('key', ['mult', 'multiplicity', 'charge'])
def test_overridden_product_state_fails_before_jobs(tmp_path, backend, key):
    qc = CompletedQC(tmp_path)
    qc.qc = backend
    qc.par['calc_kwargs'] = {key: 1}
    search = ProductSpinSelection(endpoint(point('OH', 2), point('O', 3), mult=2))
    with pytest.raises(ValueError):
        search.advance(qc, qc.par)
    assert not qc.sp_calls


def test_reaction_generator_keeps_selected_optimization_and_failure_status(tmp_path, monkeypatch):
    monkeypatch.setattr('kinbot.bimol_products.Optimize', CompletedOptimization)
    qc = CompletedQC(tmp_path)
    reaction = SimpleNamespace(irc_prod=endpoint(point('C2H6'), point('C2H4'), mult=3),
                               instance_name='ethyl_pair', final_irc_processed=True)
    generator = ReactionGenerator.__new__(ReactionGenerator)
    generator.qc, generator.par = qc, qc.par
    generator.species = SimpleNamespace(reac_obj=[reaction], reac_ts_done=[2])
    generator._record_bimol_status = Mock()
    assert generator._advance_product_spins(0, reaction)
    assert generator.species.reac_ts_done == [4]
    assert reaction.products == [opt.species for opt in reaction.prod_opt]
    assert reaction.prod_opt == reaction.product_spin_search.selected_optimizations()
    assert reaction.product_spin_selection['selected_combination']
    # A separate unresolved channel retains its completed IRC and is omitted.
    reaction2 = SimpleNamespace(irc_prod=reaction.irc_prod, instance_name='other',
                                final_irc_processed=True)
    qc.failed.add(qc.sp_calls[0])
    assert generator._advance_product_spins(0, reaction2)
    assert generator.species.reac_ts_done == [-999]
    assert reaction2.final_irc_processed
    assert generator._record_bimol_status.call_args.kwargs['final_status'] == 'product_spin_unresolved'


def real_qc(tmp_path, monkeypatch, backend='gaussian', fc_compile=False):
    monkeypatch.chdir(tmp_path)
    l0 = dict(calculator='gaussian', method='AM1')
    l1 = dict(calculator=backend, method='CAM-B3LYP', basis='6-31G*')
    if backend == 'fairchem':
        l0 = dict(calculator=backend, model_path='uma-s-1p1', device='cpu')
        l1 = dict(calculator=backend, model_path='uma-m-1p1', device='cpu', task_name='omol')
    Path('input.json').write_text(json.dumps(dict(
        barrier_threshold=100., bimol=1, structure=[[], []],
        queuing='local', error_missing_local=1, fc_compile=fc_compile,
        l0_profile=l0, l1_profile=l1)))
    monkeypatch.setattr(QuantumChemistry, 'submit_qc', Mock(return_value=1))
    return QuantumChemistry(Parameters('input.json', show_warnings=False).par)


@pytest.mark.parametrize('backend', ['gaussian', 'qchem', 'orca', 'fairchem'])
@pytest.mark.parametrize('failure', [False, True])
def test_generated_single_point_worker_records_success_and_startup_failure(
        tmp_path, monkeypatch, backend, failure):
    qc = real_qc(tmp_path, monkeypatch, backend)
    species = point('C2H4', 3)
    job = qc.qc_product_spin(species, species.geom)
    assert qc._level_for(job) == 'l1'
    captured = {}

    class Energy(Calculator):
        implemented_properties = ['energy']
        def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            self.results = {'energy': -27.}

    def calculator(code, module, name, kwargs, *args, compile_model=False):
        captured.update(kwargs)
        assert kwargs['mult' if backend != 'qchem' else 'multiplicity'] == 3
        assert kwargs['charge'] == 0
        assert not {'opt', 'freq', 'geom'} & {k.lower() for k in kwargs}
        if backend == 'qchem':
            assert kwargs['jobtype'] == 'sp'
        elif backend == 'orca':
            assert 'EnGrad' not in kwargs['orcasimpleinput']
            assert kwargs['directory'] == job
        elif backend == 'gaussian':
            assert Path(kwargs['label']).parent / kwargs['chk'] == Path(job)
        if failure:
            raise RuntimeError('test startup failure')
        return Energy()

    monkeypatch.setattr('kinbot.bimol.create_calculator', calculator)
    source = Path(job + '.py').read_text()
    exec(compile(source, job + '.py', 'exec'), {})
    row = list(qc.db.select(name=job))[-1]
    assert row.data.status == ('error' if failure else 'normal')
    assert qc.check_qc(job) == ('error' if failure else 'normal')
    assert row.data.multiplicity == 3
    np.testing.assert_array_equal(row.positions, species.geom)
    assert Path(job + '.done').is_file()
    assert captured
    # Saved jobs are not rewritten and changed coordinates/state cannot reuse them.
    Path(job + '.py').write_text('saved worker evidence')
    assert qc.qc_product_spin(species, species.geom) == job
    assert Path(job + '.py').read_text() == 'saved worker evidence'
    moved = species.geom.copy()
    moved[0, 0] += .01
    assert qc.qc_product_spin(species, moved) != job
    singlet = spin_copy(species, 1)
    assert qc.qc_product_spin(singlet, singlet.geom) != job


@pytest.mark.parametrize('formula,mult,charge', [('C2H4', 1, 0), ('C2H4', 3, 0), ('O', 2, 1)])
@pytest.mark.parametrize('compile_model', [False, True])
def test_fairchem_product_worker_uses_l1_model_and_fragment_state(
        tmp_path, monkeypatch, formula, mult, charge, compile_model):
    qc = real_qc(tmp_path, monkeypatch, 'fairchem', fc_compile=compile_model)
    species = point(formula, mult, charge)
    job = qc.qc_product_spin(species, species.geom)
    predictor = object()
    load = Mock(return_value=predictor)
    monkeypatch.setattr('kinbot.fairchem_utils.load_predictor', load)

    class Energy(Calculator):
        implemented_properties = ['energy']

        def __init__(self, model, task_name):
            super().__init__()
            assert model is predictor
            assert task_name == 'omol'

        def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            assert atoms.info == dict(charge=charge, spin=mult)
            np.testing.assert_array_equal(atoms.positions, species.geom)
            self.results = {'energy': -27.}

    monkeypatch.setitem(sys.modules, 'fairchem.core', SimpleNamespace(FAIRChemCalculator=Energy))
    exec(compile(Path(job + '.py').read_text(), job + '.py', 'exec'), {})
    load.assert_called_once_with('uma-m-1p1', 'cpu', compile_model=compile_model)
    row = qc.db.get(name=job)
    assert row.data.backend == 'fc'
    assert row.data.method == 'uma-m-1p1'
    assert row.data.basis is None
    assert row.data.charge == charge and row.data.multiplicity == mult
    assert qc.check_qc(job) == 'normal'
    assert qc.get_qc_energy(job)[1] == pytest.approx(-27. * constants.EVtoHARTREE)


def test_single_point_is_bimolecular_only(tmp_path, monkeypatch):
    qc = real_qc(tmp_path, monkeypatch)
    qc.l1.par['bimol'] = 0
    species = point('O', 3)
    with pytest.raises(ValueError, match='bimolecular'):
        qc.qc_product_spin(species, species.geom)


def test_single_point_orca_submission_preserves_full_path_without_scratch_cleanup(tmp_path, monkeypatch):
    submit = QuantumChemistry.submit_qc
    qc = real_qc(tmp_path, monkeypatch, 'orca')
    species = point('O', 3)
    job = qc.qc_product_spin(species, species.geom)
    backend = qc.l1
    backend.queuing = 'slurm'
    backend.queue_job_limit = 0
    backend.check_qc = Mock(return_value=0)
    process = Mock()
    process.communicate.return_value = (b'Submitted batch job 123\n', b'')
    monkeypatch.setattr('kinbot.qc.subprocess.Popen', Mock(return_value=process))
    assert submit(backend, job, 1) == 1
    script = Path(job + constants.qext['slurm']).read_text()
    assert f'python {job}.py' in script
    assert 'SCRATCH_DIR' not in script
    # Existing ORCA optimization workers keep their former submission behavior.
    assert submit(backend, 'ordinary_well', 1) == 1
    assert 'SCRATCH_DIR' in Path('ordinary_well' + constants.qext['slurm']).read_text()


def test_real_product_optimizer_completes_and_retains_selected_properties(tmp_path, monkeypatch):
    qc = real_qc(tmp_path, monkeypatch)
    search = ProductSpinSelection(endpoint(point('OH', 2), point('O', 3), mult=2))
    sp_writer = qc.l1.qc_product_spin

    def write(species, job, energy, zpe):
        qc.db.write(Atoms(species.atom, positions=species.geom), name=job,
            data=dict(status='normal', energy=energy / constants.EVtoHARTREE,
                      zpe=zpe, frequencies=[3600.] if species.natom > 1 else []))

    def single_point(species, geom):
        job = sp_writer(species, geom)
        write(species, job, -10. - species.mult, 0.)
        Path(job + '.done').write_text('done\n')
        return job

    def optimize(species, geom):
        job = routing_name(species) + '_well'
        # Triplet O is lower in electronic energy, but its artificial test ZPE
        # reverses the sum; the comparison must read both saved quantities.
        write(species, job, -10. - species.mult,
              3. if species.mult == 3 else 0.)
        Path(job + '.log').write_text('done\n')

    monkeypatch.setattr(qc.l1, 'qc_product_spin', single_point)
    monkeypatch.setattr(qc.l1, 'qc_opt', optimize)
    assert search.advance(qc, qc.par)
    selected = search.selected_optimizations()
    assert [opt.species.mult for opt in selected] == [2, 1]
    assert all(opt.shigh == opt.shir == 1 for opt in selected)
    assert all(opt.selected_job == opt.species.source_job for opt in selected)
    assert selected[0].species.freq == [3600.]
