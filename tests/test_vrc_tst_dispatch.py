"""VRC correction jobs remain distinct and complete before rotdPy handoff."""

import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import Mock, patch

from ase import Atoms
from ase.db import connect
import numpy as np
import pytest

from kinbot.vrc_tst_scan import VTS, _drive_vrc_dispatch
from kinbot.fragments import Fragment


@pytest.mark.parametrize('pes,execute', [(0, 1), (1, 0)])
def test_pes_workers_prepare_rotd_inputs_without_duplicate_sampling(tmp_path, monkeypatch, pes, execute):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr('kinbot.rotdpy.ensure_available', lambda: None)
    seen = []
    monkeypatch.setattr('kinbot.pes.create_rotdpy_inputs',
        lambda parameters, *args, **kwargs: seen.append(parameters['rotdpy_run']))
    well = SimpleNamespace(wellorts=1, name='parent', chemid='parent', reac_obj=[])
    par = {'pes': pes, 'pp_orient': 'geometric', 'rotdpy_run': 1,
           'vrc_tst_scan': {}, 'vrc_tst_noscan': {}}
    vts = VTS(well, par, None)
    vts.scan_reac['r'] = SimpleNamespace(products=[
        SimpleNamespace(wellorts=1, name='a'), SimpleNamespace(wellorts=1, name='b')])
    vts.calculate_correction_potentials()
    assert seen == [execute]
    assert par['rotdpy_run'] == 1


def test_homo_oriented_vrc_fails_before_submission_without_gaussian_tools(
        monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr('kinbot.anl.site._executable', lambda *args: None)
    vts = VTS(SimpleNamespace(reac_obj=[]), {
        'pp_orient': 'homo', 'rotdpy_run': 0,
        'vrc_tst_scan': {}, 'vrc_tst_noscan': {},
    }, SimpleNamespace())

    with pytest.raises(FileNotFoundError, match='formchk'):
        vts.calculate_correction_potentials()


def test_vrc_fragment_handoff_rejects_incomplete_geometry():
    product = SimpleNamespace(
        wellorts=1, name='fragment', natom=2, atom=['H', 'H'], mult=1)
    qc = SimpleNamespace(get_qc_geom=Mock(return_value=(
        -1, np.full((2, 3), np.nan), ['H', 'H'])))
    vts = VTS(SimpleNamespace(), {'pp_orient': 'geometric'}, qc)
    vts.scan_reac['r'] = SimpleNamespace(products=[product])

    with pytest.raises(RuntimeError, match='finite, complete geometry'):
        vts.save_products(['r'])


def test_homo_orientation_honors_configured_mode():
    fragment = Fragment.__new__(Fragment)
    fragment.orient = {}
    fragment.par = {'pp_orient': 'homo'}
    expected = np.asarray([[1., 0., 0.]], dtype=np.float32)
    fragment.pp_from_homo = Mock(return_value=expected)

    assert np.array_equal(fragment.get_pp_orientation(2), expected)
    fragment.pp_from_homo.assert_called_once_with(index=2)


class _FakeMolpro:
    def __init__(self, species, parameters):
        self.species = species

    def create_molpro_input(self, name='', VTS=False, sample=False,
                            vrc_level=None):
        assert VTS
        Path(f'vrctst/molpro/{name}.inp').write_text(
            f'***,{name}\n! sample={sample}\n! level={vrc_level}\n')

    def get_molpro_energy(self, key, name='', VTS=False):
        assert key == 'MYENERGY'
        assert VTS
        if not Path(f'vrctst/molpro/{name}.out').is_file():
            return 0, -1.
        return 1, (-10. if name.endswith('_fr') else -11.)


@pytest.mark.parametrize('separation', [30., 55.])
def test_noscan_uses_distinct_sampling_and_high_level_asymptotes(separation):
    with TemporaryDirectory() as temporary:
        previous = Path.cwd()
        os.chdir(temporary)
        try:
            Path('vrctst/molpro').mkdir(parents=True)
            atoms = Atoms('H2', positions=[[0., 0., 0.], [0., 0., separation]])
            connect('kinbot.db').write(
                atoms, name='vrctst/r_vts_pt_asymptote_fr',
                data={'status': 'normal'})
            well = SimpleNamespace(
                chemid='parent', geom=np.asarray([[0., 0., 0.],
                                                   [0., 0., .74]]))
            parameters = {
                'vrc_tst_scan_points': [2.5],
                'vrc_tst_asymptote_distance': separation,
                'vrc_tst_scan_molpro_key': 'MYENERGY',
                'vrc_tst_sample_method': 'caspt2(2,2)',
                'vrc_tst_sample_basis': 'vdz',
                'vrc_tst_high_method': 'mrci+q(2,2)',
                'vrc_tst_high_basis': 'vtz',
            }
            vts = VTS(well, parameters, None)
            product_a = SimpleNamespace(
                atom=['H'], geom=[[0., 0., 0.]], mult=2)
            product_b = SimpleNamespace(
                atom=['H'], geom=[[0., 0., separation]], mult=2)
            vts.scan_reac['r'] = SimpleNamespace(
                usym=[[[0]], [[0]]], equiv=[[0], [1]],
                maps=[np.asarray([0]), np.asarray([1])],
                products=[product_a, product_b],
                species=SimpleNamespace(mult=1))
            observed = set()

            def dispatch(pending):
                observed.update(pending)
                for name in pending:
                    Path(f'vrctst/molpro/{name}.out').write_text(
                        'Molpro calculation terminated\n')
                    source = Path(f'vrctst/molpro/{name}.inp')
                    Path(f'vrctst/molpro/{name}.input.sha256').write_text(
                        __import__('hashlib').sha256(
                            source.read_bytes()).hexdigest() + '\n')

            with patch('kinbot.vrc_tst_scan.Molpro', _FakeMolpro), \
                    patch.object(vts, '_dispatch_molpro_corrections', dispatch):
                vts.energies(['r'], noscan=True)

            assert observed == {
                'r_vts_pt_asymptote_fr', 'r_vts_pt_asymptote'}
            correction = json.loads(
                Path('vrctst/corr_r.json').read_text())
            assert correction['dist'] == [separation]
            assert correction['e_samp'] == [0.]
            assert correction['e_high'] == [0.]
            assert correction['e_inf_samp'] == -10.
            assert correction['e_inf_high'] == -11.
            assert correction['levels']['trusted_correction'] == {
                'method': 'mrci+q(2,2)', 'basis': 'vtz'}
        finally:
            os.chdir(previous)


def test_default_multireference_scan_publishes_only_caspt2_geometries():
    with TemporaryDirectory() as temporary:
        previous = Path.cwd()
        os.chdir(temporary)
        try:
            parent_geom = np.asarray([
                [0., 0., 0.], [0., 0., .74],
                [0., 0., 3.24], [0., 0., 3.98]])
            product_a = SimpleNamespace(
                atom=['H', 'H'], geom=parent_geom[:2].copy(), mult=1,
                bond01=np.asarray([[0, 1], [1, 0]]))
            product_b = SimpleNamespace(
                atom=['H', 'H'], geom=parent_geom[2:].copy(), mult=1,
                bond01=np.asarray([[0, 1], [1, 0]]))
            reaction = SimpleNamespace(
                scan_coo=[1, 2], maps=[np.asarray([0, 1]),
                                       np.asarray([2, 3])],
                equiv=[[1], [2]], products=[product_a, product_b],
                species=SimpleNamespace(
                    atom=['H'] * 4, geom=parent_geom.copy(), charge=0,
                    mult=1),
                irc_prod=SimpleNamespace(
                    geom=parent_geom.copy(), bondlist=[[0, 1], [2, 3]]),
                do_vdW=False)
            parameters = {
                'vrc_tst_scan_method': 'caspt2(2,2)',
                'vrc_tst_scan_basis': 'aug-cc-pVTZ',
                'vrc_tst_scan_points': [2.5],
            }
            vts = VTS(SimpleNamespace(), parameters, None)
            vts.scan_reac['r'] = reaction
            calls = []

            def optimize(reaction, geom, distance, logical_name):
                calls.append((distance, logical_name))
                return (Atoms(reaction.species.atom, positions=geom),
                        -272.1, Path('/verified/caspt2'),
                        {'status': 'executed'})

            with patch.object(
                    vts, '_molpro_multireference_scan_task', optimize):
                vts.do_scan(['r'])

            assert calls == [(2.5, 'r_vts_pt00')]
            db = connect('kinbot.db')
            for name in ('r_vts_pt00', 'r_vts_pt00_fr',
                         'r_vts_pt_asymptote',
                         'r_vts_pt_asymptote_fr'):
                row = list(db.select(name=f'vrctst/{name}', sort='-id'))[0]
                assert row.data['vrc_scan_backend'] == 'molpro'
                assert row.data['vrc_scan_method'] == 'caspt2(2,2)'
                assert row.data['vrc_scan_basis'] == 'aug-cc-pVTZ'
            relaxed = list(db.select(
                name='vrctst/r_vts_pt00', sort='-id'))[0].toatoms()
            assert np.linalg.norm(
                relaxed.positions[1] - relaxed.positions[2]) == \
                pytest.approx(2.5)
        finally:
            os.chdir(previous)


@pytest.mark.parametrize('fragment_mults', [(3, 2), (2, 1)])
def test_energy_anchor_geometry_cbs_and_nonzero_asymptotic_gap(tmp_path, monkeypatch, fragment_mults):
    """Synthetic energies exercise the actual dispatcher/assembler, not QC."""
    from kinbot import constants
    from kinbot.anl.extrapolation import two_point_cbs
    monkeypatch.chdir(tmp_path)
    Path('vrctst/molpro').mkdir(parents=True)
    distances = [2., 3., 4., 30.]
    symbols = ['C', 'S', 'H', 'H', 'H']
    database = connect('kinbot.db')
    geometries = {}
    for index, distance in enumerate(distances):
        base = 'r_vts_pt_asymptote' if index == 3 else f'r_vts_pt{index:02d}'
        geom = np.array([[0.,0.,0.],[distance,0.,0.],[distance,1.,0.],
                         [distance,-1.,0.],[0.,0.,1.]])
        frozen = geom.copy(); frozen[2,1] += .12
        for name, positions in [(base,geom),(base+'_fr',frozen)]:
            geometries[name] = positions
            database.write(Atoms(symbols, positions=positions), name='vrctst/'+name,
                           data={'status':'normal'})
    product_maps = ([0,3,4],[1,2]) if fragment_mults == (3,2) else ([0,4],[1,2,3])
    products = [SimpleNamespace(atom=[symbols[i] for i in mapping],
        geom=geom[mapping], mult=mult) for mapping,mult in zip(product_maps,fragment_mults)]
    parameters = dict(vrc_tst_scan_points=distances[:-1], vrc_tst_scan_molpro_key='MYENERGY',
        vrc_tst_correction_model='spin_splitting',
        vrc_tst_sample_method='caspt2-f12(5,4)', vrc_tst_sample_basis='avdz-f12',
        vrc_tst_state_average_roots={'2':4,'4':1},
        vrc_tst_spin_splitting={'reference_multiplicity':4,
            'gap_levels':[dict(label='cas',method='caspt2-f12(5,4)',
                basis=['avdz-f12','avtz-f12'],extrapolation_power=3.7)]})
    vts = VTS(SimpleNamespace(geom=geom), parameters, None)
    vts._surface_name = 'r'  # avoid plotting
    vts.scan_reac['r'] = SimpleNamespace(products=products, usym=[[[0]],[[1]]],
        equiv=[[0],[1]], maps=[np.array(m) for m in product_maps],
        species=SimpleNamespace(atom=symbols,charge=0,mult=2))
    values = {
        'sampling': [-9.98,-9.99,-9.995,-10.],
        'reference:0': [-9.80,-9.78,-9.77,-9.75],
        'reference:1': [-9.90,-9.88,-9.87,-9.85],
        'gap:cas:target:0': [-9.88,-9.86,-9.85,-9.83],
        'gap:cas:target:1': [-9.99,-9.965,-9.95,-9.92],
        'gap:cas:reference:0': [-9.60,-9.58,-9.57,-9.55],
        'gap:cas:reference:1': [-9.70,-9.68,-9.67,-9.65]}
    saved = {}
    class Native:
        def __init__(self, species, par): self.species = species
        def create_molpro_input(self, name, VTS, sample, vrc_level):
            index = 3 if 'asymptote' in name else int(name.split('_vts_pt')[1][:2])
            base = 'r_vts_pt_asymptote' if index == 3 else f'r_vts_pt{index:02d}'
            expected = base + ('_fr' if sample or index == 3 else '')
            assert self.species.geom == pytest.approx(geometries[expected])
            saved[name] = values[vrc_level['role']][index]
            Path(f'vrctst/molpro/{name}.inp').write_text(str(vrc_level))
        def get_molpro_energy(self, key, name, VTS):
            return (1,saved[name]) if Path(f'vrctst/molpro/{name}.out').exists() else (0,-1.)
    def dispatch(pending):
        from kinbot.input_refresh import file_sha256
        for name in pending:
            p = Path('vrctst/molpro')/name
            p.with_suffix('.out').write_text('Synthetic converged energy fixture')
            p.with_suffix('.input.sha256').write_text(file_sha256(p.with_suffix('.inp')))
    monkeypatch.setattr('kinbot.vrc_tst_scan.Molpro', Native)
    monkeypatch.setattr(vts, '_dispatch_molpro_corrections', dispatch)
    vts.energies(['r'])
    correction = json.loads(Path('vrctst/corr_r.json').read_text())
    def cbs(role):
        return np.array([two_point_cbs(a,b,upper_cardinal=3,power=3.7)
                         for a,b in zip(values[role+':0'],values[role+':1'])])
    high = cbs('reference') + cbs('gap:cas:target') - cbs('gap:cas:reference')
    assert correction['e_inf_high'] == pytest.approx(high[-1])
    assert correction['e_high'] == pytest.approx((high-high[-1])*constants.AUtoKCAL)
    assert correction['e_samp'] == pytest.approx((np.array(values['sampling'])+10.)*constants.AUtoKCAL)
    assert correction['frags_mult'] == list(fragment_mults)
    assert correction['vrc_molecule']['multiplicity'] == 2
    assert len(saved) == 28


def test_explicit_gaussian_scan_keeps_legacy_qc_route():
    geometry = np.asarray([
        [0., 0., 0.], [0., 0., .74],
        [0., 0., 3.], [0., 0., 3.74],
    ])
    reaction = SimpleNamespace(
        do_vdW=False, scan_coo=[1, 2], maps=[np.asarray([0, 1]),
                                             np.asarray([2, 3])],
        equiv=[[1], [2]],
        species=SimpleNamespace(geom=geometry.copy(), natom=4))
    qc = SimpleNamespace()
    qc.qc_vts = Mock(return_value='legacy-gaussian-job')
    qc.get_qc_geom = Mock(return_value=(0, geometry.copy()))
    qc.check_qc = Mock(return_value='normal')
    vts = VTS(SimpleNamespace(), {
        'vrc_tst_scan_method': 'ub3lyp',
        'vrc_tst_scan_points': [],
    }, qc)
    vts.scan_reac['r'] = reaction
    with patch('kinbot.vrc_tst_scan.time.sleep'):
        assert vts.do_scan(['r'], noscan=True) == ['legacy-gaussian-job']
    qc.qc_vts.assert_called_once()


def test_multireference_scan_rejects_internal_fragment_dissociation():
    reaction = SimpleNamespace(
        maps=[np.asarray([0, 1]), np.asarray([2, 3])],
        products=[
            SimpleNamespace(geom=np.asarray([[0., 0., 0.], [0., 0., .74]])),
            SimpleNamespace(geom=np.asarray([[0., 0., 0.], [0., 0., .74]])),
        ],
        irc_prod=SimpleNamespace(bondlist=[[0, 1], [2, 3]]))
    positions = np.asarray([
        [0., 0., 0.], [0., 0., .74],
        [0., 0., 3.], [0., 0., 4.2],
    ])
    broken = VTS._broken_scan_fragment_bonds(reaction, positions)
    assert len(broken) == 1
    assert broken[0][:2] == (2, 3)
    assert broken[0][2] == pytest.approx(.74)
    assert broken[0][3] == pytest.approx(1.2)


def test_spin_splitting_dispatches_all_cbs_components_and_records_provenance():
    with TemporaryDirectory() as temporary:
        previous = Path.cwd()
        os.chdir(temporary)
        try:
            Path('vrctst/molpro').mkdir(parents=True)
            atoms = Atoms('H2', positions=[[0., 0., 0.], [0., 0., 30.]])
            connect('kinbot.db').write(
                atoms, name='vrctst/r_vts_pt_asymptote_fr',
                data={'status': 'normal'})
            well = SimpleNamespace(
                chemid='parent', geom=np.asarray([[0., 0., 0.],
                                                   [0., 0., .74]]))
            parameters = {
                'vrc_tst_scan_points': [2.5],
                'vrc_tst_scan_molpro_key': 'MYENERGY',
                'vrc_tst_sample_method': 'caspt2(2,2)',
                'vrc_tst_sample_basis': 'vdz',
                'vrc_tst_high_method': 'mrci+q(2,2)',
                'vrc_tst_high_basis': 'vtz',
                'vrc_tst_correction_model': 'spin_splitting',
                # Explicit optional multi-method model: production defaults
                # now use one CASPT2 gap and a single-basis CC anchor.
                'vrc_tst_spin_splitting': {
                    'reference': {'method': 'uccsd(t)-f12b',
                                  'basis': ['avtz-f12', 'avqz-f12'],
                                  'extrapolation_power': 3.7},
                    'gap_levels': [
                        {'label': 'caspt2_f12', 'method': 'caspt2-f12(2,2)',
                         'basis': ['avtz-f12', 'avqz-f12'],
                         'extrapolation_power': 3.7, 'shift': .1, 'ipea': .25},
                        {'label': 'mrci_f12_q', 'method': 'mrci-f12+q(2,2)',
                         'basis': ['avtz-f12', 'avqz-f12'],
                         'extrapolation_power': 3.7}]},
            }
            vts = VTS(well, parameters, None)
            product_a = SimpleNamespace(
                atom=['H'], geom=[[0., 0., 0.]], mult=2)
            product_b = SimpleNamespace(
                atom=['H'], geom=[[0., 0., 30.]], mult=2)
            vts.scan_reac['r'] = SimpleNamespace(
                usym=[[[0]], [[0]]], equiv=[[0], [1]],
                maps=[np.asarray([0]), np.asarray([1])],
                products=[product_a, product_b],
                species=SimpleNamespace(mult=1))
            observed = set()

            def dispatch(pending):
                observed.update(pending)
                for name in pending:
                    Path(f'vrctst/molpro/{name}.out').write_text(
                        'Molpro calculation terminated\n')
                    source = Path(f'vrctst/molpro/{name}.inp')
                    Path(f'vrctst/molpro/{name}.input.sha256').write_text(
                        __import__('hashlib').sha256(
                            source.read_bytes()).hexdigest() + '\n')

            with patch('kinbot.vrc_tst_scan.Molpro', _FakeMolpro), \
                    patch.object(vts, '_dispatch_molpro_corrections', dispatch):
                vts.energies(['r'], noscan=True)

            assert len(observed) == 11
            assert 'r_vts_pt_asymptote_fr' in observed
            assert 'r_vts_pt_asymptote_ss_triplet_cc_1' in observed
            assert 'r_vts_pt_asymptote_ss_mrci_f12_q_triplet_1' in observed
            correction = json.loads(
                Path('vrctst/corr_r.json').read_text())
            assert correction['correction_model'] == 'spin_splitting'
            assert correction['levels']['trusted_correction']['model'] == \
                'spin_splitting'
            assert correction['spin_splitting']['triplet_cbs']['power'] == \
                pytest.approx(3.7)
            assert len(correction['native_components']) == 11
            mrci_native = correction['native_components'][
                'gap:mrci_f12_q:singlet:1'][0]
            assert mrci_native['energy_contract']['method'] == \
                'MRCI-F12+Q'
            assert mrci_native['energy_contract']['energy_variable'] == \
                'ENERGD(1)'
            assert mrci_native['energy_contract'][
                'davidson_correction'] == 'relaxed_reference'
        finally:
            os.chdir(previous)


def test_molpro_dispatch_normalizes_numpy_charge_and_multiplicity():
    with TemporaryDirectory() as temporary:
        previous = Path.cwd()
        os.chdir(temporary)
        try:
            Path('vrctst/molpro').mkdir(parents=True)
            Path('vrctst/molpro/r.inp').write_text(
                '***,VRC correction\nMolpro calculation terminated\n')
            species = SimpleNamespace(
                atom=['H', 'H'],
                geom=np.asarray([[0., 0., 0.], [0., 0., .74]]),
                charge=np.float64(0.), mult=np.int64(1))
            parameters = {
                'vrc_tst_walltime': '00:10:00',
                'single_point_ppn': 4,
                'vrc_tst_min_stack_mw': 50,
                'queue_name': 'short',
                'vrc_tst_max_nodes': 1,
            }
            vts = VTS(SimpleNamespace(), parameters, None)
            captured = {}

            def prepare_dispatch(spec, run_dir):
                from kinbot.anl.dispatch import _atoms
                _atoms(spec['molecule'])
                json.dumps(spec)
                captured['spec'] = spec
                output = run_dir / 'tasks' / 'r' / 'r.out'
                output.parent.mkdir(parents=True)
                output.write_text('Molpro calculation terminated\n')

            with patch('kinbot.anl.dispatch.prepare', prepare_dispatch), \
                    patch('kinbot.anl.dispatch.preflight'), \
                    patch('kinbot.anl.dispatch.drive_run',
                          return_value=(0, {'r': 'complete'})):
                vts._dispatch_molpro_corrections({'r': species})

            molecule = captured['spec']['molecule']
            assert type(molecule['charge']) is int
            assert type(molecule['multiplicity']) is int
            assert captured['spec']['intent'][
                'cap_walltime_to_explicit_partition'] is True
            assert Path('vrctst/molpro/r.out').read_text() == \
                'Molpro calculation terminated\n'
        finally:
            os.chdir(previous)


@pytest.mark.parametrize('fmax', [None, 2.e-4])
@pytest.mark.parametrize('fragment_retry', [False, True])
@pytest.mark.parametrize('walltime', ['01:00:00', 'partition_max'])
def test_relaxed_scan_dispatch_receives_force_and_electronic_thresholds(
        tmp_path, monkeypatch, fmax, fragment_retry, walltime):
    monkeypatch.chdir(tmp_path)
    geom = np.asarray([[0., 0., 0.], [0., 0., 2.5]])
    reaction = SimpleNamespace(
        scan_coo=[0, 1], equiv=[[0], [1]], products=[],
        species=SimpleNamespace(atom=['H', 'H'], charge=0, mult=1))
    par = dict(vrc_tst_scan_method='caspt2(2,2)', vrc_tst_scan_basis='vdz',
               vrc_tst_walltime=walltime, single_point_ppn=4,
               vrc_tst_min_stack_mw=50, queue_name='short')
    if fmax is not None:
        par['vrc_tst_scan_fmax'] = fmax
    captured = {}

    class Prepared(Exception):
        pass

    def prepare(spec, run_dir):
        captured.update(spec)
        raise Prepared

    monkeypatch.setattr('kinbot.anl.dispatch.prepare', prepare)
    with pytest.raises(Prepared):
        VTS(SimpleNamespace(), par, None)._molpro_multireference_scan_task(
            reaction, geom, 2.5, 'test_point', fragment_retry=fragment_retry)
    task = captured['tasks'][0]
    assert task['optimizer']['fmax'] == (5.e-4 if fmax is None else fmax)
    assert task['optimizer']['constraints']['fixed_bonds'] == [[0, 1]]
    kw = task['profile']['calculator_kwargs']
    assert kw['scf_convergence'] == 1.e-8
    assert kw['correlation_convergence'] == 1.e-7
    if fragment_retry:
        from kinbot.run_walltime import seconds
        assert seconds(task['resources']['walltime']) == (3600 if walltime == '01:00:00' else 7200)
        assert task['optimizer']['step_retries'] == 0
        assert task['optimizer']['steps'] == 100
        assert captured['intent']['vrc_fragment_retry'] is True
    else:
        assert task['resources']['walltime'] == walltime


def test_molpro_dispatch_surfaces_native_task_failure():
    with TemporaryDirectory() as temporary:
        previous = Path.cwd()
        os.chdir(temporary)
        try:
            Path('vrctst/molpro').mkdir(parents=True)
            Path('vrctst/molpro/r.inp').write_text('***,failed VRC point\n')
            species = SimpleNamespace(
                atom=['H', 'H'], geom=[[0., 0., 0.], [0., 0., .74]],
                charge=0., mult=1)
            parameters = {
                'vrc_tst_walltime': '00:10:00', 'single_point_ppn': 2,
                'vrc_tst_min_stack_mw': 50, 'queue_name': 'short',
                'vrc_tst_max_nodes': 1,
            }
            vts = VTS(SimpleNamespace(), parameters, None)
            state = {'tasks': {'r': {'status': 'failed'}}}

            def prepare_dispatch(_, run_dir):
                task_dir = run_dir / 'tasks' / 'r'
                task_dir.mkdir(parents=True)
                (task_dir / 'execution.json').write_text(json.dumps({
                    'status': 'failed',
                    'error': 'PMPI_Init: OFI endpoint open failed',
                }))

            with patch('kinbot.anl.dispatch.prepare', prepare_dispatch), \
                    patch('kinbot.anl.dispatch.preflight'), \
                    patch('kinbot.anl.dispatch._load',
                          side_effect=lambda run_dir: (Path(run_dir), {}, state)), \
                    patch('kinbot.anl.dispatch.drive_run',
                          return_value=(1, {'r': 'failed'})):
                with pytest.raises(RuntimeError, match='PMPI_Init.*OFI endpoint'):
                    vts._dispatch_molpro_corrections({'r': species})
        finally:
            os.chdir(previous)


def test_vrc_dispatch_retries_one_execution_failure_and_preserves_archive(
        tmp_path):
    failed = {
        'tasks': {'relaxed_scan': {
            'status': 'failed', 'attempt': 1,
            'failure_phase': 'execution'}}}
    complete = {'tasks': {'relaxed_scan': {'status': 'complete'}}}
    archive = tmp_path / 'attempts' / 'relaxed_scan' / '1'
    archive.mkdir(parents=True)
    with patch('kinbot.anl.dispatch.drive_run', side_effect=[
            (1, {'relaxed_scan': 'failed'}),
            (0, {'relaxed_scan': 'complete'})]), \
            patch('kinbot.anl.dispatch._load', return_value=(
                tmp_path, {}, failed)), \
            patch('kinbot.anl.dispatch.retry_failed',
                  return_value=archive) as retry, \
            patch('kinbot.anl.dispatch.preflight_until_scheduler_ready') \
            as preflight:
        result = _drive_vrc_dispatch(
            tmp_path, ('relaxed_scan',), interval=0.)
    assert result == (0, {'relaxed_scan': 'complete'})
    retry.assert_called_once_with(tmp_path, 'relaxed_scan')
    preflight.assert_called_once_with(tmp_path, interval=0.)


def test_vrc_dispatch_does_not_retry_validation_failure(tmp_path):
    task = tmp_path / 'tasks' / 'relaxed_scan'
    task.mkdir(parents=True)
    (task / 'execution.json').write_text(json.dumps({
        'status': 'failed', 'failure_phase': 'validation',
        'error': 'constraint was lost'}))
    failed = {'tasks': {'relaxed_scan': {
        'status': 'failed', 'attempt': 1}}}
    with patch('kinbot.anl.dispatch.drive_run', return_value=(
            1, {'relaxed_scan': 'failed'})), \
            patch('kinbot.anl.dispatch._load', return_value=(
                tmp_path, {}, failed)), \
            patch('kinbot.anl.dispatch.retry_failed') as retry:
        result = _drive_vrc_dispatch(
            tmp_path, ('relaxed_scan',), interval=0.)
    assert result == (1, {'relaxed_scan': 'failed'})
    retry.assert_not_called()


def test_vrc_dispatch_retries_each_task_at_most_once(tmp_path):
    failed = {'tasks': {'relaxed_scan': {
        'status': 'failed', 'attempt': 2,
        'failure_phase': 'scheduler'}}}
    with patch('kinbot.anl.dispatch.drive_run', return_value=(
            1, {'relaxed_scan': 'failed'})), \
            patch('kinbot.anl.dispatch._load', return_value=(
                tmp_path, {}, failed)), \
            patch('kinbot.anl.dispatch.retry_failed') as retry:
        result = _drive_vrc_dispatch(
            tmp_path, ('relaxed_scan',), interval=0.)
    assert result == (1, {'relaxed_scan': 'failed'})
    retry.assert_not_called()


def test_asymptotic_mr_geometry_feeds_fixed_pivots_and_every_frozen_point(tmp_path, monkeypatch):
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.conformers import Conformers
    from kinbot.parameters import Parameters
    from kinbot.qc import QuantumChemistry
    monkeypatch.chdir(tmp_path)
    Path('profiles.json').write_text(json.dumps({
        'theory_preset': 'b3lyp-b2plyp-anl', 'barrier_threshold': 50.}))
    profiles = Parameters('profiles.json', show_warnings=False).par
    qc = QuantumChemistry(profiles)
    geometry = np.array([[0., 0., 0.], [0., 0., .74],
                         [0., 0., 3.24], [0., 0., 3.98]])
    products = [StationaryPoint('h2', 0, 1, atom=['H', 'H'], natom=2,
                geom=geometry[part].copy())
                for part in (slice(0, 2), slice(2, 4))]
    for product in products:
        product.bond01 = np.array([[0, 1], [1, 0]])
        product.confs = Conformers(product, profiles, qc)
    def reject_service_copy(*args):
        pytest.fail('Asymptotic extraction copied a live conformer/QC controller')
    monkeypatch.setattr(Conformers, '__deepcopy__', reject_service_copy, raising=False)
    reaction = SimpleNamespace(scan_coo=[1, 2], maps=[[0, 1], [2, 3]],
        equiv=[[1], [2]], products=products, do_vdW=False,
        species=SimpleNamespace(atom=['H']*4, geom=geometry.copy(), charge=0, mult=1),
        irc_prod=SimpleNamespace(geom=geometry.copy(), bondlist=[[0, 1], [2, 3]]))
    vts = VTS(SimpleNamespace(), dict(vrc_tst_fragment_geometry='asymptotic_scan',
        vrc_tst_scan_method='caspt2(2,2)', vrc_tst_scan_basis='avdz-f12',
        vrc_tst_scan_points=[2.5, 2.7], vrc_tst_asymptote_distance=30.), None)
    vts.scan_reac['r'] = reaction
    optimized, pivots = [], []
    def optimize(reaction, geom, distance, logical_name, **kwargs):
        optimized.append((distance, logical_name))
        a = Atoms('H4', positions=geom)
        # At the asymptote, change only internal fragment bond lengths.
        if distance == 30.:
            a.positions[0] = a.positions[1] - [0, 0, .80]
            a.positions[3] = a.positions[2] + [0, 0, .82]
        return a, -2., tmp_path/logical_name, {'details': {'final_geometry_sha256': 'accepted'}}, {}
    def pivot(product, par, destination, evidence):
        pivots.append((np.linalg.norm(product.geom[1]-product.geom[0]), evidence))
        return str(destination.with_suffix('.cube'))
    monkeypatch.setattr(vts, '_molpro_multireference_scan_task', optimize)
    monkeypatch.setattr('kinbot.vrc_orbitals.prepare_fixed_pivot', pivot)
    vts.do_scan(['r'])
    assert optimized == [(30., 'r_vts_fragments_asymptote'), (2.5, 'r_vts_pt00'), (2.7, 'r_vts_pt01')]
    assert [item[0] for item in pivots] == pytest.approx([.80, .82])
    assert products[0].geom == pytest.approx(geometry[:2])  # source objects preserved
    assert products[1].geom == pytest.approx(geometry[2:])
    for original, fragment in zip(products, reaction.products):
        assert fragment.confs is original.confs
        assert fragment.confs.qc is qc
    assert reaction.fragment_orbital_cubes == ['vrctst/pivot_cubes/r/fragment_0.cube',
                                               'vrctst/pivot_cubes/r/fragment_1.cube']
    db = connect('kinbot.db')
    for point in ('pt00', 'pt01', 'pt_asymptote'):
        a = list(db.select(name=f'vrctst/r_vts_{point}_fr'))[0].toatoms()
        assert a.get_distance(0, 1) == pytest.approx(.80)
        assert a.get_distance(2, 3) == pytest.approx(.82)
