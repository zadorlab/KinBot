"""Production recipes, generated native decks and assembly must agree."""
from copy import deepcopy
from dataclasses import replace

import pytest
from ase import Atoms

from kinbot.anl.dispatch import validate_spec, _geometry_hash
from kinbot.anl.recipes import recipe
from kinbot.anl.validation import automatic_composite_spec
from kinbot.anl.results import validate_result_parser, parse_cfour_dboc
from kinbot.ase_modules.calculators.molpro import render_input


def molecule(mult=2):
    return {'symbols': ['S', 'H'], 'positions': [[0., 0., 0.], [0., 0., 1.34]],
            'charge': 0 if mult == 2 else 1, 'multiplicity': mult,
            'source': {'geometry_sha256': _geometry_hash(Atoms('SH', positions=[[0, 0, 0], [0, 0, 1.34]])), 'record_sha256': 'b'*64,
                       'l2_theory': {'calculator': 'gaussian', 'method': 'B2PLYP',
                                     'basis': 'cc-pVTZ', 'dispersion': 'GD3BJ'},
                       'energy_hartree': -398.}}


@pytest.mark.parametrize('tier,geom,dboc,dbasis', [
    ('ANL0', 'CCSD(T)', 'HF', 'D'),
    ('ANL0-F12', 'CCSD(T)-F12b', 'CCSD', 'D'),
    ('ANL1', 'CCSD(T)', 'HF', 'T'),
    ('ANL1-F12', 'CCSD(T)-F12b', 'CCSD', 'T'),
])
@pytest.mark.parametrize('augmented', [False, True])
@pytest.mark.parametrize('mult', [1, 2])
def test_all_eight_variants_generate_consistent_graphs(tier, geom, dboc, dbasis, augmented, mult):
    spec = automatic_composite_spec(molecule(mult), highest_tier=tier, augmented=augmented)
    tasks = {t['id']: t for t in spec['tasks']}
    resolved = deepcopy(spec)
    for t in resolved['tasks']:
        t['resources'].update(cores=2, memory_mb=32000, partition='test')
    validate_spec(resolved)
    req = {r.key: r for r in recipe(tier, component_version=2,
                                   augmented=augmented, multiplicity=mult).requirements}
    a = 'aug-' if augmented else ''
    f = '-F12' if tier.endswith('F12') else ''
    z = 'Q' if tier.startswith('ANL1') else 'T'
    assert tasks['l3_geometry']['profile']['method'] == geom
    assert tasks['l3_geometry']['profile']['basis'] == f'{a}cc-pV{z}Z{f}'
    assert req['reference_cbs'].settings['geometry_method'] == geom
    assert req['reference_cbs'].settings['geometry_basis'] == f'{a}cc-pV{z}Z{f}'
    if tier.startswith('ANL1'):
        assert tasks['l3_geometry']['geometry_from'] == 'l3_geometry_tz'
        assert tasks['l3_geometry_tz']['profile']['basis'] == f'{a}cc-pVTZ{f}'
        assert tasks['harmonic_tz']['geometry_from'] == 'l3_geometry_tz'
        assert tasks['harmonic_qz']['geometry_from'] == 'l3_geometry'
    h = tasks['harmonic_qz' if tier.startswith('ANL1') else 'harmonic']
    assert h['result_parser'].get('method', 'CCSD(T)') == geom
    assert req['harmonic_zpe'].method == geom
    assert tasks['cfour_dboc']['result_parser']['level'] == dboc
    assert tasks['cfour_dboc']['result_parser']['basis'] == f'{a}cc-pV{dbasis}Z'
    assert req['dboc'].method == dboc
    assert req['dboc'].basis == f'{a}cc-pV{dbasis}Z'
    assert tasks['cfour_dboc']['result_parser']['reference'] == ('RHF' if mult == 1 else 'UHF' if tier.endswith('-F12') else 'ROHF')
    cv_cardinals = ('Q', '5') if z == 'Q' else ('T', 'Q')
    for ident in (f'cv_{core}_{n.lower()}z' for core in ('ae', 'fc') for n in cv_cardinals):
        assert tasks[ident]['result_parser']['method'] == geom
        cardinal = ident.split('_')[-1][0].upper()
        assert tasks[ident]['result_parser']['basis'] == f'{a}cc-pwCV{cardinal}Z'
        assert f'basis={a}cc-pwCV{cardinal}Z' in tasks[ident]['input_template']
        assert 'geometry_fallback' not in tasks[ident]
    assert req['core_valence_cbs'].settings['native_basis_pair'] == (
        tuple(f'{a}cc-pwCV{n}Z' for n in cv_cardinals))
    assert req['core_valence_cbs'].settings['upper_cardinal'] == (5 if z == 'Q' else 4)
    rel_ids = ('rel_dkh_qz', 'rel_dkh_5z', 'rel_nonrel_qz', 'rel_nonrel_5z') \
        if z == 'Q' else ('rel_dkh', 'rel_nonrel')
    for ident in rel_ids:
        assert tasks[ident]['result_parser']['method'] == 'CCSD(T)'
        assert '-f12' not in tasks[ident]['input_template'].lower()
        basis = f'{a}cc-pwCV{ident[-2].upper()}Z-DK' if z == 'Q' else f'{a}cc-pCVTZ-DK'
        assert tasks[ident]['result_parser']['basis'] == basis
        if 'nonrel' in ident:
            assert 'set,dkroll=0' in tasks[ident]['input_template']
        assert tasks[ident]['result_parser']['explicit_hamiltonian'] is True
    for cardinal in (('qz', '5z') if z == 'Q' else ('',)):
        suffix = '_'+cardinal if cardinal else ''
        dkh = tasks['rel_dkh'+suffix]['input_template']
        nr = tasks['rel_nonrel'+suffix]['input_template']
        # Exactly matched science: independent artifact names and the
        # Hamiltonian control are the only changes between these decks.
        assert dkh.replace('rel_dkh', 'rel_pair').replace('set,dkho=2', 'HAMILTONIAN') == (
            nr.replace('rel_nonrel', 'rel_pair').replace('set,dkroll=0', 'HAMILTONIAN'))
    for ident in ('ccsdt_dz', 'ccsdt_tz', 'ccsdtq_dz', 'ccsdtq_tz', 'ccsdtqp_dz'):
        if ident in tasks:
            assert '-f12' not in tasks[ident]['input_template'].lower()
    if f:
        # Deliberately preserve the independently extrapolated conventional arm.
        ids = ('anl1_5z', 'anl1_6z') if z == 'Q' else ('anl0_qz', 'anl0_5z')
        assert all(tasks[i]['result_parser']['method'] == 'CCSD(T)' for i in ids)
        assert req['reference_cbs'].settings['component_weights'] == (0.5, 0.5)


@pytest.mark.parametrize('symbols,mult,effective', [
    (['H'], 2, 'HF'), (['H', 'H'], 1, 'CCSD-F12b')])
def test_f12_rank_exact_geometry_and_frequencies(symbols, mult, effective):
    m = molecule(mult)
    m.update(symbols=symbols, positions=[[0., 0., i*.74] for i in range(len(symbols))], charge=0)
    m['source']['geometry_sha256'] = _geometry_hash(Atoms(symbols, positions=m['positions']))
    spec = automatic_composite_spec(m, highest_tier='ANL1-F12')
    for t in spec['tasks']:
        t['resources'].update(cores=2, memory_mb=32000, partition='test')
    validate_spec(spec)
    tasks = {t['id']: t for t in spec['tasks']}
    assert tasks['l3_geometry']['profile']['calculator_kwargs']['effective_method'] == effective
    assert not any(t.startswith('cv_') for t in tasks)


def test_f12_forces_differentiate_the_precise_scaled_f12_energy():
    deck = render_input(Atoms('SH', positions=[[0, 0, 0], [0, 0, 1.34]]),
                        basis='aug-cc-pVTZ-F12', geometry_name='geometry.xyz',
                        method='CCSD(T)-F12b', mult=2)
    assert 'uccsd(t)-f12b,scale_trip=1' in deck
    assert 'show[1,f28.16],kb_geom_energy' in deck
    assert 'forces,numerical,variable=kb_geom_energy,startcmd=rhf' in deck


def test_ccsd_dboc_parser_cannot_accept_hf_input_or_summary():
    from test_anl_results import _CFOUR_CH4
    request = {'kind': 'cfour_dboc', 'file': 'cfour.out', 'level': 'CCSD'}
    with pytest.raises(ValueError):
        validate_result_parser(request, backend='cfour', outputs=['cfour.out'],
                               template='*CFOUR(CALC=SCF\nDBOC=ON)')
    with pytest.raises(ValueError, match='CCSD'):
        parse_cfour_dboc(_CFOUR_CH4, level='CCSD')
    cc = _CFOUR_CH4.replace('at Hartree-Fock level', 'at CCSD level')
    assert parse_cfour_dboc(cc, level='CCSD')['selected_level'] == 'CCSD'


def test_archived_component_version_is_not_reinterpreted():
    old = automatic_composite_spec(molecule(), highest_tier='ANL1-F12', component_version=1)
    new = automatic_composite_spec(molecule(), highest_tier='ANL1-F12')
    a, b = ({t['id']: t for t in s['tasks']} for s in (old, new))
    assert a['l3_geometry']['profile']['method'] == 'CCSD(T)'
    assert b['l3_geometry']['profile']['method'] == 'CCSD(T)-F12b'
    assert a['cfour_dboc']['result_parser']['level'] == 'HF'
    assert b['cfour_dboc']['result_parser']['level'] == 'CCSD'
    assert recipe('ANL1-F12', component_version=1).version == 1
    assert recipe('ANL1-F12', component_version=2).version == 2


@pytest.mark.parametrize('tier', ['ANL1', 'ANL1-F12'])
def test_q5_corrections_use_matching_hamiltonian_pairs_and_correct_cbs_power(monkeypatch, tier):
    from test_anl_workflow import _electronic_component
    from kinbot.anl.workflow import core_valence_task_component, scalar_relativistic_task_component
    from kinbot.anl.extrapolation import two_point_cbs
    requirements = {r.key: r for r in recipe(tier, component_version=2, augmented=True).requirements}
    parts = {}
    for cardinal, z in [('Q', 'qz'), ('5', '5z')]:
        for arm, energy in [('dkh', -40.12 if z == 'qz' else -40.223),
                            ('nonrel', -40.10 if z == 'qz' else -40.20)]:
            name = f'rel_{arm}_{z}'
            parts[name] = _electronic_component(name, energy, f'aug-cc-pwCV{cardinal}Z-DK',
                core='all-electron', relativistic='DKH2' if arm == 'dkh' else 'none')
        for core, energy in [('ae', -40.12 if z == 'qz' else -40.231),
                             ('fc', -40.10 if z == 'qz' else -40.20)]:
            name = f'cv_{core}_{z}'
            part = _electronic_component(name, energy, f'aug-cc-pwCV{cardinal}Z',
                                        core='all-electron' if core == 'ae' else 'frozen')
            if tier.endswith('-F12'):
                part = replace(part, method='CCSD(T)-F12b', settings={**part.settings,
                    'scale_trip': 1, 'f12_auxiliary': requirements['core_valence_cbs'].settings[
                        'f12_auxiliary_by_basis'][part.basis]})
            parts[name] = part
    monkeypatch.setattr('kinbot.anl.workflow.task_component', lambda _, name, **kwargs: parts[name])
    rel = scalar_relativistic_task_component('run', ('rel_dkh_qz', 'rel_dkh_5z'),
        ('rel_nonrel_qz', 'rel_nonrel_5z'), requirement=requirements['scalar_relativistic'], state_id='state-A')
    assert rel.value_hartree == pytest.approx(two_point_cbs(-.02, -.023, upper_cardinal=5, power=3.7))
    cv = core_valence_task_component('run', all_electron_lower='cv_ae_qz', all_electron_upper='cv_ae_5z',
        frozen_core_lower='cv_fc_qz', frozen_core_upper='cv_fc_5z',
        requirement=requirements['core_valence_cbs'], state_id='state-A')
    assert cv.value_hartree == pytest.approx(two_point_cbs(-.02, -.031, upper_cardinal=5, power=3.7))
    parts['rel_nonrel_5z'] = replace(parts['rel_nonrel_5z'], basis='aug-cc-pwCV5Z')
    with pytest.raises(ValueError, match='all-electron Molpro'):
        scalar_relativistic_task_component('run', ('rel_dkh_qz', 'rel_dkh_5z'),
            ('rel_nonrel_qz', 'rel_nonrel_5z'), requirement=requirements['scalar_relativistic'], state_id='state-A')


def test_ccsd_dboc_obeys_production_convergence_controls():
    spec = automatic_composite_spec(molecule(), highest_tier='ANL1-F12', cc_convergence=1e-7)
    task = next(t for t in spec['tasks'] if t['id'] == 'cfour_dboc')
    assert 'CC_CONV=7\n' in task['input_template']
    assert 'LINEQ_CONV=7\n' in task['input_template']


def test_one_electron_cv_reduction_needs_no_fitting_but_two_electrons_do(monkeypatch):
    from test_anl_workflow import _electronic_component
    from kinbot.anl.workflow import core_valence_task_component
    requirement = next(r for r in recipe('ANL1-F12', component_version=2).requirements
                       if r.key == 'core_valence_cbs')
    parts = {}
    for core in ('ae', 'fc'):
        for z in ('qz', '5z'):
            ident = f'cv_{core}_{z}'
            p = _electronic_component(ident, -10.0, 'cc-pwCV'+z.upper(),
                core='all-electron' if core == 'ae' else 'frozen')
            parts[ident] = replace(p, method='CCSD(T)-F12b', settings={**p.settings,
                'rank_exact': {'electron_count': 1, 'effective_method': 'HF'},
                'effective_program_variant': 'RHF', 'scale_trip': 1})
    monkeypatch.setattr('kinbot.anl.workflow.task_component', lambda _, name, **kw: parts[name])
    def assemble():
        return core_valence_task_component('run', all_electron_lower='cv_ae_qz',
            all_electron_upper='cv_ae_5z', frozen_core_lower='cv_fc_qz',
            frozen_core_upper='cv_fc_5z', requirement=requirement, state_id='state-A')
    assert assemble().value_hartree == 0
    for ident, p in parts.items():
        parts[ident] = replace(p, settings={**p.settings,
            'rank_exact': {'electron_count': 2, 'effective_method': 'CCSD-F12b'},
            'effective_program_variant': 'RHF-UCCSD-F12b'})
    with pytest.raises(ValueError, match='auxiliaries'):
        assemble()


@pytest.mark.parametrize('bad', ['missing', 'before_basis', 'opposite', 'conflicting'])
def test_dk_basis_requires_explicit_nonrel_control_in_deck_and_native_output(bad):
    from kinbot.anl.results import _check_explicit_hamiltonian
    task = next(t for t in automatic_composite_spec(molecule(), highest_tier='ANL1')['tasks']
                if t['id'] == 'rel_nonrel_qz')
    text = task['input_template']
    if bad == 'missing':
        text = text.replace('set,dkroll=0\n', '')
    elif bad == 'before_basis':
        text = 'set,dkroll=0\n'+text.replace('set,dkroll=0\n', '')
    elif bad == 'opposite':
        text = text.replace('set,dkroll=0', 'set,dkho=2')
    else:
        text = text.replace('set,dkroll=0', 'set,dkroll=0\nset,dkho=2')
    with pytest.raises(ValueError, match='explicit requested Hamiltonian'):
        _check_explicit_hamiltonian(text, task['result_parser'])
    with pytest.raises(ValueError):
        validate_result_parser(task['result_parser'], backend='molpro',
                               outputs=task['required_outputs'], template=text)


@pytest.mark.parametrize('heading', ['UCCSD(T)', 'UCCSD(T)-F12', 'UCCSD(T)-F12B'])
def test_f12_harmonic_reads_native_driver_label_but_requires_f12b_result(heading):
    from test_anl_results import _MOLPRO_HARMONIC
    from kinbot.anl.results import parse_molpro_harmonic
    native = (_MOLPRO_HARMONIC.replace('uccsd(t)\n', 'uccsd(t)-f12b,scale_trip=1\n')
        .replace('Triples (T) contribution ', 'Triples (T) contribution (scaled) ')
        .replace('for UCCSD(T))', f'for {heading})')
        .replace(' PROGRAM * FREQUENCIES',
                 ' RHF-UCCSD-F12b energy -40.1\n !RHF-UCCSD(T)-F12 energy -40.106123456789\n PROGRAM * FREQUENCIES'))
    result = parse_molpro_harmonic(native, basis='cc-pVTZ', reference='RHF', method='CCSD(T)-F12b')
    assert result['method'] == 'CCSD(T)-F12b'
    with pytest.raises(ValueError, match='F12b energy'):
        parse_molpro_harmonic(native.replace('!RHF-UCCSD(T)-F12 energy', '!RHF-UCCSD(T) energy'),
                             basis='cc-pVTZ', reference='RHF', method='CCSD(T)-F12b')
