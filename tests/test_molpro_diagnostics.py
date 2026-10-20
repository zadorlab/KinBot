"""Native evidence extraction and automatic runtime coverage, not QC simulation."""
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

from ase import Atoms
import pytest

from kinbot.molpro_diagnostics import (diagnostic_options, ensure_native_printing,
    native_print_directives, parse_native_output, write_native_report)
from kinbot.molpro import _vrc_multireference_method, validate_vrc_state_model

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('case,expected', [('ch3s', 10), ('ch3ss', 22)])
def test_retained_native_caspt2_outputs_have_extractable_evidence(case, expected):
    archive = ROOT / f'examples/anl/{case}_anl1_f12/orbital_reviews/caspt2_native.evidence.tar.gz'
    checked = 0
    with tarfile.open(archive) as bundle:
        for member in bundle.getmembers():
            if not member.name.endswith('.out'):
                continue
            native = bundle.extractfile(member).read().decode()
            inp = bundle.extractfile(member.name[:-4]+'.native.inp').read().decode()
            report = parse_native_output(native, inp)
            assert report['issues'] == [], (member.name, report['issues'])
            assert report['diagnostic_status'] == 'evidence_recorded'
            assert report['identity_verdict'] == 'requires_orbital_and_state_overlap_review'
            assert report['native_geometries'] and report['native_basis_evidence']
            multi = report['native_sections'][0]
            assert multi['program'] == 'MULTI'
            natural = next(t for t in multi['orbital_tables'] if t['heading'] == 'NATURAL ORBITALS (STATE AVERAGED)')
            active = [o for o in natural['orbitals'] if o['role'] == 'active']
            model = report['requested_casscf'][0]
            assert [o['index'] for o in active] == model['destination_active_indices']
            assert len(active) == (4 if case == 'ch3s' else 5)
            assert all(o['ao_coefficients'] for o in active)
            assert len(multi['state_energies']) == sum(s['roots'] for s in model['states'])
            assert multi['ci_tables']
            for row in active:
                assert native.splitlines()[row['line']-1].strip().startswith(row['index'])
            checked += 1
    assert checked == expected


def real_case():
    archive = ROOT / 'examples/anl/ch3s_anl1_f12/orbital_reviews/caspt2_native.evidence.tar.gz'
    with tarfile.open(archive) as bundle:
        member = next(m for m in bundle.getmembers() if m.name.endswith('.out'))
        return (bundle.extractfile(member).read().decode(),
                bundle.extractfile(member.name[:-4]+'.native.inp').read().decode())


def test_reports_are_native_snapshots_and_detect_missing_or_changed_partitions(tmp_path):
    text, inp = real_case()
    output = tmp_path / 'sample.out'
    source = output.with_suffix('.inp')
    output.write_text(text)
    source.write_text(inp)
    path = write_native_report(output)
    data = json.loads(path.read_text())
    assert data['sources']['output']['sha256'] == hashlib.sha256(output.read_bytes()).hexdigest()
    assert data['sources']['input']['sha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert data['diagnostic_status'] == 'evidence_recorded'
    assert data['requested_correlation_calls']
    bad = parse_native_output(text, inp.replace('closed,10', 'closed,9'))
    assert any('counts disagree' in s for s in bad['issues'])
    assert any('partitions disagree' in s for s in bad['issues'])
    bad = parse_native_output(text.replace('NATURAL ORBITALS (STATE AVERAGED)', 'omitted'), inp)
    assert any('occupation table missing' in s for s in bad['issues'])
    output.write_text(' PROGRAM * MULTI\n ? Error\n GLOBAL ERROR fehler\n')
    write_native_report(output)
    bad = json.loads(path.read_text())
    assert bad['diagnostic_status'] == 'review_required'
    assert bad['native_error_lines'] and not bad['normal_termination_seen']
    assert bad['sources']['output']['sha256'] != data['sources']['output']['sha256']
    assert not list(tmp_path.glob('*.tmp'))


def test_repeated_native_multi_calls_are_not_silently_discarded():
    text, inp = real_case()
    one = parse_native_output(text, inp)
    many = parse_native_output(text + text, inp)
    assert many['issues'] == []
    assert len(many['native_sections']) == 2*len(one['native_sections'])
    assert [s['requested_model'] for s in many['native_sections'] if s['program'] == 'MULTI'] == [0, 0]


@pytest.mark.parametrize('options', [{'orbital_coefficient_threshold': 0},
    {'ci_coefficient_threshold': float('nan')}, {'virtual_orbitals': True},
    {'virtual_orbitals': -1}, {'energy_convergence': 1e-8}])
def test_invalid_diagnostic_controls(options):
    with pytest.raises(ValueError, match='orbital_diagnostics'):
        diagnostic_options(options)


@pytest.mark.parametrize('method', ['CASPT2(5,4)', 'CASPT2-F12(5,4)',
                                   'MRCI+Q(5,4)', 'MRCI-F12+Q(5,4)'])
def test_same_print_contract_for_general_mr_methods(method):
    model = {'state_average': [{'multiplicity': 2, 'roots': 1}],
             'selected_roots': {'2': 1},
             'orbital_diagnostics': {'orbital_coefficient_threshold': 0.02,
                                     'ci_coefficient_threshold': 0.03,
                                     'virtual_orbitals': 3}}
    deck, selected = _vrc_multireference_method(method, 25, 1, state_model=model)
    assert 'gprint,orbitals,orben,civector,ref' in deck
    assert 'gthresh,thrprint=0.02,printci=0.03' in deck and 'orbprint,3' in deck
    assert 'thrden=1.d-7,thrvar=1.d-10' in deck
    if method.startswith('MRCI'):
        assert 'KBVRCSEL=energd(1)' in deck
    else:
        assert 'shift=0.1,ipea=0.25' in deck
    assert validate_vrc_state_model(model)['orbital_diagnostics'] == model['orbital_diagnostics']


def test_dispatcher_custom_decks_get_printing_without_overriding_explicit_options():
    from kinbot.anl.dispatch import _render_input
    template = ('gthresh,energy=1.d-10,thrprint=0.1,printci=0.02\n'
                '{multi;occ,2;closed,0;wf,2,1,0;}\n{mrci;}\n')
    result = _render_input(template, Atoms('H2'), {}, {'cores': 1, 'memory_mb': 1000}, 'molpro')
    assert 'gprint,orbitals,orben,civector,ref' in result
    assert 'orbprint,2;' in result
    assert result.count('thrprint=') == 1 and result.count('printci=') == 1
    assert ensure_native_printing(result) == result
    assert ensure_native_printing('{ccsd(t)-f12}\n') == '{ccsd(t)-f12}\n'


def test_legacy_custom_template_and_reader_share_evidence(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot.molpro import Molpro
    import numpy as np
    monkeypatch.chdir(tmp_path)
    Path('molpro').mkdir()
    Path('custom.tpl').write_text('basis=cc-pVDZ\n{{multi;occ,2;closed,0;wf,2,1,0;}}\n{{mrci;}}\n')
    species = SimpleNamespace(atom=['H', 'H'], geom=np.array([[0., 0., 0.], [0., 0., 0.74]]),
                              natom=2, name='h2', charge=0, mult=1)
    calculator = Molpro(species, {'single_point_template': 'custom.tpl'})
    calculator.create_molpro_input(name='h2')
    assert 'gprint,orbitals,orben,civector,ref' in Path('molpro/h2.inp').read_text()
    Path('molpro/h2.out').write_text('PROGRAM * MULTI\nGLOBAL ERROR fehler\n')
    assert calculator.get_molpro_energy('MYENERGY', name='h2') == (0, -1)
    assert json.loads(Path('molpro/h2.orbital_audit.json').read_text())['native_error_lines']


def test_solver_retries_preserve_printing_and_numerical_controls_separately():
    from kinbot.molpro_solver import apply_solver_options
    deck, _ = _vrc_multireference_method('CASPT2-F12(5,4)', 25, 1)
    text = 'gthresh,energy=1.d-10,coeff=1.d-10\n' + deck
    changed = apply_solver_options(text, {'scf': {'energy_convergence': 1e-8}})
    assert native_print_directives() in changed
    assert 'gthresh,energy=1.d-8,coeff=1.d-10' in changed
    assert apply_solver_options(changed, {'scf': {'energy_convergence': 1e-8}}) == changed


def test_sampling_failure_still_writes_native_evidence(tmp_path, monkeypatch):
    from kinbot.rotdpy_molpro import MolproSample
    from types import SimpleNamespace
    import kinbot.rotdpy_molpro as bridge
    monkeypatch.chdir(tmp_path)
    text, inp = real_case()
    Path('sample.inp').write_text(inp)
    def run(*args, **kwargs):
        Path('sample.out').write_text(text + '\nGLOBAL ERROR fehler\n')
        return SimpleNamespace(returncode=255)
    monkeypatch.setattr(bridge.subprocess, 'run', run)
    monkeypatch.setattr(bridge, 'qc_runtime_environment', lambda *a, **k: ({}, {'scratch_directory': str(tmp_path)}))
    monkeypatch.setattr(bridge, 'cleanup_qc_runtime', lambda *a: None)
    with pytest.raises(RuntimeError, match='status 255'):
        MolproSample._run_basis_job(None, SimpleNamespace(name='sample', procs=1))
    assert json.loads(Path('sample.orbital_audit.json').read_text())['native_error_lines']


def test_cc_output_does_not_create_a_cas_report(tmp_path):
    p = tmp_path / 'cc.out'
    p.write_text('PROGRAM * CCSD\nMolpro calculation terminated\n')
    assert write_native_report(p) is None
    assert not p.with_suffix('.orbital_audit.json').exists()


def test_mrci_davidson_variants_and_reference_information_stay_distinct():
    text, inp = real_case()
    text += '''
 PROGRAM * MRCI-F12
 Number of optimized states: 2  Roots: 1 2
 Number of reference states: 2  Roots: 1 2
 Reference symmetry: 1 Doublet
 !CI-F12 STATE 1.1 Energy -437.40555602
 Cluster corrected energies -437.45 (Davidson, fixed reference)
 Cluster corrected energies -437.46 (Davidson, relaxed reference)
 Cluster corrected energies -437.47 (Davidson, rotated reference)
 SETTING KBVRCSEL = -437.46 AU
 Molpro calculation terminated
'''
    report = parse_native_output(text, inp)
    section = report['native_sections'][-1]
    assert section['program'] == 'MRCI-F12'
    assert section['state_energies'][0]['energy_Eh'] == pytest.approx(-437.40555602)
    evidence = '\n'.join(r['text'] for r in section['evidence'])
    assert all(label in evidence for label in ('fixed reference', 'relaxed reference', 'rotated reference'))
    assert 'KBVRCSEL = -437.46' in evidence


@pytest.mark.parametrize('code', [0, 255])
def test_external_dispatch_writes_report_before_classifying_exit(tmp_path, monkeypatch, code):
    import kinbot.anl.dispatch as dispatch
    from types import SimpleNamespace
    text, inp = real_case()
    (tmp_path / 'job.inp').write_text(inp)
    (tmp_path / 'geometry.xyz').write_text('1\nStaging fixture\nH 0 0 0\n')
    scratch = tmp_path / 'scratch'
    scratch.mkdir()
    task = {'id': 'job', 'kind': 'external', 'backend': 'molpro',
            'resources': {'cores': 1, 'memory_mb': 8000},
            'input_name': 'job.inp', 'command': ['molpro', '{input}'],
            'required_outputs': ['job.out']}
    def run(command, **kwargs):
        local = scratch / 'work'
        assert kwargs['cwd'] == local
        assert (local / 'job.inp').read_text() == inp
        for flag, expected in [('-d', scratch), ('-I', scratch), ('-W', local)]:
            assert command[command.index(flag) + 1] == str(expected)
        (local / 'job.out').write_text(text)
        assert not (tmp_path / 'job.out').exists()
        return SimpleNamespace(returncode=code)
    monkeypatch.setattr(dispatch.subprocess, 'run', run)
    monkeypatch.setattr(dispatch, 'qc_runtime_environment',
                        lambda *a, **k: ({}, {'scratch_directory': str(scratch)}))
    monkeypatch.setattr(dispatch, 'cleanup_qc_runtime', lambda *a: None)
    monkeypatch.setattr(dispatch, '_validate_external_outputs', lambda *a: None)
    if code:
        with pytest.raises(dispatch.CalculationUnavailableError, match='255'):
            dispatch._run_external(tmp_path, {'task': task})
    else:
        details = dispatch._run_external(tmp_path, {'task': task})
        assert details['generated_files'] == ['job.orbital_audit.json']
    data = json.loads((tmp_path / 'job.orbital_audit.json').read_text())
    assert data['native_sections'][0]['active']['total'] == 4


def test_cli_reads_native_input_in_preference_to_original(tmp_path):
    import sys
    text, inp = real_case()
    (tmp_path/'a.out').write_text(text)
    (tmp_path/'a.inp').write_text('obsolete input')
    (tmp_path/'a.native.inp').write_text(inp)
    result = subprocess.run([sys.executable, '-m', 'kinbot.molpro_diagnostics', str(tmp_path)],
                            capture_output=True, text=True, check=True)
    assert 'evidence_recorded' in result.stdout
    data = json.loads((tmp_path/'a.orbital_audit.json').read_text())
    assert data['sources']['input']['path'].endswith('a.native.inp')
