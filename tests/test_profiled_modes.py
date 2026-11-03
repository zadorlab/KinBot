"""The complete routed acceptance path must honor native mode frames."""
import inspect
from pathlib import Path

import pytest

from kinbot.profiled_qc import ProfiledQuantumChemistry
from kinbot.qc import QuantumChemistry
from tests.profiled_mode_fixture import stage_profiled_modes, compare_profiled_modes


def test_profiled_public_overrides_accept_the_backend_call_contract():
    # Check all shared overrides, not just the argument that failed remotely.
    # Bind both positional and keyword invocation without calling/submitting QC.
    for name, method in vars(ProfiledQuantumChemistry).items():
        base = getattr(QuantumChemistry, name, None)
        if name.startswith('_') or not callable(method) or not callable(base):
            continue
        parameters = list(inspect.signature(base).parameters.values())
        signature = inspect.signature(method)
        positional, keywords = [], {}
        for p in parameters:
            if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD):
                positional.append(object())
            elif p.kind == p.KEYWORD_ONLY:
                keywords[p.name] = object()
        signature.bind(*positional, **keywords)
        signature.bind(**{p.name: object() for p in parameters
                          if p.kind in (p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY)})


@pytest.mark.parametrize('preset', ['b3lyp-b2plyp-anl', 'uma-b2plyp-anl'])
def test_profiled_acceptance_reads_native_frames_and_replays_without_qc(tmp_path, monkeypatch, preset):
    monkeypatch.chdir(tmp_path)
    fixture = stage_profiled_modes(tmp_path, preset)
    def forbidden(*args, **kwargs):
        pytest.fail('Accepted comparison must not submit or invalidate a calculation')
    monkeypatch.setattr(QuantumChemistry, 'submit_qc', forbidden)
    monkeypatch.setattr(QuantumChemistry, 'invalidate_qc', forbidden)
    before = {p: p.read_bytes() for p in tmp_path.iterdir() if p.is_file()}
    for _ in range(2):
        opt = compare_profiled_modes(fixture)
        assert opt.shigh == 1
        assert opt.mode_comparison['status'] == 'calculated'
        assert opt.mode_comparison['overlap'] == pytest.approx(1.)
        assert opt.selected_job == opt.species.source_job == fixture.higher
        assert opt.qc._level_for(fixture.higher) == 'l2'
        assert {p: p.read_bytes() for p in before} == before


@pytest.mark.parametrize('problem', ['missing', 'wrong_frame', 'incomplete'])
def test_profiled_modes_keep_native_evidence_checks(tmp_path, monkeypatch, problem):
    monkeypatch.chdir(tmp_path)
    fixture = stage_profiled_modes(tmp_path)
    qc = QuantumChemistry(fixture.par)
    _, geom = qc.get_qc_geom(fixture.higher, len(fixture.symbols))
    path = Path(fixture.higher+'.log')
    if problem == 'missing':
        path.unlink()
    elif problem == 'wrong_frame':
        geom[0, 0] += .1
    else:
        path.write_text(path.read_text().replace('Normal termination of Gaussian', 'unfinished Gaussian'))
    with pytest.raises((OSError, ValueError)):
        qc.read_imag_mode(fixture.higher, len(geom), geometry=geom,
                         symbols=fixture.symbols)
