"""Render actual Gaussian input, beyond checking intermediate profile dictionaries."""
import ast
from io import StringIO
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms

from kinbot.gaussian_input import route_options, merge_gaussian_kwargs
from kinbot.ase_modules.io.gaussian import write_gaussian_in
from kinbot.ase_modules.calculators.gaussian import fixed_frequency_parameters
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry


def test_parenthesized_merge_and_direct_writer_are_equivalent():
    merged = merge_gaussian_kwargs({'scf': 'xqc,conver=8'}, {'SCF': '(XQC,Conver=10)'})
    assert merged == {'scf': 'XQC,Conver=10'}
    for kwargs in (merged, {'scf': '(XQC,Conver=10)'}):
        stream = StringIO()
        write_gaussian_in(stream, Atoms('H'), method='HF', basis='STO-3G', mult=2, **kwargs)
        assert 'scf(XQC,Conver=10)' in stream.getvalue()


def test_nested_named_options_and_scalar_overrides():
    assert route_options('(ReadAtoms(C,H),MaxCycle=50)') == ['ReadAtoms(C,H)', 'MaxCycle=50']
    assert merge_gaussian_kwargs({'method': 'HF', 'mem': '1GB', 'scf': 'xqc'},
                                {'method': 'B2PLYP', 'mem': '2GB', 'scf': ''}) == {
                                    'method': 'B2PLYP', 'mem': '2GB', 'scf': ''}


@pytest.mark.parametrize('value', ['xqc,(xqc,conver=10)', '(xqc,conver=10', 'xqc,,tight', 'xqc)'])
def test_invalid_option_groups_fail_before_submission(value):
    with pytest.raises(ValueError, match='Gaussian option list'):
        route_options(value)


@pytest.mark.parametrize('method', ['ANL0', 'ANL0-F12', 'ANL1', 'ANL1-F12'])
@pytest.mark.parametrize('flags', [{}, {'ts': 1, 'step': 1, 'max_step': 1},
    {'irc': 'Forward'}, {'hir': 1}, {'hir': 1, 'rigid': 1}, {'vts': 1}, {'scan': 1}])
def test_all_anl_l2_routes_render_without_duplicate_or_nested_scf(method, flags, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps({'composite_method': method,
                                            'barrier_threshold': 1000}))
    par = Parameters('input.json').par
    qc = QuantumChemistry(par).l2
    kwargs = qc.merge_kwargs(qc.get_qc_arguments('probe', 2, 0, 25, high_level=1, **flags))
    kwargs.pop('label')
    # Sella's force and separate native Hessian calls both reach this writer.
    for options, properties in [(kwargs, ['energy', 'forces']),
                                (fixed_frequency_parameters(kwargs), ['energy'])]:
        stream = StringIO()
        write_gaussian_in(stream, Atoms('CH3S', positions=np.zeros((5, 3))),
                          properties=properties, **options)
        route = stream.getvalue().split('\n\n')[0]
        assert route.count('scf(') == 1
        assert 'scf(xqc,conver=10)' in route
        assert 'scf(xqc,(' not in route
        assert 'EmpiricalDispersion(GD3BJ)' in route
        assert 'integral(UltraFine)' in route


def test_production_high_level_worker_renders_real_route(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps({'composite_method': 'ANL1-F12', 'barrier_threshold': 1000}))
    qc = QuantumChemistry(Parameters('input.json').par)
    species = SimpleNamespace(chemid=123, name='test', mult=2, charge=0, nel=25,
        natom=5, wellorts=0, atom=['C', 'H', 'H', 'H', 'S'], geom=np.zeros((5, 3)))
    monkeypatch.setattr(qc.l2, 'submit_qc', lambda *args, **kw: 1)
    qc.qc_opt(species, species.geom, high_level=1)
    tree = ast.parse(Path('123_well_high.py').read_text())
    kwargs = next(ast.literal_eval(node.value) for node in tree.body if isinstance(node, ast.Assign)
                  and any(isinstance(t, ast.Name) and t.id == 'kwargs' for t in node.targets))
    kwargs.pop('label')
    stream = StringIO()
    write_gaussian_in(stream, Atoms('CH3S', positions=np.zeros((5, 3))), **kwargs)
    assert 'scf(xqc,conver=10)' in stream.getvalue()
    assert 'B2PLYP/cc-pVTZ' in stream.getvalue()
