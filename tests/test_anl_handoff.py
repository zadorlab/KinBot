"""Completed ANL/CBH records attach to exact reconstructed network keys."""

import hashlib
import json
from types import SimpleNamespace

import pytest

from kinbot.anl.cbh import zero_k_from_composite_record
from kinbot.anl.handoff import apply_handoff


def _composite(path, *, charge=0, multiplicity=1):
    components = {
        'electronic': {
            'key': 'electronic', 'value_hartree': -40.1,
            'review_required': False, 'source_sha256': '1' * 64,
        },
        'zpe': {
            'key': 'zpe', 'value_hartree': .1,
            'review_required': False, 'source_sha256': '2' * 64,
        },
    }
    path.write_text(json.dumps({
        'schema': 1, 'status': 'complete', 'recipe': 'ANL0-F12',
        'charge': charge, 'multiplicity': multiplicity,
        'components': components,
        'electronic_terms': [{'component': 'electronic', 'coefficient': 1}],
        'zero_point_terms': [{'component': 'zpe', 'coefficient': 1}],
        'electronic_hartree': -40.1, 'zero_point_hartree': .1,
        'zero_k_hartree': -40.0,
    }))


def _formation(path, final):
    path.write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'formation': {
            'target_smiles': final.smiles, 'rung': 0,
            'method': final.method, 'reaction_energy_0k_kj_mol': 1.0,
            'formation_0k_kj_mol': -2.0, 'references': {},
            'atct_version': 'fixture', 'atct_source_sha256': '3' * 64,
            'energy_sources': {final.smiles: final.source},
        },
    }))


def test_handoff_attaches_exact_species_and_records_hashes(tmp_path, monkeypatch):
    composite = tmp_path / 'methane.json'
    _composite(composite)
    final = zero_k_from_composite_record('C', composite)
    formation = tmp_path / 'methane-cbh.json'
    _formation(formation, final)
    root = SimpleNamespace(
        name='methane', chemid=123, smiles='C', charge=0, mult=1,
        reac_obj=[], reac_ts_done=[])
    monkeypatch.setattr(
        'kinbot.anl.handoff.routing_name', lambda species: str(species.chemid))
    manifest = tmp_path / 'handoff.json'
    manifest.write_text(json.dumps({
        'schema': 1, 'mode': 'cbh-anl',
        'species': {'123': {
            'composite': composite.name, 'formation': formation.name}},
        'transition_states': {},
    }))
    report = tmp_path / 'applied.json'
    result = apply_handoff(root, manifest, report_file=report)
    assert root.final_zero_k_energy.hartree == pytest.approx(-40.)
    assert root.formation_enthalpy_0k.formation_0k_kj_mol == pytest.approx(-2.)
    assert result['species']['123']['composite_sha256'] == hashlib.sha256(
        composite.read_bytes()).hexdigest()
    assert json.loads(report.read_text()) == result


def test_handoff_fails_closed_on_missing_network_species(tmp_path, monkeypatch):
    root = SimpleNamespace(
        name='root', chemid=1, smiles='C', charge=0, mult=1,
        reac_obj=[], reac_ts_done=[])
    monkeypatch.setattr(
        'kinbot.anl.handoff.routing_name', lambda species: str(species.chemid))
    manifest = tmp_path / 'handoff.json'
    manifest.write_text(json.dumps({
        'schema': 1, 'mode': 'anl', 'species': {},
        'transition_states': {}}))
    with pytest.raises(ValueError, match='complete network'):
        apply_handoff(root, manifest, report_file=tmp_path / 'report.json')


def test_handoff_can_attach_vpt2_to_transition_state(tmp_path, monkeypatch):
    composite = tmp_path / 'ts.json'
    _composite(composite)
    ts = SimpleNamespace(
        name='ts', chemid=3, smiles='C', charge=0, mult=1)
    product = SimpleNamespace(
        name='product', chemid=2, smiles='C', charge=0, mult=1)
    reaction = SimpleNamespace(
        instance_name='reaction', do_vdW=False, ts=ts,
        prod_opt=[SimpleNamespace(species=product)])
    root = SimpleNamespace(
        name='root', chemid=1, smiles='C', charge=0, mult=1,
        reac_obj=[reaction], reac_ts_done=[-1])
    monkeypatch.setattr(
        'kinbot.anl.handoff.routing_name', lambda species: str(species.chemid))
    attached = {}
    monkeypatch.setattr(
        'kinbot.anl.handoff.attach_task_vpt2_frequencies',
        lambda point, run, task_id, **options: (
            setattr(point, 'anl_thermochemistry_frequency_source', {
                'run': str(run), 'task_id': task_id, **options}),
            attached.setdefault(point.name, True)))
    manifest = tmp_path / 'handoff.json'
    manifest.write_text(json.dumps({
        'schema': 1, 'mode': 'anl',
        'species': {
            '1': {'composite': composite.name},
            '2': {'composite': composite.name},
        },
        'transition_states': {'reaction': {
            'composite': composite.name,
            'vpt2': {'run': 'ts-vpt2', 'task_id': 'gaussian_vpt2'},
        }},
    }))
    result = apply_handoff(
        root, manifest, report_file=tmp_path / 'report.json')
    assert attached == {'ts': True}
    assert result['transition_states']['reaction']['vpt2']['task_id'] == (
        'gaussian_vpt2')
