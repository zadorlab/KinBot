"""Installed review interfaces and arbitrary campaign inventories."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from kinbot import review_campaign as review
from kinbot.execution import atomic_json


@pytest.mark.parametrize('module', [*review.MODULES.values(), 'kinbot.review_campaign'])
def test_installed_entry_points_work_outside_checkout(tmp_path, module):
    result = subprocess.run([sys.executable, '-B', '-m', module, '--help'], cwd=tmp_path,
        env=dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1])),
        text=True, capture_output=True, timeout=30)
    assert result.returncode == 0 and 'usage:' in result.stdout


def test_arbitrary_independent_inventory_launches_all_workers_before_polling(tmp_path, monkeypatch):
    records = [dict(name='arbitrary_alpha', kind='orbitals', directory='model-comparisons'),
               dict(name='arbitrary_beta', kind='relaxations', directory='geometry-checks')]
    for row in records:
        (tmp_path/row['directory']).mkdir()
    launched = []
    class Child:
        def __init__(self, command, **kwargs):
            launched.append(command); self.pid=len(launched); self.returncode=None
            assert command[3:5] == ['-m', review.MODULES[records[self.pid-1]['kind']]]
            assert '--queue-all' in command and '--allow-node' in command
        def poll(self):
            assert len(launched) == len(records)
            self.returncode = 0
            return 0
    monkeypatch.setattr(review.subprocess, 'Popen', Child)
    monkeypatch.setattr(review.requested_orbital_recon, 'evidence', lambda p:'Recorded synthetic evidence')
    assert review.supervise(tmp_path, records, queue_all=True, allow_nodes=['available']) == {
        'arbitrary_alpha':0, 'arbitrary_beta':0}
    journal = json.loads((tmp_path/'queue_all_finished.json').read_text())
    assert journal['production_model_accepted'] is False
    assert (tmp_path/'model-comparisons/orbital_evidence.txt').exists()


@pytest.mark.parametrize('fault', ['outside','duplicate_path','duplicate_name','invalid_kind','invalid_label','empty'])
def test_invalid_inventory_is_rejected_before_launch(tmp_path, fault):
    folder = tmp_path/'case'; folder.mkdir()
    row = dict(name='one', kind='orbitals', directory='case')
    records = [row]
    if fault == 'outside': row['directory'] = '..'
    elif fault == 'duplicate_path': records.append(dict(row,name='two'))
    elif fault == 'duplicate_name': records.append(dict(row,directory='.'))
    elif fault == 'invalid_kind': row['kind'] = 'execute_arbitrary_code'
    elif fault == 'invalid_label': row['name'] = '../escape'
    else: records = []
    with pytest.raises(ValueError):
        review.inventory(tmp_path, records)


def test_general_modules_do_not_import_example_implementations():
    import ast
    from kinbot import requested_orbital_recon, vrc_relaxation_recon
    for module in (review, requested_orbital_recon, vrc_relaxation_recon):
        tree = ast.parse(Path(module.__file__).read_text())
        imports = [n.module or '' for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)]
        imports += [a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names]
        assert not any('examples' in n or 'progress_2026' in n for n in imports)
