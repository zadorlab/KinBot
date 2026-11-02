"""Opt-in complete captured handoff through the real no-kinbot CLI twice.

Set KINBOT_TEST_MESS_HANDOFF to the native captured archive. Ordinary CI uses
compact native excerpts in test_mess_mirror_capture and population_identity.
"""
import hashlib
import json
import os
from pathlib import Path
import sys
import tarfile

import pytest

from kinbot import pes
from kinbot.mess_networks import connected_models, _HEADER, split_model


SCRIPT = Path(__file__).resolve().parents[1]/'examples/anl/ch3s_anl1_f12/scouting_r8/postprocess_blodgett.sh'


def prepare(root, home, monkeypatch, jobs=(), processes=()):
    from kinbot import run_cleanup
    monkeypatch.setattr(Path, 'home', lambda: home)
    monkeypatch.setattr(sys, 'argv', ['prepare', str(root)])
    monkeypatch.setattr(run_cleanup, '_jobs', lambda _: jobs)
    monkeypatch.setattr(run_cleanup, '_processes', lambda _: processes)
    code = SCRIPT.read_text().split("<<'PY'\n", 1)[1].split('\nPY\n', 1)[0]
    exec(compile(code, str(SCRIPT), 'exec'), {})


@pytest.mark.parametrize('busy', ['jobs', 'processes', 'incomplete'])
def test_postprocessing_handoff_checks_both_cases_before_any_writes(tmp_path, monkeypatch, busy):
    root=tmp_path/'run'
    originals={}
    for case in ('ch3s','ch3ss'):
        d=root/case; (d/'123').mkdir(parents=True)
        raw=dict(me=1,run_me=0,rotdpy_mess_mode='off')
        (d/'input.json').write_text(json.dumps(raw))
        (d/'chemids').write_text('123\n')
        status='failed' if busy=='incomplete' and case=='ch3ss' else 'complete'
        (d/'123/kinbot_status.json').write_text(json.dumps(dict(status=status)))
    originals={p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
    with pytest.raises(AssertionError):
        prepare(root,tmp_path,monkeypatch,jobs=[{}] if busy=='jobs' else (),
                processes=[{}] if busy=='processes' else ())
    assert {p:p.read_bytes() for p in root.rglob('*') if p.is_file()} == originals
    assert not (tmp_path/'kinbot-evidence').exists()


def test_complete_captured_scout_handoff_without_qc(tmp_path, monkeypatch):
    archive = os.environ.get('KINBOT_TEST_MESS_HANDOFF')
    if not archive:
        pytest.skip('Set KINBOT_TEST_MESS_HANDOFF to the complete native capture')
    with tarfile.open(archive) as tar:
        tar.extractall(tmp_path, filter='data')
    manifest = json.loads((tmp_path/'manifest.json').read_text())
    assert manifest['purpose'] == 'mess_handoff_replay' and manifest['status'] == 'captured'
    for item in manifest['files']:
        source = tmp_path/item['path']
        assert hashlib.sha256(source.read_bytes()).hexdigest() == item['sha256']
    root = tmp_path/'run'
    # The earlier native QC audit supplies rxn20's matching frequency output;
    # the handoff archive itself intentionally captured no large native logs.
    evidence = json.loads((Path(__file__).parent/'reference/mess_symmetry_review.json').read_text())
    reaction = evidence['states']['rxn_20']['input_header'].split('!')[1].strip()
    workers = [d for d in (root/'ch3ss').iterdir() if d.is_dir() and reaction.startswith(d.name+'_')]
    assert len(workers) == 1
    native = workers[0]/(reaction+'_high.log')
    native_bytes = (Path(__file__).parent/'reference/mirror_even_transition.log').read_bytes()
    native.write_bytes(native_bytes)
    protected = {p: p.read_bytes() for p in root.rglob('*') if p.is_file() and (
        p.name in ('input.json', 'kinbot.db', '.kinbot_search_state.json', 'l1_screening.json')
        or p.name == p.parent.name+'.json' or p.suffix == '.mess')}

    def forbidden(*args, **kwargs):
        pytest.fail('Complete postprocessing must not launch a worker or QC')
    monkeypatch.setattr(pes, 'submit_job', forbidden)
    monkeypatch.setattr(pes, 'check_status', forbidden)
    monkeypatch.setattr(pes, 'get_l3energy', forbidden)
    monkeypatch.setattr(pes, 't1_analysis', forbidden)
    monkeypatch.setattr(pes, '_write_rotdpy_slurm_template', forbidden)
    monkeypatch.setattr('kinbot.qc.QuantumChemistry.submit_qc', forbidden)
    monkeypatch.setattr(pes.time, 'sleep', lambda _: None)
    prepare(root, tmp_path, monkeypatch)
    all_models = {}
    for source in sorted(root.glob('*/input.json')):
        par = json.loads(source.read_text())
        assert not par.get('L3_calc') and not par.get('rotdpy_run')
        assert not par.get('vrc_tst_scan') and not par.get('vrc_tst_noscan')
        replay_input = source.with_name('mess_scout.json')
        assert json.loads(replay_input.read_text()) == dict(par, correct_submerged=1, run_me=0, run_messpf=0, queuing='local')
        monkeypatch.chdir(source.parent)
        monkeypatch.setattr(sys, 'argv', ['pes', replay_input.name, 'no-kinbot'])
        previous = None
        for _ in range(2):
            pes.main()
            assert json.loads(Path('pes_status.json').read_text())['status'] == 'complete'
            networks = json.loads(Path('me/mess_networks.json').read_text())
            models = {n['stem']: Path('me', n['stem']+'.inp').read_text() for n in networks}
            assert models and all(n['status'] == 'ready' for n in networks)
            for model in models.values():
                assert len(connected_models(model)) == 1
                blocks = split_model(model)[1]
                assert len({_HEADER.search(b)[2].split()[0] for b in blocks}) == len(blocks)
                assert '{' not in model and 'Core Rotd' not in model
            if source.parent.name == 'ch3ss':
                text = ''.join(models.values())
                assert text.count('kinbot_shared_mirror_surface_weight 0.5') == 2
                import re
                named = {_HEADER.search(b)[2].split()[0]: b for b in split_model(text)[1]}
                for label in ('w_10', 'rxn_29'):
                    assert float(re.search(r'SymmetryFactor\s+(\S+)', named[label])[1]) == .5
                # Captured product labels retained an open-chain input SMILES
                # after optimization made a cyclic molecule. Rebuild display
                # labels from the actual selected geometry, never the reverse.
                for label in ('pr_6', 'pr_18'):
                    assert '1' in _HEADER.search(named[label])[2].split('!')[-1]
                assert '1' not in _HEADER.search(named['pr_11'])[2].split('!')[-1]
                assert native.read_bytes() == native_bytes
            assert all(p.read_bytes() == value for p, value in protected.items())
            if previous is not None:
                assert models == previous
            previous = models
        all_models[source.parent.name] = models
    # Repeat the preparation block too: archives previous outputs, no new QC.
    prepare(root, tmp_path, monkeypatch)
    assert all(p.read_bytes() == value for p, value in protected.items())
    assert len(all_models) == 2
    text = ''.join(model for models in all_models.values() for model in models.values())
    assert text.count('energies[kcal/mol] [41.67] raised to connected-well ground 42.34') == 1
