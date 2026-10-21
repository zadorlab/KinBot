"""Polling completed TS clusters must not rewrite unchanged calculation files."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace as NS
from unittest.mock import patch

from ase.build import molecule

from kinbot import bimol
from kinbot.stationary_pt import StationaryPoint


def prepared_cluster(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = molecule('CH4')
    point = StationaryPoint('test', 0, 1, atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    point.characterize()
    point.energy, point.zpe = -3., .01
    reactions = []
    for i in range(2):
        directory = Path(f'bimol/search/{i}')
        directory.mkdir(parents=True)
        (directory / 'result.json').write_text(json.dumps({'id': str(i)}))
        reactions.append(NS(bimol_search=True, family_name='abstraction',
            instance_name=f'test_{i}', instance=[0, 1, 2], ts=copy.deepcopy(point),
            stereopath_id='same-site', search_parent='parent',
            search_candidate={'directory': str(directory)},
            l0_candidate_id=str(i), l0_ts_id=str(i),
            l0_validation={'endpoints': [
                {'state': role, 'configured_identity': {'status': 'assigned', 'id': role}}
                for role in ('reactant', 'product')]}))
    reactions[1].ts.geom[1, 0] += .01
    for key in ('unique_candidates', 'validated_l0'):
        directory = Path('bimol') / key
        directory.mkdir()
        (directory / 'index.json').write_text(json.dumps({key: [{'id': str(i)} for i in range(2)]}))
        for i in range(2):
            (directory / str(i)).mkdir()
            (directory / str(i) / 'metadata.json').write_text(json.dumps({'id': str(i)}))
    return NS(reac_obj=reactions, reac_ts_done=[6, 6]), {
        'multi_conf_tst': 0, 'bimol_energy_tolerance': .25, 'bimol_ts_rmsd': .1}


def test_cluster_batches_index_writes_and_skips_unchanged_polls(tmp_path, monkeypatch):
    species, par = prepared_cluster(tmp_path, monkeypatch)
    with patch('kinbot.bimol._json', wraps=bimol._json) as write:
        bimol.release_final_ts_clusters(species, par)
    assert sum(Path(c.args[0]).name == 'index.json' for c in write.call_args_list) == 2
    species.reac_ts_done[0] = 2
    with patch('kinbot.bimol._json', wraps=bimol._json) as write:
        bimol.release_final_ts_clusters(species, par)
    assert sum(Path(c.args[0]).name == 'index.json' for c in write.call_args_list) == 2
    before = {p: p.read_bytes() for p in Path('bimol').rglob('*.json')}
    with patch('kinbot.bimol._json') as write, patch.object(Path, 'read_text') as read:
        for _ in range(20):
            bimol.release_final_ts_clusters(species, par)
    write.assert_not_called()
    read.assert_not_called()
    assert before == {p: p.read_bytes() for p in before}
    species.reac_ts_done[0] = -1
    bimol.release_final_ts_clusters(species, par)
    for key in ('unique_candidates', 'validated_l0'):
        records = json.loads((Path('bimol') / key / 'index.json').read_text())[key]
        for record in records:
            archive = record['final_conformer_archive']
            assert archive[0]['final_irc_workflow_status'] == -1
            assert archive == record['lifecycle']['final_conformer_archive']


def test_repeated_status_and_unknown_id_do_not_write_but_real_change_does(tmp_path, monkeypatch):
    species, _ = prepared_cluster(tmp_path, monkeypatch)
    reaction = species.reac_obj[0]
    bimol.update_bimol_status(reaction, final_status='final_irc_valid')
    with patch('kinbot.bimol._json') as write:
        bimol.update_bimol_status(reaction, final_status='final_irc_valid')
        bimol._update_index('bimol/validated_l0/index.json', 'validated_l0', 'missing', {'status': 1})
    write.assert_not_called()
    bimol.update_bimol_status(reaction, final_status='final_exported')
    saved = json.loads(Path('bimol/search/0/result.json').read_text())
    assert saved['final_status'] == saved['lifecycle']['final_status'] == 'final_exported'


def test_failed_representative_still_releases_an_alternative(tmp_path, monkeypatch):
    species, par = prepared_cluster(tmp_path, monkeypatch)
    bimol.release_final_ts_clusters(species, par)
    species.reac_ts_done[0] = -999
    bimol.release_final_ts_clusters(species, par)
    assert species.reac_ts_done == [-999, 7]
    assert species.reac_obj[1].final_irc_attempted
    data = json.loads(Path('bimol/validated_l0/index.json').read_text())['validated_l0'][1]
    assert data['final_status'] == 'retrying_final_irc'
    assert data['final_conformer_archive'][0]['final_irc_workflow_status'] == -999
