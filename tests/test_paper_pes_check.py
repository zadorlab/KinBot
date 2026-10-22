import json

from ase import Atoms
from ase.db import connect

from kinbot import paper_pes_check as check


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def fixture(tmp_path):
    case = tmp_path/'ch3ss'
    put(case/'input.json', {'smiles': 'CS[S]', 'mult': 2})
    put(case/'pes_exploration.json', dict(phase='exploring', initial_worker='400',
        token='current', epoch=1, workers={'200': {'status': 'exploring'}}))
    put(case/'pes_l2_selection.json', dict(status='selected', token='current', epoch=1,
        workers={'200': {'selected': ['200_htransfer'], 'omitted': {}}}))
    labels = {'200': '[CH2]SS', '300': 'SC[S]', '400': 'CS[S]',
              '11': 'C=S', '12': '[SH]', '21': '[CH]=S', '22': 'S',
              '31': 'S=CS', '32': '[H]', '41': '[CH3]', '42': '[S][S]',
              '500': '[CH]=[SH]S'}
    for worker in ('200', '300', '400'):
        w = case/worker
        put(w/(worker+'.json'), dict(smiles=labels[worker], mult=2))
        put(w/'accepted_l2_network.json', {'network': {'route_to_smiles': labels}})
        (w/'kinbot_monitor.out').write_text('')
    (case/'200/kinbot_monitor.out').write_text('-1 23 200_htransfer 400\n')
    (case/'400/kinbot_monitor.out').write_text('-1 7 400_hom_sci_1_2 41 42\n')
    return case


def index(report, name):
    return {r['paper']: r for r in report[name]}


def test_connectivity_uses_h_attachment_not_formula_or_bond_order():
    assert len({check.endpoint(v[1]) for k, v in check.NODES.items() if k != 'W1'}) == 7
    assert check.connectivity('SC[S]') == check.connectivity('[H]SC([H])([H])[S]')
    assert check.connectivity('SC[S]') != check.connectivity('[CH2]SS')
    assert check.connectivity('[S][S]') == check.connectivity('S=S')


def test_read_only_reverse_paths_and_extra_discoveries(tmp_path):
    case = fixture(tmp_path)
    (case/'300/kinbot_monitor.out').write_text('-1 23 300_other 500\n')
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in case.rglob('*') if p.is_file()}
    result = check.snapshot(tmp_path)
    assert not result['errors']
    nodes, edges = index(result, 'paper_nodes'), index(result, 'paper_connections')
    assert nodes['W4']['status'] == 'accepted_connectivity'
    assert edges['TS4']['status'] == 'accepted_endpoint_and_path_type_match'
    assert edges['TS4']['evidence'][0]['selection'] == 'retained'
    assert edges['scission']['status'] == 'accepted_endpoint_and_path_type_match'
    assert [r['route'] for r in result['additional_accepted_species']] == ['500']
    assert [r['name'] for r in result['additional_accepted_connections']] == ['300_other']
    assert before == {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in case.rglob('*') if p.is_file()}
    assert 'read-only' in check.render(result)


def test_separated_fragments_do_not_prove_complex_or_tight_ts(tmp_path):
    case = fixture(tmp_path)
    (case/'200/kinbot_monitor.out').write_text('-1 7 200_hom_sci_1_2 11 12\n')
    result = check.snapshot(case)
    nodes, edges = index(result, 'paper_nodes'), index(result, 'paper_connections')
    assert nodes['P1']['status'] == 'accepted_connectivity'
    assert nodes['W1']['status'] == 'not_found_in_accepted_complex_records'
    assert edges['TS2']['status'] == 'accepted_endpoints; different_path_type'
    assert not edges['TS2']['W1_complex_recorded']
    assert edges['entrance']['status'] == 'not_found'
    assert 'net endpoint match does not establish' in check.render(result)


def test_complex_record_requires_current_acceptance_and_matching_fragments(tmp_path):
    case = fixture(tmp_path)
    w = case/'200'
    path = w/'accepted_l2_network.json'
    saved = json.loads(path.read_text())
    saved['network']['reactions'] = [dict(name='200_dissociation', reactant='200',
        products=['12', '11'], product_complex='vdw')]
    put(path, saved)
    (w/'kinbot_monitor.out').write_text('-1 23 200_dissociation 11 12\n')
    result = check.snapshot(case)
    assert index(result, 'paper_nodes')['W1']['routing_keys'] == ['vdw']
    assert 'orientation_unverified' in index(result, 'paper_nodes')['W1']['status']
    (w/'kinbot_monitor.out').write_text('-999 23 200_dissociation\n')
    result = check.snapshot(case)
    assert index(result, 'paper_nodes')['W1']['status'] == 'not_found_in_accepted_complex_records'
    assert index(result, 'paper_connections')['TS2']['status'] == 'rejected_candidate_only'


def test_pending_endpoint_hints_are_not_discoveries(tmp_path):
    case = fixture(tmp_path)
    (case/'300/kinbot_monitor.out').write_text('1 16 300_dissociation\n')
    (case/'300/kinbot.log').write_text('Based on the end of IRC, reaction 300_dissociation leads to products [31, 32]\n')
    result = check.snapshot(case)
    assert index(result, 'paper_nodes')['P3']['status'] == 'candidate_endpoint_in_progress'
    assert index(result, 'paper_connections')['TS5']['status'] == 'candidate_in_progress'
    assert not result['additional_accepted_connections']


def test_stale_selection_and_unknown_products_are_visible(tmp_path):
    case = fixture(tmp_path)
    sel = case/'pes_l2_selection.json'
    data = json.loads(sel.read_text())
    data['token'] = 'previous'
    put(sel, data)
    (case/'300/kinbot_monitor.out').write_text('-1 23 300_unknown 999\n')
    result = check.snapshot(case)
    assert not result['selection_current']
    assert index(result, 'paper_connections')['TS4']['evidence'][0]['selection'] == 'not_yet_selected'
    assert [r['name'] for r in result['unresolved_endpoints']] == ['300_unknown']
    assert not result['additional_accepted_connections']


def test_missing_current_monitor_does_not_resurrect_snapshot_success(tmp_path):
    case = fixture(tmp_path)
    path = case/'200/kinbot_monitor.out'
    path.unlink()
    result = check.snapshot(case)
    assert len(result['errors']) == 1
    assert result['errors'][0]['file'] == str(path)
    assert index(result, 'paper_connections')['TS4']['status'] == 'not_found_in_current_records'


def test_duplicate_searches_remain_evidence_not_unique_saddles(tmp_path):
    case = fixture(tmp_path)
    (case/'400/kinbot_monitor.out').write_text('-1 23 400_reverse 200\n-1 23 400_alias 200\n')
    result = check.snapshot(case)
    edge = index(result, 'paper_connections')['TS4']
    assert len(edge['evidence']) == 3
    assert result['search_counts']['accepted'] == 3
    assert not result['additional_accepted_connections']


def test_native_product_cache_is_read_only_and_latest_failure_wins(tmp_path):
    case = fixture(tmp_path)
    worker = case/'300'
    (worker/'kinbot_monitor.out').write_text('-1 23 300_new_products 999 31\n')
    database = worker/'kinbot.db'
    db = connect(database)
    db.write(Atoms('H', positions=[[0., 0., 0.]]), name='999_well_high', data={'status': 'normal'})
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in case.rglob('*') if p.is_file()}
    report = check.snapshot(case)
    assert not report['errors']
    assert index(report, 'paper_connections')['TS5']['status'] == 'accepted_endpoint_and_path_type_match'
    assert report['cache_identity_evidence']['999']['source_row_id'] == 1
    assert before == {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in case.rglob('*') if p.is_file()}
    db.write(Atoms('H', positions=[[0., 0., 0.]]), name='999_well_high', data={'status': 'error'})
    report = check.snapshot(case)
    assert index(report, 'paper_connections')['TS5']['status'] == 'not_found_in_current_records'
    assert report['unresolved_endpoints'][0]['name'] == '300_new_products'


def test_additional_unfinished_wells_and_candidates_are_separate(tmp_path):
    case = fixture(tmp_path)
    put(case/'500/500.json', {'smiles': '[CH]=[SH]S', 'mult': 2})
    (case/'500/kinbot_monitor.out').write_text('1 16 500_candidate\n')
    (case/'500/kinbot.log').write_text('reaction 500_candidate leads to products [400]\n')
    report = check.snapshot(case)
    assert [r['route'] for r in report['additional_pending_wells']] == ['500']
    assert [r['name'] for r in report['additional_in_progress_connections']] == ['500_candidate']
    assert not report['additional_accepted_species']
