"""Prepare supplemental checks from a stopped, independently reviewed run.

No source result is edited, cancelled, or relabelled. All choices are in the
plan JSON; normal KinBot preparation/execution supplies the actual decks.
"""
import argparse
from copy import deepcopy
import fcntl
import json
from pathlib import Path
import shutil
import tempfile

from kinbot.execution import atomic_json, file_sha256
from kinbot import molpro_reconnaissance as recon, run_cleanup
from kinbot.run_format import ensure_current_run
from run_confirmation import validated_record


def source_snapshot(source, review):
    record = validated_record(source)
    expected = {s['title']: s for s in review['systems']}
    if set(expected) != {s['title'] for s in record['systems']}:
        raise ValueError('Source chains differ from the reviewed export.')
    systems, proof = {}, []
    for system in record['systems']:
        manifest = json.loads(Path(system['manifest']).read_text())
        cfg = json.loads(Path(system['configuration']).read_text())
        old = {j['name']: j for j in expected[system['title']]['jobs']}
        if set(old) != {j['name'] for j in manifest['jobs']}:
            raise ValueError('Source points differ from the reviewed export.')
        jobs = {}
        for job in manifest['jobs']:
            prior = old[job['name']]
            d = Path(job['directory'])
            inp, out = d/(job['name']+'.inp'), d/(job['name']+'.out')
            if job['status'] != prior['status'] or file_sha256(inp) != prior['input_sha256']:
                raise ValueError('Source status/input changed since review: '+str(d))
            digest = file_sha256(out) if out.exists() else None
            if digest != prior['output_sha256']:
                raise ValueError('Source output changed since review: '+str(d))
            if job['status'] == 'complete':
                execution = json.loads((d/'execution.json').read_text())
                if execution.get('status') != 'complete' or execution.get('returncode') != 0:
                    raise ValueError('Completed source has invalid execution: '+str(d))
                if file_sha256(job['wfu']) != prior['wfu_sha256']:
                    raise ValueError('Completed source checkpoint changed: '+str(d))
            jobs[job['anchor']] = job
            proof.append({'chain': system['title'], 'name': job['name'],
                          'status': job['status'], 'input_sha256': prior['input_sha256'],
                          'output_sha256': digest, 'wfu_sha256': prior.get('wfu_sha256')})
        systems[system['title']] = (cfg, jobs)
    return record, systems, proof


def prepare(source, root, plan_path, review_path):
    source, root = Path(source).resolve(), Path(root).resolve()
    if source == root or source in root.parents or root in source.parents:
        raise ValueError('Recovery must use a separate sibling run directory.')
    plan, review = json.loads(Path(plan_path).read_text()), json.loads(Path(review_path).read_text())
    with (source/'.confirmation.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        live = run_cleanup.inspect(source)
        if live['jobs'] or live['processes']:
            raise RuntimeError('Source run is still live; preserve it and wait for completion.')
        record, systems, proof = source_snapshot(source, review)
        provenance = {'source_directory': str(source), 'plan': plan, 'review_sha256': file_sha256(review_path),
                      'source_catalog_sha256': record['catalog_sha256'], 'source_points': proof}
        if root.exists():
            existing = validated_record(root)
            if json.loads((root/'recovery_provenance.json').read_text()) != provenance:
                raise ValueError('Recovery definition changed; use a fresh directory.')
            live = run_cleanup.inspect(root)
            if live['jobs'] or live['processes']:
                raise RuntimeError('Recovery already running; use monitor, not another controller.')
            if existing['status'] not in ('prepared_review_required', 'interrupted'):
                raise ValueError('Recovery already finished or failed; review its results before another retry.')
            for s in existing['systems']:
                for j in json.loads(Path(s['manifest']).read_text())['jobs']:
                    if file_sha256(Path(j['directory'])/(j['name']+'.inp')) != j['input_sha256']:
                        raise ValueError('Recovery input changed.')
            return existing
        prepared = []
        for item in plan['chains']:
            cfg, jobs = systems[item['source_chain']]
            cfg = deepcopy(cfg)
            anchors = {a['name']: a for a in cfg['anchors']}
            cfg['title'] = item['title']
            cfg['anchors'] = [deepcopy(anchors[name]) for name in item['points']]
            for a in cfg['anchors']:
                a['orbital_chain'] = cfg['title']
            cfg['geometry_angstrom'] = cfg['anchors'][0]['geometry_angstrom']
            cfg['probe_definition']['state_model']['correlation_options'].update(item['correlation_options'])
            if 'correlation_checks' in item:
                cfg['correlation_checks'] = deepcopy(item['correlation_checks'])
            if item.get('seed_from_completed'):
                parent = jobs[item['seed_from_completed']]
                if parent['status'] != 'complete':
                    raise ValueError('Recovery seed must be a completed point.')
                seed = {'path': parent['wfu'], 'sha256': parent['wfu_sha256'],
                        'basis': parent['basis'], 'project_orbitals': False}
                # This WFU already contains the selected CAS. Do not apply the
                # old full-valence seed partition or orbital permutation again.
                space = cfg['probe_definition']
                space['orbital_partition'] = {
                    'source_active_electrons': space['active_electrons'],
                    'source_active_orbitals': space['active_orbitals'],
                    'close_lowest': 0, 'discard_highest': 0, 'orbital_swaps': []}
            else:
                seed = deepcopy(next(iter(cfg['initial_records'].values())))
            if file_sha256(seed['path']) != seed['sha256']:
                raise ValueError('Recovery seed hash mismatch.')
            first_name = cfg['bases'][0]['name']+'__'+item['points'][0]
            cfg['initial_records'] = {first_name: seed}
            cfg['confirmation_provenance'].update(recovery=item, source_run=str(source))
            # Validate every candidate before making the destination run.
            with tempfile.TemporaryDirectory() as tmp:
                p = Path(tmp)/'input.json'
                atomic_json(p, cfg)
                data = recon._load(p)
                previous = seed
                for anchor in cfg['anchors']:
                    name = cfg['bases'][0]['name']+'__'+anchor['name']
                    recon.validate_input_blocks(recon._deck(data, cfg['bases'][0]['basis'], anchor, name, previous))
                    previous = name
            prepared.append((cfg, seed))
        root.mkdir(parents=True)
        ensure_current_run(root, create=True)
        atomic_json(root/'recovery_provenance.json', provenance)
        atomic_json(root/'catalog.json', json.loads((source/'catalog.json').read_text()))
        entries = []
        for cfg, seed in prepared:
            dest = root/'seeds'/(cfg['title']+'.wfu')
            dest.parent.mkdir(exist_ok=True)
            shutil.copy2(seed['path'], dest)
            if file_sha256(dest) != seed['sha256']:
                raise ValueError('Copied seed failed verification.')
            seed['path'] = str(dest)
            p = root/'inputs'/(cfg['title']+'.json')
            p.parent.mkdir(exist_ok=True)
            atomic_json(p, cfg)
            manifest = recon.prepare(p, root/cfg['title'])
            m = json.loads(manifest.read_text())
            m['slurm_excluded_nodes'] = record['execution']['excluded_nodes']
            m['startup_retries'] = record['execution']['startup_retries']
            m['solver_options'] = json.loads((source/'catalog.json').read_text())['settings']['solver_options']
            atomic_json(manifest, m)
            entries.append({'title': cfg['title'], 'configuration': str(p),
                            'configuration_sha256': file_sha256(p), 'manifest': str(manifest), 'jobs': len(m['jobs'])})
        result = {'schema': 1, 'status': 'prepared_review_required', 'production_approved': False,
                  'catalog_sha256': file_sha256(root/'catalog.json'), 'systems': entries,
                  'execution': {**record['execution'], 'max_jobs': plan['max_jobs']}}
        atomic_json(root/'confirmation.json', result)
        return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('root', type=Path)
    p.add_argument('--plan', type=Path, default=Path(__file__).with_name('recovery_plan.json'))
    p.add_argument('--review', type=Path, default=Path(__file__).with_name('native_r2_review.json'))
    a = p.parse_args()
    result = prepare(a.source, a.root, a.plan, a.review)
    print(json.dumps({'status': result['status'], 'points': sum(s['jobs'] for s in result['systems']),
                      'chains': len(result['systems']), 'source_preserved': str(a.source), 'recovery': str(a.root)}, indent=2))
