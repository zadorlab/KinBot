"""Stage reverse abstraction searches from retained, native PES product pairs.

These independent discovery workers do not claim to be extra PES wells. Their
validated TSs and both IRC endpoints are exported for entrance-complex matching
before a new combined MESS model is authorized.
"""
import argparse
from copy import deepcopy
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import sqlite3
import subprocess
import sys
import time

import numpy as np

from kinbot import constants, l1_screening
from kinbot.execution import atomic_json, file_sha256
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.species_routing import routing_name


def _database(worker):
    from ase.db import connect
    path = (worker/'kinbot.db').resolve(strict=True)
    db = connect(str(path), use_lock_file=False)
    db.connection = sqlite3.connect(path.as_uri()+'?mode=ro', uri=True)
    db.change_count = 0
    return db


def _fragment(worker, name, par):
    point = l1_screening.read_points(worker, [name+'_well'],
                                    l1_screening.theory(par), geometry=True)[0]
    db = _database(worker)
    try:
        row = next(db.select(name='stereochemistry/'+name+'_well', sort='-id', limit=1), None)
        if row is None:
            raise ValueError(f'{worker}/{name}: no native fragment electronic-state reference.')
        ref = dict(row.data)
    finally:
        db.connection.close()
    raw = ref['raw_input']
    fragment = StationaryPoint(name, raw['charge'], raw['multiplicity'],
        atom=point['symbols'], geom=np.array(point['positions']))
    fragment.characterize()
    from kinbot.electronic_states import apply_electronic_state_model
    apply_electronic_state_model(fragment, par)
    if routing_name(fragment) != name:
        raise ValueError(f'{worker}/{name}: optimized fragment differs from its requested state.')
    structure = [value for symbol, position in zip(raw['atoms'], raw['geometry_angstrom'])
                 for value in [symbol, *position]]
    return fragment, point, structure


def plan(case):
    """Read-only, deduplicated inventory; neither paper labels nor run names."""
    case = Path(case).resolve(strict=True)
    ref = l1_screening.reference(case)
    entries = {}
    for worker_name in (case/'chemids').read_text().splitlines():
        worker = case/worker_name
        summary = worker/('summary_'+worker_name+'.out')
        if not summary.is_file():
            raise ValueError(f'{summary}: finish or capture this worker before staging reverse discovery.')
        par = Parameters(worker/(worker_name+'.json')).par
        for line in summary.read_text().splitlines():
            fields = line.split()
            if not fields or fields[0] != 'SUCCESS':
                continue
            products = fields[3:-2] if fields[-1].startswith('vdW') else fields[3:]
            if len(products) != 2:
                continue
            products = sorted(products)
            key = tuple([*products, par['mult']])
            if key in entries:
                entries[key]['sources'].append(dict(worker=worker_name, reaction=fields[2]))
                continue
            frags, points, structures = zip(*[_fragment(worker, name, par) for name in products])
            if sum(f.charge for f in frags) != par['charge']:
                raise ValueError('Reverse fragments do not conserve the parent charge.')
            spin = par['mult']-1
            if spin not in range(abs(frags[0].mult-frags[1].mult),
                                 frags[0].mult+frags[1].mult-1, 2):
                raise ValueError('Reverse fragments cannot couple to the source PES multiplicity.')
            # Enumerate both directions using the real bimolecular finder.
            from kinbot.reaction_finder_bimol import ReactionFinderBimol
            combined = StationaryPoint('pair', [f.charge for f in frags], par['mult'],
                structure=deepcopy(list(structures)), fragA=frags[0], fragB=frags[1])
            combined.characterize()
            combined.chemid = '_'.join(products)
            combined.name = combined.chemid
            settings = deepcopy(par)
            settings.update(bimol=1, bimol_families=['abstraction'], skip_families=[],
                charge=[int(f.charge) for f in frags],
                mult=[int(f.mult) for f in frags]+[int(par['mult'])],
                structure=list(structures), smiles='', stereo_reference=None,
                pes=0, me=0, run_me=0, L3_calc=0, rotdpy_run=0, composite_method='',
                vrc_tst_scan={}, vrc_tst_noscan={}, anl_handoff_manifest='',
                thermochemistry_model_output='', bimol_result_only=True,
                homolytic_saddle_search=False, profiled_l1_reuse_only=0,
                pes_phase_policy='legacy', pes_barrier_reference=None,
                bimol_pes_reference=ref, l1_pes_energy_ceiling=par.get('l1_pes_energy_ceiling'),
                input_file='', do_clean=False, selected_reactions={},
                skip_reactions=[], reaction_search=1)
            selection = l1_screening.evaluate(ref,l1_screening.theory(par),list(points),
                par.get('l1_pes_energy_ceiling'),charge=par['charge'])
            if selection['decision'] != 'retain':
                raise ValueError(f'{products}: source asymptote is not retained by the current L1 policy.')
            ceiling = par.get('l1_pes_energy_ceiling')
            if ceiling is None:
                raise ValueError('Reverse PES staging needs an explicit l1_pes_energy_ceiling; '
                                 'no new energy limit is assumed from an arbitrary source worker.')
            settings.update(barrier_threshold=max(0.,ceiling-selection['relative_kcal_mol']),
                            barrier_threshold_L2=None)
            ReactionFinderBimol(combined,settings,None).find_reactions()
            digest = hashlib.sha256(json.dumps(key).encode()).hexdigest()[:20]
            entries[key] = dict(id=digest, products=products, source=str(worker),
                sources=[dict(worker=worker_name,reaction=fields[2])],
                reactions=combined.reac_name, settings=settings,
                l1_asymptote=selection,
                source_summary_sha256=file_sha256(summary))
    return dict(schema=1, case=str(case), entries=list(entries.values()))


def _seed_fragment_qc(source, target, products):
    """Copy only fragment observations/artifacts; retain provenance, not jobs."""
    from ase.db import connect
    from kinbot.run_format import ensure_current_run, MARKER
    ensure_current_run(source)
    shutil.copy2(source/MARKER,target/MARKER)
    def database_generation():
        return {p.name:file_sha256(p) for p in
                (source/'kinbot.db',source/'kinbot.db-wal') if p.is_file()}
    generation = database_generation()
    def selected(name):
        return any(Path(name.removeprefix('stereochemistry/')).name.startswith(p+'_')
                   for p in products)
    source_db = _database(source)
    rows=[]
    try:
        rows=list(source_db.select())
        with connect(str(target/'kinbot.db')) as dest:
            for row in rows:
                if selected(row.get('name','')):
                    dest.write(row.toatoms(), key_value_pairs=row.key_value_pairs,
                               data=dict(row.data))
    finally:
        source_db.connection.close()
    files=[]
    for sub in ('','conf','hir','perm','perm/conf','perm/hir'):
        directory=source/sub
        if not directory.is_dir():continue
        for path in directory.iterdir():
            if not path.is_file() or path.is_symlink() or not selected(path.name):continue
            # Executable/scheduler ownership belongs to the old directory.
            if path.suffix in ('.py','.sbatch','.slurm','.pbs') or path.name.endswith('.runtime.json'):
                continue
            rel=path.relative_to(source)
            (target/rel).parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(path,target/rel)
            digest=file_sha256(target/rel)
            if digest != file_sha256(path):
                raise RuntimeError(f'{path}: changed during reverse-search staging; source work must be quiescent.')
            files.append(dict(path=str(rel),sha256=digest))
    route=json.loads((source/'.kinbot_theory_jobs.json').read_text())
    route['jobs']={name:{**entry,'job_id':None} for name,entry in route['jobs'].items()
                   if selected(name)}
    atomic_json(target/'.kinbot_theory_jobs.json',route)
    if generation != database_generation():
        raise RuntimeError('Source fragment database changed during staging; retry preparation after its owner finishes.')
    atomic_json(target/'fragment_seed.json',dict(source=str(source),files=files,
        source_database_sha256=generation['kinbot.db'],source_database_generation=generation,products=products))


def prepare(case):
    case=Path(case).resolve(strict=True)
    base=case/'reverse_abstractions'
    base.mkdir(exist_ok=True)
    with (base/'prepare.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        marker=base/'manifest.json'
        data=json.loads(marker.read_text()) if marker.exists() else plan(case)
        for entry in data['entries']:
            if entry.get('input') and file_sha256(base/entry['input'])!=entry['input_sha256']:
                raise ValueError('Staged reverse input changed; refusing an implicit revision.')
        if data.get('status') in ('staged_for_discovery_not_integrated','integrated'):return data
        data['status']='preparing'
        atomic_json(marker,data)
        for entry in data['entries']:
            if entry.get('input'):continue
            if not entry['reactions']:
                entry['status']='no_abstraction_motif'
                continue
            directory=base/entry['id']
            temporary=base/(entry['id']+'.staging')
            settings=deepcopy(entry['settings'])
            settings['input_file']=str(directory/'input.json')
            if directory.exists():
                # Crash after atomic rename but before manifest publication.
                if json.loads((directory/'input.json').read_text())!=settings or not (directory/'fragment_seed.json').is_file():
                    raise RuntimeError(f'{directory}: unrecognized staging contents.')
            else:
                if temporary.exists():
                    import uuid
                    os.rename(temporary,base/(entry['id']+'.interrupted.'+uuid.uuid4().hex))
                temporary.mkdir()
                _seed_fragment_qc(Path(entry['source']),temporary,entry['products'])
                atomic_json(temporary/'input.json',settings)
                os.rename(temporary,directory)
            entry.pop('settings')
            entry.update(status='staged',input=entry['id']+'/input.json',
                         input_sha256=file_sha256(directory/'input.json'))
            atomic_json(marker,data)
        data['status']='staged_for_discovery_not_integrated'
        atomic_json(marker,data)
        return data


def inputs(case):
    base=Path(case)/'reverse_abstractions'
    marker=base/'manifest.json'
    if not marker.exists():return []
    data=json.loads(marker.read_text())
    if data.get('status') not in ('staged_for_discovery_not_integrated','integrated'):
        raise RuntimeError('Reverse abstraction staging is incomplete; rerun its preparation first.')
    paths=[(base/e['input']).resolve(strict=True) for e in data['entries'] if e.get('input')]
    if any(not p.is_relative_to(base.resolve()) for p in paths):
        raise ValueError('Reverse input escapes its case.')
    return paths


def ensure_integrated(case):
    marker=Path(case)/'reverse_abstractions'/'manifest.json'
    if not marker.exists():
        return
    manifest = json.loads(marker.read_text())
    if manifest.get('status') == 'integrated':
        return
    # An exhausted reverse search is not a missing chemical reaction. When
    # every requested motif has a current terminal nonfinding, there is no new
    # saddle to integrate. Do not require somebody to invent one or edit the
    # manifest to let the original PES proceed.
    from kinbot.search_state import NAME, validate_input
    if manifest.get('status') == 'staged_for_discovery_not_integrated':
        clear = True
        for entry in manifest['entries']:
            if not entry.get('input'):
                clear &= entry.get('status') == 'no_abstraction_motif' and not entry['reactions']
                continue
            inp = (marker.parent/entry['input']).resolve(strict=True)
            if not inp.is_relative_to(marker.parent.resolve()):
                raise ValueError('Reverse input escapes its case.')
            directory = inp.parent
            paths = [directory/name for name in
                     ('kinbot_status.json', 'reverse_search_result.json', NAME)]
            if not all(p.is_file() for p in paths):
                clear = False
                continue
            controller, result, state = [json.loads(p.read_text()) for p in paths]
            decisions = {r['reaction']: r['status'] for r in result['reactions']}
            if (controller.get('status') != 'complete' or controller.get('returncode') != 0
                    or result.get('status') != 'discovery_complete_needs_integration'
                    or not set(entry['reactions']).issubset(decisions)
                    or any(value != -999 for value in decisions.values())
                    or set(state.get('reactions', {})) != set(decisions)
                    or any(r.get('status') != -999 or r.get('reassess')
                           for r in state.get('reactions', {}).values())):
                clear = False
                continue
            validate_input(inp)  # A stale result cannot settle a new revision.
        if clear:
            return
    raise RuntimeError('Reverse abstraction discovery is unfinished or has accepted results '
        'requiring entrance-well matching and integration before replacing the combined MESS '
        'model. See reverse_abstractions/manifest.json and reverse_search_result.json; old MESS '
        'is stale. Completed searches with no retained reaction do not require integration '
        'or an automatic retry.')


def launch(case, code):
    """The ordinary refresh driver owns these controllers too."""
    from kinbot.run_cleanup import _processes
    for path in inputs(case):
        marker=path.parent/'launch.json'
        if marker.exists() or _processes(path.parent):continue
        env=dict(os.environ, PYTHONPATH=str(code),
            OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
        with (path.parent/'driver.log').open('ab') as log:
            child=subprocess.Popen([sys.executable,'-u','-B','-m','kinbot.kb',str(path)],
                cwd=path.parent,env=env,stdin=subprocess.DEVNULL,
                stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        atomic_json(marker,dict(pid=child.pid,host=socket.gethostname(),input=str(path),code=str(code)))


def characterize_entrance(endpoint, reactants, par, qc):
    """Obtain matched L1 depth and L2 data for a bound entrance complex.

    Native failures are recorded, not converted to barrierless verdicts or
    automatically retried. Reentry adopts the same named QC observations.
    """
    from kinbot.l1_endpoint import ensure
    from kinbot.optimize import Optimize
    while True:
        checked = ensure(qc, endpoint, par, family='abstraction')
        if checked['status'] != 'waiting':
            break
        time.sleep(1)
    if checked['status'] != 'complete':
        return dict(status='unresolved', l1_check=checked)
    point = checked['point']
    fragments = l1_screening.read_points(Path.cwd(),
        [routing_name(p)+'_well' for p in reactants], l1_screening.theory(par))
    depth = (sum(p['energy_hartree']+p['zpe_hartree'] for p in fragments)
             -point['energy_hartree']-point['zpe_hartree'])*constants.AUtoKCAL
    result = dict(l1_check=checked, fragment_points=fragments, depth_kcal_mol=depth)
    if depth <= par['vdW_detection']:
        return dict(result, status='no_retained_entrance_complex',
                    reason='Matched L1 depth does not exceed vdW_detection; no artificial MESS well created.')
    endpoint.energy = point['energy_hartree']
    endpoint.zpe = point['zpe_hartree']
    endpoint.freq = point['frequencies']
    optimization = Optimize(endpoint, par, qc, just_high=True)
    while optimization.shigh not in (1, -999):
        optimization.do_optimization()
        if optimization.shigh not in (1, -999):
            time.sleep(1)
    if optimization.shigh == -999:
        return dict(result, status='unresolved', reason='Entrance L2 optimization failed its native validation.')
    return dict(result, status='characterized', l2=_point(optimization.species))


def _point(p):
    result = dict(name=p.name,symbols=list(map(str,p.atom)),positions=np.asarray(p.geom).tolist(),
        charge=int(p.charge),multiplicity=int(p.mult),energy_hartree=float(p.energy),
        zpe_hartree=float(p.zpe),frequencies=list(map(float,p.freq)))
    if (not np.isfinite(result['positions']).all()
            or not all(np.isfinite(result[key]) for key in ('energy_hartree','zpe_hartree'))
            or not np.isfinite(result['frequencies']).all()):
        raise ValueError('Nonfinite reverse-search stationary-point data.')
    return result


def export_result(species, par, qc):
    from kinbot.reaction_path import reaction_path_id
    records=[]
    for i,obj in enumerate(species.reac_obj):
        item=dict(reaction=obj.instance_name,status=int(species.reac_ts_done[i]))
        if species.reac_ts_done[i]==-1:
            if not hasattr(obj,'irc_reactant'):
                raise ValueError('Accepted reverse abstraction has no retained entrance IRC endpoint.')
            # Product objects can be renamed to species IDs and optimized in
            # place. Read the native endpoints by their preserved job IDs,
            # never by the later product name or its L2 coordinates.
            jobs = obj.irc_endpoint_jobs
            points = {side:l1_screening.read_points(Path.cwd(), [job],
                l1_screening.theory(par), geometry=True, require_zpe=False)[0]
                for side, job in jobs.items()}
            entrance = StationaryPoint(jobs['reactant'], species.charge, species.mult,
                atom=points['reactant']['symbols'], geom=np.asarray(points['reactant']['positions'], dtype=float))
            entrance.characterize()
            item.update(ts_l2=_point(obj.ts_opt.species),products_l2=[_point(o.species) for o in obj.prod_opt],
                stereopath=reaction_path_id(obj),
                entrance_l1_job=jobs['reactant'], entrance_l1=points['reactant'],
                product_l1_job=jobs['product'])
            item['entrance'] = characterize_entrance(entrance,
                (species.fragA,species.fragB),par,qc)
            item['product_complex_l2'] = (_point(obj.irc_prod_opt.species)
                if obj.do_vdW and obj.irc_prod_opt.shigh==1 else None)
            item['product_endpoint_l1'] = points['product']
        records.append(item)
    atomic_json('reverse_search_result.json',dict(schema=1,status='discovery_complete_needs_integration',
        source_reference=par.get('bimol_pes_reference'),reactions=records,
        reactants_l2=[_point(species.fragA),_point(species.fragB)],
        database_sha256=file_sha256('kinbot.db'),
        mess_ready=False,reason='Match and characterize entrance complex before adding a well-to-product barrier.'))


def status(case):
    """Read controller/result receipts without querying or changing QC jobs."""
    base=Path(case)/'reverse_abstractions'
    data=json.loads((base/'manifest.json').read_text())
    result=[]
    for entry in data['entries']:
        row=dict(id=entry['id'],products=entry['products'],searches=len(entry['reactions']))
        if not entry.get('input'):
            row['status']=entry.get('status','not_staged')
        else:
            directory=(base/entry['input']).parent
            native=directory/'reverse_search_result.json'
            controller=directory/'kinbot_status.json'
            if native.is_file():
                report=json.loads(native.read_text())
                row.update(status=report['status'],result=str(native),
                    accepted=sum(r['status']==-1 for r in report['reactions']),
                    entrance_statuses=[r['entrance']['status'] for r in report['reactions'] if r['status']==-1])
            elif controller.is_file():
                report=json.loads(controller.read_text())
                row.update(status=report.get('status','unknown'),controller=str(controller))
            else:
                row['status']='launched_without_status' if (directory/'launch.json').exists() else 'staged'
        result.append(row)
    return dict(case=str(Path(case).resolve()),status=data['status'],pairs=result)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',type=Path)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--prepare',action='store_true')
    mode.add_argument('--status',action='store_true',help='Read discovery progress and native result locations.')
    args=parser.parse_args(argv)
    if args.status:
        print(json.dumps(status(args.case),indent=2))
        return
    data=prepare(args.case) if args.prepare else plan(args.case)
    print(json.dumps(dict(case=data['case'],pairs=len(data['entries']),
        searches=sum(len(e['reactions']) for e in data['entries']),
        pairs_with_abstraction=sum(bool(e['reactions']) for e in data['entries']),
        status=data.get('status','read_only_plan')),indent=2))


if __name__=='__main__':main()
