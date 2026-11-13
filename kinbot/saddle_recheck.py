"""Bounded, independently restartable saddle/IRC checks at an explicit theory.

Use the ordinary QC, reaction-guess and IRC owners. Results are observations
for later PES integration. A target label never determines the observed path.
"""
import argparse
from contextlib import contextmanager
from dataclasses import asdict
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import time
from types import SimpleNamespace

import numpy as np
from ase import Atoms

from kinbot.execution import atomic_json, file_sha256
from kinbot.frequencies import frequency_validation_error
from kinbot.gaussian_modes import read_normal_modes
from kinbot.irc import IRC
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.reac_family import carry_out_reaction
from kinbot.reaction_guess import alternative_name
from kinbot.reactions.reac_R_Addition_MultipleBond import RAdditionMultipleBond
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import require_supported_identity
from kinbot.theory import resolve_profiles


@contextmanager
def in_directory(directory):
    previous = Path.cwd()
    os.chdir(directory)
    try:
        yield
    finally:
        os.chdir(previous)


def point(record, name):
    species = StationaryPoint(name, record['charge'], record['multiplicity'],
        atom=np.array(record['symbols']), geom=np.array(record['positions'], float))
    species.characterize()
    return species


def validate(manifest):
    if manifest.get('schema') != 1 or not manifest.get('cases'):
        raise ValueError('Saddle recheck needs schema 1 and explicit cases.')
    names = []
    for case in manifest['cases']:
        name = case['name']
        if not re.fullmatch(r'[A-Za-z0-9_-]+', name):
            raise ValueError('Case names must be single portable directory names.')
        names.append(name)
        if case['mode'] not in ('optimize', 'reuse', 'family'):
            raise ValueError('Unknown saddle recheck mode.')
        seed = case['seed']
        geometry = np.asarray(seed['positions'], float)
        if geometry.shape != (len(seed['symbols']), 3) or not np.isfinite(geometry).all():
            raise ValueError('A seed must have a complete finite geometry.')
        if type(seed['charge']) is not int or type(seed['multiplicity']) is not int or seed['multiplicity'] < 1:
            raise ValueError('Each seed needs an explicit integer charge and multiplicity.')
        if case['mode'] == 'reuse' and not case.get('native'):
            raise ValueError('Saddle reuse requires hash-bound native evidence.')
        if case['mode'] == 'family':
            from kinbot.reaction_guess import validate as validate_guess
            validate_guess({case['source_reaction']: case['guess']})
            motif = case['motif']
            if len(motif) != 3 or len(set(motif)) != 3 or any(
                    type(i) is not int or i < 0 or i >= len(geometry) for i in motif):
                raise ValueError('Invalid three-atom family motif.')
    if len(set(names)) != len(names):
        raise ValueError('Duplicate recheck case name.')
    return manifest


def stage(manifest_path, root, source_root):
    manifest_path, root = Path(manifest_path).resolve(), Path(root).resolve()
    source_root = Path(source_root).resolve()
    manifest = validate(json.loads(manifest_path.read_text()))
    root.mkdir(parents=True, exist_ok=True)
    frozen = root/'campaign.json'
    if frozen.exists() and json.loads(frozen.read_text()) != manifest:
        raise ValueError('Existing campaign has different scientific inputs. Use a new directory.')
    atomic_json(frozen, manifest)
    for case in manifest['cases']:
        folder = root/case['name']
        folder.mkdir(exist_ok=True)
        case_path = folder/'case.json'
        if case_path.exists() and json.loads(case_path.read_text()) != case:
            raise ValueError(f'{folder}: existing seed definition differs.')
        atomic_json(case_path, case)
        parameters = dict(manifest['parameters'])
        if case['mode'] == 'family':
            parameters['reaction_search_guesses'] = {case['source_reaction']: case['guess']}
        path = folder/'input.json'
        if path.exists() and json.loads(path.read_text()) != parameters:
            raise ValueError(f'{folder}: existing input differs.')
        atomic_json(path, parameters)
        # Native reuse is a copy into this independent diagnostic calculation.
        # Source calculations and their databases are never modified.
        for suffix, receipt in case.get('native', {}).get('files', {}).items():
            if suffix not in ('.com', '.log'):
                raise ValueError('Native saddle evidence accepts .com and .log only.')
            source = (source_root/receipt['path']).resolve()
            if not source.is_relative_to(source_root) or file_sha256(source) != receipt['sha256']:
                raise ValueError(f'{source}: native source hash differs.')
            destination = folder/('source'+suffix)
            if destination.exists():
                if file_sha256(destination) != receipt['sha256']:
                    raise ValueError(f'{destination}: copied native evidence differs.')
            else:
                shutil.copy2(source, destination)
    return root


def adopt(case, qc, parameters):
    """Reuse one verified, fixed-theory saddle with its native Hessian evidence."""
    from kinbot.l1_screening import normalized_profile
    native, seed = case['native'], case['seed']
    target = normalized_profile(asdict(resolve_profiles(parameters)['l1']))
    if normalized_profile(native['profile']) != target:
        raise ValueError('Native saddle theory differs from the requested check theory.')
    for suffix, receipt in native['files'].items():
        if file_sha256('source'+suffix) != receipt['sha256']:
            raise ValueError('Copied source evidence changed.')
    modes = read_normal_modes('source.log', seed['symbols'], seed['positions'])
    data = native['data']
    reason = frequency_validation_error(modes['frequencies_cm1'], 1, parameters['imagfreq_threshold'])
    if reason or data['status'] != 'normal' or not np.isfinite([data['energy'], data['zpe']]).all():
        raise ValueError(reason or 'Invalid native saddle result.')
    if not np.allclose(modes['frequencies_cm1'], data['frequencies'], atol=.01, rtol=0):
        raise ValueError('Native frequencies differ from the selected source row.')
    if not qc.db.count(name='saddle'):
        qc.db.write(Atoms(seed['symbols'], seed['positions']), name='saddle',
                    data=dict(data, saddle_recheck_source=native))
    if not Path('saddle.log').exists():
        shutil.copy2('source.log', 'saddle.log')
    atomic_json('saddle.reused.json', dict(source=native, native_modes=modes))
    # This is an explicit import receipt, with no claim of a new QC execution.
    Path('saddle_sella.log').write_text('Adopted hash-verified native saddle; see saddle.reused.json\ndone\n')


def observation(qc, job, seed):
    rows = list(qc.db.select(name=job))
    if not rows:
        return None
    row = rows[-1]
    data = {key:row.data[key] for key in ('status','error','energy','zpe') if key in row.data}
    if 'frequencies' in row.data:
        data['frequencies'] = list(map(float,row.data['frequencies']))
    result = dict(job=job, row_id=row.id, symbols=list(row.symbols),
                  positions=row.positions.tolist(), data=data)
    try:
        species = point(dict(seed, symbols=result['symbols'], positions=result['positions']), job)
        result['identity'] = require_supported_identity(species)
        result['smiles'] = species.smiles
    except ValueError as error:
        result['identity_error'] = str(error)
    if data.get('frequencies'):
        try:
            result['native_modes'] = read_normal_modes(job+'.log',row.symbols,row.positions)
        except (OSError,ValueError) as error:
            result['native_mode_error'] = str(error)
    return result


def tick(case, qc, parameters):
    state_path = Path('recheck.json')
    state = json.loads(state_path.read_text()) if state_path.exists() else dict(schema=1, step=0)
    if state.get('complete'):
        return state
    species = point(case['seed'], 'saddle')
    job = 'saddle'
    if case['mode'] == 'family':
        job = alternative_name(case['source_reaction'], case['guess'])
        reaction = RAdditionMultipleBond(species, qc, parameters, case['motif'], job)
        from kinbot.search_state import SearchState
        # Ordinary PES discovery has already submitted a parent optimization.
        # These checks can start directly with a path, so freeze the explicit
        # profile before recording its first submission checkpoint.
        if hasattr(qc,'_save_routes') and not Path('.kinbot_theory_jobs.json').exists():
            qc._save_routes()
        species.reac_obj = [reaction]
        species.reac_type = ['R_Addition_MultipleBond']
        species.reac_ts_done = [0]
        species.reac_step = [state['step']]
        species.reac_scan_energy = [[]]
        checkpoint = SearchState(SimpleNamespace(species=species,qc=qc,par=parameters,inp='input.json'))
        step = species.reac_step[0]
        if step <= reaction.max_step:
            if species.reac_ts_done[0] == -999 or qc.check_qc(job) == 'error':
                state.update(complete=True, outcome='saddle_search_not_retained',
                             saddle=observation(qc,job,case['seed']))
                species.reac_ts_done[0] = -999
                checkpoint.save(0)
                atomic_json(state_path,state)
                return state
            step = carry_out_reaction(reaction, step, parameters['qc_command'])
            state['step'] = step
            species.reac_step[0] = step
            species.reac_ts_done[0] = -999 if step < 0 else 1 if step > reaction.max_step else 0
            checkpoint.save(0)
            if step < 0:
                state.update(complete=True, outcome='saddle_search_not_retained')
            atomic_json(state_path, state)
            return state
    elif case['mode'] == 'reuse' and not Path('saddle.reused.json').exists():
        adopt(case, qc, parameters)
    status = qc.check_qc(job)
    if status == 0 and case['mode'] == 'optimize':
        qc.qc_opt_ts(species, species.geom)
        status = 'running'
    if status in (0, 'running'):
        state.update(stage='saddle', qc_status=status)
    elif status == 'error':
        state.update(complete=True, outcome='saddle_not_validated', saddle=observation(qc,job,case['seed']))
    else:
        saddle = observation(qc, job, case['seed'])
        reason = frequency_validation_error(saddle['data'].get('frequencies', []), 1,
                                            parameters['imagfreq_threshold'])
        if reason:
            state.update(complete=True, outcome='saddle_not_validated', reason=reason, saddle=saddle)
        else:
            reaction = SimpleNamespace(instance_name=job, qc=qc, species=species)
            IRC(reaction, parameters).do_irc_calculations()
            statuses = {d:qc.check_qc(job+'_IRC_'+d) for d in ('F','R')}
            state.update(stage='irc', saddle=saddle, irc_statuses=statuses)
            if all(s in ('normal','error') for s in statuses.values()):
                state['endpoints'] = {d:observation(qc,job+'_IRC_'+d+'_prod',case['seed']) for d in statuses}
                state['irc_paths'] = {d:json.loads(Path(job+'_IRC_'+d+'.irc_path.json').read_text())
                                      for d in statuses if Path(job+'_IRC_'+d+'.irc_path.json').exists()}
                valid = (all(p and p['data'].get('status')=='normal' for p in state['endpoints'].values())
                         and all(statuses[d]=='normal' and state['irc_paths'].get(d,{}).get('converged') for d in statuses))
                state.update(complete=True, outcome='irc_endpoints_available' if valid else 'irc_incomplete')
    atomic_json(state_path, state)
    return state


def run(root, interval=10):
    root = Path(root).resolve()
    manifest = validate(json.loads((root/'campaign.json').read_text()))
    with (root/'controller.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print(f'Existing controller owns {root}; continuing with that owner.')
            return 0
        from kinbot.driver_status import run_recorded
        def execute():
            contexts = {}
            for case in manifest['cases']:
                with in_directory(root/case['name']):
                    if json.loads(Path('case.json').read_text()) != case:
                        raise ValueError('Seed definition changed after staging.')
                    expected = dict(manifest['parameters'])
                    if case['mode']=='family':
                        expected['reaction_search_guesses'] = {case['source_reaction']:case['guess']}
                    if json.loads(Path('input.json').read_text()) != expected:
                        raise ValueError('Effective recheck input changed after staging.')
                    parameters = Parameters('input.json',show_warnings=False).par
                    contexts[case['name']] = (QuantumChemistry(parameters), parameters)
            while True:
                results = {}
                for case in manifest['cases']:
                    with in_directory(root/case['name']):
                        results[case['name']] = tick(case, *contexts[case['name']])
                atomic_json(root/'results.json', results)
                if all(r.get('complete') for r in results.values()):
                    return 0
                time.sleep(interval)
        with in_directory(root):
            return run_recorded(execute, 'controller.json')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    stage_parser = sub.add_parser('stage')
    stage_parser.add_argument('manifest',type=Path)
    stage_parser.add_argument('root',type=Path)
    stage_parser.add_argument('--source-root',type=Path,required=True)
    for name in ('run','status'):
        child=sub.add_parser(name)
        child.add_argument('root',type=Path)
        if name=='status':
            child.add_argument('--details',action='store_true')
    args=parser.parse_args(argv)
    if args.action=='stage':
        print(stage(args.manifest,args.root,args.source_root))
    elif args.action=='run':
        run(args.root)
    else:
        manifest=validate(json.loads((args.root/'campaign.json').read_text()))
        for case in manifest['cases']:
            path=args.root/case['name']/'recheck.json'
            record=json.loads(path.read_text()) if path.exists() else {'stage':'not_started'}
            display=record if args.details else {k:record[k] for k in
                ('stage','step','outcome','complete','irc_statuses') if k in record}
            print(case['name'],json.dumps(display,sort_keys=True))


if __name__=='__main__':
    main()
