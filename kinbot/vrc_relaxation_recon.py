"""General finite relaxation panels from explicit orbital selections.

These graphs never publish production models. Only the separation is fixed;
internal fragment rearrangements remain visible in final.xyz and the report.
"""
import argparse
from contextlib import ExitStack
from copy import deepcopy
import fcntl
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from ase.io import read as atoms_read
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, atomic_text, file_sha256
from kinbot import molpro_reconnaissance as recon, reaction_decisions, run_policy
from kinbot.stationary_pt import StationaryPoint
from kinbot.electronic_states import _canonical_smiles
from kinbot.molpro import validate_vrc_state_model
from kinbot.parameters import Parameters
from kinbot.vrc_graph import select_graph, scientific_spec


def read(path):
    return json.loads(Path(path).read_text())


def native_evidence(manifest):
    """Calculation identity, excluding mutable scheduling/derived audit fields.

    Input/output/WFU hashes bind the native calculation; geometry, model and
    seed provenance are also explicit. Neither exclusions nor scheduler status
    identify a new calculation. This is not a scientific acceptance verdict.
    """
    fields = ('name', 'directory', 'basis', 'anchor', 'input_sha256', 'template_sha256',
              'probe_active_space', 'probe_method', 'state_model', 'root_states',
              'symbols', 'geometry_angstrom', 'distance_angstrom', 'spin_context',
              'channel', 'orientation', 'correlation_checks', 'orbital_predecessor',
              'initial_record', 'initial_wfu_sha256', 'wfu', 'wfu_sha256',
              'output_sha256', 'reference_review_contract_sha256')
    return dict(schema=manifest.get('schema'), source=manifest.get('source'),
                electron_count=manifest.get('electron_count'),
                jobs=[{k: j.get(k) for k in fields} for j in manifest['jobs']])


def _panel_evidence(specs, output, evidence_snapshots):
    """Read current manifests once and find hash-bound historical snapshots."""
    paths = {str(Path(read(Path(s['intent']['relaxation_recon']['source_screen'])/
                          'review_launch.json')['systems'][0]['manifest']).resolve())
             for s in specs}
    current = {p: Path(p).read_bytes() for p in paths}
    snapshots = {}
    for p in [*evidence_snapshots, *sorted((output/'.native_evidence').glob('*.json'))]:
        data = Path(p).read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        if Path(p).parent == output/'.native_evidence' and Path(p).stem != digest:
            raise ValueError('Saved native manifest snapshot changed: '+str(p))
        snapshots[digest] = data
    for s in specs:
        evidence = s['intent']['relaxation_recon']['evidence']
        for path in paths & evidence.keys():
            if hashlib.sha256(current[path]).hexdigest() != evidence[path]:
                raise ValueError('Native manifest changed during panel validation; repeat at an idle boundary.')
    snapshots.update({hashlib.sha256(data).hexdigest(): data for data in current.values()})
    return current, snapshots


def _same_panel(saved, wanted, current, snapshots):
    before, after = scientific_spec(saved), scientific_spec(wanted)
    old = before['intent']['relaxation_recon']['evidence']
    new = after['intent']['relaxation_recon']['evidence']
    for path, data in current.items():
        if path not in old or path not in new or old[path] == new[path]:
            continue
        historical = snapshots.get(old[path])
        if historical is None:
            raise ValueError('Missing hash-matched native manifest snapshot: '+path)
        if native_evidence(json.loads(historical)) != native_evidence(json.loads(data)):
            raise ValueError('Native calculation evidence changed: '+path)
        # Only this verified manifest hash may differ. All other provenance,
        # geometry, theory, state definitions and constraints remain exact.
        old[path] = new[path]
    return before == after


def specifications(source, request, selection, case_input, *, distances_angstrom=None):
    source, request, selection, case_input = [Path(p).resolve(strict=True) for p in
                                             (source,request,selection,case_input)]
    record = read(source/'review_launch.json')
    for name,digest in record['files'].items():
        if file_sha256(source/name) != digest:
            raise ValueError('Saved orbital evidence changed: '+name)
    current, old, chosen = read(request), read(source/'request.json'), read(selection)
    par = Parameters(str(case_input)).par
    if (current.get('status') != 'needs_review' or chosen.get('reviewed') is not True
            or reaction_decisions.rejected(request.parents[3],current['reaction'])):
        raise ValueError('Panel requires an active review request and explicit orbital selection.')
    def value(key, obj):
        v=obj.get(key)
        if key=='products':
            return sorted((_canonical_smiles(p['smiles']),p['charge'],p['multiplicity']) for p in v)
        return v
    for key in ('parent','reaction','symbols','bond','charge','multiplicity','products'):
        if value(key,current)!=value(key,old) or value(key,current)!=value(key,chosen['request_match']):
            raise ValueError('Orbital panel request mismatch: '+key)
    molecules=[]
    for data in (old,current):
        p=StationaryPoint('probe',data['charge'],data['multiplicity'],
                         atom=np.array(data['symbols']),geom=np.array(data['positions']))
        p.characterize();molecules.append(p)
    if not np.array_equal(molecules[0].bond01,molecules[1].bond01):
        raise ValueError('Parent connectivity/atom order differs from the orbital screen.')
    manifest_path=Path(record['systems'][0]['manifest']).resolve()
    manifest_bytes=manifest_path.read_bytes()
    manifest=json.loads(manifest_bytes);jobs={j['name']:j for j in manifest['jobs']}
    required=chosen.get('required_complete_jobs',list(jobs))
    for name in required:
        job=jobs[name]
        if (job['status']!='complete' or job.get('validation_errors')
                or file_sha256(Path(job['directory'])/(name+'.inp'))!=job['input_sha256']
                or file_sha256(job['wfu'])!=job['wfu_sha256']
                or file_sha256(Path(job['directory'])/(name+'.out'))!=job['output_sha256']):
            raise ValueError('Native orbital screen is incomplete or changed: '+name)
    config=read(record['systems'][0]['configuration'])
    reference=current['positions'];reference_species=molecules[1]
    if (source/'geometry_source.json').exists():
        geometry=read(source/'geometry_source.json')
        if geometry['network']['name']!=current['reaction'] or geometry['symbols']!=current['symbols']:
            raise ValueError('Product-side orbital provenance differs from this request.')
        reference=geometry['positions']
        reference_species=StationaryPoint('product',current['charge'],current['multiplicity'],
            atom=np.array(current['symbols']),geom=np.array(reference))
        reference_species.characterize()
    config.update(geometry_angstrom=reference,multiplicity=current['multiplicity'],kinbot_input=str(case_input))
    anchors=recon._coordinates(config,source/'input.json') if 'coordinates' in config else config['anchors']
    if distances_angstrom is not None:
        if (not isinstance(distances_angstrom,list) or not distances_angstrom
                or any(type(d) not in (int,float) or not np.isfinite(d) or d<=0
                       for d in distances_angstrom)
                or len(set(distances_angstrom))!=len(distances_angstrom)):
            raise ValueError('Relaxation distances must be distinct positive finite numbers.')
        coordinates={(tuple(a['bond']),tuple(a['moving_atoms'])) for a in anchors}
        if len(coordinates)!=1 or any(a.get('_rotation') for a in anchors):
            raise ValueError('A distance override requires one unrotated separation coordinate.')
        anchors=[dict(deepcopy(anchors[0]),distance_angstrom=d) for d in sorted(distances_angstrom)]
    basis=par['vrc_tst_scan_basis'];seed=chosen['basis_seeds'][basis];job=jobs[seed['job']]
    if job['basis']!=basis or seed['job'] not in required:
        raise ValueError('Panel needs an explicitly verified basis-matched seed.')
    model=deepcopy(chosen['state_model'])
    model['orbital_selection']=dict(mode='molpro_record',path=job['wfu'],sha256=job['wfu_sha256'],
        record=seed['record'],orbital_type='natural',align_active_orbitals=False,active_indices=seed['active_indices'])
    active=chosen['active_space']
    model=validate_vrc_state_model(model,active_electrons=active['active_electrons'],
        active_orbitals=active['active_orbitals'],nelectron=manifest['electron_count'],atom_count=len(current['symbols']))
    method=chosen['parameters']['vrc_tst_scan_method']
    specs=[]
    for anchor in anchors:
        pair=[i-1 for i in anchor['bond']]
        moving=[i-1 for i in anchor['moving_atoms']]
        maps=[sorted(set(range(len(current['symbols'])))-set(moving)),moving]
        proof=dict(reaction=current['reaction'],source_screen=str(source),request=str(request),
            distance_angstrom=anchor['distance_angstrom'],maps=maps,coordinate=pair,
            reference_positions=reference,reference_bonds=reference_species.bond01.tolist(),
            evidence={str(request):file_sha256(request),str(selection):file_sha256(selection),
                      str(manifest_path):hashlib.sha256(manifest_bytes).hexdigest()})
        task=dict(id='relaxed_scan',kind='ase_optimize',geometry_from='initial',geometry_output='final.xyz',
            profile=dict(calculator='molpro',method=method,basis=basis,command='molpro',optimizer='sella',
                calculator_kwargs=dict(state_model=model,caspt2_shift=par['vrc_tst_scan_shift'],
                    caspt2_ipea=par['vrc_tst_scan_ipea'],gradient_mode=par['vrc_tst_scan_gradient'],
                    scf_convergence=par['vrc_tst_scf_convergence'],correlation_convergence=par['vrc_tst_correlation_convergence'])),
            optimizer=dict(fmax=par['vrc_tst_scan_fmax'],steps=50,step_retries=0,order=0,
                sella_kwargs=dict(internal=True),constraints=dict(fixed_bonds=[pair])),
            resources=dict(cores=8,memory_mb=104000,min_stack_mw=1000,partition=par['queue_name'],walltime='04:00:00'))
        spec=dict(schema=1,name='relaxed-orbital-recon',molecule=dict(symbols=current['symbols'],
            positions=recon._anchor_geometry(config,anchor).tolist(),charge=current['charge'],multiplicity=current['multiplicity']),
            limits=dict(max_nodes=1),tasks=[task],intent=dict(relaxation_recon=proof,production_model_accepted=False,
                vrc_fragment_retry=True,cap_walltime_to_explicit_partition=True))
        dispatch.validate_spec(spec);specs.append(spec)
    return specs


def prepare(production,plan,output, *, validate_only=False, evidence_snapshots=()):
    production,plan=Path(production).resolve(strict=True),Path(plan).resolve(strict=True)
    output=Path(output).resolve()
    if output.is_relative_to(production):raise ValueError('Keep diagnostic panels outside production.')
    policy=run_policy.read(production) or {};slots=policy.get('orbital_review_reserved_slots',0)
    if slots<1:raise ValueError('No review slots reserved.')
    entries=read(plan)['panels']; specs=[]
    for entry in entries:
        specs.extend(specifications(production.parent/entry['source'],production/entry['request'],
            plan.parent/entry['selection'],production/entry['case_input'],
            distances_angstrom=entry.get('distances_angstrom')))
    with ExitStack() as locks:
        if not validate_only:
            output.mkdir(parents=True,exist_ok=True)
            lock=locks.enter_context((output/'.prepare.lock').open('a'))
            fcntl.flock(lock,fcntl.LOCK_EX)
        current,snapshots=_panel_evidence(specs,output,evidence_snapshots)
        receipt=output/'panel.json'
        graphs=[]
        if receipt.exists():
            saved=read(receipt)
            graphs=saved['graphs']
            if (saved['plan_sha256']!=file_sha256(plan)
                    or saved['production']!=str(production)
                    or len(graphs)!=len(specs) or len(set(graphs))!=len(graphs)
                    or any(Path(g).resolve().parent!=output for g in graphs)
                    or not all(_same_panel(dispatch._load(g)[1],s,current,snapshots)
                               for g,s in zip(graphs,specs))):
                raise ValueError('Existing panel plan or evidence changed; inspect before preparing another panel.')
        elif any(output.glob('probe_*/workflow.json')):
            raise ValueError('Incomplete panel preparation; inspect saved graphs before replacing owners.')
        if validate_only:
            return dict(graphs=graphs,production_model_accepted=False)
        # Preserve historical bytes and existing graph/workflow/job identities.
        # A legacy panel may use a hash-matched candidate source snapshot.
        for digest,data in snapshots.items():
            path=output/'.native_evidence'/(digest+'.json')
            path.parent.mkdir(exist_ok=True)
            if not path.exists():
                atomic_text(path,data.decode('utf-8'))
        if not (output/'queue_all_policy.json').exists():
            run_policy.update(output,max_jobs=min(3,slots),exclude=policy['excluded_nodes'],manage_pending=True)
        if not receipt.exists():
            graphs=[str(select_graph(output,'probe',spec)) for spec in specs]
        atomic_json(receipt,dict(schema=1,production=str(production),plan_sha256=file_sha256(plan),graphs=graphs))
    return dict(graphs=graphs,production_model_accepted=False)


def report(output):
    result=[]
    for graph in read(Path(output)/'panel.json')['graphs']:
        g,spec,state=dispatch._load(graph);proof=spec['intent']['relaxation_recon']
        row=dict(graph=graph,reaction=proof['reaction'],distance=proof['distance_angstrom'],task=state['tasks']['relaxed_scan'])
        final=g/'tasks/relaxed_scan/final.xyz'
        if final.exists():
            execution=read(final.parent/'execution.json')
            if execution.get('status')!='executed':
                row['unaccepted_geometry']=str(final)
                result.append(row)
                continue
            if execution.get('artifacts',{}).get('final.xyz')!=file_sha256(final):
                raise ValueError('Diagnostic final geometry changed.')
            xyz=atoms_read(final).positions;initial=np.array(proof['reference_positions']);bonds=proof['reference_bonds']
            pairs=[]
            for i in range(len(xyz)):
                for j in range(i+1,len(xyz)):
                    if any(i in m and j in m for m in proof['maps']):
                        before=float(np.linalg.norm(initial[i]-initial[j]));after=float(np.linalg.norm(xyz[i]-xyz[j]))
                        pairs.append(dict(atoms=[i+1,j+1],initial=before,final=after,bond=bool(bonds[i][j]),
                                          stretched=bool(bonds[i][j] and after>1.2*before)))
            row['internal_pairs']=pairs
        result.append(row)
    return dict(panels=result,production_model_accepted=False)


def run(output, *, queue_all=False, allow_nodes=()):
    output=Path(output).resolve(strict=True);record=read(output/'panel.json')
    with ExitStack() as locks:
        for p in [output/'.run.lock',*[Path(g)/'drive.lock' for g in record['graphs']]]:
            stream=locks.enter_context(p.open('a'));fcntl.flock(stream,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for graph in record['graphs']:
            dispatch.preflight_until_scheduler_ready(graph)
        if queue_all:
            allowed=set(run_policy._nodes(allow_nodes))
            marker=output/'queue_all_policy.json'
            previous=set(read(marker).get('allowed_nodes',[])) if marker.exists() else set()
            policy=run_policy.read(record['production'])
            run_policy.update(output,max_jobs=len(record['graphs']),
                exclude=set(policy['excluded_nodes'])-allowed,allow=allowed-previous,manage_pending=True)
            atomic_json(marker,dict(allowed_nodes=sorted(allowed),queue_all=True))
            print(f'Queueing all {len(record["graphs"])} prepared relaxed probes; live jobs are adopted.',flush=True)
        elif allow_nodes:
            raise ValueError('Node overrides require explicit --queue-all review scheduling.')
        while True:
            if queue_all:
                slots=run_policy.read(output)['max_jobs']
            else:
                policy=run_policy.read(record['production']);slots=policy.get('orbital_review_reserved_slots',0)
                run_policy.update(output,max_jobs=min(3,slots),exclude=policy['excluded_nodes'],manage_pending=True)
            try:
                states=[dispatch.advance(g,submit=slots>0) for g in record['graphs']]
            except dispatch.SchedulerUnavailableError as error:
                print('Scheduler unavailable; retaining diagnostic jobs:',error,flush=True)
                time.sleep(20)
                continue
            if all(e['status'] in ('complete','failed','unavailable','blocked','skipped')
                   for s in states for e in s['tasks'].values()):break
            time.sleep(20)
    result=report(output);atomic_json(output/'report.json',result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);commands=p.add_subparsers(dest='command',required=True)
    prep=commands.add_parser('prepare');prep.add_argument('production');prep.add_argument('plan');prep.add_argument('output')
    for name in ('run','report'):
        sub=commands.add_parser(name);sub.add_argument('output')
        if name=='run':
            sub.add_argument('--queue-all',action='store_true')
            sub.add_argument('--allow-node',action='append',default=[])
    a=p.parse_args();print(json.dumps(prepare(a.production,a.plan,a.output) if a.command=='prepare'
                                    else run(a.output,queue_all=a.queue_all,allow_nodes=a.allow_node)
                                    if a.command=='run' else report(a.output),indent=2))


if __name__=='__main__':main()
