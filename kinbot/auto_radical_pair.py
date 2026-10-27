"""Opt-in, native-verified single-state CAS(2,2) for saturated radical pairs.

The policy is deliberately finite and restrictive. Unsupported chemistry and
failed scientific tests produce a persistent review request, never a guessed
production partition. The molecular example is not part of this module.
"""
from copy import deepcopy, copy
from pathlib import Path
from types import SimpleNamespace
import json
import logging
import re
import subprocess
import time

import networkx as nx
import numpy as np
from ase.data import atomic_numbers

from kinbot.execution import atomic_json, file_sha256
from kinbot.constants import AUtoKCAL
from kinbot.vrc_reference_review import digest
from kinbot import molpro_reconnaissance as recon
from kinbot import radical_pair_evidence as evidence

logger = logging.getLogger('KinBot')
POLICY = dict(version=1, maximum_valence_orbitals=14,
    potential_tolerance_kcal_mol=.5, angular_tolerance_kcal_mol=.5,
    return_tolerance_kcal_mol=.02, asymptote_tolerance_kcal_mol=.05,
    minimum_subspace_projection=.90, minimum_return_projection=.99,
    minimum_ci_overlap=.99, maximum_ionic_probability=.02,
    maximum_fragment_charge_error=.02, maximum_t1=.03,
    maximum_reference_density_distance=.20, maximum_extra_fractionality=.10,
    maximum_allocation_retries=2)


class ReviewRequired(ValueError):
    pass


def _require(value, message):
    if not value:
        raise ReviewRequired(message)


def _graph(symbols, bonds):
    graph = nx.from_numpy_array(np.asarray(bonds))
    nx.set_node_attributes(graph, dict(enumerate(map(str, symbols))), 'element')
    return graph


def context(parameters, reaction):
    """Check applicability against the accepted products, not their names."""
    from kinbot.molpro import full_valence_active_space
    from kinbot.vrc_correction import spin_splitting_settings, sampling_settings
    p = reaction.species
    symbols = list(map(str, p.atom)); xyz = np.asarray(p.geom, dtype=float)
    bond = list(map(int, reaction.scan_coo))
    _require(int(p.mult) == 1 and int(p.charge) == 0,
             'Automatic radical-pair policy requires a neutral singlet parent.')
    _require(getattr(reaction, 'family_name', 'hom_sci') in ('hom_sci', 'hs') and
             '_hom_sci_' in reaction.instance_name,
             'Automatic radical-pair policy applies only to accepted homolytic scission.')
    _require(set(symbols) <= {'C', 'H'} and len(bond) == 2 and len(set(bond)) == 2
             and all(0 <= i < len(symbols) for i in bond) and 'C' in [symbols[i] for i in bond],
             'Automatic policy requires a C-C or C-H bond in a saturated hydrocarbon.')
    graph = _graph(symbols, p.bond)
    _require(nx.is_tree(graph) and all(v['weight'] == 1 for *_, v in graph.edges(data=True))
             and graph.has_edge(*bond), 'Only acyclic, single-bond hydrocarbon graphs are supported.')
    _require(all(graph.degree(i) == (4 if s == 'C' else 1) for i, s in enumerate(symbols)),
             'Parent valence is not a closed-shell saturated hydrocarbon.')
    graph.remove_edge(*bond)
    components = [set(nx.node_connected_component(graph, i)) for i in bond]
    products = list(reaction.products)
    _require(len(products) == 2 and all(int(x.mult) == 2 and int(x.charge) == 0 for x in products),
             'Accepted products must be two neutral doublet radicals.')
    mappings = []; used = set()
    for product in products:
        pg = _graph(product.atom, product.bond)
        match = None
        for index, component in enumerate(components):
            if index in used:
                continue
            gm = nx.algorithms.isomorphism.GraphMatcher(graph.subgraph(component), pg,
                node_match=lambda a, b: a['element'] == b['element'],
                edge_match=lambda a, b: a['weight'] == b['weight'])
            if gm.is_isomorphic():
                inverse = {v: k for k, v in gm.mapping.items()}
                match = [inverse[i] for i in range(len(product.atom))]; used.add(index); break
        _require(match is not None, 'Accepted product connectivity does not match the severed bond.')
        mappings.append(match)
    _require(xyz.shape == (len(symbols), 3) and np.isfinite(xyz).all(), 'Invalid parent geometry.')
    _require(not parameters.get('vrc_tst_state_average_roots') and
             not parameters.get('vrc_tst_selected_roots') and not parameters.get('vrc_tst_state_model')
             and not parameters.get('vrc_tst_state_model_file'),
             'Explicit electronic-state/orbital overrides require their own review.')
    spin = spin_splitting_settings(parameters, parent_multiplicity=1, fragment_multiplicities=[2, 2])
    sample = sampling_settings(parameters)
    _require(spin['target_multiplicity'] == 1 and spin['reference_multiplicity'] == 3
             and spin['reference']['method'].lower() == 'uccsd(t)-f12b'
             and len(spin['gap_levels']) == 1, 'Automatic policy requires one CASPT2-F12 gap and triplet CC reference.')
    gap = spin['gap_levels'][0]
    _require(gap['method'].lower() == parameters['vrc_tst_scan_method'].lower() ==
             sample['method'].lower() == 'caspt2-f12(2,2)', 'Automatic policy requires CASPT2-F12(2,2) throughout.')
    _require(isinstance(gap['basis'], list) and len(gap['basis']) == 2 and
             gap['basis'] == spin['reference']['basis'], 'MR gap and CC reference must use the same CBS basis pair.')
    shift, ipea = gap.get('shift', .1), gap.get('ipea', .25)
    _require(shift == parameters.get('vrc_tst_scan_shift', .1) and
             ipea == parameters.get('vrc_tst_scan_ipea', .25) and shift >= 0 and ipea >= 0,
             'Scan and gap shift/IPEA controls must match.')
    bases = list(dict.fromkeys([*(sample['basis'] if isinstance(sample['basis'], list) else [sample['basis']]),
                               parameters['vrc_tst_scan_basis'], *gap['basis']]))
    full = full_valence_active_space(symbols, 0, include_light_atoms=True)
    _require(full['active_orbitals'] <= POLICY['maximum_valence_orbitals'],
             'Full valence verification exceeds the automatic policy budget; a reviewed model is required.')
    distances = sorted(set(map(float, parameters['vrc_tst_scan_points'])))
    far = float(parameters.get('vrc_tst_asymptote_distance', 30.))
    _require(len(distances) >= 3 and distances[0] >= 1.8 and distances[-1] < far and far >= 25.,
             'Automatic verification requires scan R>=1.8 A and a separate asymptote >=25 A.')
    # Fixed finite coverage: ends, interior quartiles, and two separated limits.
    radial = sorted(set([distances[0],min(distances,key=lambda x:abs(x-distances[0]-.5)),
                         min(distances,key=lambda x:abs(x-distances[0]-1.5)),
                         *[distances[int((len(distances)-1)*f)] for f in (.25,.5,.75)],
                         distances[-1], far]))
    return dict(policy=deepcopy(POLICY), reaction=reaction.instance_name, symbols=symbols,
                geometry=xyz.tolist(), bond=bond, mappings=mappings, full_space=full,
                electrons=sum(atomic_numbers[s] for s in symbols), bases=bases,
                spin=spin, sampling=sample, radial=radial, distances=distances, far=far,
                scan_basis=parameters['vrc_tst_scan_basis'], shift=shift, ipea=ipea,
                settings={k: v for k, v in parameters.items() if k.startswith('vrc_tst_')})


def move(c, xyz, distance):
    xyz = np.array(xyz, dtype=float, copy=True); a, b = c['bond']
    delta = xyz[b]-xyz[a]; norm = np.linalg.norm(delta)
    _require(norm > .1, 'Degenerate reaction coordinate.')
    group = next(m for m in c['mappings'] if b in m)
    xyz[group] += delta/norm*(distance-norm)
    return xyz


def _base(c, p, title, geometry, symbols=None, mult=1):
    return dict(schema=1, title=title, symbols=symbols or c['symbols'],
        geometry_angstrom=np.asarray(geometry).tolist(), charge=0, multiplicity=mult,
        bases=[dict(name=f'b{i}', basis=b) for i,b in enumerate(c['bases'])],
        ppn=int(p['single_point_ppn']), partition=p['queue_name'],
        queue_job_limit=int(p['queue_job_limit']), orbital_print_virtuals=20,
        orbital_print_threshold=.05, scf_convergence=1.e-10,
        correlation_convergence=1.e-7, scf_options={'max_iterations':200})


def _write_immutable(path, value):
    if path.exists():
        _require(json.loads(path.read_text()) == value, 'Immutable verification definition changed: '+str(path))
    else:
        atomic_json(path, value)


def _run_sources(root, sources, parameters):
    from kinbot import run_policy
    manifests = []
    for name, source in sources.items():
        path = root/(name+'.json'); _write_immutable(path, source)
        manifest=root/name/'reconnaissance.json'
        if not manifest.is_file():
            manifest=recon.prepare(path,root/name)
        saved=json.loads(manifest.read_text())
        _require(Path(saved['source']).resolve()==path.resolve(),'Verification source identity changed.')
        shared=run_policy.find(manifest.parent)
        shared=str(shared) if shared is not None else None
        # Complete pool attachment after a crash between preparation and this
        # write, but never reassign an attempted/live panel to another pool.
        untouched=(saved.get('status')=='prepared' and not saved.get('excluded_jobs') and
                   all(not j.get('job_id') and not j.get('submission_status') and
                       not j.get('attempts') for j in saved['jobs']))
        if shared and not saved.get('shared_policy_file') and untouched:
            saved['shared_policy_file']=shared
            atomic_json(manifest,saved)
        _require(saved.get('shared_policy_file')==shared,'Verification resource-pool identity changed.')
        for j in saved['jobs']:
            prepared=Path(j['directory'])/(j['name']+'.inp')
            _require(prepared.is_file() and file_sha256(prepared)==j['input_sha256'],
                     'Prepared native input changed; preserving the attempt for review.')
        manifests.append(manifest)
    logger.info('Automatic orbital verification: %s, %s probes, shared job cap %s',
                root,sum(len(json.loads(m.read_text())['jobs']) for m in manifests),parameters['queue_job_limit'])
    # run_many reconciles live allocations and retains bounded retry histories.
    # A scheduler query outage is not a scientific failure or a retry attempt.
    while True:
        try:
            recon.run_many(manifests, max_jobs=int(parameters['queue_job_limit']), interval=20.,
                           keep_going=True, allocation_retries=POLICY['maximum_allocation_retries'])
            break
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            command = error.cmd
            if not command or Path(command[0]).name != 'squeue':
                raise
            logger.warning('Automatic orbital verification is waiting for Slurm: %s', error)
            time.sleep(20.)
    result = {}
    for name, manifest in zip(sources, manifests):
        report = recon.audit(manifest)
        _require(all(j['status'] == 'complete' and not j.get('validation_errors') for j in report['jobs']),
                 'Native verification panel is incomplete: '+str(manifest))
        result[name] = report['jobs']
    return result


def _xml(job, kind='natural'):
    item = job['orbital_xml'][kind]
    _require(file_sha256(item['path']) == item['sha256'], 'Native orbital export changed.')
    return evidence.read_orbitals(item['path'])


def _output(job):
    path = Path(job['directory'])/(job['name']+'.out')
    _require(file_sha256(path) == job['output_sha256'], 'Native output changed.')
    return path.read_text(errors='replace')


def scf_phase(root, c, p, far_geometry):
    sources = {}
    for name, mapping, mult in [('dimer',list(range(len(c['symbols']))),3),
                               *[(f'fragment{i}',m,2) for i,m in enumerate(c['mappings'])]]:
        s = _base(c,p,name,np.asarray(far_geometry)[mapping], [c['symbols'][i] for i in mapping],mult)
        s.update(orbital_character_only=True, anchors=[{'name':'far','geometry_angstrom':s['geometry_angstrom']}])
        sources[name] = s
    jobs = _run_sources(root,sources,p)
    seeds, observations = {}, []
    for index, basis in enumerate(c['bases']):
        d = next(j for j in jobs['dimer'] if j['basis']==basis)
        a = _xml(d,'scf'); fragments = [_xml(next(j for j in jobs[f'fragment{i}'] if j['basis']==basis),'scf')
                                      for i in range(2)]
        somo, populations = evidence.fragment_embedding(a,fragments,c['mappings'])
        indices = np.flatnonzero(abs(a.occupations-1.) < 1.e-5)
        _require(len(indices)==2,'Separated dimer ROHF must contain exactly two SOMOs.')
        eigen = np.linalg.svd(somo[:,indices], compute_uv=False)**2
        _require(min(eigen) >= POLICY['minimum_subspace_projection'], 'Dimer SOMOs do not span both fragment radicals.')
        counts = [sum(atomic_numbers[c['symbols'][i]] for i in m) for m in c['mappings']]
        _require(max(abs(np.array(populations)-counts)) <= POLICY['maximum_fragment_charge_error'],
                 'Separated ROHF does not describe the two neutral fragments.')
        directory = Path(d['directory']); inp=directory/(d['name']+'.inp'); out=directory/(d['name']+'.out')
        seed=dict(path=d['wfu'],sha256=d['wfu_sha256'],basis=basis,active_indices=(indices+1).tolist(),
                  source_input=str(inp),source_input_sha256=file_sha256(inp),
                  source_output=str(out),source_output_sha256=file_sha256(out))
        from kinbot.recon_scf_seed import validate
        validate(seed,c['electrons'],2,2)
        seeds[basis]=seed
        observations.append(dict(basis=basis,fragment_populations=populations,somo_projection=eigen.tolist()))
    atomic_json(root/'fragment_identity.json',dict(observations=observations))
    return seeds,jobs


def _state(mult, selection=None):
    m=dict(state_average=[dict(multiplicity=mult,roots=1,weights=[1.])],
           selected_roots={str(mult):1},canonical_orbitals=True)
    if selection is not None:m['orbital_selection']=selection
    return m


def panel_source(c,p,name,mult,anchors,seeds=None,full=False):
    s=_base(c,p,name,anchors[0]['geometry_angstrom'],mult=mult)
    space=c['full_space'] if full else dict(active_electrons=2,active_orbitals=2)
    model=_state(mult,space.get('orbital_selection'))
    s.update(probe_active_space=True,probe_method='casscf',scf_multiplicity=3,
        anchors=anchors, probe_definition=dict(active_electrons=space['active_electrons'],
        active_orbitals=space['active_orbitals'],state_model=model,
        chemical_note='All-atom valence comparator' if full else 'Two independently identified fragment SOMOs'),
        casscf_options=dict(max_iterations=100,solver='so_sci',energy_convergence=1.e-9))
    e,o=space['active_electrons'],space['active_orbitals']
    s['correlation_checks']=[dict(method=f'caspt2-f12({e},{o})',shift=c['shift'],ipea=c['ipea'])]
    if not full:
        s.update(ci_basis='csf',ci_print_threshold=1.e-8)
        s['casscf_options']['solver']='second_order'
        s['probe_definition']['orbital_source']='scf'
        s['initial_records']={}
        seen=set()
        for a in anchors:
            chain=a.get('orbital_chain')
            if not chain or chain not in seen:
                for i,basis in enumerate(c['bases']):s['initial_records'][f'b{i}__'+a['name']]=deepcopy(seeds[basis])
            if chain:seen.add(chain)
        s['correlation_checks'] += [dict(method='caspt2-f12(2,2)',shift=c['shift']+.1,ipea=c['ipea']),
                                     dict(method='mrci-f12+q(2,2)')]
        if mult==3:
            s['correlation_checks'].append(dict(method='uccsd(t)-f12b',multiplicity=3))
            s['reference_identity_check']=True
    return s


def _selection(job, mult):
    return dict(orbital_selection=dict(mode='molpro_record',path=job['wfu'],sha256=job['wfu_sha256'],
        record='2140.2',orbital_type='natural',align_active_orbitals=True), **_state(mult))


def relax(root,c,p,reaction,pilots):
    """Use the production constrained optimizer; keep all its native attempts."""
    from kinbot.vrc_tst_scan import VTS
    q=copy(reaction); q.maps=deepcopy(c['mappings']); q.equiv=[[c['bond'][0]],[c['bond'][1]]]
    q.irc_prod=SimpleNamespace(bondlist=list(_graph(c['symbols'],reaction.species.bond).edges()))
    pp=deepcopy(p); pp['vrc_tst_active_space']=dict(mode='selected',active_electrons=2,active_orbitals=2)
    driver=VTS(reaction.species,pp,None)
    job=next(j for j in pilots['pilot1'] if j['basis']==c['scan_basis'])
    model=_selection(job,1); geometry=move(c,c['geometry'],c['far'])
    result={}; receipts=[]
    for i,distance in enumerate(reversed(c['radial'])):
        geometry=move(c,geometry,distance)
        atoms,energy,graph,execution,model=driver._molpro_multireference_scan_task(q,geometry,distance,
            'auto_'+digest(c)[:16]+f'_relax_{i}',state_model_override=model)
        geometry=atoms.positions
        result[f'r{i}']=dict(geometry_angstrom=geometry.tolist(),distance_angstrom=distance)
        receipts.append(dict(graph=str(graph),execution_sha256=file_sha256(graph/'tasks/relaxed_scan/execution.json'),
                             final_sha256=file_sha256(graph/'tasks/relaxed_scan/final.xyz')))
    atomic_json(root/'relaxations.json',dict(points=result,receipts=receipts))
    return result


def anchors(c,points):
    """Relaxed radial path, angular neighbors, reverse path, independent starts."""
    forward=[]; angle_ref={}
    for name,point in points.items():
        forward.append(dict(name=name,orbital_chain='radial',geometry_angstrom=point['geometry_angstrom']))
    # Three rotations at each of the three closest separations, plus matching
    # far geometry. Rotation about the reactive atom keeps R exactly fixed.
    near=sorted(points,key=lambda k:points[k]['distance_angstrom'])[:3]
    a,b=c['bond']; moving=max(c['mappings'],key=len)
    pivot=a if a in moving else b
    for cut,(axis_id,angle) in enumerate(((0,60.),(1,45.),(2,90.))):
        for name in ['r0',*reversed(near)]:
            xyz=np.array(points[name]['geometry_angstrom']); u=xyz[b]-xyz[a];u/=np.linalg.norm(u)
            v=np.eye(3)[np.argmin(abs(u))]; v=np.cross(u,v);v/=np.linalg.norm(v)
            axis=(u,v,np.cross(u,v))[axis_id];t=np.deg2rad(angle); x=xyz[moving]-xyz[pivot]
            xyz[moving]=xyz[pivot]+x*np.cos(t)+np.cross(axis,x)*np.sin(t)+np.outer(x@axis,axis)*(1-np.cos(t))
            key=f'a{cut}_{name}'; points[key]=dict(geometry_angstrom=xyz.tolist(),distance_angstrom=points[name]['distance_angstrom'])
            angle_ref[key]=name
            forward.append(dict(name=key,orbital_chain=f'angle{cut}',geometry_angstrom=xyz.tolist()))
    ordered=[];pairs=[]
    for chain in dict.fromkeys(a['orbital_chain'] for a in forward):
        rows=[a for a in forward if a['orbital_chain']==chain];ordered.extend(rows)
        for a in reversed(rows):
            key='back_'+a['name'];points[key]=deepcopy(points[a['name']]);pairs.append((a['name'],key))
            ordered.append(dict(a,name=key))
        for a in (rows[0],rows[-1]):
            key='ind_'+a['name'];points[key]=deepcopy(points[a['name']]);pairs.append((a['name'],key))
            ordered.append(dict(a,name=key,orbital_chain=key))
    return ordered,angle_ref,pairs


def _energy(job,index):
    value=job['correlated_energies_hartree'][index]['energy_hartree']
    _require(np.isfinite(value),'Nonfinite native correlated energy.')
    return value


def _job_index(jobs):
    return {(j['basis'],j['anchor']):j for j in jobs}


def _files(jobs):
    files={}
    for group in jobs.values():
        for j in group:
            d=Path(j['directory']); receipt=d/'execution.json'; ex=json.loads(receipt.read_text())
            for path in [receipt,Path(ex['native_input']),d/(j['name']+'.inp'),d/(j['name']+'.out'),Path(j['wfu']),
                         *[Path(x['path']) for x in j['orbital_xml'].values()],
                         *[d/name for name in ex.get('reference_identity_matrices',{})]]:
                files[str(path)]=file_sha256(path)
    return files


def assess(c,points,angle_ref,pairs,jobs,scf_jobs):
    """Apply the fixed policy to native results; absence never means a pass."""
    from kinbot.orbital_validation import compare_potentials,assess_path
    from kinbot.molpro_reference_identity import report as reference_report
    indices={name:_job_index(group) for name,group in jobs.items()}
    checks=[];details={};nc=(c['electrons']-2)//2;active=list(range(nc,nc+2))
    def check(name,passed,value):
        checks.append(dict(test=name,passed=bool(passed),observed=value))
    from kinbot.molpro_diagnostics import parse_native_output
    for name,group in jobs.items():
        mult=int(name[-1]);space=c['full_space'] if name.startswith('full') else dict(active_electrons=2,active_orbitals=2)
        for job in group:
            text=_output(job)
            parsed=parse_native_output(text)
            multi=[s for s in parsed['native_sections'] if s['program'].upper() in ('MULTI','MCSCF')]
            _require(len(multi)==1,'Verification requires one native MULTI calculation per probe.')
            native=multi[0];states=native.get('states',[])
            _require((native.get('active') or {}).get('total')==space['active_orbitals'] and
                     len(states)==1 and states[0].get('active_electrons')==space['active_electrons'] and
                     states[0].get('roots')==1 and states[0].get('spin','').lower()=={1:'singlet',3:'triplet'}[mult],
                     'Native CAS dimensions or electronic state differ from the automatic contract.')
    for basis in c['bases']:
        compact={};expanded={};alternate={};independent={}
        for point in points:
            geometry=[[s,*xyz] for s,xyz in zip(c['symbols'],points[point]['geometry_angstrom'])]
            for collection,index,model in ((compact,0,'compact'),(expanded,0,'full'),
                                            (alternate,1,'compact'),(independent,2,'compact')):
                collection[point]=dict(geometry_angstrom=geometry,energies_hartree={str(mult):
                    _energy(indices[f'{model}{mult}'][(basis,point)],index) for mult in (1,3)})
        pp={k:dict(v,geometry_angstrom=compact[k]['geometry_angstrom']) for k,v in points.items()}
        for label,comparison in [('full_valence',expanded),('shift',alternate),('mrci',independent)]:
            result=compare_potentials(pp,compact,comparison,asymptote='r0',angular_reference=angle_ref,
                                      tolerance_kcal_mol=POLICY['angular_tolerance_kcal_mol'])
            error=max(v['maximum_absolute_difference_kcal_mol'] for v in result['raw_comparisons'])
            check(basis+'/'+label,error<=POLICY['potential_tolerance_kcal_mol'] and
                  result['status']=='angular_test_passed',error)
            details[basis+'/'+label]=result
        path=assess_path(pp,compact,path=[k for k in points if re.fullmatch(r'r\d+',k)],asymptote='r0',
                         return_pairs=pairs,return_tolerance_kcal_mol=POLICY['return_tolerance_kcal_mol'])
        check(basis+'/return_energy',path['status']=='return_test_passed',path);details[basis+'/path']=path
        expanded_path=assess_path(pp,expanded,path=[k for k in points if re.fullmatch(r'r\d+',k)],asymptote='r0',
            return_pairs=pairs,return_tolerance_kcal_mol=POLICY['return_tolerance_kcal_mol'])
        check(basis+'/reference_return_energy',expanded_path['status']=='return_test_passed',expanded_path)
        far_compact=[compact['r0']['energies_hartree'][str(m)] for m in (1,3)]
        far_full=[expanded['r0']['energies_hartree'][str(m)] for m in (1,3)]
        check(basis+'/spin_degeneracy',max(abs(np.diff(far_compact)[0]),abs(np.diff(far_full)[0]))*AUtoKCAL
              <=POLICY['asymptote_tolerance_kcal_mol'],dict(compact=far_compact,full=far_full))
        for mult in (1,3):
            nearfar='r1'
            check(basis+f'/plateau{mult}',abs(compact[nearfar]['energies_hartree'][str(mult)]-
                  compact['r0']['energies_hartree'][str(mult)])*AUtoKCAL<=POLICY['asymptote_tolerance_kcal_mol'],
                  (compact[nearfar]['energies_hartree'][str(mult)]-compact['r0']['energies_hartree'][str(mult)])*AUtoKCAL)
            for point in points:
                job=indices[f'compact{mult}'][(basis,point)];nat=_xml(job)
                ci=evidence.natural_ci(_output(job),mult)
                ci_density=(2*evidence.singlet_matrix(ci)@evidence.singlet_matrix(ci).T if mult==1 else np.eye(2))
                _require(np.max(abs(ci_density-np.diag(nat.occupations[active])))<=2.e-5,
                         'Transformed CI and saved natural-orbital density disagree.')
                fulljob=indices[f'full{mult}'][(basis,point)];full=_xml(fulljob)
                closed=(c['electrons']-c['full_space']['active_electrons'])//2
                fv=list(range(closed,closed+c['full_space']['active_orbitals']))
                proj=evidence.projection(full,nat,fv,active)
                check(basis+f'/{point}/contains{mult}',min(proj)>=POLICY['minimum_subspace_projection'],proj)
                # Fractional density outside the candidate subspace, in the
                # reference NO basis. Closed shells give n(2-n)=0 exactly.
                t=evidence.transform(full,nat)[active,:]
                leakage=float(np.sum(full.occupations*(2-full.occupations)*(1-np.sum(t*t,axis=0))))
                check(basis+f'/{point}/extra_correlation{mult}',leakage<=POLICY['maximum_extra_fractionality'],leakage)
                if mult==3:
                    identity=reference_report(job['directory'],job['name'],job['probe_active_space'],c['electrons'],3)
                    t1=re.findall(r'T1 diagnostic:\s*([0-9.EeDd+-]+)',_output(job))
                    _require(len(t1)==1,'Expected one native CC T1 diagnostic.')
                    t1=float(t1[0].replace('D','E'))
                    ok=(min(identity['rohf_somo_active_projection_eigenvalues'])>=POLICY['minimum_subspace_projection']
                        and identity['closest_density_root']==1 and identity['reference_states'][0]['density_distance']
                        <=POLICY['maximum_reference_density_distance'] and t1<=POLICY['maximum_t1'])
                    check(basis+f'/{point}/cc_reference',ok,dict(t1=t1,identity=identity))
                if point=='r0':
                    fragment_rows=[_xml(next(j for j in scf_jobs[f'fragment{i}'] if j['basis']==basis),'scf') for i in range(2)]
                    somo,pop=evidence.fragment_embedding(nat,fragment_rows,c['mappings'])
                    projection=np.linalg.svd(somo[:,active],compute_uv=False)**2
                    check(basis+f'/fragment_somo{mult}',min(projection)>=POLICY['minimum_subspace_projection'],projection.tolist())
                    target=[sum(atomic_numbers[c['symbols'][i]] for i in m) for m in c['mappings']]
                    check(basis+f'/fragment_charge{mult}',max(abs(np.array(pop)-target))<=POLICY['maximum_fragment_charge_error'],pop)
                    if mult==1:
                        ionic=evidence.ionic_probability(ci,somo[:,active])
                        check(basis+'/covalent_singlet',ionic<=POLICY['maximum_ionic_probability'],ionic)
                    monomer_energy=sum(next(j for j in scf_jobs[f'fragment{i}'] if j['basis']==basis)['energy_hartree']
                                        for i in range(2))
                    error=(job['energy_hartree']-monomer_energy)*AUtoKCAL
                    check(basis+f'/matched_fragment_energy{mult}',abs(error)<=POLICY['asymptote_tolerance_kcal_mol'],error)
            for a,b in pairs:
                ja,jb=(indices[f'compact{mult}'][(basis,k)] for k in (a,b))
                complete_transform=evidence.transform(_xml(ja),_xml(jb))
                t=complete_transform[np.ix_(active,active)]
                proj=np.linalg.svd(t,compute_uv=False)**2
                overlap=evidence.wavefunction_overlap(evidence.natural_ci(_output(ja),mult),
                    evidence.natural_ci(_output(jb),mult),complete_transform,nc,mult)
                check(basis+f'/{a}/{b}/state{mult}',min(proj)>=POLICY['minimum_return_projection'] and
                      overlap>=POLICY['minimum_ci_overlap'],dict(projection=proj.tolist(),ci_overlap=overlap))
        # CC must agree with an independent MR reference on the triplet
        # interaction potential, rather than only at a constant energy offset.
        for point in points:
            cc=(_energy(indices['compact3'][(basis,point)],3)-_energy(indices['compact3'][(basis,'r0')],3))*AUtoKCAL
            mr=(independent[point]['energies_hartree']['3']-independent['r0']['energies_hartree']['3'])*AUtoKCAL
            check(basis+'/'+point+'/cc_mr_potential',abs(cc-mr)<=POLICY['potential_tolerance_kcal_mol'],cc-mr)
    # Use the actual configured CBS formula as well as checking each basis.
    # A small per-basis error can be amplified by extrapolation.
    from kinbot.anl.extrapolation import two_point_cbs
    from kinbot.vrc_correction import _basis_cardinal
    gap=c['spin']['gap_levels'][0];lo,hi=gap['basis']
    def cbs(model,mult,point,index=0,level=gap):
        return two_point_cbs(_energy(indices[f'{model}{mult}'][(lo,point)],index),
            _energy(indices[f'{model}{mult}'][(hi,point)],index),upper_cardinal=_basis_cardinal(hi),
            power=level.get('extrapolation_power',3.7))
    cbs_observations={}
    for label,model,index in [('compact','compact',0),('full','full',0),('shift','compact',1),('mrci','compact',2)]:
        cbs_observations[label]={point:dict(geometry_angstrom=[[s,*v] for s,v in zip(c['symbols'],p['geometry_angstrom'])],
            energies_hartree={str(m):cbs(model,m,point,index) for m in (1,3)}) for point,p in points.items()}
    pp={k:dict(v,geometry_angstrom=cbs_observations['compact'][k]['geometry_angstrom']) for k,v in points.items()}
    for label in ('full','shift','mrci'):
        result=compare_potentials(pp,cbs_observations['compact'],cbs_observations[label],asymptote='r0',
            angular_reference=angle_ref,tolerance_kcal_mol=POLICY['angular_tolerance_kcal_mol'])
        error=max(v['maximum_absolute_difference_kcal_mol'] for v in result['raw_comparisons'])
        check('CBS/'+label,error<=POLICY['potential_tolerance_kcal_mol'] and result['status']=='angular_test_passed',error)
        details['CBS/'+label]=result
    for point in points:
        reference=c['spin']['reference']
        cc=(cbs('compact',3,point,3,reference)-cbs('compact',3,'r0',3,reference))*AUtoKCAL
        mr=(cbs('compact',3,point,2)-cbs('compact',3,'r0',2))*AUtoKCAL
        check('CBS/'+point+'/cc_mr',abs(cc-mr)<=POLICY['potential_tolerance_kcal_mol'],cc-mr)
    effective={point:dict(geometry_angstrom=cbs_observations['full'][point]['geometry_angstrom'],
        energies_hartree={'1':cbs('compact',3,point,3,c['spin']['reference'])+
                         cbs('compact',1,point)-cbs('compact',3,point)}) for point in points}
    reference={point:dict(geometry_angstrom=row['geometry_angstrom'],energies_hartree={'1':row['energies_hartree']['1']})
               for point,row in cbs_observations['full'].items()}
    final=compare_potentials(pp,effective,reference,asymptote='r0',angular_reference=angle_ref,
                            tolerance_kcal_mol=POLICY['angular_tolerance_kcal_mol'])
    maximum=max(v['maximum_absolute_difference_kcal_mol'] for v in final['raw_comparisons'])
    check('CBS/reconstructed_singlet',maximum<=POLICY['potential_tolerance_kcal_mol'] and
          final['status']=='angular_test_passed',maximum)
    details['CBS/reconstructed_singlet']=final
    sample=c['sampling']
    from kinbot.anl.extrapolation import two_point_cbs
    sampling={}
    for point in points:
        if isinstance(sample['basis'],str):
            energy=_energy(indices['compact1'][(sample['basis'],point)],0)
        else:
            a,b=sample['basis'];energy=two_point_cbs(_energy(indices['compact1'][(a,point)],0),
                _energy(indices['compact1'][(b,point)],0),upper_cardinal=_basis_cardinal(b),
                power=sample.get('extrapolation_power',3.7))
        sampling[point]=dict(geometry_angstrom=reference[point]['geometry_angstrom'],energies_hartree={'1':energy})
    angular=compare_potentials(pp,sampling,reference,asymptote='r0',angular_reference=angle_ref,
                              tolerance_kcal_mol=POLICY['angular_tolerance_kcal_mol'])
    check('sampling/corrected_angular_residual',angular['status']=='angular_test_passed',angular['angular_residuals'])
    details['sampling/corrected_angular_residual']=angular
    return dict(schema=1,policy=c['policy'],contract_sha256=digest(c),checks=checks,details=details,
                status='accepted' if checks and all(v['passed'] for v in checks) else 'needs_review',
                scope='Finite radial/angular, relaxed-geometry, return and independent-start electronic-model verification; '
                      'production sampling, kinetics and numerical convergence retain their own gates.',
                files={**_files(jobs),**_files(scf_jobs)})


def surface(c,jobs,verification):
    children={}
    for mult in (1,3):
        children[str(mult)]=dict(schema=2,state_average_policy='reviewed',
            molecule=dict(symbols=c['symbols'],charge=0),active_space=dict(active_electrons=2,active_orbitals=2),
            basis_models={basis:_selection(next(j for j in jobs[f'compact{mult}'] if j['basis']==basis and j['anchor']=='r0'),mult)
                          for basis in c['bases']})
    parameters=dict(vrc_tst_active_space=dict(mode='selected',active_electrons=2,active_orbitals=2,
                    source='Native automatic radical-pair policy v1: '+str(verification)),
                    vrc_tst_state_model=dict(schema=3,spin_models=children),vrc_tst_state_model_file='')
    from kinbot.vrc_multistate import validate_surface_model
    return validate_surface_model(dict(schema=1,reviewed=True,
        review_note='Passed bounded native radical-pair verification; evidence: '+str(verification),
        channels=[dict(bond=[i+1 for i in c['bond']],reaction=c['reaction'],combination='single',
        bottlenecks=[dict(name='single',coordinate='reactive',distances_angstrom=c['distances'],parameters=parameters,
                         states=[dict(multiplicity=1,root=1,weight=1.,energy_offset_cm_inverse=0.,
                                      energy_reference_note='Nondegenerate singlet; triplet is the correction reference, not an additional reactive surface.')])])]))


def _cached(root,c):
    result=root/'verification.json';model=root/'surface.json'
    if not result.is_file():return None
    saved=json.loads(result.read_text())
    _require(saved.get('contract_sha256')==digest(c),'Verification contract changed.')
    if saved['status']!='accepted':
        raise ReviewRequired('Native automatic verification requires review: '+str(result))
    _require(saved.get('files') and saved.get('selected_surface'),'Accepted surface evidence is incomplete.')
    for path,sha in saved['files'].items():
        _require(Path(path).is_file() and file_sha256(path)==sha,'Accepted native evidence changed: '+path)
    if not model.is_file():
        _publish(root,saved)
    loaded=json.loads(model.read_text())
    without_binding=deepcopy(loaded)
    for channel in without_binding['channels']:
        for stage in channel['bottlenecks']:
            stage['parameters'].pop('vrc_tst_reference_review',None)
    _require(without_binding==saved['selected_surface'],'Accepted surface changed after verification.')
    return loaded


def _publish(root,report):
    """Finish a recoverable publication after the durable scientific decision."""
    model=deepcopy(report['selected_surface']);reference=root/'reference_review.json'
    _write_immutable(reference,dict(kind='automatic_radical_pair_v1',verification=dict(
        path=str(root/'verification.json'),sha256=file_sha256(root/'verification.json'))))
    model['channels'][0]['bottlenecks'][0]['parameters']['vrc_tst_reference_review']=dict(
        path=str(reference),sha256=file_sha256(reference))
    _write_immutable(root/'surface.json',model)
    return model


def validate_receipt(parameters,saved):
    """Native automatic evidence is a separate, fully bound reference review."""
    from kinbot.vrc_reference_review import contract
    witness=saved.get('verification',{})
    path=Path(witness.get('path',''))
    _require(path.is_file() and file_sha256(path)==witness.get('sha256'),'Automatic verification receipt changed.')
    report=json.loads(path.read_text())
    _require(report.get('status')=='accepted' and report.get('policy')==POLICY and
             report.get('checks') and all(x.get('passed') is True for x in report['checks']),
             'Automatic native evidence did not pass the current policy.')
    _require(report.get('production_contract')==json.loads(json.dumps(contract(parameters))),
             'Automatic verification belongs to another production model or electronic contract.')
    _require(report.get('files'),'Automatic verification has no native artifacts.')
    for name,sha in report['files'].items():
        _require(Path(name).is_file() and file_sha256(name)==sha,'Automatic native evidence changed: '+name)


def ensure(parameters,reaction):
    """Called by the owning VRC driver before production fragment/scan jobs."""
    from kinbot.vrc_multistate import SurfaceModelRequired
    try:
        c=context(parameters,reaction)
    except ReviewRequired as error:
        raise SurfaceModelRequired(str(error)) from error
    root=Path('vrctst/auto_radical_pair')/digest(c)
    root=root.resolve();root.mkdir(parents=True,exist_ok=True)
    try:
        with recon._manifest_lock(root):
            _write_immutable(root/'contract.json',c)
            model=_cached(root,c)
            if model is None:
                atomic_json(root/'status.json',dict(phase='fragment_scf',contract_sha256=digest(c)))
                far=move(c,c['geometry'],c['far'])
                seeds,initial_scf=scf_phase(root/'initial_scf',c,parameters,far)
                pilot_anchor=[dict(name='far',geometry_angstrom=far.tolist(),orbital_chain='pilot')]
                pilots=_run_sources(root/'pilots',{f'pilot{m}':panel_source(c,parameters,f'pilot{m}',m,pilot_anchor,seeds) for m in (1,3)},parameters)
                atomic_json(root/'status.json',dict(phase='relaxation',contract_sha256=digest(c)))
                points=relax(root,c,parameters,reaction,pilots)
                seeds,scf=scf_phase(root/'relaxed_scf',c,parameters,points['r0']['geometry_angstrom'])
                panel,angles,pairs=anchors(c,points)
                _write_immutable(root/'coverage.json',dict(points=points,angles=angles,return_pairs=pairs))
                sources={f'{name}{m}':panel_source(c,parameters,f'{name}{m}',m,panel,seeds,full=name=='full')
                         for name in ('compact','full') for m in (1,3)}
                atomic_json(root/'status.json',dict(phase='verification_panel',contract_sha256=digest(c)))
                jobs=_run_sources(root/'panel',sources,parameters)
                report=assess(c,points,angles,pairs,jobs,scf)
                report['files'].update(_files(initial_scf))
                report['files'].update(_files(pilots))
                for name in ('contract.json','coverage.json','relaxations.json'):
                    report['files'][str(root/name)]=file_sha256(root/name)
                relaxation=json.loads((root/'relaxations.json').read_text())
                for r in relaxation['receipts']:
                    d=Path(r['graph'])/'tasks/relaxed_scan'
                    report['files'][str(d/'execution.json')]=r['execution_sha256']
                    report['files'][str(d/'final.xyz')]=r['final_sha256']
                if report['status']=='accepted':
                    model=surface(c,jobs,root/'verification.json')
                    from kinbot.vrc_reference_review import contract
                    stage=model['channels'][0]['bottlenecks'][0]
                    from kinbot.vrc_multistate import state_parameters
                    leaf=state_parameters(parameters,stage,stage['states'][0],[2,2])
                    report['production_contract']=json.loads(json.dumps(contract(leaf)))
                    report['selected_surface']=deepcopy(model)
                atomic_json(root/'verification.json',report)
                _require(report['status']=='accepted','Native verification failed; inspect '+str(root/'verification.json'))
                model=_publish(root,report)
            from kinbot.vrc_reference_review import validate
            from kinbot.vrc_multistate import state_parameters
            stage=model['channels'][0]['bottlenecks'][0]
            validate(state_parameters(parameters,stage,stage['states'][0],[2,2]),root)
            atomic_json(root/'status.json',dict(phase='accepted',contract_sha256=digest(c)))
            return model['channels'][0],dict(source=str(root/'surface.json'),sha256=file_sha256(root/'surface.json'),
                                             verification=str(root/'verification.json'))
    except (ValueError,recon.ProbeFailure,KeyError) as error:
        if not (root/'verification.json').exists():
            atomic_json(root/'verification.json',dict(status='needs_review',contract_sha256=digest(c),
                        policy=c['policy'],reason=str(error)))
        atomic_json(root/'status.json',dict(phase='needs_review',reason=str(error),contract_sha256=digest(c)))
        raise SurfaceModelRequired(str(error)+'; evidence: '+str(root)) from error
    except RuntimeError as error:
        if str(error).startswith('Another reconnaissance controller'):
            raise SurfaceModelRequired('Automatic verification is owned by a live controller: '+str(root)) from error
        atomic_json(root/'status.json',dict(phase='reconciliation_required',reason=str(error),contract_sha256=digest(c)))
        raise SurfaceModelRequired('Automatic verification needs scheduler/attempt reconciliation: '+str(root)+'; '+str(error)) from error
