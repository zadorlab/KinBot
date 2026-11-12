"""Read native campaign records and print a portable text review with hashes."""
import argparse
import importlib.util
import json
from pathlib import Path

from kinbot.execution import file_sha256
from kinbot.molpro_diagnostics import parse_native_output
from kinbot.radical_pair_evidence import read_orbitals, projection

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('completion_campaign',HERE/'campaign.py')
campaign=importlib.util.module_from_spec(spec);spec.loader.exec_module(campaign)


def emit(kind,value):
    print(kind+' '+json.dumps(value,sort_keys=True,allow_nan=False),flush=True)


def report(root):
    root=Path(root).resolve(strict=True)
    emit('CAMPAIGN',dict(directory=str(root),plan_sha256=file_sha256(root/'plan.json'),production_model_accepted=False))
    receipt=campaign.read(root/'campaign.json')
    complete=campaign.read(root/'completion.json') if (root/'completion.json').exists() else {}
    emit('COMPLETION',complete)
    manifests=sorted(set(receipt['manifests']+complete.get('manifests',[])+[str(p) for p in root.glob('*/h*/native/reconnaissance.json')]))
    selected_screen_titles=set()
    bases=campaign.read(root/'plan.json')['bases']
    for path in root.glob('h*_followup_trials.json'):
        prefix=path.name.split('_')[0]
        for name in campaign.read(path)['selected']:
            selected_screen_titles.update(f'{prefix}_{name}_{b.replace("-","_")}' for b in bases)
    indexed={};counts={};issues=[]
    for filename in manifests:
        path=Path(filename);data=campaign.read(path);cfg=campaign.read(path.parents[1]/'input.json')
        title=cfg['title']
        emit('CONFIGURATION',dict(title=title,input_sha256=file_sha256(path.parents[1]/'input.json'),
             definition=cfg.get('probe_definition'),initial_records=cfg.get('initial_records'),anchors=cfg['anchors']))
        for job in data['jobs']:
            status=job.get('status','waiting');counts[status]=counts.get(status,0)+1
            ident=[title,job['basis'],job['anchor']]
            emit('POINT_STATUS',dict(identity=ident,status=status,job_id=job.get('job_id'),validation_errors=job.get('validation_errors')))
            if status!='complete':continue
            directory=Path(job['directory']);out=directory/(job['name']+'.out');inp=out.with_suffix('.inp')
            changed=[]
            for p,key in ((out,'output_sha256'),(inp,'input_sha256'),(Path(job['wfu']),'wfu_sha256')):
                if not p.is_file() or file_sha256(p)!=job.get(key):changed.append(str(p))
            if changed:
                issues.append(dict(identity=ident,changed_native_files=changed));continue
            parsed=parse_native_output(out.read_text(errors='replace'),inp.read_text())
            row=dict(identity=ident,geometry=job['geometry_angstrom'],roots=job['root_states'],
                     active_space=job.get('probe_active_space'),correlation=job.get('correlated_energies_hartree'),
                     input_sha256=job['input_sha256'],output_sha256=job['output_sha256'],wfu=job['wfu'],wfu_sha256=job['wfu_sha256'],
                     normal_termination=parsed['normal_termination_seen'],native_issues=parsed['issues'],native_errors=parsed['native_error_lines'])
            emit('NATIVE_POINT',row)
            indexed[tuple(ident)]=job
            # Full-path, reference and binding jobs carry native active character and CI evidence.
            # Endpoint screens retain all energies; their chosen representatives occur in these records.
            if path.is_relative_to(root/'screen') and title not in selected_screen_titles:continue
            sections=[]
            for section in parsed['native_sections']:
                entry={k:v for k,v in section.items() if k not in ('orbital_tables','evidence')}
                closed=section.get('closed',{}).get('total',0);active=section.get('active',{}).get('total',0)
                entry['orbital_tables']=[dict(heading=t['heading'],line=t['line'],orbitals=[o for o in t['orbitals']
                    if max(1,closed-2)<=int(o['index'].split('.')[0])<=closed+active+3]) for t in section.get('orbital_tables',[])]
                entry['evidence']=section.get('evidence',[])
                sections.append(entry)
            emit('NATIVE_CHARACTER',dict(identity=ident,sections=sections))
    # Reused screen points are included under their full-path identity for return checks.
    for path in sorted(root.glob('h*_reused.json')):
        name=path.name.removesuffix('_reused.json')
        for basis in campaign.read(root/'plan.json')['bases']:
            if name.endswith('_'+basis):
                title=name[:-(len(basis)+1)]+'_full_'+basis.replace('-','_')
                for anchor,job in campaign.read(path).items():indexed[(title,basis,anchor)]=job
    for key,right in indexed.items():
        title,basis,anchor=key
        original={'back_near':'near','back6':'r6','back_far':'far','back_in':'in'}.get(anchor)
        if original is None or (title,basis,original) not in indexed:continue
        left=indexed[(title,basis,original)]
        try:
            a=read_orbitals(Path(left['directory'])/(left['name']+'.natural.orbitals.xml'))
            b=read_orbitals(Path(right['directory'])/(right['name']+'.natural.orbitals.xml'))
            space=right['probe_active_space'];n=sum(campaign.atomic_numbers[s] for s in right['symbols'])
            closed=(n-space['active_electrons'])//2;indices=list(range(closed,closed+space['active_orbitals']))
            values=projection(a,b,indices,indices)
            emit('RETURN_CHECK',dict(identity=key,reference=original,active_subspace_squared_singular_values=values,
                  difference_kcal_mol=campaign.KCAL*(campaign.energy(right)-campaign.energy(left)),
                  left_xml_sha256=file_sha256(Path(left['directory'])/(left['name']+'.natural.orbitals.xml')),
                  right_xml_sha256=file_sha256(Path(right['directory'])/(right['name']+'.natural.orbitals.xml'))))
        except (ValueError,OSError,KeyError) as error:
            emit('RETURN_CHECK_UNAVAILABLE',dict(identity=key,reason=str(error)))
    for path in sorted(root.glob('h*_followup_trials.json')):
        emit('TRIAL_ALLOCATION',campaign.read(path))
    for path in sorted(root.glob('h*_partitions.json')):
        data=campaign.read(path)
        emit('HYDROGEN_IDENTITY',dict(file=path.name,sources={b:s['partition'] for b,s in data['sources'].items()}))
    emit('SUMMARY',dict(counts=counts,integrity_issues=issues,production_model_accepted=False))
    return 2 if issues else 0


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('directory',type=Path)
    raise SystemExit(report(parser.parse_args().directory))
