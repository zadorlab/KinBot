"""Reconcile existing choices without changing or accepting any surface model."""
from collections import Counter
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]


def read(path):
    return json.loads(path.read_text())


def proof(path, relative_to):
    return dict(path=str(path.relative_to(relative_to)),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def installed(channel):
    result = []
    for bottleneck in channel['bottlenecks']:
        p = bottleneck['parameters']
        model = p['vrc_tst_state_model']
        result.append(dict(name=bottleneck['name'], active_space=p['vrc_tst_active_space'],
                           state_model=model, reactive_states=bottleneck['states'],
                           coordinate=bottleneck.get('coordinate')))
    return result


def build(export):
    inventory = export/'runs/ch3s_ch3ss_orbital_reviews_20261027/requests.20261027T033613Z.json'
    captured = read(inventory)
    if captured['issues']:
        raise ValueError(captured['issues'])
    recipes = {}
    old = HERE.parent/'r7_compact_union_20261025'
    for entry in read(old/'publication_plan.json')['selections']:
        path = old/entry['recipe']; data = read(path)
        name = data['request_match']['reaction']
        if name in recipes:
            raise ValueError('Duplicate selection: '+name)
        recipes[name] = (path, data, 'prepared_previously', entry['source'])
    for path in sorted((HERE/'selections').glob('*.json')):
        data = read(path); name = data['request_match']['reaction']
        if name in recipes:
            raise ValueError('Do not replace a prior selection: '+name)
        recipes[name] = (path, data, 'prepared_from_detached_panel', data.get('source_model'))
    rows = []
    for item in captured['requests']:
        q = item['request']; name = q['reaction']
        row = dict(case=item['case'], reaction=name, request=q,
                   request_sha256=item['sha256'], installed_channels=item['installed_channels'])
        if q['status'] == 'model_accepted':
            if not item['installed_channels'] or name in recipes:
                raise ValueError('Inconsistent installed selection: '+name)
            row.update(disposition='carry_forward_installed',
                       models=[m for c in item['installed_channels'] for m in installed(c)])
        elif name in recipes:
            path, data, status, source = recipes.pop(name)
            for key, value in data['request_match'].items():
                if q.get(key) != value:
                    raise ValueError(f'{name}: changed request field {key}')
            row.update(disposition=status, selection=proof(path, REPO),
                       source_review=source, models=[dict(active_space=data['active_space'],
                       state_model=data['state_model'], reactive_states=data['states'],
                       basis_seeds=data['basis_seeds'], rationale=data['review_note'])])
        else:
            row.update(disposition='unresolved', models=[])
        rows.append(row)
    if recipes:
        raise ValueError('Prepared choices absent from inventory: '+str(sorted(recipes)))
    inputs = []
    for case, parent in (('ch3s', '470831050000000000002'),
                         ('ch3ss', '791471490990000000002')):
        path = export/'production'/case/'input.json'; data = read(path)
        inputs.append(proof(path, export))
        for c in data['vrc_tst_surface_model']['channels']:
            name = c.get('reaction', parent+'_hom_sci_'+'_'.join(map(str,c['bond'])))
            rows.append(dict(case=case, reaction=name, disposition='carry_forward_initial',
                             installed_channels=[c], models=installed(c),
                             evidence=f'examples/anl/{case}_anl1_f12/orbital_reviews/ACTIVE_SPACE.md'))
    names = [r['reaction'] for r in rows]
    if len(set(names)) != len(names) or len(rows) != 34:
        raise ValueError('Known 34-channel inventory is not reconciled.')
    return dict(schema=1, status='incomplete_not_launch_ready',
                scope='Choice inventory, not a publication or a restart-readiness certificate.',
                source_inventory=proof(inventory, export), case_inputs=inputs,
                counts=dict(Counter(r['disposition'] for r in rows)),
                channels=sorted(rows,key=lambda r:(r['case'],r['reaction'])))


def render(data):
    lines = ['# Known channel choices for the r8 handoff', '',
             'This inventory preserves the installed and previously prepared decisions. '
             'It does not approve unresolved cases, revive a rejected reaction, or override L1 pruning. '
             'The full JSON retains the exact roots, weights, checkpoint hashes and original rationale.', '',
             '| Reaction | Disposition | CAS | Orbital averaging | Reactive roots | Evidence |',
             '| --- | --- | --- | --- | --- | --- |']
    for row in data['channels']:
        cas=[]; averages=[]; flux=[]
        for model in row['models']:
            s=model['active_space'];cas.append(f"({s['active_electrons']},{s['active_orbitals']})")
            state=model['state_model']
            basis=state.get('basis_models')
            states=list(basis.values()) if basis else [state]
            averages.append('/'.join(sorted(set(', '.join(
                f"{x['roots']}{ {2:'D',4:'Q'}.get(x['multiplicity'],str(x['multiplicity']))}"
                for x in b['state_average']) for b in states))))
            flux.append(', '.join(f"{ {2:'D',4:'Q'}.get(x['multiplicity'],str(x['multiplicity']))}{x['root']}"
                                  for x in model['reactive_states']))
        evidence=row.get('selection',{}).get('path',row.get('evidence','Captured installed model / request inventory'))
        lines.append(f"| `{row['reaction']}` | {row['disposition']} | {'; '.join(cas) or 'pending'} | "
                     f"{'; '.join(averages)} | {'; '.join(flux)} | {evidence} |")
    lines.extend(['', '## Retained explanations', ''])
    for row in data['channels']:
        if not row['models']:continue
        lines.extend([f"### {row['reaction']}", ''])
        for m in row['models']:
            lines.extend([m.get('rationale',m['active_space'].get('source','See linked original review.')), ''])
    return '\n'.join(lines)+'\n'


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('export',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();data=build(a.export.resolve(strict=True))
    a.output.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    a.output.with_suffix('.md').write_text(render(data))
    print(json.dumps(data['counts'],sort_keys=True))
