"""Reconcile the all-channel native export without changing calculation files."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re

from kinbot.molpro_diagnostics import parse_native_output

HERE = Path(__file__).resolve().parent
ARCHIVE = '15b6ef48db20e9aa97b901a3e15759003728bc3bea585865cffcef4c90f3a256'


def read(path):
    return json.loads(path.read_text())


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def audit(root, manifest_path):
    receipt = read(root/'VERIFIED.json')
    if receipt['archive_sha256'] != ARCHIVE or receipt['issues']:
        raise ValueError('Require the verified all-channel export.')
    manifest = read(manifest_path)
    index = {f['member']: f for f in manifest['files']}

    def evidence(path):
        member = str(path.relative_to(root))
        expected = index[member]
        if digest(path) != expected['sha256']:
            raise ValueError('Changed native evidence: '+member)
        return dict(path=member, sha256=expected['sha256'])

    def native(out):
        inp = out.with_suffix('.inp')
        if not inp.exists():
            return dict(output=evidence(out), missing_input=True)
        text = out.read_text(errors='replace')
        deck = inp.read_text()
        parsed = parse_native_output(text, deck)
        sections = [s for s in parsed['native_sections'] if s['program']=='MULTI']
        first = sections[0] if sections else {}
        basis = re.search(r'^basis\s*=\s*([^\n]+)', deck, re.M | re.I)
        return dict(output=evidence(out), input=evidence(inp),
            basis=basis.group(1) if basis else None,
            normal=parsed['normal_termination_seen'], errors=parsed['native_error_lines'],
            issues=parsed['issues'], requested=parsed['requested_casscf'],
            correlation=parsed.get('requested_correlation_calls'),
            geometries=parsed.get('native_geometries'),
            states=first.get('state_energies'), ci=first.get('ci_tables'),
            orbitals=first.get('orbital_tables', [])[-1:] or [],
            multiple_geometry_sections=len(sections),
            energies=re.findall(r'RS2-F12/[^\n]*energy=\s*([-\d.]+)',text))

    ledger = read(HERE/'prior_selection_ledger.json')
    rows=[]
    for channel in ledger['channels']:
        name=channel['reaction']; parent=name.split('_',1)[0]
        folder=root/'production'/channel['case']/parent/'vrctst'
        row=dict(reaction=name, disposition=channel['disposition'],
                 r8_selection=channel['r8_selection'], points=[], graphs=[])
        for out in sorted((folder/'molpro').glob(name+'__*.out')):
            row['points'].append(native(out))
        for p in sorted((folder/'molpro_scan').glob(name+'__*/state.json')):
            state=read(p); workflow=read(p.parent/'workflow.json')
            task=p.parent/'tasks/relaxed_scan'
            outputs=sorted(task.glob('relaxed_scan_step_*.out'))
            row['graphs'].append(dict(state=evidence(p),workflow=evidence(p.parent/'workflow.json'),
                distance=workflow['intent'].get('vrc_scan_distance_angstrom'),
                statuses={k:v.get('status') for k,v in state['tasks'].items()},
                last_native=native(outputs[-1]) if outputs else None))
        row['summary']=dict(native_single_points=len(row['points']),
            normal_single_points=sum(p.get('normal',False) and not p.get('errors') for p in row['points']),
            point_bases=dict(Counter(p.get('basis') for p in row['points'])),
            geometry_statuses=dict(Counter(g['statuses'].get('relaxed_scan') for g in row['graphs'])))
        rows.append(row)
        print(name, row['summary'], flush=True)
    return dict(schema=1, archive_sha256=ARCHIVE, verification=receipt,
                channels=rows, limitations=[
                    'Normal termination is not scientific acceptance.',
                    'Geometry-gradient outputs contain displaced sections: reported MULTI character is the first central section.',
                    'Complete graphs are historical records; no live scheduler was queried.',
                    'Existing radial paths do not substitute for matched expansion or finite angular tests.'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('manifest',type=Path)
    args=p.parse_args();result=audit(args.root,args.manifest)
    with gzip.open(HERE/'production_native_audit.json.gz','wt') as f:
        json.dump(result,f)
