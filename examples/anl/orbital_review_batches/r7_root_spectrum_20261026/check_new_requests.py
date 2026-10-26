"""Read-only validation of five supplied requests and one shared endpoint.

Refuse changed requests/endpoints. This shares orbital evidence only, never
merges reaction paths or accepts production models. Repeat before publication.
"""
import argparse
import json
from pathlib import Path

from kinbot import requested_orbital_recon as review
from kinbot.execution import file_sha256


def check(production, plan=None):
    production=Path(production).expanduser().resolve(strict=True)
    plan=Path(plan or Path(__file__).resolve().parent)
    evidence=review._read(plan/'new_requests_evidence.json')
    results={}
    for item in evidence['requests']:
        expected=item['request']
        recipe=next(review._read(p) for p in sorted(plan.glob('new_*.json'))
                    if review._read(p).get('request_match',{}).get('reaction')==expected['reaction'])
        path=production/'ch3ss'/expected['parent']/'vrctst/orbital_review'/expected['reaction']/'request.json'
        request=review._read(path)
        if request.get('status')!='needs_review' or any(request.get(k)!=v for k,v in recipe['request_match'].items()):
            raise ValueError('Live request differs from supplied evidence: '+str(path))
        proof=dict(request=str(path),request_sha256=file_sha256(path))
        if item['geometry']:
            xyz,anchors,source=review.product_complex_geometry(path,request,recipe['geometry_source'])
            if xyz!=item['geometry']['positions']:
                raise ValueError('Live product endpoint differs from supplied evidence: '+expected['reaction'])
            proof.update(geometry_source=source,anchors=anchors)
        results[expected['reaction']]=proof
    for alias in review._read(plan/'aliases.json')['aliases']:
        p=review._read(plan/alias['primary_recipe'])['request_match']['reaction']
        a=review._read(plan/alias['alias_recipe'])['request_match']['reaction']
        left=results[p];right=results[a]
        for key in ('symbols','positions','fragments'):
            if left['geometry_source'][key]!=right['geometry_source'][key]:
                raise ValueError('Alias no longer has the same mapped geometry: '+a)
        if left['anchors']!=right['anchors']:
            raise ValueError('Alias no longer has the same scan anchors: '+a)
        right['orbital_screen_shared_with']=p
    return dict(status='matched',requests=results,
                note='Read-only request/geometry validation; no reaction equivalence or orbital selection implied.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('production',type=Path)
    args=parser.parse_args()
    print(json.dumps(check(args.production),indent=2))
