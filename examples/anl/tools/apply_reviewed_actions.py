#!/usr/bin/env python3
"""Apply an external review manifest through KinBot's live control APIs.

Preview by default. No sbatch, cancellation, controller replacement or ledger
reset occurs here. Run-specific identities belong in the supplied manifest.
"""
import argparse
import json
from pathlib import Path

from kinbot import live_control, rotdpy_campaign
from kinbot.execution import file_sha256


def inside(root, relative):
    path = (root / relative).resolve(strict=True)
    if not path.is_relative_to(root):
        raise ValueError(f'Path outside run: {relative}')
    return path


def prepare(root, manifest):
    if manifest.get('schema') != 1:
        raise ValueError('Unknown review manifest schema.')
    reviews, retries = [], []
    for item in manifest['reviews']:
        inp = inside(root, item['input'])
        if file_sha256(inp) != item['input_sha256']:
            raise ValueError(f'{inp}: reviewed input changed.')
        saved = json.loads(rotdpy_campaign.receipt(inp).read_text())
        current = rotdpy_campaign.validate(inp, saved)
        if current['correction']['source_sha256'] != item['correction_sha256']:
            raise ValueError(f'{inp}: reviewed correction changed.')
        correction = inside(root, current['correction']['source'])
        data = json.loads(correction.read_text())
        for components in data['native_components'].values():
            for component in components:
                for ext, key in (('.inp', 'input_sha256'), ('.out', 'output_sha256')):
                    native = inside(root, correction.parent/'molpro'/(component['job']+ext))
                    if file_sha256(native) != component[key]:
                        raise ValueError(f'{native}: reviewed native evidence changed.')
        reviews.append((inp, item['note'], current))
    for item in manifest['retries']:
        target = inside(root, item['target'])
        current = live_control.generation(target, 'task')
        if any(current.get(k) != v for k,v in item['expected'].items()):
            print('RETAIN changed generation', item['target'], current)
            continue
        if current['status'] not in ('failed', 'unavailable'):
            print('RETAIN current status', item['target'], current)
            continue
        pending = [r for _,r in live_control.requests(target.parent, kinds={'task'})
                   if r['target'] == str(target)]
        if pending:
            print('ALREADY QUEUED', item['target'], [r['id'] for r in pending])
            continue
        retries.append((target, item))
    return reviews, retries


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('manifest', type=Path)
    p.add_argument('--apply', action='store_true')
    args = p.parse_args()
    root = args.root.resolve(strict=True)
    reviews, retries = prepare(root, json.loads(args.manifest.read_text()))
    print(f'Reviewed actions: {len(reviews)} correction receipts, {len(retries)} task retry requests.')
    for inp, note, expected in reviews:
        if args.apply:
            rotdpy_campaign.accept(inp, note, expected=expected)
        print('ACCEPTED' if args.apply else 'REVIEW READY', inp.name)
    for target, item in retries:
        if args.apply:
            r = live_control.request(root, target, reason=item['reason'], expected=item['expected'])
            print('QUEUED FOR OWNER', r['id'], target.relative_to(root))
        else:
            print('RETRY READY', target.relative_to(root), item['expected']['job_id'])


if __name__ == '__main__':
    main()
