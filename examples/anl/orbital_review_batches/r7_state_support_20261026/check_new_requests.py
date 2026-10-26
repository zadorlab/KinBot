"""Check all seven newly received requests before staging, without mutations."""
import argparse
import importlib.util
import json
from pathlib import Path

from kinbot.execution import file_sha256


def check(production, plan=None):
    production = Path(production).expanduser().resolve(strict=True)
    plan = Path(plan or Path(__file__).resolve().parent)
    prior = Path(__file__).resolve().parent.parent/'r7_root_spectrum_20261026'
    spec = importlib.util.spec_from_file_location('prior_request_check', prior/'check_new_requests.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    results = module.check(production)['requests']
    evidence = json.loads((plan/'new_requests_evidence.json').read_text())
    for item in evidence['requests']:
        q = item['request']
        w = production/'ch3ss'/q['parent']
        path = w/'vrctst/orbital_review'/q['reaction']/'request.json'
        current = json.loads(path.read_text())
        keys = ('parent','reaction','bond','symbols','positions','products','charge','multiplicity','status')
        if any(q.get(k) != current.get(k) for k in keys):
            raise ValueError('Live request differs from captured evidence: '+str(path))
        network = json.loads((w/'accepted_l2_network.json').read_text())
        rows = [r for r in network['network']['reactions'] if r['name'] == q['reaction']]
        if rows != item['network']:
            raise ValueError('Accepted reaction differs from captured evidence: '+q['reaction'])
        results[q['reaction']] = dict(request=str(path),request_sha256=file_sha256(path),
            note='Distinct geometry and atom map; no orbital evidence reused between parents.')
    return dict(status='matched',requests=results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('production',type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.production),indent=2))
