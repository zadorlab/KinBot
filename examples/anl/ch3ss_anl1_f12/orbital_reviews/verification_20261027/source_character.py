"""Extract traceable starting-orbital evidence; coefficients are not populations."""
import argparse
import gzip
import json
from pathlib import Path

from kinbot.execution import file_sha256
from kinbot.molpro_diagnostics import parse_native_output

HERE = Path(__file__).resolve().parent


def audit(root):
    plan = json.loads((HERE / 'plan.json').read_text())
    entries = []
    for case in plan['cases']:
        for basis, source in case['sources'].items():
            selected = sorted({i for v in case['variants'] for i in v['indices'][basis]})
            sources = [('near', source['output'])]
            if 'far_sources' in case:
                sources.append(('far', case['far_sources'][basis]['files']['.out']))
            for location, member in sources:
                output = root / member
                inp = output.with_suffix('.inp')
                expected = plan['source_files'][member]['sha256']
                if file_sha256(output) != expected:
                    raise ValueError('Changed source: ' + str(output))
                parsed = parse_native_output(output.read_text(), inp.read_text())
                sections = [s for s in parsed['native_sections'] if s['program'] == 'MULTI']
                tables = [t for s in sections for t in s['orbital_tables']
                          if 'NATURAL ORBITALS' in t['heading']]
                if not tables:
                    raise ValueError('Missing native natural-orbital table: ' + member)
                orbitals = {int(o['index'].split('.')[0]): o for o in tables[-1]['orbitals']}
                entries.append(dict(
                    case_id=case['id'], basis=basis, location=location,
                    output_member=member, output_sha256=expected,
                    input_sha256=file_sha256(inp), symbols=source['symbols'],
                    coordinate=source['coordinate'], heading=tables[-1]['heading'],
                    table_line=tables[-1]['line'],
                    requested_indices=selected,
                    native_orbitals=[orbitals[i] for i in selected if i in orbitals],
                    unprinted_indices=[i for i in selected if i not in orbitals]))
    return dict(
        schema=1, archive_sha256=plan['archive_sha256'], entries=entries,
        interpretation='Native printed AO coefficients support qualitative character assignment. '
                       'They are not orthogonalized orbital populations or proof of optimized retention. '
                       'An unprinted external trial direction requires the full XML and new native result.',
        production_model_accepted=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('export_root', type=Path)
    p.add_argument('--output', type=Path, default=HERE / 'source_character.json.gz')
    a = p.parse_args()
    result = audit(a.export_root)
    with gzip.open(a.output, 'wt') as f:
        json.dump(result, f)
    print(json.dumps(dict(report=str(a.output), source_tables=len(result['entries']),
                          unprinted_indices=sum(len(e['unprinted_indices']) for e in result['entries']),
                          production_model_accepted=False)))
