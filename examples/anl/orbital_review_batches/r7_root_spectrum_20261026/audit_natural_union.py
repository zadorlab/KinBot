"""Reproduce the completed natural-union audit from two verified exports.

Arguments are the extracted compact-continuation and natural-union archives.
Only reads files; writes JSON to stdout. WFU integrity is checked separately by
the export verifier. All MO comparisons use identical geometries and AO bases.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re

from kinbot.molpro_diagnostics import parse_native_output


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('older_export',type=Path)
    parser.add_argument('newer_export',type=Path)
    args=parser.parse_args()
    old_root=args.older_export/'production/orbital_recon'
    new_root=args.newer_export/'production/orbital_recon'
    script=Path(__file__).resolve().parent.parent/'r7_compact_union_20261025/compare_orbitals.py'
    spec=importlib.util.spec_from_file_location('orbital_comparison',script)
    comparison=importlib.util.module_from_spec(spec);spec.loader.exec_module(comparison)
    result=dict(schema=1, jobs=[], comparisons=[],
        method='Native MULTI results and complete same-geometry MO change of basis; '
               'AO coefficients are qualitative character, not populations.')
    for b in ('adz','atz'):
        for anchor in ('exit_l2','exit_r_30p0'):
            name=b+'__'+anchor
            rel=Path('natural_union_20261025/csh2_sh_exit/csh2_sh_exit_natural_union')/name
            d=new_root/rel
            native=(d/(name+'.out')).read_text()
            audit=parse_native_output(native)
            section=audit['native_sections'][-1]
            tables=[t for t in section['orbital_tables'] if 'NATURAL' in t['heading']]
            orbitals=[]
            for o in tables[-1]['orbitals']:
                if 16<=int(o['index'].split('.')[0])<=22:
                    o=dict(o)
                    o['ao_coefficients']=sorted(o.get('ao_coefficients',[]),
                        key=lambda a:abs(a['coefficient']),reverse=True)[:6]
                    orbitals.append(o)
            x=json.loads((d/'execution.json').read_text())
            hashes={f:hashlib.sha256((d/f).read_bytes()).hexdigest() for f in
                    (name+'.out',name+'.native.inp',name+'.natural.orbitals.xml','execution.json')}
            result['jobs'].append(dict(path=str(rel),sha256=hashes,
                execution={k:x.get(k) for k in ('status','slurm_job_id','input_transform','returncode')},
                normal_termination=audit['normal_termination_seen'],
                state_energies=section['state_energies'],orbitals=orbitals,
                ci_tables=section['ci_tables'],
                real_times_seconds=[float(v) for v in re.findall(r'REAL TIME\s*\*\s*([\d.]+) SEC',native)]))
        name=b+'__exit_l2'
        left=old_root/'compact_resolution_20261025/csh2_sh_exit_resolution/csh2_sh_exit_resolution_union'/name/(name+'.natural.orbitals.xml')
        right=new_root/'natural_union_20261025/csh2_sh_exit/csh2_sh_exit_natural_union'/name/(name+'.natural.orbitals.xml')
        c=comparison.compare(Path('/'),left,right,{
            'intended_union_vs_optimized_active':(list(range(17,23)),list(range(17,23))),
            'replace_new_17_with_inactive_SH_16':([16,18,19,20,21,22],list(range(17,23)))})
        c['sha256']={'older/'+str(left.relative_to(old_root)):c['sha256'][str(left)],
                     'newer/'+str(right.relative_to(new_root)):c['sha256'][str(right)]}
        c['left']='older/'+str(left.relative_to(old_root))
        c['right']='newer/'+str(right.relative_to(new_root))
        result['comparisons'].append(c)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
