"""Reproduce fixed-orbital results from verified extracted native evidence.

Arguments: extracted compact-continuation archive, extracted root-spectrum archive.
Read-only. Check source output/XML hashes and unchanged-geometry MO spans; never
select a model or infer fragment ground-state ordering from a fixed-orbital gap.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re

from kinbot.molpro_diagnostics import parse_native_output


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def audit(previous, current):
    plan=Path(__file__).resolve().parent
    compare_script=plan.parent/'r7_compact_union_20261025/compare_orbitals.py'
    spec=importlib.util.spec_from_file_location('compare_native',compare_script)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    new=current/'production/orbital_recon/root_spectrum_20261026/csh2_sh_exit/csh2_sh_exit_fixed_spectrum'
    old=previous/'production/orbital_recon/compact_resolution_20261025/csh2_sh_exit_resolution/csh2_sh_exit_resolution_union'
    manifest=json.loads((new/'reconnaissance.json').read_text());results=[]
    source_jobs={j['name']:j for j in json.loads((old/'reconnaissance.json').read_text())['jobs']}
    for j in manifest['jobs']:
        name=j['name'];out=new/name/(name+'.out');xml=new/name/(name+'.natural.orbitals.xml')
        if digest(out)!=j['output_sha256'] or digest(xml)!=j['orbital_xml']['natural']['sha256']:
            raise ValueError('Changed native output or XML: '+name)
        source_out=old/name/(name+'.out');source_xml=old/name/(name+'.natural.orbitals.xml')
        source_job=source_jobs[name]
        if (digest(source_out)!=source_job['output_sha256'] or
                digest(source_xml)!=source_job['orbital_xml']['natural']['sha256']):
            raise ValueError('Changed source output or XML: '+name)
        native=parse_native_output(out.read_text());section=next(s for s in native['native_sections'] if s['program']=='MULTI')
        orbitals=next(t['orbitals'] for t in section['orbital_tables'] if 'NATURAL' in t['heading'])
        # Absolute paths are accepted by pathlib's join; both XML files are
        # independently verified in the exported archives before this audit.
        cmp=m.compare(Path('/'),source_xml,xml,{
            'source_union_vs_fixed_active':(list(range(17,23)),list(range(17,23)))})
        cmp.pop('right_orbitals_in_left_basis')
        # Make the saved evidence portable across archive extraction locations.
        cmp['left']=str(source_xml.relative_to(previous))
        cmp['right']=str(xml.relative_to(current))
        cmp['sha256']={cmp['left']:digest(source_xml),cmp['right']:digest(xml)}
        energies=j['state_energies_hartree'];e={(x['multiplicity'],x['root']):x['energy_hartree'] for x in energies}
        results.append(dict(job=name,basis=j['basis'],status=j['status'],
            source_output_sha256=digest(source_out),native_output_sha256=digest(out),
            natural_xml_sha256=digest(xml),wfu_sha256=j['wfu_sha256'],
            native_state_energies=section['state_energies'],reported_state_energies=energies,
            gap_D2_D1_cm_inverse=(e[2,2]-e[2,1])*219474.6313632,
            gap_D3_D1_cm_inverse=(e[2,3]-e[2,1])*219474.6313632,
            gap_Q1_D3_cm_inverse=(e[4,1]-e[2,3])*219474.6313632,
            active_orbitals=[o for o in orbitals if o.get('role')=='active'],
            ci_tables=section['ci_tables'],subspace_comparison=cmp,
            real_times_seconds=[float(x) for x in re.findall(r'REAL TIME\s*\*\s*([\d.]+)\s*SEC',out.read_text())]))
    return dict(scope='Four fixed-orbital CASCI diagnostics; no optimized production model or isolated-fragment spin ordering inferred.',jobs=results)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('compact_archive',type=Path);parser.add_argument('spectrum_archive',type=Path)
    args=parser.parse_args()
    print(json.dumps(audit(args.compact_archive,args.spectrum_archive),indent=2))
