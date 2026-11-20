from pathlib import Path
import json,gzip,hashlib,re
from math import dist
from kinbot.constants import AUtoKCAL
R=Path('local_runs/r8_bindings_20261119');P=Path('examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/review_20261119_bindings')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
ledger=read(P/'binding_decisions.json');raw=Path('examples/anl/tmp/r9-orbital-final.yjxXq5MV.txt')
points={};returns={};summary=None
for line in raw.open():
 key,_,s=line.partition(' ')
 if key not in ('NATIVE_POINT','RETURN_CHECK','SUMMARY'):continue
 a=json.loads(s)
 if key=='NATIVE_POINT':points[tuple(a['identity'])]=a
 elif key=='RETURN_CHECK':returns[tuple(a['identity'])]=a
 else:summary=a
assert len(points)==578 and summary['production_model_accepted'] is False
recomputed=0
# Recompute each listed value from original exported energies, with an explicit
# configuration-name join only for full campaigns seeded from earlier r30.
for m in ledger['h_loss_models']:
 title=m['title'];basis='avdz-f12' if 'avdz' in title else 'avtz-f12';base_title=title.replace('_full_','_')
 def energy(point):
  key=(title,basis,point)
  if key not in points:key=(base_title,basis,point)
  p=points[key];assert p['normal_termination'] and not p['native_errors'] and not p['native_issues']
  energies=[v['energy_hartree'] for v in p['correlation'] if v['root']==1 and v['ipea']==.25 and v['shift']==.1 and v['native_method']=='CASPT2-F12']
  assert len(energies)==1
  return energies[0]
 ref=energy('r30')
 for pt,want in m['energies_from_own_r30_kcal'].items():
  assert abs((energy(pt)-ref)*AUtoKCAL-want)<1e-10;recomputed+=1
 r=[v for k,v in returns.items() if k[0]==title]
 assert m['return_min_squared_singular_value']==min((x for v in r for x in v.get('active_subspace_squared_singular_values',[])),default=None)
f=read(P/'h_loss_fragment_binding.json')
for v in [dict(fragment=f['target_geometry'],distances=f['target_distances'])]+f['checks']:
 g=v['fragment'];s=g['symbols'];x=g['positions'];i=s.index('C');ss=[j for j,v in enumerate(s) if v=='S']
 assert sorted(dist(x[i],x[j]) for j in ss)==v['distances']['C_S_angstrom'] or max(abs(a-b) for a,b in zip(sorted(dist(x[i],x[j]) for j in ss),v['distances']['C_S_angstrom']))<1e-14
 assert abs(dist(x[ss[0]],x[ss[1]])-v['distances']['S_S_angstrom'])<1e-14
n=json.loads(gzip.decompress((P/'network_binding_evidence.json.gz').read_bytes()))
for r in n['reverse']:
 log=Path(r['native_log']);assert sha(log)==r['native_log_sha256']
 # Read all native frequency rows; the final negative mode must occur explicitly.
 freq=[float(x) for line in log.read_text(errors='replace').splitlines() if 'Frequencies --' in line for x in line.split('--',1)[1].split()]
 assert min(abs(x-r['imaginary_cm1']) for x in freq)<1e-6
 assert sum(v<0 for v in r['native_frequencies'])==1
bound=gzip.decompress((P/'ch3ss_provenance_bound.inp.gz').read_bytes()).decode();source=Path(n['input'])
phys=lambda text:' '.join(line.split('!',1)[0] for line in text.splitlines()).split()
assert phys(bound)==phys(source.read_text())
receipt=read(source.parent/'execution.json');assert receipt['returncode']==0 and receipt['input_sha256']==sha(source)
report=dict(raw_orbital_export_sha256=sha(raw),native_points=578,recomputed_relative_energies=recomputed,reverse_native_logs_verified=3,fragment_comparisons_verified=len(f['checks']),physical_MESS_tokens_identical=True,prior_native_receipt=str(source.parent/'execution.json'),prior_native_receipt_sha256=sha(source.parent/'execution.json'),production_acceptance=False)
(P/'independent_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
