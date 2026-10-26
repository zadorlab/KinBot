"""Read-only native evidence for the four fixed-orbital state diagnostics."""
import argparse
import importlib.util
import json
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('production', type=Path)
    args = parser.parse_args()
    script = Path(__file__).resolve().parent.parent/'r7_compact_union_20261025/compare_orbitals.py'
    spec = importlib.util.spec_from_file_location('orbital_comparison', script)
    comparison = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(comparison)
    root = args.production.expanduser().resolve()/'orbital_recon'
    saved = root/'root_spectrum_20261026/csh2_sh_exit/review_launch.json'
    record = json.loads(saved.read_text())
    for system in record['systems']:
        manifest = json.loads(Path(system['manifest']).read_text())
        for job in manifest['jobs']:
            directory = Path(job['directory']); name = job['name']
            print('JOB', json.dumps({k:job.get(k) for k in (
                'name', 'status', 'validation_errors', 'basis', 'state_energies_hartree',
                'initial_record', 'orbital_predecessor', 'initial_wfu_sha256', 'wfu_sha256')}, sort_keys=True))
            for filename in ('execution.json', name+'.native.inp'):
                p = directory/filename
                if p.is_file():
                    print('FILE', p, '\n'+p.read_text())
            p = directory/(name+'.out')
            if p.is_file():
                text = p.read_text(errors='replace')
                start = re.search(r'(?im)^\s*PROGRAM\s*\*\s*MULTI\b', text)
                print('NATIVE_MULTI', str(p), '\n'+(text[start.start():] if start else 'NO MULTI SECTION'))
            if job.get('status') == 'complete':
                b = name.split('__')[0]
                old = Path('compact_resolution_20261025/csh2_sh_exit_resolution/csh2_sh_exit_resolution_union')/name/(name+'.natural.orbitals.xml')
                new = directory.relative_to(root)/(name+'.natural.orbitals.xml')
                print('COMPARISON', b, json.dumps(comparison.compare(root, old, new, {
                    'source_union_vs_fixed_CASCI_active': (list(range(17,23)), list(range(17,23)))}), sort_keys=True))


if __name__ == '__main__':
    main()
