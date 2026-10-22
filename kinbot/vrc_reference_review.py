"""Explicit, native-evidence review gate for new spin-splitting selections.

Density comparisons are diagnostics, never an automatic chemical approval.
Legacy unchanged models can continue; newly accepted spin-splitting models
must carry this review. Direct MR potentials do not need a CC reference.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from kinbot.execution import atomic_json, file_sha256


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def contract(parameters):
    from kinbot.vrc_correction import expected_levels, spin_splitting_settings
    spin = spin_splitting_settings(parameters, parent_multiplicity=int(parameters.get('mult', 1)))
    return dict(levels=expected_levels(parameters, settings=spin),
                target_root=parameters.get('vrc_tst_selected_roots', {}).get(
                    str(spin['target_multiplicity']), 1),
                scan_method=parameters.get('vrc_tst_scan_method'),
                shift=parameters.get('vrc_tst_scan_shift', .1),
                ipea=parameters.get('vrc_tst_scan_ipea', .25),
                asymptote_angstrom=parameters.get('vrc_tst_asymptote_distance', 30.))


def inspect(parameters, checks):
    """Re-read native results and bind both identity diagnostics and CC output."""
    from kinbot.molpro_reference_identity import report
    from kinbot.vrc_correction import spin_splitting_settings
    scientific = json.loads(json.dumps(contract(parameters)))
    token = digest(scientific)
    spin = spin_splitting_settings(parameters, parent_multiplicity=int(parameters.get('mult', 1)))
    reference = spin['reference_multiplicity']
    bases = spin['reference']['basis']
    bases = set(bases if isinstance(bases, list) else [bases])
    for level in spin['gap_levels']:
        bases.update(level['basis'] if isinstance(level['basis'], list) else [level['basis']])
    coverage = set()
    evidence = []
    for check in checks:
        path = Path(check['manifest']).expanduser().resolve(strict=True)
        manifest = json.loads(path.read_text())
        source = Path(manifest['source'])
        config = json.loads(source.read_text())
        if config.get('reference_review_contract_sha256') != token:
            raise ValueError('Reference probe was not staged for this exact reviewed contract.')
        jobs = [j for j in manifest['jobs'] if j['name'] == check['job']]
        if len(jobs) != 1:
            raise ValueError('Reference check must identify one native job.')
        job = jobs[0]
        if job.get('status') != 'complete' or job.get('validation_errors'):
            raise ValueError('Reference review requires completed, validated MR and CC calculations.')
        correlated = job.get('correlated_energies_hartree', [])
        cc = [e for e in correlated if e['multiplicity'] == reference
              and e['method'].lower() == spin['reference']['method'].lower()]
        mr = [e for e in correlated if e['multiplicity'] == reference
              and e['method'].lower() != spin['reference']['method'].lower()]
        if len(cc) != 1 or not mr:
            raise ValueError('Reference check lacks native CC and MR reference energies.')
        directory = Path(job['directory'])
        result = report(directory, job['name'], job['probe_active_space'],
                        manifest['electron_count'], reference)
        if not result.get('reference_states'):
            raise ValueError('Reference check lacks state-specific native density analysis.')
        output = directory/(job['name']+'.out')
        if file_sha256(output) != job['output_sha256']:
            raise ValueError('Reference native output changed.')
        receipt = directory/'execution.json'
        execution = json.loads(receipt.read_text())
        native = Path(execution['native_input'])
        if file_sha256(native) != execution['native_input_sha256']:
            raise ValueError('Executed native reference input changed.')
        if execution.get('reference_review_contract_sha256') != token:
            raise ValueError('Native execution did not record the staged reference review contract.')
        distance = job.get('distance_angstrom')
        if distance is None:
            distance = job.get('reference_review_distance_angstrom')
        role = ('asymptote' if distance == scientific['asymptote_angstrom']
                else 'interior' if isinstance(distance, (float,int))
                and 0 < distance < scientific['asymptote_angstrom'] else None)
        if role is None:
            raise ValueError('Reference check needs a finite interior or the configured asymptote.')
        coverage.add((job['basis'], role))
        files = {str(p):file_sha256(p) for p in (source, output, receipt, native)}
        files.update({str(directory/name):sha for name,sha in execution['reference_identity_matrices'].items()})
        evidence.append(dict(manifest=str(path), job=job['name'], job_sha256=digest(job),
                             basis=job['basis'], role=role, files=files, analysis=result))
    if any((basis, role) not in coverage for basis in bases for role in ('interior','asymptote')):
        raise ValueError('Reference review needs native interior and asymptote checks in every CC/MR basis.')
    return dict(schema=1, contract=scientific, contract_sha256=token, checks=evidence)


def validate(parameters, base):
    from kinbot.vrc_correction import correction_model, apply_vrc_accuracy_preset
    # Raw refresh inputs can inherit spin splitting from the accuracy preset.
    # Resolve that choice before deciding whether the native review is needed.
    if 'vrc_tst_correction_model' not in parameters:
        parameters = apply_vrc_accuracy_preset(deepcopy(parameters), parameters.keys())
    if correction_model(parameters) != 'spin_splitting':
        return
    binding = parameters.get('vrc_tst_reference_review')
    if not isinstance(binding, dict) or set(binding) != {'path', 'sha256'}:
        raise ValueError('New spin-splitting selections require a native CC/MR reference review '
                         '(vrc_tst_reference_review); orbital occupations alone do not approve the reference.')
    path = Path(binding['path']).expanduser()
    if not path.is_absolute():
        path = Path(base)/path
    if file_sha256(path) != binding['sha256']:
        raise ValueError('CC/MR reference review receipt changed.')
    saved = json.loads(path.read_text())
    fields = ('state_pairing', 'cc_quality', 'degenerate_subspaces', 'orbital_continuity')
    if (saved.get('approved') is not True or
            any(not str(saved.get('findings', {}).get(k, '')).strip() for k in fields)):
        raise ValueError('Explicit reference-state, CC-quality, degeneracy and continuity findings are required.')
    fresh = inspect(parameters, saved.get('checks', []))
    if any(saved.get(k) != v for k,v in fresh.items()):
        raise ValueError('CC/MR review is stale or belongs to a different native/scientific contract.')


def validate_additions(previous, current, parameters, base):
    """Grandfather only byte-equivalent existing channels, never new choices."""
    old = (previous or {}).get('channels', [])
    for channel in current['channels']:
        if channel in old:
            continue
        for stage in channel['bottlenecks']:
            for state in stage['states']:
                leaf = {**parameters, **deepcopy(stage['parameters']), 'mult':state['multiplicity']}
                leaf['vrc_tst_selected_roots'] = {**leaf.get('vrc_tst_selected_roots', {}),
                                                str(state['multiplicity']):state['root']}
                validate(leaf, base)


def main():
    import argparse
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='contract':
        parser=argparse.ArgumentParser(description='Fingerprint the exact proposed scalar spin-splitting model.')
        parser.add_argument('parameters',type=Path)
        args=parser.parse_args(sys.argv[2:])
        scientific=contract(json.loads(args.parameters.read_text()))
        print(json.dumps(dict(contract=scientific,reference_review_contract_sha256=digest(scientific)),indent=2))
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('parameters', type=Path, help='Resolved scalar leaf parameters, including selected orbital model.')
    parser.add_argument('checks', type=Path, help='JSON list of manifest/job pairs.')
    parser.add_argument('output', type=Path)
    parser.add_argument('--findings', type=Path, help='Explicit reviewed findings; omitting this produces an unapproved draft.')
    args = parser.parse_args()
    parameters = json.loads(args.parameters.read_text())
    result = inspect(parameters, json.loads(args.checks.read_text()))
    result['approved'] = args.findings is not None
    result['findings'] = json.loads(args.findings.read_text()) if args.findings else {}
    if args.findings and any(not str(result['findings'].get(k,'')).strip() for k in
            ('state_pairing','cc_quality','degenerate_subspaces','orbital_continuity')):
        raise ValueError('All four scientific review findings are required for approval.')
    if args.output.exists():
        raise FileExistsError('Use a new review receipt path to retain prior evidence.')
    atomic_json(args.output, result)
    print(json.dumps(dict(path=str(args.output.resolve()), sha256=file_sha256(args.output))))



def validate_input_change(previous, current, base):
    """Check changed/new contracts in an explicit full controller refresh."""
    from kinbot.vrc_multistate import normalize_inline_model_paths, load_surface_model
    a,b=deepcopy(previous),deepcopy(current)
    normalize_inline_model_paths(a,base);normalize_inline_model_paths(b,base)
    def surfaces(p):
        result=dict(p.get('vrc_tst_surface_models',{}))
        loaded=load_surface_model(dict(p,input_file=str(Path(base)/'input.json')))
        if loaded: result['__initial__']=loaded[0]
        return result
    old,new=surfaces(a),surfaces(b)
    keys=('vrc_tst_correction_model','vrc_tst_accuracy_preset','vrc_tst_active_space',
          'vrc_tst_state_model','vrc_tst_state_model_file','vrc_tst_spin_splitting',
          'vrc_tst_scan_shift','vrc_tst_scan_ipea','vrc_tst_asymptote_distance',
          'vrc_tst_state_average_roots','vrc_tst_selected_roots')
    changed=any(a.get(k)!=b.get(k) for k in keys)
    for key,model in new.items():
        validate_additions(None if changed else old.get(key),model,b,base)
    if not new and changed:
        validate(b,base)


if __name__ == '__main__':
    main()
