"""Audit planned question coverage, never certify uncomputed native outcomes."""
import argparse
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('verification_panel',HERE/'panel.py')
panel=importlib.util.module_from_spec(spec);spec.loader.exec_module(panel)

QUESTIONS={
    'minimum_chemical_space': 'Compare the named compact/control spaces, native CI and fragment orbital character; identify which added directions matter to the interaction potential.',
    'optimized_orbital_retention': 'Inspect complete native orbitals and same-geometry subspaces; reject substitution of a required fragment direction.',
    'fragment_and_state_correlation': 'Inspect native CI/root energies and all supported multiplicities at separation against the requested fragment states; no silent root or product relabeling.',
    'path_branch_dependence': 'Compare near/back_near and r6/back6 energies, CI and orbital spans in both bases.',
    'separated_reference': 'Compare 20 and 30 angstrom and use each model/orientation own asymptote.',
    'angular_dependence': 'Compare two independent angular configurations with their corresponding separated references.',
    'relaxation_dependence': 'Inspect near/mid candidate relaxations, unchanged fragment connectivity and all spaces on the same final geometries.',
    'basis_dependence': 'Compare matched aDZ-F12/aTZ-F12 interaction energies and D/T extrapolations.',
}


def audit():
    plan=panel.read(HERE/'plan.json');coverage=plan['coverage'];cases=plan['cases']
    if len(coverage)!=34 or len({r['reaction'] for r in coverage})!=34:
        raise ValueError('Known-channel inventory is incomplete or duplicated.')
    included={r['reaction'] for r in coverage if r['action']!='excluded'}
    inherited={r['reaction'] for r in coverage if r['action']=='existing_native_review'}
    planned={r for c in cases for r in c['reactions']}
    if included!=inherited|planned or inherited&planned:
        raise ValueError('Missing or repeated inherited/planned channel.')
    if len(included)!=29 or len(inherited)!=3 or len(cases)!=25:
        raise ValueError('Unexpected scope; reassess the coverage explicitly.')
    if plan['protocol']['max_jobs']!=23 or panel.MAX_JOBS!=23:
        raise ValueError('Scheduling cap disagrees with the authorized handoff.')
    rows=[]
    for c in cases:
        if len(c['variants'])<2 or not c['reactive_states']:
            raise ValueError('Missing control or reactive-state definition.')
        support={s['multiplicity']:s['roots'] for s in c['state_model']['state_average']}
        for state in c['reactive_states']:
            if not 0<state['root']<=support[state['multiplicity']]:
                raise ValueError('Reactive root outside supported manifold.')
        for basis in plan['protocol']['bases']:
            s=c['sources'][basis]
            names={a['name'] for a in panel.anchors(s)}
            expected={'near','mid','r6','r20','r30','back6','back_near','angle_1','angle_far_1','angle_2','angle_far_2'}
            if names!=expected:raise ValueError('Missing radial/angular/return control.')
            for v in c['variants']:
                panel.partition(s['electrons'],s['source_space'],v['electrons'],v['indices'][basis])
            for member in (s['wfu'],s['output'],s['manifest'],s['configuration']):
                if member not in plan['source_files']:raise ValueError('Unbound source evidence.')
        extra=[]
        if c['id']==27:
            if set(c.get('far_sources',{}))!=set(plan['protocol']['bases']):
                raise ValueError('Missing independent far branches for CH+HSSH.')
            extra.append('Independent inward starts from both far sulfur occupied combinations, plus frontier-only and occupied-union controls.')
        if c['h_exit']:
            extra.append('Pair-expansion CAS and independent MRCI+Q control; external trial orbital character must be established from native results.')
        rows.append(dict(case_id=c['id'],reactions=c['reactions'],
            trials=[dict(name=v['name'],electrons=v['electrons'],indices=v['indices'],purpose=v['purpose']) for v in c['variants']],
            question_checks=QUESTIONS,additional_checks=extra,
            evidence_status='planned_not_native_passed',promotion_ready=False))
    return dict(schema=1,status='known_question_coverage_checked_native_outcomes_pending',
        plan_sha256=panel.file_sha256(HERE/'plan.json'),archive_sha256=plan['archive_sha256'],
        counts=dict(known_channels=34,in_scope=29,excluded=5,inherited_native_reviews=3,
                    channel_identities_with_new_tests=26,distinct_verification_cases=25,max_concurrent_jobs=23),
        inherited_channels=[r for r in coverage if r['action']=='existing_native_review'],
        cases=rows,immediate_promotion_certified=False,
        pending_decisions=[
            'Native numerical and chemical acceptance for each candidate/control; normal termination does not imply acceptance.',
            'Matched-ZPE L1 shallow-well decisions for 13 channels are separate from orbital review. All remain included.',
            'Final basis-specific output/WFU hashes, optimized active indices, roots and weights must be bound into the selected R8 models after review.',
            'The run importer, recovery handoff and R8 input/preflight are separate readiness requirements.'],
        limits=[
            'This is a coverage audit, not a prediction that every calculation or candidate will succeed.',
            'Unexpected state switches, substituted/diffuse trial orbitals or material control differences remain failures requiring scientific resolution.',
            'No numerical tolerance or approved model was relaxed to make this audit pass.'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'pre_submission_audit.json')
    a=p.parse_args();report=audit();a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(report=str(a.output),status=report['status'],counts=report['counts'],
                         immediate_promotion_certified=False),indent=2))
