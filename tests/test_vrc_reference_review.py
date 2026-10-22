from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from kinbot import vrc_reference_review as review
from kinbot.execution import atomic_json, file_sha256


@pytest.fixture
def native_review(tmp_path):
    # An exactly represented five-electron quartet in a CAS(3,4).
    # These small matrices test analysis, provenance and acceptance separately.
    par=dict(mult=2, vrc_tst_correction_model='spin_splitting',
        vrc_tst_active_space='explicit',
        vrc_tst_scan_method='caspt2-f12(3,4)', vrc_tst_sample_method='caspt2-f12(3,4)',
        vrc_tst_sample_basis='avdz-f12', vrc_tst_high_method='caspt2-f12(3,4)',
        vrc_tst_high_basis=['avdz-f12','avtz-f12'],
        vrc_tst_spin_splitting=dict(target_multiplicity=2,reference_multiplicity=4))
    token=review.digest(review.contract(par))
    config=tmp_path/'probe.json';atomic_json(config,dict(reference_review_contract_sha256=token))
    occ=np.diag([2.,1.,1.,1.,0.,0.,0.]);eye=np.eye(7)
    checks=[]
    for basis in ('avdz-f12','avtz-f12'):
        for distance in (3.,30.):
            d=tmp_path/f'{basis}-{distance}';d.mkdir()
            matrices={}
            for name,value in dict(s=eye,cas=eye,hf=eye,dhf=occ,dr1=occ).items():
                p=d/f'job.identity_{name}.out'
                p.write_text('BEGIN_DATA\n'+'\n'.join(' '.join(map(str,row)) for row in value)+'\nEND_DATA\n')
                matrices[p.name]=file_sha256(p)
            native=d/'job.native.inp';native.write_text('native reference input fixture')
            atomic_json(d/'execution.json',dict(reference_identity_matrices=matrices,
                native_input=str(native),native_input_sha256=file_sha256(native),
                reference_review_contract_sha256=token))
            output=d/'job.out';output.write_text('native CC/MR fixture')
            job=dict(name='job',status='complete',directory=str(d),basis=basis,
                distance_angstrom=distance,output_sha256=file_sha256(output),
                probe_active_space=dict(active_electrons=3,active_orbitals=4),
                correlated_energies_hartree=[dict(method='uccsd(t)-f12b',multiplicity=4,root=1,energy_hartree=-4.),
                    dict(method='caspt2-f12(3,4)',multiplicity=4,root=1,energy_hartree=-3.9)])
            manifest=d/'reconnaissance.json'
            atomic_json(manifest,dict(source=str(config),jobs=[job],electron_count=5))
            checks.append(dict(manifest=str(manifest),job='job'))
    return par,checks,tmp_path


def approve(par,checks,path):
    receipt=review.inspect(par,checks)
    receipt.update(approved=True,findings={k:'Explicit reviewed fixture evidence' for k in
        ('state_pairing','cc_quality','degenerate_subspaces','orbital_continuity')})
    atomic_json(path,receipt)
    return dict(par,vrc_tst_reference_review=dict(path=str(path),sha256=file_sha256(path)))


def test_reference_requires_native_analysis_and_explicit_review(native_review):
    par,checks,tmp=native_review
    with pytest.raises(ValueError,match='native CC/MR reference review'): review.validate(par,tmp)
    # A direct potential does not require or consume these CC diagnostics.
    review.validate(dict(par,vrc_tst_correction_model='direct'),tmp)
    bound=approve(par,checks,tmp/'review.json')
    review.validate(bound,tmp)
    d=Path(checks[0]['manifest']).parent
    (d/'job.identity_dr1.out').write_text('changed')
    with pytest.raises(ValueError,match='export changed'): review.validate(bound,tmp)


def test_reference_rejects_missing_basis_failure_and_changed_contract(native_review):
    par,checks,tmp=native_review
    with pytest.raises(ValueError,match='every CC/MR basis'): review.inspect(par,checks[:2])
    with pytest.raises(ValueError,match='exact reviewed contract'):
        review.inspect(dict(par,vrc_tst_scan_shift=.2),checks)
    p=Path(checks[0]['manifest']);m=json.loads(p.read_text());m['jobs'][0]['status']='incomplete';atomic_json(p,m)
    with pytest.raises(ValueError,match='completed, validated'): review.inspect(par,checks)


def test_reference_preserves_unchanged_legacy_models_but_gates_additions():
    from tests.test_input_refresh import surface
    model=surface();par=dict(vrc_tst_correction_model='spin_splitting')
    review.validate_additions(model,deepcopy(model),par,Path.cwd())
    with pytest.raises(ValueError,match='native CC/MR reference review'):
        review.validate_additions(None,model,par,Path.cwd())


def test_preset_cannot_bypass_reference_review(tmp_path):
    par = dict(vrc_tst_accuracy_preset='anl1-f12', vrc_tst_active_space='explicit')
    with pytest.raises(ValueError, match='native CC/MR reference review'):
        review.validate(par, tmp_path)
    review.validate(dict(par, vrc_tst_correction_model='direct'), tmp_path)


def test_every_example_leaf_retains_scan_model_and_changes_only_energy_assembly(tmp_path, monkeypatch):
    from kinbot.vrc_restart import revised_input, validate_input
    from kinbot.vrc_correction import correction_components, resolved_state_model, spin_splitting_settings
    from kinbot.vrc_multistate import state_parameters
    from kinbot.parameters import Parameters
    import hashlib
    root=Path(__file__).parents[1]
    monkeypatch.setattr('kinbot.anl.site._executable',lambda *a:'native-tool')
    for name in ('ch3s','ch3ss'):
        payload=json.loads((root/f'examples/anl/{name}_anl1_f12/{name}.json').read_text())
        def seeds(value):
            if isinstance(value,dict):
                if value.get('mode')=='molpro_record':
                    p=tmp_path/(hashlib.sha256(value['path'].encode()).hexdigest()+'.wfu')
                    p.write_bytes(b'render-only reviewed orbital fixture')
                    value.update(path=str(p),sha256=file_sha256(p))
                for v in value.values():seeds(v)
            elif isinstance(value,list):
                for v in value:seeds(v)
        seeds(payload)
        payload['vrc_tst_correction_model']='spin_splitting'
        for channel in payload['vrc_tst_surface_model']['channels']:
            for stage in channel['bottlenecks']:
                stage['parameters']['vrc_tst_correction_model']='spin_splitting'
        path=tmp_path/f'{name}.json';atomic_json(path,payload)
        validate_input(path,{'direct_caspt2_cbs':True})
        original=Parameters(str(path),show_warnings=False).par
        direct=revised_input(original,{'direct_caspt2_cbs':True})
        assert direct['vrc_tst_scan_points']==original['vrc_tst_scan_points']
        assert direct['rotdpy_dist']==original['rotdpy_dist']
        for a,b in zip(original['vrc_tst_surface_model']['channels'],direct['vrc_tst_surface_model']['channels']):
            for sa,sb in zip(a['bottlenecks'],b['bottlenecks']):
                for state in sa['states']:
                    before=state_parameters(original,sa,state,[2,3])
                    after=state_parameters(direct,sb,state,[2,3])
                    spin=spin_splitting_settings(before,parent_multiplicity=state['multiplicity'],fragment_multiplicities=[2,3])
                    for basis in ('avdz-f12','avtz-f12'):
                        old=resolved_state_model(before,None,0,spin,basis=basis)
                        new=resolved_state_model(after,None,0,None,basis=basis)
                        assert old==new
                    terms=correction_components(after)
                    assert [c['role'] for c in terms]==['sampling','trusted:0','trusted:1']
                    assert all('ccsd' not in c['method'] for c in terms)
