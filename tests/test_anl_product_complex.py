"""Noncovalent binding uses harmonic differences; VPT2 thermochemistry is retained."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot.anl import automatic
from kinbot.anl.validation import automatic_composite_spec
from tests.test_anl_automatic import _species


def complex_plan():
    root = _species('parent', 'CC=O', ['C','C','O','H','H','H','H'])
    a = _species('methane', 'C', ['C','H','H','H','H'])
    b = _species('co', '[C-]#[O+]', ['C','O'])
    point = _species('complex', 'C.[C-]#[O+]', root.atom)
    root.reac_obj = [SimpleNamespace(instance_name='exit', family_name='r12_insertion_R',
        products=[a,b], ts=_species('ts','',root.atom), do_vdW=True,
        irc_prod_opt=SimpleNamespace(species=point))]
    root.reac_ts_done = [-1]; root.reac_type = ['r12_insertion_R']
    return automatic.build_network_request(root, dict(composite_method='ANL0-F12',
        anl_size_policy=automatic.DEFAULT_SIZE_POLICY, atct_version='1.220'))


def test_complex_is_not_a_cbh_target_but_retains_all_electronic_anl_tasks():
    plan = complex_plan()
    entry = plan['network']['product_complexes']['complex']['record']
    assert 'complex' not in plan['cbh']['targets']
    assert all(e['smiles'] != entry['smiles'] for e in plan['cbh']['calculated_states'])
    ordinary = automatic_composite_spec(entry['molecule'], highest_tier='ANL0-F12')
    assert ordinary['intent']['stationary_point'] == 'minimum'
    assert {'gaussian_vpt2','harmonic','l3_geometry'} <= {t['id'] for t in ordinary['tasks']}
    assert next(t for t in ordinary['tasks'] if t['id']=='l3_geometry')['optimizer']['order']==0


def test_disconnected_complex_assembles_harmonic_binding_profile(tmp_path,monkeypatch):
    entry = complex_plan()['network']['product_complexes']['complex']['record']
    spec = dict(name='anl-automatic-post-l2',molecule=entry['molecule'],tasks=[],
                intent={})
    monkeypatch.setattr(automatic,'_load',lambda p:(tmp_path,spec,dict(tasks={})))
    profiles = automatic.assemble_available_profiles(tmp_path,entry,tmp_path/'profiles',
        state_id='complex')
    assert set(profiles['thermochemistry'])==set(profiles['barrier'])
    for label,path in profiles['barrier'].items():
        assert 'harmonic-ZPE' in label
        record=json.loads(Path(path).read_text())
        assert 'vpt2_correction' not in record['components']
        assert record['zero_k_hartree']==pytest.approx(entry['energy_hartree']+entry['harmonic_zpe_hartree'])


def test_binding_excludes_anharmonic_profiles_without_removing_them(tmp_path,monkeypatch):
    entry = complex_plan()['network']['product_complexes']['complex']['record']
    harmonic = 'L2:B2PLYP-D3BJ/cc-pVTZ:harmonic-ZPE'
    anharmonic = 'L2:B2PLYP-D3BJ/cc-pVTZ:anharmonic-ZPE'
    profiles = dict(barrier={harmonic:'harmonic.json'},
                    thermochemistry={harmonic:'harmonic.json',anharmonic:'vpt2.json'})
    energies={'C':{harmonic:-40.,anharmonic:-41.},
              '[C-]#[O+]':{harmonic:-113.,anharmonic:-114.}}
    before=deepcopy((profiles,energies))
    reads=[]
    def load(smiles,path):
        reads.append(path)
        return -153.01
    monkeypatch.setattr(automatic,'zero_k_from_composite_record',load)
    method,values=automatic._complex_binding_method(entry,
        [dict(smiles=s) for s in energies],profiles,energies)
    assert method==harmonic and reads==['harmonic.json']
    assert values[method]-sum(v[method] for v in energies.values())==pytest.approx(-.01)
    assert (profiles,energies)==before
    del energies['C'][harmonic]
    with pytest.raises(ValueError,match='common accepted harmonic'):
        automatic._complex_binding_method(entry,[dict(smiles=s) for s in energies],profiles,energies)


def test_computed_vpt2_remains_available_but_binding_uses_harmonic_difference(tmp_path,monkeypatch):
    from collections import Counter
    from kinbot.anl.cbh import CBHReaction, FormationEnthalpy, solve_formation_enthalpy, HARTREE_TO_KJ_MOL
    network=complex_plan()['network']
    entry=network['product_complexes']['complex']['record']
    fragments=[r for r in network['stable_species'] if r['routing_name'] in ('methane','co')]
    profiles={};energies={}
    for i,record in enumerate([entry,*fragments]):
        spec=dict(name='anl-automatic-post-l2',molecule=record['molecule'],
                  tasks=[dict(id='gaussian_vpt2')],intent={})
        monkeypatch.setattr(automatic,'_load',lambda p,s=spec:(tmp_path,s,
            dict(tasks=dict(gaussian_vpt2=dict(status='complete')))))
        components=tuple(automatic._l2_source_component(record,key=k,quantity='zpe',
            value=v,method='B2PLYP-D3BJ',basis='cc-pVTZ',state_id=record['routing_name'])
            for k,v in [('harmonic_zpe',record['harmonic_zpe_hartree']),('vpt2_correction',-.001*(i+1))])
        monkeypatch.setattr(automatic,'_vpt2_zpe_components',lambda *a,c=components,**k:c)
        p=automatic.assemble_available_profiles(tmp_path,record,tmp_path/record['routing_name'],state_id=record['routing_name'])
        profiles[record['smiles']]=p
        energies[record['smiles']]={label:automatic.zero_k_from_composite_record(record['smiles'],path)
            for label,path in p['thermochemistry'].items()}
        assert any('anharmonic' in label for label in energies[record['smiles']])
    method,values=automatic._complex_binding_method(entry,fragments,profiles[entry['smiles']],energies)
    selected={r['smiles']:energies[r['smiles']][method] for r in fragments}
    selected[entry['smiles']]=values[method]
    refs={r['smiles']:FormationEnthalpy(r['smiles'],0,'VPT2-inclusive reference',0.,-10.-i,{},
          'test','a'*64,{},uncertainty_kj_mol=0.) for i,r in enumerate(fragments)}
    reaction=CBHReaction(0,entry['smiles'],{entry['smiles']:-1,**dict(Counter(r['smiles'] for r in fragments))},
        {r['smiles']:dict(Counter(r['molecule']['symbols'])) for r in [entry,*fragments]})
    result=solve_formation_enthalpy(reaction,selected,reference_formations=refs,binding_energy=True)
    binding=entry['energy_hartree']+entry['harmonic_zpe_hartree']-sum(
        r['energy_hartree']+r['harmonic_zpe_hartree'] for r in fragments)
    assert result.formation_0k_kj_mol==pytest.approx(-21.+binding*HARTREE_TO_KJ_MOL)
    assert result.references==refs and result.construction=='binding_energy'
    assert result.method==method and 'harmonic-ZPE' in method
