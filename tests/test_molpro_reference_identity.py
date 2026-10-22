import json

import numpy as np
import pytest

from kinbot.molpro_reference_identity import compare, matrix, report
from kinbot.execution import atomic_json, file_sha256


@pytest.mark.parametrize('multiplicity', [1, 2, 3, 4, 5])
@pytest.mark.parametrize('selector', ['spin', 'ms2', 'saved_dm'])
@pytest.mark.parametrize('reference_first', [False, True])
def test_individual_densities_saved_to_explicit_record_and_legacy_repaired(multiplicity, selector, reference_first):
    from kinbot.molpro_reference_identity import (natural_cards, export_deck,
        repair_legacy_selectors, density_sets, validate_density_saves)
    entries = [dict(multiplicity=multiplicity,roots=2),dict(multiplicity=multiplicity+2,roots=3)]
    if not reference_first: entries.reverse()
    model = {'state_average':entries}; ms2=multiplicity-1
    cards = natural_cards(model,multiplicity)
    multi = '{multi;\n'+''.join(f"wf,25,1,{e['multiplicity']-1}\nstate,{e['roots']}\n" for e in entries)
    multi += cards+'canonical,2141.2\nnatorb,2140.2\n}\n'
    markers = ''.join(f'reconr{i}=energy({i})\nshow,reconr{i}\n' for i in range(1,6))
    correlation = '{rs2-f12;wf,25,1,3;}\n{rhf;wf,25,1,3;}\n'
    expected = multi+markers+export_deck('probe',model,multiplicity,part='mr')+correlation+export_deck('probe',model,multiplicity,part='hf')
    legacy = multi+markers+correlation+'{matrop\n'
    matrices=[('s','s'),('cas','orb,2140.2,natural'),('hf','orb,2150.2'),('dhf','den,2150.2')]
    matrices += [(f'dr{r}',f'den,3000.2,state={r}.1,ms2={ms2}') for r in (1,2)]
    legacy += ''.join(f'load,{k},{v}\nwrite,{k},probe.identity_{k}.out,new,scientific\n' for k,v in matrices)+'}\n'
    if selector != 'saved_dm':
        legacy_cards=''.join(f'natorb,{3000+r}.2,state={r}.1,spin={ms2}\n' for r in (1,2))
        legacy=legacy.replace(cards,legacy_cards)
        for r in (1,2):
            legacy=legacy.replace(f'den,3000.2,state={r}.1,ms2=',f'den,{3000+r}.2,state={r}.1,{selector}=')
    unrelated='\nload,user,den,2140.2,state=1.1,spin=1\n'
    fixed,repairs=repair_legacy_selectors(legacy+unrelated)
    assert fixed.replace('\n\n','\n') == (expected+unrelated).replace('\n\n','\n')
    assert repairs[-1]['root_sets'] == density_sets(model,multiplicity)
    assert fixed.count('orbital,3000.2') == 1
    assert repair_legacy_selectors(fixed) == (fixed,[])
    assert repair_legacy_selectors(unrelated) == (unrelated,[])
    output=' PROGRAM * MULTI (Direct Multiconfiguration SCF)\n'
    roots=[r for e in entries for r in range(1,e['roots']+1)]
    output+=''.join(f' First order charge density matrix for state {r}.1 saved on record 3000.2 (density set {i})\n' for i,r in enumerate(roots,1))
    output+=' PROGRAM * MATROP\n'
    assert validate_density_saves(fixed,output)['native_save_sequence_verified']
    with pytest.raises(ValueError,match='save sequence'):
        validate_density_saves(fixed,output.replace('density set 1','density set 9'))
    with pytest.raises(ValueError,match='selection'):
        validate_density_saves(fixed.replace('set=1','set=9') if reference_first else fixed.replace('set=4','set=9'),output)


def test_native_density_matching_detects_wrong_root_and_is_rotation_invariant(tmp_path):
    # Four occupied spatial orbitals, five electrons, quartet. Reference has
    # its three SOMOs inside the CAS; root 2, not labelled root 1, matches it.
    rng = np.random.default_rng(19)
    q, _ = np.linalg.qr(rng.normal(size=(7, 7)))
    s = np.diag(np.linspace(.7, 1.6, 7))
    hf = np.diag(1/np.sqrt(np.diag(s)))@q
    occ = np.array([2.,1,1,1,0,0,0])
    alternate = np.array([2.,1,1,0,1,0,0])
    a = {'s':s, 'hf':hf, 'cas':hf, 'dhf':hf@np.diag(occ)@hf.T,
         'dr1':hf@np.diag(alternate)@hf.T, 'dr2':hf@np.diag(occ)@hf.T}
    result = compare(a,electrons=5,active_electrons=3,active_orbitals=4,multiplicity=4)
    assert result['closest_density_root'] == 2
    assert result['reference_states'][0]['density_distance'] < 1.e-12
    assert result['rohf_somo_active_projection_eigenvalues'] == pytest.approx([1,1,1])
    # State-average rotations and signs must not look like orbital flips.
    u, _ = np.linalg.qr(rng.normal(size=(4,4)))
    a['cas'] = hf.copy(); a['cas'][:,1:5] = hf[:,1:5]@u
    rotated = compare(a,electrons=5,active_electrons=3,active_orbitals=4,multiplicity=4)
    assert rotated['rohf_somo_active_projection_eigenvalues'] == pytest.approx([1,1,1])
    paths = {}
    for key, value in a.items():
        p = tmp_path/f'job.identity_{key}.out'
        p.write_text('# MATROP\nBEGIN_DATA\n'+'\n'.join(' '.join(f'{v:.16E}' for v in row) for row in value)+'\nEND_DATA\n')
        assert matrix(p) == pytest.approx(value)
        paths[p.name] = file_sha256(p)
    atomic_json(tmp_path/'execution.json',dict(reference_identity_matrices=paths))
    assert report(tmp_path,'job',dict(active_electrons=3,active_orbitals=4),5,4)['closest_density_root'] == 2
    p.write_text('changed')
    with pytest.raises(ValueError,match='changed'): report(tmp_path,'job',dict(active_electrons=3,active_orbitals=4),5,4)
    a['hf'] = a['hf'].T
    with pytest.raises(ValueError,match='orthonormality'): compare(a,electrons=5,active_electrons=3,active_orbitals=4,multiplicity=4)
