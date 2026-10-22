"""Opt-in native ROHF/CASSCF reference diagnostics for finite VRC probes.

These compare one-particle densities and orbital subspaces, not many-electron
wavefunction overlaps. They cannot certify a correlated state or a full scan.
"""
import json
from pathlib import Path
import re

import numpy as np

from kinbot.execution import file_sha256


def state_records(model, multiplicity):
    roots = next(e['roots'] for e in model['state_average']
                 if e['multiplicity'] == multiplicity)
    # MULTI DM saves all individual state densities to ONE explicitly named
    # ORBITAL record. NATORB records are not separate per-root density stores.
    return {r: '3000.2' for r in range(1, roots + 1)}


def natural_cards(model, multiplicity):
    """Save individual AO densities without changing SA orbitals/weights.

    Explicit ORBITAL avoids dependence on which NATORB/CANONICAL card occurs
    first. Keep the production SA canonical/natural exports on 2141/2140.
    """
    state_records(model, multiplicity)  # Require the reference spin manifold.
    if sum(e['roots'] for e in model['state_average']) + 1 > 50:
        raise ValueError('MULTI supports at most 50 saved density matrices.')
    return 'orbital,3000.2\ndm\n'


def density_sets(model, multiplicity):
    """C1 MULTI DM sets in declared WF/STATE order, before the SA density."""
    offset = 0
    for entry in model['state_average']:
        if entry['multiplicity'] == multiplicity:
            return {root: offset+root for root in range(1, entry['roots']+1)}
        offset += entry['roots']
    raise ValueError('Reference multiplicity absent from the state average.')


def export_deck(stem, model, multiplicity, *, part='all'):
    """Export MR immediately after MULTI; export HF immediately after ROHF.

    Mixed-spin MULTI in Molpro 2024 can save the individual densities but fail
    their MS2 lookup when the requested spin is not the last WF block. Select
    the explicitly saved set instead. validate_density_saves checks the native
    save sequence; an unknown sequence is never silently accepted.
    """
    matrices, comments = [], ''
    if part in ('all', 'mr'):
        layout = ','.join(f"{e['multiplicity']}:{e['roots']}" for e in model['state_average'])
        comments = f'! KinBot reference density layout: {layout}; reference={multiplicity}\n'
        matrices = [('s', 's'), ('cas', 'orb,2140.2,natural')]
        matrices += [(f'dr{root}', f'den,3000.2,set={iset}')
                     for root, iset in density_sets(model, multiplicity).items()]
    if part in ('all', 'hf'):
        matrices += [('hf', 'orb,2150.2'), ('dhf', 'den,2150.2')]
    if not matrices:
        raise ValueError('Unknown reference export part.')
    return comments + '{matrop\n' + ''.join(
        f'load,{name},{source}\nwrite,{name},{stem}.identity_{name}.out,new,scientific\n'
        for name, source in matrices) + '}\n'


def validate_density_saves(native, output):
    """Bind numbered density selection to the observed MULTI save sequence."""
    contract = re.findall(r'! KinBot reference density layout: ([\d:,]+); reference=(\d+)', native)
    if not contract:
        return None  # Old completed attempts keep their original evidence.
    if len(contract) != 1:
        raise ValueError('Ambiguous reference density layout.')
    layout, reference = contract[0]
    entries = [dict(multiplicity=int(m), roots=int(n)) for m, n in
               (pair.split(':') for pair in layout.split(','))]
    blocks = re.split(r'(?im)^\s*PROGRAM\s*\*\s*', output)
    multi = [b for b in blocks if re.match(r'MULTI\b', b, re.I)]
    if len(multi) != 1:
        raise ValueError('Cannot verify a unique native MULTI density save.')
    saves = re.findall(r'First order charge density matrix for state\s+(\d+\.\d+)\s+'
                       r'saved on record\s+3000\.2\s+\(density set\s+(\d+)\)', multi[0], re.I)
    expected = [(f'{root}.1', str(i)) for i, root in enumerate(
        (r for e in entries for r in range(1, e['roots']+1)), 1)]
    if saves != expected:
        raise ValueError('Native MULTI density save sequence differs from WF/STATE layout.')
    sets = density_sets({'state_average': entries}, int(reference))
    actual = re.findall(r'(?im)^load,dr(\d+),den,3000\.2,set=(\d+)\s*$', native)
    if actual != [(str(r), str(s)) for r, s in sets.items()]:
        raise ValueError('Native density selection differs from declared reference roots.')
    return dict(layout=entries, reference_multiplicity=int(reference),
                root_sets=sets, native_save_sequence_verified=True)


def repair_legacy_selectors(text):
    """Repair both generations of our failed per-root density export decks.

    Keep prepared inputs immutable; the executed input and receipt record the
    correction. General user density selectors are never rewritten.
    """
    repairs = []
    pattern = r'(?im)^load,dr(\d+),den,(\d+\.2),state=(\d+\.1),(?:spin|ms2)=(\d+)[ \t]*$'
    for match in list(re.finditer(pattern, text)):
        root, old_record, state, ms2 = match.groups()
        card = f'natorb,{old_record},state={state},spin={ms2}'
        blocks = [b for b in re.finditer(r'(?is)\{multi\b[^{}]*\}', text)
                  if re.search(r'(?im)^'+re.escape(card)+r'[ \t]*$', b[0])]
        if old_record != f'{3000+int(root)}.2' or state != root+'.1' or len(blocks) != 1:
            continue
        block = blocks[0][0]
        if re.search(r'(?im)^orbital,', block) and 'orbital,3000.2\ndm\n' not in block:
            raise ValueError('Unexpected diagnostic MULTI output record; inspect before repair.')
        replacement = '' if 'orbital,3000.2\ndm\n' in block else 'orbital,3000.2\ndm'
        updated = re.sub(r'(?im)^'+re.escape(card)+r'[ \t]*$', replacement, block)
        text = text.replace(block, updated, 1)
        text = text.replace(match[0], f'load,dr{root},den,3000.2,state={state},ms2={ms2}', 1)
        repairs.append(dict(record='3000.2', previous_record=old_record,
                            state=state, ms2=int(ms2), saved_individual_densities=True))
    # Migrate both old per-root decks (now normalized above) and the preceding
    # DM/MS2 generation. Only our complete diagnostic MATROP block is moved.
    multis = list(re.finditer(r'(?is)\{multi\b[^{}]*\}', text))
    exports = [b for b in re.finditer(r'(?is)\{matrop\b[^{}]*\}', text)
               if re.search(r'(?im)^load,dr\d+,den,3000\.2,state=', b[0])]
    if not exports:
        return text, repairs
    if len(multis) != 1 or len(exports) != 1:
        raise ValueError('Cannot identify a unique diagnostic MULTI/export block.')
    block, multi = exports[0][0], multis[0][0]
    if 'orbital,3000.2\ndm' not in multi:
        raise ValueError('Individual MULTI density saving is missing.')
    entries = [dict(multiplicity=int(spin)+1, roots=int(roots)) for spin, roots in
               re.findall(r'(?im)^wf,\d+,1,(\d+)\s*\nstate,(\d+)\s*$', multi)]
    selectors = re.findall(r'(?im)^load,dr(\d+),den,3000\.2,state=(\d+)\.1,ms2=(\d+)\s*$', block)
    stems = set(re.findall(r'(?im)^write,\w+,([^,\n]+)\.identity_\w+\.out,new,scientific\s*$', block))
    if not entries or not selectors or len(stems) != 1 or len({s[2] for s in selectors}) != 1:
        raise ValueError('Cannot bind diagnostic densities to the MULTI state layout.')
    mult = int(selectors[0][2])+1
    model = {'state_average': entries}
    if [(int(a), int(b)) for a,b,_ in selectors] != [(r,r) for r in density_sets(model,mult)]:
        raise ValueError('Diagnostic density roots differ from the MULTI reference manifold.')
    stem = stems.pop()
    # Refuse unrelated operations in this block instead of dropping them.
    lines = [s.strip() for s in block.splitlines()[1:-1] if s.strip()]
    if any(not re.fullmatch(r'(?:load,(?:s,s|cas,orb,2140\.2,natural|hf,orb,2150\.2|dhf,den,2150\.2|'
                           r'dr\d+,den,3000\.2,state=\d+\.1,ms2=\d+)|'
                           r'write,(?:s|cas|hf|dhf|dr\d+),'+re.escape(stem)+
                           r'\.identity_\w+\.out,new,scientific)', s) for s in lines):
        raise ValueError('Unexpected operation in diagnostic export block.')
    text = text.replace(block, export_deck(stem,model,mult,part='hf').rstrip(), 1)
    # Preserve ENERGY arrays before any diagnostic command, at the end of the
    # generated root-marker sequence and before RS2/MRCI/ROHF can reuse records.
    markers = list(re.finditer(r'(?im)^show,reconr\d+\s*$', text))
    if not markers or markers[-1].end() > text.index('{matrop'):
        raise ValueError('Cannot locate the post-MULTI root-energy boundary.')
    at = markers[-1].end()
    text = text[:at]+'\n'+export_deck(stem,model,mult,part='mr')+text[at:]
    repairs.append(dict(numbered_density_sets=True, exported_before_correlation=True,
                        reference_multiplicity=mult, root_sets=density_sets(model,mult)))
    return text, repairs


def matrix(path):
    text = Path(path).read_text()
    blocks = re.findall(r'(?is)BEGIN_DATA\s*(.*?)\s*END_DATA', text)
    if len(blocks) != 1:
        raise ValueError(f'Expected one MATROP matrix: {path}')
    rows = [line for line in blocks[0].splitlines()
            if line.strip() and not line.lstrip().startswith(('#', '*', '!'))]
    values = [float(v.replace('D', 'E').replace('d', 'e'))
              for line in rows for v in line.replace(',', ' ').split()]
    size = int(np.sqrt(len(values)))
    if not size or size*size != len(values) or not np.isfinite(values).all():
        raise ValueError(f'Invalid square MATROP matrix: {path}')
    return np.array(values).reshape(size, size)


def compare(matrices, *, electrons, active_electrons, active_orbitals, multiplicity):
    s, cas, hf, dhf = [matrices[k] for k in ('s', 'cas', 'hf', 'dhf')]
    if any(v.shape != s.shape for v in matrices.values()):
        raise ValueError('Native reference matrices have different AO dimensions.')
    # Detect transpose/format/basis mistakes rather than assigning states from
    # a numerically meaningless overlap. All probes use C1 and one AO basis.
    orth = max(float(np.max(abs(c.T@s@c-np.eye(len(s))))) for c in (cas, hf))
    if orth > 2.e-4:
        raise ValueError(f'Native MO orthonormality failed ({orth:g}).')
    nc = (electrons-active_electrons)//2
    nh = (electrons-multiplicity+1)//2
    somo = slice(nh, nh+multiplicity-1)
    overlap = hf.T@s@cas
    singular = np.linalg.svd(overlap[somo, nc:nc+active_orbitals], compute_uv=False)
    hfd = hf.T@s@dhf@s@hf
    expected = np.zeros(len(s)); expected[:nh] = 2.; expected[somo] = 1.
    if (abs(np.trace(hfd)-electrons) > 2.e-3
            or np.max(abs(hfd-np.diag(expected))) > 2.e-3):
        raise ValueError('Exported ROHF density does not match its declared occupation.')
    states = []
    for key in sorted(k for k in matrices if re.fullmatch(r'dr\d+', k)):
        density = hf.T@s@matrices[key]@s@hf
        if (abs(np.trace(density)-electrons) > 2.e-3
                or np.max(abs(density-density.T)) > 2.e-4
                or np.linalg.eigvalsh(density).min() < -2.e-3
                or np.linalg.eigvalsh(density).max() > 2.002):
            raise ValueError(f'Invalid state-specific density: {key}')
        states.append(dict(root=int(key[2:]),
            density_distance=float(np.linalg.norm(density-hfd)),
            rohf_somo_occupations=np.diag(density)[somo].tolist()))
    if not states:
        raise ValueError('No state-specific MR reference densities.')
    states.sort(key=lambda v: v['density_distance'])
    return dict(mo_orthonormality_error=orth,
        rohf_somo_active_projection_eigenvalues=(singular**2).tolist(),
        reference_states=states, closest_density_root=states[0]['root'],
        verdict='review_required',
        limitation='Density similarity is not a many-electron overlap; degenerate states can mix. '
                   'CC amplitudes/diagnostics and MR CI character must also be reviewed.')


def report(directory, stem, space, electrons, multiplicity):
    directory = Path(directory)
    receipt = json.loads((directory/'execution.json').read_text())
    files = receipt.get('reference_identity_matrices', {})
    if not files:
        return {'verdict': 'unavailable', 'reason': 'No native reference matrices recorded.'}
    matrices = {}
    for name, digest in files.items():
        path = directory/name
        if path.parent != directory or file_sha256(path) != digest:
            raise ValueError('Native reference identity export changed: '+name)
        key = name.removeprefix(stem+'.identity_').removesuffix('.out')
        matrices[key] = matrix(path)
    return compare(matrices, electrons=electrons, active_electrons=space['active_electrons'],
                   active_orbitals=space['active_orbitals'], multiplicity=multiplicity)
