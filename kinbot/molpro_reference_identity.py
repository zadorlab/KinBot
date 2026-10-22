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
    return {r: f'{3000 + r}.2' for r in range(1, roots + 1)}


def natural_cards(model, multiplicity):
    # MULTI NATORB's SPIN is 2S; DENSITY/LOAD's SPIN is S. Do not share
    # that keyword between these two interfaces (Molpro manual).
    return ''.join(f'natorb,{record},state={root}.1,spin={multiplicity-1}\n'
                   for root, record in state_records(model, multiplicity).items())


def export_deck(stem, model, multiplicity):
    """Run immediately after ROHF, before CC can replace its density."""
    matrices = [('s', 's'), ('cas', 'orb,2140.2,natural'),
                ('hf', 'orb,2150.2'), ('dhf', 'den,2150.2')]
    matrices += [(f'dr{root}', f'den,{record},state={root}.1,ms2={multiplicity-1}')
                 for root, record in state_records(model, multiplicity).items()]
    return '{matrop\n' + ''.join(
        f'load,{name},{source}\nwrite,{name},{stem}.identity_{name}.out,new,scientific\n'
        for name, source in matrices) + '}\n'


def repair_legacy_selectors(text):
    """Repair only our paired NATORB/identity-export cards, not user densities.

    Keep prepared inputs immutable; the executed input and receipt record the
    correction. NATORB still uses its documented SPIN=2S convention.
    """
    repairs = []
    pattern = r'(?im)^load,dr(\d+),den,(\d+\.2),state=(\d+\.1),spin=(\d+)\s*$'
    def replace(match):
        root, record, state, ms2 = match.groups()
        if (record != f'{3000+int(root)}.2' or state != root+'.1'
                or not re.search(r'(?im)^natorb,'+re.escape(record)+',state='
                                 +re.escape(state)+',spin='+ms2+r'\s*$', text)):
            return match[0]
        repairs.append(dict(record=record, state=state, ms2=int(ms2)))
        return f'load,dr{root},den,{record},state={state},ms2={ms2}'
    return re.sub(pattern, replace, text), repairs


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
