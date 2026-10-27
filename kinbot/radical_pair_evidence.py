"""Native orbital and CAS(2,2) evidence for the restricted radical-pair policy.

MO transforms below are only defined for identical nuclei and AO bases. They
are not cross-geometry overlaps, and no AO coefficient is called a population.
"""
from dataclasses import dataclass
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np

NS = {'m': 'http://www.molpro.net/schema/molpro-output',
      'c': 'http://www.xml-cml.org/schema'}
XLINK = '{http://www.w3.org/1999/xlink}href'


def numbers(text):
    return [float(v.replace('D', 'E').replace('d', 'e')) for v in text.split()]


@dataclass
class Orbitals:
    coefficients: np.ndarray
    occupations: np.ndarray
    symbols: list
    positions: np.ndarray
    atom_basis: list
    atom_columns: list


def read_orbitals(path):
    root = ET.parse(path).getroot()
    sets = root.findall('.//m:orbitals', NS)
    bases = root.findall('.//m:basisSet', NS)
    atoms = root.findall('.//c:atom', NS)
    if len(sets) != 1 or len(bases) != 1 or not atoms:
        raise ValueError('One complete native orbital set and AO basis are required: '+str(path))
    rows = list(sets[0]); basis = bases[0]
    c = np.array([numbers(row.text or '') for row in rows])
    n = np.array([float(row.get('occupation', 'nan')) for row in rows])
    if (c.ndim != 2 or len(c) != int(basis.get('length')) or
            c.shape[1] != int(basis.get('cartesianLength')) or
            not np.isfinite(c).all() or not np.isfinite(n).all()):
        raise ValueError('Incomplete/nonfinite native MO coefficient matrix.')
    groups = {g.get('id'): g for g in basis.findall('m:basisGroup', NS)}
    assigned = {a.get('id'): [] for a in atoms}
    for assoc in basis.findall('m:association', NS):
        ag = re.findall(r"@id='([^']+)'", assoc.find('m:atoms', NS).get(XLINK))
        bg = re.findall(r"@id='([^']+)'", assoc.find('m:bases', NS).get(XLINK))
        ag = [v for v in ag if v in assigned]; bg = [v for v in bg if v in groups]
        if not ag or not bg:
            raise ValueError('Unrecognized native AO-to-atom association.')
        for atom in ag:
            assigned[atom].extend(bg)
    atom_basis, atom_columns, offset = [], [], 0
    for atom in atoms:
        descriptions = []; count = 0
        for key in assigned[atom.get('id')]:
            g = groups[key]; lo, hi = int(g.get('minL')), int(g.get('maxL'))
            if lo != hi:
                raise ValueError('Mixed angular-momentum basis groups are unsupported.')
            exponents = numbers(g.find('m:basisExponents', NS).text)
            contractions = [numbers(x.text) for x in g.findall('m:basisContraction', NS)]
            if not contractions or any(len(v) != len(exponents) for v in contractions):
                raise ValueError('Incomplete native contraction definition.')
            descriptions.append((lo, tuple(exponents), tuple(map(tuple, contractions))))
            count += len(contractions)*(lo+1)*(lo+2)//2
        atom_basis.append(descriptions)
        atom_columns.append(list(range(offset, offset+count))); offset += count
    if offset != c.shape[1]:
        raise ValueError('Native atom/basis associations do not span the Cartesian coefficients.')
    xyz = np.array([[float(a.get(k)) for k in ('x3', 'y3', 'z3')] for a in atoms])
    return Orbitals(c, n, [a.get('elementType') for a in atoms], xyz, atom_basis, atom_columns)


def transform(left, right):
    """Return <right_i|left_j> from two complete, same-geometry MO sets."""
    if (left.symbols != right.symbols or left.atom_basis != right.atom_basis or
            left.coefficients.shape != right.coefficients.shape or
            not np.allclose(left.positions, right.positions, rtol=0, atol=2.e-8)):
        raise ValueError('MO comparison requires the same atom order, geometry and AO basis.')
    a, b = left.coefficients, right.coefficients
    t = b @ np.linalg.pinv(a, rcond=1.e-12)
    if (np.linalg.norm(t@a-b)/max(np.linalg.norm(b), 1.e-30) > 1.e-7 or
            np.max(abs(t@t.T-np.eye(len(t)))) > 2.e-5):
        raise ValueError('Native MO change of basis failed reconstruction/orthonormality.')
    return t


def projection(left, right, left_indices, right_indices):
    block = transform(left, right)[np.ix_(right_indices, left_indices)]
    return (np.linalg.svd(block, compute_uv=False)**2).tolist()


def fragment_embedding(dimer, fragments, mappings):
    """Embed independently computed fragment MOs in the identical dimer AO basis."""
    vectors, somo = [], []
    for frag, mapping in zip(fragments, mappings):
        if len(mapping) != len(frag.symbols):
            raise ValueError('Fragment atom map has the wrong length.')
        embedded = np.zeros((len(frag.coefficients), dimer.coefficients.shape[1]))
        for i, j in enumerate(mapping):
            if (frag.symbols[i] != dimer.symbols[j] or frag.atom_basis[i] != dimer.atom_basis[j]
                    or not np.allclose(frag.positions[i], dimer.positions[j], rtol=0, atol=2.e-8)):
                raise ValueError('Independent fragment and dimer AO basis/geometry differ.')
            embedded[:, dimer.atom_columns[j]] = frag.coefficients[:, frag.atom_columns[i]]
        ids = np.flatnonzero(abs(frag.occupations-1.) < 1.e-5)
        if len(ids) != 1:
            raise ValueError('Each neutral doublet must have exactly one native ROHF SOMO.')
        t = embedded @ np.linalg.pinv(dimer.coefficients, rcond=1.e-12)
        if np.linalg.norm(t@dimer.coefficients-embedded)/np.linalg.norm(embedded) > 1.e-7:
            raise ValueError('Fragment orbitals do not lie in the dimer physical AO subspace.')
        if np.max(abs(t@t.T-np.eye(len(t)))) > 2.e-5:
            raise ValueError('Embedded fragment MOs are not orthonormal in the dimer metric.')
        vectors.append(t); somo.append(t[ids[0]])
    s = np.array(somo)
    if np.max(abs(s@s.T-np.eye(2))) > 2.e-5:
        raise ValueError('Separated fragment SOMOs have appreciable interfragment overlap.')
    populations = [float(np.sum(t*t*dimer.occupations[None, :])) for t in vectors]
    return s, populations


def natural_ci(text, multiplicity):
    """Read the transformed single-state CAS(2,2) CI vector, not a later MRCI table."""
    sections = re.split(r'(?im)^\s*PROGRAM \*', text)
    multis = [s for s in sections if re.match(r'\s*MULTI\b', s)]
    if len(multis) != 1:
        raise ValueError('Expected exactly one native MULTI section.')
    tail = re.split(r'(?i)Natural orbital dump[^\n]*2140\.2[^\n]*', multis[0])
    if len(tail) != 2 or 'Transform the CI vector' not in tail[1]:
        raise ValueError('Missing CI vector explicitly transformed to the saved natural orbitals.')
    block = tail[1].split('Energy:', 1)[0]
    entries = re.findall(r'^\s*([20ab]{2})\s+([-+0-9.EeDd]+)\s*$', block, re.M)
    ci = {k: float(v.replace('D', 'E').replace('d', 'e')) for k, v in entries}
    allowed = {'20', '02', 'ab'} if multiplicity == 1 else {'aa'}
    if not ci or len(ci) != len(entries) or set(ci)-allowed:
        raise ValueError('Unrecognized single-state CAS(2,2) spin-adapted CI expansion.')
    norm = sum(v*v for v in ci.values())
    if abs(norm-1.) > 2.e-5:
        raise ValueError('Printed natural-orbital CI vector is incomplete or unnormalised.')
    return ci


def singlet_matrix(ci):
    return np.array([[ci.get('20', 0.), ci.get('ab', 0.)/np.sqrt(2)],
                     [ci.get('ab', 0.)/np.sqrt(2), ci.get('02', 0.)]])


def ionic_probability(ci, fragment_to_active):
    """Double occupancy of either independent SOMO at the separated limit."""
    w = fragment_to_active @ singlet_matrix(ci) @ fragment_to_active.T
    return float(w[0, 0]**2+w[1, 1]**2)


def ci_overlap(left_ci, right_ci, active_transform, multiplicity):
    if multiplicity == 3:
        return float(np.linalg.det(active_transform)**2)
    a = singlet_matrix(left_ci)
    b = active_transform.T @ singlet_matrix(right_ci) @ active_transform
    return float(np.sum(a*b)**2)


def wavefunction_overlap(left_ci, right_ci, mo_transform, closed, multiplicity):
    """Exact determinant overlap for two real CAS(2,2) states, including inactive MOs.

    MO transform is <right|left>. Singlet CSFs expand in alpha/beta ordered
    determinants; the two open-shell determinants have the same amplitude.
    This includes closed-active mixing instead of multiplying isolated overlaps.
    """
    def determinants(ci):
        core=list(range(closed));a,b=closed,closed+1
        if multiplicity==3:
            return [(core+[a,b],core,ci['aa'])]
        return [(core+[a],core+[a],ci.get('20',0.)),
                (core+[b],core+[b],ci.get('02',0.)),
                (core+[a],core+[b],ci.get('ab',0.)/np.sqrt(2)),
                (core+[b],core+[a],ci.get('ab',0.)/np.sqrt(2))]
    total=0.
    for la,lb,lc in determinants(left_ci):
        for ra,rb,rc in determinants(right_ci):
            total+=lc*rc*np.linalg.det(mo_transform[np.ix_(ra,la)])*np.linalg.det(mo_transform[np.ix_(rb,lb)])
    return float(total*total)
