"""Reproduce same-geometry orbital comparisons from the transferred native XML.

Usage: python compare_orbitals.py /path/to/extracted/archive > comparison.json
Read-only; this neither edits orbital files nor accepts a production model.
"""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

NS = {'m': 'http://www.molpro.net/schema/molpro-output',
      'c': 'http://www.xml-cml.org/schema'}


def read(path):
    root = ET.parse(path).getroot()
    sets = root.findall('.//m:orbitals', NS)
    if len(sets) != 1:
        raise ValueError(f'Expected exactly one orbital set: {path}')
    orbitals = sets[0]
    matrix = np.array([[float(v) for v in o.text.split()] for o in orbitals])
    basis = root.find('.//m:basisSet', NS)
    atoms = [dict(a.attrib) for a in root.findall('.//c:atom', NS)]
    if len(matrix) != int(basis.attrib['length']):
        raise ValueError(f'Complete MO basis required: {path}')
    molecule_id = root.find('m:molecule', NS).attrib['id']
    # Per-job XML molecule IDs are bookkeeping. Retain the full basis-to-atom
    # associations, normalizing only that verified local molecule identifier.
    basis_definition = ET.tostring(basis, encoding='unicode').replace(
        "molecule[@id='"+molecule_id+"']", "molecule[@id='CURRENT']")
    return matrix, atoms, basis_definition


def compare(root, left, right, indices):
    a, atoms_a, basis_a = read(root/left)
    b, atoms_b, basis_b = read(root/right)
    if atoms_a != atoms_b or basis_a != basis_b or a.shape != b.shape:
        raise ValueError('Comparison requires identical atom order, geometry and AO basis.')
    # Both complete MO sets are orthonormal in the same AO metric. Therefore
    # b = T a implies T_ij = <b_i|a_j>. The rectangular Cartesian representation
    # has redundant components; its pseudoinverse gives this unique MO transform.
    transform = b @ np.linalg.pinv(a)
    residual = float(np.linalg.norm(transform @ a-b)/np.linalg.norm(b))
    unitary = float(np.max(abs(transform @ transform.T-np.eye(len(transform)))))
    if residual > 1.e-9 or unitary > 1.e-8:
        raise ValueError('XML coefficient sets do not give an orthonormal MO transform.')
    result = dict(left=str(left), right=str(right),
                  sha256={str(p): hashlib.sha256((root/p).read_bytes()).hexdigest()
                          for p in (left, right)},
                  relative_reconstruction_error=residual, max_orthogonality_error=unitary,
                  subspaces={}, right_orbitals_in_left_basis={})
    for label, (i, j) in indices.items():
        block = transform[np.ix_(np.array(i)-1, np.array(j)-1)]
        result['subspaces'][label] = dict(right_indices=i, left_indices=j,
            overlap_singular_values=np.linalg.svd(block, compute_uv=False).tolist())
    for i in range(17, 23):
        order = np.argsort(abs(transform[i-1]))[-6:][::-1]
        result['right_orbitals_in_left_basis'][str(i)] = [
            dict(index=int(j)+1, squared_overlap=float(transform[i-1, j]**2)) for j in order]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('extracted_archive', type=Path)
    args = parser.parse_args()
    root = args.extracted_archive/'production/orbital_recon'
    results = []
    for b in ('adz', 'atz'):
        left = Path(f'compact_resolution_20261025/csh2_sh_exit_resolution/'
                    f'csh2_sh_exit_resolution_union/{b}__exit_l2/{b}__exit_l2.natural.orbitals.xml')
        right = Path(f'compact_continuation_20261025/csh2_sh_exit_resolution/'
                     f'csh2_sh_exit_resolution_inward/{b}__in_near/{b}__in_near.natural.orbitals.xml')
        results.append(compare(root, left, right, {
            'active_five': (list(range(18, 23)), list(range(18, 23))),
            'include_inactive_17': (list(range(17, 23)), list(range(17, 23)))}))
    print(json.dumps(dict(method='Complete same-geometry orthonormal MO change of basis; '
                                'no occupation cutoff, AO-population proxy or cross-geometry overlap.',
                          comparisons=results), indent=2))


if __name__ == '__main__':
    main()
