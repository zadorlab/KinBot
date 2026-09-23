"""Rehydrate portable, read-only methanol observations for contract tests."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from kinbot import frequencies, symmetry
from kinbot.calculation import geometry_reference
from kinbot.hindered_rotors import HIR
from kinbot.stationary_pt import StationaryPoint


def methanol_data(run='MeOH_H_rotor'):
    return json.loads((Path(__file__).parent / 'reference/methanol_counting.json').read_text())['runs'][run]


def saved_point(data):
    p = StationaryPoint(data['source_job'], data['charge'], data['multiplicity'],
                        atom=data['atoms'], geom=np.array(data['geometry_angstrom']),
                        wellorts=data['wellorts'])
    # characterize(bond_mx=...) skips bond perception; it does not assign the
    # supplied matrix. Restore the saved TS union before deriving atom IDs.
    p.bond = np.array(data['bond'])
    p.bond01 = (p.bond > 0).astype(int)
    p.bonds = [np.array(b) for b in data['bonds']]
    p.rads = [np.array(r) for r in data['rads']]
    if 'reac_bond' in data:
        p.reac_bond = np.array(data['reac_bond'])
    p.find_cycle()
    p.characterize(bond_mx=p.bond)
    p.dihed = copy.deepcopy(data['dihedrals'])
    p.source_job, p.source_row_id = data['source_job'], data['source_row_id']
    p.energy, p.zpe = data['electronic_energy_hartree'], data['zpe_hartree']
    p.freq = list(data['raw_frequencies_cm-1'])
    p.reduced_freqs = list(p.freq)
    p.hess = data.get('hessian_native_eV_angstrom-2', data.get('hessian_hartree_bohr-2', []))
    p.conformer_representation = 'selected conformer'
    if 'optical_reference' in data:
        p.optical_reference = copy.deepcopy(data['optical_reference'])
    symmetry.calculate_symmetry(p)
    if data['scans']:
        h = HIR(p, SimpleNamespace(qc='fc'),
                {'nrotation': 12, 'plot_hir_profiles': False, 'rotor_0_test': 1})
        p.hir = h
        # Fixture provenance is explicit: the extraction verified the saved
        # selected row, lifecycle rotor and actual reference scan geometry.
        h.scan_reference = geometry_reference(p)
        h.scan_reference['backend'] = h.qc.qc
        h.scan_reference['dihedrals'] = copy.deepcopy(p.dihed)
        h.scan_reference['provenance'] = 'saved selected row and lifecycle input, reconstructed for this fixture only'
        h.rigid_scan = False
        for index, scan in enumerate(data['scans']):
            h.hir_status.append([0 if point.get('accepted', point['status'] == 'normal')
                                 else 1 for point in scan['points']])
            h.hir_energies.append([point['electronic_energy_hartree'] for point in scan['points']])
            h.hir_geoms.append([point['geometry_angstrom'] for point in scan['points']])
            h.scan_jobs.append([point['source_job'] for point in scan['points']])
            h.point_observations.append(copy.deepcopy(scan['points']))
            assert h.fourier_fit('fixture', np.arange(12) * np.pi / 6, index)
        _, p.reduced_freqs = frequencies.get_frequencies(p, p.hess, p.geom)
    return p
