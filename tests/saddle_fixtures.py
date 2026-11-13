"""Synthetic mode evidence for symbolic MESS mirror/ensemble fixtures."""
import re
import numpy as np
from kinbot.mess_networks import split_model, _HEADER
from kinbot.saddle_identity import annotate


def bind_mode(species):
    """Give a synthetic saddle one documented unstable breathing coordinate."""
    from kinbot import constants
    from kinbot.calculation import array_fingerprint
    mass = np.repeat([constants.exact_mass[str(s)] for s in species.atom],3)
    mode = (np.asarray(species.geom)-np.asarray(species.geom).mean(axis=0)).ravel()*np.sqrt(mass)
    mode /= np.linalg.norm(mode)
    species.hess = (.01*np.eye(len(mode))-.1*np.outer(mode,mode))*np.sqrt(np.outer(mass,mass))
    species.source_job = species.name
    species.source_row_id = 1
    species.optical_hessian_reference = dict(source_job=species.name,source_row_id=1,
        atoms=list(map(str,species.atom)),geometry_sha256=array_fingerprint(species.geom),
        hessian_sha256=array_fingerprint(species.hess),hessian_massweighted=False,
        hessian_unit='hartree / bohr^2')


def bind_member_modes(species):
    from copy import copy
    from dataclasses import replace
    from kinbot.conformer_records import retain
    records = []
    for record in species.conformer_inventory:
        point = copy(species)
        point.name = record.source_job
        point.geom = np.asarray(record.geometry)
        bind_mode(point)
        records.append(replace(record, electronic_energy_hartree=record.zero_energy_hartree-species.zpe,
            zpe_hartree=species.zpe, hessian=point.hess.tolist(),
            hessian_reference=point.optical_hessian_reference))
    retain(species,records,[r.index for r in records])


def add_saddle_evidence(text):
    prefix, blocks, end = split_model(text)
    result = []
    for block in blocks:
        header = _HEADER.search(block)
        if header[1] == 'Barrier':
            energy = re.findall(r'ZeroEnergy\[kcal/mol\]\s+([-+\d.eE]+)', block)
            xyz = re.search(r'Geometry\[angstrom\]\s+(\d+)\s*\n', block)
            if xyz:
                rows = [line.split() for line in block[xyz.end():].splitlines()[:int(xyz[1])]]
                symbols = [row[0] for row in rows]
                positions = np.array([row[1:4] for row in rows], float)
            else:
                symbols = ['C','H','F','Cl']
                positions = np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
                edge = re.search(r'! kinbot_mirror_endpoints (\S+) (\S+) (\S+) (\S+)', block)
                if edge and sorted((edge[1],edge[3])) > sorted((edge[2],edge[4])):
                    positions *= [-1,1,1]
            mode = positions-positions.mean(axis=0)
            if not np.linalg.norm(mode):
                mode[:] = 1.
            record = dict(symbols=symbols, positions=positions.tolist(), mode=mode.tolist(),
                          isotopes=[0]*len(symbols), charge=0, multiplicity=1,
                          frequencies=[-100.,200.,300.], energy_hartree=-100.+float(energy[0])/627.509
                          if energy else -100., source_job=header[2].split()[0])
            block = annotate(block,record)
        result.append(block)
    return prefix+''.join(result)+end
