"""Import the completed B–W3 saddle into the existing two-direction IRC owner.

Run on the machine with the original dispatch output. Staging is idempotent;
the existing saddle, checkpoint and task state are read-only inputs.
"""
from dataclasses import asdict
from pathlib import Path
import argparse
import json
import re

import numpy as np
from ase.io import read
from ase.units import Hartree
from kinbot import saddle_recheck
from kinbot.anl.results import parse_gaussian_harmonic
from kinbot.execution import atomic_json, file_sha256
from kinbot.gaussian_modes import read_normal_modes
from kinbot.parameters import Parameters
from kinbot.theory import resolve_profiles

HERE=Path(__file__).resolve().parent
SOURCE_HASH='5f2c324dfb57601c41aa6744b0b171ca0144fb25c98c216c16a56066177a677f'
EXPECTED_FREQUENCIES=[-154.9890,71.4266,274.1081,324.5331,443.9839,1003.2018,
                      1015.3656,1091.2340,1501.8277,2720.3157,3104.4891,3195.4147]

def build_manifest(source, *, expected_hash=None):
    source=Path(source).resolve()
    expected_hash=expected_hash or SOURCE_HASH
    log=source/'tasks/saddle/saddle.log'
    if file_sha256(log)!=expected_hash:
        raise ValueError('Saddle log differs from the reviewed native result.')
    receipt=json.loads((source/'tasks/saddle/execution.json').read_text())
    if receipt.get('status')!='executed':
        raise ValueError('The source calculation has not executed successfully.')
    workflow=json.loads((source/'workflow.json').read_text())
    if (workflow['molecule']['charge'],workflow['molecule']['multiplicity'])!=(0,2):
        raise ValueError('The source electronic state differs from the reviewed doublet.')
    text=log.read_text()
    harmonic=parse_gaussian_harmonic(text,method='B2PLYP',basis='cc-pVTZ',
        dispersion='GD3BJ',natoms=6,stationary_point='transition_state')
    if not np.allclose(harmonic['frequencies_cm_inverse'],EXPECTED_FREQUENCIES,atol=.0001,rtol=0):
        raise ValueError('The spectrum differs from the reviewed saddle.')
    atoms=read(log,format='gaussian-out')
    modes=read_normal_modes(log,atoms.symbols,atoms.positions)
    energy=atoms.get_potential_energy()  # ASE reads the double-hybrid E2 total.
    zpe=harmonic['zpe']['hartree']
    totals=re.findall(r'Sum of electronic and zero-point Energies=\s*([-+\d.]+)',text)
    if not totals or abs(energy/Hartree+zpe-float(totals[-1]))>2.e-6:
        raise ValueError('Double-hybrid energy does not match the native thermochemistry sum.')
    parameters=json.loads((HERE/'path_parameters.json').read_text())
    profile=asdict(resolve_profiles(Parameters(str(HERE/'path_parameters.json'),show_warnings=False).par)['l1'])
    profile['optimizer']='gaussian'  # Native Opt=TS, not a new Sella optimization.
    case=dict(name='b_w3',mode='reuse',seed=dict(symbols=list(atoms.symbols),
        positions=atoms.positions.tolist(),charge=0,multiplicity=2),
        native=dict(profile=profile,source_optimizer='Gaussian Opt=(TS,CalcFC,Tight)',
            data=dict(status='normal',energy=energy,zpe=zpe,frequencies=modes['frequencies_cm1']),
            files={suffix:dict(path='tasks/saddle/saddle'+suffix,
                   sha256=file_sha256(source/('tasks/saddle/saddle'+suffix))) for suffix in ('.com','.log')}),
        purpose='Identify both endpoints of the newly established B2PLYP saddle; no endpoint identity is assumed.')
    return dict(schema=1,parameters=parameters,cases=[case],
        source_job_id=receipt.get('slurm_job_id'),source_log_sha256=expected_hash,
        intended_execution='paths-only; compare optimized endpoints with existing positive-frequency minima')

def stage(source, destination):
    source=Path(source).resolve();destination=Path(destination).resolve()
    manifest=build_manifest(source)
    if destination==source or destination.is_relative_to(source):
        raise ValueError('Use a separate sibling path campaign; preserve the completed source.')
    frozen=destination/'campaign.json'
    if frozen.exists() and json.loads(frozen.read_text())!=manifest:
        raise ValueError('Existing path campaign differs; refusing to overwrite it.')
    destination.mkdir(parents=True,exist_ok=True)
    request=destination/'requested_campaign.json'
    if request.exists() and json.loads(request.read_text())!=manifest:
        raise ValueError('Existing request differs; refusing to overwrite it.')
    atomic_json(request,manifest)
    saddle_recheck.stage(request,destination,source)
    print(destination)
    return destination

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('destination',type=Path)
    args=parser.parse_args()
    stage(args.source,args.destination)
