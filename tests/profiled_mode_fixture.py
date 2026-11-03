"""Simulated completed jobs backed by a captured native Gaussian spectrum.

The two job labels deliberately reuse one spectrum to test rotation/routing;
they are not independent L1/L2 calculations or new scientific evidence.
"""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from ase import Atoms

from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.optimize import Optimize
from kinbot.gaussian_modes import read_normal_modes
from kinbot.reaction_path import set_endpoint_populations
from tests.test_native_mess_symmetry import LOG, point


def stage_profiled_modes(directory, preset='b3lyp-b2plyp-anl'):
    """Caller owns cwd; create only disposable test evidence in directory."""
    directory = Path(directory)
    inp = directory/'modes.json'
    inp.write_text(json.dumps(dict(theory_preset=preset, fc_model_path='/test/uma.pt',
        queuing='local', conformer_search=1, multi_conf_tst=0, high_level=1,
        barrier_threshold=50.)))
    par = Parameters(str(inp), show_warnings=False).par
    qc = QuantumChemistry(par)
    p = point('rxn_20')
    native = read_normal_modes(LOG, p.atom, p.geom)
    # An explicitly routed L2 job without the conventional suffix checks that
    # restart uses its recorded backend, not only filename heuristics.
    lower, higher = 'native_lower', 'selected_refinement'
    rot = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    for job, level, geom in ((lower, 'l1', p.geom.copy()),
                             (higher, 'l2', p.geom @ rot + 5.)):
        qc._record(job, level)
        backend = qc.backends[level]
        backend.db.write(Atoms(list(p.atom), positions=geom), name=job,
            data=dict(status='normal', energy=-100., zpe=.01,
                      frequencies=native['frequencies_cm1']))
        (directory/(job+'.log')).write_bytes(LOG.read_bytes())
        # Simulated owning-worker completion, separate from the captured log.
        (directory/(job+'_sella.log')).write_text('simulated completion\ndone\n')
        if backend.frequency_mode != 'native_hessian':
            mode = native['modes'][int(np.argmin(native['frequencies_cm1']))]
            (directory/(job+'_vib.xyz')).write_text(str(p.natom)+'\nsimulated ASE artifact\n'+''.join(
                str(s)+' '+' '.join(map(str, [*g, *v]))+'\n'
                for s, g, v in zip(p.atom, p.geom, mode)))
    return SimpleNamespace(par=par, lower=lower, higher=higher, symbols=list(p.atom))


def compare_profiled_modes(fixture):
    """Use fresh router, real result readers, acceptance and record adoption."""
    qc = QuantumChemistry(fixture.par)
    ts = point('rxn_20')
    # The captured MESS header identifies w_8 -> w_3. TS populations come
    # from the endpoints, not from treating a saddle graph as a stable well.
    set_endpoint_populations(ts, [point('w_8')], [point('w_3')])
    opt = Optimize(ts, fixture.par, qc, just_high=True)
    opt.l1_reference_job = fixture.lower
    opt.compare_structures(job=fixture.higher)
    return opt
