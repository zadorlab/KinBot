"""Bounded response to repeated Cartesian minimum steps at Sella's radius floor."""
from collections import deque
from pathlib import Path
import json

import numpy as np

from kinbot.execution import atomic_json


class OptimizerCycleError(RuntimeError):
    """The evaluated minimum remains cycling after bounded radius reductions."""


class FloorCycle:
    def __init__(self, optimizer, fmax, report):
        self.optimizer = optimizer
        self.fmax = fmax
        self.report = Path(report) if report else None
        self.points = deque(maxlen=12)
        self.events = []
        self.last_step = None
        if self.report and self.report.is_file():
            saved = json.loads(self.report.read_text())
            self.events = saved['events']
            if self.events:
                # A resumed worker inherits the exhausted numerical settings.
                # A fresh seed may still converge, but it receives no new
                # allowance to repeat the same radius reductions.
                optimizer.delta_min = min(optimizer.delta_min, self.events[-1]['new_floor'])

    def __call__(self):
        opt = self.optimizer
        if opt.nsteps == self.last_step:
            return
        self.last_step = opt.nsteps
        # Use Sella's already evaluated point. Reading atoms could evaluate a
        # different last Hessian probe or trigger additional QC.
        current = opt.pes.curr
        if current.get('g') is None or current.get('x') is None:
            return
        force = float(np.linalg.norm(np.asarray(current['g']).reshape(-1, 3), axis=1).max())
        if (force <= self.fmax or opt.delta > opt.delta_min * 1.01
                or not np.isfinite(force)):
            self.points.clear()
            return
        self.points.append(dict(step=opt.nsteps, x=np.array(current['x'], copy=True),
                                gradient=np.array(current['g'], copy=True), force=force,
                                energy=float(current['f']), rho=float(opt.rho)))
        if len(self.points) < self.points.maxlen:
            return
        points = list(self.points)
        # Require three repeated cycles (period 2--4), persistent unconverged
        # forces and several rejected model predictions. Mere slow descent
        # or step exhaustion is insufficient.
        repeated = any(all(np.linalg.norm(points[i]['gradient'] - points[i-period]['gradient'])
                           <= .05 * np.linalg.norm(points[i]['gradient'])
                           for i in range(period, len(points)))
                       for period in (2, 3, 4))
        reversals = sum(np.dot(a['gradient'], b['gradient']) <
                        -.8 * np.linalg.norm(a['gradient']) * np.linalg.norm(b['gradient'])
                        for a, b in zip(points, points[1:])) >= 3
        poor_model = sum(p['rho'] < 0 for p in points) >= 3
        stalled = min(p['force'] for p in points[6:]) >= .9 * min(p['force'] for p in points[:6])
        if not (repeated and reversals and poor_model and stalled):
            return
        exhausted = sum(e['action'] == 'reduce_radius_floor' for e in self.events) >= 2
        event = dict(old_floor=opt.delta_min, new_floor=opt.delta_min if exhausted else opt.delta_min / 10,
                     action='exhausted' if exhausted else 'reduce_radius_floor',
                     points=[dict(step=p['step'], positions=p['x'].tolist(),
                                  force_eV_angstrom=p['force'],
                                  gradient_eV_angstrom=p['gradient'].tolist(),
                                  energy_eV=p['energy'], rho=p['rho'])
                             for p in points])
        self.events.append(event)
        if self.report:
            atomic_json(self.report, dict(events=self.events, fmax=self.fmax,
                        status=event['action'], hessian_difference_step=opt.eta))
        self.points.clear()
        if exhausted:
            raise OptimizerCycleError('Persistent Cartesian minimum cycling after two '
                                      'trust-radius floor reductions. Inspect optimizer cycle evidence.')
        # The finite-difference Hessian step and convergence threshold retain
        # their independent meanings. Only the lower trust-radius bound moves.
        opt.delta_min = event['new_floor']
        opt.delta = min(opt.delta, opt.delta_min)


def attach_cycle_monitor(optimizer, fmax, report):
    from sella import Sella
    from sella.peswrapper import PES
    if (not isinstance(optimizer, Sella) or optimizer.ord != 0
            or type(optimizer.pes) is not PES or optimizer.atoms.constraints
            or getattr(optimizer, 'user_constraints', None) is not None
            or any(optimizer.pes.cons.internals.get(key)
                   for key in ('bonds', 'angles', 'dihedrals', 'other'))):
        return
    monitor = getattr(optimizer, '_kinbot_floor_cycle', None)
    if monitor is None:
        path = str(report) + '.cycle.json' if report else None
        monitor = FloorCycle(optimizer, fmax, path)
        optimizer.attach(monitor)
        optimizer._kinbot_floor_cycle = monitor
    monitor.fmax = fmax
