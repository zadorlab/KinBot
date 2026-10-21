"""Bounded ordinary-worker continuation from a verified evaluated geometry.

Only unconstrained Sella/BFGS workers with literal initial Atoms qualify.
No frequency rejection, SCF failure or arbitrary exception is reclassified.
"""
import ast
import json
from pathlib import Path
import re

import numpy as np
from ase import Atoms
from ase.io import read

from kinbot.execution import file_sha256
from kinbot.optimization_budget import step_retries


def evidence(qc, job, job_id):
    from kinbot.timeout_retry import _matching
    log,worker=Path(job+'_sella.log'),Path(job+'.py')
    if not log.is_file() or not worker.is_file():
        return None
    with log.open('rb') as f:
        f.seek(max(0,log.stat().st_size-16384));lines=f.read().decode(errors='replace').splitlines()
    failures=[s for s in lines if s.startswith('Worker failure:')]
    if not failures:
        return None
    failure=failures[-1]
    exhausted=re.fullmatch(r'Worker failure: Optimization did not converge within (\d+) steps\.',failure)
    ode=failure=='Worker failure: Optimization failed: RuntimeError: Geometry update ODE is taking too long to converge!'
    if not exhausted and not ode:
        return None
    tree=ast.parse(worker.read_text())
    initial=None; solvers=set()
    for node in ast.walk(tree):
        if not isinstance(node,ast.Assign) or len(node.targets)!=1 or not isinstance(node.targets[0],ast.Name):
            continue
        if not isinstance(node.value,ast.Call) or not isinstance(node.value.func,ast.Name):
            continue
        if node.targets[0].id=='mol' and node.value.func.id=='Atoms':
            try:
                values={k.arg:ast.literal_eval(k.value) for k in node.value.keywords}
                initial=Atoms(symbols=values['symbols'],positions=values['positions'])
            except (ValueError,KeyError,TypeError):
                return None
        if node.targets[0].id=='opt' and node.value.func.id in ('Sella','BFGS'):
            if any(k.arg=='constraints' for k in node.value.keywords):
                return None
            solvers.add(node.value.func.id)
    if initial is None or not solvers or (ode and 'Sella' not in solvers):
        return None
    trajectory=Path(job+'.traj')
    if not trajectory.is_file():
        return None
    rows=[s.split() for s in lines if s.startswith(('Sella ','BFGS:'))]
    if not rows:
        return None
    row=rows[-1]
    if exhausted and int(row[1])<int(exhausted[1]):
        return None
    # Sella logs its evaluated trajectory index; later Hessian probes are not
    # necessarily the accepted point. Older logs without that column use last.
    index=int(row[-1]) if row[0]=='Sella' and len(row)==9 else -1
    first,last=read(trajectory,index=0),read(trajectory,index=index)
    if (first.get_chemical_symbols()!=initial.get_chemical_symbols()
            or last.get_chemical_symbols()!=initial.get_chemical_symbols()
            or not np.allclose(first.positions,initial.positions,atol=1e-6,rtol=0)
            or not np.isfinite(last.positions).all()
            or not np.isfinite(last.get_forces()).all()
            or not np.isfinite(last.get_potential_energy())
            or abs(last.get_potential_energy()-float(row[3]))>2e-6):
        return None
    budget_file=Path(job+'.optimization_budget.json')
    used=json.loads(budget_file.read_text()).get('chunks_used') if budget_file.exists() else 1
    if type(used) is not int or used<0:
        return None
    remaining=1+step_retries()-used
    if exhausted and remaining<=0:
        return None
    allocation=Path('.kinbot_timeout_allocations')/(job_id+'.json')
    if not allocation.is_file():
        return None
    native=json.loads(allocation.read_text())
    if not _matching(native,Path.cwd(),job_id) or not native.get('TimeLimit') or not native.get('Partition'):
        return None
    return dict(job_id=job_id,source='optimizer_ode' if ode else 'optimizer_steps',
        previous_walltime=native['TimeLimit'],partition=native['Partition'],
        worker_sha256=file_sha256(worker),trajectory_sha256=file_sha256(trajectory),
        positions=last.positions.tolist(),chunks_used=used,chunks=max(1,remaining),
        native_sha256=file_sha256(log))


def render(source, proof, job):
    """Change only the verified seed, step budget, and (ODE only) coordinates."""
    class Continuation(ast.NodeTransformer):
        def visit_Assign(self,node):
            node=self.generic_visit(node)
            if (len(node.targets)==1 and isinstance(node.targets[0],ast.Name)
                    and node.targets[0].id=='mol' and isinstance(node.value,ast.Call)
                    and isinstance(node.value.func,ast.Name) and node.value.func.id=='Atoms'):
                for k in node.value.keywords:
                    if k.arg=='positions':
                        k.value=ast.parse(repr(proof['positions']),mode='eval').body
            if (proof['source']=='optimizer_ode' and len(node.targets)==1
                    and isinstance(node.targets[0],ast.Name) and node.targets[0].id=='sella_kwargs'):
                return [node,*ast.parse("sella_kwargs['internal'] = False").body]
            return node
        def visit_Call(self,node):
            node=self.generic_visit(node)
            if (isinstance(node.func,ast.Attribute) and isinstance(node.func.value,ast.Name)
                    and node.func.value.id=='opt' and node.func.attr=='run'):
                node.func=ast.Name(id='run_optimizer',ctx=ast.Load());node.args.insert(0,ast.Name(id='opt',ctx=ast.Load()))
            if isinstance(node.func,ast.Name) and node.func.id=='run_optimizer':
                node.keywords=[k for k in node.keywords if k.arg not in ('chunks','used','report')]
                for k,v in dict(chunks=proof['chunks'],used=proof['chunks_used'],report=job+'.optimization_budget.json').items():
                    node.keywords.append(ast.keyword(arg=k,value=ast.Constant(value=v)))
            return node
    tree=Continuation().visit(ast.parse(source))
    tree.body.insert(0,ast.ImportFrom(module='kinbot.optimization_budget',names=[ast.alias(name='run_optimizer')],level=0))
    ast.fix_missing_locations(tree)
    return ast.unparse(tree)+'\n'
