"""Immutable accepted geometry points for interrupted ASE optimizations.

This is geometry continuation, not a native force/frequency or Hessian restart.
The next calculator evaluation and all final scientific checks still run.
"""
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import uuid

import numpy as np
from ase import Atoms
from ase.io import read

from kinbot.execution import atomic_json, file_sha256

NAME = '.kinbot_checkpoints'


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def publish(path, data):
    """Flush complete checkpoint bytes before replacing the visible pointer."""
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=None
    try:
        with tempfile.NamedTemporaryFile(mode='w',dir=path.parent,delete=False) as stream:
            temporary=Path(stream.name)
            stream.write(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')
            stream.flush();os.fsync(stream.fileno())
        os.replace(temporary,path)
        fd=os.open(path.parent,os.O_RDONLY|getattr(os,'O_DIRECTORY',0))
        try:os.fsync(fd)
        finally:os.close(fd)
    finally:
        if temporary is not None and temporary.exists():temporary.unlink()


def restore_attempt_start(directory, record, atoms):
    """Verify the separately bound restart point when checking a trajectory."""
    path=directory/'optimization_resume.json'
    if not path.exists():return
    receipt=json.loads(path.read_text());seed=receipt.get('checkpoint')
    if not seed:return
    expected=dict(record_sha256=digest(record),starting_symbols=atoms.get_chemical_symbols(),
                  starting_positions=atoms.positions.tolist())
    source=(directory/seed['generation']).resolve(strict=True)
    if (receipt.get('contract')!=expected or not source.is_relative_to((directory/NAME).resolve())
            or file_sha256(source)!=seed['sha256']):
        raise RuntimeError('Optimizer attempt starting checkpoint changed.')
    data=json.loads(source.read_text())
    if data.get('contract')!=expected or data['positions']!=seed['positions']:
        raise RuntimeError('Optimizer attempt start differs from checkpoint.')
    atoms.positions=data['positions']


def stage_previous(directory, previous):
    """Copy a bound accepted point into an identical restaged task.

    The complete earlier attempt stays in its archive. Resource/scientific/task
    changes do not inherit this checkpoint through the identical-task path.
    """
    previous = Path(previous)
    old = previous/'task.json'
    current = directory/'task.json'
    if not old.is_file() or digest(json.loads(old.read_text())) != digest(json.loads(current.read_text())):
        return
    for pointer in (previous/NAME/'optimization').glob('*/latest.json'):
        if pointer.resolve() != pointer or not re.fullmatch(r'[a-f0-9]{64}',pointer.parent.name):
            raise RuntimeError('Archived optimizer checkpoint path is invalid.')
        receipt = json.loads(pointer.read_text())
        name = receipt.get('generation', '')
        if not re.fullmatch(r'point-[a-f0-9]{32}\.json',name):
            raise RuntimeError('Archived optimizer checkpoint generation is invalid.')
        source = pointer.parent/name
        data = source.read_bytes()
        if source.is_symlink() or hashlib.sha256(data).hexdigest() != receipt.get('sha256'):
            raise RuntimeError('Archived optimizer checkpoint changed before restaging.')
        target = directory/NAME/'optimization'/pointer.parent.name
        target.mkdir(parents=True,exist_ok=True)
        publish(target/name,json.loads(data))
        # Preserve the hash of the newly serialized, identical point.
        publish(target/'latest.json',dict(generation=name,sha256=file_sha256(target/name)))
        publish(target/'import.json',dict(source=str(source),source_sha256=receipt['sha256'],
                                         imported_sha256=file_sha256(target/name)))


@contextmanager
def ownership(directory):
    root=Path(directory)/NAME;root.mkdir(exist_ok=True)
    with (root/'optimization.lock').open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as e:
            raise RuntimeError('Another process owns this optimization; no files were changed.') from e
        yield


class Progress:
    def __init__(self, directory, record, atoms, *, used_chunks, chunk_steps, max_chunks):
        self.directory=Path(directory)
        self.contract=dict(record_sha256=digest(record),
                           starting_symbols=atoms.get_chemical_symbols(),
                           starting_positions=atoms.positions.tolist())
        self.root=self.directory/NAME/'optimization'/digest(self.contract)
        self.root.mkdir(parents=True,exist_ok=True)
        self.chunk_steps=chunk_steps
        self.limit=max_chunks*chunk_steps
        self.steps=used_chunks*chunk_steps
        self.base=self.steps
        self.loaded=None
        pointer=self.root/'latest.json'
        if pointer.exists():
            receipt=json.loads(pointer.read_text());name=receipt.get('generation','')
            if not re.fullmatch(r'point-[a-f0-9]{32}\.json',name):
                raise RuntimeError('Invalid optimizer checkpoint generation.')
            path=self.root/name
            if path.is_symlink() or file_sha256(path)!=receipt.get('sha256'):
                raise RuntimeError('Optimizer checkpoint changed; prior work retained.')
            data=json.loads(path.read_text())
            if (data.get('schema')!=1 or data.get('contract')!=self.contract
                    or type(data.get('steps_used')) is not int
                    or not self.steps <= data['steps_used'] <= self.limit
                    or data.get('step_limit')!=self.limit):
                raise RuntimeError('Optimizer checkpoint is incompatible with the task or step budget.')
            point=Atoms(data['symbols'],positions=data['positions'])
            forces=np.asarray(data['forces_ev_angstrom'])
            if (point.get_chemical_symbols()!=atoms.get_chemical_symbols()
                    or forces.shape!=point.positions.shape
                    or not np.isfinite(point.positions).all() or not np.isfinite(forces).all()
                    or not np.isfinite(data['energy_ev'])):
                raise RuntimeError('Invalid evaluated optimizer checkpoint.')
            atoms.positions=point.positions
            self.loaded=dict(generation=str(path.relative_to(self.directory)),sha256=receipt['sha256'],
                             positions=data['positions'],steps_used=data['steps_used'])
            self.steps=self.base=data['steps_used']

    def archive_previous_attempt(self, task_id):
        """Retain current partial optimizer/force files before counters restart."""
        paths=[]
        for p in self.directory.iterdir():
            if p.is_symlink():continue
            if p.is_file() and (p.name in ('optimization.traj','optimization.log','optimization_budget.json',
                                          'optimization_resume.json')
                    or p.name.startswith(task_id+'_step_')):
                paths.append(p)
        if not paths:return None
        archive=self.root/'attempts'/uuid.uuid4().hex;archive.mkdir(parents=True)
        # Short hashed WFU names may not have the step prefix. Copy them before
        # later evaluations can reuse their names, retaining staged seeds in place.
        from kinbot.qc_scratch import references
        extras=set()
        for p in paths:
            if p.suffix=='.inp':
                extras.update(self.directory/n for n in references(p.read_text(),'molpro'))
        files={}
        for p in sorted(extras-set(paths)):
            if p.is_file() and not p.is_symlink():
                shutil.copy2(p,archive/p.name);files[p.name]=file_sha256(archive/p.name)
        for p in paths:
            p.rename(archive/p.name);files[p.name]=file_sha256(archive/p.name)
        atomic_json(archive/'manifest.json',dict(schema=1,files=files,checkpoint=self.loaded))
        return str(archive.relative_to(self.directory))

    def attach(self, optimizer, atoms):
        # Called once per optimizer instance. A minimum-recovery refinement can
        # construct another optimizer, but cannot reset the global step count.
        self.base=self.steps
        def save():
            point=atoms
            if optimizer.__class__.__name__=='Sella':
                # Sella's trajectory also includes Hessian probes. Its log
                # explicitly identifies the current accepted trajectory frame.
                rows=[s.split() for s in (self.directory/'optimization.log').read_text().splitlines()
                      if s.startswith('Sella ')]
                if not rows or len(rows[-1])!=9 or int(rows[-1][1])!=optimizer.nsteps:
                    raise RuntimeError('Cannot identify the accepted Sella checkpoint frame.')
                row=rows[-1];point=read(self.directory/'optimization.traj',index=int(row[-1]))
                if abs(point.get_potential_energy()-float(row[3]))>2.e-6:
                    raise RuntimeError('Sella accepted checkpoint energy differs from its log.')
            energy=float(point.get_potential_energy());forces=point.get_forces(apply_constraint=False)
            if (point.get_chemical_symbols()!=self.contract['starting_symbols']
                    or not np.isfinite(point.positions).all() or not np.isfinite(energy)
                    or not np.isfinite(forces).all()):
                raise RuntimeError('Optimizer did not provide a finite evaluated checkpoint.')
            steps=self.base+optimizer.nsteps
            if not self.steps <= steps <= self.limit:
                raise RuntimeError('Optimizer checkpoint exceeded its original step budget.')
            data=dict(schema=1,contract=self.contract,symbols=point.get_chemical_symbols(),
                      positions=point.positions.tolist(),energy_ev=energy,forces_ev_angstrom=forces.tolist(),
                      steps_used=steps,step_limit=self.limit,optimizer=optimizer.__class__.__name__,
                      slurm_job_id=os.environ.get('SLURM_JOB_ID'),
                      slurm_restart_count=os.environ.get('SLURM_RESTART_COUNT','0'),
                      created_utc=datetime.now(timezone.utc).isoformat(),
                      restart_scope='geometry_only_recompute_force_and_optimizer_hessian')
            path=self.root/('point-'+uuid.uuid4().hex+'.json')
            publish(path,data)
            publish(self.root/'latest.json',dict(generation=path.name,sha256=file_sha256(path)))
            self.steps=steps
        optimizer.attach(save,interval=1)

    def remaining(self):
        return max(0,self.limit-self.steps)

    def receipt(self, archive):
        atomic_json(self.directory/'optimization_resume.json',dict(schema=1,contract=self.contract,
                    checkpoint=self.loaded,archive=archive,steps_used=self.steps,step_limit=self.limit,
                    restart_scope='geometry_only_recompute_force_and_optimizer_hessian'))
