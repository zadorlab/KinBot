"""Bounded, immutable publication of native restart candidates from local scratch.

Byte integrity and stable source files do not prove native restart applicability.
Backends must still validate method, geometry, state, records and restart phase.
"""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import stat
import time
import uuid

from kinbot.execution import file_sha256
from kinbot.anl.optimization_progress import publish

DEFAULTS = {'qc_checkpoint_interval': 0, 'qc_checkpoint_max_mb': 16384,
            'qc_checkpoint_generations': 2}
NAMED = {'fort.16', 'fort.17', 'JOBARC', 'JAINDX', 'OPTARC', 'FCMINT',
         'FCMFINAL', 'GRD', 'HESSIAN'}
SUFFIXES = {'.wfu', '.chk', '.rwf', '.gbw'}


def validate(values):
    for name,default in DEFAULTS.items():
        value=values.get(name,default)
        low,high = ((0,86400) if name=='qc_checkpoint_interval' else
                    (1,1048576) if name=='qc_checkpoint_max_mb' else (2,10))
        if type(value) is not int or not low <= value <= high:
            raise ValueError(f'{name} must be an integer from {low} to {high}.')


def configure(values):
    validate(values)
    for key in DEFAULTS:
        if key in values:os.environ['KINBOT_'+key.upper()]=str(values[key])


def settings():
    result={k:int(os.environ.get('KINBOT_'+k.upper(),v)) for k,v in DEFAULTS.items()}
    validate(result)
    return result


def signature(path):
    s=path.stat(follow_symlinks=False)
    if not stat.S_ISREG(s.st_mode):
        raise OSError(f'Native checkpoint is not a regular file: {path}')
    return (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns)


def _copy(source, target, expected, stop):
    """Never follow a replaced symlink or chase a growing scratch file."""
    flags=os.O_RDONLY | getattr(os,'O_NOFOLLOW',0) | getattr(os,'O_NONBLOCK',0)
    fd=os.open(source,flags)
    sha=hashlib.sha256()
    with os.fdopen(fd,'rb') as incoming, target.open('xb') as outgoing:
        s=os.fstat(incoming.fileno())
        actual=(s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns)
        if not stat.S_ISREG(s.st_mode) or actual!=expected:
            raise OSError(f'Native checkpoint replaced before copy: {source}')
        remaining=expected[2]
        while remaining:
            if stop is not None and stop.is_set():
                raise InterruptedError('Checkpoint copy stopped before publication.')
            block=incoming.read(min(1024**2,remaining))
            if not block:
                raise OSError(f'Native checkpoint shrank during copy: {source}')
            remaining-=len(block);sha.update(block);outgoing.write(block)
        if incoming.read(1):
            raise OSError(f'Native checkpoint grew during copy: {source}')
        outgoing.flush();os.fsync(outgoing.fileno())
    if signature(source)!=expected or file_sha256(target)!=sha.hexdigest():
        raise OSError(f'Native checkpoint changed during copy: {source}')
    return sha.hexdigest()


def _flush_tree(root):
    # Flush directory entries too: a durable pointer must not reference files
    # whose names were only present in the filesystem's volatile cache.
    for parent,_,files in os.walk(root,topdown=False):
        for name in files:
            with (Path(parent)/name).open('rb') as stream:os.fsync(stream.fileno())
        fd=os.open(parent,os.O_RDONLY|getattr(os,'O_DIRECTORY',0))
        try:os.fsync(fd)
        finally:os.close(fd)


class Publisher:
    def __init__(self,directory,work,runtime,backend,label,inputs,*,config=None):
        self.config=settings() if config is None else {**DEFAULTS,**config}
        validate(self.config)
        if not label or Path(label).name!=label or label in ('.','..'):
            raise ValueError('Native checkpoint label must be a basename.')
        if any(Path(name).name!=name or name in ('.','..') for name in inputs):
            raise ValueError('Native checkpoint inputs must be basenames.')
        self.interval=self.config['qc_checkpoint_interval']
        self.directory=Path(directory);self.work=Path(work)
        self.label=label
        self.scratch=Path(runtime['scratch_directory']).resolve()
        self.runtime=runtime
        self.attempt=uuid.uuid4().hex
        self.next_at=time.monotonic()+self.interval
        self.failures=0
        self.root=self.directory/'.kinbot_checkpoints'/'native'/label/self.attempt
        self.contract={}
        if self.interval:
            self.contract=dict(backend=backend,label=label,
                inputs={name:file_sha256(self.work/name) for name in sorted(inputs) if (self.work/name).is_file()},
                task_sha256=file_sha256(self.directory/'task.json') if (self.directory/'task.json').is_file() else None,
                slurm_job_id=os.environ.get('SLURM_JOB_ID'),
                slurm_restart_count=os.environ.get('SLURM_RESTART_COUNT','0'),
                scratch_directory=str(self.scratch),host=socket.gethostname())

    def _candidates(self):
        paths=[]
        for parent,dirs,files in os.walk(self.scratch,followlinks=False):
            parent=Path(parent)
            dirs[:]=[name for name in dirs if not (parent/name).is_symlink()
                     and len(parent.relative_to(self.scratch).parts)<2]
            for name in files:
                p=parent/name
                if p.name in NAMED or p.suffix.lower() in SUFFIXES:
                    if p.is_symlink() or p.resolve()!=p:
                        raise OSError(f'Native checkpoint is a symlink: {p}')
                    if signature(p)[2]:paths.append(p)
        return sorted(paths)

    def capture(self,stop=None):
        paths=self._candidates()
        if not paths:return None
        initial={p:signature(p) for p in paths}
        total=sum(s[2] for s in initial.values())
        if total>self.config['qc_checkpoint_max_mb']*1024**2:
            raise OSError(f'Native checkpoint set ({total} bytes) exceeds configured copy ceiling.')
        previous=self.root/'latest.json'
        self._prune()
        if previous.exists():
            prior=json.loads(previous.read_text())
            # Unchanged native bytes need no new shared-filesystem copy. The
            # signature check is an optimization, not a restoration proof.
            if prior.get('source_signatures')=={str(p):list(s) for p,s in initial.items()}:
                return prior
        token='generation-'+uuid.uuid4().hex
        temporary=self.root/('.'+token+'.preparing');temporary.mkdir(parents=True)
        entries=[]
        try:
            for source in paths:
                target=temporary/'files'/source.relative_to(self.scratch)
                target.parent.mkdir(parents=True,exist_ok=True)
                digest=_copy(source,target,initial[source],stop)
                entries.append(dict(source=str(source),path=str(target.relative_to(temporary)),
                                    bytes=initial[source][2],sha256=digest))
            # Keep the exact input text and bounded diagnostics with the binary
            # files. Output tails are observations, not completion/phase proofs.
            diagnostic=temporary/'diagnostics';diagnostic.mkdir()
            for name in self.contract['inputs']:
                source=self.work/name
                if source.suffix.lower() not in SUFFIXES and source.name not in NAMED:
                    if file_sha256(source)!=self.contract['inputs'][name]:
                        raise OSError(f'Native input changed: {name}')
                    shutil.copy2(source,diagnostic/name)
                    if file_sha256(diagnostic/name)!=self.contract['inputs'][name]:
                        raise OSError(f'Native input changed during copy: {name}')
            for p in self.work.iterdir():
                if p.is_file() and not p.is_symlink() and p.suffix in ('.out','.log','.err','.stdout','.stderr'):
                    with p.open('rb') as stream:
                        stream.seek(max(0,p.stat().st_size-1024**2))
                        (diagnostic/(p.name+'.tail')).write_bytes(stream.read(1024**2))
            if self._candidates()!=paths or any(signature(p)!=s for p,s in initial.items()):
                raise OSError('Native checkpoint set changed before generation publication.')
            if stop is not None and stop.is_set():
                raise InterruptedError('Checkpoint copy stopped before publication.')
            manifest=dict(schema=1,contract=self.contract,files=entries,bytes=total,
                          diagnostics={p.name:file_sha256(p) for p in diagnostic.iterdir()},
                          created_at=datetime.now(timezone.utc).isoformat(),
                          status='restart_candidate_requires_native_validation')
            _flush_tree(temporary)
            publish(temporary/'manifest.json',manifest)
            destination=self.root/token;temporary.rename(destination)
            pointer=dict(generation=token,sha256=file_sha256(destination/'manifest.json'),
                         source_signatures={str(p):list(s) for p,s in initial.items()})
            publish(previous,pointer)
            self._prune()
            self.runtime['native_checkpoint']=dict(directory=str(destination),**manifest)
            return pointer
        finally:
            if temporary.exists():shutil.rmtree(temporary)

    def _prune(self):
        """Bound even orphan generations left by interrupted pointer publication."""
        pointer=self.root/'latest.json'
        current=json.loads(pointer.read_text()).get('generation') if pointer.exists() else None
        generations=sorted(self.root.glob('generation-*'),key=lambda p:p.stat().st_mtime_ns,reverse=True)
        keep={p for p in generations if p.name==current}
        for p in generations:
            if len(keep)<self.config['qc_checkpoint_generations']:keep.add(p)
        for p in generations:
            if p not in keep and p.is_dir() and not p.is_symlink():shutil.rmtree(p)

    def poll(self,*,force=False,stop=None):
        if not self.interval or (not force and time.monotonic()<self.next_at):return
        try:
            result=self.capture(stop)
        except (OSError,ValueError,KeyError) as error:
            self.failures+=1
            self.runtime['native_checkpoint_error']=str(error)
            self.next_at=time.monotonic()+min(3600,max(30,self.interval)*2**min(self.failures,6))
            # Keep retrying after storage recovers, retaining prior generations.
            self._status()
            return None
        self.failures=0;self.runtime.pop('native_checkpoint_error',None)
        self.next_at=time.monotonic()+self.interval
        self._status()
        return result

    def _status(self):
        try:
            publish(self.directory/(self.label+'.checkpoint.json'),dict(schema=1,
                attempt=self.attempt,checkpoint=self.runtime.get('native_checkpoint'),
                error=self.runtime.get('native_checkpoint_error'),consecutive_failures=self.failures,
                updated_at=datetime.now(timezone.utc).isoformat()))
        except OSError:
            pass  # A disconnected shared filesystem cannot receive its own diagnostic.
