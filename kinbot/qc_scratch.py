"""Compute-node scratch selection and bounded publication of native QC results.

Capacity is not locality. In particular, neither HOME nor the calculation
directory is a scratch fallback. All paths are checked on the executing host.
"""
from contextlib import contextmanager
import ctypes
import errno
import json
import os
from pathlib import Path
import pwd
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import threading


LOCAL_FILESYSTEMS = {'ext2', 'ext3', 'ext4', 'xfs', 'btrfs', 'zfs', 'ufs',
                     'apfs', 'hfs', 'hfsplus', 'overlay'}


def configure(parameters):
    """Carry single-input site preferences to all nested workers via Slurm."""
    if 'qc_scratch_fallbacks' in parameters:
        paths = parameters['qc_scratch_fallbacks']
        if (not isinstance(paths, list) or
                any(not isinstance(p, str) or not Path(p).expanduser().is_absolute()
                    or os.pathsep in p for p in paths)):
            raise ValueError('qc_scratch_fallbacks must be a list of absolute local paths')
        os.environ['KINBOT_QC_SCRATCH_FALLBACKS'] = os.pathsep.join(paths)
    if 'qc_scratch_min_mb' in parameters:
        value = parameters['qc_scratch_min_mb']
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError('qc_scratch_min_mb must be a positive integer')
        os.environ['KINBOT_QC_SCRATCH_MIN_MB'] = str(value)


def account_scratch():
    """Use the actual account, not an inherited USER/HOME from another job."""
    return Path('/scratch') / pwd.getpwuid(os.getuid()).pw_name


def filesystem(path):
    """Describe the containing mount; unknown/network/RAM mounts fail closed."""
    path = Path(path).resolve()
    if sys.platform == 'linux':
        def unescape(value):
            return re.sub(r'\\([0-7]{3})', lambda m: chr(int(m[1], 8)), value)
        mounts = []
        for line in Path('/proc/self/mountinfo').read_text().splitlines():
            left, right = line.split(' - ', 1)
            fields, detail = left.split(), right.split()
            mount = Path(unescape(fields[4]))
            if path == mount or mount in path.parents:
                mounts.append((len(mount.parts), str(mount), detail[0], unescape(detail[1])))
        if not mounts:
            raise OSError('Cannot identify scratch filesystem')
        _, mount, kind, source = max(mounts)
    elif sys.platform == 'darwin':
        existing = path
        while not existing.exists():
            existing = existing.parent
        class StatFS(ctypes.Structure):
            _fields_ = [('bsize', ctypes.c_uint32), ('iosize', ctypes.c_int32),
                        ('counts', ctypes.c_uint64 * 5), ('fsid', ctypes.c_int32 * 2),
                        ('owner', ctypes.c_uint32), ('type', ctypes.c_uint32),
                        ('flags', ctypes.c_uint32), ('subtype', ctypes.c_uint32),
                        ('fstype', ctypes.c_char * 16), ('mount', ctypes.c_char * 1024),
                        ('source', ctypes.c_char * 1024), ('reserved', ctypes.c_uint32 * 8)]
        info = StatFS()
        library = ctypes.CDLL(None, use_errno=True)
        if library.statfs(os.fsencode(existing), ctypes.byref(info)):
            error = ctypes.get_errno()
            raise OSError(error, os.strerror(error))
        kind, mount, source = (value.decode() for value in (info.fstype, info.mount, info.source))
    else:
        raise OSError('Scratch filesystem verification is unsupported on this platform')
    if kind not in LOCAL_FILESYSTEMS:
        raise OSError(f'Scratch must be disk-backed local storage; found {kind} at {mount}')
    return {'mount': mount, 'filesystem': kind, 'device': source}


def candidates(env, backend):
    preferred = account_scratch()
    # Never create /scratch on the root disk when the site has no such mount.
    if preferred.parent.is_dir():
        yield 'account_scratch', preferred
    keys = ['KINBOT_QC_SCRATCH', f'KINBOT_{backend.upper()}_SCRATCH']
    for key in dict.fromkeys(keys):
        if env.get(key):
            yield key, Path(env[key]).expanduser()
    for key in ('KINBOT_QC_SCRATCH_FALLBACKS', f'KINBOT_{backend.upper()}_SCRATCH_FALLBACKS'):
        if env.get(key):
            for value in env[key].split(os.pathsep):
                if not value.strip():
                    raise ValueError(f'{key} contains an empty path')
                yield key, Path(value).expanduser()
    keys = ['SLURM_TMPDIR', 'SCRATCH', 'TMPDIR']
    if backend in ('gaussian', 'gauss'):
        keys.append('GAUSS_SCRDIR')
    for key in keys:
        if env.get(key):
            yield key, Path(env[key]).expanduser()
    yield 'system_tmp', Path('/tmp')


def allocate(env, backend, minimum_mb, *, probe, work_directory=None, avoid_roots=()):
    nodes = env.get('SLURM_NNODES', env.get('SLURM_JOB_NUM_NODES', '1'))
    if int(nodes) != 1:
        raise ValueError('Node-local QC workspaces require a single-node allocation')
    if isinstance(minimum_mb, bool) or not isinstance(minimum_mb, int) or minimum_mb < 1:
        raise ValueError('QC scratch minimum must be a positive integer')
    avoided = list(avoid_roots)
    for key in ('KINBOT_QC_AVOID_SCRATCH_ROOTS', f'KINBOT_{backend.upper()}_AVOID_SCRATCH_ROOTS'):
        inherited = json.loads(env.get(key, '[]'))
        if not isinstance(inherited, list) or any(not isinstance(p, str) or not p for p in inherited):
            raise ValueError(f'{key} must contain a JSON list of paths')
        avoided.extend(inherited)
    for key in ('KINBOT_QC_SCRATCH_MIN_MB', f'KINBOT_{backend.upper()}_SCRATCH_MIN_MB'):
        if env.get(key):
            if not re.fullmatch(r'[1-9][0-9]*', env[key]):
                raise ValueError(f'{key} must be a positive integer')
            minimum_mb = max(minimum_mb, int(env[key]))
    roots = list(candidates(env, backend))  # Validate all configured alternatives first.
    homes = {Path(pwd.getpwuid(os.getuid()).pw_dir), Path(env.get('HOME', Path.home()))}
    for home in list(homes):
        try:
            homes.add(home.resolve())
        except OSError:
            pass  # Broken HOME must not prevent choosing working local storage.
    rejected, seen, errors = [], set(), []
    for source, candidate in roots:
        scratch = None
        try:
            # Reject obvious home paths before touching a possibly broken mount.
            if any(candidate == h or h in candidate.parents for h in homes):
                raise OSError('Home directories cannot be used for QC scratch')
            root = candidate.resolve()
            if any(root == Path(previous).expanduser().resolve() for previous in avoided):
                raise OSError('Scratch root had an explicit I/O failure in the preceding attempt')
            if root in seen:
                continue
            seen.add(root)
            if any(root == h or h in root.parents for h in homes):
                raise OSError('Scratch resolves into a home directory')
            if work_directory and root == Path(work_directory).resolve():
                raise OSError('The result directory cannot be used for QC scratch')
            mount = filesystem(root)
            root.mkdir(parents=True, exist_ok=True)
            free = shutil.disk_usage(root).free // 1024**2
            if free < minimum_mb:
                raise OSError(f'requires {minimum_mb} MB; {free} MB free')
            scratch = Path(tempfile.mkdtemp(prefix=f'kinbot-{backend}-', dir=root))
            probe(scratch)
            return scratch, {
                'scratch_directory': str(scratch), 'scratch_root': str(root),
                'scratch_source': source, 'scratch_host': socket.gethostname(),
                'scratch_available_mb': free, 'scratch_minimum_mb': minimum_mb,
                'scratch_io_checked': True, 'scratch_mount': mount,
                'scratch_rejected_candidates': rejected,
            }
        except OSError as error:
            if scratch is not None:
                shutil.rmtree(scratch, ignore_errors=True)
            rejected.append(f'{source}={candidate}: {error}')
            errors.append(error.errno)
    code = next((e for e in errors if e in (errno.EIO, errno.ESTALE, errno.ETIMEDOUT)), errno.ENOSPC)
    raise OSError(code, 'No verified local QC scratch is usable: ' + '; '.join(rejected))


def _basename(name):
    if not name or Path(name).name != name or name in ('.', '..'):
        raise ValueError(f'Native QC file must be a relative basename: {name!r}')
    return name


def references(text, backend):
    """Only explicitly referenced restart/input files are staged, never all WFUs."""
    names = set()
    if backend == 'molpro':
        for name in re.findall(r'(?im)(?:^|;)\s*file\s*,\s*\d+\s*,\s*([^,;\s]+)', text):
            names.add(_basename(name.strip('\"\'')))
        for name in re.findall(r'(?im)(?:^|;)\s*(?:include\s*,|geometry\s*=)\s*([\w./-]+)', text):
            names.add(_basename(name))
        for name in re.findall(r'(?im)(?:^|;)\s*put\s*,\s*\w+\s*,\s*([^,;\s]+)', text):
            _basename(name.strip('\"\''))
    elif backend in ('gaussian', 'gauss'):
        for name in re.findall(r'(?im)^\s*%(?:oldchk|chk)\s*=\s*(\S+)', text):
            name = _basename(name)
            names.add(name if Path(name).suffix else name + '.chk')
        for name in re.findall(r'(?im)^\s*%(?:rwf|int|d2e|skr)\s*=\s*(\S+)', text):
            for field in name.split(','):
                if field and not re.fullmatch(r'[-+]?\d+(?:[KMGTP]B|[KMGTP]W)?', field, re.I):
                    _basename(field)
    return names


def molpro_command(command, runtime, work):
    """Override site/default repositories, including inherited -W and -I."""
    clean, skip = [], False
    options = ('-d', '-W', '-I', '--directory', '--wavefunction-file-repository',
               '--main-file-repository')
    for arg in command:
        if skip:
            skip = False
            continue
        if arg in options:
            skip = True
        elif any(arg.startswith(key + '=') for key in options if key.startswith('--')):
            continue
        elif any(arg.startswith(key) and len(arg) > 2 for key in ('-d', '-W', '-I')):
            continue
        else:
            clean.append(arg)
    return [clean[0], '-d', runtime['scratch_directory'], '-W', str(work),
            '-I', runtime['scratch_directory'], *clean[1:]]


@contextmanager
def workspace(directory, runtime, *, backend, inputs, outputs=(), label='native'):
    """Execute in local storage; publish results/checkpoints, not integral files.

    A small live tail remains visible to the existing monitor. If publication
    fails, local evidence is retained and no successful completion is returned.
    """
    directory = Path(directory).resolve()
    work = Path(runtime['scratch_directory']) / 'work'
    work.mkdir()
    runtime['native_work_directory'] = str(work)
    staged = set(inputs)
    for name in inputs:
        name = _basename(name)
        path = directory / name
        if path.suffix.lower() in ('.inp', '.com'):
            staged.update(references(path.read_text(), backend))
    for name in staged:
        source = directory / _basename(name)
        if source.is_file():
            shutil.copy2(source, work / name)
        elif name in inputs:
            raise FileNotFoundError(source)
    runtime_file = directory / f'{_basename(label)}.runtime.json'
    from kinbot.execution import atomic_json
    atomic_json(runtime_file, runtime)
    stop = threading.Event()
    live = directory / f'{label}.live.log'
    def mirror():
        while not stop.wait(30):
            try:
                paths = sorted((p for p in work.iterdir() if p.suffix in ('.out', '.log', '.err', '.stderr')),
                               key=lambda p: p.stat().st_mtime)[-3:]
                chunks = []
                for path in paths:
                    with path.open('rb') as stream:
                        stream.seek(max(0, path.stat().st_size - 16384))
                        chunks.append(path.name + '\n' + stream.read().decode(errors='replace'))
                if chunks:
                    live.write_text('\n'.join(chunks))
            except OSError:
                return  # Do not repeatedly hammer an unavailable result filesystem.
    thread = threading.Thread(target=mirror, daemon=True)
    thread.start()
    try:
        yield work
    finally:
        stop.set()
        thread.join(timeout=1)
        try:
            explicit = set(outputs) | staged
            suffixes = {'.out', '.log', '.err', '.stdout', '.stderr', '.xml', '.xyz',
                        '.molden', '.wfu', '.chk', '.fchk', '.grad'}
            for path in work.iterdir():
                if path.is_symlink():
                    raise OSError(f'Refusing to publish a native symlink: {path}')
                if path.is_file() and (path.name in explicit or path.suffix in suffixes
                                       or path.name.endswith('.wfu.failed')):
                    target = directory / path.name
                    temporary = directory / f'.{path.name}.publish.{os.getpid()}'
                    shutil.copy2(path, temporary)
                    os.replace(temporary, target)
            if live.exists():
                live.unlink()
        except OSError as error:
            runtime['retain_scratch'] = True
            raise RuntimeError(f'QC result publication failed; local files retained on '
                               f'{socket.gethostname()} at {work}: {error}') from error


def main(argv=None):
    """Native shell entry point; preparation on a login node selects no scratch."""
    import argparse
    from kinbot.anl.runtime import qc_runtime_environment, cleanup_qc_runtime
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--backend', required=True)
    parser.add_argument('--input')
    parser.add_argument('--probe', action='store_true', help='Check local scratch on this host; launch no QC')
    parser.add_argument('--min-mb', type=int, default=4096)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if args.probe:
        import json
        from kinbot.anl.runtime import _check_scratch_io
        _, runtime = allocate(dict(os.environ), args.backend, args.min_mb, probe=_check_scratch_io)
        print(json.dumps(runtime, indent=2))
        cleanup_qc_runtime(runtime)
        return 0
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or not args.input:
        parser.error('--input and a native command after -- are required unless using --probe')
    path = Path(args.input).resolve(strict=True)
    env, runtime = qc_runtime_environment(command[0], args.backend, work_directory=path.parent)
    try:
        with workspace(path.parent, runtime, backend=args.backend,
                       inputs=[path.name], label=path.stem) as local:
            command = [*command, path.name]
            if args.backend == 'molpro':
                command = molpro_command(command, runtime, local)
            with (local / f'{path.stem}.stdout').open('wb') as out, \
                    (local / f'{path.stem}.stderr').open('wb') as err:
                result = subprocess.run(command, cwd=local, env=env, stdout=out, stderr=err)
        return result.returncode
    finally:
        cleanup_qc_runtime(runtime)


if __name__ == '__main__':
    sys.exit(main())
