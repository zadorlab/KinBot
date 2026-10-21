"""Journaled, opt-in home-cache cleanup between stopped and fresh campaigns.

Never recursively clean home, package caches, or node-local scratch from a
login node. A scheduler error is not evidence of an empty queue. Reviewed
orbital seeds must already have been copied/verified by the caller.
"""
import errno
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import socket
import stat
import subprocess
import uuid

from kinbot.execution import atomic_json
from kinbot.qc_scratch import account_identity


def validate(options):
    allowed = {'empty_home_wfu', 'remove_verified_home_scratch'}
    if not isinstance(options, dict) or set(options) - allowed:
        raise ValueError('Unknown cache cleanup options.')
    if any(type(v) is not bool for v in options.values()):
        raise ValueError('Cache cleanup options must be booleans.')


def _identity(path):
    st = path.lstat()
    if not stat.S_ISDIR(st.st_mode) or st.st_uid != os.getuid():
        raise RuntimeError(f'Cleanup requires an owned, non-symlink directory: {path}')
    return {'device': st.st_dev, 'inode': st.st_ino, 'mode': stat.S_IMODE(st.st_mode)}


def _overlap(a, b):
    return a == b or a in b.parents or b in a.parents


def _mounts():
    # Includes bind mounts, which can have the same st_dev as their parent.
    text = Path('/proc/self/mountinfo').read_text()
    return [Path(re.sub(r'\\([0-7]{3})', lambda m: chr(int(m[1], 8)), s.split()[4]))
            for s in text.splitlines()]


def _check_tree(path, protected):
    _identity(path)
    if any(_overlap(path, p.resolve()) for p in protected):
        raise RuntimeError(f'Cleanup target overlaps protected run/source data: {path}')
    if any(m == path or path in m.parents for m in _mounts()):
        raise RuntimeError(f'Cleanup target contains a mount point: {path}')


def assert_idle(targets):
    """Check all account allocations, and local QC/controllers/open cache files.

    This is not a reservation against a user independently starting new work.
    Call again after moving a target, before removing its retired contents.
    """
    account = account_identity()['name']
    q = subprocess.run(['squeue', '-h', '-u', account, '-o', '%i|%T'],
                       capture_output=True, text=True, check=True, timeout=20)
    if q.stdout.strip():
        raise RuntimeError('Cache cleanup waits for ALL account allocations, including pending '
                           'and completing jobs. Unrelated jobs were not cancelled. '
                           + ' | '.join(q.stdout.splitlines()[:8]))
    proc = Path('/proc')
    if not proc.is_dir():
        raise RuntimeError('Cache cleanup must run on the Linux controller host.')
    targets = [p.resolve() for p in targets]
    for p in proc.iterdir():
        if not p.name.isdigit() or int(p.name) == os.getpid():
            continue
        try:
            if p.stat().st_uid != os.getuid():
                continue
            args = (p/'cmdline').read_bytes().decode(errors='replace').split('\0')
            command = ' '.join(args)
            native = any(re.fullmatch(r'(?:molpro(?:\.exe)?|g16|l\d+\.exe|dmrcc|mrcc|xcfour)',
                                     Path(a).name) for a in args if a)
            controller = any(a in ('kinbot.pes', 'kinbot.kb', 'kinbot.anl.dispatch',
                                    'kinbot.run_refresh', 'kinbot.run_recover') for a in args)
            if native or controller:
                raise RuntimeError(f'Cache cleanup waits for local QC/controller PID {p.name}: {command[:180]}')
            for f in [p/'cwd', *(p/'fd').iterdir()]:
                try:
                    target = Path(os.readlink(f).removesuffix(' (deleted)'))
                except FileNotFoundError:
                    continue  # A process may close a descriptor during inspection.
                if target.is_absolute() and any(target == t or t in target.parents for t in targets):
                    raise RuntimeError(f'Cache cleanup waits for PID {p.name}, which uses {target}')
        except OSError as error:
            if error.errno not in (errno.ENOENT, errno.ESRCH):
                raise RuntimeError(f'Cannot establish cache inactivity for PID {p.name}: {error}') from error


def _scratch_candidates(archive, home, skipped):
    """Only direct children of the old home scratch root named in r6 evidence."""
    cache = home/'.kinbot_scratch'
    if not cache.exists() and not cache.is_symlink():
        return []
    _identity(cache)
    found = {}
    host = socket.gethostname().split('.')[0]
    count = 0
    for parent, dirs, files in os.walk(archive, followlinks=False):
        dirs[:] = [n for n in dirs if not n.startswith('.') and not Path(parent, n).is_symlink()]
        count += len(dirs) + len(files)
        if count > 500000:
            raise RuntimeError('Scratch inventory limit reached; no cleanup started.')
        for name in files:
            if not name.endswith(('.runtime.json', '.failure.json')):
                continue
            p = Path(parent, name)
            if p.is_symlink() or p.stat().st_size > 4*1024*1024:
                continue
            m = json.loads(p.read_text())
            runtime = m.get('runtime', m)
            raw = runtime.get('scratch_directory')
            if not raw:
                continue
            target = Path(raw)
            if target.parent != cache or not re.fullmatch(r'kinbot-[A-Za-z0-9_.-]+', target.name):
                continue
            # Per-node paths must not be interpreted on the login node. Even
            # historical home paths with ambiguous host provenance are retained.
            if runtime.get('scratch_host', m.get('host', '')).split('.')[0] != host:
                skipped[str(target)] = 'different or unknown execution host'
                continue
            if runtime.get('retain_scratch'):
                skipped[str(target)] = 'runtime requests retention'
                continue
            execution = p.parent/'execution.json'
            if execution.is_file():
                e = json.loads(execution.read_text())
                if 'publication failed' in str(e.get('error', '')).lower():
                    skipped[str(target)] = 'unpublished native evidence'
                    continue
            found[str(target)] = target
    for p in cache.iterdir():
        if str(p) not in found and str(p) not in skipped:
            skipped[str(p)] = 'not identified by archived runtime records'
    return [p for n, p in sorted(found.items()) if n not in skipped
            and (p.exists() or p.is_symlink())]


def cleanup(archive, destination, options, *, protected=()):
    """Resume a single cleanup transaction; never clear a replacement directory.

    Moving each target to its journaled sibling first means an interruption
    cannot cause the next invocation to delete newly created cache contents.
    """
    validate(options)
    archive = Path(archive).resolve(strict=True)
    destination = Path(destination).resolve(strict=True)
    home = Path.home().resolve(strict=True)
    protected = [archive, destination, *map(Path, protected)]
    path = destination/'cache_cleanup.json'
    if not any(options.values()):
        return {'status': 'disabled'}
    if not shutil.rmtree.avoids_symlink_attacks:
        raise RuntimeError('This Python lacks descriptor-safe recursive deletion; cleanup refused.')
    # Serialize cleanup across campaigns, not just invocations of this handoff.
    fd = os.open(home/'.kinbot_cache_cleanup.lock', os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, 'r+') as lock:
        if os.fstat(lock.fileno()).st_uid != os.getuid():
            raise RuntimeError('Cache cleanup lock is not owned by this account.')
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if path.exists():
            journal = json.loads(path.read_text())
            if (journal['archive'] != str(archive) or journal['home'] != str(home)
                    or journal['options'] != options):
                raise ValueError('Cache cleanup journal differs from requested cleanup.')
            if journal['status'] == 'complete':
                return journal
        else:
            skipped = {}
            targets = _scratch_candidates(archive, home, skipped) if options.get(
                'remove_verified_home_scratch') else []
            wfu = home/'wfu'
            if options.get('empty_home_wfu') and (wfu.exists() or wfu.is_symlink()):
                targets.insert(0, wfu)
            journal = dict(schema=1, status='planned', archive=str(archive), home=str(home),
                           options=options, skipped=skipped, entries=[])
            for p in targets:
                _check_tree(p, protected)
                identity = _identity(p)
                quarantine = p.with_name('.kinbot-retired-cache-'+uuid.uuid4().hex)
                journal['entries'].append(dict(path=str(p), retired=str(quarantine),
                    identity=identity, recreate=p == wfu, status='planned'))
            atomic_json(path, journal)
        targets = [Path(e[k]) for e in journal['entries'] for k in ('path', 'retired')]
        assert_idle(targets)
        for entry in journal['entries']:
            if entry['status'] == 'removed':
                continue
            source, retired = Path(entry['path']), Path(entry['retired'])
            if retired.exists() or retired.is_symlink():
                if _identity(retired) != entry['identity']:
                    raise RuntimeError(f'Retired cache identity changed: {retired}')
            elif entry['status'] == 'planned':
                if _identity(source) != entry['identity']:
                    raise RuntimeError(f'Cache identity changed since planning: {source}')
                _check_tree(source, protected)
                source.rename(retired)
            elif entry['status'] != 'removing':
                raise RuntimeError(f'Missing retired cache: {retired}')
            # A crash after rename is recoverable via the recorded inode. Do
            # not touch any new source directory, now or on later invocations.
            if entry['recreate'] and not source.exists():
                source.mkdir(mode=entry['identity']['mode'])
            entry['status'] = 'removing'
            atomic_json(path, journal)
            if retired.exists():
                _check_tree(retired, protected)
                assert_idle(targets)
                print(f'Removing retired cache: {source}', flush=True)
                shutil.rmtree(retired)
            entry['status'] = 'removed'
            atomic_json(path, journal)
        journal['status'] = 'complete'
        atomic_json(path, journal)
        return journal
