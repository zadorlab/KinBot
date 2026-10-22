"""Bounded, work-conserving turns between waiting PES/branch submissions.

Called only under run_policy's admission lock. Tickets never own allocations;
Slurm and the existing admission ledger remain the capacity authority.
"""
import hashlib
import json
import os
from pathlib import Path
import socket
import time

from kinbot.execution import atomic_json

LEASE_SECONDS = 120.


def group(root, directory):
    scope = root
    branch = 'search'
    for parent in (directory, *directory.parents):
        if not parent.is_relative_to(root):
            break
        if parent.name == 'vrctst':
            branch = 'vrc'
        elif (parent/'network_request.json').is_file() or parent.name == 'anl_automatic':
            branch = 'anl'
        if (parent/'pes_exploration.json').is_file():
            scope = parent
            break
    return str(scope.relative_to(root)) + ':' + branch


def take_turn(root, directory, command, vacancies):
    """Register/renew interest, returning whether this caller may submit now."""
    path = root/'.kinbot_policy_waiters.json'
    state = json.loads(path.read_text()) if path.exists() else dict(schema=1, waiting={}, served={})
    if state.get('schema') != 1:
        raise ValueError('Unsupported shared-pool waiter schema.')
    now, host = time.time(), socket.gethostname()
    waiting, served = state['waiting'], state['served']
    for token, item in list(waiting.items()):
        stale = now-item['heartbeat'] > LEASE_SECONDS or now < item['heartbeat']
        if not stale and item['host'] == host:
            try:
                os.kill(item['pid'], 0)
            except ProcessLookupError:
                stale = True
            except PermissionError:
                pass
        if stale:
            del waiting[token]
    token = hashlib.sha256(json.dumps([str(directory), list(command)]).encode()).hexdigest()
    item = waiting.setdefault(token, dict(group=group(root,directory), first=now,
                                         directory=str(directory), command=list(command)))
    item.update(heartbeat=now, pid=os.getpid(), host=host)
    groups = {}
    for key, value in sorted(waiting.items(), key=lambda kv:(kv[1]['first'], kv[0])):
        groups.setdefault(value['group'], []).append(key)
    ordered = sorted(groups, key=lambda g:(served.get(g, 0), waiting[groups[g][0]]['first'], g))
    eligible = []
    while ordered and len(eligible) < vacancies:
        for g in ordered:
            eligible.append(groups[g].pop(0))
            if len(eligible) >= vacancies:
                break
        ordered = [g for g in ordered if groups[g]]
    allowed = token in eligible
    if allowed:
        state['turn'] = state.get('turn', 0) + 1
        served[item['group']] = state['turn']
        del waiting[token]
    atomic_json(path, state)
    return allowed
