#!/usr/bin/env python3
"""Read-only owner/queue/native-tail capture for final evidence reconciliation.

No KinBot dispatcher is imported; no task, policy or controller state is changed.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import socket
import subprocess


def read(path):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError) as error:
        return {'read_error': str(error)}


def tail(path, limit=5000):
    try:
        with path.open('rb') as stream:
            stream.seek(max(0, path.stat().st_size-limit))
            return stream.read(limit).decode(errors='replace')
    except OSError as error:
        return str(error)


def command(args):
    try:
        result = subprocess.run(args, text=True, capture_output=True, timeout=15, check=False)
        return dict(args=args, returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
    except (OSError, subprocess.TimeoutExpired) as error:
        return dict(args=args, error=str(error))


def snapshot(root):
    native=root/'native';remaining=root/'remaining_l2'
    result=dict(captured_utc=datetime.now(timezone.utc).isoformat(),host=socket.gethostname(),root=str(root))
    result['saved']={str(p.relative_to(root)):read(p) for p in
        [native/'completion.json',native/'native_owner.json',remaining/'controller.json',remaining/'results.json']}
    processes=[]
    proc=Path('/proc')
    for p in proc.iterdir() if proc.is_dir() else []:
        if not p.name.isdecimal():continue
        try:
            if p.stat().st_uid!=os.getuid():continue
            args=(p/'cmdline').read_bytes().decode(errors='replace').replace('\0',' ').strip()
            if str(root) not in args:continue
            processes.append(dict(pid=int(p.name),args=args,wchan=tail(p/'wchan',256),
                status=tail(p/'status',3500),stack=tail(p/'stack',1500)))
        except OSError:continue
    result['matching_processes']=processes
    result['queue']=command(['squeue','--me','-h','-o','%i|%T|%M|%N|%Z|%R'])
    result['nodes']=command(['sinfo','-N','-h','-p','month-long-cpu','-o','%N|%T|%c|%m'])
    for parent in [root,*root.parents]:
        p=parent/'.kinbot_resource_policy.json'
        if p.is_file():
            result['resource_policy']=dict(path=str(p),policy=read(p),
                admissions=read(parent/'.kinbot_policy_admissions.json'))
            break
    result['owner_log_tails']={p.name:tail(p) for p in sorted(root.glob('l2.*.log'))[-3:]}
    jobs=[]; failures={}
    for graph,task in [('cc_si_W1','refine'),('sp_paper_ts1_F','qz')]:
        d=native/graph/'tasks'/task;execution=read(d/'execution.json')
        job=str(execution.get('slurm_job_id',''))
        if job.isdecimal():jobs.append(job)
        failures[graph+'/'+task]=dict(execution=execution,files={p.name:tail(p) for p in d.iterdir()
            if p.is_file() and p.suffix in ('.out','.log','.err','.stderr')}) if d.is_dir() else {'missing':str(d)}
    result['incomplete_native_tails']=failures
    if jobs:result['accounting']=command(['sacct','-j',','.join(jobs),'--noheader','--parsable2',
        '--format=JobID,State,ExitCode,DerivedExitCode,NodeList,Reason,MaxRSS'])
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('validation_root',type=Path)
    args=parser.parse_args()
    print(json.dumps(snapshot(args.validation_root.resolve()),indent=2))
