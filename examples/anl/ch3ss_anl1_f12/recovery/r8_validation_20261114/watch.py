"""Read-only panel/Slurm snapshot; Ctrl+C stops this watcher only."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import getpass
import json
from pathlib import Path
import socket
import subprocess
import time


def read(path):
    try: return json.loads(path.read_text())
    except (OSError, ValueError): return {}


def native_reason(directory, ident):
    """Compact native cause, rather than the launcher's abort traceback."""
    text = ''
    for suffix in ('.log', '.out'):
        p = directory/(ident+suffix)
        if not p.is_file(): continue
        with p.open('rb') as stream:
            stream.seek(max(0,p.stat().st_size-16000))
            text += stream.read().decode(errors='replace')
    for marker, reason in [('RxnCyc: Invalid step size found!', 'invalid_DVV_time_step'),
                           ('RHF not converged', 'ROHF_not_converged'),
                           ('The problem occurs in ccsd', 'CC_not_converged'),
                           ('Maximum number of corrector steps', 'IRC_corrector_limit')]:
        if marker in text: return reason
    return 'inspect_native_output'


def snapshot(root):
    root = Path(root).resolve()
    print(datetime.now(timezone.utc).isoformat(), socket.gethostname())
    owner = read(root/'native_owner.json'); pid = owner.get('pid')
    live = False
    if isinstance(pid, int) and pid > 0:
        try:
            args = (Path('/proc')/str(pid)/'cmdline').read_bytes().split(b'\0')
            live = (str(root).encode() in args and b'run' in args
                    and any(a.endswith((b'/resume.py', b'/campaign.py')) for a in args))
        except OSError: pass
    print('OWNER', 'process present' if live else 'process absent/unverified',
          'pid', pid, 'saved status', owner.get('status'), owner.get('error',''))
    progress = root/'progress.json'
    if progress.is_file():
        p = read(progress)
        print('PROGRESS age_seconds', round(time.time()-progress.stat().st_mtime),
              'finished', p.get('execution_finished'), 'infrastructure', p.get('infrastructure',{}))
    jobs = {}; available = False
    try:
        query = subprocess.run(['squeue','-h','-u',getpass.getuser(),'-o','%i|%T|%M|%Z'],
                               text=True,capture_output=True,timeout=15,check=True)
        for line in query.stdout.splitlines():
            fields = line.strip().split('|',3)
            if len(fields)==4: jobs[fields[0]] = fields[1:]
        available = True
    except (OSError,subprocess.SubprocessError) as error:
        print('QUEUE unavailable:', error)
    names = read(root/'stage.json').get('workflows',[])
    plan = read(root/'reconciliation_v1.json')
    extra = dict(plan.get('curvature',{}))
    for filename in ('reconciliation_repairs_v2.json', 'reconciliation_electronic_v3.json', 'reconciliation_departure_v3.json'):
        repair = read(root/filename)
        plan['irc'] = plan.get('irc',[]) + repair.get('irc',[])
        extra.update(repair.get('curvature',{}))
    pilots = {(p['graph'],p['task']['id']) for p in plan.get('irc',[])}
    counts = Counter()
    for name in [*names, *extra]:
        state = read(root/name/'state.json')
        for ident,entry in state.get('tasks',{}).items():
            if name.startswith('irc_') and ident!='hessian': counts[entry.get('status','unknown')] += 1
            if (name,ident) not in pilots and name not in extra: continue
            job = str(entry.get('job_id',''))
            execution = read(root/name/'tasks'/ident/'execution.json')
            native = execution.get('status') if str(execution.get('slurm_job_id'))==job else None
            scheduler = jobs.get(job,['absent from queue' if available else 'unknown'])[0]
            print(name,ident,'state='+str(entry.get('status')), 'job='+job,
                  'queue='+scheduler,'native='+str(native))
            if native == 'failed':
                print('  cause='+native_reason(root/name/'tasks'/ident,ident))
    print('IRC task states',dict(counts))
    if available:
        # The native directory is under RUN/validation_DATE/native.
        run = root.parent.parent
        running = [(job,*fields) for job,fields in jobs.items() if Path(fields[2]).is_relative_to(run)]
        print('ALL R8 allocations',len(running))
        for row in running: print('JOB',*row)
    held = sum(len(x) for x in read(progress).get('held_for_review',{}).values())
    print('HELD unsubmitted tasks',held,'scientific approval: false')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--interval',type=int,default=0,help='Repeat every N seconds (minimum 30).')
    args = parser.parse_args()
    if args.interval and args.interval < 30: parser.error('Use an interval of at least 30 seconds.')
    try:
        while True:
            snapshot(args.root)
            if not args.interval: return
            time.sleep(args.interval)
    except KeyboardInterrupt: pass


if __name__ == '__main__': main()
