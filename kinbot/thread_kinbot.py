import sys,os
import subprocess
import time
from pathlib import Path

_processes = {}

def run_threads(jobs, name, max_running=10, pes=0):
    """
    This method runs all the tests, instead of using the python threading
    it submits python runs to the head node
    
    jobs is a dictionary of paths to the directories where the input files are (keys)
    and the names of the json input files (values)
    
    max_running is the maximum threads that are allowed to run simultaneously 
    """

    if type(max_running) is not int or max_running < 1:
        raise ValueError('max_running must be a positive integer.')
    running = []
    finished = []
    pids = {}
    while len(finished) < len(jobs):
        
        while len(running) < max_running and len(running) + len(finished) < len(jobs):
            # start a new job
            job = sorted(jobs.keys())[len(running) + len(finished)]
            pid = submit_job(job, jobs[job], pes)
            pids[job] = pid
            running.append(job)
        
        #check if a thread is done
        for job in running:
            if not check_status(job,pids[job]):
                finished.append(job)
        
        #remove the finished threads
        for job in finished: 
            if job in running:
                running.remove(job)
        
        with open('%s_thread_summary.txt'%name,'w+') as f:
            f.write('Total\t\t%i\n'%len(jobs))
            f.write('Running\t\t%i\n'%len(running))
            f.write('Finished\t%i\n\n'%len(finished))
            
            
            f.write('Running:\n')
            for job in running:
                f.write('\t%s\n'%job)
                
            f.write('\nFinished:\n')
            for job in finished: 
                f.write('\t%s\n'%job)

        time.sleep(1)

def check_status(job,pid):
    process = _processes[pid]
    code = process.poll()
    if code is None:
        return 1
    if code:
        raise RuntimeError(f'KinBot in {job} exited with status {code}; '
                           'inspect thread_kinbot.out and thread_kinbot.err.')
    return 0


def submit_job(job, inpfile, pes):
    """
    Submit a kinbot run using subprocess and return the pid
    """
    command = [sys.executable, '-m', 'kinbot.pes' if pes else 'kinbot.kb', inpfile]
    directory = Path(job).expanduser()
    # PIPEs that nobody drains can deadlock a verbose worker.
    with (directory/'thread_kinbot.out').open('ab') as output, \
            (directory/'thread_kinbot.err').open('ab') as errors:
        process = subprocess.Popen(command, cwd=directory, stdout=output,
                                   stdin=subprocess.DEVNULL, stderr=errors)
    _processes[process.pid] = process
    pid = process.pid
    return pid 
