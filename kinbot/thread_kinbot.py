import sys,os
import subprocess
import threading, time

def run_threads(jobs, name, max_running=10, pes=0):
    """
    This method runs all the tests, instead of using the python threading
    it submits python runs to the head node
    
    jobs is a dictionary of paths to the directories where the input files are (keys)
    and the names of the json input files (values)
    
    max_running is the maximum threads that are allowed to run simultaneously 
    """

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
    command = ['ps','-u','root','-N','-o','pid,s,user,%cpu,%mem,etime,args']
    process = subprocess.Popen(command,shell=False,stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    out,err = process.communicate()
    lines = out.split('\n')
    for line in lines:
        if len(line)> 0:
            if '%i'%pid == line.split()[0]:
                return 1
    return 0


def submit_job(job, inpfile, pes):
    """
    Submit a kinbot run using subprocess and return the pid
    """
    if pes:
        command = ["python","/home/jzador/KinBot/kinbot/pes.py",inpfile,"&"]
    else:
        command = ["python","/home/jzador/KinBot/kinbot/kb.py",inpfile,"&"]
    process = subprocess.Popen(command,cwd = os.path.expanduser(job), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    time.sleep(1)
    pid = process.pid
    return pid 

