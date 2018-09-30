###################################################
##                                               ##
## This file is part of the KinBot code v2.0     ##
##                                               ##
## The contents are covered by the terms of the  ##
## BSD 3-clause license included in the LICENSE  ##
## file, found at the root.                      ##
##                                               ##
## Copyright 2018 National Technology &          ##
## Engineering Solutions of Sandia, LLC (NTESS). ##
## Under the terms of Contract DE-NA0003525 with ##
## NTESS, the U.S. Government retains certain    ##
## rights to this software.                      ##
##                                               ##
## Authors:                                      ##
##   Judit Zador                                 ##
##   Ruben Van de Vijver                         ##
##                                               ##
###################################################
import sys,os
import subprocess
import threading, time
import kinbot

def run_threads(jobs, name, max_running = 10):
    """
    This method runs all the tests, instead of using the python threading
    it submits python runs to the head node
    
    jobs is a list of paths to the directories where the input files are
    max_running is the maximum threads that are allowed to run simultaneously 
    """

    running = []
    finished = []
    pids = {}
    while len(finished) < len(jobs):
        
        while len(running) < max_running and len(running) + len(finished) < len(jobs):
            # start a new job
            job = jobs[len(running) + len(finished)]
            pid = submit_job(job)
            pids[job] = pid
            
            #t = threading.Thread(name=job,target=run_kinbot,args=(job,))
            #t.start()
            running.append(job)
        
        #check if a thread is done
        for job in running:
            if not check_status(job,pids[job]):
                finished.append(job)
        
        #remove the finished threads
        for job in finished: 
            if job in running:
                running.remove(job)
        
        f = open('%s_thread_summary.txt'%name,'w+')
        f.write('Total\t\t%i\n'%len(jobs))
        f.write('Running\t\t%i\n'%len(running))
        f.write('Finished\t%i\n\n'%len(finished))
        
        for job in finished: 
            f.write('\t%s\n'%job)
        
        f.close()

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


def submit_job(job):
    """
    Submit a kinbot run usung subprocess and return the pid
    """
    command = ["python","kinbot.py","%s"%job[2:],"&"]
    process = subprocess.Popen(command,stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    time.sleep(1)
    pid = process.pid
    return pid 

