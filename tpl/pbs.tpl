#! /bin/bash -f
#PBS -N {name}
#PBS -l nodes=1:ppn={ppn}
#PBS -l walltime=11:59:59
#PBS -o {dir}/$PBS_JOBNAME.stdout
#PBS -e {dir}/$PBS_JOBNAME.err
#PBS -m n

## {queue_name}

