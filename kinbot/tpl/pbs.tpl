#! /bin/bash -f
#PBS -N {name}
#PBS -l nodes=1:ppn={ppn}
#PBS -q {queue_name}
#PBS -o {errdir}/$PBS_JOBNAME.stdout
#PBS -e {errdir}/$PBS_JOBNAME.err
#PBS -m n

export OMP_NUM_THREADS={ppn}

