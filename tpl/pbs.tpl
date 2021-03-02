#! /bin/bash -f
#PBS -N {name}
#PBS -l nodes=1:ppn={ppn}
#PBS -q {queue_name}
#PBS -o {dir}/$PBS_JOBNAME.stdout
#PBS -e {dir}/$PBS_JOBNAME.err
#PBS -m n

export OMP_NUM_THREADS={ppn}

