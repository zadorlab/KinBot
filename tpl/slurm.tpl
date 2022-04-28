#! /bin/bash

#SBATCH -N 1
#SBATCH -c {ppn}
#SBATCH -q {queue_name}
#SBATCH -o {errdir}/{name}.stdout
#SBATCH -e {errdir}/{name}.err
{slurm_feature}

export OMP_NUM_THREADS={ppn}
