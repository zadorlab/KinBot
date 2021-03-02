#! /bin/bash

#SBATCH -N 1
#SBATCH -c {ppn}
#SBATCH -q {queue_name}
#SBATCH -o {dir}/{name}.stdout
#SBATCH -e {dir}/{name}.err
{slurm_feature}

export OMP_NUM_THREADS={ppn}
