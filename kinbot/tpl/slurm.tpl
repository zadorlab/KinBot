#! /bin/bash

#SBATCH -N 1
#SBATCH -c {ppn}
#SBATCH --partition={queue_name}
#SBATCH -o {errdir}/{name}.stdout
#SBATCH -e {errdir}/{name}.err
#SBATCH --mem=10gb
{slurm_feature}
{slurm_exclude_directive}
{walltime_directive}

export OMP_NUM_THREADS={ppn}
