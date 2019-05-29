#! /bin/bash

#SBATCH -N 1
#SBATCH -c {ppn}
#SBATCH -t 05:00:00
#SBATCH -o {dir}/{name}.stdout
#SBATCH -e {dir}/{name}.err
{slurm_feature}

# {queue_name}
