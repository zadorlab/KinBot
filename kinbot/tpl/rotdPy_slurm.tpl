#!/usr/bin/env bash
#SBATCH --job-name=rotd-s{surf_id}-f{face_id}-n{samp_id}
#SBATCH --nodes=1
#SBATCH --ntasks={procs}
#SBATCH --mem={mem}M
#SBATCH --time=@WALLTIME@
@PARTITION_DIRECTIVE@
@EXCLUSIVE_DIRECTIVE@
#SBATCH --output=surf{surf_id}_face{face_id}_samp{samp_id}.stdout
#SBATCH --error=surf{surf_id}_face{face_id}_samp{samp_id}.err

set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
# Every generated rotdPy sample is a single-node Molpro job.  Preserve a site
# override; otherwise keep Intel MPI off unavailable PSM3/OFI network devices.
export I_MPI_FABRICS="${I_MPI_FABRICS:-shm}"
export MPLCONFIGDIR="$PWD/.matplotlib"
mkdir -p "$MPLCONFIGDIR"
@PYTHON@ surf{surf_id}_face{face_id}_samp{samp_id}.py
