ml purge

newgrp - ggaussian
module load Gaussian/g16_E.01-intel-2018a

module load Python/3.6.4-intel-2018a

## purge all loaded modules
module purge

## change to the UGent gaussian user group and load modules
newgrp - ggaussian
module load Gaussian/g16_E.01-intel-2018a
module load Python/3.6.4-intel-2018a

## store directory name of the submit dir
ORIGDIR=$SLURM_SUBMIT_DIR/

## make a directory on the VSC scratch
WORKDIR=$VSC_SCRATCH_NODE/$SLURM_JOB_ID
mkdir -p $WORKDIR/hir
mkdir -p $WORKDIR/conf

## start the Gaussian calculation
cd $WORKDIR
cp $ORIGDIR/{python_file} $WORKDIR/{python_file}
cp $ORIGDIR/{name}.chk $WORKDIR/{name}.chk

python {python_file} {arguments}
formchk {name}.chk

## copy results
cp {name}.log $ORIGDIR/{name}.log
cp {name}_prod.log $ORIGDIR/{name}_prod.log
cp {name}.com $ORIGDIR/{name}.com
cp {name}_prod.com $ORIGDIR/{name}_prod.com
cp {name}.chk $ORIGDIR/{name}.chk
cp {name}.fchk $ORIGDIR/{name}.fchk
cd $ORIGDIR
rm -rf $WORKDIR

