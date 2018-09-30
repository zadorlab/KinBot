#! /bin/bash -f
#PBS -N _name_
#PBS -l nodes=1:ppn=_ppn_
#PBS -q medium
#PBS -o perm/$PBS_JOBNAME.stdout
#PBS -e perm/$PBS_JOBNAME.err
#PBS -m n

cd ${PBS_O_WORKDIR}
_code_ _name_._exti_ > _name_._exto_
echo "done" >> _name_._exto_

