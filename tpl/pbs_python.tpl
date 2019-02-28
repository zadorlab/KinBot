
newgrp - ggaussian
module load Gaussian/g16_E.01-intel-2018a

module load Python/3.6.4-intel-2018a
cd ${{PBS_O_WORKDIR}}
python {python_file} {arguments}

formchk {name}.chk

