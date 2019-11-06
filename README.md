[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/zadorlab/KinBot)

# KinBot: Automated reaction pathway search for gas-phase molecules with C, H, O, N and S atoms

<p>
    <img src="graphics/kinbot_logo_V2.png" width="220" height="240" />
</p>

## Description
This repository contains the KinBot code version 2.0,
a tool for automatically searching for reactions on the potential energy surface.

If you are using this tool in scientific publications, please reference this git repo and the following publication:

Ruben Van de Vijver, Judit Zádor

KinBot: Automated stationary point search on potential energy surfaces

Computer Physics Communication, 2019 (in press currently, please check if there is a page number assigned)

https://doi.org/10.1016/j.cpc.2019.106947

## How to Install
Make sure all dependencies are correctly installed (see the manual in the docs/ directory for a list of dependencies). 

Clone the project to the place where you want to run it. Make sure you switch to the latest version, e.g., 2.0.1:

    git branch 2.0.1

You can find the latest stable version's tag if you click on the Branch button above on this page.

In your local space go into the KinBot/ directory. Run the following:

    python setup.py build
    python setup.py install
    
If you do not have admin priveleges, you might have to run

    python setup.py build
    python setup.py install --user
    
Moreover, if you plan to modify the code, you need to install it as:

    python setup.py build
    python setup.py develop --user

Alternatively, you can also install the latest stable version (however, this is currently way behind compared to what is on GitHub, so currently the github installation is preferred) using conda:

    conda install -c rdkit -c openbabel -c zadorlab kinbot


## How to Run
To run KinBot (which will only explore one well), make an input file (e.g. input.json) and run:

    kinbot input.json

To run a full PES search, make an input file (e.g. input.json) and run:

    pes input.json

You can find additional command line arguments in the manual. 

## Documentation
See [wiki](https://github.com/zadorlab/KinBot/wiki).

## List of files in this project
See [list](https://github.com/zadorlab/KinBot/wiki/KinBot-file-structure).

## Authors
* Judit Zador (jzador@sandia.gov)
* Ruben Van de Vijver (Ruben.VandeVijver@UGent.be)

## Website
https://kinbot.sandia.gov

## Acknowledgement
This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE organizations, the Office of Science and the National Nuclear Security Administration, responsible for the planning and preparation of a capable exascale ecosystem including software, applications, hardware, advanced system engineering, and early test bed platforms to support the nation's exascale computing imperative. RVdV was also supported by the AITSTME project as part of the Predictive Theory and Modeling component of the Materials Genome Initiative. Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and Engineering Solutions of Sandia, LLC., a wholly owned subsidiary of Honeywell International, Inc., for the U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525. 
