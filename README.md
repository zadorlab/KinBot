# KinBot: Automated reaction pathway search for gas-phase molecules with C, H, O, N and S atoms

## Description
This repository contains the KinBot code version 2.0,
a tool for automatically searching for reactions on the potential energy surface.

## How to Install
Make sure all dependencies are correctly installed (see the manual in the docs/ directory for a list of dependencies). 

Clone the project on your machine and go into the KinBot/ directory. Run the following:

    python setup.py build
    python setup.py install
    
Alternatively, you can also install the latest stable version (which might be slightly behind compared to what is on GitHub) using conda:

    conda install -c rdkit -c openbabel -zadorlab kinbot


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
This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE organizations, the Office of Science and the National Nuclear Security Administration, responsible for the planning and preparation of a capable exascale ecosystem including software, applications, hardware, advanced system engineering, and early test bed platforms to support the nation's exascale computing imperative. RVdV was also supported by the AITSTME project as part of the Predictive Theory and Modeling component of the Materials Genome Initiative. Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and Engineering Solutions of Sandia, LLC., a wholly owned subsidiary of Honeywell International, Inc., for the U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525. The views expressed in the article do not necessarily represent the views of the U.S. Department of Energy or the United States Government. 
