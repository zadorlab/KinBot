# KinBot: Automated reaction pathway search for gas-phase molecules with C, H, O and S atoms

## Description
This repository contains the KinBot code version 2.0,
a tool for automatically searching for reactions on the potential energy surface.

## How to Install
Clone the project on your machine and go into the KinBot/ directory. Run the following:

    python setup.py build
    python setup.py install

## How to Run
To run KinBot (which will only explore one well), make an input file (e.g. input.json) and run:

    kinbot input.json

To run a full PES search, make an input file (e.g. input.json) and run:

    pes input.json

You can find additional command line arguments in the manual. 

## Documentation
See the manual in the docs/ directory. 

## Authors
* Judit Zador (jzador@sandia.gov)
* Ruben Van de Vijver (vandevijver.ruben@gmail.com)

