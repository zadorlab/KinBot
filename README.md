[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/zadorlab/KinBot)

# KinBot: Automated reaction pathway search for gas-phase molecules

<p>
    <img src="graphics/kinbot_logo_V2.png" width="220" height="240" />
</p>

## Description
This repository contains the KinBot code version 2.0,
a tool for automatically searching for reactions on the potential energy surface.

If you are using this tool in scientific publications, please reference this git repo and the following publication:

Ruben Van de Vijver, Judit Zádor: KinBot: Automated stationary point search on potential energy surfaces, Computer Physics Communication, 2019, 106947 https://doi.org/10.1016/j.cpc.2019.106947

## How to Install
Make sure all dependencies are correctly installed. The dependencies are lister here https://github.com/zadorlab/KinBot/wiki/Setting-Up-KinBot-on-Your-System 

Clone the project to the place where you want to run it. Make sure you switch to the latest version, e.g., 2.0.5:

    git branch 2.0.5

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

## Papers using KinBot

* Van de Vijver, R., Zádor, J.: KinBot: Automated stationary point search on potential energy surfaces. Computer Physics Communications, 2020 248
* Satya Prakash Joshi, Prasenjit Seal, Timo Theodor Pekkanen, Raimo Sakari Timonen and Arkke Johannes Eskola: Direct Kinetic Measurements and Master Equation Modelling of the Unimolecular Decomposition of Resonantly-Stabilized CH2CHCHC(O)OCH3 Radical and an Upper Limit Determination for CH2CHCHC(O)OCH3+O2 Reaction. Z. Phys. Chem., 2020, https://doi.org/10.1515/zpch-2020-1612
* Katherine S.Lockwood, Nicole J.Labbe: Insights on keto-hydroperoxide formation from O2 addition to the beta-tetrahydrofuran radical. Proceedings of the Combustion Institute, Volume 38, Issue 1, 2021, Pages 533-541 https://doi.org/10.1016/j.proci.2020.06.357
* Leonid Sheps, Amanda L. Dewyer, Maria Demireva, and Judit Zádor: Quantitative Detection of Products and Radical Intermediates in Low-Temperature Oxidation of Cyclopentane. J. Phys. Chem. A 2021, 125, 20, 4467–4479, https://doi.org/10.1021/acs.jpca.1c02001

Old version:
* Van de Vijver, R., Van Geem, K. M., Marin, G. B., Zádor, J.: Decomposition and isomerization of 1-pentanol radicals and the pyrolysis of 1-pentanol. Combustion and Flame, 2018 196 500-514.
* Grambow, C. A., Jamal, A., Li, Y.-P., Green, W. H., Zádor, J., Suleimanov, Y. V.: New unimolecular reaction pathways of a g-ketohydroperoxide from combined application of automated reaction discovery methods. Journal of the American Chemical Society, 2018 140 1035-1048.
* Rotavera, B., Savee, J. D., Antonov, I. O., Caravan, R. L., Sheps, L., Osborn, D. L., Zádor, J., Taatjes, C. A.: Influence of oxygenation in cyclic hydrocarbons on chain-termination reactions from R + O2: tetrahydropyran and cyclohexane. Proceedings of the Combustion Institute, 2017 36 597-606.
* Antonov, I. O., Zádor, J., Rotavera, B., Papajak, E., Osborn, D. L., Taatjes, C. A., Sheps, L.: Pressure-Dependent Competition among Reaction Pathways from First- and Second-O2 Additions in the Low-Temperature Oxidation of Tetrahydrofuran. Journal of Physical Chemistry A 2016, 120 6582-6595.
* Antonov, I. O., Kwok, J., Zádor, J., Sheps, L.: OH + 2-butene: A combined experimental and theoretical study in the 300-800 K temperature range. Journal of Physical Chemistry A, 2015 119 7742-7752.
* Zádor, J., Miller, J.A.: Adventures on the C3H5O potential energy surface: OH + propyne, OH + allene and related reactions. Proceedings of the Combustion Institute, 2015 35 181-188.
* Rotavera, B., Zádor, J., Welz, O., Sheps, L., Scheer, A.M., Savee, J.D., Ali, M.A., Lee, T.S., Simmons, B.A., Osborn, D.L., Violi, A., Taatjes, C.A.: Photoionization mass spectrometric measurements of initial reaction pathways in low-temperature oxidation of 2,5-dimethylhexane. Journal of Physical Chemistry A, 2014 44 10188-10200.

## Acknowledgement
This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE organizations, the Office of Science and the National Nuclear Security Administration, responsible for the planning and preparation of a capable exascale ecosystem including software, applications, hardware, advanced system engineering, and early test bed platforms to support the nation's exascale computing imperative. RVdV was also supported by the AITSTME project as part of the Predictive Theory and Modeling component of the Materials Genome Initiative. Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and Engineering Solutions of Sandia, LLC., a wholly owned subsidiary of Honeywell International, Inc., for the U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525. 
