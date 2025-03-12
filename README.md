[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/zadorlab/KinBot)

# KinBot: Automated Reaction Kinetics of Gas-Phase Organic Species over Multiwell Potential Energy Surfaces

<p>
    <img src="https://raw.githubusercontent.com/zadorlab/KinBot/master/graphics/kinbot_logo_V2.png" width="220" height="240" />
</p>

## Description
This repository contains the KinBot code version 2.2.1,
a tool for automatically searching for reactions on the potential energy surface.

If you are using this tool in scientific publications, please reference the following publications:

* Ruben Van de Vijver, Judit Zádor: KinBot: _Automated stationary point search on potential energy surfaces_, Comp. Phys. Comm., **2019**, 248, 106947. https://doi.org/10.1016/j.cpc.2019.106947
```
@article{Vijver2020,
   author = {Van de Vijver, Ruben and Z\'ador, Judit},
   title = {KinBot: Automated stationary point search on potential energy surfaces},
   journal = {Comput. Phys. Commun.},
   volume = {248},
   pages = {106947},
   year = {2020},
   type = {Journal Article}
}
```
* Judit Zádor, Carles Martí, Ruben Van de Vijver, Sommer L. Johansen, Yoona Yang, Hope A. Michelsen, Habib N. Najm: _Automated reaction kinetics of gas-phase organic species over multiwell potential energy surfaces_, J. Phys. Chem. A, **2023**, 127, 565–588. https://doi.org/10.1021/acs.jpca.2c06558

```
@article{Zador2022,
   author = {Z\'ador, Judit and Mart\'i, Carles and Van de Vijver, Ruben and Johansen, Sommer L. and Yang, Yoona and Michelsen, Hope A. and Najm, Habib N.},
   title = {Automated reaction kinetics of gas-phase organic species over multiwell potential energy surfaces},
   journal = {J. Phys. Chem. A},
   volume = {127},
   pages = {565-588},
   year = {2023},
   type = {Journal Article}
}
```

We appreciate if you send us the DOI of your published paper that used KinBot, so we can feature it here below.

## How to Install

KinBot can be installed both in three different ways, from the PyPI index (`pip install`), from the conda-forge repo (`conda install`) or by cloning this github repo and then install it locally.

### PyPI

    pip install kinbot

> **Note**
>  KinBot only works with Python >= 3.10.

### conda-forge

    conda install -c conda-forge kinbot

### From Github

If you want to have the very last version of KinBot without waiting for a 
release or you want to modify KinBot acccording to your needs you can clone the project 
from github:

    git clone git@github.com:zadorlab/KinBot.git

and then, from within the KinBot directory produced after cloning, type:

    pip install -e .
 
> **Note**
> If you want to modify KinBot yourself it's better to fork the project 
> into your own repository and then clone it.

## How to Run
To run a single-well exploration of KinBot, make an input file (e.g. input.json) and run:

    kinbot input.json

To run a full PES search, make an input file (e.g. input.json) and run:

    pes input.json

You can find additional command line arguments in the manual. 

## Documentation
See the [wiki](https://github.com/zadorlab/KinBot/wiki) for keywords, and our [tutorial](https://hackmd.io/@jzador/kinbot_workshop_2023#/) for a more hands-on introduction to the code.

## List of files in this project
See [list](https://github.com/zadorlab/KinBot/wiki/KinBot-file-structure).

## Authors
* Judit Zádor (jzador@sandia.gov)
* Ruben Van de Vijver 
* Amanda Dewyer
* Carles Martí
* Clément Soulié (csoulie@sandia.gov)

## Papers using KinBot
1. Hansen, N., Gaiser, N., Bierkandt, T., Oßwald, P., Köhler, M., Zádor, J., Hemberger, P.: _Identification of dihydropentalenes as products of the molecular-weight growth reaction of cyclopentadienyl plus propargyl_. J. Phys. Chem. A, **2025** 129 1714-1725. https://doi.org/10.1021/acs.jpca.4c06549
2. Almeida, T. G., Martí, C., Kurtén, T., Zádor, J., Johansen, S. L.: _Theoretical analysis of the OH-Initiated atmospheric oxidation reactions of imidazole_. Phys. Chem. Chem. Phys., **2024** 26 23570-23587. https://doi.org/10.1039/D4CP02103G
3. Yuan, E. C.-Y., Kumar, A., Guan, X., Hermes, E. D., Rosen, A. S., Zádor, J., Head-Gordon, T., Blau, S. M.: _Analytical ab initio Hessian from a Deep Learning Potential for Transition State Optimization_. Nat. Comm., **2024** 15 8865. https://doi.org/10.1038/s41467-024-52481-5
4. Doner, A. C., Zádor, J., Rotavera, B.: _Stereoisomer-dependent rate coefficients and reaction mechanisms of 2-ethyloxetanylperoxy radicals_. Proc. Combust. Inst., **2024**, 40, 105578. https://doi.org/10.1016/j.proci.2024.105578
5. Hansen, N. A, Price, T. D., Filardi, L. R., Gurses, S. M., Zhou, W., Hansen, N., Osborn, D. L. Zádor, J., Kronawitter, C. X.: _The photoionization of methoxymethanol: Fingerprinting a reactive C2 oxygenate in a complex reactive mixture_. J. Chem. Phys., **2024**, 160, 124306. https://doi.org/10.1063/5.0197827
6. Martí, C., Devereux, C., Najm, H. N., Zádor, J.: _Evaluation of rate coefficients in the gas-phase using a machine learned potential_. J. Phys. Chem. A, **2024**, 128, 1958–1971. https://doi.org/10.1021/acs.jpca.3c07872
7. Lang, J., Foley, C. D., Thawoos, S., Behzadfar, A., Liu, Y., Zádor, J., Suits, A. G.: _Reaction dynamics of S(3P) with 1,3-butadiene and isoprene: Crossed beam scattering, low temperature flow experiments, and high-level electronic structure calculations_. Farad. Discuss., **2024**, 251, 550-572. https://doi.org/10.1039/D4FD00009A
8. Wang, D., Tian, Z.-Y., Zheng, Z.-H., Li, W., Wu, L.-N., Kuang, J.-J., Yang, J.-Z.: _Experimental and modeling study of the n, n-dimethylformamide pyrolysis at atmospheric pressure_. Combust. Flame, **2024**, 260, 113240. https://doi.org/10.1016/j.combustflame.2023.113240
9. Doner, A. C., Zádor, J., Rotavera, B.: _Unimolecular reactions of 2,4-dimethyloxetanyl radicals._ J. Phys. Chem A, **2023**, 127, 2591–2600 https://doi.org/10.1021/acs.jpca.2c08290
10. Li, H., Lang, J., Foley, C. D., Zádor, J., Suits, A. G.: _Sulfur (3P) reaction with conjugated dienes gives cyclization to thiophenes under single collision conditions._ J. Phys. Chem. Letters, **2023**, 14, 7611–7617. https://doi.org/10.1021/acs.jpclett.3c01953
11. Martí, C., Michelsen, H. A., Najm, H. N., Zádor, J.: _Comprehensive kinetics on the C7H7 potential energy surface under combustion conditions._ J. Phys. Chem. A, **2023**, 127, 1941–1959. https://pubs.acs.org/doi/full/10.1021/acs.jpca.2c08035
12. Zádor, J, Martí, C., Van de Vijver, R., Johansen, S. L., Yang, Y., Michelsen, H. A., Najm, H. N.: _Automated reaction kinetics of gas-phase organic species over multiwell potential energy surfaces._ J. Phys. Chem. A, **2023**, 127, 565–588. https://doi.org/10.1021/acs.jpca.2c06558
13. Lockwood, K. S., Ahmed, S. F., Huq, N. A., Stutzman, S. C., Foust, T. D., Labbe, N. J.: _Advances in predictive chemistry enable a multi-scale rational design approach for biofuels with advantaged properties_ Sustainable Energy Fuels, **2022**, 6, 5371-5383. https://doi.org/10.1039/D2SE00773H
14. Takahashi, L., Yoshida, S., Fujima, J., Oikawa, H., Takahashi, K.: _Unveiling the reaction pathways of hydrocarbons via experiments, computations and data science._ Phys. Chem. Chem. Phys., **2022**, 24, 29841-29849. https://pubs.rsc.org/en/content/articlelanding/2022/CP/D2CP04499D
15. Doner, A. C., Zádor, J., Rotavera, B.: _Stereoisomer-dependent unimolecular kinetics of 2,4-dimethyloxetane peroxy radicals._ Faraday Discuss., **2022**, 238, 295-319. https://doi.org/10.1039/D2FD00029F
16. Ramasesha, K., Savee, J. D., Zádor, J., Osborn, D. L.: _A New Pathway for Intersystem Crossing: Unexpected Products in the O(3P) + Cyclopentene Reaction._ J. Phys. Chem. A, **2021**, 125 9785-9801. https://doi.org/10.1021/acs.jpca.1c05817
17. Rogers, C. O, Lockwood, K. S., Nguyen, Q. L. D., Labbe, N. J.: _Diol isomer revealed as a source of methyl ketene from propionic acid unimolecular decomposition._ Int. J. Chem. Kinet., **2021**, 53, 1272–1284. https://doi.org/10.1002/kin.21532
18. Lockwood, K. S., Labbe, N. J.: _Insights on keto-hydroperoxide formation from O2 addition to the beta-tetrahydrofuran radical._ Proceedings of the Combustion Institute, **2021**, 38, 1, 533. https://doi.org/10.1016/j.proci.2020.06.357
19. Sheps, L., Dewyer, A. L., Demireva, M., and Zádor, J.: _Quantitative Detection of Products and Radical Intermediates in Low-Temperature Oxidation of Cyclopentane._ J. Phys. Chem. A **2021**, 125, 20, 4467. https://doi.org/10.1021/acs.jpca.1c02001
20. Zhang, J., Vermeire, F., Van de Vijver, R., Herbinet, O.; Battin-Leclerc, F., Reyniers, M.-F., Van Geem, K. M.: _Detailed experimental and kinetic modeling study of 3-carene pyrolysis._ Int. J. Chem. Kinet., **2020**, 52, 785-795. https://doi.org/10.1002/kin.21400
21. Van de Vijver, R., Zádor, J.: _KinBot: Automated stationary point search on potential energy surfaces._ Computer Physics Communications, **2020**, 248, 106947. https://doi.org/10.1016/j.cpc.2019.106947
22. Joshi, S. P., Seal, P., Pekkanen, T. T., Timonen, R. S., Eskola, A. J.: _Direct Kinetic Measurements and Master Equation Modelling of the Unimolecular Decomposition of Resonantly-Stabilized CH2CHCHC(O)OCH3 Radical and an Upper Limit Determination for CH2CHCHC(O)OCH3+O2 Reaction._ Z. Phys. Chem., **2020**, 234, 1251. https://doi.org/10.1515/zpch-2020-1612

Older Version of KinBot:
1. Van de Vijver, R., Van Geem, K. M., Marin, G. B., Zádor, J.: _Decomposition and isomerization of 1-pentanol radicals and the pyrolysis of 1-pentanol._ Combustion and Flame, **2018,** 196, 500. https://doi.org/10.1016/j.combustflame.2018.05.011
2. Grambow, C. A., Jamal, A., Li, Y.-P., Green, W. H., Zádor, J., Suleimanov, Y. V.: _Unimolecular reaction pathways of a g-ketohydroperoxide from combined application of automated reaction discovery methods._ J. Am. Chem. Soc., 2018, 140, 1035. https://doi.org/10.1021/jacs.7b11009
3. Rotavera, B., Savee, J. D., Antonov, I. O., Caravan, R. L., Sheps, L., Osborn, D. L., Zádor, J., Taatjes, C. A.: _Influence of oxygenation in cyclic hydrocarbons on chain-termination reactions from R + O2: tetrahydropyran and cyclohexane._ Proceedings of the Combustion Institute, **2017,** 36, 597. https://doi.org/10.1016/j.proci.2016.05.020
4. Antonov, I. O., Zádor, J., Rotavera, B., Papajak, E., Osborn, D. L., Taatjes, C. A., Sheps, L.: _Pressure-Dependent Competition among Reaction Pathways from First- and Second-O2 Additions in the Low-Temperature Oxidation of Tetrahydrofuran._ J. Phys. Chem. A, **2016,** 120 6582. https://doi.org/10.1021/acs.jpca.6b05411
5. Antonov, I. O., Kwok, J., Zádor, J., Sheps, L.: OH + 2-butene: A combined experimental and theoretical study in the 300-800 K temperature range. J. Phys. Chem. A, **2015,** 119, 7742. https://doi.org/10.1021/acs.jpca.5b01012
6. Zádor, J., Miller, J.A.: _Adventures on the C3H5O potential energy surface: OH + propyne, OH + allene and related reactions._ Proceedings of the Combustion Institute, **2015,** 35, 181. https://doi.org/10.1016/j.proci.2014.05.103
7. Rotavera, B., Zádor, J., Welz, O., Sheps, L., Scheer, A.M., Savee, J.D., Ali, M.A., Lee, T.S., Simmons, B.A., Osborn, D.L., Violi, A., Taatjes, C.A.: _Photoionization mass spectrometric measurements of initial reaction pathways in low-temperature oxidation of 2,5-dimethylhexane._ J. Phys. Chem. A, **2014,** 44, 10188. https://doi.org/10.1021/jp507811d

## Acknowledgement
This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE organizations, the Office of Science and the National Nuclear Security Administration, responsible for the planning and preparation of a capable exascale ecosystem including software, applications, hardware, advanced system engineering, and early test bed platforms to support the nation's exascale computing imperative. RVdV was also supported by the AITSTME project as part of the Predictive Theory and Modeling component of the Materials Genome Initiative. Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and Engineering Solutions of Sandia, LLC., a wholly owned subsidiary of Honeywell International, Inc., for the U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525. 
