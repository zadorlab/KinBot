[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/zadorlab/KinBot)

# KinBot: Automated Reaction Kinetics of Gas-Phase Organic Species over Multiwell Potential Energy Surfaces

<p>
    <img src="https://raw.githubusercontent.com/zadorlab/KinBot/master/graphics/kinbot_logo_V2.png" width="220" height="240" />
</p>

## Description
This repository contains the KinBot code version 2.4.1,
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

KinBot can be installed in three different ways: from the PyPI index (`pip install`), from the conda-forge repo (`conda install`), or by cloning this github repo and installing it locally.

### PyPI

    pip install kinbot

> **Note**
>  KinBot only works with Python >= 3.11.

RDKit >= 2026.3.5 is a required dependency. KinBot records its version and
stereochemistry settings in the log and `.kinbot_run.json`. Restart a calculation
with the same RDKit version and settings. Calculations from older KinBot versions
are incompatible with this release; start them in new directories.

To use the FAIRChem/UMA machine-learning interatomic potential backend
(`qc = 'fc'`), install the optional `fc` dependencies (this also pulls in
`fairchem-core`, `torch`, and a compatible `numpy`/`scipy`):

    pip install kinbot[fc]

The optional `plot` extra (`matplotlib`) enables the plotting
utilities:

    pip install kinbot[plot]

### conda-forge

    conda install -c conda-forge kinbot

### From Github

If you want to have the very last version of KinBot without waiting for a 
release or you want to modify KinBot according to your needs you can clone the project 
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

SMILES inputs use RDKit with a fixed seed to generate the starting geometry.
The same SMILES and RDKit version give the same starting stereoisomer and
calculation name. Explicit tetrahedral and double-bond SMILES assignments are
retained. Unspecified stereochemistry selects one reproducible starting
structure; it does not request all stereoisomers. Supply coordinates to select
a particular cumulene orientation. Coordinate inputs are not changed.

### Product complexes in MESS

`me_skip_vdW` controls the final MESS network in both a single-well search and
a PES calculation. The default, `0`, includes accepted product complexes as
wells, with barrierless exits to the separated products. Set `me_skip_vdW: 1`
to omit these wells and exits and connect each inner transition state directly
to the separated products.

This option does not change the search or optimization of complexes;
`vdW_detection` still determines whether a complex is bound. The inner
transition state's Eckart tunneling depths keep the complex energy reference
with either setting. Both writers use one lowest-energy complex per set of
product stereoisomers, as in the existing PES model. For MC-TST, the selected
complex structure defines the tunneling reference, while its lowest retained
conformer defines the MESS well ground energy.

Omitting a complex removes its collisional stabilization and its explicit
capture/redissociation competition. The two settings need not give the same
rates. `correct_submerged` only raises barriers to connected bound-well ground
energies; it does not raise an inner TS merely because it lies below the
separated fragments. This option does not add support for a MESS barrier
between two sets of separated reactants and products.

### Preliminary conformer searches (L0)

For wells, `l0_conformer_search: 1` enables a preliminary, inexpensive
conformer search before L1 when `conformer_search: 1` is enabled. It remains
disabled by default and does not enable preliminary transition-state conformer
searches. The existing QC backend, command, calculation filenames and defaults
remain unchanged.

| Input | Default | Previous name (still accepted) |
| --- | --- | --- |
| `l0_conformer_search` | `0` | `semi_emp_conformer_search` |
| `l0_method` | `am1` | `semi_emp_method` |
| `l0_max_dihed` | `5` | `max_dihed_semi_emp` |
| `l0_random_conf` | `500` | `random_conf_semi_emp` |
| `l0_conformer_threshold` | `5` kcal/mol | `semi_emp_confomer_threshold` |

`l0_max_dihed` sets the maximum number of dihedral angles for exhaustive
sampling. Beyond that number, `l0_random_conf` sets the number of randomly
sampled conformers. `l0_conformer_threshold` is the energy window above the
lowest L0 conformer for retaining seeds for L1.

The previous `L0_conformer_search` spelling is also accepted. Old names emit a
deprecation warning. Conflicting values for an alias and its replacement are
rejected, including the historically misspelled `semi_emp_confomer_threshold`.

### Independent L0, L1 and L2 methods

`l0_profile`, `l1_profile` and `l2_profile` select the calculator separately
for each level. These reuse the profile names and fields from Hrishikesh Ram's
`composite` branch. For example, this input uses a small UMA model for preliminary
exploration, a larger UMA model for discovery, and Gaussian for refinement:

```json
{
  "l0_profile": {
    "calculator": "fairchem", "model_path": "uma-s-1p1",
    "task_name": "omol", "device": "cuda"
  },
  "l1_profile": {
    "calculator": "fairchem", "model_path": "uma-m-1p1",
    "task_name": "omol", "device": "cuda"
  },
  "l2_profile": {
    "calculator": "gaussian", "method": "M062X", "basis": "6-31G(d)",
    "command": "g16", "optimizer": "native",
    "calculator_kwargs": {"integral": "UltraFine"}
  },
  "conformer_search": 1, "l0_conformer_search": 1, "high_level": 1
}
```

This is a configuration fragment to combine with the molecule's usual input.
`model_path` accepts a local checkpoint or a supported pretrained model name.
For Gaussian AM1 exploration use `l0_profile` with `calculator: "gaussian"`,
`method: "AM1"`, `basis: ""` and the desired Gaussian `command`.

L0 covers preliminary conformers and ring preparation.
L1 covers ordinary well/conformer optimization, reaction discovery and its IRC.
L2 is enabled by `high_level: 1`; it refines stationary points. Hindered rotors
use the selected stationary point's level. The preliminary conformer switch
remains restricted to wells.

Explicit profile fields override the corresponding legacy settings for the
selected optimization level. Adding `l1: "gauss"` alone does not change the
reaction search: Gaussian constrained TS steps still use AM1; Q-Chem steps
still use B3LYP/STO-3G. The final TS optimization uses L1. Scans and barrierless
saddles retain `scan_method`/`scan_basis` and `barrierless_saddle_*` settings.
These options also apply with profiles. Explicit `calculator_kwargs` override
generated calculator keywords. A FairChem L1 profile uses its UMA model
throughout the search. L0 does not control the constrained TS steps.
`l0`, `l1` and `l2` also accept calculator aliases such as `gauss`, `fc` and `uma`.
Ordinary inputs without profiles retain their previous settings, including the
Gaussian command used by preliminary searches. Flat `l0_qc`,
`l0_qc_command`, `l0_basis` and `l0_calc_kwargs` inputs remain accepted;
`l0_method` is a fallback beneath `l0_profile.method`.

`optimizer` can be `native` for a backend's own optimizer or `sella`; `bfgs`
selects BFGS for wells, while saddles and constrained searches continue to use
Sella. FairChem uses its existing Sella/BFGS workers. `frequency_mode: "auto"`
retains the backend's usual frequency calculation, and `ase_forces` uses finite
differences with Sella/BFGS. Gaussian also supports `native_hessian` after Sella;
Q-Chem profiles require `optimizer: "native"` and support native Hessians.
Sella/BFGS with Q-Chem is rejected because those workers do not yet request
force-only calculations. Ordinary inputs without profiles are unchanged.
NWChem profiles are rejected because their calculation workers are not supported.
Q-Chem can optimize ordinary conformers, but it cannot prepare ring conformers.
For calculations that need ring preparation, select Gaussian or FairChem at L0;
Q-Chem can still supply L1 or L2. An attempted Q-Chem ring preparation stops with
this explanation before KinBot writes or submits the ring job.

MESS and PES describe the resolved settings. An explicit profile `label` is
retained. ANL presets and composite-recipe selection are not included; specify
the methods in the profiles themselves. An empty `theory_preset` field remains
accepted for previously written inputs, but a named preset is rejected.
`frequency_mode` selects how harmonic frequencies are calculated; it does not
request vibrational perturbation theory or anharmonic corrections.

Each generated worker uses its assigned model, command and calculator keywords.
The profile route file, `.kinbot_theory_jobs.json`, records the methods before
the first job and records each job's level for reading results and restarting.
KinBot refuses to add profiles to existing calculations without a matching
theory record. It also refuses to restart a profiled calculation after profiles
are removed. The existing `.kinbot_run.json` also remembers profile use, so
removing a custom `profiled_route_file` setting cannot hide it. Calculations
started before this flag was introduced receive it when restarted with their
original profiles and route-file setting; keep those settings for that restart.
PES checks the same theory records before replacing species inputs or reading
results, including with `no-kinbot`. These checks compare profile settings;
they do not add method tracking to older calculations made without profiles.
The comparison also includes ORCA's `orcasimpleinput` and `orcablocks`, the
NNPES `nn_model`, Gaussian integration-grid and electronic-guess settings,
and the existing scan-method settings used by the selected calculator.
Records from earlier versions that omitted these settings are rejected even
with unchanged input, because they cannot establish how the saved energies
were calculated. Keep those results unchanged and use a fresh directory for
calculations with the updated code.

Hindered-rotor scans use L1 when `high_level` is 0 and L2 when it is 1. Changing
`high_level` after rotor jobs have been recorded requires a fresh directory;
KinBot refuses to reuse scans from the other level, including during PES
postprocessing. Changing this option before any rotor jobs were recorded is
still allowed.
Keep the scientific settings and route-file location fixed when resuming a run;
use a fresh directory when changing methods, models, optimizers or frequency
treatment. The executable `command` and CPU/GPU `device` can change on restart.
They are recorded separately and do not identify the calculation method.
Earlier profile records used different TS-search defaults and lack a complete
record of these choices; start a fresh directory for those runs.
CPU, memory and queue limits
remain the existing global settings. L3 retains its existing settings.


## Documentation
See the [wiki](https://github.com/zadorlab/KinBot/wiki) for keywords, and our [tutorial](https://hackmd.io/@jzador/ry1DSsEyyx#/) for a more hands-on introduction to the code.

## List of files in this project
See [list](https://github.com/zadorlab/KinBot/wiki/KinBot-file-structure).

## Authors
* Judit Zádor (jzador@sandia.gov)
* Ruben Van de Vijver 
* Amanda Dewyer
* Carles Martí
* Clément Soulié (csoulie@sandia.gov)
* Luka Dockx

## Papers using KinBot
1.	Kendall, D. T., Zádor J.: _Benchmarking the UMA foundation interatomic potential for gas-phase chemical kinetics_ J. Phys. Chem. A, **2026**, 130, 4450-4468. https://pubs.acs.org/doi/10.1021/acs.jpca.6c01748
2.	Kjaersgaard, A., Zádor J., Zwier, T. S., Shiels, O.J., Sheps, L.: _Measurement of the CH3SCH2O2 -> CH2SCH2OOH rate coefficient and mass spectrometric characterization of hydroperoxymethyl thioformate (HPMTF) at T = 400 – 550 K._ J. Phys. Chem. A, **2025**, 129, 10452-10462. https://doi.org/10.1021/acs.jpca.5c06279
3. Shiels, O. J., Kjaersgaard, A., Zádor, J., Sheps, L.: _Low-temperature autooxidation of cyclopentanone and 3-pentanone: The critical role of competing radical chain-branching and chain inhibiting pathways._ Phys. Chem. Chem. Phys., **2025**, 27, 23813-23825. https://doi.org/10.1039/D5CP03403E
4. Osborn, D. L., Samanta, B. R., Soulié, C., Reisler, H., Zádor, J.: _Chemistry of sugar formation in the gas phase: Following the activated aldehyde_. J. Am. Chem. Soc., **2025** 147, 32468–32479. https://doi.org/10.1021/jacs.5c05357
5. Hansen, N., Gaiser, N., Bierkandt, T., Oßwald, P., Köhler, M., Zádor, J., Hemberger, P.: _Identification of dihydropentalenes as products of the molecular-weight growth reaction of cyclopentadienyl plus propargyl_. J. Phys. Chem. A, **2025**, 129, 1714-1725. https://doi.org/10.1021/acs.jpca.4c06549
6. Wang, Z.-M., Wang, D., Mansy, A. E., Tian, Z.-Y.: _Pyrolysis and kinetic modeling investigation of 1-methoxy-2-propanol_. Combust. Flame, **2025**, 271, 113804. https://doi.org/10.1016/j.combustflame.2024.113804
7. Almeida, T. G., Martí, C., Kurtén, T., Zádor, J., Johansen, S. L.: _Theoretical analysis of the OH-Initiated atmospheric oxidation reactions of imidazole_. Phys. Chem. Chem. Phys., **2024**, 26, 23570-23587. https://doi.org/10.1039/D4CP02103G
8. Yuan, E. C.-Y., Kumar, A., Guan, X., Hermes, E. D., Rosen, A. S., Zádor, J., Head-Gordon, T., Blau, S. M.: _Analytical ab initio Hessian from a Deep Learning Potential for Transition State Optimization_. Nat. Comm., **2024**, 15, 8865. https://doi.org/10.1038/s41467-024-52481-5
9. Doner, A. C., Zádor, J., Rotavera, B.: _Stereoisomer-dependent rate coefficients and reaction mechanisms of 2-ethyloxetanylperoxy radicals_. Proc. Combust. Inst., **2024**, 40, 105578. https://doi.org/10.1016/j.proci.2024.105578
10. Hansen, N. A, Price, T. D., Filardi, L. R., Gurses, S. M., Zhou, W., Hansen, N., Osborn, D. L. Zádor, J., Kronawitter, C. X.: _The photoionization of methoxymethanol: Fingerprinting a reactive C2 oxygenate in a complex reactive mixture_. J. Chem. Phys., **2024**, 160, 124306. https://doi.org/10.1063/5.0197827
11. Martí, C., Devereux, C., Najm, H. N., Zádor, J.: _Evaluation of rate coefficients in the gas-phase using a machine learned potential_. J. Phys. Chem. A, **2024**, 128, 1958–1971. https://doi.org/10.1021/acs.jpca.3c07872
12. Lang, J., Foley, C. D., Thawoos, S., Behzadfar, A., Liu, Y., Zádor, J., Suits, A. G.: _Reaction dynamics of S(3P) with 1,3-butadiene and isoprene: Crossed beam scattering, low temperature flow experiments, and high-level electronic structure calculations_. Farad. Discuss., **2024**, 251, 550-572. https://doi.org/10.1039/D4FD00009A
13. Wang, D., Tian, Z.-Y., Zheng, Z.-H., Li, W., Wu, L.-N., Kuang, J.-J., Yang, J.-Z.: _Experimental and modeling study of the n, n-dimethylformamide pyrolysis at atmospheric pressure_. Combust. Flame, **2024**, 260, 113240. https://doi.org/10.1016/j.combustflame.2023.113240
14. Wang, D. Dai, W.-K., Tian, Z.-Y., Liu, Y.-W., Wang, Z.-M., Wu, L.-N., Kuang, J.-J., Wang ,Q.-P., Xu, Q., Wang, Z.-D.: _Theoretical and experimental study on the pyrolysis of N-methylpyrrolidone_. J Anal Appl Pyrolysis, **2024**, 183, 106751. https://doi.org/10.1016/j.jaap.2024.106751
15. Doner, A. C., Zádor, J., Rotavera, B.: _Unimolecular reactions of 2,4-dimethyloxetanyl radicals._ J. Phys. Chem A, **2023**, 127, 2591–2600 https://doi.org/10.1021/acs.jpca.2c08290
16. Li, H., Lang, J., Foley, C. D., Zádor, J., Suits, A. G.: _Sulfur (3P) reaction with conjugated dienes gives cyclization to thiophenes under single collision conditions._ J. Phys. Chem. Letters, **2023**, 14, 7611–7617. https://doi.org/10.1021/acs.jpclett.3c01953
17. Martí, C., Michelsen, H. A., Najm, H. N., Zádor, J.: _Comprehensive kinetics on the C7H7 potential energy surface under combustion conditions._ J. Phys. Chem. A, **2023**, 127, 1941–1959. https://pubs.acs.org/doi/full/10.1021/acs.jpca.2c08035
18. Zádor, J, Martí, C., Van de Vijver, R., Johansen, S. L., Yang, Y., Michelsen, H. A., Najm, H. N.: _Automated reaction kinetics of gas-phase organic species over multiwell potential energy surfaces._ J. Phys. Chem. A, **2023**, 127, 565–588. https://doi.org/10.1021/acs.jpca.2c06558
19. Lockwood, K. S., Ahmed, S. F., Huq, N. A., Stutzman, S. C., Foust, T. D., Labbe, N. J.: _Advances in predictive chemistry enable a multi-scale rational design approach for biofuels with advantaged properties_ Sustainable Energy Fuels, **2022**, 6, 5371-5383. https://doi.org/10.1039/D2SE00773H
20. Takahashi, L., Yoshida, S., Fujima, J., Oikawa, H., Takahashi, K.: _Unveiling the reaction pathways of hydrocarbons via experiments, computations and data science._ Phys. Chem. Chem. Phys., **2022**, 24, 29841-29849. https://pubs.rsc.org/en/content/articlelanding/2022/CP/D2CP04499D
21. Doner, A. C., Zádor, J., Rotavera, B.: _Stereoisomer-dependent unimolecular kinetics of 2,4-dimethyloxetane peroxy radicals._ Faraday Discuss., **2022**, 238, 295-319. https://doi.org/10.1039/D2FD00029F
22. Ramasesha, K., Savee, J. D., Zádor, J., Osborn, D. L.: _A New Pathway for Intersystem Crossing: Unexpected Products in the O(3P) + Cyclopentene Reaction._ J. Phys. Chem. A, **2021**, 125 9785-9801. https://doi.org/10.1021/acs.jpca.1c05817
23. Rogers, C. O, Lockwood, K. S., Nguyen, Q. L. D., Labbe, N. J.: _Diol isomer revealed as a source of methyl ketene from propionic acid unimolecular decomposition._ Int. J. Chem. Kinet., **2021**, 53, 1272–1284. https://doi.org/10.1002/kin.21532
24. Lockwood, K. S., Labbe, N. J.: _Insights on keto-hydroperoxide formation from O2 addition to the beta-tetrahydrofuran radical._ Proceedings of the Combustion Institute, **2021**, 38, 1, 533. https://doi.org/10.1016/j.proci.2020.06.357
25. Sheps, L., Dewyer, A. L., Demireva, M., and Zádor, J.: _Quantitative Detection of Products and Radical Intermediates in Low-Temperature Oxidation of Cyclopentane._ J. Phys. Chem. A **2021**, 125, 20, 4467. https://doi.org/10.1021/acs.jpca.1c02001
26. Zhang, J., Vermeire, F., Van de Vijver, R., Herbinet, O.; Battin-Leclerc, F., Reyniers, M.-F., Van Geem, K. M.: _Detailed experimental and kinetic modeling study of 3-carene pyrolysis._ Int. J. Chem. Kinet., **2020**, 52, 785-795. https://doi.org/10.1002/kin.21400
27. Van de Vijver, R., Zádor, J.: _KinBot: Automated stationary point search on potential energy surfaces._ Computer Physics Communications, **2020**, 248, 106947. https://doi.org/10.1016/j.cpc.2019.106947
28. Joshi, S. P., Seal, P., Pekkanen, T. T., Timonen, R. S., Eskola, A. J.: _Direct Kinetic Measurements and Master Equation Modelling of the Unimolecular Decomposition of Resonantly-Stabilized CH2CHCHC(O)OCH3 Radical and an Upper Limit Determination for CH2CHCHC(O)OCH3+O2 Reaction._ Z. Phys. Chem., **2020**, 234, 1251. https://doi.org/10.1515/zpch-2020-1612

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
