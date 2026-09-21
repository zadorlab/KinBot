# Profiled ethane HPC validation

This external-site test exercises one automatic KinBot reaction workflow and
the portable non-MRCC ANL interfaces:

1. FairChem UMA L1 initial optimization, conformer search, ethane C-C
   homolytic reaction discovery, and product optimization.
2. B2PLYP-D3(BJ)/cc-pVTZ L2 stationary-point refinement through ASE/Sella and
   L2 hindered-rotor scans.
3. A fresh, accepted L2 parent geometry exported from `kinbot.db` into the
   exclusive-node dispatcher.
4. Molpro CCSD(T)/cc-pVTZ ASE/Sella geometry, followed by concurrent Molpro
   harmonic, F12/TZ, F12/QZ, and CCSD(T)/DZ jobs, CFOUR DBOC, and a
   frequency-only Gaussian VPT2 calculation.
5. Native-output hash verification, method-aware parsing, and F12 CBS
   extrapolation.

The final audit status is intentionally
`interface_complete_recipe_incomplete`. MRCC is disabled, and the branch does
not yet have a pinned ANL1-F12 equation or production providers for every
core-valence, scalar-relativistic, spin-orbit, and higher-order component.
Therefore this test must not publish an ANL1-F12 energy or heat of formation.

Run `run.sh PARTITION MAX_CONCURRENT_L3_NODES` after installing FairChem,
obtaining access to the UMA model, loading the licensed QC programs, and
activating the KinBot environment. The script is restartable: KinBot reuses
its database and the ANL dispatcher reuses its immutable workflow state.
