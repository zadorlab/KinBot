# Profiled ethane HPC validation

This external-site test exercises one automatic KinBot reaction workflow and
the portable non-MRCC ANL interfaces:

1. FairChem UMA L1 initial optimization, conformer search, ethane C-C
   homolytic reaction discovery, and product optimization.
2. B2PLYP-D3(BJ)/cc-pVTZ L2 stationary-point refinement through ASE/Sella and
   L2 hindered-rotor scans.
3. The accepted C-C homolysis is treated as a barrierless channel. Gaussian
   prepares the VRC asymptote, two exclusive-node Molpro correction jobs run
   at the configured VRC sampling/high levels, and KinBot writes the runnable
   rotdPy input, runs its reduced Slurm/Molpro sampling, and validates the
   hashed surface-flux and number-of-states outputs. No stationary saddle or
   imaginary frequency is required.
4. A hard gate verifies the accepted reaction, both methyl stoichiometric
   entries, all four requested parent rotor points, the VRC correction record,
   and a completed, hash-verified rotdPy result.
5. A fresh, accepted L2 parent geometry exported from `kinbot.db` into the
   exclusive-node dispatcher.
6. Molpro CCSD(T)/cc-pVTZ ASE/Sella geometry, followed by concurrent Molpro
   harmonic, F12/TZ, F12/QZ, and CCSD(T)/DZ jobs, CFOUR DBOC, and a
   frequency-only Gaussian VPT2 calculation.
7. Native-output hash verification, method-aware parsing, and F12 CBS
   extrapolation.

The final audit status is intentionally
`interface_complete_recipe_incomplete`. MRCC is disabled, and the branch does
not yet have a pinned ANL1-F12 equation or production providers for every
core-valence, scalar-relativistic, spin-orbit, and higher-order component.
Therefore this test must not publish an ANL1-F12 energy or heat of formation.

Run `run.sh PARTITION MAX_CONCURRENT_NODES [FAIRCHEM_MODEL]` after
installing FairChem, obtaining access to the UMA model, loading the licensed
QC programs, initializing and installing the private `external/ROTD_py`
submodule, and activating the KinBot environment. `FAIRCHEM_MODEL` may be a
registered name or, preferably for offline compute nodes, the absolute path
to a previously downloaded checkpoint. It can also be supplied through
`KINBOT_FAIRCHEM_MODEL`. The script is restartable: KinBot reuses its database
and the ANL dispatcher reuses its immutable workflow state.

The methyl product still receives its ordinary minimum Hessian because MESS
needs fragment partition functions. The barrierless reaction itself never
receives a transition-state Hessian. The generated input contains one
dividing-surface distance, small temperature/energy/angular grids, at most
eight requested samples, one Molpro sampling process, and one queued sampling
job. Its sampling level is `caspt2(2,2)/vdz`; the asymptotic correction also
evaluates the configured `caspt2(2,2)/avtz` high level. This is an execution
and interface smoke test, not a production-converged VRC calculation.
