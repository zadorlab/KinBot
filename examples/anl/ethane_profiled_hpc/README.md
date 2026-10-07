# Profiled ethane HPC validation

This external-site test exercises one automatic KinBot reaction workflow and
the portable ANL interfaces:

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
6. Molpro CCSD(T)/cc-pVTZ ASE/Sella geometry, followed by one globally
   throttled fan-out containing Molpro harmonic, F12, core-valence, and
   relativistic jobs; direct-MRCC higher-order jobs; CFOUR DBOC; and a
   frequency-only Gaussian VPT2 calculation on its matching L2 geometry.
7. Native-output hash verification, method-aware parsing, correction
   assembly, and one-run ANL audit.

Fresh runs create `anl_composite`. The affordable default is the complete
ANL0/ANL0-F12 interface set. Set `KINBOT_ANL_SCOPE=anl1` to include the
ANL1-only CCSDT(Q)/cc-pVTZ and CCSDTQ(P)/cc-pVDZ jobs; these may require days.
An existing `anl_interface` directory is treated as immutable legacy state
and continues through its original interface-only audit.

Run `run.sh PARTITION MAX_CONCURRENT_NODES [FAIRCHEM_MODEL]` after
installing FairChem, obtaining access to the UMA model, loading the licensed
QC programs, initializing and installing the private `external/ROTD_py`
submodule, and activating the KinBot environment. `FAIRCHEM_MODEL` may be a
registered name or, preferably for offline compute nodes, the absolute path
to a previously downloaded checkpoint. It can also be supplied through
`KINBOT_FAIRCHEM_MODEL`. `KINBOT_ANL_SCOPE` accepts `anl0` or `anl1`. The
script is restartable: KinBot reuses its database and the ANL dispatcher
reuses its immutable workflow state.

The methyl product still receives its ordinary minimum Hessian because MESS
needs fragment partition functions. The barrierless reaction itself never
receives a transition-state Hessian. The generated input contains one
dividing-surface distance, small temperature/energy/angular grids, at most
eight requested samples, one Molpro sampling process, and one queued sampling
job. Its sampling level is `caspt2(2,2)/vdz`; the asymptotic correction also
evaluates the configured `caspt2(2,2)/avtz` high level. This is an execution
and interface smoke test, not a production-converged VRC calculation.

## Production methyl + methyl VRC continuation

After the reduced end-to-end test has completed, run
`run_vrc_production.sh PARTITION 8 [FAIRCHEM_MODEL]` against the same run
directory. It preserves the reduced ROTD_py directory, performs a real
multipoint correction scan, evaluates the ROTD samples at
CASPT2(2e,2o)/cc-pVDZ, evaluates the trusted correction at
MRCI+Q(2e,2o)/cc-pVTZ using Molpro's Davidson-corrected `ENERGD`, and samples
24 dividing-surface constructions from 2.3 through 4.6 angstrom with the
0.85 dynamical correction and the normal temperature, energy,
angular-momentum, and Monte Carlo convergence grids. At most eight exclusive
sampling jobs or eight exclusive correction jobs run concurrently.

Every correction output is tied to its generated input hash, so the earlier
CASPT2/aug-cc-pVTZ asymptote cannot be reused as the new MRCI+Q correction.
The final `rotdpy_production_gate.json` is written only after all surfaces
report converged Monte Carlo samples and their numeric ROTD/MESS outputs pass
hash verification. Set `rotdpy_mess_mode` to `production` in the later MESS
run to prevent the reduced smoke result from entering a kinetics model.
`monitor_vrc_production.sh` reports the active versioned correction graph,
the correction levels and point count, Monte Carlo products, production gate,
queue, and recent KinBot messages without changing the run.
