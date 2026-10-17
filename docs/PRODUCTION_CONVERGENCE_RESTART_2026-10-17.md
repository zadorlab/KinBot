# Fresh two-PES restart with practical convergence thresholds

The r6 launch was interrupted by cluster trouble. The settings below are
unchanged; the current [recovery command](PRODUCTION_HANDOFF_RECOVERY_2026-10-17.md)
archives r6 and starts fresh r7 runs with repeat-safe handoff checks.

The user superseded the saved-results resume with a full clean restart of both
examples on 17 October. No calculated results are shared between the two PES or
copied from the retired attempt. The reviewed orbital WFU files remain input
data used to initialize the selected active spaces; they are not cached energies.

## Numerical contract

| Stage | Requested setting | Native meaning |
| --- | --- | --- |
| Gaussian L1/L2, including CBH references and L2 VPT2 | SCF 1e-8 | `Conver=8`: RMS density change 1e-8 and maximum change 1e-6; not an energy-error bound. |
| ANL Molpro SCF and VRC SCF | 1e-8 | Global energy threshold; SCF and local correlated thresholds are rendered separately. |
| ANL CC and high-spin VRC CC anchor | 1e-7 | Molpro local CC energy threshold 1e-7 Eh; squared-amplitude threshold no looser than its documented 1e-10 default. MRCC `sctol=8,cctol=7`. The CFOUR HF-DBOC task explicitly receives `SCF_CONV=8`; retained legacy CFOUR CC templates use `CC_CONV=7,LINEQ_CONV=7`. These backends have different convergence metrics. |
| CASPT2/MRCI | Energy 1e-7 Eh; coefficient norm 1e-5 | Literal `THRDEN=1.d-7,THRVAR=1.d-10`; THRVAR is the squared norm. The energy threshold is relaxed from 1e-8; the coefficient threshold is retained, one order tighter than Molpro's documented global COEFF default of 1e-4. |
| CASSCF orbitals | Energy 1e-8 Eh; orbital gradient 1e-5 | Existing reviewed MULTI settings retained. These are separate from the CASPT2/MRCI iterative solution. |
| Final L1/L2/L3 geometries | Maximum force norm 5e-4 eV/angstrom | About 9.72e-6 Eh/bohr. Sella uses a force norm, not Gaussian Berny's complete Tight force/displacement criteria. |
| Molpro constrained VRC relaxation | Maximum projected force norm 5e-4 eV/angstrom | Exposed as `vrc_tst_scan_fmax`, replacing a fixed 0.03 value. The reactive-distance constraint remains enforced. |

The single JSON in each example explicitly specifies these production values.
`anl_scf_convergence` reaches both levels of the DFT preset; a deliberate
per-level `calculator_kwargs.scf` override remains supported. Reference inputs
inherit the parent's parameters. The ANL graph carries the SCF/CC/geometry
settings into every automatic task. The six basis-specific VRC models specify
their correlated energy threshold explicitly. The shared CASPT2/MRCI renderer
also defaults to 1e-7 energy and 1e-5 coefficient norm for new inputs without
overrides. Historical reconnaissance inputs/outputs retain their actual settings.

These are reasonable numerical tolerances for subchemical thermochemistry,
not a certification of total chemical accuracy or a rate-error bound. A 1e-7 Eh
energy change is about 0.000063 kcal/mol, but a stopping threshold is not a bound
on the remaining energy error. Soft modes, finite-difference derivatives,
model/basis errors and electronic-state choices remain relevant. Native
frequency validation, bounded minimum recovery and unresolved-SCF failure
rejection remain active. The 1% RotdPy target is sampling precision.

Sources: [Gaussian SCF documentation, vendor mirror](https://www.conflex.co.jp/gaussian_support/scf.php),
[Gaussian optimization documentation](https://www.conflex.co.jp/gaussian_support/opt.php),
[Molpro CCSD thresholds](https://www.molpro.net/manual/doku.php?id=the_closed_shell_ccsd_program),
[Molpro MRCI options](https://www.molpro.net/manual/doku.php?id=the_mrci_program),
[Molpro global thresholds](https://www.molpro.net/manual/doku.php?id=program_control#global_thresholds_gthresh),
[Molpro MULTI convergence](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#convergence_thresholds).

## Latest failures and repairs

The CH3SS SH reference L2 job first failed XQC, then its CDIIS retry stalled at
RMSDP=8.95e-10 and MaxDP=1.75e-8. Gaussian printed a nonconvergence warning
and then normal termination. The reader correctly rejected it at Conver10.
Those residuals are below the newly requested Conver8 limits; this supports
the relaxation but is not proof that every fresh SCF will converge. The actual
source is the user's native report for
`r5_dft/ch3ss/791471490990000000002/anl_automatic/reference_l2/state_20691a101e49c36b/330330000000000000002_well_high.log`.
The failed CDIIS calculation rules out blaming XQC alone. The working atomic
symmetry treatment remains in place. No untested molecular symmetry restriction
or change of spin/reference type was added.

The prior resume helper incorrectly archived the entire VRC directory,
including successful fragment checkpoint/output files, while keeping their
database rows. This caused waits on old job IDs and missing HOMO-pivot
checkpoints. The helper now preserves fragment handoffs, with checked same-run
restoration for artifacts moved by that old helper. The chosen **fresh launch
does not call this restore/resume path**.

The composite thread now publishes its complete/failed branch state immediately,
even while VRC is still active. This prevents the dashboard from reporting a
known failed CBH branch as running until VRC returns.

The threshold audit also found that a local `thresh,energy` substitution matched
inside `gthresh,energy`, overwriting the global SCF value with the CC value.
It now excludes the global card. L3 gradient and VRC CC-anchor renderers now
have explicit local CC thresholds instead of inheriting the SCF energy limit.
Molpro squared-amplitude controls are kept distinct from energy controls.
The HF-DBOC template previously omitted an explicit SCF limit; it now receives
the same selected SCF convergence as the other automatic ANL tasks.

## Clean restart

After pulling `composite`, `examples/anl/tools/restart_on_blodgett.sh` prepares
and preflights the fresh parent
`/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261017_r6_dft_scf8`.
Only after successful preflight does it stop the two r5 PES controllers/jobs,
verify they are inactive, and rename their complete common parent to a dated
`.stopped.*` sibling. Existing external repair archives remain intact.
It then starts both fresh drivers detached on the login host. The old attempt's
logs, databases, failed markers and calculated models are absent from the new
run. It neither deletes shared scratch nor touches unrelated jobs.

Each PES retains its 10-allocation cap, all eleven bad-node exclusions, the
DFT preset, selected orbital indices/roots, IPEA 0.25, SHIFT 0.1, aDZ/aTZ 3.7
gap extrapolation, experimental electronic levels/SO corrections, ordinary single-bottleneck VRC
and 1% flux-error target.
MRCI is not enabled by changing its general numerical defaults.

The parent-directory monitor still covers both runs. End-to-end ANL,
RotdPy and MESS completion remains pending native execution. Local tests
exercise actual profile resolution, reference inheritance, all ANL backends'
rendered controls, CAS/MRCI defaults and overrides, archive integrity and branch
reporting. They do not execute licensed Gaussian or Molpro locally.

Local verification: the full suite passed **1930 tests and 424 subtests**
(zero failures/skips, two warnings), including seven native
RotdPy/MESS/MESSPF checks. After the later barrier/temperature input edits,
**53 tests and 11 subtests** passed for preparation, input propagation,
reference/profile routing and VRC/PES integration. These counts are separate,
not additive. [Commands, logs and source hashes](../examples/anl/tools/convergence_restart_validation.json)
record the tested scope. Shell syntax and `git diff --check` also passed.

## Barrier screening and temperature range

The later user instruction replaces the 1000 kcal/mol search cutoff with
`barrier_threshold_L2=400` and `barrier_threshold_add=10`. KinBot therefore
screens L1 candidates at 410 kcal/mol and applies the 400 kcal/mol L2 test to
converged non-homolytic TSs. These are ZPE-corrected energies, not free-energy
barriers. The existing homolysis endpoint allowance remains in effect.
PES child thresholds account for the product-well energy so that the energy
ceiling is referenced to the initial well. No reaction family is disabled.

This is a conservative example choice, not a new universal cutoff algorithm:
at 5000 K, RT=9.9360 kcal/mol and exp(-400/RT)=3.284e-18. Even this number is
not a rigorous rate bound: the [transition-state-theory expression](https://goldbook.iupac.org/terms/view/T06470)
also contains entropy and a frequency prefactor. At zero activation entropy,
the Eyring estimate for a 400 kcal/mol activation enthalpy is 3.42e-4 s^-1;
using 200 kcal/mol would instead give about 1.89e5 s^-1. That is why the cutoff
was not reduced to a common low-temperature value. Failed TS searches can
still occur before a barrier can be evaluated. Final fluxes determine importance.

Each single input explicitly sets MESS `TemperatureList` to 300–5000 K at
100 K intervals and `rotdpy_temperature_grid=[300,100,1,48]` to match it.
The explicit RotdPy energy grid `[0,10,1.05,190]` retains its default values
and ends near 4018 kcal/mol. `ExcessEnergyOverTemperature=30` is unchanged;
`ModelEnergyLimit=800` leaves room above 400 kcal/mol plus 30 RT at 5000 K
(about 698 kcal/mol). The [MESS manual, section 5.3](https://tcg.cse.anl.gov/papr/codes/mess/messmanual_3_23_16.pdf)
defines these as distinct matrix/model limits. Generic KinBot barrier,
temperature and MESS energy defaults are unchanged. No second settings file
or molecule-specific threshold code is introduced.

The numerical range extension does not validate a complete high-temperature
electronic model. In particular, [NIST](https://physics.nist.gov/PhysRefData/Handbook/Tables/sulfurtable5_a.htm)
places sulfur's separate 1D term at 9238.609 cm^-1. Its equilibrium electronic
weight at 5000 K is about 4% of the total ground-3P-term weight. These inputs
describe the reviewed ground-term reaction surfaces and fine structure;
they do not add excited-term reaction surfaces or change the SO centroid by
mixing different terms. Consequently the 1% Monte Carlo target cannot be
interpreted as 1% physical accuracy over the newly extended temperature range.
