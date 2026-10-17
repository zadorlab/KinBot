# DFT preset restart and automatic recovery of rejected minima

## Native evidence and decision

The CH3S L2 worker in
`/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016_r2/ch3s/470831050000000000002`
converged to small forces, but its independent Gaussian Hessian contained
**−587.3340 cm⁻¹**. Rejecting that point as a minimum was correct. The worker
previously hid the reason in an error-only result; the new worker records the
actual spectrum rejection or exception in stderr, the Sella log and its result.

The Sella and native frequency energies were −437.92436235217 and
−437.92436243561 Eh, respectively. The user then ran an isolated directional
check at unchanged B2PLYP-D3(BJ)/cc-pVTZ, with independent checkpoint copies.
Its native report is `/home/hram/kinbot-runs/l2-curvature.eJZhr1/report.txt`;
the user reported an empty `error.txt`. The following are the supplied report
values, not an independently downloaded native archive:

| Displacement (Å, full Cartesian norm) | Energy (Eh) | Projected force (eV/Å) |
|---:|---:|---:|
| 0 | −437.92436243565 | 0.00039646464957125704 |
| −0.01 | −437.92436520235 | −0.015047182302865945 |
| +0.01 | −437.92436547877 | 0.01572179531092174 |
| −0.02 | −437.92437292005 | −0.025860918893797505 |
| +0.02 | −437.92437339153 | 0.026221427610951958 |

The energy/force curvatures are −1.58093/−1.53845 eV/Å² at 0.01 Å and
−1.45855/−1.30206 eV/Å² at 0.02 Å. Both signs lower the energy; both ways of
estimating curvature are negative. This establishes a downhill direction on
the calculated surface, rather than a frequency-parser artifact. It does not
establish the ultimate electronic cause of the instability or prove a global
minimum can be reached with a bounded local retry.

The reference ⟨S²⟩ was 0.7569 (pure doublet: 0.75). That observation and the
curvature data do not justify changing all calculations to restricted-open-shell
B2PLYP. Menon and Radom found that opposing HF and MP2 spin-contamination
effects can partly cancel in unrestricted double hybrids; RO is not a universal
repair. Their benchmark is supporting context, not validation of this sulfur
surface. [J. Phys. Chem. A 112, 13225–13230 (2008), DOI 10.1021/jp803064k](https://pubmed.ncbi.nlm.nih.gov/18759419/).
The L2 recovery retains the electronic method, multiplicity, dispersion and
basis. CAS, root and VRC settings are unchanged. The later user-approved
**L1** change and full restart are separate decisions, described below.

## Second failure and selected preset

CH3SS passed its initial L2 calculation but rejected the CH3S product
`470831050000000000002_well` at L1, inside its `791471490990000000002`
worker directory. The user supplied the Sella log and error pickle: the UMA
spectrum began at **−53.0623707468 cm⁻¹**, beyond the existing −50 cm⁻¹
criterion. The remaining eight modes were positive. The error was a real
frequency rejection, not a Gaussian crash or evidence of MP2 spin contamination.
It blocked the S–S channel. Changing the running input would not retroactively
recover that rejected reaction.

The user requested replacing UMA with inexpensive DFT and cleanly restarting
both PES. The reusable `theory_preset="b3lyp-b2plyp-anl"` now selects:

| Level | Method | Optimization/frequencies |
| --- | --- | --- |
| L1 | Gaussian B3LYP-D3(BJ)/cc-pVDZ | Sella, then native Gaussian Hessian |
| L2 | Gaussian B2PLYP-D3(BJ)/cc-pVTZ | Sella, then native Gaussian Hessian |

Gaussian selects the unrestricted form for these open-shell multiplicities;
closed-shell singlets are not forced unrestricted. Both retain XQC and
UltraFine. The original UMA presets remain available and their implicit ANL
selection is unchanged. Neither this preset nor the recovery contains sulfur,
species-ID or bond-specific conditions. The Gaussian manual lists B3LYP and
B2PLYP dispersion parameters and DFT derivatives: [Gaussian DFT reference](https://conflex.co.jp/gaussian_support/dft.php).
This establishes supported methods, not successful convergence for these PES.

A generated-worker audit caught additional profile propagation defects:
legacy AM1 preliminary TS/ring searches could override the DFT profile;
legacy scan/reference defaults could combine MP2 with DFT dispersion;
reaction/IRC writers could omit the resolved profile's dispersion settings;
and the reaction-family final Sella TS writer lacked the native frequency-mode
format argument. Those paths now retain explicit profiles and use their
correct frequency policy. The old non-profiled preliminary-method behavior is
preserved. Explicit scan overrides still apply. The legacy `_well_mp2` job
name is unchanged, but its input uses the resolved scan-reference method.

The original invalid SCF-list repair remains in place. The new tests inspect
actual generated well, scan-reference, conformer, ring, constrained TS, final
TS and IRC workers, as well as rendered native inputs. They do not certify
Gaussian's runtime behavior.

## General production behavior

The feature is in `kinbot/minimum_recovery.py` and the ordinary generated
Sella minimum worker. It contains no element, molecule, bond or species-ID
dispatch. It currently supports **Gaussian + external Sella/BFGS optimization
+ native Gaussian Hessians**, for wells and minimum conformers. It does not
apply minimum descent to transition states or constrained VRC scans.

1. Optimize normally and evaluate a fresh fixed-geometry Hessian.
2. If a minimum fails the ordinary frequency criterion and has a mode below
   `-imagfreq_threshold`, read its lowest imaginary mode. Check atom order,
   finite values and geometry agreement; rotate the native mode into the
   optimizer frame using a proper rotation. No atom permutation or reflection
   is applied.
3. Retain the rejected `.com`, `.log`, `.chk`, `.traj` and `_sella.log` in a
   unique `<label>_minimum_recovery/attempt_N_*/` directory. Evaluate both
   displacement signs with separate copies of the same rejected checkpoint.
4. Select the lower trial only if it lowers the energy by more than 10⁻⁶ eV.
   A native error or two non-descending trials stops recovery. This tiny
   energy guard prevents a numerical tie from being called descent; it is not
   an accuracy target or a change to frequency acceptance.
5. Reinitialize the optimizer at the selected geometry, retaining its trial
   checkpoint, then evaluate another Hessian. Only the normal optimization
   and frequency checks can accept the resulting point. Existing downstream
   species/connectivity checks still apply.

Ordinary JSON controls (defaults shown):

```json
{
  "minimum_recovery_attempts": 2,
  "minimum_recovery_displacement": 0.02,
  "imagfreq_threshold": 50.0
}
```

`minimum_recovery_attempts=0` disables recovery. The full Cartesian displacement
norm is in Å, not Å per atom; attempt two uses twice the first displacement.
Each attempt adds two force/energy calculations and another optimization and
Hessian. It stays within the same QC allocation. A clean start uses this path
without edited starting geometries or a separate manual displacement step.
Parser errors, nonfinite data and optimizer failures do not trigger mode descent.

The accepted/error pickle remains the driver handoff. A small `.result.json`
also exposes status and recovery count to the read-only monitor. Each recovery
has `recovery.json` with both geometries/energies, selection, displacement and
failure reason. The monitor shows those reports without treating retained
rejected Hessians as fresh failures. It still includes both prepared PES trees.

The conformer generator also now retains the basename checkpoint for Gaussian
Sella workers, matching its `Guess=Read` usage and the calculator's working
directory. Native Gaussian conformer checkpoint policy otherwise stays intact.
The templates accept the previous formatter's fields as well: new controls
travel through a private entry removed from Sella options before optimizer
construction. An already running controller can render the updated template
with the default controls; it does not need new formatting arguments. A
previously generated worker script is not rewritten by pulling this change.

## Restart and retry policy

Use [the production restart script](../examples/anl/tools/restart_on_blodgett.sh)
after pulling `composite`. It prepares and verifies `_r3_dft`, including six
orbital seeds, resolved input hashes, both-level Gaussian previews, dependency
checks and small native RotdPy/MESS tests **before** stopping old work. It then
cancels only jobs/processes belonging to both `_r2` PES trees, waits for them
to disappear, and renames their shared parent to `.stopped.<UTC timestamp>`.
Unrelated calculations and shared scratch are untouched. Fresh detached
controllers run on the login host; only QC/sampling allocations use Slurm.
One monitor of `_r3_dft` follows both PES. Archive paths and retirement
inventories are written into the replacement directory.

| Condition | Production behavior |
| --- | --- |
| SCF difficulty | Existing Gaussian `XQC` fallback, at unchanged method/spin. |
| Gaussian minimum converged but significant imaginary frequency | Up to two mode-descent attempts, both signs tested, checkpoint isolation, new optimizer and Hessian each time. |
| Non-descending trials, missing/mismatched mode, native error, exhausted budget | Record the cause and fail the point; no method substitution or relaxed frequency acceptance. |
| Sella step limit/parser failure/nonfinite result | Explicit error result; not automatically treated as an imaginary-mode problem. |
| Intermediate constrained TS path optimizer exception | Existing finite-input-geometry fallback is retained; a final TS still needs its own optimization/frequency/IRC checks. |
| Recognized Slurm-job I/O failure | Existing allocation-local guard defaults to two retries (`KINBOT_SLURM_IO_RETRIES`); scoped artifacts and attempt journals are retained. It does not classify ordinary QC nonconvergence as I/O. |
| Hindered-rotor recovery | Existing finite retry budget (`hir_recovery_attempts=2`, step multiplier 2); independent point validation remains required. |
| Rotd sampling job retry | Existing `rotdpy_max_retries=2`; failed/convergence-incomplete flux is not accepted into MESS. |
| Failed launch preflight or incomplete teardown | Stop the transition before replacement submission; show the cause. |

A clean start uses the recovery without manual geometry edits. Native Gaussian
acceptance remains to be demonstrated on Blodgett. Do not turn a normal native
termination or an empty queue into an end-to-end success claim.

An **optional** isolated check helper is retained at
`examples/anl/tools/check_minimum_recovery.py`. It accepts an input
JSON, an earlier generated worker and a new check directory, plus `--level l1`
or `--level l2`. It reads the original undisplaced geometry, charge and
multiplicity with AST literal parsing and uses the ordinary production job
writer. It never executes the old Python script or injects check results into
an existing PES. This check is not a prerequisite for the requested full restart.

## Local validation scope

Tests execute generated well, conformer and TS workers with substituted QC
and optimizer calculations, but real native-output energy/frequency/mode readers,
frame alignment, checkpoint branching and pickle publication. They cover descent,
fresh-Hessian acceptance, repeated rejection, uphill trials, QC exceptions,
disabled recovery, malformed settings, nonfinite spectra, rotated coordinates,
reflection/geometry mismatches, and conformer subdirectory checkpoint placement.
The native check helper is tested for unchanged starting geometry, unchanged
source files, normal generator use and duplicate refusal. The monitor test
includes two PES trees with simultaneous jobs and a nested recovery check.
These tests do not run Gaussian locally or certify native convergence.

## Final input audit amendments

The user requested VeryTight geometry accuracy and a single production input.
Both PES inputs now explicitly use `sella_fmax=0.0001` eV/angstrom, and the
new DFT preset supplies that default when not overridden. Automatic L3 had
inherited the validation graph's loose 0.03 value; `anl_geometry_fmax=0.0001`
is now a general validated parameter propagated through stable-species and TS
plans to every automatic L3 optimization. This is a maximum-force norm at the
VeryTight scale, not Gaussian's complete Berny stopping test. Gaussian SCF
uses XQC,Conver=10 with UltraFine; ANL SCF/CC controls are 1e-10/1e-8.
No tolerance is made less strict by a minimum-recovery attempt.

All channel/orbital/state controls are inline in each PES JSON, with explicit
indices and CAS(5,4)/(7,5) labels. Relative seed paths are resolved before
workers change directory. The previous checkpoint continuation dropped the
explicit alignment flag; that defect is fixed, and continued records retain
natural-orbital type, updated packed indices, hash and alignment choice.
The general recon installer writes inline selections too. Obsolete duplicate
surface/state JSON files are removed; scientific native outputs remain intact.
RotdPy's requested statistical flux target is 1 percent, with the explicit
finite sampling caps retained and convergence checked rather than assumed.
