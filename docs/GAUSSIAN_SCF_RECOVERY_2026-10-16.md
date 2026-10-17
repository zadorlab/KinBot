# Gaussian atomic SCF and convergence validation

## Latest native result and correction (22:16 PDT)

The initial CDIIS retry was **not sufficient**. In `_r4_dft`, the sulfur
`320000000000000000003_well.log` printed `Convergence criterion not met`,
RMSDP=3.86e-6 and MaxDP=3.22e-5, yet exited normally. KinBot's reader accepted
the printed energy and the retry was incorrectly recorded as complete.
The exit-code-only recovery test missed this case. This was a KinBot
acceptance bug, not evidence that the requested threshold was satisfied.

Two isolated tests at `/home/hram/kinbot-runs/sulfur-scf.vm5DG2` separated
the issues. `Symmetry=None` stalled at RMSDP=1.21e-9 and failed in Link 508.
Forcing `Symmetry=(PG=D2h)` instead failed before SCF with
`DN error #2 in IniRep` in Link 202; that diagnostic input was inappropriate
for this atom/build and is not used in the fix.

The corrected test at
`/home/hram/kinbot-runs/sulfur-auto-symmetry.XJjD7K/sulfur.log` used:

```text
#P UB3LYP/cc-pVDZ EmpiricalDispersion=GD3BJ Integral=UltraFine
SCF=(XQC,Conver=10,NoVarAcc,NoIncFock,MaxCycle=256,Symm)
```

With **no molecular Symmetry keyword**, Gaussian 16 B.01 detected OH with
D2h as its largest Abelian subgroup. Charge 0/triplet sulfur converged in
9 cycles: RMSDP=1.26e-11, MaxDP=1.32e-10, energy=-398.125080451 Eh,
S-squared=2.0022 before annihilation, normal termination and exit 0.
The requested RMS/max limits were 1e-10/1e-8. These values are transcribed
from the user's native report; full logs remain at the cluster paths above.
No change to spin, functional, basis, grid or convergence threshold was needed.
The result supports a symmetry/conditioning explanation, but does not reveal
Gaussian's exact internal source of the unconstrained iteration failure.

## General atom handling and strict reader

For **any one-atom Gaussian input**, the writer now removes the molecular
symmetry override, requests `SCF=Symm,NoVarAcc,NoIncFock`, and supplies
MaxCycle=256 only when no cycle limit was specified. It retains the requested
SCF algorithm, convergence, guess, electronic method, charge, multiplicity,
basis, grid and coordinates. Molecular inputs are unchanged. Users needing
their own atomic symmetry policy can set `"atomic_scf_symmetry": false` in
Gaussian calculator kwargs (globally or in a theory profile). This boolean
is a KinBot control and never appears in a Gaussian route. The same writer
is used by worker execution and input previews.

Gaussian documents automatic point-group detection and separately documents
SCF=Symm as preserving initial occupations within each Abelian irrep:
[Symmetry reference](https://www.conflex.co.jp/gaussian_support/symmetry.php),
[SCF reference](https://www.conflex.co.jp/gaussian_support/scf.php).
This does not claim that every possible atomic state/initial occupation is
the ground state; user-specified states still require appropriate review.

The reader now rejects an unresolved native SCF nonconvergence warning even
if `SCF Done` and normal termination follow it. An explicit subsequent QC
convergence message can resolve the warning, as in a successful XQC fallback.
A later geometry or linked calculation cannot hide an earlier unresolved
failure. The calculator clears cached results before reading; a rejected retry
is recorded as failed and cannot publish energy/forces as successful results.

## Whether r4 needs restarting

The 22:17:52 PDT monitor shows CH3S C-S homolysis **accepted**, with the later
sulfur conformer SCF explicitly converged and sulfur L2 normally terminated
after 3 SCF cycles. No rejected CH3S channels are shown. CH3SS initial L2
frequencies are still running. The selected conformer replaces the original
product properties, and `Optimize.compare_structures` loads the accepted L2
record through `load_calculation_record`. Thus this snapshot does not require
discarding the two running PES. It does not certify future channel outcomes,
the complete native L2 density history, or completed ANL/RotdPy/MESS work.
The original unconverged L1 file remains historical evidence, not a converged
result. The reported L1 energies differ by about 8e-9 Eh.

Pulling the patch leaves jobs untouched; **new** Gaussian worker processes
use it. Already running Python processes retain imported code. The current
CH3SS worker is molecular and its symmetry policy is unchanged. Optional
`restart_on_blodgett.sh` now preflights `_r5_dft`, stops/archives both `_r4_dft`
PES trees, and launches fresh drivers. It is not required merely to update
the monitor, nor does pulling alone rewrite old result/database records.

## Native observation

The r3 DFT production attempt reached initial L1/L2 minima in both examples.
CH3S C–S homolysis was subsequently rejected because its triplet sulfur product
`320000000000000000003_well` failed. Its native input requested
B3LYP-D3(BJ)/cc-pVDZ, charge 0, multiplicity 3, UltraFine integration and
`SCF=(XQC,Conver=10)`. The log ended at 21:16:34 PDT on 16 October 2026 with
`No lower point found -- run aborted.` and an error termination in `l508.exe`.
The repeated KinBot error lines were polling of this one failed generation.
This was an electronic SCF line-search failure on an atom, not a geometry
optimization or CAS failure.

The nearly flat line search suggests numerical conditioning/precision may
matter. The supplied output does not establish the precise internal cause,
spin contamination, or an incorrect spin state. It does not justify changing
the spin, method or accepting a nonconverged energy.

## General automatic retry

The KinBot Gaussian calculator now retries **once per calculation invocation**
when a failed process has both a native Link 508 error termination and either
the reported line-search abort or a convergence-failure diagnostic. It uses
`SCF=(Conver=10,CDIIS,NoVarAcc,NoIncFock,MaxCycle=256)` for these inputs.
Any explicitly larger cycle limit is retained. The caller's density threshold,
method, basis, dispersion, grid, geometry, charge, multiplicity and guess policy
are retained; there is no fallback electronic level or looser acceptance rule.
Unrelated syntax, I/O and other-link failures do not invoke this retry.
With the corrected native reader, failure of the retry propagates rather than
accepting an unconverged result. Subsequent calculation invocations start with the original
requested parameters.

Gaussian documents CDIIS as an alternative to QC, with damping implied;
NoVarAcc uses full integral accuracy from the outset, and NoIncFock disables
incremental Fock construction. Conver=N is a density convergence criterion,
not an energy-change criterion. These documented options motivate this
algorithm/precision retry; they do not guarantee convergence for this atom.
Source: [Gaussian SCF keyword reference hosted by distributor CONFLEX](https://www.conflex.co.jp/gaussian_support/scf.php).

The default is enabled for the Gaussian calculator, without molecule-specific
logic. A user can disable it with `"scf_recovery": false` in Gaussian
`calc_kwargs` or the relevant theory profile's calculator kwargs. This is a
KinBot calculator control, never a Gaussian route keyword.

Before retrying, KinBot copies the failed input, native log and available
checkpoint into `restart_gaussian_scf/<label>.<unique>/`. The adjacent
`<label>.scf_recovery.json` records the original/retry SCF settings, archive
paths and hashes, and retry outcome. The combined monitor displays that record
and excludes the archived failed log from current-job diagnostics.

## Historical r4 deployment

The user selected a clean restart of **both** examples rather than preserving
running jobs. The earlier `examples/anl/tools/restart_on_blodgett.sh`
prepared and preflighted `_r4_dft`, stopped controllers and allocations under both
`_r3_dft` PES directories, verified that they had stopped, and moved the old
parent to a dated `.stopped.*` archive. Both fresh controllers ran detached on
the login host. Existing per-PES allocation caps, exclusions, scientific input
parameters and binary orbital checkpoints are unchanged. No unrelated jobs or
shared scratch are removed. The parent-directory monitor still covers both.

An already running Python worker cannot acquire this fix by pulling alone.
New Gaussian worker processes import it automatically; restarting both ensures
that both drivers also reconsider previously rejected channels.

Local regression tests cover trigger discrimination, bounded success/failure,
disable control, preserved thresholds and artifacts, actual native-route
rendering through the calculator, monitor visibility, and scoped two-PES
retirement/launch helpers. Additional regression cases exercise the real output
reader with the native false-success excerpt, successful XQC, cleared caches,
failed retry metadata, and general atom/diatomic/override input routing.
Native atom success is the isolated test above, not a production-completion claim.

Validation of the atom/reader changes: **98 passed in 12.19 seconds**.
Related worker, minimum-recovery, selected-record, theory-profile and production
tests: **95 passed in 11.80 seconds**. Shell syntax and `git diff --check` pass.
Reproduce the focused suite with:

```bash
.venv/bin/python -m pytest tests/test_gaussian_scf_recovery.py tests/test_gaussian_runtime.py tests/test_gaussian_input_routes.py tests/test_run_monitor.py tests/test_caspt2_handoff.py tests/test_run_cleanup.py -q --disable-warnings --tb=short
bash -n examples/anl/tools/restart_on_blodgett.sh
.venv/bin/python -m pytest tests/test_sella_worker_results.py tests/test_minimum_recovery.py tests/test_calculation_records.py tests/test_theory_profiles.py tests/test_caspt2_production.py -q --disable-warnings --tb=short
```
