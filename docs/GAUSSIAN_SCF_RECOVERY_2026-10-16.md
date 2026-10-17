# Gaussian QC-SCF recovery and r4 restart

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
Failure of the retry propagates normally rather than looping or accepting an
unconverged result. Subsequent calculation invocations start with the original
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

## Deployment and clean restart

The user selected a clean restart of **both** examples rather than preserving
running jobs. After pulling, `examples/anl/tools/restart_on_blodgett.sh`
prepares and preflights `_r4_dft`, stops controllers and allocations under both
`_r3_dft` PES directories, verifies that they have stopped, and moves the old
parent to a dated `.stopped.*` archive. Both fresh controllers run detached on
the login host. Existing per-PES allocation caps, exclusions, scientific input
parameters and binary orbital checkpoints are unchanged. No unrelated jobs or
shared scratch are removed. The parent-directory monitor still covers both.

An already running Python worker cannot acquire this fix by pulling alone.
New Gaussian worker processes import it automatically; restarting both ensures
that both drivers also reconsider previously rejected channels.

Local regression tests cover trigger discrimination, bounded success/failure,
disable control, preserved thresholds and artifacts, actual native-route
rendering through the calculator, monitor visibility, and scoped two-PES
retirement/launch helpers. Native Gaussian execution of this retry remains to
be verified on Blodgett. Local tests are not a full production-completion claim.

Validation: **90 passed in 11.97 seconds**; `bash -n` also passed for the
restart script. Reproduce the focused suite with:

```bash
.venv/bin/python -m pytest tests/test_gaussian_scf_recovery.py tests/test_gaussian_runtime.py tests/test_gaussian_input_routes.py tests/test_run_monitor.py tests/test_caspt2_handoff.py tests/test_run_cleanup.py -q --disable-warnings --tb=short
bash -n examples/anl/tools/restart_on_blodgett.sh
```
