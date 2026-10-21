# CASPT2-only KinBot production handoff

**Current handoff, 21 October: fresh r7 with ANL components v2.** Use
[`fresh_run.py`](fresh_run.py) and the editable [`fresh_run.json`](fresh_run.json).
The [complete component matrix, audit limits, restart and monitor commands](../../../docs/ANL_COMPONENTS_V2_AND_R7.md)
supersede the older r5/r6/recovery scripts below. The new helper preflights
before stopping r6, archives its evidence, uses pinned code, shares one node
pool, and never repeats a recorded driver launch. This campaign configuration
also opts into emptying `~/wfu` and removing verified old home scratch after
the account's jobs stop; see the cleanup scope and safeguards in the handoff.

## Historical handoffs (do not use these commands for the new r7)

**Interrupted cluster handoff, 17 October:** use
`bash examples/anl/tools/recover_on_blodgett.sh` after pulling `composite`.
This supersedes the r6 launch command: it preflights fresh r7 inputs, stops and
archives r6 (and r5 if still present), then starts the two r7 drivers. It runs
detached and prints a unique progress-log path immediately. The details and
retry behavior are in the [recovery handoff](../../../docs/PRODUCTION_HANDOFF_RECOVERY_2026-10-17.md).

The handoff preparation is complete. The two end-to-end PES calculations
have **not** been run locally or certified in advance. The user handles SSH
and launches them on Blodgett; subsequent failures/results will be reviewed
through completion. Example/configuration cleanup is complete; purging historical native evidence
remains deferred until production finishes.

Latest local verification: **1930 tests + 424 subtests** passed, followed by
**53 tests + 11 subtests** after the 5000 K screening/input change. The
[current validation record](convergence_restart_validation.json) includes
commands, native smoke-test coverage and limits. Passing these checks does
not certify the remaining native Gaussian/Molpro production calculations.

**Repair update, 16 October 2026:** the first launch at 18:31 PDT (`c340cc89`)
passed the five small backend preflight tests and archived nine old campaigns,
but both initial L2 jobs failed at 18:35 with a malformed Gaussian SCF option
list. The input merger/writer is repaired; PES now records worker failure
instead of attempting postprocessing. A clean replacement and detailed monitor
are supplied below. See the [failure evidence and audit](../../../docs/GAUSSIAN_PRODUCTION_REPAIR_2026-10-16.md).

**Current restart update, 17 October:** `_r5_dft` reached CBH-ANL and the VRC
relaxed scan, then failed in interfaces documented in the
[downstream repair record](../../../docs/PRODUCTION_DOWNSTREAM_REPAIR_2026-10-17.md).
The repairs place WFU loading before the requested Molpro geometry in every
production renderer, separate Gaussian VPT2 framework labels from its symmetry
cap, bypass VPT2 symmetry processing, and preserve requested CBH reference
identity across equivalent bond-order perception. After a subsequent SH SCF
failure and a resume bug that moved fragment checkpoints, the user selected a
**full fresh restart**, with no calculated results reused between the PES or
from r5. The [convergence and restart record](../../../docs/PRODUCTION_CONVERGENCE_RESTART_2026-10-17.md)
documents the repairs and the new numerical settings. Native confirmation of these repairs
and completion of ANL/RotdPy/MESS remain pending. The previous 1000 kcal/mol cutoff did not cause these interface failures.
The new user-selected screening policy uses 400 kcal/mol at L2, with a
10 kcal/mol L1 allowance, for temperatures through 5000 K.

This repair passed **403 regression tests and five native RotdPy/MESS/MESSPF
checks**. [Commands, source hashes and limits](downstream_repair_validation.json)
are retained. The new Gaussian/Molpro repair inputs still need native Blodgett
confirmation; this is not a claim of completed kinetics or absence of all bugs.

**Previous atomic-SCF repair:** `_r3_dft` passed the initial minima, then rejected
CH3S C–S homolysis when the sulfur atom's Gaussian QC-SCF line search aborted.
The [general bounded SCF retry](../../../docs/GAUSSIAN_SCF_RECOVERY_2026-10-16.md)
retains Conver10 and the requested theory and adds a CDIIS/full-integral-accuracy
retry. The user selected fresh `_r4_dft` runs for both PES, archiving `_r3_dft`.
The initial sulfur L1 retry subsequently proved unconverged despite normal
termination. The reader now rejects that false success. A native isolated
atomic-symmetry test genuinely converged in 9 cycles at Conver10; general
single-atom handling now uses that policy, with `atomic_scf_symmetry=false`
available in Gaussian kwargs. Molecular settings are unchanged.

The 22:17 PDT `_r4_dft` snapshot showed CH3S C-S **accepted**, a later converged
sulfur conformer, and normal sulfur L2 completion. The user subsequently chose
clean `_r5_dft` runs overnight. Keep the earlier unconverged L1 log as history.

**Earlier minimum repair:** the `_r2` CH3S L2 minimum failed its Hessian;
independent displacements confirmed negative curvature. CH3SS passed its
initial L2 job but later rejected the CH3S UMA product at L1. The new
`b3lyp-b2plyp-anl` preset uses B3LYP-D3(BJ)/cc-pVDZ at L1 and retains
B2PLYP-D3(BJ)/cc-pVTZ at L2. Both Gaussian levels have bounded minimum
recovery with fresh Hessian validation. The [current evidence and retry plan](../../../docs/MINIMUM_RECOVERY_2026-10-16.md)
records what is fixed and what still requires native execution.
Unconstrained automatic ANL L3 minima now also have bounded curvature/descent
recovery before dependent energies. This adds 6N force evaluations per check;
the independent native harmonic check remains required. All recovery controls
are explicit in the two inputs, with no changes to spin or electronic levels.

The current inputs are **one JSON per PES**: all channel, state, solver and
orbital-index choices are inline. The obsolete separate surface/state JSON
files have been removed. Only the six hash-checked binary orbital checkpoints
are staged as data. Methods explicitly read CAS(5,4) or CAS(7,5). RotdPy's
statistical target in the base example JSON is **1%**. The r7 campaign's
`fresh_run.json` overrides it to **5%** and sets 47 relaxed-scan points
(2–8 angstrom by 0.2, 8.5–12 by 0.5, 13–20 by 1). Its ROTD dividing surfaces
and separate 30-angstrom asymptotic reference are unchanged. Final L1/L2,
automatic L3 and constrained VRC force
thresholds are **0.0005 eV/angstrom**, with SCF Conver8 / Molpro 1e-8 and
CC/CASPT2/MRCI energy convergence 1e-7. Molpro coefficient tolerances remain
separate: CC squared-amplitude and CAS/MRCI `THRVAR` are 1e-10. The reviewed
CASSCF orbital energy/gradient tolerances remain 1e-8/1e-5.
Preflight renders both Gaussian profiles and every declared CASPT2
channel/root/basis before stopping the old run.

Every CAS/MRCI evaluation now retains [native orbital/state diagnostics and an
extracted report](../../../docs/MOLPRO_NATIVE_ORBITAL_DIAGNOSTICS.md) beside the
output, including RotdPy samples. It records occupations, printed AO coefficients,
CI/state information and hashes; it does not infer no-swapping from dimensions.

Previous full regression: **1,856 tests and 424 subtests passed**, zero failures
or skips, including native backend checks. The [validation record](dft_restart_validation.json)
retains the command, source hashes and compressed test report. This is local
verification; native Gaussian/Molpro recovery and production remain to be run.

## Accepted inputs and scope

| Channel | CAS | Active natural orbitals | SA spatial roots per spin (D/Q) |
| --- | --- | --- | --- |
| CH3S C–S | (5,4) | 11–14 | 3/3 |
| CH3SS C–S | (7,5) | 18–22: S2 pi/pi-star + methyl radical | 1/1 |
| CH3SS S–S | (7,5) | 18–22: CH3S sulfur pair + departing S 3p | 6/6 |

These one-based indices refer to six supplied, hash-bound, basis-specific
natural-orbital records. The two CH3SS channels have different subspaces.
The retained sparse radial/angular/return evidence supports these minimum
spaces; it does not prove that no uncomputed geometry can change state character.

Production uses ordinary **single-bottleneck VRC**, CASPT2-F12 sampling at
aDZ-F12, a constrained relaxed CASPT2-F12/aDZ-F12 target-root scan, and
aDZ/aTZ spin gaps extrapolated with **exponent 3.7** and anchored to ground
high-spin ROHF-UCCSD(T)-F12b also extrapolated aDZ/aTZ with exponent 3.7.
The CC reference and gap method/basis/exponent settings are independent input
controls. MR asymptotic supermolecule geometries supply the rigid fragments;
fixed-geometry auxiliary L2 orbitals supply HOMO pivot directions. See
[the targeted handoff](../../../docs/VRC_MR_FRAGMENT_ORBITALS.md).
**IPEA=0.25, real SHIFT=0.1 Eh** apply
to all new CASPT2 stages. MRCI is excluded. The L2 search cutoff is 400 kcal/mol (L1: 410); MESS/RotdPy cover 300–5000 K.
MESS ModelEnergyLimit is 800 kcal/mol, separate from the search cutoff.
ANL1-F12/CBH thermochemistry is a separate, higher-level workflow.

The user subsequently verified all 27 H/C/S basis-library lookups on Blodgett
Molpro 2024.1, including augmented QZ-F12 and 5Z-F12. No production basis
change is required. See the [library-check record](../../../docs/MOLPRO_BASIS_LIBRARY_CHECK_2026-10-16.md)
for the exact coverage and its distinction from F12 execution validation.

Every CASPT2 root uses the full declared same-spin reference list and common
weights. The general `energy_selection: {"4": "lowest"}` option evaluates
all declared quartet roots for the correction and selects their minimum:
it does not assume CASSCF root 1 is always the lowest correlated root.
All roots must succeed. It does not provide universal root-character tracking,
add XMS mixing, or guarantee a correct ground-state CC solution. Numerical
geometry gradients remain bound to explicitly selected target roots.

## Documentation and evidence

- [Full workflow, expected candidate chemistry, resources and commands](../ch3s_anl1_f12/orbital_reviews/RUNNING.md).
- [Every explicit parameter, nested control, consequential default and master/branch provenance](../ch3s_anl1_f12/orbital_reviews/PARAMETERS.md).
- [Maintained 12-item status checklist](../../../docs/ORBITAL_RECON_DELIVERABLES.md).
- [Orbital occupations, chemical character and native evidence](../ch3ss_anl1_f12/orbital_reviews/joint_history/docs/SULFUR_ACTIVE_ORBITAL_EVIDENCE.md).
- [Matched shift-0.1 native acceptance](../ch3ss_anl1_f12/orbital_reviews/joint_history/docs/CASPT2_SHIFT01_NATIVE_REVIEW_2026-10-16.md).
- [General orbital-reconnaissance procedure](../../../docs/ORBITAL_REVIEW_WORKFLOW.md).
- [Deferred cleanup inventory](../../../docs/REPOSITORY_CLEANUP_PLAN.md).

Local verification: **571 tests + 31 subtests**, plus **five separate native
RotdPy/MESS/MESSPF checks**. The latter retain native files in
[the backend evidence package](../ch3ss_anl1_f12/orbital_reviews/joint_history/production_checks/backend_smoke_review.json).
The new minimum-selector wrapper is regression-tested; its individual CASPT2
root calls have native confirmation, but this wrapper itself has not yet
executed in Molpro. No full sulfur VRC flux or MESS network is claimed complete.

The subsequent full production-repair audit passed **1,748 tests and 424
subtests**, with **zero failures or skips** and native RotdPy/MESS/MESSPF
checks enabled. These counts describe separate runs, not additive coverage.
The [repair validation record](../ch3ss_anl1_f12/orbital_reviews/joint_history/production_checks/repair_validation.json)
records the command, source hashes and coverage limits.
Cluster preflight now also renders 72 Gaussian input previews from the actual
production profiles. Native Gaussian execution remains a Blodgett check.

Atomic sulfur, CH3S and SH fine-structure levels are forwarded to MESS/PF using its existing
electronic-level capability. The KinBot JSON forwarding and explicit term-
centroid energy hook were added on this branch; neither invents a VRC SOC
potential. No inner/outer stages or automatic equilibration are enabled.
S, SH and CH3S request that correction once in ANL: -195.756111, -188.4175
and -127.5 cm^-1, respectively; CH3S uses the approved effective approximation. The
[species review](../../../docs/SULFUR_ELECTRONIC_LEVELS_2026-10-17.md) records
experimental sources, the complete discovered inventory, the CH3S JT caveat,
and separate energy-zero/partition-function uses. Transport parameters remain
provisional example estimates.

## Fresh two-PES restart after the downstream repairs

The helper prepares `_20261017_r6_dft_scf8`, checks the rendered inputs and
runs the small native RotdPy/MESS/MESSPF tests before stopping anything. It
then stops only both r5 PES trees, waits for their jobs/processes to disappear,
and archives the complete r5 parent to a dated `.stopped.*` sibling. Both fresh
controllers detach on the login host. The new run contains no old database,
failure markers, energies or thermochemistry models. The reviewed WFU files
remain orbital-selection inputs. No calculated results are shared across PES.
Earlier external repair archives remain intact. Full reports are saved;
terminal output summarizes the fresh paths and launch PIDs.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/tools/restart_on_blodgett.sh
BASH
```

Each PES has one active worker and a 10-allocation cap: at most 20 allocations
across the two runs. All eleven known bad nodes are excluded:
`node[1,9-10,13-14,16,21,23-24,30,33]`. Actual scheduler exclusions remain a
runtime check. Preparation saves complete resolved JSON and input hashes.
The helper refuses duplicate launch. Failed preflight or incomplete retirement
prevents launch; unrelated jobs and shared scratch remain untouched.
`production_retirement.*.json` records the scoped inventory and archive;
`launch_report.json` records the new controllers. Old native evidence remains
accessible in the archive. All reactions are reconsidered in the fresh search.

Both automatic branches write `anl_automatic/electronic_state_audit.json`.
It covers PES species, every planned CBH rung/reference, and terminal ATcT
identities, with input sources and any unreviewed open-shell states. The
[review](../../../docs/SULFUR_ELECTRONIC_LEVELS_2026-10-17.md) also records
the requested future experimental SO database; present data remain explicit
in each single input JSON.

The old saved-results `resume_on_blodgett.sh` is not used for this restart.

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python -m kinbot.run_monitor /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261017_r6_dft_scf8'
```

The parent directory monitors both PES. Both overviews appear first; routine
result records and dependency waits are summarized. For a full report, redirect
`--all --details` to a file rather than capturing a long watch screen.

For live native/controller log text instead of snapshots, run the same Python
monitor command with `--follow`. `--all` expands reaction rows and log coverage;
`--details` expands dispatch records. The
[monitor guide](../../../docs/RUN_MONITOR.md) explains status meanings and bounded reads.

After launch, review optimized structures/IRC, accepted channels, ANL component
and formation-enthalpy closure, relaxed corrections/asymptotes, root/subspace
behavior, Monte Carlo convergence, grid-boundary minima and final MESS rates.
The Monte Carlo target (5% for r7) is statistical precision, not a total
kinetic-error bound or confidence level.

## Example layout (current)

The `caspt2_production` and acetaldehyde example folders have been removed.
Each sulfur example has exactly one root JSON plus an `orbital_reviews/`
subdirectory. That subdirectory contains its own runtime WFU files, CASPT2
justification, parameter reference and native evidence. No archive extraction
or state-model JSON is required to run the input. The shared tools here only
prepare, check, stop/archive, launch and monitor ordinary KinBot runs.
The original mixed historical exports are retained once under CH3SS
`orbital_reviews/joint_history`, with links from both reviews; neither input
reads them. [layout_manifest.json](layout_manifest.json) records relocated paths.

See the current [DFT restart validation](dft_restart_validation.json).
The full local suite was rerun after this layout change; the final result is below. Native Gaussian/Molpro
production completion cannot be certified by these local tests.

## Final local validation for this restart

The final reorganized tree passed **1,818 tests and 424 subtests**, with
zero failures or skips and two warnings, including five native RotdPy/MESS/
MESSPF checks. [Commands, hashes, logs and coverage limits](dft_restart_validation.json)
are retained. Each input was also tested after copying only that JSON and
its own WFU data, without any other example or configuration file.
Native Gaussian/Molpro production execution remains a Blodgett check.
