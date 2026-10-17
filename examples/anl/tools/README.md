# CASPT2-only KinBot production handoff

The handoff preparation is complete. The two end-to-end PES calculations
have **not** been run locally or certified in advance. The user handles SSH
and launches them on Blodgett; subsequent failures/results will be reviewed
through completion. Example/configuration cleanup is complete; purging historical native evidence
remains deferred until production finishes.

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
identity across equivalent bond-order perception. The user selected saved-results
resume in `_r5_dft`, with an external archive of obsolete attempts and fresh
controller logs, using the commands below. Native confirmation of these repairs
and completion of ANL/RotdPy/MESS remain pending. The barrier cutoff remains
1000; it did not cause these failures and does not prevent failed TS searches.

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
statistical target is **1%**. Final L1/L2 and automatic L3 force thresholds
are **0.0001 eV/angstrom**, with SCF Conver10 / Molpro 1e-10 and CC 1e-8.
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
aDZ-F12, a constrained relaxed CASPT2-F12/aTZ-F12 target-root scan, and
aDZ/aTZ spin gaps extrapolated with **exponent 3.7** and anchored to ground
high-spin ROHF-UCCSD(T)-F12b/aTZ-F12. **IPEA=0.25, real SHIFT=0.1 Eh** apply
to all new CASPT2 stages. MRCI is excluded. The search cutoff is 1000 kcal/mol.
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

## Resume accepted results after the downstream repairs

The helper checks the saved `_r5_dft` run, including rendered inputs and small
native backend tests. Both PES must be inactive; otherwise it leaves their
jobs untouched and refuses migration. It retains accepted calculations and
genuine reaction rejections, archives obsolete attempts/controller logs in a
dated `.repair_archive.*` sibling, updates the experimental levels and SO
corrections, archives obsolete derived thermochemistry/MESS files,
and resumes both controllers detached on the login host. SH's successful VPT2
is reparsed without QC. Failed H2S VPT2 and the invalid first VRC scans are
restaged. HSS's accepted L2 reference is retained with checked identity handling.
Full reports are saved; terminal output summarizes the archive and launch PIDs.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/tools/resume_on_blodgett.sh
BASH
```

Each PES has one active worker and a 10-allocation cap: at most 20 allocations
across the two runs. All eleven known bad nodes are excluded:
`node[1,9-10,13-14,16,21,23-24,30,33]`. Actual scheduler exclusions remain a
runtime check. Preparation saves complete resolved JSON and input hashes.
The helper refuses duplicate launch. A preflight, active-job or unreviewed-failure
check prevents launch; unrelated jobs and shared scratch remain untouched.
The archive manifest records each original path. `resume_report.json` and
`resume_provenance.json` record the migration. Accepted database/native records
remain in place. Small failed-job markers retain negative cache decisions;
the monitor identifies hash-matched historical outcomes separately. Old failure
history remains accessible in the archive. Reconstructed rejection decisions
can still appear in new logs; no scientific result is relabeled as successful.

`restart_on_blodgett.sh` remains an optional **fresh r6 replacement** command;
it is not the saved-results restart selected here.

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python -m kinbot.run_monitor /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016_r5_dft'
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
The 1% Monte Carlo target is statistical precision, not a total kinetic-error bound.

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
