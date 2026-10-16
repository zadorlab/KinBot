# One-time Molpro confirmation handoff

The local selection is complete: **CH3S C–S CAS(5,4), CH3SS C–S CAS(7,5),
CH3SS S–S CAS(7,5)**. The two CH3SS selections contain different orbitals.
These are the minimum chemically supported spaces for corrected VRC sampling,
pending native Molpro confirmation. Full KinBot production has not run yet.

## Start on Blodgett

On the existing `composite` checkout, with its existing virtual environment:

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/orbital_reviews/molpro_confirmation/launch_on_blodgett.sh
```

No archive extraction is needed after pulling. The launcher validates the three
existing reviewed WFU seeds and generates all new inputs before stopping the
old campaign. It preserves native results in place, writes a dated inventory
and metadata backup, stops only the superseded campaign's controller and jobs,
then launches the new confirmation. It refuses to overwrite an existing new
run. A scheduler outage or possible reusable completed point stops the
transition for review rather than submitting duplicates.

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py monitor /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016'
```

The new run has **39 jobs in nine chains**: 30 CASPT2-F12 points across aDZ-F12
and aTZ-F12, plus nine aTZ-F12 MRCI-F12+Q points. They check 30 Å asymptotes,
representative radial/angular behavior, and exact returns. Each chain starts
with its far-point gate. No new geometry optimization or production correction
scan is included. I/O/startup/allocation retries are bounded and recorded;
scientific nonconvergence is not accepted as success.

Eight ranks request 104000 MB/node, below the reported 126000 MB limit. Known
bad nodes are excluded in the catalogue. CASSCF maxit is capped at Molpro's
999 limit; SCF/correlation limits are 2000. Correlation energy/coefficient
criteria are 1e-8/1e-5, giving THRVAR=1e-10. The existing relaxed-Davidson
`CLUSTER=0; ENERGD(1)` convention and `CABS_SINGLES=1` are unchanged. Read the
[package instructions](examples/anl/orbital_reviews/molpro_confirmation/README.md)
for resume, export, root/state conventions and native acceptance checks.

## Twelve deliverables

| # | Deliverable | Status at handoff |
| --- | --- | --- |
| 1 | Literature-grounded minimum-space selection | Local decision complete; native sampling-basis cost pending |
| 2 | ORCA chemistry, radial/angular and return checks | Complete; 47 curated native points |
| 3 | Minimal aDZ/aTZ Molpro confirmation | Prepared and locally tested; run next |
| 4 | General JSON orbital/state/solver/resource controls | Implemented and locally tested; native use pending |
| 5 | General orbital reconnaissance procedure | Revised, documented and locally tested |
| 6 | Both production KinBot inputs | CAS counts configured; install native orbital/surface selections after review |
| 7 | User-run production with continued bug-fix support | Follows Molpro review; no end-to-end result claimed |
| 8 | Generality and integration audit | Complete locally: 1,622 tests plus 391 subtests pass |
| 9 | Reproducible documentation and native evidence | Archived, checksummed and independently replayed |
| 10 | Full-directory cleanup inventory | Recorded; cleanup deferred until production completes |
| 11 | Consolidated handoff | Molpro commands ready; final production handoff follows acceptance |
| 12 | Preserve/reuse/retire old Blodgett work | Scoped adapter tested; executes with launch above |

The old Blodgett campaign is a standalone example study. General KinBot recon
can prepare and report equivalent declared checks; it supports the chemical
decision rather than automatically proving a unique active space. After native
review we install the basis-specific orbital/state models and proceed with the
ordinary KinBot inputs, supporting fixes through CBH-ANL, RotdPy and MESS.

The maintained [detailed checklist](docs/ORBITAL_RECON_DELIVERABLES.md),
[scientific decision](docs/SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md),
[general procedure](docs/ORBITAL_REVIEW_WORKFLOW.md), and
[cleanup inventory](docs/REPOSITORY_CLEANUP_PLAN.md) record the scope and caveats.

## Retained artifacts and verification

- `molpro_confirmation_transfer_20261016.tar.gz` at the repository root is
  committed with its SHA-256 sidecar. It contains the runtime source, handoff
  package, documentation, test evidence and compact ORCA decision evidence,
  including native text inputs/outputs for the 47 curated points. Its internal
  `transfer_manifest.json` hashes every payload file. Extraction is optional
  and should use a fresh directory when inspecting this snapshot.
- The full native archive, including orbital binaries/coefficients, seeds and
  analysis dependencies, is retained locally at
  `local_runs/qchem_casscf_20261015/reproducibility/orca_minimum_space_final_20261016.tar.gz`.
  Its SHA-256 is
  `13762a1124092834313eee3bf463383f9c5658535351715dc4efdf981feca523`.
  This 311 MiB artifact is not duplicated in Git. Its manifest/checksum and
  verification report are included in the compact transfer bundle.
- Full-archive verification checked 1,469 files, then used isolated extracted
  analysis source to reproduce all 47 curated observations exactly, recompute
  return/active-span diagnostics and plots, and generate fresh replay inputs.
  No new electronic-structure calculation was claimed by this replay.
- The full local regression passed **1,622 tests and 391 subtests**. Shell
  syntax, bundle member hashes and the launcher/retirement scenarios were
  checked. Native Molpro behavior and current cluster health remain to be
  established by this confirmation; no zero-error guarantee is made.

The compact bundle intentionally omits orbital binaries and is not a replacement
for the full evidence archive when repeating orbital-overlap analysis. Existing
Blodgett WFU seeds are reused only after exact hash checks. No old outputs,
failed-attempt evidence or scratch directories were deleted during handoff.
