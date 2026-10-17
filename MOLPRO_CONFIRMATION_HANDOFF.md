# Historical Molpro confirmation handoff

**Superseded for submission:** Molpro 2023 reproduced the S–S MRCI DSYGV
failure. The user closed MRCI, and the matched CASPT2 shift-0.1 check passed.
Use the [CASPT2 production handoff](CASPT2_PRODUCTION_HANDOFF.md) for the next
launch. The commands below are retained as historical diagnostic provenance;
do not rerun these MRCI probes as production prerequisites.

## S–S MRCI build comparison after failed shift trial (historical)

The minimum spaces remain **CH3S C–S (5,4), CH3SS C–S (7,5), CH3SS S–S
(7,5)**. The CASPT2 confirmation and four simultaneous-root CH3S MRCI points
pass native review. The S–S far point has now completed according to the
controller. Center failed at the first simultaneous MRCI iteration with
`dsygv failed`; tilted is waiting for center's checkpoint. The stopped
seven-point run has **five validated, one failed, one waiting**.

The 0.1 Eh denominator-shift retry **also failed** at the same first-iteration
solve (job 52908243). Its tilted point has not run. The next submission is
one center point with the installed **molpro23** module, keeping the failed
trial's electronic input and original starting seed identical. It is a
controlled comparison of executables/builds, not a known cure or a decision
to downgrade production. All six roots per spin, CABS=3, SFIX, tolerances,
0.1 Eh shifts and relaxed-Davidson selectors are retained. See the
[native diagnosis and documentation review](docs/MOLPRO_MRCI_DSYGV_RECOVERY_2026-10-16.md).
No public release note found in this review establishes a fix for this error.

The user has selected **IPEA=0.25, CASPT2 real SHIFT=0.1 Eh, and MRCI
Davidson iteration denominator shifts of 0.1 Eh** for the eventual production
KinBot inputs. The failed retry used the requested MRCI shifts.
Completed CASPT2 confirmation remains documented at SHIFT=0.3 Eh. This is a
recorded production choice, not a change to existing runs or a claim of
validated shift insensitivity. The [decision record](docs/SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md)
lists the three reviewed orbital sets and the remaining native state audit.

## Submit the one-point build comparison

The user handles SSH. On Blodgett:

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/orbital_reviews/molpro_confirmation/probe_mrci_build_on_blodgett.sh
```

Source:
`/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_solver_retry`.
Destination:
`/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_molpro23_probe`.
The script cancels/deletes nothing, leaves other campaigns alone and refuses a
live source or duplicate destination controller. Do not run the historical
full-restart, recovery, CABS-gate, original multistate or shift-retry scripts again.

The script loads molpro23 in its own process, verifies a 2023 version, and
checks the generated batch setup also selects that version. It verifies the
source failure and copies the pristine original seed, never the failed center
WFU. Exact prepared-input and seed hashes must match the failed shift trial.
`build_comparison.json` records these checks and executable information.
All previous results remain intact; the completed S–S far spectrum still
needs independent native review.

The [plan JSON](examples/anl/orbital_reviews/molpro_confirmation/mrci_build_probe_plan.json)
changes no electronic controls and contains no tilted point. A successful
center would justify review and matched-build far/tilted validation before
production use. A checkpoint/input incompatibility in the older build is
distinct from reproduction of DSYGV. The test does not establish a code bug
or isolate a math-library defect. If DSYGV recurs, the next requirement is
LAPACK INFO/matrix diagnostics from a diagnostic build or Molpro support,
rather than another parameter sweep.

Resources remain **eight ranks, 1000 MW/rank and 104000 MB/node**, below the
reported 126000 MB node limit. Concurrency is one. Existing exclusions and
bounded infrastructure retries are inherited. CASSCF MAXIT=999; correlation
MAXIT/MAXITI=2000, MAXVI=60, input MAXDAV=20, THRDEN=1e-8, THRVAR=1e-10.
The failed six-root native call printed an effective MAXDAV=24. No default
settings for unrelated KinBot calculations are changed.

## Monitor and export

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py monitor /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_molpro23_probe --all --roots'
```

Expect one center point with two simultaneous calls (doublet and quartet),
twelve selected energies in total. Passing iteration 1 is progress, not final
acceptance. This is not a fresh active-space survey.

Export the stopped source now to retain the completed S–S far calculation and
failed center context. This is independent of the new retry:

```bash
cd /home/hram/KinBot
.venv/bin/python examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py collect \
  /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate \
  --output /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate.evidence.tar.gz
```

The shift retry has stopped. Export it with:

```bash
cd /home/hram/KinBot
.venv/bin/python examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py collect \
  /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_solver_retry \
  --output /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_solver_retry.evidence.tar.gz
```

After the build probe stops, collect its evidence as well:

```bash
cd /home/hram/KinBot
.venv/bin/python examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py collect \
  /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_molpro23_probe \
  --output /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_molpro23_probe.evidence.tar.gz
```

The user transfers the archives to `examples/anl/orbital_reviews`. Collection
refuses to overwrite an existing archive or collect a live/uncertain run.
Original historical failures remain in their monitors and are not new failures
of this trial.

## Deliverable status

| # | Deliverable | Current status |
|---|---|---|
| 1 | Literature-grounded minimum spaces and cost | Choices retained; native sampling timings documented |
| 2 | Local ORCA chemical/radial/angular/return checks | Complete; 47 curated points retained |
| 3 | Native aDZ/aTZ Molpro transfer | CASPT2 and four CH3S MRCI points reviewed; S–S far completed per controller, center failed twice and tilted unrun |
| 4 | General JSON controls | Simultaneous MRCI confirmed for CH3S; optional denominator controls added for the explicit numerical trial |
| 5 | General orbital reconnaissance | Procedure documented; chemistry remains in example data |
| 6 | Complete production inputs | CAS counts configured; final accepted stage models still to install |
| 7 | User-run end-to-end support | Remains in scope through CBH-ANL/RotdPy/MESS completion |
| 8 | Generality/integration audit | 25 focused tests pass; exact native-deck/seed preservation and module-path switching tested; no local Molpro execution |
| 9 | Reproducible evidence | Existing exports and CH3S audit retained; source failure and trial rationale documented |
| 10 | Repository cleanup | Inventory retained; cleanup deferred until end-to-end runs finish |
| 11 | Handoff | One center build comparison ready; no confirmed solver fix or production readiness |
| 12 | Preserve earlier work | All old results retained; no cancellation or deletion |

The [maintained checklist](docs/ORBITAL_RECON_DELIVERABLES.md) is the full scope.
The [general procedure](docs/ORBITAL_REVIEW_WORKFLOW.md) documents the controls.
This handoff uses **git pull**, not the older repository-root transfer tarball;
that historical tarball does not contain this retry. Do not extract it over the
checkout. Native acceptance requires state/orbital and interaction review,
not just successful controller completion.
