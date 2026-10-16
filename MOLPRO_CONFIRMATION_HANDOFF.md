# Molpro confirmation: targeted S–S MRCI solver retry

The minimum spaces remain **CH3S C–S (5,4), CH3SS C–S (7,5), CH3SS S–S
(7,5)**. The CASPT2 confirmation and four simultaneous-root CH3S MRCI points
pass native review. The S–S far point has now completed according to the
controller. Center failed at the first simultaneous MRCI iteration with
`dsygv failed`; tilted is waiting for center's checkpoint. The stopped
seven-point run has **five validated, one failed, one waiting**.

The next submission is a **controlled numerical trial**, not a guaranteed
repair: center → tilted with Davidson denominator shifts of 0.1 Eh for pairs,
singles and internals. All six roots per spin, CABS=3, SFIX, tolerances and
relaxed-Davidson selectors are retained. The completed far checkpoint/energy
is reused. The [native diagnosis and documentation review](docs/MOLPRO_MRCI_DSYGV_RECOVERY_2026-10-16.md)
explain the trial, the missing LAPACK INFO code and why imaginary-shift advice
for perturbative CCD/DCD F12 corrections is a different context.

## Start the two-point trial

The user handles SSH. On Blodgett:

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/orbital_reviews/molpro_confirmation/retry_mrci_solver_on_blodgett.sh
```

Source:
`/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate`.
Destination:
`/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_solver_retry`.
The script cancels/deletes nothing, leaves other campaigns alone and refuses a
live source or duplicate destination controller. Do not run the historical
full-restart, recovery, CABS-gate or original multistate launch scripts again.

Before preparation it rechecks native completion and checkpoint consistency
on temporary manifest copies, verifies the diagnosed center error, then binds
the source input/output/WFU hashes into the new run's provenance. All five
completed source results remain intact. The selected far WFU is copied without
reapplying its former full-valence orbital permutation. The retained
`source_completion_snapshot.json` is a completion check, not approval of the
far state spectrum. Its full export still needs independent local review.

The [plan JSON](examples/anl/orbital_reviews/molpro_confirmation/mrci_preconditioner_plan.json)
changes only `correlation_options.denominator_shifts`. It renders
`shiftp=0.1,shifts=0.1,shifti=0.1` on each MRCI command. The numerical target
is the same MRCI eigenproblem, so this does not require repeating the completed
far calculation. A failed center leaves tilted unsubmitted; it cannot be
silently accepted or bypassed.

Resources remain **eight ranks, 1000 MW/rank and 104000 MB/node**, below the
reported 126000 MB node limit. Concurrency is one. Existing exclusions and
bounded infrastructure retries are inherited. CASSCF MAXIT=999; correlation
MAXIT/MAXITI=2000, MAXVI=60, input MAXDAV=20, THRDEN=1e-8, THRVAR=1e-10.
The failed six-root native call printed an effective MAXDAV=24. No default
settings for unrelated KinBot calculations are changed.

## Monitor and export

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py monitor /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_solver_retry --all --roots'
```

Initially expect two points: center active/queued and tilted waiting. Each
point has two simultaneous calls (doublet and quartet), twelve selected
energies in total. Passing iteration 1 is progress, not final acceptance.

Export the stopped source now to retain the completed S–S far calculation and
failed center context. This is independent of the new retry:

```bash
cd /home/hram/KinBot
.venv/bin/python examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py collect \
  /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate \
  --output /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate.evidence.tar.gz
```

After the retry stops, export it with:

```bash
cd /home/hram/KinBot
.venv/bin/python examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py collect \
  /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_solver_retry \
  --output /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate_solver_retry.evidence.tar.gz
```

The user transfers both archives to `examples/anl/orbital_reviews`. Collection
refuses to overwrite an existing archive or collect a live/uncertain run.
Original historical failures remain in their monitors and are not new failures
of this trial.

## Deliverable status

| # | Deliverable | Current status |
|---|---|---|
| 1 | Literature-grounded minimum spaces and cost | Choices retained; native sampling timings documented |
| 2 | Local ORCA chemical/radial/angular/return checks | Complete; 47 curated points retained |
| 3 | Native aDZ/aTZ Molpro transfer | CASPT2 and four CH3S MRCI points reviewed; S–S far completed per controller, center retry/tilted pending |
| 4 | General JSON controls | Simultaneous MRCI confirmed for CH3S; optional denominator controls added for the explicit numerical trial |
| 5 | General orbital reconnaissance | Procedure documented; chemistry remains in example data |
| 6 | Complete production inputs | CAS counts configured; final accepted stage models still to install |
| 7 | User-run end-to-end support | Remains in scope through CBH-ANL/RotdPy/MESS completion |
| 8 | Generality/integration audit | Local rendering/provenance tests pass; native trial remains untested locally |
| 9 | Reproducible evidence | Existing exports and CH3S audit retained; source failure and trial rationale documented |
| 10 | Repository cleanup | Inventory retained; cleanup deferred until end-to-end runs finish |
| 11 | Handoff | Two-point solver trial ready; production readiness remains pending |
| 12 | Preserve earlier work | All old results retained; no cancellation or deletion |

The [maintained checklist](docs/ORBITAL_RECON_DELIVERABLES.md) is the full scope.
The [general procedure](docs/ORBITAL_REVIEW_WORKFLOW.md) documents the controls.
This handoff uses **git pull**, not the older repository-root transfer tarball;
that historical tarball does not contain this retry. Do not extract it over the
checkout. Native acceptance requires state/orbital and interaction review,
not just successful controller completion.
