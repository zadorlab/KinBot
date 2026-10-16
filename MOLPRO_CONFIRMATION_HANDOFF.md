# Molpro confirmation: root-protected MRCI supplement

The new exports have been independently reviewed: **20 shared-reference
CASPT2 points and both CABS=3 asymptote gates pass**. The minimum active spaces
remain **CH3S C–S (5,4), CH3SS C–S (7,5), CH3SS S–S (7,5)**. The successful
original CH3SS C–S checks are reused.

The tilted CH3S MAXDAV20 MRCI retry converged numerically, but separate root
calls reached effectively duplicate states. Increasing the Davidson subspace
resolved stagnation without protecting distinct roots. The new input solves
all declared roots of each spin simultaneously. It retains relaxed Davidson,
the active spaces, geometries and tolerances. The
[native review and exact evidence](docs/MOLPRO_RECOVERY_REVIEW_2026-10-16.md)
explain the diagnosis, CABS results and acceptance conditions. **The new
simultaneous Molpro calculations have not yet been executed; full production
is not yet approved.**

## Start the seven-point supplement

The user handles SSH. On Blodgett:

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/orbital_reviews/molpro_confirmation/resolve_mrci_on_blodgett.sh
```

The source is the original stopped `_r2` run. Its inputs, outputs and completed
checkpoints are checked against the retained native review before preparation.
The destination is a separate sibling:
`/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate`.
Existing results and any independent running jobs are preserved. The script
cancels/deletes nothing and refuses to start a duplicate controller. Do not
rerun the historical full-restart, recovery or CABS-gate launch scripts.

The [explicit plan](examples/anl/orbital_reviews/molpro_confirmation/multistate_resolution_plan.json)
contains three independent chains:

| Chain | Points | Purpose |
|---|---|---|
| CH3S aTZ MRCI simultaneous, CABS=1 | tilted | Isolate root protection against the completed MAXDAV20 point |
| CH3S aTZ MRCI simultaneous, CABS=3 | far, center, tilted | Matched coupled-CABS interactions and tilted CABS control |
| CH3SS S–S aTZ MRCI simultaneous, CABS=3 | far, center, tilted | Complete the blocked reference checks with protected roots |

There are three or six roots per spin, solved together, with energy-ordered
`ENERGD(i)` selectors. These are correlated eigenstate ranks, not an assertion
that CASSCF root labels retain their character through crossings. The repeated
far points validate the new eigensolve and provide matched dissociation zeros.
No CASPT2, IPEA-series, expanded-CAS or geometry-optimization jobs are submitted.

Eight ranks use 1000 MW/rank; the reservation remains **104000 MB/node**, below
the reported 126000 MB node limit. Concurrency is three jobs. All known node
exclusions and bounded infrastructure retries are retained. Simultaneous roots
can cost more memory/time; native resource adequacy remains to be observed.
CASSCF MAXIT=999; correlation MAXIT/MAXITI=2000, MAXVI=60, MAXDAV=20,
THRDEN=1e-8, THRVAR=1e-10. CABS=3 is explicit in the two matched chains;
the general default is unchanged.

## Monitor and export

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py monitor /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate --roots'
```

The monitor now groups simultaneous roots into the correct same-spin native
call. During a solve, several energy markers can remain unprinted until the
whole call returns; that alone does not imply a stall. `--all --roots` also
lists completed points and waiting predecessors. These commands are read-only.

After the supplement stops:

```bash
/home/hram/KinBot/.venv/bin/python \
  examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py collect \
  /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate \
  --output /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_multistate.evidence.tar.gz
```

Transfer that archive to `examples/anl/orbital_reviews` for native review.
An existing archive is not overwritten. Earlier historical failures remain in
their original monitors; they are not new failures of this supplement.

## Deliverable status

| # | Deliverable | Current status |
|---|---|---|
| 1 | Literature-grounded minimum spaces and cost | Choices retained; native sampling timings now documented |
| 2 | Local ORCA chemical/radial/angular/return checks | Complete; 47 curated points retained |
| 3 | Native aDZ/aTZ Molpro transfer | CASPT2 recovery and CABS3 gates reviewed; MRCI simultaneous roots/interactions pending |
| 4 | General JSON controls | Explicit simultaneous MRCI option added; array selectors and gradient slot tested |
| 5 | General orbital reconnaissance | Procedure documented; chemistry remains in example data |
| 6 | Complete production inputs | CAS counts configured; accepted final stage models still to install |
| 7 | User-run end-to-end support | Remains in scope, including bug fixes through CBH-ANL/RotdPy/MESS completion |
| 8 | Generality/integration audit | 400 targeted tests pass, including all seven new transformed input paths; no local native Molpro |
| 9 | Reproducible evidence | Both new exports verified; native reviews, orbital/energy findings and analysis source retained |
| 10 | Repository cleanup | Inventory retained; cleanup deferred until end-to-end runs finish |
| 11 | Handoff | This MRCI-only supplement is ready to submit; production is not yet approved |
| 12 | Preserve earlier work | All old results retained; no cancellation or deletion |

The [maintained checklist](docs/ORBITAL_RECON_DELIVERABLES.md) is the full scope.
The [general procedure](docs/ORBITAL_REVIEW_WORKFLOW.md) describes the controls
for other systems. The scripts here are example adapters to that general code.

The optional repository-root `molpro_confirmation_transfer_20261016.tar.gz`
and SHA sidecar are refreshed with this change. Pulling is sufficient; do not
extract the snapshot over a running checkout. Full native exports remain in
`orbital_reviews` outside the source bundle; their extracted copies and per-file
hashes are identified in the review. Production readiness requires native state
and interaction review, not merely a successful controller status.
