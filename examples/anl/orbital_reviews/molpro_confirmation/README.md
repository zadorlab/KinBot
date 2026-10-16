# Minimal native transfer confirmation

**Current result:** 20 shared-reference CASPT2 recovery points and both CABS=3
asymptote gates pass independent native review. The tilted MAXDAV20 MRCI retry
completed numerically but contains duplicate-root evidence. The new seven-point
MRCI-only supplement uses simultaneous same-spin roots and matched CABS=3
far/center/tilted energies. Start with the
[current handoff](../../../../MOLPRO_CONFIRMATION_HANDOFF.md) and
[native review](../../../../docs/MOLPRO_RECOVERY_REVIEW_2026-10-16.md).
Run `resolve_mrci_on_blodgett.sh`; do not relaunch historical recovery/full
restart/gate scripts. All previous results are preserved. The selected spaces
remain (5,4), (7,5), (7,5); full KinBot production is still pending final native
reference acceptance and model installation.

`run_confirmation.py monitor RUN --roots` reports correlated calls and emitted
energy markers, including grouped simultaneous roots; `--all` adds successful
points and waiting predecessors. These read-only displays do not substitute
for scientific review. The older calculation definitions below remain a
reproducibility record; their independent-root/CABS policies are superseded
where the new plan explicitly says so.

This is a fixed-geometry backend check. It does not repeat optimizations or the
production relaxed/frozen correction scans. Nuclear coordinates are the exact
local comparison coordinates, including the C3v/Cs reference cuts and the
deliberately tilted checks. Electronic calculations use C1, retaining the full
declared state manifold. The input records identify the source configuration
and its hash for every geometry.

## Calculations and conventions

| Calculation | Points per channel | Basis | Native jobs |
| --- | --- | --- | --- |
| CASPT2-F12 | 30 A, interaction center, bonded-side point, exact return to center, tilted center | aDZ-F12 and aTZ-F12 | 30 |
| MRCI-F12+Q | 30 A, center, tilted center | aTZ-F12 | 9 |

There are nine independent orbital-propagation chains and **39 jobs total**.
Every requested spin/root is computed at each point. The first point gates its
dependent points; there is no duplicate far-point pilot. All chemistry and site
settings live in `geometry_and_seed_catalog.json`; the scripts are adapters to
the general KinBot preparation, native executor, retry guard and monitor.

The spin-free CASSCF averages use three roots per spin for CH3S C–S, one per
spin for CH3SS C–S, and six per spin for CH3SS S–S. Doublet/quartet block weights
are equal. The latter six arise from two CH3S electronic components times three
S(3P) components, not six degenerate sulfur states. These optimization weights
are separate from physical VRC statistical weights.

CASPT2 uses real shift 0.3 Eh and IPEA 0.25, following the user's requested
baseline. Existing local ORCA screening used IPEA zero; its matched sensitivity
checks were cancelled at the user's request. Do not attribute an ORCA/Molpro difference solely to F12
while these conventions differ. These are explicit proposed example settings,
not inherited preset defaults.
MRCI reads the natural orbitals and uses the relaxed-reference Davidson energy
selector `ENERGD(1)` with explicit `OPTION,CLUSTER=0`. Molpro defines this as
`ENERGD0(1)`; fixed and rotated alternatives are `ENERGD1(1)` and `ENERGD2(1)`.
The assigned energy is read **after** the full MRCI command returns, including
CABS singles. Do not parse the first/last printed cluster-correction line or
accept a value printed before a later fatal error. See the
[Molpro energy-variable definitions](https://www.molpro.net/manual/doku.php?id=the_mrci_program#cluster_corrections_for_multi-state_mrci).

The original catalogue retains **CABS_SINGLES=1**, the existing singles-CI default.
The old (9,7) S–S far-point job converged its MRCI equations, printed Davidson
values, then failed with `Orthogonalization problem` in `mrf12_singles.F90`.
It is not a completed result. This error is not cured by increasing iteration
limits. The documented perturbational alternative is 2 (3 is the coupled
variant). The minimum-space `_r2` run reproduced the error. The recovery
trial of 2 then failed in both asymptotes with `strange..` in the same routine;
matching reference counts and converged MRCI were confirmed in both. The
separate `cabs3_gate_plan.json` now tests 3 at those two geometries only, for
all requested roots, before any further scans. Both gates have now passed; matched interacting energies and simultaneous
root protection are the remaining checks.
The earlier decision to leave the default alone has been superseded by this
native evidence and the user's request for a resolution. The official minimum-space far-point
gates must establish actual success. If they reproduce this failure, retain
the results and test a declared alternative before installing production
models. The old failure used an expanded CAS and multiple reference states;
the repaired recipe uses minimum spaces and its explicit full-reference policy, so
failure or success cannot be assumed from that old result.
[Molpro CABS options](https://www.molpro.net/manual/doku.php?id=the_mrci_program#explicitly_correlated_mrci-f12). All roots are retained in the report. A fixed selected-root
placeholder in the catalogue is not proof that this root remains the lowest
correlated state at every geometry.
`settings.reference_policies` explicitly sets the correlation reference treatment
for each method; it is separate from the CASSCF state average. The current
original recipe computes each correlated root separately. Original CASPT2 uses `selected_root`;
MRCI-F12 uses `state_average`, with all declared same-spin references and an
explicit matching NSTATI. Lone excited-root references triggered an installed
F12 solver restriction and are now rejected during generation. Review duplicate or exchanged
root character in the native results; neither a root number nor `REFSTATE`
alone establishes root tracking. The reference policy and sparse MRCI point
list are editable catalogue data, as are the SCF limits.

The supplemental `recovery_plan.json` changes CASPT2 to `state_average` for
the CH3S and S–S manifolds while keeping IPEA 0.25/shift 0.3. Existing completed
controls show that this removes the ~1.175 kcal/mol CH3S artifact. The change
is rechecked on the five current geometries in both bases. C–S in CH3SS has
one reference per spin and needs no repeat. The original catalogue and old
energies remain historical evidence, not the newly selected reference policy.

The original reviewed aTZ-F12 natural-orbital records are hash-bound. Initial
`START,...,OVL` projection and reoptimization produce the target-geometry model
in each requested basis. **A projected guess is not an approved aDZ-F12 orbital
model.** Review its newly optimized native orbitals and asymptote before using
that destination-basis WFU in production. Source-active swaps occur once, before
closing/discarding orbitals; subsequent points propagate the target CAS without
repeating the source swaps. Both absolute source Molpro indices and their
conversion to source-active indices are recorded.

Foreign-basis seeds are staged on FILE,3; the native input writes a fresh
FILE,2, avoiding the source/destination dimension collision. Each initial
or propagated RHF explicitly sets spin and determinant occupations. See the
[feature-by-feature manual audit and native failure record](../../../../docs/MOLPRO_CONFIRMATION_INPUT_AUDIT_2026-10-16.md).

To replace the failed confirmation while retaining all its results, run:

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/orbital_reviews/molpro_confirmation/restart_on_blodgett.sh
```

The replacement is `minimum_active_space_molpro_confirmation_20261016_r2`.
This script is for the confirmation itself; do not reuse the original
`launch_on_blodgett.sh` to restart it. The old directory is preserved without
deleting evidence. All 39 points are rerun consistently under the corrected
protocol; this is separate from later production KinBot calculations.

Numerical settings are SCF E=1e-8 Eh/maxit=2000; CASSCF E=1e-8 Eh,
gradient=1e-5/maxit=999; correlation E=1e-8 Eh/coefficient=1e-5,
maxit=maxiti=2000, reference-vector capacity=60. The coefficient criterion
renders Molpro THRVAR=1e-10, not 1e-20. Eight ranks with 1000 MW each request
104000 MB per node. This fits the reported 126000 MB node limit and leaves
40000 MB above the nominal 64000 MB rank stacks; actual peak use must still be
observed. Current exclusions and bounded startup/I/O/allocation retries are
editable catalogue data. Inspect current node health before a future launch.

## Prepare and run on Blodgett

The user transfers the updated repository/scripts and handles SSH. The three
catalogued source WFU paths must exist on Blodgett with their recorded hashes.
Preparation verifies all three before creating the run directory and submits
nothing. Use a **new directory**, preserving the earlier studies.

First launch from the updated `composite` branch (the user handles SSH):

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/orbital_reviews/molpro_confirmation/launch_on_blodgett.sh
```

The repository-root `molpro_confirmation_transfer_20261016.tar.gz` is a
checksummed transfer snapshot. After pulling, use the checked-in scripts
directly; no archive extraction or reinstall is needed for this checkout's
existing virtual environment. See [the root handoff record](../../../../MOLPRO_CONFIRMATION_HANDOFF.md).

The script prepares and validates the new decks **before** stopping old work.
It then writes a dated preservation report beside the old campaign, stops its
controller and only jobs with working directories in that old tree, and audits
the retained files. The old native outputs/WFU/retry history remain in place.
It checks the three required seed hashes and scheduler access before sending
signals. A scheduler outage, timeout, or possible matching completed point
stops the launch for review; it does not assume an empty queue or waste a valid
result. The report is `OLD_ROOT.handoff.TIMESTAMP/evidence.json`.

The most recent supplied snapshot had 59 validated, two failed and seven live
expanded-space S–S MRCI jobs. Those seven are superseded by the minimum-space
checks. Live IDs are discovered again at execution, not hardcoded from the
snapshot. No other campaign or user's allocation is cancelled. The full-manifest
68 waiting entries belonged to the old pilot-only study and are not another
68 required jobs.

The default new directory is
`/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016`.
The launcher refuses an existing new directory before stopping old work.
Override either directory with its first/second argument (old/new). If a
prepared launch stopped during transition, finish/review that transition and
then use the explicit resume command below; do not regenerate inputs.

```bash
KB_PY=/home/hram/KinBot/.venv/bin/python
KB_CONFIRM=/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016
KB_TOOLS=/home/hram/KinBot/examples/anl/orbital_reviews/molpro_confirmation
"$KB_PY" "$KB_TOOLS/run_confirmation.py" monitor "$KB_CONFIRM"
# Resume an interrupted confirmation controller only:
nohup "$KB_PY" -u "$KB_TOOLS/run_confirmation.py" run "$KB_CONFIRM" \
  >> "$KB_CONFIRM/controller.log" 2>&1 < /dev/null &
```

The monitor includes actual Slurm exclusions, native stages/iterations, file
ages, retry counts and current error text. Repeat its command under `watch -n
10` if desired. Resume an interrupted controller by rerunning the same `run`
command; a lock prevents two controllers. The executor validates preserved
completed jobs and waits for existing allocations rather than resubmitting
them. Killing the controller alone preserves native jobs. Terminal scientific
failures require diagnosis; this adapter does not relabel them as convergence.

After the controller and its native jobs have stopped, export the evidence:

```bash
"$KB_PY" "$KB_TOOLS/run_confirmation.py" collect "$KB_CONFIRM" \
  --output "${KB_CONFIRM}.evidence.tar.gz"
```

Collection re-audits the native files, requires a successful scheduler query
with no live allocations for this run, and refuses to overwrite an archive.
It includes original reviewed seeds and resulting WFU files, so it may be
large. `export_inventory.json` records every archived file's SHA-256; the
archive has its own `.sha256` sidecar. An incomplete calculation can be exported
for diagnosis and is still labelled incomplete. The user transfers this archive
for review; the scripts do not use SSH.

## Required review and evidence retention

The completed manifests contain state-labelled CASSCF and correlated energies,
native convergence checks, occupations, checkpoint hashes and orbital XML paths.
Retain `catalog.json`, `confirmation.json`, all input snapshots, manifests,
`execution.json`, native `.inp`/`.out`, XML, WFU, launcher logs and retry histories.
Record Molpro version and the loaded site setup with the result archive.

Review the following before installing basis-specific selections:

1. All requested roots converge, correct frozen-core counts are printed
   (6 for CH3S, 11 for CH3SS), and native method/shift/energy selectors agree.
2. The intended orbitals remain active in each basis; no core or unrelated
   valence orbital replaces them. Rotations within the same active span are
   allowed. Inspect projected starting guesses as well as optimized outputs.
3. Each spin-free asymptote has the appropriate fragment manifold; preserve
   genuine CH3S fragment splitting. Use one common dissociation zero per
   declared model rather than separately zeroing every excited root.
4. Return energies and orbital character reproduce the forward calculation;
   compare relative radial/angular energies and spin gaps with the retained
   ORCA data. Identify correlated root reordering explicitly. Sparse points
   do not establish smoothness at every later Monte Carlo geometry.
5. Compare the sparse MRCI-F12+Q and CASPT2-F12 spin gaps and review the physical
   surface assignment before enabling the averaged-gap correction. Do not
   mistake equal orbital-average weights for surface populations.

The selected WFU is recorded with `kinbot.molpro_reconnaissance select` and
installed with `kinbot.orbital_review install-selection`. The physical surface
model requires its separate state/weight review. Native
acceptance remains a reviewed decision; merely finishing these jobs does not
certify production kinetics or guarantee error-free cluster execution.
