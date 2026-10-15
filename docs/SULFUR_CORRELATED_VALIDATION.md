# CH3S/CH3SS: correlated validation and completion of the examples

The [combined CASSCF reference panel](SULFUR_ORBITAL_REVIEW_2026-10-15.md)
is complete. The user requested all tests needed for high accuracy at minimum
cost, followed by completed CBH/ANL–ROTD–MESS examples. The next native batch
is a finite CASPT2-F12 comparison, not another reconnaissance expansion.
No smaller space or cheaper production basis is approved by its preparation.

Current milestone: the broad CASSCF orbital/initialization panel is complete;
the first correlated candidate comparison is partly complete and needs the
incomplete-point restart below. Sparse independent-method checks, final
inner/outer state-average choices, basis/cost validation and native production
are still outstanding. The available outputs do not yet justify selecting a
final production CAS for either example. Scheduler interruptions and iteration
limits are not evidence that every candidate is chemically unsuitable.

**Current next commands:** use
[the finite follow-up batches](SULFUR_REMAINING_CHECKS_2026-10-15.md).
They include the incomplete-point repair below, asymptotic-degeneracy
diagnostics, independent-method comparisons, and matched DZ/TZ tests. Do not
launch the older standalone repair in parallel with that combined follow-up.

## First native batch

The acceptance target is **accurate interaction potentials and spin gaps, with
stable state character and smooth, affordable ROTD sampling**. Neither
smoothness nor agreement with the largest CASSCF space proves accuracy.
Chemical completeness and correlated reference adequacy are separate checks.
A CAS(2,2) shortcut is not inferred from the presence of two radical fragments.

| Channel | Candidate CAS spaces | Larger diagnostic space | Radial stencil (Å) |
| --- | --- | --- | --- |
| CH3S C–S | (5,4), (9,6) | (13,8) | 2.9, 3.1, 3.3 |
| CH3SS C–S | (9,7), (13,9) | (19,12) | 3.1, 3.3, 3.5 |
| CH3SS S–S | (9,7), (13,9) | (19,12) | 3.5, 3.7, 3.9 |

There are **39 independent single-geometry jobs**, with a shared cap of 20:
each of the six smaller candidates gets three neighboring radial points, a
30 Å radial reference, and a 60° tilted point at the central distance. Each
larger diagnostic space gets only the center, far and tilted points initially.
This limits its cost while providing common reference comparisons and a
first orientation-sensitivity check. It is a finite accuracy/continuity/cost
screen, not complete coverage of either kinetic bottleneck.

All chemical hypotheses and seed choices are editable data in
[`validation_candidates.json`](../examples/anl/orbital_reviews/validation_candidates.json).
The preparation script reads the requested sampling method, basis, shift,
IPEA and scan grid from each linked KinBot input. It takes three neighboring
scan points near the reviewed intermediate-distance seed. For the present
inputs these resolve to CASPT2-F12/aug-cc-pVTZ-F12, shift 0.3 Eh and IPEA
0.25 Eh. Counts in preset method strings are replaced by each candidate's
counts. These shifts are method parameters, not convergence tolerances.

Each job reoptimizes a verified completed TZ-F12 natural-orbital record,
then calculates every root in the existing channel average: 3 doublets +
3 quartets, 1 + 1, and 6 + 6, respectively (**260 correlated root evaluations**).
This holds state averaging fixed across the initial space comparison; these
averages are diagnostic and are not assumed to be the final kinetic averages.
The stencil uses the correction scan's atom-pair coordinate, not a COM
surface coordinate or a predicted bottleneck location.

`prepare-partitions` accepts explicit source directories, target counts,
`close_lowest` and `discard_highest`. These refer to the source active natural
orbitals in their reviewed order. Closing orbitals makes them doubly occupied
in the candidate CI, **but their shapes still optimize in CASSCF**; this is
not Molpro's frozen-core approximation. The source/target electron and orbital
counts must agree. Each point requires an explicit verified WFU at the same basis. Optional
`target_distance_angstrom` moves the configured fragment along its existing
bond coordinate, retaining its declared orientation, and reoptimizes from that
seed. Source and target geometries remain distinct in provenance. There is no
implicit cross-basis transfer. Changing the target CAS
does not change or relabel the source energies. A changed state average is
an explicit separate model, and is supported by the generic partition plan.
Partitioned starts omit DIAB against the old, differently sized active space.

New optimizations use SO-SCI with explicit CASSCF energy target 1e-8 Eh and
the native gradient default. Exact CI is enabled only up to eight orbitals;
larger spaces use iterative CI. Global SCF/correlation thresholds retain their
source settings. The final reviewed selection can carry these validated
`casscf_options` into production, including the explicit energy tolerance;
selecting a WFU alone no longer has to lose its numerical controls. Native convergence and the requested active counts must be
confirmed in the output. The successful reference tolerance pilot supports
this starting choice, not a promise that every candidate will converge.

A correlated-recon bug was fixed before this launch: the old renderer removed
canonical-orbital output even when CASPT2 checks were requested. A reused WFU
could then supply stale canonical orbitals. The new MULTI step writes current
pseudo-canonical orbitals to 2141.2 and natural orbitals to 2140.2; each
CASPT2 check explicitly selects 2141.2, and each MRCI check selects 2140.2.
This follows Molpro's [CASPT2 orbital recommendation, including its IPEA caveat](https://www.molpro.net/manual/doku.php?id=multireference_rayleigh_schroedinger_perturbation_theory).
Explicit CI/RS2/MRCI nonconvergence also invalidates the result, even if the
program exits normally and prints energies. Shared production energy readers
now reject those failure messages as well.

## Native selector failure and recovery

The first batch exported at 13:22 on 15 October launched 20 of 39 jobs.
Eighteen finished CASSCF and all requested RS2-F12 calls (84 correlated root
energies). The generated `KBVRCSEL=energy(2)` read an undefined second element,
so Molpro assigned zero without energy units and KinBot rejected the results.
This was a shared implementation error, not evidence of failed CAS choices.
Two jobs (52907669 and 52907677) stopped mid-calculation; their batch stderr
reports an existing allocation retry journal. That is consistent with a Slurm
requeue, which retains its job ID. The export does not establish why Slurm
restarted those jobs. Nineteen jobs had not been submitted.

The old code incorrectly transferred the CCSD-F12 A/B convention to Molpro's
MR programs. **CASPT2-F12 and MRCI-F12 have no corresponding A/B choice here.**
For `STATE,1,root`, `ENERGY(1)` and `ENERGD(1)` select the sole calculated
state, including an excited root. `+Q` remains relaxed-reference Davidson
with `OPTION,CLUSTER=0`. The native method is now named accurately in contracts
and presets. Old MR `-f12b` names are explicitly recorded compatibility aliases;
`-f12a` requests are rejected. CCSD(T)-F12b is unchanged. See Molpro's
[RS2 state selection](https://www.molpro.net/manual/doku.php?id=multireference_rayleigh_schroedinger_perturbation_theory)
and [MRCI variables and F12 options](https://www.molpro.net/manual/doku.php?id=the_mrci_program).

Recovery validates each printed `!RSPT2 STATE` energy against the exact native
input and output program, electron count, active partition, multiplicity, root,
level shift and IPEA. Normal execution, convergence, input hashes, XML exports
and checkpoint checks remain required. It records recovery provenance and the
native method, leaving the original outputs and prepared inputs untouched.
Only this identifiable legacy selector defect permits recovery; generic
undefined variables still fail production validation. New executions of old
prepared decks correct the selector in the separately hashed native input.

The Slurm retry guard now distinguishes scheduler restart counts under one
allocation lock, archives interrupted products, verifies and restores original
seeds, and retains each restart's journal. Duplicate execution of the same
restart count is still rejected. This follows Slurm's
[requeue and SLURM_RESTART_COUNT semantics](https://slurm.schedmd.com/sbatch.html).
Collection now also includes batch scripts and the controller log.

### Preliminary scientific diagnostic from the recovered data

For CH3S at 30 Å, the three roots within each spin block span about
0.00015 kcal/mol at CASSCF, but **1.078, 1.177 and 1.176 kcal/mol** at CASPT2-F12
for CAS(5,4), (9,6) and (13,8), respectively. These are raw root-index spreads,
not a proven assignment of every correlated state to a particular fragment
term. They flag an unresolved asymptotic/state-treatment problem; increasing
CAS size alone has not removed it. The diagnostic panel remains useful, but
no candidate can be accepted solely on smooth curves or relative agreement.
The later MRCI+Q, root-character, rotational/term consistency and shift/state-
treatment checks must resolve this before a production surface is selected.
Do not hide the splitting by averaging or applying unvalidated offsets.

Recovery was replayed against all 18 native outputs (84 labelled energies),
and all 40 prepared/native input hashes from the 20 attempted jobs matched.
The local audit replay used explicit synthetic checkpoint placeholders because
binary WFU files are not in the export; Blodgett must verify the real files.
The full local suite passed 1,438 tests and 391 subtests; a final focused run
passed 126 tests after tightening recovery's successful-return-code check.
The only full-suite warnings were the existing deprecated `imp` import and
unwritable external Sella cache. No new native Molpro job ran locally.

## Second failure snapshot (14:48 on 15 October)

The later status paste and the subsequently supplied native excerpts establish
separate failure classes. The local full archive still predates this attempt.

- Pilots 10 and 21 have unfinished execution records after leaving the queue.
  The cause is not established by `status="running"` alone.
- Pilot 26's RS2-F12 step rejects its MCSCF orbital record as unconverged.
  MULTI was still improving at its 50-iteration limit: energy change
  1.6e-7 Eh and orbital gradient 2.6984e-4, against targets 1e-8 Eh and 1e-5.
  This is actual nonconvergence; `IGNORE_ERROR` would invalidate the comparison.
- Pilot 36 converged CASSCF in eight iterations and computed the first RS2-F12
  root. The second root exhausted 50 reference-CI iterations and explicitly
  requested more `MAXITI`. The new general `correlation_options` field raises
  that budget without loosening the requested accuracy.
- Several other pilots have no native output. The previous validator expanded
  this into many missing-orbital/root/occupation errors. Those derivative
  messages do not establish a chemistry failure. Audit now surfaces launch
  evidence, and execution records cover failures during seed/input staging.

The reference-CI failure phrase is now rejected by both reconnaissance and
production convergence validation, including when the process exits normally.
These changes do not establish the cause of the unfinished or output-free
jobs. Pilot 27 has no execution record, guard journal or batch streams. Slurm
accounting is disabled on this installation, so the scheduler cause remains
unknown. The restart below selects only incomplete points, retaining the
original verified seeds and chemistry. It raises CASSCF's ceiling to 150 and
reference CI's to 200; converged jobs still stop at their existing tolerances.
It does not claim these budgets guarantee native convergence. Keep completed
results and every original checkpoint.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite
KB_PY=/home/hram/KinBot/.venv/bin/python
RUN=/home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_20261015
"$KB_PY" -m kinbot.orbital_campaign stop-controller "$RUN"
"$KB_PY" -m kinbot.orbital_campaign collect "$RUN" \
    /home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_latest_20261015.tar.gz
BASH
```

Collection now also queries `sacct` for this campaign's recorded allocations
and prior attempts, including job states, exit codes/signals and restart counts.
The result is included as `scheduler_accounting.json`. Unavailable accounting
is recorded without preventing export of native evidence. See the
[Slurm accounting fields](https://slurm.schedmd.com/sacct.html).

Then, in a local Mac terminal:

```bash
bash <<'BASH'
set -euo pipefail
cd /Users/hrishikeshram/KinBot
scp hram@blodgett.ca.sandia.gov:/home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_latest_20261015.tar.gz \
    examples/anl/orbital_reviews/
BASH
```

## Restart incomplete points on Blodgett

Keep the existing campaign and all five reference/repair campaign directories;
the original binary WFUs are still required. The archive intentionally omits
WFUs, so the full checkpoint audit happens on Blodgett. The generic
`prepare-incomplete` command re-audits the finite single-geometry pilots and
copies only incomplete configurations into a fresh campaign. Completed points
remain in the original campaign; their paths/hashes are recorded in the new
plan. It preserves the original verified seeds and numerical controls except
the explicitly raised iteration ceilings. It never uses interrupted products
as starting guesses. No molecular names enter this selection logic.
All exits are inside a child shell and preserve the SSH session.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite
KB_PY=/home/hram/KinBot/.venv/bin/python
OLD=/home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_20261015
RUN=/home/hram/kinbot-runs/ch3s_ch3ss_potential_repair_20261015
test ! -e "$RUN"
"$KB_PY" -m pip install -e . --no-build-isolation --no-deps
"$KB_PY" -m kinbot.orbital_campaign stop-controller "$OLD"
"$KB_PY" -m kinbot.run_cleanup "$OLD" --stop --timeout 300
"$KB_PY" -m kinbot.orbital_campaign collect "$OLD" \
    /home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_latest_20261015.tar.gz
"$KB_PY" -m kinbot.orbital_review prepare-incomplete "$OLD" "$RUN" \
    --casscf-max-iterations 150 --reference-max-iterations 200
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
export KINBOT_MOLPRO_SCRATCH_MIN_MB=65536
mkdir -p "$KINBOT_MOLPRO_SCRATCH"
"$KB_PY" -m kinbot.orbital_campaign preflight "$RUN"
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$RUN" \
    --pilot-only --probe-method casscf --max-jobs 20 --max-refinements 0 \
    >> "$RUN/controller.log" 2>&1 < /dev/null &
printf '%s\n' "$!" > "$RUN/controller.pid"
"$KB_PY" -m kinbot.orbital_campaign progress "$RUN"
tail -n 20 "$RUN/controller.log"
squeue -u hram -o '%.18i %.10T %.48j %.10M %Z'
BASH
```

Monitor:

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python -m kinbot.orbital_campaign progress /home/hram/kinbot-runs/ch3s_ch3ss_potential_repair_20261015'
```

The expected successful status is `selected_pilot_complete`, with every newly
prepared pilot complete. Each system also has a duplicate full-panel manifest that
remains unrun. Do not drop `--pilot-only`: doing so repeats those geometries.
Elapsed times cover the reoptimization and all requested correlated roots;
they are not per-root timings. Slurm availability determines actual concurrency.

After the controller has stopped and its jobs have left Slurm, export even if
some comparisons failed; preserved native diagnostics will guide the next step:

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
RUN=/home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_20261015
STAMP=$(date +%Y%m%d_%H%M%S)
REPAIRED=/home/hram/kinbot-runs/ch3s_ch3ss_potential_repair_20261015
"$KB_PY" -m kinbot.orbital_review report "$REPAIRED" "$REPAIRED/review_report_$STAMP.json" \
    --additional-campaign "$RUN"
"$KB_PY" -m kinbot.orbital_campaign collect "$REPAIRED" \
    /home/hram/kinbot-runs/ch3s_ch3ss_potential_repair_20261015.tar.gz
BASH
```

`collect` includes the native evidence and the report just written above.
Copy the archive from a local Mac terminal:

```bash
bash <<'BASH'
set -euo pipefail
cd /Users/hrishikeshram/KinBot
scp hram@blodgett.ca.sandia.gov:/home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_latest_20261015.tar.gz \
    hram@blodgett.ca.sandia.gov:/home/hram/kinbot-runs/ch3s_ch3ss_potential_repair_20261015.tar.gz \
    examples/anl/orbital_reviews/
BASH
```

## Remaining decisions and acceptance checks

1. Compare candidate-minus-reference **interaction energies** and spin gaps,
   matching geometry and state character; total-energy offsets alone are not
   an active-space error estimate. Inspect intruders, reference weights, root
   order and occupations. The common center/far and tilted points test reference agreement and
   anisotropy; the radial stencils expose local branch changes and slope
   changes. None of these alone establishes complete production PES accuracy.
   When a sampling correction is proposed, judge its final corrected potential
   as well as the raw sampling error. A predominantly radial bias may be
   correctable; an unvalidated orientation-dependent residual is not removed
   merely by fitting the radial scan.
   `potential_profiles` groups single-point systems by model, basis, channel,
   region, orientation and state average, with labelled relative energies,
   secants and timings. Compare small/medium/large spaces at matching points,
   retaining root-character checks; do not rank them by total energies.
2. Extend spaces that pass both accuracy and continuity checks to the remaining
   relevant radial/tilted anchors,
   and test reduced inner state averages independently of the space change.
   Keep the fragment degeneracies needed for outer surfaces. A smaller outer
   CH3SS S–S CAS(7,5) is only worth an additional comparison if the measured
   savings justify it; the compressed high-state failure of its projection
   does not by itself decide an outer-only model.
3. Benchmark MRCI-F12+Q against CASPT2-F12 on the surviving spaces, including
   shift/IPEA sensitivity where indicated. Test ROHF-UCCSD(T)-F12b for the
   high-spin anchor and the production spin-splitting correction. Do not
   replace a multireference low-spin dissociation by a single-reference CC
   calculation. Preserve the intentional RHF/ROHF-UCC convention.
4. Compare DZ-F12 sampling with TZ-F12 corrections at common geometries and
   state averages. Generate basis-specific reference orbitals before changing
   the basis. Use a small QZ-F12 check if needed to resolve the TZ error;
   neither routine CBS nor a full QZ recon grid is the default. Check the
   orientation dependence left after the radial correction, not just the
   radial basis difference. Atomic spin–orbit data remain experimental.
5. Review and install the chosen channel/inner/outer spaces with
   `orbital_review install-selection`, together with explicit electronic
   weights/offsets and equilibration assumptions in the surface model. Preserve
   the CBH/ANL, ROTD and MESS controls in the example inputs. No vdW minimum
   is required to route an accepted `hom_sci` channel through two variational
   bottlenecks.
6. Run both native examples end to end, verifying CBH/ANL correction evidence,
   converged ROTD statistical sampling, inner/outer flux handling and MESS
   rates. As an initial stringent screen, investigate interaction/gap changes
   above 0.02 kcal/mol; final acceptance should require rate stability within
   about 2% across the requested temperature grid, with Monte Carlo uncertainty
   small enough to resolve that change. These are working numerical targets,
   not an absolute physical error claim. Basis economies in VRC do not by
   themselves justify removing independent ANL thermochemical corrections.

The native correlated comparisons must be completed and the final end-to-end
runs performed on Blodgett. Local preparation/parser tests cannot replace them.
