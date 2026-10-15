# CH3S/CH3SS: correlated validation and completion of the examples

The [combined CASSCF reference panel](SULFUR_ORBITAL_REVIEW_2026-10-15.md)
is complete. The user requested all tests needed for high accuracy at minimum
cost, followed by completed CBH/ANL–ROTD–MESS examples. The next native batch
is a finite CASPT2-F12b comparison, not another reconnaissance expansion.
No smaller space or cheaper production basis is approved by its preparation.

## First native batch

The acceptance target is **accurate interaction potentials and spin gaps, with
stable state character and smooth, affordable ROTD sampling**. Neither
smoothness nor agreement with the largest CASSCF space proves accuracy.
Chemical completeness and correlated reference adequacy are separate checks.
A CAS(2,2) shortcut is not inferred from the presence of two radical fragments.

| Channel | Candidate CAS spaces | Larger diagnostic space | Radial stencil (Å) |
| --- | --- | --- | --- |
| CH3S C–S | (5,4), (9,6) | (13,8) | 3.1, 3.3, 3.5 |
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
inputs these resolve to CASPT2-F12b/aug-cc-pVTZ-F12, shift 0.3 Eh and IPEA
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

## Launch on Blodgett

Keep the five reference/repair campaigns at their current paths: their binary
WFUs are the inputs. The following block runs the 39 new pilot calculations, at most 20
concurrently. It does not restart the old failed attempts. Rerunning the block
stops only this new controller and adopts its live allocations and completed
jobs. It never automatically retries a scientifically failed calculation.
`set -e` is confined to the child shell, so failure preserves the SSH session.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite
KB_PY=/home/hram/KinBot/.venv/bin/python
RUN=/home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_20261015
PLAN=/home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_20261015.plan.json
"$KB_PY" -m pip install -e . --no-build-isolation --no-deps
if [ ! -f "$RUN/campaign.json" ]; then
    if [ ! -f "$PLAN" ]; then
        "$KB_PY" examples/anl/orbital_reviews/prepare_validation.py "$PLAN"
    fi
    "$KB_PY" -m kinbot.orbital_review prepare-partitions "$PLAN" "$RUN"
fi
"$KB_PY" -m kinbot.orbital_campaign stop-controller "$RUN"
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
watch -n 10 '/home/hram/KinBot/.venv/bin/python -m kinbot.orbital_campaign progress /home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_20261015'
```

The expected successful status is `selected_pilot_complete`, with **39/39
pilots complete**. Each system also has a duplicate full-panel manifest that
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
"$KB_PY" -m kinbot.orbital_review report "$RUN" "$RUN/review_report_$STAMP.json"
"$KB_PY" -m kinbot.orbital_campaign collect "$RUN" \
    /home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_20261015.tar.gz
BASH
```

`collect` includes the native evidence and the report just written above.
Copy the archive from a local Mac terminal:

```bash
bash <<'BASH'
set -euo pipefail
cd /Users/hrishikeshram/KinBot
scp hram@blodgett.ca.sandia.gov:/home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_20261015.tar.gz \
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
3. Benchmark MRCI-F12b+Q against CASPT2-F12b on the surviving spaces, including
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

The native correlated comparisons and the final end-to-end runs remain to be
performed on Blodgett. Local preparation/parser tests cannot replace them.
