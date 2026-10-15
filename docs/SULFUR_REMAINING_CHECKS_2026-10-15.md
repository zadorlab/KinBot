# Finite CH3S/CH3SS follow-up: asymptotes, references and sampling cost

This replaces the standalone incomplete-point repair commands in
[the first correlated-panel notes](SULFUR_CORRELATED_VALIDATION.md). The broad
CASSCF reference panel is complete. Final CAS choices and native
CBH/ANL–ROTD–MESS production remain unvalidated.

## The 30 Å splitting

Same total spin does not imply equal energy: roots can correlate to different
fragment electronic terms. However, the three intended CH3 + S(3P) components
within each doublet/quartet manifold should approach degeneracy in these
spin-free calculations. There is no spin–orbit operator in these decks.

In the archived CH3S CAS(5,4) far pilot (`02_pilot/avt_f12__cs_r_30`), the native
doublet CASSCF energies are -437.068974957096, -437.068974955967 and
-437.068974720089 Eh: a span of about 0.000149 kcal/mol. Native RS2-F12 gives
-437.394347860051, -437.394347860791 and -437.396065040407 Eh: about
1.078 kcal/mol. These are native state energies, not the old marker-selector
bug. The larger CAS(13,8) also has a gap of approximately 1.18 kcal/mol.
That is not credible evidence of a physical long-range splitting.

The leading hypothesis is an imbalance in the approximate correlation
treatment of a degenerate manifold. Molpro's `STATE,1,root` uses a state-specific
density/Fock operator; `REFSTATE` can instead define a shared reference-averaged
Fock operator. Molpro also documents that its IPEA treatment is not exactly
orbital invariant. Pseudocanonical orbitals are already being used. These
facts motivate the controls below; they do not establish which contribution
caused this particular gap. C1 symmetry permits arbitrary representations of
the degenerate manifold; it is not a physical source of splitting.
See the [Molpro RS2 manual](https://www.molpro.net/manual/doku.php?id=multireference_rayleigh_schroedinger_perturbation_theory).

Do not label this gap experimental spin–orbit splitting, or erase it by
independent root offsets before understanding it. CH3SS S–S also requires
assignment of the CH3S fragment states: its six roots per spin need not all
represent the same fragment electronic term.

## Prepared work

All chemistry, spaces, fragment controls and source choices are example data
in [validation_candidates.json](../examples/anl/orbital_reviews/validation_candidates.json)
and [followup_recipe.json](../examples/anl/orbital_reviews/followup_recipe.json).
[prepare_followup.py](../examples/anl/orbital_reviews/prepare_followup.py) composes
the general campaign and partition tools. No sulfur names, paths or CAS counts
are added to shared KinBot runtime logic.

**Batch A: 47 new independent pilots plus only unfinished old pilots.**

* 24 common-geometry comparisons: center/far for every small, medium and large
  candidate, plus tilted points for the six smaller candidates. Compare
  CASPT2-F12 with IPEA=0, both selected-root and shared-reference treatments,
  against MRCI-F12+Q with the complete same-spin reference manifold.
* At the nine far points, additionally compare the original IPEA=0.25,
  shared-reference IPEA=0.25, and a smaller real shift (0.1 rather than 0.3 Eh,
  with IPEA=0). These are diagnostic Hamiltonian choices, not tolerance changes.
* Six separation/orientation controls: the smallest candidate in each channel
  at 60 Å and at a tilted 30 Å geometry, with the five CASPT2 treatments.
* Two isolated sulfur controls: CAS(4,3) targeting 3p, and CAS(6,4) including
  3s, averaging all three triplet components, with CASPT2 variants and MRCI+Q.
* Fifteen fresh AVAS/CASSCF DZ-F12 references: center, left, right, far and
  tilted for each channel. They use no TZ checkpoint.

The original successful pilots stay in their existing directories. Unfinished
ones reuse their original verified seeds in the new campaign, with up to
150 CASSCF macroiterations and 200 internal/reference-CI iterations. New MR
comparisons also allow 150 correlated iterations. Failed native convergence
remains a failure; `IGNORE_ERROR` is never used. Accounting being disabled does
not block preparation or export, and cannot establish why an old job vanished.

**Batch B: 87 independent pilots, after Batch A.**

* 39 matched DZ-F12 candidate points, with original IPEA=0.25 and both
  selected-root/shared-reference IPEA=0 CASPT2 comparisons.
* 24 TZ-F12 points completing the zero-IPEA comparisons and supplying matched
  high-spin anchors. Fifteen TZ comparisons already supplied by Batch A are
  not repeated. Where a new CC anchor is needed, its short CASSCF prelude and
  CASPT2 checks are repeated in the same job to retain a self-contained record.
* A restricted-reference UCCSD(T)-F12b high-spin anchor at all five small-space
  points per channel, in both bases. This supplies one ground state of the
  requested spin, not the full excited manifold. Inspect its reference
  character; it is not an unconditional exact benchmark.
* 24 reduced-inner-average checks: the two smaller spaces of CH3S C–S and
  CH3SS S–S, with one or two roots per spin at center/far/tilted geometries.
  Compare CASPT2 and MRCI+Q. These averages are inner-region candidates;
  incomplete asymptotic manifolds must not be silently installed for the outer
  bottleneck. CH3SS C–S already uses one root per spin in this comparison.

Both batches use a **single campaign-wide cap of 20 allocations** and no
adaptive expansion. `--pilot-only` is essential: these single-geometry
systems have unused duplicate full-panel manifests. Successful completion is
`selected_pilot_complete`, even though those duplicates remain unrun.
Method checks within a job reuse its CASSCF orbitals; job elapsed time includes
all those methods/roots. Use the native program timings to estimate sampling
cost rather than interpreting total job time as one potential evaluation.

## Run Batch A now on Blodgett

This stops only the old validation campaign and preserves its evidence and
checkpoint sources. Do not remove old reference directories. Failures exit
the child shell, leaving the SSH connection open. Use a fresh output path;
the block refuses to overwrite an existing follow-up.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite
KB_PY=/home/hram/KinBot/.venv/bin/python
OLD=/home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_20261015
RUN=/home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015
test ! -e "$RUN"
"$KB_PY" -m pip install -e . --no-build-isolation --no-deps
"$KB_PY" -m kinbot.orbital_campaign stop-controller "$OLD"
"$KB_PY" -m kinbot.run_cleanup "$OLD" --stop --timeout 300
"$KB_PY" -m kinbot.orbital_campaign collect "$OLD" \
    /home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_latest_20261015.tar.gz
"$KB_PY" examples/anl/orbital_reviews/prepare_followup.py references "$OLD" "$RUN"
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
export KINBOT_MOLPRO_SCRATCH_MIN_MB=65536
mkdir -p "$KINBOT_MOLPRO_SCRATCH"
"$KB_PY" -m kinbot.orbital_campaign preflight "$RUN"
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$RUN" \
    --pilot-only --probe-method casscf --max-jobs 20 --max-refinements 0 \
    >> "$RUN/controller.log" 2>&1 < /dev/null &
printf '%s\n' "$!" > "$RUN/controller.pid"
"$KB_PY" -m kinbot.orbital_campaign progress "$RUN"
BASH
```

Monitor in another terminal (Ctrl-C stops only `watch`):

```bash
cd /home/hram/KinBot
watch -n 10 '/home/hram/KinBot/.venv/bin/python -m kinbot.orbital_campaign progress /home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015; squeue -u hram -o "%.18i %.10T %.48j %.10M %R"'
```

## Run Batch B after A completes

The status check prevents a second 20-job campaign from starting alongside A.
Preparation verifies each completed DZ reference and uses only its own basis
checkpoint. These are comparison calculations, not an approval of IPEA=0 or
of any production CAS. If A reveals a systematic failure of both zero-IPEA
treatments, review that evidence before spending the second batch.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
REF=/home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015
RUN=/home/hram/kinbot-runs/ch3s_ch3ss_sampling_checks_20261015
test ! -e "$RUN"
"$KB_PY" - "$REF" <<'PY'
import json, sys
from pathlib import Path
c = json.loads((Path(sys.argv[1]) / 'campaign.json').read_text())
if c['status'] != 'selected_pilot_complete':
    raise SystemExit('Batch A has not completed; inspect/export it before starting Batch B.')
PY
"$KB_PY" examples/anl/orbital_reviews/prepare_followup.py potentials "$REF" "$RUN"
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
export KINBOT_MOLPRO_SCRATCH_MIN_MB=65536
mkdir -p "$KINBOT_MOLPRO_SCRATCH"
"$KB_PY" -m kinbot.orbital_campaign preflight "$RUN"
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$RUN" \
    --pilot-only --probe-method casscf --max-jobs 20 --max-refinements 0 \
    >> "$RUN/controller.log" 2>&1 < /dev/null &
printf '%s\n' "$!" > "$RUN/controller.pid"
"$KB_PY" -m kinbot.orbital_campaign progress "$RUN"
BASH
```

```bash
cd /home/hram/KinBot
watch -n 10 '/home/hram/KinBot/.venv/bin/python -m kinbot.orbital_campaign progress /home/hram/kinbot-runs/ch3s_ch3ss_sampling_checks_20261015; squeue -u hram -o "%.18i %.10T %.48j %.10M %R"'
```

## Export after jobs and controllers stop

This also works after A alone; B is included only if it exists. Native failures
remain in the archive. Disabled `sacct` is recorded, not treated as an empty
history or an export failure. The exports omit large WFUs, which stay on Blodgett.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
OLD=/home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_20261015
REF=/home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015
POT=/home/hram/kinbot-runs/ch3s_ch3ss_sampling_checks_20261015
STAMP=$(date +%Y%m%d_%H%M%S)
"$KB_PY" -m kinbot.orbital_review report "$REF" "$REF/review_report_$STAMP.json" \
    --additional-campaign "$OLD"
"$KB_PY" -m kinbot.orbital_campaign collect "$REF" \
    /home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015.tar.gz
if [ -f "$POT/campaign.json" ]; then
    "$KB_PY" -m kinbot.orbital_review report "$POT" "$POT/review_report_$STAMP.json" \
        --additional-campaign "$REF" --additional-campaign "$OLD"
    "$KB_PY" -m kinbot.orbital_campaign collect "$POT" \
        /home/hram/kinbot-runs/ch3s_ch3ss_sampling_checks_20261015.tar.gz
fi
BASH
```

Copy A and the latest original results from a **local Mac terminal**:

```bash
bash <<'BASH'
set -euo pipefail
cd /Users/hrishikeshram/KinBot
scp hram@blodgett.ca.sandia.gov:/home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_latest_20261015.tar.gz \
    hram@blodgett.ca.sandia.gov:/home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015.tar.gz \
    examples/anl/orbital_reviews/
BASH
```

After B's export, copy it similarly:

```bash
bash <<'BASH'
set -euo pipefail
cd /Users/hrishikeshram/KinBot
scp hram@blodgett.ca.sandia.gov:/home/hram/kinbot-runs/ch3s_ch3ss_sampling_checks_20261015.tar.gz \
    examples/anl/orbital_reviews/
BASH
```

## What these batches decide, and what still depends on the results

1. Identify whether the asymptotic gap changes with distance, orientation,
   IPEA/reference policy, isolated-fragment treatment or correlation method.
   Do not infer that every same-spin root must be degenerate without fragment
   term assignments. Shared-reference CASPT2 is not XMS-CASPT2.
2. Compare interaction energies and spin gaps at matching states/geometries,
   including angular residuals after a radial correction. A larger CAS or a
   smoother curve alone cannot establish accuracy. Inspect retained orbitals
   and CI character after reoptimization rather than equating ranks to chemistry.
3. Use the stencils, changes in state averaging, basis comparisons and native
   timings to choose promising sampling/correction combinations. The combined
   report retains all method/reference labels. Profiles do not silently merge
   different campaigns or state averages; comparison of A/B points must retain
   their explicit provenance and match the full method labels.
4. Only surviving combinations warrant additional coverage of inner/outer
   anchors and, if the measured DZ/TZ difference or reference disagreement
   requires it, sparse QZ-F12 checks. There is no routine CBS or full QZ scan.
   A finite preselected batch cannot guarantee that no further chemistry
   checks will be needed; these results decide their necessity and scope.
5. The user then approves channel/stage-specific spaces and methods. Install
   the reviewed selections and run the two native end-to-end examples,
   checking ROTD statistical convergence and sensitivity of the MESS rates.
   Local preparation tests do not establish native Molpro success or scientific
   sufficiency, and neither batch automatically approves a production model.
