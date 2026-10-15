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

## Resume Batch A after the node21 startup failures

The retained records for jobs 52907719 and 52907720 both identify node21,
`FAILED`, `ExitCode=0:53`, and zero/one second of runtime. Their directories
contain only the prepared input and batch script. This is consistent with
Slurm's batch I/O initialization failure before the shell/guard starts, not
Molpro nonconvergence. Slurm uses `ESLURMD_IO_ERROR=4021`; its reported signal
portion can be 53. See the upstream
[error definitions](https://github.com/SchedMD/slurm/blob/master/slurm/slurm_errno.h)
and [step manager](https://github.com/SchedMD/slurm/blob/master/src/slurmd/slurmstepd/mgr.c).
The exact filesystem/stdio problem requires node-side diagnostics; the
provided records do not establish a quota or permission cause. Job 52907713
was separately requeued and was running on node16; it should be preserved.

The controller now retains `scontrol` evidence and automatically handles this
narrow startup failure with bounded retries and campaign-scoped host
exclusions. This one-time block additionally authorizes retry of the **11
reported empty attempts**, if their scheduler records have already expired.
It skips live jobs, changed IDs, runtime evidence, and contradictory scheduler
records. This explicit incident repair is distinct from the conservative
automatic policy. Successful results, active workers, input chemistry and
checkpoint sources are preserved. Do not recreate or clean the campaign.

The block stops only the controller, updates code, then resumes the same
campaign with a cap of 20. `set -e` applies to the child shell, not the SSH
login shell. Repeating it does not retry replacement job IDs.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
RUN=/home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015
"$KB_PY" -m kinbot.orbital_campaign stop-controller "$RUN"
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite
"$KB_PY" -m pip install -e . --no-build-isolation --no-deps
"$KB_PY" - "$RUN" <<'PY'
import json, subprocess, sys
from pathlib import Path
from kinbot import orbital_campaign as campaign, molpro_reconnaissance as recon
from kinbot.execution import atomic_json, slurm_job_record
root = Path(sys.argv[1])
reported = {str(i) for i in (52907719,52907720,52907721,52907722,52907723,
                           52907724,52907725,52907727,52907728,52907729,52907730)}
with campaign._lock(root):
    data = json.loads((root / 'campaign.json').read_text())
    live = set(subprocess.run(['squeue','-h','-o','%A'], check=True,
               capture_output=True, text=True, timeout=30).stdout.split())
    paths = {Path(s[k]) for s in data['systems'] for k in ('pilot_manifest','manifest')}
    for path in sorted(paths):
        m = json.loads(path.read_text())
        # This host is from the supplied incident evidence, not a code default.
        m['slurm_excluded_nodes'] = sorted(set(m.get('slurm_excluded_nodes', [])) | {'node21'})
        atomic_json(path, m)
        for j in m['jobs']:
            jid = str(j.get('job_id', ''))
            if jid not in reported or jid in live or j.get('status') == 'complete':
                continue
            d, name = Path(j['directory']), j['name']
            artifacts = [d/'execution.json', d/(name+'.native.inp'), d/(name+'.out')]
            streams = [d/(name+s) for s in ('.stdout','.stderr','.launcher.stdout','.launcher.stderr')]
            if (any(p.exists() for p in artifacts) or
                any(p.exists() and p.stat().st_size for p in streams) or
                any((d/'.kinbot_slurm_io_retries').glob('*/'+jid))):
                print('Preserving execution evidence:', jid, d)
                continue
            record = slurm_job_record(jid)
            if record and (record.get('JobId') != jid or record.get('WorkDir') != str(d)
                           or record.get('JobState') != 'FAILED' or record.get('ExitCode') != '0:53'):
                print('Preserving allocation with different scheduler evidence:', jid, record)
                continue
            if len(m['jobs']) != 1:
                raise SystemExit('Expected a single-probe incident manifest: '+str(path))
            j['manual_startup_recovery'] = {'job_id': jid, 'scheduler_record': record,
                'reason': 'Explicit retry of reported empty pre-start attempt; native files absent.'}
            atomic_json(path, m)
            print(recon.retry_failed(path, include_interrupted=True))
PY
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
export KINBOT_MOLPRO_SCRATCH_MIN_MB=65536
mkdir -p "$KINBOT_MOLPRO_SCRATCH"
"$KB_PY" -m kinbot.orbital_campaign preflight "$RUN"
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$RUN" \
    --pilot-only --probe-method casscf --max-jobs 20 --max-refinements 0 \
    >> "$RUN/controller.log" 2>&1 < /dev/null &
printf '%s\n' "$!" > "$RUN/controller.pid"
"$KB_PY" -m kinbot.orbital_campaign progress "$RUN" --failures-only
BASH
```

Compact monitoring (Ctrl-C stops only `watch`):

```bash
cd /home/hram/KinBot
watch -n 10 '/home/hram/KinBot/.venv/bin/python -m kinbot.orbital_campaign progress /home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015 --failures-only; squeue -u hram -o "%.18i %.10T %.48j %.10M %R"; tail -n 8 /home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015/controller.log'
```

## Resume Batch A and add the CH3 + S2 CAS(7,5) comparison

Use this block for the latest interrupted-worker failures. It stops only the
controller, keeps healthy Slurm allocations running, reuses validated completed
probes, and adopts five new CAS(7,5) probes under the same 20-job cap. Do not delete
the earlier review campaign: its reviewed WFU files are required seeds.

The added points are the matched center/left/right/tilted stencil and the 30 Å
reference. They compare CASPT2-F12 IPEA/reference treatments and MRCI-F12+Q at TZ.
The example recipe retains the methyl radical orbital and both S2 pi/pi* pairs;
it freezes the sigma bond and excludes sigma*. The far checkpoint requires an
explicit orbital exchange before partitioning. Checkpoint hashes bind that choice
to the reviewed source. This is a candidate, not a production CAS approval.

The request is **104000M per allocation**, below node2's reported **126000M**,
leaving 22000M outside the allocation. Eight Molpro ranks retain 1000 MW each;
the request includes program overhead and headroom. This is an allocation limit,
not proof of the program's peak memory use. Slurm will only schedule on nodes
meeting the request. The environment override also applies to older prepared
80000M scripts; each new submission records it in the manifest. Existing live
allocations keep their original resource limits.

Confirmed native I/O failures are retried automatically. An interrupted worker
with no final record is retried after the publication grace period. There are at
most two controller retries, preserving every attempt; the existing batch I/O
guard also has a bounded retry allowance. Explicit convergence failures and known
OOM, timeout or cancellation remain failures. Independent probes continue, and
the campaign reports failure at the end if any required probe remains invalid.
Retained `scontrol` records are captured without requiring `sacct`.

All shell settings below are confined to the child `bash`; an error leaves the
SSH login shell connected. The block is repeatable. No healthy job is cancelled.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
RUN=/home/hram/kinbot-runs/ch3s_ch3ss_reference_checks_20261015
"$KB_PY" -m kinbot.orbital_campaign stop-controller "$RUN"
env -u LD_LIBRARY_PATH -u LD_PRELOAD git -c submodule.recurse=false pull --ff-only origin composite
"$KB_PY" -m pip install -e . --no-build-isolation --no-deps
"$KB_PY" examples/anl/orbital_reviews/prepare_cas75.py "$RUN/cas7_5_addition"
"$KB_PY" -m kinbot.orbital_campaign append "$RUN" "$RUN/cas7_5_addition"
unset SBATCH_MEM_PER_CPU SBATCH_MEM_PER_GPU
export SBATCH_MEM_PER_NODE=104000
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
export KINBOT_MOLPRO_SCRATCH_MIN_MB=65536
mkdir -p "$KINBOT_MOLPRO_SCRATCH"
"$KB_PY" -m kinbot.orbital_campaign preflight "$RUN"
"$KB_PY" - "$RUN" <<'PY'
import json, subprocess, sys
from pathlib import Path
c = json.loads((Path(sys.argv[1]) / 'campaign.json').read_text())
m = json.loads(Path(c['systems'][0]['pilot_manifest']).read_text())
j = m['jobs'][0]
subprocess.run(['sbatch', '--test-only', '--mem=104000M', j['name'] + '.sbatch'],
               cwd=j['directory'], check=True)
PY
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$RUN" \
    --pilot-only --probe-method casscf --max-jobs 20 --max-refinements 0 \
    --keep-going --allocation-retries 2 \
    >> "$RUN/controller.log" 2>&1 < /dev/null &
printf '%s\n' "$!" > "$RUN/controller.pid"
"$KB_PY" -m kinbot.orbital_campaign progress "$RUN" --failures-only
BASH
```

Monitor with the compact `watch` command above. Export with the existing block
below; the parent archive includes the nested addition and its review provenance.
Do not start a separate controller inside `cas7_5_addition`.

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
unset SBATCH_MEM_PER_CPU SBATCH_MEM_PER_GPU
export SBATCH_MEM_PER_NODE=104000
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
export KINBOT_MOLPRO_SCRATCH_MIN_MB=65536
mkdir -p "$KINBOT_MOLPRO_SCRATCH"
"$KB_PY" -m kinbot.orbital_campaign preflight "$RUN"
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$RUN" \
    --pilot-only --probe-method casscf --max-jobs 20 --max-refinements 0 \
    --keep-going --allocation-retries 2 \
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
"$KB_PY" -m kinbot.orbital_campaign collect "$OLD" \
    /home/hram/kinbot-runs/ch3s_ch3ss_potential_validation_latest_20261015.tar.gz
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
   The [PCCP 2007 guidance](ORBITAL_REVIEW_WORKFLOW.md#what-accuracy-means-for-this-workflow)
   also motivates the CH3 + S2 CAS(7,5) pi-manifold addition below. Its native
   results remain unvalidated; the base A/B recipes do not automatically
   install it or supply its DZ basis selection. If it survives the TZ checks,
   verify the corresponding DZ orbital characters before testing DZ sampling.
   CAS(9,7) must not be presented as a demonstrated minimum-cost choice merely
   because it was the smallest candidate in the original recipe.
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
