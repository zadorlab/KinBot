# October 22 geometry/VPT2 repair handoff

Follow-up: the native C2V,Loose attempt also failed. The extended
[general VPT2 procedure](VPT2_RECOVERY.md) now shares a complete strategy history,
including archived operator attempts, and considers an untried PG=C1 cap after
the native-group alternatives. Load it with one code refresh using the command
below; subsequent eligible retries are automatic. This is a recovery-policy
update, not a claim that the unresolved native calculation has succeeded.

Latest correction: the subsequent refinement refresh failed before replacing
controllers because an independently running ROTD campaign was classified as
an unknown child. The old `controller_refresh.json` still described the earlier
successful refresh. The queued requests were not consumed or lost.

Refresh now retains verified campaign owners, their descendants and their
pool-owned allocations. It checks input hashes and process start identities,
and coordinates with the pool submission lock before stopping PES/KB owners.
Unregistered children and allocations still prevent a handoff.

Native `hsdel` / singular-expansion failures now have at most two automatic
retries in generated Molpro geometry tasks: a changed primary P-space, then a
smaller internal expansion workspace and tighter energy/coefficient convergence.
The active space, requested roots, weights, NSTATI, CASPT2 shift and IPEA remain
unchanged. The failed attempt is archived; successful energies and gradients
must still pass the original validators. Unknown errors and a user-specified
P-space do not trigger this automatic ladder. Deliberately bounded fragment
probes retain their explicit retry budgets.

The affected production scan is CH3S -> CH3 + S, doublet root 2, point 16
(5.2 A), not C + SH3. Root 1 already has a resolved 48-point scan. Cached L2
fragment energies in the supplied audit give a provisional harmonic D0 of
69.00945 kcal/mol, not a completed root-2 variational barrier.

The previous 16:08 background controller refresh did launch the new PES/KB
owners in the supplied 16:29 monitor. Existing QC was retained and two independent
ROTD campaigns launched. Their then-current 0/116 counts did not yet demonstrate
sampling progress or completion.

The following changes are ready for one code refresh. Future requests through
`run_control refine` and the native VPT2 recovery need no controller restart.

- Geometry/frequency refinement follows the dependency graph at TZ or QZ, for
  minima or saddles. Completed dependent single points are archived and
  recomputed at the replacement geometry; unrelated calculations are retained.
- The selected geometry uses fmax=0.0001 eV/angstrom. Its numerical derivatives
  and dependent harmonic Hessians use SCF/CC convergence of 1e-11, preserving
  any stricter setting. The original stationary-point acceptance remains.
- VPT2 recovery reads Gaussian's native framework assignments and automatically
  derives a bounded retry, including a Loose detection retry after corroborated
  initial/mass-framework disagreement. It has no molecule/group lookup table.
- An explicit timeout retry can use saved allocation ownership plus its matching
  TIME LIMIT diagnostic after Slurm purges the allocation and no worker receipt
  was written. Deliberately bounded VRC probes still do not retry endlessly.
- Reference examples have moved outside the installed package into `validation/`.

## Scientific interpretation of the current TS

The retained reaction is CH2SH -> H2 + HSC, discovered under the
`r12_insertion_R_5_1_4` search name. The latest logs report successful IRC endpoints,
L2 identity RMSD 0.114739 A and mode correlation 1.00. The optional final product
complex was omitted because its depth was -0.365586 kcal/mol. That does not reject
the retained TS/product channel. Within the configured L1/L2 acceptance protocol,
this is an accepted discovery; the search-family name is not its reaction identity.

L2 row 376 has one imaginary mode (-57.8313 cm^-1), with two other soft modes
at 50.5039 and 109.2887 cm^-1. The TZ geometry reached fmax=0.000457302 eV/A,
below its original 0.0005 target, and its energy matches the harmonic calculation
to about 9e-11 hartree. The geometry passed the force criterion. The TZ Hessian
has three imaginary modes (237.78i, 167.53i, 77.93i), so it does not certify a
first-order saddle. The remote H atom was already present in the accepted L2
geometry. These data do not distinguish numerical curvature error from a real
change of stationary-point order at L3. The paired retry tests that question;
it does not presume the answer or label the channel a fake reaction.

## Refresh with the already queued retries

Run on Blodgett. The TS refinement and 12-hour fragment retry are already
queued; do not publish either request again. Requests are
consumed by the owning graph when admitted; they do not bypass PES admission
or shared capacity. Only affected allocations may be cancelled by this explicit
request; the refresh itself retains native QC allocations.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_PY="$HOME/KinBot/.venv/bin/python"
KB_LOG="$KB_RUN/controller_refresh_$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$KB_PY" -u -m kinbot.run_refresh "$KB_RUN" --apply \
  > "$KB_LOG" 2>&1 < /dev/null &
printf 'Refresh PID: %s\nLog: %s\n' "$!" "$KB_LOG"
BASH
```

Once the new refresh reports `phase: launched`, VPT2 and the native singularity
recovery are automatic when their owning graphs are polled.
It need not receive a hardcoded C2V retry command. The current native evidence
is exactly the diagnostic accepted by the new general rule. Neither retry has
been executed remotely from the development checkout; native convergence is
still to be verified.

The read-only status commands are:

```bash
"$HOME/KinBot/.venv/bin/python" -m kinbot.run_control show \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
"$HOME/KinBot/.venv/bin/python" -m kinbot.run_monitor \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
```

Do not republish an already queued request. If its receipt says `blocked`
because a task generation changed, inspect the new state and preview a new
request. A queued receipt is not a submitted calculation; an applied receipt
records archival/reset, and the subsequent task state supplies the new job ID.

## Explicit budget for the timed-out fragment attempt (already queued)

Request 897941f10f974b7a94e2fbbef40b49dd already queues one retry of the
recorded failed point with 12 hours and at most 300 optimizer steps, retaining
its original force tolerance and scientific model. The old attempt had only
2 hours and stopped at step 98, before convergence. This starts from its
original declared seed; it does not claim a valid converged branch. A still
active or already successful attempt is retained. The owning scan may finish
its other points before consuming the request.

```bash
"$HOME/KinBot/.venv/bin/python" - <<'PY'
from pathlib import Path
import json
from kinbot import live_control
r = Path.home()/'kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7'
paths = list(r.glob('ch3s/*/vrctst/molpro_scan/*vts_pt07_fragment_retry*/tasks/relaxed_scan'))
if len(paths) != 1:
    raise SystemExit(f'Expected one saved attempt, found {len(paths)}; no request created.')
d = paths[0]
s = json.loads((d.parent.parent/'state.json').read_text())['tasks']['relaxed_scan']
if s.get('status') not in ('failed', 'unavailable'):
    raise SystemExit(f'Retaining current attempt: {s.get("status")}; no request created.')
print(json.dumps(live_control.request(r, d,
    reason='One explicit longer retry after confirmed two-hour timeout before convergence',
    fallback={'schema': 1, 'resources': {'walltime': '0-12:00:00'},
              'optimizer': {'steps': 300}}), indent=2))
PY
```

Connectivity gaps at other scan points remain scientific review items. No
successful outcome, stationary-point order, or complete VRC scan is manufactured
by these changes. Historical L1/L2 search rejections are separate from required
production-component failures and are not all resubmitted by this handoff.

## Read-only paper comparison

Run after pulling the checkout. This command changes no run file, submits
nothing and prints the comparison plus coordinates, frequencies and energies
needed for the remaining paper validation. Paste the output for the comparison.
Paper labels/energies are explicit reference data outside `kinbot/`.

```bash
cd "$HOME/KinBot"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -B -m validation.pes_compare \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/ch3ss" \
  --reference validation/data/d4cp04873c_fig6.json --evidence
```

The supplied 17:18 read-only snapshot has accepted W2 and W4 connectivity,
the W4-W2/TS4 endpoint and path-type match, and W4-P4 scission. W3 has been
discovered, with current acceptance unavailable. W2-P1 is a candidate in
progress; this does not establish the paper's W1 complex or TS2. W1, P2, P3,
and TS1/2/3/5/6 are not established in the current accepted records. The two
accepted TS4 search names have essentially identical energies/frequencies and
must not be counted as two distinct physical mechanisms from their names.

Additional accepted connections outside Figure 6 are W4 -> CH3S + S and
W4 -> [H][C]([H])=[S]([H])[S]. The latter is an additional H-shift isomer.
Further candidates are still in progress. There were 5 accepted searches,
33 rejections, 141 searches in progress and one unavailable worker monitor.

Provisional harmonic L2 E0 values relative to cached CH2S + SH (kcal/mol):

| Point | Current L2 | Paper Figure 6 |
|---|---:|---:|
| W2 | -18.03 | -18.46 |
| W3 (acceptance pending) | -22.85 | -22.82 |
| W4 | -43.88 | -43.02 |
| TS4 | 7.63 | 7.38 |
| P4 | -0.28 | 3.57 |

These values use the supplied native electronic energies plus harmonic ZPE,
with P1 = -835.934570239890 hartree. They are not final ANL/L3 values; the
paper uses a different composite protocol. P4 differs by 3.85 kcal/mol and
requires checking at the final tier. The additional isomer is at +11.01,
its connecting TS at +26.85, and CH3S + S at +30.52 kcal/mol on this same
provisional zero. Product cache rows without saved calculator parameters do
not establish spin/theory metadata; the audit now prints unknown rather than
inventing a default singlet. No result is fitted to the paper.

An unfinished/unmapped snapshot cannot prove a channel absent. Endpoint
matches do not establish exact saddle/conformer identity, particularly W1.
