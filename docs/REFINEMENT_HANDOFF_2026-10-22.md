# October 22 geometry/VPT2 repair handoff

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

## Refresh and request the paired retry

Run on Blodgett. The distinct refinement request kind is ignored by older
controllers, so it can be queued before the single code refresh. Requests are
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
"$KB_PY" -m kinbot.run_control refine "$KB_RUN" \
  --target ch3s/471030810040000000002/anl_automatic/states/ts_0450141d9a891e9c/composite/tasks/harmonic_tz \
  --fmax 0.0001 --scf 1e-11 --cc 1e-11 --cancel-running \
  --reason 'Paired geometry/Hessian refinement after TZ stationary-point validation failure' --apply
nohup "$KB_PY" -u -m kinbot.run_refresh "$KB_RUN" --apply \
  > "$KB_RUN/refinement_controller_refresh.log" 2>&1 < /dev/null &
printf 'Refresh PID: %s\nLog: %s\n' "$!" "$KB_RUN/refinement_controller_refresh.log"
BASH
```

Once the refresh reports `phase: launched`, the new VPT2 recovery is automatic.
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

## Explicit budget for the timed-out fragment attempt

After the refreshed controllers are loaded, this queues one retry of the
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

An unfinished/unmapped snapshot cannot prove a channel absent. Endpoint matches
are reported separately from paper saddle/conformer identity, especially the
W1 van der Waals complex. A full numerical comparison still awaits this output.
