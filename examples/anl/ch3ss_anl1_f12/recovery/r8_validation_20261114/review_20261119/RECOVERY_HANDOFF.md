# R8 necessary paper-reconciliation calculations

The [decision ledger](CALCULATION_DECISIONS.md) now selects **four B2PLYP
scans and one missing reverse IRC**. Existing forward paths, fragment energies,
CC results and orbital assets are reused. The user's empty queue/process check
allows these owners to be started; both retain their existing lock/reconciliation
protections against a concurrent or repeated invocation.

## Update, stage and launch on Blodgett

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
KB_VALIDATION="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114"
KB_CLOSURE="$KB_VALIDATION/paper_closure_20261119"
KB_TOOLS="$PWD/examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114"
.venv/bin/python -B "$KB_TOOLS/paper_closure.py" stage "$KB_CLOSURE"
nohup "$PWD/.venv/bin/python" -B -u "$KB_TOOLS/paper_closure.py" run "$KB_CLOSURE" \
  >> "$KB_CLOSURE/owner.log" 2>&1 < /dev/null &
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.saddle_recheck run \
  "$KB_VALIDATION/remaining_l2" --only ch3ss_rxn_1_1 --paths-only \
  >> "$KB_VALIDATION/remaining_l2/selected_owner.log" 2>&1 < /dev/null &
BASH
```

`stage` submits nothing. The two `run` commands request the five initial
allocations. Their shared policy preserves the exclusion list and capacity
limit. The scans use eight cores per exclusive node. An existing owner or
allocation is adopted, not replaced. Finished records are preserved on repeat.

The selected IRC run stops before endpoint-frequency submission so existing
spectra can be matched first. Its `controller.selected.json` receipt describes
only this selected run; it does not mark the 44-case campaign complete.

## Read-only progress

```bash
squeue --me -h -o '%i|%T|%j|%Z'
KB_VALIDATION="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114"
cat "$KB_VALIDATION/paper_closure_20261119/progress.json"
tail -n 25 "$KB_VALIDATION/paper_closure_20261119/owner.log"
tail -n 25 "$KB_VALIDATION/remaining_l2/selected_owner.log"
.venv/bin/python -B -m kinbot.saddle_recheck status "$KB_VALIDATION/remaining_l2"
```

Geometries and native output: `paper_closure_20261119/<scan>/tasks/scan/scan.log`.
The original selected IRC geometry/path files remain in
`remaining_l2/ch3ss_rxn_1_1/`; its endpoint structures are also in `kinbot.db`.

## Capture completed evidence

```bash
cd "$HOME/KinBot" || exit
.venv/bin/python -B \
  examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/collect.py collect \
  "$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114"
```

The collector now includes the new panel's plan, native inputs/outputs,
execution/progress receipts and derived scan observations. It snapshots the
remaining-L2 databases consistently and verifies archive hashes. It writes an
export under `~/kinbot-evidence` and does not alter or submit calculations.

## What completes the scientific decision

The scan owner extracts **correlated B2PLYP energies**, converged constrained
geometries, coordinate progress and candidate maxima. Native job completion is
separate from a complete constrained trace. The next review maps every point's
H attachments, spin/state continuity and inner/outer branch. A continuous
maximum requires an unconstrained saddle plus two-sided path verification;
a switch or failed trace is diagnosed before any further job is selected.

After that topology is established, freeze the energy convention, distinct
versus duplicate paths and population factors. Preserve literal SI and
energy-consistent Eckart variants: two SI W1-side depths differ by 3 kcal/mol.
Corrected kinetic/scope sensitivity and exact orbital/state asset binding
remain required before R9. MESS is still paused.

## Verification and scope

Local tests exercise the real dispatcher and selected owner with simulated
Slurm: staging, capacity deferral, active allocations with arriving outputs,
interrupted owners, restart, repeated invocation, locks, input mutation and
scheduler outages. A separate replay of the captured `ch3ss_rxn_1_1` files uses
the actual QC status reader and requests exactly the missing reverse path.
The test counts and artifact hashes are in [verification.json](verification.json).
Native completion of these new calculations is not yet verified.

No CC retry, whole-panel sampling restart or new orbital screen is part of this
handoff. The published MAXIT delimiter repair remains shared software
maintenance; it does not require rerunning an optional diagnostic.
