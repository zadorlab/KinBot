# R8 necessity-limited handoff

Use [CALCULATION_DECISIONS.md](CALCULATION_DECISIONS.md) for the current
scientific decisions and calculation admission. **No new QC task is selected.**
The earlier conditional TS1-F CC retry command has been withdrawn from this
handoff: repairing a software defect does not establish a scientific need to
rerun its affected calculation.

## Update and capture the remaining operational evidence

Run on Blodgett from the existing `composite` checkout:

```bash
cd "$HOME/KinBot" || exit
git pull --ff-only origin composite
.venv/bin/python -B \
  examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/finalization_snapshot.py \
  "$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114"
```

The snapshot reads the current owners, matching process wait states, queue,
node policy, L2 owner-log tails, and scheduler/native evidence for two truncated
failures. It submits, refreshes, cancels and edits no run. The update changes
repository files only; it launches no jobs. Missing executables or inaccessible
process details are reported in the capture.

The 13 completed path/endpoint pairs in the remaining-L2 panel must be reused.
Neither a stale owner record nor a missing endpoint-frequency record alone
justifies repeating their parent paths. The read-only capture is needed to
choose an owner action without duplicate submissions; it is not another
scientific calculation.

## Retry decisions

- No CC TS2/TS3 reoptimization, gradients, Hessians or paths are selected.
- No CC TS5 path or W1 refinement retry is selected.
- No compact `12_1` CC retry or `12_2` duplicate departure trial is selected.
- No whole 44-case restart, `resume_sampling.sh`, or panel-wide refresh is selected.
- Additional QC requires the exact model-relevant gap and reuse check in the
  decision ledger. R9 remains unapproved while its retained model gaps are open.

The shared MAXIT fix in `adf3dcba9` remains valid software maintenance: preserve
Molpro's comma delimiter and permit a one-time, same-limit repair only for
matching native input-error/retry provenance. Its 266-test local regression
covers actual dispatcher orchestration with a simulated scheduler, archival,
repeated invocation, live-work preservation, probe decisions and sampling.
No repaired native CC execution has been claimed or requested here.
