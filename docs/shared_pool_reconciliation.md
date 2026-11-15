# Releasing finished QC allocation reservations

A shared pool counts both live Slurm allocations and submissions whose outcome
has not yet been reconciled. A short job can finish and leave Slurm's retained
job records before the next admission poll. Ordinary QC, including HIR and IRC
jobs, now has the same evidence-based reconciliation as dispatcher and orbital
reconnaissance jobs.

For an ordinary QC claim, KinBot requires all of the following:

- A successful current queue query, with the allocation absent.
- An explicit `Invalid job id specified` response from `scontrol`; an outage or
  an empty response is insufficient.
- A known submitted allocation and an exactly matching saved submission receipt.
- The generated QC wrapper's identity, guard token and output-marker arguments.
- A terminal wrapper execution record for that allocation's latest recorded
  requeue, or an allocation-tagged exit marker when no wrapper journal exists.

A failed calculation can release its allocation without accepting its result or
requesting a retry. A bare legacy exit code, a missing job ID, a mismatched
receipt, or a running/malformed requeue record keeps the reservation. Old
allocation journals remain usable even if a newer calculation has overwritten
the logical job's current exit marker. No recursive checkpoint search is used.

New controllers reconcile these records during normal admission. For an already
running controller, a separate process can reconcile the shared ledger under
the same admission lock:

```bash
python -B -m kinbot.run_policy reconcile /path/to/run
```

This command writes the reconciled ledger and evidence journal, so it is **not
read-only**. It does not submit, cancel, restart or modify any calculation or
controller. Existing controllers see the freed capacity on their next admission
poll. A fresh controller process is needed to load the new automatic behavior;
repeating the command safely handles further completed ordinary-QC claims from
controllers that were started with older code. Unconfirmed claims remain held.

The reconciliation journals are `.kinbot_policy_reconciled/<token>.json`; they
retain the original claim, scheduler diagnostic and hashed terminal evidence.

Tests exercise the actual shared admission and dispatcher with simulated Slurm:
eight live allocations plus twelve old QC claims admit exactly twelve further
tasks under a twenty-node cap. Repeated invocation preserves live allocations,
native files and previous results, with no duplicate submission. This is local
orchestration validation, not a verified remote scheduler handoff.
