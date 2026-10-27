# Outage preservation and recovery boundaries

`run_outage_audit` captures evidence without starting an owner, dispatcher or
recovery routine. `checkpoint_preserve` copies named native checkpoint files
without accepting them as results or altering their source. Neither command
restarts chemistry, cancels allocations, changes scientific settings, or prunes
the PES/CBH network.

## Capture before deciding on recovery

```bash
python -B -m kinbot.run_outage_audit RUN \
  --output-parent /persistent/path/outage-evidence --job OLD_JOB_ID
```

Repeat `--job` for up to 20 important old allocations. The capture reads the
current queue, duplicate accounting records for named/current scoped jobs, and
`scontrol` for named jobs. An unavailable scheduler remains **unknown**, never
an empty queue. Inputs, task ledgers, accepted networks, review requests,
runtime pointers, small trajectories and output tails are archived with hashes.
Large checkpoints, databases and sampling pickles are inventoried only. Pickles
are never loaded. Budget omissions, read failures and changes during capture
are explicit. This is not an atomic filesystem snapshot or a completed-result
validation. The source tree is not written; use an output directory outside it.

## Preserve node-local checkpoints before releasing allocations

Use this only after stopping submission owners and sending STOP to the run's
running allocations, with a receipt recording their IDs. Slurm STOP signalling
does not release the allocations and may leave their queue state RUNNING.
Pending jobs have no running native process to freeze. Ordinary cancellation
may trigger site cleanup; it does not guarantee Python `finally` blocks run.

```bash
python -B -m kinbot.checkpoint_preserve /path/to/audit/report.json \
  --output-parent /persistent/path/checkpoints \
  --allocations /path/to/outage_freeze.RECEIPT.json
```

The controller-host command requires locally discovered run processes to be
stopped and jobs still to be RUNNING under the recorded root. It starts at most
two concurrent `srun --jobid=... --overlap` copy steps in **existing** allocations,
one per node. It does not request new allocations or resume stopped chemistry.
Site restrictions can prevent those steps; failed steps and missing nodes remain
unresolved and need site-assisted copying. Never infer success just because
the launcher exits: review `steps.json` and every `preservation.json`.

The node command checks host, ownership, directory scope, symlinks, disk space
and copy size. It saves `.wfu`, `.chk`, `.rwf`, MRCC `fort.16`/`fort.17`, selected
CFOUR archives, inputs, geometry/XML files, and bounded diagnostic tails. Large
integral files are deliberately outside this checkpoint collection. Each copied
file has a SHA-256 verified against its destination and before/after source
metadata. Stable bytes do not establish an internally complete checkpoint.
The default copy ceiling is 128 GiB per node; failures/omissions are reported.

The allocation orchestrator also preserves named checkpoint artifacts from
owned sibling `kinbot-*` scratch directories. These may include older attempts
whose runtime pointer was overwritten, or unrelated KinBot work. They are
quarantined with **unknown attribution** and cannot be reused until matched to
the correct run, geometry, basis, state and native attempt. Source files are
never changed or removed. Direct node invocation only includes these siblings
when `--include-unattributed` is requested.

Historical workspaces on nodes with no surviving allocation are not reachable
by this mechanism. Give administrators the scratch host/path inventory and
request retention of those workspaces and any pre-outage filesystem snapshots.

## What is and is not currently restartable

| Calculation | Retained progress and required checks |
| --- | --- |
| Completed L1/L2/PES work | Preserve databases, accepted networks, inputs and owner generation records. Reconcile with task/native provenance; do not repeat accepted exploration wholesale. |
| Completed ANL components | Existing dispatch validation checks staged inputs, geometry and execution artifacts. Native result parsers remain required. A successful native output without a wrapper receipt needs explicit reconciliation. |
| ASE/Sella L3 geometry | A valid evaluated trajectory frame or a provenance-matched native force input/output can supply a continuation geometry. The existing automatic continuation handles **step exhaustion**, not arbitrary mid-step outages. It starts a new optimizer; it does not restore its internal Hessian/history. |
| Automatic Slurm requeue of ASE geometry | The existing task runner can restart from `geometry.xyz`/an explicitly bound continuation, reuse numbered calculator filenames, and reopen `optimization.traj`. Requeue is **not** evidence of native continuation. Preserve remaining old files before any new task attempt. |
| Molpro numerical Hessian | Molpro supports restart using saved numerical-frequency task records (default `5500.2`, with Hessian/frequency records). An arbitrary orbital-only WFU is insufficient. KinBot does not yet implement a validated general outage restart for this task. |
| MRCC high-order single point | Existing `anl.dispatch resume-mrcc` can stage `rest=1` from retained amplitudes after inactive-job, input and dependency checks. Interrupted node-local amplitudes must first be recovered and validated. It is not safe to apply to an active/requeued allocation. |
| Gaussian VPT2 / other native jobs | Preserve checkpoint plus relevant native workspace artifacts. A `.chk` can be a seed without containing resumable progress for the interrupted phase. No universal VPT2 continuation is established by this audit. |
| VRC scans and orbital reviews | Retain per-point native results, basis-specific orbital records, hashes and manifests. Controller `VALIDATED` means execution checks, not scientific approval of active character or roots. |
| RotdPy | Preserve sampling payloads, surface inputs and attempt journals. Validate completed samples and reconcile live allocations through the owning campaign before resuming; do not delete/reinitialize samples to recover the owner. |

Sources: [Molpro frequencies](https://www.molpro.net/manual/doku.php?id=harmonic_vibrational_frequencies_frequencies),
[MRCC manual](https://www.mrcc.hu/MRCC/manual/pdf/manual.pdf),
[Slurm scancel](https://slurm.schedmd.com/scancel.html),
[Slurm srun](https://slurm.schedmd.com/srun.html). Native support also depends on
the installed program version and which records survived.

## Migration to a new run

Keep the old run at its original path, with accepted results and failed attempts
intact. A new run needs an explicit import ledger for reused results, checkpoint
continuations and unresolved tasks. Copying a directory or rewriting a status to
`complete` is not a validated migration. Existing inter-PES graph sharing checks
identity and has lifecycle constraints; it is not a general partial-graph import.

Before starting the new run, prevent unreviewed Slurm requeues from overwriting
prior attempts, finish the owner/restart reconciliation, and test repeated
orchestration with live-job preservation and sampling adoption. Restarted native
calculations must then demonstrate that the intended checkpoint was read.

A proposed 150 kcal/mol cutoff should use each PES's own initial reference and
consistent L2 electronic energies plus ZPVE. Distinguish TS energies, well
energies and sums of product energies. Keep missing/uncertain cases for review;
the cutoff is a pruning proposal, not proof of negligible kinetics at every
temperature/energy. CBH references have different stoichiometries: remove them
only when no retained thermochemistry target depends on them. Final pruning
requires exploration closure; incomplete exploration permits provisional lists.
