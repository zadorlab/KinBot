# Shared Slurm resource pools and live controls

KinBot can share one allocation budget across multiple PES controllers and
their ordinary QC, CBH/ANL, VRC, ROTD_py sampling, and MESS/MESMER jobs.
An opt-in `.kinbot_resource_policy.json` in a common ancestor directory
defines the pool. No account, molecule, partition, node name, or site capacity
is hardcoded into the implementation. Standalone orbital-reconnaissance
campaign tools retain their separate campaign controls; do not run those
independent submitters inside an active shared-pool tree.

Without a policy, existing input defaults and per-branch limits remain in
effect. With a policy, its live limit replaces the ordinary QC, dispatcher,
ANL branch, and ROTD job-count limits. The frozen scientific plans retain
their original fields and hashes, so adding a pool does not invalidate
accepted geometries, single points, or samples. These older job-count fields
are no longer reservations in shared mode. Resource controls do not change
cores, memory, scratch selection, convergence thresholds, or chemistry.

## Enable a pool

For a single new run, the existing JSON input can contain:

```json
"queue_live_policy": {
  "max_jobs": 12,
  "excluded_nodes": ["site_bad_node"]
}
```

The default is `{}` (disabled). This initializes a policy alongside the input,
unless an ancestor policy already governs the run. Input defaults never
overwrite subsequent live edits. This option requires `queuing: "slurm"`.

For multiple runs, initialize their **common ancestor** before launching or
refreshing them:

```bash
python -m kinbot.run_policy set /path/to/campaign --max-jobs 23 \
  --exclude site_bad_node,reserved_node
```

Prepared campaign inputs and known worker/dispatcher node exclusions are
preserved when enabling the pool. Explicit node names are required; bracketed
Slurm hostlist syntax is not accepted by this CLI. Each new allocation uses
one exclusive node. Retained allocations consume their reported node count.
Nested pools are rejected because they could split a supposedly shared budget.

All cooperating controllers must run code with shared admission support.
Creating the policy cannot hot-patch an older running Python process. For an
existing salvage run, install the code and perform the existing controller
refresh once. Future policy edits do not need a refresh.

## Change the running pool

```bash
python -m kinbot.run_policy set /path/to/campaign --max-jobs 18
python -m kinbot.run_policy set /path/to/campaign --exclude nodeA,nodeB
python -m kinbot.run_policy set /path/to/campaign --allow nodeA
python -m kinbot.run_policy set /path/to/campaign --max-jobs 0
python -m kinbot.run_policy show /path/to/campaign
```

`--exclude` adds exclusions; `--allow` removes them, including exclusions
learned before the edit. A newly diagnosed node failure can exclude that node
again. These flags change the pool's eligibility, not the administrator's
Slurm node state or partition access. `max_jobs=0` pauses new admissions;
lowering a limit lets existing allocations finish. Changes are read on the
next admission check (normally within the controller polling interval).

Exclusions apply to **new submissions**. Running and already pending jobs keep
their submitted settings. Use the existing `run_refresh --replace-pending`
procedure only when pending jobs must be replaced: it holds and verifies a
pending job before cancelling it and retains a job that started in the
meantime. Running calculations are never moved by a policy edit.

## Conservative admission and dependencies

* Each submitted job reserves capacity, including PENDING, RUNNING, and
  COMPLETING allocations. This prevents simultaneously starting pending jobs
  from exceeding the budget.
* KinBot tasks waiting on geometry or energy dependencies have **not** been
  submitted and reserve no capacity. Their accepted parent results are needed
  to generate and validate downstream QC inputs. This feature does not
  pre-submit a speculative dependency queue to Slurm.
* Other branches and PES runs can use every unreserved slot. There are no
  fixed half-and-half reservations. Submission is opportunistic; there is no
  guarantee of equal shares or an exact global start order.
* A process-shared file lock covers occupancy lookup, the durable admission
  claim, and the `sbatch` call. A thread lock also serializes threads within a
  controller. The filesystem hosting the pool must support coherent POSIX
  advisory locks and atomic renames, as required by the existing graph locks.
* Live retained jobs under the pool directory count even if submitted before
  the policy existed. Symlink aliases are resolved when comparing workdirs.
* Scheduler query failure pauses new submissions. A timed-out or unrecognized
  submission response retains an uncertain claim; it is not blindly retried.
  An acknowledged job absent from the first queue poll remains reserved until
  queue visibility or terminal scheduler evidence reconciles it. Unresolved
  claims are visible with `show` and require investigation, not ledger deletion.

Policy revisions are saved in `.kinbot_policy_history/`, and accepted
submissions record the policy revision and exclusions in
`.kinbot_policy_submissions/`. The pending-admission journal is
`.kinbot_policy_admissions.json`. The monitor shows the shared limit instead
of stale per-case caps. The read-only audit checks new jobs against their
submission revision and distinguishes retained pre-policy allocations.

Admission supplies `--nodes=1`, `--exclusive`, and `--exclude` on the command
line while retaining the scientific scripts. Slurm command-line options take
precedence over environment and script options. See the official
[sbatch documentation](https://slurm.schedmd.com/sbatch.html).

## Blodgett handoff, 2026-10-19

The supplied snapshot contains 29 usable nodes. Excluding the seven
drained/down/unknown nodes plus reserving `node9,node16,node21,node22,node24,node34`
leaves at most **23** for the campaign. Five reserve nodes have existing
exclusions or recent scratch-failure evidence; `node22` is the extra reserve.
This leaves those six outside this campaign; it cannot reserve them against
other users. Slurm continues to enforce node state, partition, and site policy.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
git pull --ff-only origin composite
KB_PY=/home/hram/KinBot/.venv/bin/python
KB_RUN=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage
"$KB_PY" -m kinbot.run_policy set "$KB_RUN" --max-jobs 23 \
  --exclude node1,node9,node10,node13,node14,node16,node21,node22,node23,node24,node30,node33,node34
"$KB_PY" -m kinbot.run_policy show "$KB_RUN"
"$KB_PY" -u -m kinbot.run_refresh "$KB_RUN" --apply --replace-pending
"$KB_PY" -m kinbot.run_policy show "$KB_RUN"
BASH
```

This archives/replaces only the controllers and verified pending attempts.
It preserves the run directory, accepted results, and running allocations.
If a local helper is active, refresh deliberately refuses to stop it; repeat
the same command once that short helper exits. If scheduler identity or hold
status cannot be verified, it stops rather than risk duplicate submissions.

Keep the existing monitor command pointed at the same run root. Scratch
remains the verified node-local policy, preferring `/scratch/<actual account>`;
the shared-pool feature does not relax it. Existing bounded scratch/I/O retry
rules remain in effect. An arbitrary solver failure is not reclassified as
recoverable merely to use more nodes.
