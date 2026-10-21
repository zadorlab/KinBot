# Outer PES L2 selection before ANL

The previous audit patch exposed evidence but left a scheduling error: each
well could schedule ANL for every accepted search before the outer PES chose
which pathways to retain. This change moves that choice before ANL dispatch.

## Ordering and scientific scope

1. Every explored well publishes its accepted L2 network with the current
   controller token and exploration epoch, before acknowledging readiness.
2. After all discovered wells finish exploration, the outer PES applies its
   existing `select_summary_reaction` rule to the combined L2 evidence. It
   compares absolute L2 E+ZPE, recognizes reverse endpoints, preserves distinct
   stereochemical pathway classes, and applies the existing saddle/placeholder
   and exact-energy complex tie rules. This is the same kinetic selection
   policy used by final PES assembly, not a new geometry-equivalence test.
3. `pes_l2_selection.json` records selected and omitted observations. ANL uses
   only the selected saddles, stable products, retained complexes, and their
   recursively required CBH references. Exact-identity aliases can consume an
   existing canonical graph even when that owner's observation is selected out.
4. VRC receives the full accepted network and retains its independent execution
   and active-space review path. Reviewed VRC channels now start after their
   own well completes exploration, without waiting for the other wells. VRC
   review does not occupy an exploration slot or block ANL after L2 selection.
5. Handoff, worker MESS inputs and final summaries use the selected network.
   Omitted successful searches are labelled `SELECTED_OUT`, not `FAILED`.
   Final ANL auditing records both the selected plan and original accepted QC
   catalog. Raw reaction monitors continue to report the exploration outcome.

Missing or inconsistent L2 evidence blocks new ANL planning rather than falling
back to scheduling every raw search. The independent VRC branch can still run.
Ordinary non-composite KinBot and explicitly unmanaged/legacy execution retain
their existing behavior. No molecule, account or cluster is hardcoded.

## Preservation and repair

`network_request.json` remains the accepted scientific catalog. No native QC
graph is renamed or deleted merely because the kinetic selection omits it.
Selected work adopts existing dispatcher state and job IDs; it does not create
a new graph for a resource change. Previously running, now unselected QC may
finish in its existing directory, but does not spawn further ANL work through
this selected driver. Such a completion is not automatically proof of scientific
acceptance; normal validation applies if that graph is later consumed.

A new exploration epoch recomputes the global L2 selection. Derived profiles,
formations and handoff files are archived when the selected set changes, with
an interruption-safe selection journal. QC graphs and native files stay in
place and can be reused if a later selection needs them again. This projection
does not overwrite an old graph with changed geometry inputs.

The explicit stopped-controller `--deduplicate-anl` repair now also compares
previously reserved PES graph owners. Previously those owners were skipped.
It prefers retained completed work, verifies each alias directly against its
final owner (no transitive tolerance assumption), and uses the existing
scheduler-verification/cancellation/archive procedure. Normal live owner
registration remains immutable to avoid reciprocal ownership races.

A borrower no longer waits for the canonical owner's *whole* CBH network to
finish. After releasing its private graph locks, it acquires the canonical
graph's exclusive drive lock and advances/adopts that graph. An unrelated
reference failure or an owner selected out of kinetics cannot strand it.

## Updating the current run

Pulling alone updates monitoring; pinned running controllers require the
controlled refresh to use the new scheduling behavior. The dry run reads the
current manifests and queue before any cancellation decision. The apply path
stops/replaces controllers, retains unrelated allocations, and cancels/archives
only graph duplicates supported by its scientific and scheduler checks.

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
/home/hram/KinBot/.venv/bin/python -u -m kinbot.run_refresh \
  /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage \
  --deduplicate-anl
# After the dry-run report completes successfully:
/home/hram/KinBot/.venv/bin/python -u -m kinbot.run_refresh \
  /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage \
  --deduplicate-anl --apply
```

The current run is still exploring. Its replacement controllers reconstruct
accepted chemistry through the existing caches, run reviewed VRC for completed
wells, finish remaining L1/L2 work, and then select the full L2 PES before
starting new ANL work. Existing scratch,
exclusivity, resource pool, exclusion, walltime and recovery settings are not
changed by this update. New VRC channels still require reviewed models.

## Independent VRC gate follow-up

The first selection patch left the global exploration wait ahead of both
branches. This unnecessarily held VRC. The wait now applies to ANL planning
and the final join, while VRC keeps the main thread needed by its signal
handlers. Even a one-slot input can keep ANL waiting independently of a VRC
review hold; scheduler admission still enforces the existing QC resource cap.

Live reaction/network edits increment the exploration epoch and wait for all
active VRC/ANL branch threads to acknowledge it before changing reaction
objects or plans. No allocation is cancelled by this coordination. Accepted
L2 snapshots and epoch acknowledgements are published together under the
owner lock. The monitor reports exploration slots and unstarted wells, making
the separate `simultaneous_kinbot` limit visible.

For an already reconciled run on `ebb3c3ff`, enable this follow-up with
`git pull --ff-only origin composite` and the ordinary
`python -u -m kinbot.run_refresh RUN --apply`. This refresh archives controller
records and replaces the controllers while retaining existing QC allocations
and results. Do not repeat migration, rejected-reaction, deduplication or
pending-replacement flags merely to enable the scheduling change. It does not
increase `simultaneous_kinbot`: the next well remains gated by that explicit
input setting while the current well's last search is running.

## Validation and limits

Tests include actual local subprocesses crossing the exploration/selection
barrier, reverse discoveries, stereochemical alternatives, strict handoff
coverage, filtering of unused CBH work, unchanged saved job IDs/native files,
canonical graph adoption, repeated selection/reopening, and reconciliation of
already registered graph owners. These tests use fixture QC artifacts; they
do not constitute a new Gaussian/Molpro/MRCC/CFOUR run on Blodgett.

Validation on October 21:

- Independent VRC gate follow-up: 211 passed, 1 skipped across phase, branch,
  selection, input refresh, monitoring, controller refresh, live control and
  network revision tests. Real subprocess tests cover VRC before the next
  well finishes and ANL completing during a VRC hold, with one-slot and
  multi-slot inputs. These are orchestration tests with fixture QC callbacks.
- Focused selection, refresh, ownership, handoff and compatibility suite:
  253 passed, 1 skipped.
- Full repository suite: 2,502 passed, 7 skipped, 11 failed. The same 11 failures
  were reproduced separately from an unmodified archive of pre-change commit
  `8e57a0f7`: CFOUR standalone-result fixture, two sulfur-example expectations,
  three CASPT2 production expectations, three older MRCI probe tests, orbital
  campaign batch test, and species-routing HOMO test. This is not a clean full
  suite, and this patch does not claim to fix those failures.
- Final selection-policy/monitor checks, including distinguishing old pinned
  controllers and stale epochs: 72 passed.

The combined assembly test follows the real plan, profile, recursive CBH and
handoff code with fixture electronic results. It covers both selected and
legacy networks, repeat execution, and preservation of native graph records.

This fixes the scheduling gap. It does not establish the cause of the three
reported CH3SS stage-four rejections without their native acceptance diagnostics,
or guarantee that a scientifically difficult calculation will converge. It
also cannot raise an existing Slurm allocation's limit without scheduler
permission. The existing verified numerical/timeout recovery policies remain
in force when their graphs are required and advanced.
