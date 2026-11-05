# Restart handoff corrections, 4 November 2026

These changes address coordination and result handoff defects. They do not
change electronic structure methods, stereochemical acceptance, frequency
validation, rotor eligibility, or the homolysis evidence required for PST export.

- Conformer selection publishes the complete accepted calculation to its `_low`
  alias, including Gaussian/Sella completion and native Hessian files and source
  row provenance. Historical aliases lacking those files are rebuilt from the
  original conformer calculations. Failed reads never supply zero placeholder
  coordinates to the stereochemistry classifier.
- HIR-requested refinement retires the old completion receipt, pickle handoff,
  result record and optimization budget along with native outputs. Historical
  or interrupted invalidations can finish retiring pre-invalidation artifacts;
  active jobs and artifacts from a newer generation remain untouched. The
  invalidation must be explicit: this does not retry arbitrary failed saddle
  searches. Accepted-search replay may execute that previously authorized HIR
  refinement, but still refuses unrequested replacement QC.
- JSON publishers use unique same-directory temporary files, flush and fsync
  before rename, verify written bytes, and sync the directory. PES coordination
  retains a second copy of the same epoch under its existing lock. Corrupt
  primary bytes are retained. If both copies are unreadable, workers do not
  invent an empty exploration; native PES startup can reconstruct coordination
  under its exclusive driver lock. These measures cannot guarantee writes on
  a faulty filesystem or prove that NFS caused a historical corruption.
- Restart acknowledgment validates installed reviewed models only. It never
  invokes automatic radical-pair verification on identity-only fragment stubs.
  Verification runs later under the VRC owner with actual accepted products.
  Bookkeeping/storage/parse errors record `execution_error`, not orbital rejection.
  An old unassessed record containing only an identifiable JSON-parser failure
  is archived before resuming. Actual failed scientific checks remain held.
- Controller refresh supports a conventional single-PES launch without campaign
  receipts when its unique local input matches the retained accepted snapshot.
  It preserves dispatcher allocations, automatic orbital probes, profiled QC
  and sampling. A refresh inside a Slurm allocation refuses before mutation;
  batch-owned controllers restart by resubmitting their original PES batch job.
  Refresh does not create replacement Slurm controller allocations.

Validation uses local synthetic native records and simulated scheduler/process
boundaries, including complete refresh, interrupted staging, source adoption and
repeated invocation. This is not a native ethane calculation or a remote restart
verification. The private `composite_zador` branch and its three unpublished
patches were not available locally and are not represented as reviewed here.

The scoped regression suite passed 643 tests and 56 subtests, with one skipped
test. It was run from a clean archive of the published base with only this
change applied, excluding unrelated working-tree changes. An older dispatcher
test fixture was updated to include the workers inventory required by a real
coordination record; its review-epoch and allocation-preservation assertions
remain unchanged.

For the CH3S/CH3SS scout, use the ordinary in-place `run_refresh --wait-idle 600
--apply` after pulling the published commit. Do not delete scientific holds or
reopen failed saddle searches wholesale. The reported finite dissociation paths
and unsuccessful saddle attempts still require their native scientific evidence;
these handoff corrections do not establish barrierlessness.
