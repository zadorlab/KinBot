# Gaussian dependency checkpoint handoff

Dispatcher Gaussian tasks may name a completed dependency's declared `.chk`
artifact using an absolute `%oldchk` path. KinBot verifies the workflow,
dependency status, staged inputs and successful execution's artifact hashes,
then copies that checkpoint into the child's local scratch workspace. A copy
hash check closes the source-read race. Arbitrary absolute paths, symlinks,
unaccepted producers and changed checkpoints remain errors.

Only the executed copy of the input receives a local checkpoint filename. The
original prepared deck, workflow, geometry, theory and state definitions remain
unchanged. The executed deck (`INPUT.com.native`), imported checkpoint and source
hashes are retained in execution evidence. The parent checkpoint is never opened
for writing by the child. Normal native success and scientific validation gates
still apply.

For existing failures with exactly the former pre-execution path rejection:

```bash
python -B -m kinbot.anl.checkpoint_inputs /path/to/graph-or-panel
python -B -m kinbot.anl.checkpoint_inputs /path/to/graph-or-panel --apply
```

The first command previews. The second queues generation-bound requests through
existing live controls. It does not edit the owner's state or submit/cancel jobs.
The owner archives the failed attempt and affected blocked descendants before
restaging. There is one repair per task; repeating the command preserves pending
requests, replacement generations and completed work. Native numerical failures
are not selected. The live-control root must already be enabled.

The frozen example panel's original owner did not poll live controls. Its new
`examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/resume.py` does, while
retaining `campaign.py` and its scientific fingerprints. Its `refresh` command
uses existing PID-safe pause/recheck/resume machinery and the pool admission
lock. Incomplete staging/submission or a scheduler outage prevents replacement.
Only the verified local panel owner is retired; all QC allocations are adopted.
Controller handoff receipts and previous native attempts remain on disk. The
separate saddle-recheck controller and PES controllers are outside this scope.

Tests execute a stub native program through the actual worker and scratch path,
then exercise live requests, archival, replacement-owner submission/adoption,
publication interruption, queue outages, stale requests, bounded repetition and
handoff races. These are local simulated-site tests, not a successful Gaussian
calculation or a verified Blodgett handoff.
