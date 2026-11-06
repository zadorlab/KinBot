# Read-only native evidence for unresolved homolysis checks

`python -B -m kinbot.run_seed_audit RUN --homolysis-native --output-parent EXPORTS
--max-mb 16384` captures complete files from ordinary L1/L2 workers with unresolved
`homolysis_checks.json` entries. EXPORTS must be outside RUN.

Use `--ordinary-native` instead to include every ordinary PES worker, staged
reverse worker and combined MESS directory, even when no homolysis note exists
there. This includes complete HIR/conformer outputs, selected minimum evidence,
thermochemistry models and shared-work receipts needed to check population
counting and cross-worker reuse. The same read-only and completeness checks
apply. A terminal inconclusive homolysis entry is an evidence note, not an
export hold; the capture does not change that decision.

It includes native inputs and outputs, Gaussian/Sella logs, frequencies/Hessians,
current working checkpoints, trajectories, IRCs, fragment calculations, conformers, alternative
tight-saddle evidence, complete databases with sidecars, search and submission
histories, and existing ordinary retry archives. It includes all ordinary native
files in the affected worker because a different search can already contain the
relevant saddle or fragment evidence. It does not reset, ingest, submit, cancel,
retry, accept or reject calculations, and does not execute copied Python/pickle
files. Slurm accounting is queried in bounded batches for recorded
scan/TS/IRC job IDs; missing accounting is not a failed-calculation verdict.

ANL, VRC/RotdPy, raw scratch trees, periodic `.kinbot_checkpoints` generation
stores and raw RWF/WFU files are outside this L1/L2
capture. This is evidence for diagnosis, not a complete restart backup. The
manifest states this explicitly. Files larger than 256 MiB or exceeding the total
budget, unreadable files, symlinks and source changes are reported as capture
issues, never silently presented as complete evidence. A live run can therefore
yield an `incomplete` capture; it is not stopped to obtain a snapshot.

The export prints progress immediately, creates `summary.txt`,
`homolysis_inventory.json`, `scheduler.json`, `manifest.json` and a
`seed-evidence.tar.gz` archive, verifies the captured file hashes inside the
archive, and prints its SHA256. The summary preserves recorded search states;
it does not turn a naturally unsuccessful saddle search into a retry request.

Local regression checks cover full native-file inclusion and unchanged source
contents, archived failures and alternative saddle evidence, changes during a
live capture, and batched read-only accounting including unavailable responses.

## Large exports and terminal disconnects

Use the detached mode from the login host instead of keeping a foreground
collector attached to SSH:

```bash
python -u -B -m kinbot.run_seed_audit RUN \
  --ordinary-native --mess-handoff --max-mb 16384 \
  --output-parent EXPORTS --detach
```

This returns immediately with a `seed-capture-*` job directory, collector PID,
log path and status path. It starts **only** a read-only local collector. The
collector has its own process session, redirects all terminal descriptors, and
ignores SIGHUP. It survives terminal logout/hangup; this does not promise
survival of a host reboot, SIGKILL, or a site's policy that kills background
processes. No QC jobs or controllers are started, cancelled or changed.

On a new login, inspect the latest capture without needing old shell variables:

```bash
python -B -m kinbot.run_seed_audit --status EXPORTS
```

The same command accepts the exact `seed-capture-*` directory to inspect an older
capture. It prints saved status and the last log lines. `complete` means capture
and archive verification finished; `incomplete` means source/copy issues are
listed in the manifest; `failed` identifies a collector error. A saved `running`
state with an old timestamp is not proof that the collector remains alive.
`result.json` contains the archive location and SHA256. These are capture
verdicts, not scientific acceptance of any calculation.

Inventory, copying, source rechecks, compression and archive verification now
report phase changes and throttled progress (about every five seconds while
files advance). A blocked filesystem call can still delay the next update.
Overlapping case/worker scopes are walked once per inventory. A second inventory
and all source-hash rechecks remain necessary to detect live changes. Compression
uses gzip level 1, matching the existing orbital exporter. Complete native files
and all checksums are retained. The archive stays `.partial` until verification
passes; no partial archive is advertised as the final export.

Do not repeatedly launch captures while one is running: each invocation creates
a separate directory. Earlier exports and interrupted partials are retained.
Detached mode does not resume an older interrupted foreground capture. To keep
an inactive run's diagnostic snapshot stable, let capture finish before refreshing
its controllers.

## Periodic recovery generations and budget checks

The November 5 export exposed a collector defect: recursive selection included
`.kinbot_checkpoints/native/JOB/ATTEMPT/generation-*`. These periodic recovery
copies sorted before current native evidence, exhausted the 16 GiB budget and
generated millions of omission records. They are now excluded **before directory
traversal**, including checkpoint stores nested in HIR and reverse workers.
Their contents are not deleted or modified. Current working checkpoints,
checkpoint receipts, native outputs, databases, HIR/conformer files and archived
SCF/timeout/failed attempts remain selected. This audit is not a recovery backup.
`selection.json` and the final manifest list the excluded trees explicitly.

Native capture checks the selected total size before copying. If it exceeds
`--max-mb`, the collector stops with `selection.json` (counts, bytes and the
largest files), without copying a partial run or creating an archive. Directory
and file inventory limits also stop explicitly before copying, rather than
creating millions of omission records. Source mutation and hash checks still
apply after preflight. For an incomplete export, `inconclusive_checks` is null;
`captured_inconclusive_records` only counts records actually present in the
partial capture. Zero captured records never certifies no outstanding checks.

This correction changes only evidence collection. Updating the checkout is
sufficient to use it; no controller refresh or QC retry is required.
