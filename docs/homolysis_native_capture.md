# Read-only native evidence for unresolved homolysis checks

`python -B -m kinbot.run_seed_audit RUN --homolysis-native --output-parent EXPORTS
--max-mb 16384` captures complete files from ordinary L1/L2 workers with unresolved
`homolysis_checks.json` entries. EXPORTS must be outside RUN.

It includes native inputs and outputs, Gaussian/Sella logs, frequencies/Hessians,
checkpoints, trajectories, IRCs, fragment calculations, conformers, alternative
tight-saddle evidence, complete databases with sidecars, search and submission
histories, and existing ordinary retry archives. It includes all ordinary native
files in the affected worker because a different search can already contain the
relevant saddle or fragment evidence. It does not reset, ingest, submit, cancel,
retry, accept or reject calculations, and does not execute copied Python/pickle
files. Slurm accounting is queried in bounded batches for recorded
scan/TS/IRC job IDs; missing accounting is not a failed-calculation verdict.

ANL, VRC/RotdPy, raw scratch trees and raw RWF/WFU files are outside this L1/L2
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
