# Read-only QC decision evidence

Run `python -B -m kinbot.qc_decision_audit RUN --output EMPTY_DIRECTORY` to
capture evidence for acceptance and rejected-search review. The output must
be outside the calculation run. The command does not execute saved workers,
construct a QC driver, query Slurm, retry jobs or change search decisions.

The export includes:

- Search outcomes, theory/input settings, L1 pruning and endpoint-check records.
- Consistent SQLite backups including committed WAL data and full row histories.
  The text report displays the latest three rows per calculation; it never
  substitutes an older success for a newer failed row.
- Current result/exit/SCF/optimization-budget receipts, minimum-recovery attempt
  reports, and original input/output paths identified by SCF recovery receipts.
- Generated Python workers and Gaussian inputs, plus bounded native/controller
  log tails. Offsets, hashes and truncation are recorded explicitly.
- Generated PES-wide MESS inputs and the relevant installed source files at
  audit time. These files do not establish which code an earlier process loaded.

`summary.txt` lists search outcomes and consistency questions. `evidence.txt`
contains the structured report for transfer or pasting. The compressed archive
also includes the database snapshots, captured source files and native tails.

The preliminary checks replay finite geometry/energy/ZPE and signed-frequency
requirements for expected accepted TS and stable-species cache labels. Passing
these checks is **not** native scientific validation. A missing/error row can
instead reflect an alternate selected conformer, reuse mapping or changing
live run; it is a question to trace, not an automatic acceptance bug.

Rejection categories retain their evidence and do not establish that a physical
channel is absent. SCF failure, exhausted optimizer budgets, wrong stationary
point order or wrong IRC endpoints do not by themselves justify a retry. A
useful recovery proposal needs applicability and attempt-history review; only
an actual native retry can establish that it succeeds.

Each SQLite backup is consistent independently. The overall live-run export is
not an atomic snapshot; detected file/search changes are listed. Native tails
and older controller diagnostics can be incomplete. This audit cannot prove
that there are no acceptance bugs or decide counterfactual retry success.
