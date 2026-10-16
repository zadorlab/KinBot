# Deferred repository cleanup

**Inventory only. Do not execute cleanup until CH3S and CH3SS have completed
the requested end-to-end KinBot runs on Blodgett.** No files were removed or
moved for this inventory. Native inputs, outputs and orbital evidence must
remain recoverable and referenced from the final conclusions.

The [machine-readable inventory](REPOSITORY_CLEANUP_INVENTORY_2026-10-15.json)
covers every top-level entry except Git internals, does not follow symlinks,
and lists files of at least 5 MB outside the environment/package caches.
It found no read errors. Sizes are a snapshot while ORCA is running, not a
promise of reclaimable space. Parent Git tracking does not classify files
inside the `external/ROTD_py` submodule.

## Keep in the repository

| Material | Reason and action |
| --- | --- |
| `kinbot/`, `tests/`, packaging metadata, license, `graphics/` | Production source, regression coverage and distribution assets. Newly added `orbital_geometry.py`, `orbital_validation.py` and tests need review, not removal because they are untracked. |
| `external/ROTD_py`, `.gitmodules` | Required integration/submodule. Preserve its recorded revision and local changes. It is not disposable generated output. |
| General workflow docs, current example READMEs and JSON recipes | Retain a concise current entry point plus the exact configuration users need. Site settings and chemical choices belong in data rather than shared implementation. |
| `examples/anl/orbital_reviews/molpro_confirmation/` | Keep the finished compact confirmation recipes, provenance index and instructions. Large generated jobs and checkpoints belong in the external evidence archive. |
| Reusable preparation/analysis scripts | Keep general utilities with tests and documentation. Review useful local ORCA adapters for promotion into the general workflow; do not delete the only reproducible implementation. |
| Small native test fixtures | Keep deliberately chosen minimal regression fixtures; avoid replacing them with entire production output trees. |

## Archive outside normal source control, preserving references

| Material | Planned disposition and dependency |
| --- | --- |
| `local_runs/qchem_casscf_20261015/` (about 5 GiB at inventory) | Build and verify a production-evidence archive. Retain accepted native `.inp`/`.out`, parsed results, geometries, basis definitions, executable versions, solver/method settings, state/orbital diagnostics, required `.gbw`/`.wfu`/orbital exports, hashes and replay tools. Keep a compact manifest and stable archive references in the repo. |
| Superseded local campaigns and failed attempts | Archive separately from accepted evidence, with reason/status and links explaining any conclusion or recovery. Preserve original files where needed to reproduce root/solver/initialization issues. |
| `local_runs/.../preparation_history/` | Retain as history, including the withdrawn absolute-energy/full-space selection interpretation. Link the current selection record so history cannot be mistaken for current instructions. |
| `local_runs/.../reproducibility/*.tar.gz` | Compare manifest/hash contents with the final archive before identifying redundant snapshots. Do not assume two archives with similar names contain identical evidence. |
| `examples/anl/orbital_reviews/ch3s_ch3ss_orbital_review_20261015.tar.gz`, `ch3s_ch3ss_orbital_diagnostics_20261015.tar.gz`, `ch3s_ch3ss_potential_validation_20261015.tar.gz` | Retain original Blodgett exports in the evidence archive. Replace bulky working copies only after checked archive extraction and repaired documentation links. |
| `examples/anl/ch3s_anl1_f12/ch3s_tzf12_review_coz08qx7.tar.gz` | Review overlap with the other native exports. Archive rather than embed a temporary-named export beside production input. |
| `examples/anl/orbital_reviews/results/` and `combined_results/` | Preserve reproducible summaries/plots with provenance; clearly label them historical if superseded. Keep only the useful compact current plots in the main docs. |
| `docs/archive/`, dated readiness/restart/checkpoint documents | Consolidate the current instructions, retain historical records under a clear archive entry point, and check every statement against the final evidence. Do not silently erase previous numerical results or caveats. |

The source directory contains several preparation/restart helpers that were
created to recover this standalone study. Review each for reuse or archival;
the final user workflow must not depend on undocumented one-off steps. The
finished Molpro confirmation and general orbital-recon procedure replace the
old sequence of ad hoc restart instructions as the current entry points.

## Remove only after evidence and reproducibility checks

| Material | Condition before removal |
| --- | --- |
| ORCA transformed integrals and transient `*.tmp` files | Native process and all child processes have exited; required restart/orbital records were identified explicitly; accepted native outputs/results are archived. Never apply a broad suffix deletion to an active run. |
| Controller locks, PID files, transient logs | Controllers are inactive; archive execution/retry logs that establish provenance first. A stale PID is not authority to kill a reused process ID. |
| `__pycache__/`, `.pytest_cache/`, `.mpl-cache/` | Reproducible caches, already ignored. These can be regenerated after the requested end-to-end completion. |
| Duplicate generated decks/plots and superseded exports | Demonstrated duplicate content or a verified replay path exists, all documentation links resolve, and no current configuration references them. |
| Empty `conf/` | Check whether packaging or workflows expect it; remove only if unused. |

## Keep local, outside version control

| Material | Action |
| --- | --- |
| `.venv/` (about 1.05 GiB) | Keep while workflows use it. Export package/version information for reproducibility. Relocate/recreate only after the examples finish; it is already ignored. |
| `.conda-pkgs/` (about 1.09 GiB) | Package cache, already ignored. Review environment links and use the environment manager's supported cleanup after successful environment recreation; do not manually prune during validation. |
| Root `.kinbot_run.json` | Generated RDKit/stereochemistry run metadata. Retain it with the run it belongs to, then remove the misplaced root copy if no longer referenced. Do not confuse it with a source configuration. |
| Future generated native runs and large evidence exports | Add narrowly scoped ignore rules after deciding the archive/index layout. Do not ignore every JSON/output suffix: some are intentional source recipes or regression fixtures. |

## Completion procedure

1. Confirm both production workflows have completed and review their final
   thermochemistry, corrections, Monte Carlo uncertainty and MESS outputs.
2. Freeze the final evidence index, decisions, environment records and archive
   hashes; verify extraction and replay of the numerical analyses.
3. Check that retained documentation points to accessible native evidence,
   including rejected attempts needed to explain a choice.
4. Make a concrete keep/archive/remove manifest from this inventory, checking
   active processes and configuration references again at that time.
5. Carry out the agreed cleanup and rerun the appropriate source/docs checks.

This plan addresses deliverable 10 in
[the approved checklist](ORBITAL_RECON_DELIVERABLES.md). The final inventory
must be refreshed after native production because new evidence will exist.
