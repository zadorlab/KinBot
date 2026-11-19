# Evidence reviewed on 19 November 2026

Read the [scientific review](../../../../../../docs/R8_FINAL_EVIDENCE_REVIEW_20261119.md)
and [operational handoff](RECOVERY_HANDOFF.md).

| Artifact | Contents |
|---|---|
| `CALCULATION_DECISIONS.md` | Current B2PLYP judgments, retired retries and necessity gates |
| `prior_l2_connectivity.json` | Earlier converged Sella path and endpoint-frequency receipts for reuse |
| `task_inventory.tsv` | All 245 native tasks, status, job ID and native cause classes |
| `full_network_inventory.json`, `.tsv` | All 105 preserved full-network connections, exact endpoint identities and representations |
| `screening_energy_audit.json` | Every one of the 24 L1 exclusions, corresponding L2 energies or explicit missing result |
| `cc_comparison.json` | Same-convention T/Q plus DZ ZPE energies and atom-permuted proper-rotation SI geometry comparisons |
| `paper_si_checks.json` | SI energy ledger and all four Eckart depth checks |
| `path_geometry_checks.json` | Raw versus postoptimized geometries, distances and prior-geometry comparisons |
| `remaining_l2_database_audit.json` | Read-only extraction of latest path/endpoint database rows for all 44 cases |
| `remaining_endpoint_matches.json` | Proper-rotation, element-preserving permutation comparisons to prior geometries |
| `verification.json` | Snapshot identity, separate verification checks and regression result |
| `evidence_summary.png` | Native TS5 separation/recapture, TS2 gradients and remaining L2 coverage |

The received archive and raw logs remain in the local evidence store rather
than being duplicated here. On the review workstation, the extracted files are
under `local_runs/r8_validation_review_20261119/evidence`; 57 XYZ inspection
artifacts are under the sibling `geometries` directory. Raw native geometries
on Blodgett remain under `validation_20261114/native/<graph>/tasks/<task>`;
the separate remaining-L2 geometries are stored in each case's `kinbot.db`.
XYZ comments distinguish path endpoints, optimized structures and unvalidated
database endpoints. File names preserve graph/task/section provenance.

Network `scientific_status` values are requirements for binding a serialized
edge to accepted native evidence. They do not discard the qualified positive
evidence described in the scientific review or certify an edge from its label.
