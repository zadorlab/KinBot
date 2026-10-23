# General Gaussian VPT2 framework recovery

All Gaussian VPT2 tasks use the same dispatcher procedure, selected by backend,
parser, native diagnostic and attempt history. No molecule, production path or
paper result determines a retry. The working directive is recorded in
[`AGENTS.md`](../AGENTS.md).

`Symmetry=None` and `Symmetry=(PG=C1)` are different processing strategies.
Gaussian documents None as bypassing point-group assignment and symmetry
processing; PG limits the symmetry used. NoSymm preserves group detection while
preventing reorientation. Loose relaxes initial symmetry detection. See the
[Gaussian symmetry manual hosted by its distributor](https://www.conflex.co.jp/gaussian_support/symmetry.php).
None and native-reported PG retries were already supported; the missing piece
was a shared no-repeat history including the general C1 alternative after the
native-group strategies fail. A previous success motivates an applicable retry,
not a guarantee that Gaussian will complete it for every geometry.

The procedure requires Gaussian's explicit framework-group inconsistency and
error termination. A segmentation trace alone is insufficient. Native success
is reconciled through the existing result parser rather than resubmitted.

For nonlinear geometries, the planner considers:

1. Symmetry=None, if untried.
2. The new group reported by Gaussian after an empty framework under None, or
   after a reassignment whose old group matches the current requested group.
   At most two distinct non-C1 groups are attempted.
3. One Loose retry when the initial full point group and the mass-based framework
   explicitly corroborate the disagreement with the requested group.
4. PG=C1, if still untried. This is the general minimal-symmetry cap, not a
   molecule-specific assignment.

The planner reevaluates applicability after each native failure. No setting is
repeated. It reads current receipts and same-geometry/theory task archives,
including operator-requested retries, and migrates the earlier recovery receipts.
Linear geometries retain the established NoSymm handling instead of the generic
None/C1 fallback. Unsupported user symmetry controls require review.

Each retry retains geometry, method, basis, charge, multiplicity, frequency
options, and result validation. The graph owner verifies allocation and artifact
identity, archives the failed attempt, and stages its replacement. Completed and
running siblings are preserved. An existing manual fallback keeps its original
scientific binding. Unknown errors or exhaustion leave the task failed for
review; they do not manufacture an anharmonic correction or silently accept a
harmonic substitute.

## Loading the procedure

Existing Python controllers retain their loaded code. One controller refresh
loads this extension; subsequent applicable retries are automatic and require
neither another refresh nor a separate species-specific request:

```bash
python -m kinbot.run_refresh /path/to/run --apply
```

Use the run's Python environment and updated checkout. The default refresh
retains native allocations, pending allocations and verified sampling owners;
do not add `--replace-pending` for this update. Already queued live controls
remain owned by their drivers.

The `recovery` field in `workflow.json` records the chosen setting, observed
framework mismatch and complete `attempted_symmetry_settings` history. A staged
or submitted retry is not evidence of native success; completion still requires
Gaussian's result to pass the unchanged parser.

Regression coverage exercises varied native point groups, the full sequence
through exhaustion, legacy and manual histories, linear geometries, route-only
edits, allocation fences, interrupted archival and preservation of other tasks.
Native Gaussian convergence must be checked on the calculation host.
