# Cross-worker HIR and conformer calculation reuse

New Gaussian/Sella `hir/` and `conf/` calculations within the same shared
resource-policy root now select a single producer. This includes separate PES
workers and reverse-search workers. No species identifiers or example paths
are special-cased. Existing local attempts remain authoritative; activation
does not enroll, cancel, rewrite, or invalidate old running/completed/failed jobs.

A match requires the same normalized worker program (theory, spin, charge,
optimizer, frequency controls, constraints and isotope definitions), runtime,
QC command, environment and custom queue template. Resource counts and output
labels are not scientific differences. Ordered atoms must match and the maximum
Cartesian difference after a **proper rotation** must be at most 1e-6 Angstrom.
This tolerance only admits coordinate serialization noise. No reflection,
atom permutation, conformer equivalence or population merging is inferred.
External fields, embedding, externally supplied guesses and general basis files
are outside this contract and continue through normal submission. A missed cache
match never rejects a valid calculation.

The first worker owns the Slurm allocation; other workers retain durable consumer
receipts without claiming its job ID. Registry locks serialize selection and are
released before waiting for shared-pool admission. The owning driver's normal
failure/recovery handling remains authoritative. A final unsuccessful source is
an unsuccessful attempt, not proof that a reaction channel does not exist.
Changed source definitions or independent replacement allocations release the
association, allowing normal independent handling; they do not supply an old
result under a new definition.

Successful adoption requires the matching generation's zero exit marker, the
completed optimizer handoff, normal scientific result, finite geometry/energy,
and native Gaussian normal termination. Conformers additionally need finite
frequencies and ZPE. Exit zero by itself is insufficient. Immutable bundles retain
native inputs, outputs, checkpoints, optimizer trajectory and the complete result
record, with SHA256 verification. The native coordinate frame, Hessian and mode
vectors stay together. Consumer database publication comes last and is repeatable
after interruption. There is no synthetic Gaussian output or consumer Slurm ID.

`.kinbot_qc_reuse/` holds the pool's entries and bundles. Each logical job has a
`.qc_reuse_producer.json` or `.qc_reuse_consumer.json` receipt. The consumer row
records the original job, allocation, source row (when available) and shared token.
Existing controller refresh preserves these records. If interruption occurs
between submission and local registration, consumers consult durable shared-pool
admission receipts. An ambiguous admission stays pending; it never authorizes a
second allocation. A producer interrupted before submission resumes through its
owning controller. No independent takeover is attempted.

## Captured regression evidence

The reported HIR jobs 52954924 and 52954919 have identical B2PLYP/cc-pVTZ/GD3BJ
settings and the same fixed dihedral. Their proper-rotation RMSD is
1.3779413778534112e-7 Angstrom; maximum displacement is 2.1744087189192038e-7
Angstrom. They are duplicate starting calculations. Both exit markers were zero;
the supplied comparison did not contain their scientific result rows, so this is
not a separate scientific acceptance verdict.

Conformer jobs 52954921 and 52954918 share B3LYP/cc-pVDZ/GD3BJ settings but differ
in atom order and starting conformation. Exhaustive element-preserving mappings
give a minimum proper-rotation RMSD of 0.8801379054790381 Angstrom. They must remain
separate searches. Their reported result receipts were normal; final-structure
equivalence was not included in this input comparison.

`tests/data/qc_reuse_pairs.json` preserves the supplied comparison evidence.
Tests render the actual worker templates, exercise concurrent producer selection,
normal/error handoffs, immutable-copy verification, interruption/repeated adoption,
profile routing, actual submission orchestration with a simulated scheduler, and
the full controller refresh with retained native jobs and sampling. These are
local lifecycle/regression tests, not newly run Gaussian calculations on Blodgett.
