# Cross-worker minimum, L2 and HIR calculation reuse

New Gaussian/Sella `hir/`, `conf/` and minimum `_well_high` calculations within the same shared
resource-policy root now select a single producer. This includes separate PES
workers and reverse-search workers. No species identifiers or example paths
are special-cased. Existing local attempts remain authoritative; activation
does not enroll, cancel, rewrite, or invalidate old running/completed/failed jobs.

A match requires the same normalized worker program (theory, spin, charge,
optimizer, frequency controls, constraints and isotope definitions), runtime,
QC command, environment and custom queue template. Resource counts and output
labels are not scientific differences. Element-preserving atom permutations are
tested with a **proper rotation**, never a reflection. Constrained dihedral atoms
must map under the very same permutation. TS requests retain the previous
same-order rule. The usual maximum coordinate difference is 1e-6 Angstrom,
which only admits coordinate serialization noise.
External fields, embedding, externally supplied guesses and general basis files
are outside this contract and continue through normal submission. A missed cache
match never rejects a valid calculation.

## Independently converged minima

The selected L1 result can now establish that independent conformer searches
ended at the same minimum, so their new L2 refinements share one producer. The
selected L2 results similarly support sharing new HIR points. Evidence requires:

- A normal selected database record, its original optimization script, and a
  complete native Gaussian frequency spectrum bound to that same geometry.
- The same calculation program/theory, charge and multiplicity, a positive
  minimum spectrum, and the same atomic masses under the proposed mapping.
- Maximum mapped coordinate difference at most 0.0005 Angstrom, electronic
  energy difference at most 1e-6 Hartree, ZPE difference at most 1e-5 Hartree,
  and corresponding sorted frequencies within 1 cm^-1. These are identity
  checks within QC resolution, not relaxed QC acceptance/convergence criteria.

Nondefault isotope models currently opt out because these worker templates do
not carry a complete isotope contract. Missing evidence or different minima
leave the normal independent calculation in place. The original L1 results are
retained. No molecule name, connectivity identifier or formula alone establishes
minimum identity; wells, channels and population factors are not merged here.

Atom ordering can also change which terminal atoms KinBot chooses for the same
rotor axis. Before any local scan is staged, verified identical minima agree on
one mapped set of terminal atoms, preserving the local axis set and rotor
numbering. Their point geometry and full constrained-calculation definition
must still match. This prevents falsely aliasing different constrained scans.
Any existing local scan script, allocation, artifact or database row prevents
this change; running and previous scans keep their definitions. Different grids,
constraints or relaxation settings are not silently declared interchangeable.

`*.minimum_reuse.json` binds the source evidence to the exact request script.
`hir/*.minimum_mapping.json` records the original/mapped dihedrals and measured
minimum comparison. The pool's `minima/` registry retains native reference logs
and scripts. Matching and producer selection are serialized across controllers.

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
record, with SHA256 verification. Native outputs remain unmodified. For a
permuted consumer, the database geometry uses the producer's Cartesian frame
with the consumer's atom order. The Hessian follows the same 3-by-3 block
permutation, H_consumer = P H_native P^T. Native Hessian/mode readers consult the
recorded native frame before applying that mapping; formatted-checkpoint
conversion is a local derivative, not another QC calculation. Conventional
result publication retains this provenance too. HIR ingestion, mode projection,
symmetry and thermal-model construction run through their normal local drivers
using the mapped observations, rather than copying population factors.
Consumer database publication comes last and is repeatable
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

Additional regressions run the profiled L1-result -> L2 -> HIR -> projection and
publication path with independent, permuted/rotated minima and simulated native
handoffs. They check one L2 allocation and one set of HIR allocations, invariant
quadratic Hessian energies, mapped native vectors and equal projected spectra.
Full refresh tests retain the live source and unrelated sampling, adopt its
completion after restart, and repeat the handoff without another allocation.
Receipts created by the previous reuse implementation remain readable after
refresh. These tests do not assert that the two reported conformer jobs actually
converged to the same minimum; that still depends on their final native results.
