# Recovery and PES-wide composite reuse — 2026-10-21

This update addresses the supplied CH3S/CH3SS failure evidence. It does not
promise that every candidate saddle exists or that every electronic calculation
will converge. The code paths are generic; the cluster/run paths below are only
handoff arguments.

## Duplicate work and master behavior

Master's within-folder reaction filtering does not by itself prevent separate
well controllers from calculating the same saddle in opposite directions.
Traditional L3 work after outer-PES consolidation avoids that duplication.
Composite per-well DAGs therefore need shared ownership before L3 submission.

`anl_ts_reuse: {}` now enables conservative whole-graph sharing within one PES
for stable species, calculated CBH references, transition states and accepted
product complexes. It remains opt-in, with the old parameter name and
`--deduplicate-ts` CLI alias retained. `--deduplicate-anl` describes the expanded
scope. There is no sharing between independent PES directories.

* Saddles must have the same unordered endpoints/path class, state, L2 theory,
  requested composite settings, energy, ZPE and frequency spectrum, and match
  under element-preserving atom permutations and a proper rotation. Sharing
  products alone is insufficient. Reflected distinct chiral structures and
  isotope-bearing records are not automatically merged.
* Stable/CBH states additionally need matching chemical identity and spin.
  Different vdW discovery names can share a graph only after matching fragment
  identities and accepted complex structures. Distinct complex conformers remain
  separate. This shares QC; it does not invent kinetic degeneracies or remove
  distinct reaction paths.
* Existing comparisons use RMSD 0.0005 angstrom, electronic energy 1e-6 hartree,
  ZPE 1e-5 hartree and frequencies 1 cm^-1. Each positive comparison records its
  measurements and atom permutation. Ambiguous comparisons retain separate work.
* A serialized per-PES registry reserves even unique graph owners before graph
  creation. Later wells cannot demote a reserved owner. Initial migration prefers
  more completed tasks, then more submitted tasks. A duplicate graph is archived
  whole; components from independently optimized geometries are not spliced.
* Reverse discoveries finish their private graphs before waiting for/adopting
  shared graphs. Graph drive locks serialize adoption. Reference promotion and
  append-only network expansion preserve earlier ownership and graph revisions.
* Canonical graph paths reach profile assembly, VPT2, CBH reference models and
  final handoff auditing. A local consumed-reuse record preserves the actual
  mapping even when the global registry later grows.

These checks deliberately can miss duplicates when evidence is insufficient.
They do not consolidate independent ordinary L1/L2 searches or their reference
optimizations. The expensive composite calculation is the shared unit.

## Numerical repairs

* The NFS shared manifest lock now uses a readable/writable descriptor. One
  automatic migration retry recognizes the exact pre-QC `worker_load` EBADF
  traceback. It archives that failed attempt, verifies allocation termination,
  and retains unrelated successful/running work. Other EBADF errors are not
  assumed equivalent.
* Open-shell Gaussian diatomics have a bounded SCF fallback: detected symmetry,
  `SCF=Symm`, a fresh Huckel guess and CDIIS, preserving method, basis, grid,
  charge, multiplicity and requested density convergence. The supplied SH test
  converged in 18 cycles, while its fresh no-symmetry control failed. This
  addresses the observed symmetry-sensitive density oscillation rather than
  accepting the unconverged energy. Explicit nonstandard guesses/occupation
  controls are not replaced. `diatomic_scf_symmetry: false` opts out.
* The native SH force table was already in input coordinates despite Gaussian's
  preceding standard-orientation geometry. Recovery validates geometry without
  rotating those forces a second time. Polyatomic optimization and the existing
  VPT2 framework-symmetry recovery are unchanged.
* Optimizers now have two additional step chunks by default. Same-allocation
  continuation preserves the optimizer and learned Hessian. Failures from QC or
  invalid curvature do not become convergence successes.
* Old exhausted ANL/VRC optimizations can resume from a verified evaluated point
  within that budget. Sella's logged accepted trajectory index is used, avoiding
  later Hessian-probe frames. Fixed scan distances remain bound to the original
  staged geometry. A new allocation rebuilds its optimizer/Hessian; it is not a
  serialized optimizer restart. Scientific graph identity is unchanged.
* Ordinary unconstrained Sella/BFGS step failures get a verified saved-geometry
  continuation, using remaining chunks. The exact internal-coordinate ODE
  failure can retry from that point in Cartesian coordinates. Native theory,
  saddle order, force targets and frequency validation remain intact.
* The old failed SH worker can receive one migration attempt with the corrected
  calculator, independently of its timeout allowance. Ordinary retries require
  exact allocation/worker/seed provenance and use generation fencing and durable
  submission receipts. Missing evidence causes a review stop, never a guessed
  resubmission.
* Reference inputs no longer inherit the parent's stereochemical geometry or
  PES barrier reference. A live reference input is not overwritten during
  adoption. The CBH whole-molecule reference retains the requested spin/charge,
  fixing the supplied `[CH]` state conflict.

The two supplied TS failures with no qualifying imaginary mode remain rejected.
A larger optimization budget does not waive this scientific acceptance check.
Step exhaustion can still be a real failure after the bounded extensions.

The generic policy knob is `queue_live_policy.optimization_step_retries`
(integer 0–10, default 2). A live shared policy can be changed without replacing
controllers:

```bash
python -m kinbot.run_policy set /absolute/run --optimization-step-retries 2
```

Zero restores one chunk for future optimizer calls. Already running optimizer
calls and old pinned workers keep their existing behavior. Confirmed Slurm
walltime retries remain separate, retain their existing limits and increase the
allocation within the partition maximum. No native MRCC/CFOUR amplitude restart
is newly promised by this change.

## Barrier migration and retained policies

The supplied saved inputs still had a 1000 kcal/mol L1 threshold and no L2
threshold. The handoff explicitly sets a 400 kcal/mol L2 ceiling relative to the
initial well's electronic energy plus harmonic ZPE. Each later well derives its
own local ceiling from that same reference, avoiding accumulated offsets. L1
uses the existing screening margin. This is a conservative search ceiling, not
a proof of negligible kinetics at 5000 K. Defaults for other runs are unchanged.

Migration refuses before stopping controllers if an existing accepted TS is
above the proposed absolute ceiling: deleting accepted chemistry requires an
explicit network-pruning decision. It does not silently remove saved branches.

The existing scratch discovery/validation, exclusive allocations, shared pool
and exclusions, CPU/memory, walltime/timeout policy, MR active spaces, aDZ-F12
scan geometry, extrapolated correction settings, SO data and ROTD settings are
retained. Future unreviewed VRC channels remain held for review while independent
CBH work proceeds. Only the recognized numerical failures and proven duplicate
ANL graphs are selected; no blanket retry of scientific rejections is requested.

## Current in-place handoff

Run on the controller host. This replaces the Python controllers, adopts
running/pending ordinary L1/L2, nested CBH-reference and dispatcher allocations,
archives proven duplicate graphs after their owned jobs stop, and resumes the
same directory. Canonical completed/running work and existing VRC points stay
in place. Do not add `--replace-pending` or another VRC basis migration.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
git pull --ff-only origin composite
/home/hram/KinBot/.venv/bin/python -u -m kinbot.run_refresh \
  /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage \
  --deduplicate-anl --barrier-threshold-l2 400 --apply
BASH
```

Omitting `--apply` previews the current selection without stopping controllers
or changing jobs. The preview may create/validate a code snapshot. The actual
selection is recomputed after controller quiescence to catch completed/new jobs.
Repeated identical launched handoffs are no-ops; interrupted transactions must
be resumed with the same code/options. Ownership uncertainty or an unidentified
local child prevents the handoff rather than permitting duplicate submission.

```bash
watch -n 10 /home/hram/KinBot/.venv/bin/python -m kinbot.run_monitor \
  /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage
```

Check the reported retained IDs and duplicate cancellation selection against the
current queue. The restart is not confirmed on Blodgett until its new report is
available. Numerical recoveries occur when the relevant controller next revisits
that task; the full exploration barrier still controls new downstream work.

## Verification record

The final targeted regression run passed **592 tests**, with 1 skipped and
10 subtests passed. It covers dispatcher and ordinary-worker recovery, controller
replacement, shared ownership and adoption, final profile/VPT2 bookkeeping,
reference inputs, barrier migration, live policy, scratch, scheduler and worker
result validation. `git diff --check` also passed.

The final full suite reported **2468 passed, 11 failed, 7 skipped**, with
429 subtests passed. The remaining failures reproduce the previously identified
baseline issues: the external CFOUR fixture lacks geometry; sulfur-example and
CASPT2 tests retain old basis/reference expectations; MR orbital-probe mocks
expect the old working paths; and a species-routing test substitutes a
`StationaryPoint` without the fragment's `orbital_cube` attribute. This is not a
claim that the entire repository suite passes.

The Gaussian SH symmetry and force-frame evidence came from the supplied native
cluster check. No new Gaussian, Molpro, MRCC or CFOUR production calculation was
executed by this local validation. Cluster launch, task adoption and native
convergence still need confirmation from the resulting run records.
