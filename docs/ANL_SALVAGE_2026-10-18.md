# Recovering the October 18 orphaned ANL runs

The export `r5-salvage-20261018-111955.zip` shows a geometry fallback policy,
not a submission race. Job 52909613 optimizes the QZ geometry of
`ts_857b1afac1911a45`; that graph's 6Z energy was not staged. Jobs 52909609
(`ts_f866e4cbbd8d38b5`) and 52909574 (`state_684eb53c4637f96a`) belong to
other structures. Their TZ geometries failed, their QZ geometries were
skipped, and their 6Z inputs used the explicitly declared L2 fallback.

## Required geometry and retry policy

Automatic ANL now defaults to `anl_require_l3_geometry: 1` and
`anl_l3_max_retries: 1`. These are general JSON input controls. A requested
L3 optimization gets one clean retry after an execution/scheduler failure;
each failed attempt is archived. Thresholds, method, basis and starting
geometry are not relaxed. Minimum recovery inside an optimization remains
separately controlled by `anl_minimum_recovery`. Preflight, parser and
provenance failures do not qualify for a blind retry.

Recipe single points wait for their declared L3 geometry. If the required
optimization still fails, no L2 geometry is silently substituted and no
completed ANL profile is assembled. Setting `anl_require_l3_geometry: 0`
explicitly restores the older recorded geometry fallback policy. A state
assigned only L2 or SP@L2 by the size policy still uses that explicitly
assigned tier; this option does not change the size policy.

VPT2 is a frequency-only calculation on the accepted matching L2 geometry
by design. Its valid results are retained. A diagnostic explicitly named
SP@L2 is also distinct from an ANL recipe single point.

In this export, 23 of 30 complete task records are eligible for retention;
seven completed single points used the fallback and are archived. Two more
native harmonic TZ outputs can be recovered by the corrected parser:
Molpro printed their imaginary modes in a separate unsigned table. The
modes are -1640.46 and -1640.44 cm-1, with eight real modes each. Acceptance
on Blodgett still verifies the full artifact hashes and geometry provenance.
The two incompatible live 6Z jobs are stopped only after a fresh Slurm
ownership/directory check. The live QZ optimization is retained if it is
still running. The tool discovers current job IDs; it does not hardcode
these historical IDs.

## Recovery and paths

`python -m kinbot.run_salvage OLD NEW --launch` works on the original Linux
controller host with the controllers stopped. It requires successful Slurm
queries, matching submitted-job records, input checksums and controller
locks. It never uses a user-wide kill or cancellation.

The directory is atomically exchanged with a symlink using Linux
`renameat2(RENAME_EXCHANGE)`. The actual tree moves to NEW, and OLD becomes
an alias to NEW. This keeps the same directory inodes, current working
directories and absolute paths usable by native jobs already running.
If a shared filesystem does not support that exchange, OLD stays physically
in place and NEW becomes its alias instead. No gap in the live pathname is
introduced. `salvage.json` records the layout and physical root.
**Keep both pathnames while this run is in use.** Copying a changing
native directory and changing its driver independently would risk lost
writes and duplicate submissions.

Failed/replaced task directories, old launch inputs and controller logs are
stored in a dated sibling `OLD.salvage_archive.*`. Accepted native artifacts
remain unchanged. The new drivers use a hash-checked source snapshot under
`NEW/.salvage_code`, preventing later repository pulls from changing their
Python modules midway through the run. The same recovery command is
restartable after a scheduler outage or interrupted preparation. An unknown
live controller or ambiguous submission stops it rather than launching a
duplicate. Both PES controllers run as detached login-node processes.

The first October 18 recovery stopped after cancellation and directory
exchange, before launching either replacement controller. Its graph repair
mistook a retained `submitted` record for a live allocation: the dispatcher
normally retains that record for 60 seconds after queue exit while waiting
for `execution.json`, which a cancelled job may never write. Recovery now
checks Slurm again before archiving an invalid-geometry task in that state.
Only confirmed-absent allocations can be retired; a live allocation or
failed scheduler query still blocks archival. Valid jobs keep the normal
publication wait. The archived state retains the original job ID, and the
graph's `salvage.json` records `retired_jobs`. Rerun the same command below
to resume the partial r5/r6 layout; do not remove either pathname or run a
general cleanup. Regression coverage exercises this interrupted sequence,
preservation of the valid live task's hashes, and launching both controllers
once after preparation completes.

The old input's scientific settings and node exclusions are retained,
except for the explicitly requested SCF/CC controls, mandatory L3 policy,
and Molpro pivot backend. Accepted earlier results with tighter thresholds
are kept only if they use the correct geometry. Recovery does not prune an
already explored reaction network or share references between the two PES
runs. Their existing five-node ANL and five-node VRC branch budgets remain;
adopted jobs count toward those budgets.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
git pull --ff-only origin composite
KB_OLD=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016_r5_dft
KB_NEW=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage
/home/hram/KinBot/.venv/bin/python -u -m kinbot.run_salvage "$KB_OLD" "$KB_NEW" --launch
BASH
```

Monitor both cases, including adopted jobs whose Slurm working directory
still has the old spelling:

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python -m kinbot.run_monitor /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage --lines 3'
```

Use `--details` for full metadata or `--follow` instead of `watch` to follow
new native/controller log lines. The monitor resolves aliases when matching
Slurm directories. An empty queue is never a success verdict.

## Other repairs and their validation limits

- The geometry audit found a real CC threshold translation bug. Generated
  `THRESH,COEFF=1e-10` printed **THRVAR=1e-20** in Molpro 2024.1; older
  defaulted jobs also printed 1e-18. Neither is the intended 1e-10 variance.
  New and restaged CC commands use **`THRDEN=1e-7,THRVAR=1e-10`** directly,
  including conventional/F12 single points and numerical geometry/frequency
  evaluations. This preserves the intended amplitude criterion and the
  requested energy tolerance; no `NOCHECK` is used. Already accepted results
  and the useful live QZ job retain their original inputs. The native output
  establishes the extra squaring; Molpro documents the direct options in
  [CCSD](https://www.molpro.net/manual/doku.php?id=the_closed_shell_ccsd_program)
  and [the open-shell interface](https://www.molpro.net/manual/doku.php?id=open-shell_coupled_cluster_theories).
  An [official Molpro response](https://www.molpro.net/pipermail/molpro-user/2004-August/000976.html)
  also demonstrates the squaring of `THRESH,COEFF`.
- Molpro can generate the fragment HOMO cube without Gaussian checkpoints,
  `formchk` or `cubegen`. Select `vrc_tst_pivot_backend: "molpro"`,
  `vrc_tst_pivot_method: "UKS/B3LYP3"`, and
  `vrc_tst_pivot_basis: "cc-pVDZ"`. These auxiliary orbitals orient VRC
  pivots; they do not replace the CASPT2-F12 interaction energies. The native
  input/output, final geometry, forces and cube are retained and checked.
  Molpro's B3LYP3 uses the VWN3 convention associated with Gaussian B3LYP;
  Molpro's B3LYP name denotes VWN5. See the
  [Molpro DFT manual](https://www.molpro.net/manual/doku.php?id=the_density_functional_program),
  [CUBE manual](https://www.molpro.net/manual/doku.php?id=dump_density_or_orbital_values_cube),
  and [OPTG manual](https://www.molpro.net/manual/doku.php?id=geometry_optimization_optg).
- Nonlinear VPT2 retains `Symmetry=None`, including the existing workaround.
  Only geometrically linear VPT2 inputs use `NoSymm`, which retains the
  linear point-group classification needed by Link 716. The native H2
  failure explicitly said `Molecule is non-linear but NTrRo=5`.
- A separate, bounded recovery handles the observed H2S Link 717 failure
  with an **empty** old framework group under `Symmetry=None`. It retries
  with the new group printed by Gaussian (`CS` in this output), retaining
  the same geometry and theory. Other named-group mismatches retain the
  established `Symmetry=None` recovery. See the licensed distributor's
  [Gaussian symmetry documentation](https://www.conflex.co.jp/gaussian_support/symmetry.php).
- The seven interrupted Molpro tasks in the export lack a specific native
  error sentence: MPI abort/termination is insufficient to identify an
  electronic or hardware cause. Their archived diagnostics are retained;
  the recovery retries them with the requested numerical settings. It does
  not claim that their underlying native failure has been established.

The five failed TZ geometry tasks were SH (`state_20691a101e49c36b`,
force step 12), CH4 (`state_684eb53c4637f96a`, step 5), CH3S
(`state_9b7047433cae3eb9`, step 1), `C=[SH]`
(`state_074fcc5ef5188035`, step 4), and TS `ts_f866e4cbbd8d38b5`
(step 2). All aborted in a native force evaluation, rather than reaching
the geometry iteration limit. SH's BFGS bond optimization had converged to
0.000029 eV/angstrom; the abort happened during the subsequent finite
difference curvature check. The rejected CH3SS SH **L2** reference is a
different failure: its old SCF density plateau was 8.95e-10 against a
1e-10 requirement. Its clean retry uses the requested 1e-8 SCF setting.
Earlier accepted L2 geometries are retained; no imaginary mode is ignored
to turn a failed minimum check into success.

The failed L3 outputs end after successful reference-point energies or
partway through a displaced-point calculation. Their MPI aborts do not
identify OOM, disk exhaustion, or an electronic cause. One previous Slurm
guard attempted an I/O retry; its classification alone is not proof of the
final failure's cause. New force/external failures preserve bounded native
and scratch rank-log tails plus scratch capacity in `*.failure.json`
before private scratch cleanup. This improves the evidence available if a
bounded retry still fails. The unnecessarily strict CC variance is fixed,
but it cannot honestly be blamed for every truncated MPI failure.

Regression tests cover geometry invalidation, retention of live task hashes,
no duplicate submission, interrupted migration rollback, bounded L3 retry,
VPT2 routing, separate imaginary tables, both cube layouts and pivot result
validation, CC threshold migration and failure-diagnostic retention. The
Molpro pivot and scoped VPT2 retries still require native
cluster execution; local tests cannot establish success of licensed programs.

## Resource interpretation

The three reported Slurm allocations reserve exclusive 96-core nodes with
126000 MB RAM, but their native commands use **8 Molpro MPI ranks**, one
OpenMP thread each. Their stack is **1473 million eight-byte words per
rank**, about **11.8 GB/rank** or **94.3 GB total**, with remaining memory
for overhead. `MIN_MEMORY=0` requests all node memory; it is not zero use.

Increasing to 12 ranks under the same budget would give about 915 MW/rank,
below the declared 1024 MW minimum. More ranks can increase distributed
integral/scratch traffic and are not automatically faster. Keep the current
rank count until a representative benchmark supports changing it. `sstat`
MaxRSS is a per-task maximum, not necessarily aggregate node use. Advancing
CPU counters alone cannot distinguish useful computation from MPI polling.
The exported scratch paths are under `/home/hram/.cache/kinbot/molpro`;
check their filesystem topology before attributing slow integral sorting
to memory or proposing a node-local scratch change. See
[Molpro memory](https://www.molpro.net/manual/doku.php?id=general_program_structure)
and [Slurm sstat](https://slurm.schedmd.com/sstat.html).
