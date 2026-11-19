# Native reconciliation update — 16 November 2026

**19 November necessity revision:** [current decisions](review_20261119/CALCULATION_DECISIONS.md)
select four B2PLYP inner-path scans and one missing retained reverse IRC.
Failed optional CC checks do not block the available B2PLYP judgments. Earlier
retry/release commands below are historical.

**Current status, 19 November:** the native panel is finished and reviewed;
scientific closure remains open. Follow the
[evidence review and closure sequence](../../../../../docs/R8_FINAL_EVIDENCE_REVIEW_20261119.md).
The sampling/refresh commands below document earlier waves, not a new launch
recommendation. Preserve their decisions and attempt history.

## Finalization objective — 19 November 2026

The completed panel must resolve the paper agreements/disagreements and justify
the `_r9` model, following the current
[R9 acceptance criteria](../../../../../docs/r9_readiness.md). Native completion
and collection do not close that review. The latest export has been reviewed;
the separate 44-case L2 campaign remains incomplete in the captured evidence.

In particular, audit every candidate dissociation's actual bottlenecks before
carrying its scouting PST classification into R9. Beta-scission connectivity,
an optimized fragment endpoint, or an inherited `nobar` label cannot establish
a barrierless treatment. Resolve conventional inner saddles and loose outer
exits separately, then reassess affected branching/pruning and the generated
MESS network. The final deliverables are the evidence-linked paper verdicts,
corrected network, exact R9 channel/treatment/reuse plan and remaining blockers.

## Reviewed sampling and pool utilization

The supplied inventory shows a 20-job shared ceiling, 15 excluded nodes, and
21 eligible physical nodes (14 idle, six allocated, one mixed). The duplicated
partition rows refer to the same hardware. One tracked submission belongs to
the current live TS3 probe; it is retained by normal admission reconciliation.

The original blanket IRC hold is the main source of idle capacity. All 24 held
tasks are original, unattempted `forward_5`/`reverse_5` branches. Their parent
harmonic spectra each contain one imaginary vibrational mode. The policy for
this wave is explicit:

| Evidence in the same-direction step-10 log | Decision for staged step-5 |
|---|---|
| Failed at the native 100-cycle corrector limit | One half-step test, same spatial extent and convergence settings |
| Multi-point path followed by fragment recapture during endpoint optimization | One step-sensitivity test; inspect raw path separately from optimized endpoint |
| One-point gradient stop (both rxn_12_2 directions) | Retain hold; use the running departure investigation |
| Another failure, live source, changed evidence, repeated settings or altered task | Hold or refuse preparation with the specific reason |

This permits up to 22 new tasks, subject to the complete native logs on Blodgett.
The decision binds the source attempts, accepted Hessian, candidate input,
geometry and task definition. It changes admission only. The owning driver
retains all existing jobs and backfills free slots, up to the shared ceiling.
`kinbot.run_capacity` derives that ceiling from unique healthy nodes that meet
the staged jobs' cores/memory after exclusions. Its atomic policy update
preserves exclusions and retry controls, and refuses a concurrent policy change.
`sampling_capacity.json` records the applied inventory; the existing policy
history records the before/after values. A repeated handoff is idempotent.

The half-step tests retain MaxCycle=100 and the original L2 model. Numerical
frequency force tails are not used as independent stationary-gradient evidence.
Path connectivity, stationary-point identity, CC benchmarks and MESS publication
remain subject to the existing scientific review.

## Follow-up after the positive-time probes

The supplied native evidence is recorded in `native_followup_20261116.json`.

| Observation | Scientific interpretation | Next bounded action |
|---|---|---|
| All four positive-time DVV jobs executed, each stopping at one point; reaction coordinate 0.002–0.004 | The input error is fixed; these traces still do not establish connectivity | One rxn_12_2 forward departure diagnostic, capped at 64 points |
| TS2 tight central gradient: maximum component 0.027348837, norm 0.043571151 Eh/bohr; original frequency gradient norm 0.043580 | The rejected CC geometry is nonstationary; its five imaginary modes cannot establish a higher-order stationary point | Preserve the two running Hessian checks; jointly review their step dependence and the original optimization stop before selecting a geometry refinement |
| TS2 original QSD ended at iteration 21, final printed GRADMAX 0.03042668; requested tolerance 1e-5 | `END OF GEOMETRY OPTIMIZATION` did not mean convergence | Shared harmonic parsing now rejects the native final-frequency nonzero-gradient warning |
| TS3 shifted ROHF stagnated at 60 iterations; native executable is Molpro 2024.1 | A second unsuccessful electronic solve supplies no curvature result | One electronic-only quadratic ROHF trial, reusing the established orbital-reconnaissance determinant renderer |

Gaussian's [Overlay 1 controls](https://gaussian.com/overlay1/#iop_(1/108))
document a separate DVV gradient stopping test. The new departure task retains
the accepted parent Hessian, direction, 0.10 fs time step, electronic settings
and force–velocity angle test; `IOp(1/108=-1)` disables the gradient stopping test
for this **64-point diagnostic only**. It exports a raw log/checkpoint for review,
without endpoint optimization or scientific acceptance. Reaching the point cap
is not proof of a minimum. Its trace must demonstrate departure from the TS
region before longer paths or additional branches are selected. This explicitly
changes a stopping test, rather than tightening VeryTight again or extending the
corrector iteration ceiling.

The TS3 trial uses `rhf,so,maxit=60,energy=1.d-12` with explicit electron/spin
counts and C1 occupations matching the failed source: `wf,41,1,1;occ,21;closed,20`.
The existing ROHF-UCCSD(T)-F12b(T*) model, basis, geometry and thresholds remain.
Molpro's [SCF documentation](https://www.molpro.net/manual/doku.php?id=the_scf_program)
supports this solver for the observed version. The task exports orbitals and
electronic energies for comparison against the intended reference. It has no
derivative children: convergence alone cannot identify the intended orbital
solution. This recipe requires native version/occupation evidence and both
exhausted default/shifted attempts; it cannot apply recursively.

Two immutable additions, `reconciliation_electronic_v3.json` and
`reconciliation_departure_v3.json`, bind the new tasks to prior decisions and
source hashes. The existing owner applies them under its lock. The running TS1
pc150, TS2 Hessians, TS5 paths and any other allocated work continue to be
polled/adopted; unsent original IRC grid tasks remain held. Repeating the handoff
reuses the same decisions and task IDs. No new paper–R8 qualitative disagreement
has been validated by these outputs.

After pulling the published commit, run:

```bash
bash examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/resume_followup_probes.sh
```

Read-only status and evidence review:

```bash
.venv/bin/python -B examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/watch.py \
  "$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114/native" --interval 60
```

```bash
.venv/bin/python -B examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/reconcile.py review \
  "$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114/native"
```

`followup_validation.json` records local tests and the separate verification
pass. Simulated-site tests cover the complete combined decision/owner/worker
handoff, positive and failed native receipts, retained source hashes and live
jobs, held and already-allocated sampling, repeated invocation, wrong ancestry,
changed geometry/state, and transient submission/staging/archive/scheduler
races. These tests validate orchestration; the new native probes must still run
on Blodgett.

The objective remains a defensible comparison of the retained R8 connections,
barrier ordering and exclusions against the paper. Finishing jobs is not that
comparison. This update adds a finite diagnostic wave and leaves scientific
approval false. The frozen generator, source evidence and requested graphs
are unchanged. Original calculations remain available, including failures.

## What the completed paths actually establish

The user supplied job-tagged native excerpts and log hashes. The transcribed
geometries and their comparisons are in `native_review_20261116.json`;
the full logs remain on Blodgett. These are not locally executed calculations.

* **Native rxn_12_2:** both directions stopped at one point, with Gaussian's
  `PES minimum detected` message and gradient magnitudes 0.0001286 and
  0.0001123. The parent Hessian has one imaginary frequency, -178.3817 cm^-1.
  The final optimized endpoints have twelve positive vibrational modes, but
  that optimization followed the very short IRC. The forward S–S distance
  changed from 3.041 to 2.080 Å during optimization. Its final geometry matches
  the previously observed W2 conformer `paper_ts2_R` to 0.000104 Å RMSD.
  The reverse final geometry matches the previous `native_ch3ss_rxn_12_2_R`
  complex-B endpoint to 0.000064 Å. These are useful basin observations;
  a continuous W2/complex-B connection still needs independent path evidence.
  This does not validate the old low-barrier W2/W3 assignment.
* **Paper TS3:** both completed reverse branches also contain only one IRC
  point before optimization. Their positive-frequency final minima do not
  establish a full, step-converged path. The two directions and the original
  reactive mode still need joint review.
* **Paper TS5:** the 22-point reverse IRC ends with separated H. H(3) is
  3.136 Å from C and 3.665 Å from S(4), while the following optimization forms
  a 1.346 Å S(4)–H(3) bond. **The earlier inference of an H-transfer disagreement
  based on that optimized endpoint is withdrawn.** The raw path is consistent
  with the paper's H-loss direction; the detached fragment, forward direction,
  path-step dependence and matched asymptote energies still require validation.
  Its final captured minimum matches the earlier `paper_ts5_R` within 0.000112 Å.

RMSDs above minimize over permutations within each element and proper rotations,
without reflection. They compare printed geometries, not reactive trajectories.
Forward/reverse signs depend on the native imaginary-mode phase; compare the
mapped pair of endpoints rather than assuming that direction names agree across
programs or Hessians.
The production helper's covalent-radius component screen only flags possible
recapture; it does not assign a chemical channel or fragment electronic state.

## Finite native checks in this update

Five raw Gaussian path probes use the original accepted parent checkpoint,
charge, spin, B2PLYP-D3(BJ)/cc-pVTZ, grid and SCF settings. **Correction:** DVV
uses time steps, so the earlier claim of preserving spatial step/path length
was incorrect. The original DVV inputs carried negative EulerPC spatial steps
and were rejected before propagation. The repaired time-step budget is described
below. The probes have new task IDs, preserve the original tasks, and omit the
unconstrained endpoint optimization:

| Source branch | One additional calculation | Reason |
|---|---|---|
| rxn_12_2 forward_10 and reverse_10 | DVV, VeryTight | Test the one-point small-gradient stops with an independent integrator and tighter endpoint criteria |
| paper_ts2 forward_5 | DVV, VeryTight | EulerPC corrector oscillates above tolerance even at the smaller step |
| paper_ts1 forward_5 | EulerPC, MaxCycle=150, unchanged Tight | Recent correction sizes contract to 0.000012 against 0.000010; extension requires that trend in the full native log |
| paper_ts5 reverse_10 | DVV, VeryTight | Retain a raw independent path through the fragment exit, without subsequent recapture |

This is a pilot, not an assertion that DVV is converged or superior. A successful
pilot must be compared against the original path, then both directions and
step-size dependence evaluated before broader adoption. Plateaued correctors
never qualify for the iteration extension. Each recipe/source pair is fixed
once; a failed pilot ends for review. No convergence criterion is relaxed.
The documented [Gaussian IRC algorithms and criteria](https://www.conflex.co.jp/gaussian_support/irc.php)
support GradientOnly with DVV and distinguish Tight/VeryTight from iteration limits.

**TS2/TS3 coupled-cluster audit:** the native tables really contain respectively
five imaginary modes (428.54, 408.26, 214.84, 152.08, 5.16 cm^-1) and two
(191.15, 38.59 cm^-1). This is not a sign/parser error. Each rejected saved
geometry gets three fixed-geometry diagnostics: central numerical gradient at
0.005 bohr and numerical Hessians at 0.010/0.005 bohr, with SCF 1e-12 and CC
energy/squared-amplitude thresholds 1e-11. The ROHF-UCCSD(T)-F12b(T*) method,
frozen core, basis and molecular state remain unchanged. These six jobs cannot
export an accepted geometry or unblock the failed refinement's descendants.
See Molpro's [gradient](https://www.molpro.net/manual/doku.php?id=energy_gradients)
and [frequency](https://www.molpro.net/manual/doku.php?id=harmonic_vibrational_frequencies_frequencies)
definitions. Persistent negative curvatures alone do not establish a stationary
higher-order saddle: the gradient and electronic continuity must also be checked.

**W1 and native rxn_12_1 remain unresolved.** Their failed CC calculations print
T1/D1 of approximately 0.046/0.163 and 0.080/0.343. These values from unsuccessful
iterations are warning evidence, not certified diagnostics of a converged state.
Do not use their energies as benchmark barriers, suppress convergence checks,
or rerun them with larger CC iteration ceilings. Reference/orbital continuity
and the suitability of a single-reference treatment need review before selecting
another method or seed. The valid existing TS5 CC calculation and its dependent
work are preserved. Existing T/Q benchmarks remain available for comparison.

The 44 retained-L2 checks, 24 L1 exclusion checks, other panel work and ordinary
R8 jobs are preserved. Newly submitted original independent IRC grid tasks are
held for review (24 were staged in the supplied snapshot); existing allocations
are polled and adopted. The five probes and six curvature jobs share the existing
admission limits. No result is promoted to MESS or a validated PES edge.

## Observed probe failures and bounded second wave

The Blodgett owner did resume after the checkpoint-inspection fix. The user's
native output identifies three causes; launcher segmentation/MPI abort messages
alone had obscured them. These failed probes supply no new connectivity or
curvature result and cannot justify a disagreement with the paper.

* Jobs 53245741–53245744: Gaussian 16 B.01 read the accepted Hessian but rejected
  `StepSize=-5/-10` as negative **0.05/0.10 fs** DVV time steps, before the first
  propagated point. This was our input-generation defect. New task IDs ending
  `_positive_dt_v2` use `StepSize=5/10` (positive time), retaining the direction,
  VeryTight criterion and other controls. Maximum steps remain 1200/600; these
  are finite integration budgets, not equal spatial path lengths. The helper
  requires the explicit native invalid-step message and absence of propagated
  points, and re-verifies the accepted parent checkpoint. It cannot retry a
  subsequently failed positive-time task. Fresh DVV generation is also fixed.
* Jobs 53245756–53245758 (TS2): the tighter initial CC solve reached iteration
  39 before any derivative. The printed summed amplitude variances decrease
  from 3.1e-13 at iteration 30 to 1.87e-14 at 39. Printed energy changes round
  to zero but do not demonstrate the 1e-11 energy criterion. One retry uses the
  shared CC MAXIT procedure with a 100-iteration cap, gated on a small,
  contracting final native history in **each** failed audit output. This is a
  bounded test of continuing convergence; it does not declare it converged.
* Jobs 53245759–53245761 (TS3): default ROHF oscillated through 60 iterations,
  with density changes around 3.9e-5, and CC correctly refused the orbitals.
  One trial applies `{rhf,maxit=60;shift,-1.0,-0.5}`. It changes the convergence
  algorithm without increasing its limit, importing unverified orbitals or
  changing charge, spin, basis, geometry or method. This is Molpro's documented
  level-shift remedy for difficult conventional ROHF. The existing SO/SO-SCI
  reconnaissance strategies were reviewed; this finite trial uses the
  version-independent conventional solver. Root identity still requires
  comparison of the native orbitals and energies, even if it converges.

See the Molpro [SCF convergence guidance](https://www.molpro.net/manual/doku.php?id=the_scf_program)
and [CC convergence definitions](https://www.molpro.net/manual/doku.php?id=the_closed_shell_ccsd_program).
Neither `NOCHECK`, `IGNORE_ERROR` nor relaxed thresholds is used.

Each new curvature graph first runs an `electronic` calculation at the same
rejected saved geometry and exports a Molden file. Only successful electronic
execution releases the three fixed-geometry derivative tasks. Every solve
continues to enforce convergence. Each new task is attempted once; no automatic
cycle escalation is enabled. A successful gate is execution evidence, not
proof of orbital continuity or scientific acceptance.

`reconciliation_v1.json` remains byte-for-byte unchanged. A separate checksum-
bound `reconciliation_repairs_v2.json` records all ten failed source jobs and
their native evidence. The owner appends the four repaired DVV tasks through
the archival journal and stages two new electronic/derivative graphs. Old
attempts and the original scientific holds stay intact. The live TS1 pc150,
TS5 forward/reverse and screening allocations are polled/adopted unchanged.

## Update, resume, and inspect

Use the published commit specified in the chat, then run from `~/KinBot`:

```bash
bash examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/resume_probe_repairs.sh
```

The script prepares the hash-bound repair decision from the failed native
attempts and refreshes only the verified
panel controller. It neither cancels QC jobs nor resets old attempts. The owning
driver applies append-only graph changes through the existing archival journal.
A different source generation or ambiguous/live source allocation blocks the
decision instead of guessing. Repeating the command reuses the same decision.

Read-only, once (add `--interval 60` for a watcher; Ctrl+C stops the watcher):

```bash
.venv/bin/python -B examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/watch.py \
  "$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114/native"
```

Detailed raw-path/curvature review, also read-only:

```bash
.venv/bin/python -B examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/reconcile.py review \
  "$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114/native"
```

Do not use the old `campaign.py run` launcher for this wave: it predates the
decision's hold on unsent IRC grid tasks. Use `resume_probe_repairs.sh` for this
second wave. `resume_reconciliation.sh` remains the first-wave entry point.

## Validation boundary

The second-wave regression tests reproduce the negative-time rejection through
the real worker, then execute/adopt the repaired positive-time task with a
simulated Gaussian. They also test failed electronic gates blocking all
derivatives, successful gate execution releasing derivatives through the owner,
unchanged source files, repeated decisions/refreshes, source drift, held work,
and submission/staging/archive/scheduler races. The updated watcher prints
compact native causes and includes both decisions; collection includes the
new graphs and Molden evidence. The recorded test count is in
`reconciliation_validation.json`. Native convergence of the second wave remains
unverified until Blodgett executes it.

The first Blodgett handoff of `52e2069f3` stopped during decision preparation
with `Checkpoint import filename already exists`. It had enabled panel live
controls but had not refreshed the owner or submitted this diagnostic wave.
Source attempts correctly retained their imported checkpoints; review
incorrectly called the fresh-worker staging check. Dependency inspection now
verifies the accepted parent without treating retained source files as new
staging destinations. Existing attempt files stay in place, and the worker's
collision check remains strict. The updated tests first execute failed and
successful source workers, then exercise decision preparation, owner adoption,
probe execution, repeated starts and refresh races with those files retained.

Local tests exercise actual staging, dependency-checkpoint imports, archival
transactions, source/job verification, native-worker receipts, owner decision
adoption, live-allocation preservation, sampling holds/adoption, repeated starts,
interrupted publication/preparation and transient refresh races. Scheduler and
native executables are simulated. Publication is a tested software handoff,
**not verified Blodgett execution or a successful native convergence remedy**.
Native pilot results and complete path/mode/energy reconciliation remain to do.
The targeted suite passed 193 tests both in the working tree and in a clean
export of base commit `52e2069f3` with only this correction overlaid. The shell
handoff also passed `bash -n`. See `reconciliation_validation.json` for the
tested sources and command.
