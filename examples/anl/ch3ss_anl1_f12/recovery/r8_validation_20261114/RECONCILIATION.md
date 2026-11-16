# Native reconciliation update — 16 November 2026

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
charge, spin, B2PLYP-D3(BJ)/cc-pVTZ, grid, SCF settings, step size and path-length
budget. They have new task IDs, preserve the original tasks, and omit the
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

## Update, resume, and inspect

Use the published commit specified in the chat, then run from `~/KinBot`:

```bash
bash examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/resume_reconciliation.sh
```

The script enables controls only at the native panel, prepares a hash-bound
decision from the current native attempts, and refreshes only the verified
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
decision's hold on unsent IRC grid tasks. Use `resume_reconciliation.sh`.

## Validation boundary

Local tests exercise actual staging, dependency-checkpoint imports, archival
transactions, source/job verification, native-worker receipts, owner decision
adoption, live-allocation preservation, sampling holds/adoption, repeated starts,
interrupted publication/preparation and transient refresh races. Scheduler and
native executables are simulated. Publication is a tested software handoff,
**not verified Blodgett execution or a successful native convergence remedy**.
Native pilot results and complete path/mode/energy reconciliation remain to do.
The targeted suite passed 187 tests both in the working tree and in a clean
export of base commit `fa8bed91c` with only this update overlaid. The shell
handoff also passed `bash -n`. See `reconciliation_validation.json` for the
tested sources and command.
