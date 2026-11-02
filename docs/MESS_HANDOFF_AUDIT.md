# Complete L1/L2/PST MESS handoff replay

The 2026-11-01 capture contains the completed CH3S and CH3SS scout workers,
their databases, accepted/rejected search histories, worker MESS fragments,
inputs and PES handoff records. Archive SHA-256:
`9a65dbb987f4bfaf81eb96ec25c0163b7fafd3338da52a7ed76673a6ccd3dc0a`.
All 9,177 captured files passed their embedded hashes. The original capture is
retained; replay runs use separate local copies.

## General fixes

1. **Mirror completion and homolysis fallbacks.** A reflected, classified
   saddle can connect the same populations as an existing `hom_sci` PST
   fallback. Final assembly now reapplies the existing saddle-over-fallback
   rule after reflection. Only explicitly identified pure-PST homolysis
   placeholders qualify; unrelated parallel paths, ROTD models and distinct
   classified saddles are preserved. A comment identifies each superseded
   placeholder. Older worker fragments are annotated in memory, without QC
   resubmission or rewriting the worker files.
2. **Native geometry guard.** MESS's default minimum interatomic separation
   rejects the valid short H2 geometry in the captured CH3S input. The common
   MESS header now explicitly uses `AtomDistanceMin[angstrom] 0.32`, matching
   KinBot's existing partition-function headers. No coordinates or QC
   acceptance thresholds change.
3. **Backend independence.** With ROTD disabled, postprocessing no longer
   creates a ROTD setup or inspects its scheduler. With L3 disabled, it no
   longer reads L3 energies, copies L3 directories, writes L3 submission files
   or performs Molpro T1 analysis. `no-kinbot` also avoids querying a process
   status for its synthetic worker IDs. The ordinary owning driver still
   renders the complete PES and MESS inputs.
4. **Explicit scouting approximation.** Both example scout inputs enable
   the existing `correct_submerged=1` option. When requested, final assembly
   records every below-well barrier correction in the MESS input. The captured
   CH3SS model has one such adjustment: `rxn_37`, 41.67 to 42.34 kcal/mol.
   Eckart tunneling is removed at the floor. This is a statistical-model
   approximation, not a corrected QC result; original energies remain intact.
   CH3S requires no such adjustment. The default for other runs is unchanged.

The earlier population-identity fix also remains active: disconnected vdW
complexes sharing a formula/legacy identifier retain distinct MESS populations
when their partner graphs differ. Bimolecular exits do not connect otherwise
disconnected well networks.

## Validation and limits

The complete captured pair passed the real `kinbot.pes ... no-kinbot` path
twice. Final counts are **3 wells, 6 product sets, 10 barriers** for CH3S and
**17 wells, 18 product sets, 91 barriers** for CH3SS. The latter includes four
explicit complexes in addition to the thirteen explored wells.

Protected original inputs, databases, search states, screening decisions and
worker MESS fragments are byte-identical before/after replay. Tests prohibit
QC submission, worker launches, process-status polling and disabled L3/ROTD
backends. Repeated output is identical. The optional postprocessing helper is
also tested for live-work/incomplete-worker refusal before writes and repeated
preparation. Simulated-site refresh tests exercise transient staging,
submission recovery, live job preservation, sampling adoption and idempotence.

The isolated release validation finished with **405 tests passed, one
platform-specific test skipped, and 69 subtests passed**. Native MESS is
included in this run; it is not a mock of the solver.

Native MESS source archive SHA-256:
`911085be5322c7de2686ac789d0f3a2516f70b971c19ee88f55d542f6a45c908`.
Local executable SHA-256:
`735853406902b8d3807a49e4d1e435417a29355043ba40078147a304631ae47b`.
The local build uses the archived MESS source with its numerical dependencies;
the build provenance records a C++ include portability adjustment, not a
kinetic-model change.

CH3S completes the full requested 300–2000 K, 0.01–760 torr grid. The complete
CH3SS model completes a 1000 K, 760 torr native calculation and writes the full
0–200 kcal/mol microcanonical grid. Full CH3SS runs at the original energy
grain exceeded the local 900-second diagnostic limit; that is not a model
validation failure or a reason to rerun QC.

Some double-precision pressure-dependent tables contain materially negative
off-diagonal rates. A normal MESS exit does **not** certify those tables.
They must not be used for pruning without a numerical convergence/reduction
check. Native microcanonical state counts/rates and high-pressure capture
coefficients are separate outputs; collisionless branching derived from them
must be labeled with its assumptions, not presented as pressure-dependent
photolysis yields. A 193 nm photon corresponds to 148.142152 kcal/mol, not to
a unique temperature. The default PST coefficient is a scouting assumption,
not a quantitatively fitted capture potential.

## Reproduce locally

The large source capture stays outside Git. Compact native excerpts with
hashes are ordinary regression fixtures in `tests/reference/`. With the full
capture and a native executable available, run:

```bash
KINBOT_TEST_MESS_HANDOFF=/absolute/path/seed-evidence.tar.gz \
KINBOT_TEST_MESS=/absolute/path/mess_driver \
  .venv/bin/python -m pytest -q \
  tests/test_mess_handoff_replay.py tests/test_mess_mirror_capture.py \
  tests/test_mess_native_scout.py tests/test_run_refresh.py
```

Local replay and native execution do not require a further Blodgett refresh.
For an optional remote rebuild only, the archived, input-only helper is
`examples/anl/ch3s_anl1_f12/scouting_r8/postprocess_blodgett.sh`. It verifies
both cases are idle/complete and uses derived `mess_scout.json` files; it does
not change original run inputs, submit QC, or run MESS.
