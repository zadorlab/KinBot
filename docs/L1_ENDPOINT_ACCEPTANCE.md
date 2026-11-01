# L1 endpoint frequency acceptance

The endpoint frequency check uses the completed frequency calculation's own
electronic energy, ZPE and signed frequencies together. It verifies the recorded
L1 theory generation, selected atoms and geometry, normal QC completion, finite
energy/ZPE, and the existing minimum-frequency policy. It does not require an
independent endpoint force evaluation to reproduce the frequency calculation's
energy within two microhartree.

That former cross-calculation tolerance was not an SCF convergence or electronic
state test. In particular, a separate native Hessian can use different symmetry
bookkeeping. A different symmetry keyword alone does not invalidate it. The
energy difference is now recorded as a non-blocking diagnostic in
`l1_endpoint_checks.json`, alongside the source and selected frequency row IDs.
The selected energy and ZPE always come from the same frequency observation;
the old endpoint energy is not spliced into the new thermochemical result.

This does not accept missing spectra, invalid stationary-point order, changed
geometry, nonfinite energies, failed QC or unverified L1 theory generations.
It does not change SCF convergence settings or choose another electronic-state
model. Normal termination alone is not sufficient for accepting a reaction.

## Existing runs and monitoring

Refresh the owning controllers with the updated code using `kinbot.run_refresh`.
No reaction reset, database edit, repeated endpoint QC, or threshold override is
needed for a completed frequency observation blocked only by the old energy
comparison. The owner reads the existing result, preserves the old unresolved
entry in journal history, and continues the configured L1 depth/energy screening.
Completed/failed native jobs remain protected by the existing idempotent
`qc_freq` adoption logic. Unrelated accepted and rejected searches remain intact.

The read-only run monitor now prints pending L1 endpoint checks with their
reasons. An unresolved scientific check is labelled as blocked, not as a failed
reaction or an assertion that another calculation would help. Endpoint journals
belonging to already accepted/rejected searches are not displayed as live holds.

## Native regression evidence

`tests/reference/l1_endpoint_frequency_offsets.json` contains three exported
endpoint/frequency pairs, native input text, selected output lines, database
observations and input/output hashes. Its source is the read-only export
`r8-qc-decisions.IXLLO8dm`, archive SHA256
`4080c334d5c5de1b929d7345ed4311e83d7e0d8fd03238e1c6362fce3f95db0a`.

Both members of each pair terminated normally. Their input geometries, method,
basis, charge, multiplicity, SCF settings, integration grid and dispersion match;
the source uses `Symm(None)` and the Hessian uses `Symmetry(PG=C1)`. Printed spin
diagnostics agree. Native energies, frequencies and ZPE were checked against
the exported database. The respective energy differences are 2.365, 2.529 and
2.079 microhartree (0.0014841, 0.0015870 and 0.0013046 kcal/mol).

Local tests replay these observations without QC, retain database bytes/history,
and exercise repeated controller refresh plus owner adoption. The reaction-loop
test covers restart reconstruction and both sides of the 0.5 kcal/mol depth
threshold, including a nonzero source/frequency energy difference. Native QC
continuation on Blodgett remains distinct from local simulated-site validation.

## Related consistency checks audited

The following checks serve different purposes and were retained:

- `calculation.publish_optimization_result` and the MESS conformer source check
  compare an accepted observation to its own stored row, preventing accidental
  replacement by a later result.
- `qc_optimizer_recovery`, `anl.optimization_recovery`,
  `anl.legacy_optimization`, and `anl.optimization_progress` compare an optimizer
  trajectory step to that same step's printed log before restoring it. Their
  tolerances accommodate log rounding, not independent SCF evaluations.
- Native result readers and VRC energy reuse compare stored energies with the
  specific native output they claim as provenance.
- Theory-route fingerprints bind work to its generation. They are not a test
  of equality between a force job's symmetry keyword and its separate Hessian's.
- The hindered-rotor zero-angle check concerns consistency of a scan with its
  reference potential, uses 0.1 kcal/mol, and can fall back to harmonic modes.
  These three endpoint observations do not establish a defect in that check.

No other demonstrable false rejection of these normally completed endpoint
results was found in this audit. This is not a claim that every QC failure or
every possible acceptance criterion across all backends has been validated.

## Accepted-reaction restart regression

The subsequent native restart exposed a separate reconstruction bug:
`_endpoint_l1_check` skipped accepted searches. That skipped observation reuse,
not just re-screening. The old raw IRC row often has no ZPE, so replay recomputed
the complex depth using zero endpoint ZPE and tried to optimize an invented
L2 complex. The no-new-QC replay guard correctly stopped that request, but the
owning worker and PES then exited before reaching unfinished searches.

Accepted replay now reads the same validated endpoint and fragment observations
as the original search. `qc_freq` adopts an existing frequency job; the guard
against missing/replacement QC is retained. Previously rejected searches and
accepted absolute-energy screening decisions are unchanged. No reset or new
endpoint/high-level job is required to repair these six cases.

`tests/reference/l1_accepted_replay.json` records the six reported examples,
source row IDs, endpoint/fragment observations, native frequency-output hashes
and relevant native lines from the same exported archive. Native energy, ZPE
and frequencies were checked against the database. Matched depths are
0.169140, -0.115794, 0.181514, 0.050335, -0.002648 and 0.221489 kcal/mol;
the erroneous zero-endpoint-ZPE depths were 13.013–17.695 kcal/mol. All six
therefore remain below the configured 0.5 kcal/mol complex threshold.

The old test expecting accepted replay to skip this check was removed. New
regressions use all six exported cases and run the real reaction loop through
accepted reconstruction, complex classification, L2 dispatch and final complex
assessment. Simulated cached QC exercises the actual fixed-frequency adoption
and submission guard: shallow complexes deliberately have no `_prod_high`
result, and retained complexes have one. Repeated interruption/restart preserves
the terminal search record and leaves database and screening records unchanged.
Refresh tests cover both active and accepted endpoint adoption after relaunch.
These are local replay/simulated-site checks, not a claim that the new remote
controllers have already finished.

## MESS artifact handoff after completed exploration

The subsequent monitor showed all three CH3S and thirteen CH3SS workers complete.
CH3SS failed only when assembling the PES-wide MESS input: the selected
`12_shift_S_R_3_2_4_IRC_F_prod_0000.mess` was absent. In the same exported evidence
described above, worker `791871560500020000002-sb6e34ded5baaf03b` has accepted
forward/reverse TS observations 759 and 762 with identical B2PLYP total energies
of -835.88114038804 Hartree and ZPE 0.028943 Hartree. Their complex observations
745 and 754 also agree, at -835.9628768711201 and 0.028476 Hartree respectively.
All four associated Gaussian outputs terminate normally; the B2PLYP totals
(including the perturbative contribution) and ZPE agree with the database.

The worker retained the first equal-energy route's product artifacts, while
PES assembly retained the later complex route. Workers now export all accepted
route artifacts in PES mode, including their optical validation. Final PES
selection and standalone lowest-barrier selection are unchanged. Product exit
artifacts additionally include the reaction in their bounded-length filename,
so two routes to the same fragments cannot overwrite each other's deferred
PST/ROTD context. Legacy shared product files remain readable with the existing
context validation. No geometry, energy, frequency, search decision, or QC
attempt is changed by artifact generation.

The regressions call the actual worker writer, summary writer, database energy
lookup, PES selection, MESS assembly/finalization, and network splitter. They
cover equal/unequal TS energies, reversed order, distinct paths, single/multiple
conformer modes, PST-only and automatic modes, retained/omitted complexes, and
repeated generation. A simulated controller refresh preserves native files and
live allocations, then the owning writer regenerates the missing model and
completes assembly. These are local tests; successful remote assembly still
requires the normal pull/refresh handoff and its completion record.

## Disconnected-complex identity in MESS assembly

After the artifact fix, the 2026-11-01 monitor reaches mirror completion but
reports `A configured MESS population has two different names`. Two distinct
disconnected complexes share legacy ChemID `791230350000000000002`:

- `791671450790040000002_r12_insertion_R_2_1_5_IRC_R_prod_high`
  (row 793; L1 row 621) represents CH2···HSS.
- `791871560500020000002-sb6e34ded5baaf03b_12_shift_S_R_3_2_4_IRC_F_prod_high`
  (row 754; L1 row 645) represents HCS···H2S.

Their canonical molecular graph identities differ. Using only the legacy
numerical routing name in MESS population annotations aliased the two wells.
Merging them or bypassing the duplicate-name check would join different
fragment partitions incorrectly. Disconnected population annotations now also
include the existing full graph/charge/multiplicity/configuration identity.
The same identity is used on barrier endpoints and reflected populations.
Connected species names, QC names, optical factors, energies, geometries,
frequencies and accepted search decisions are unchanged. The duplicate guard
remains and now identifies the conflicting key and names.

`tests/reference/mess_complex_identity_collision.json` preserves the source
archive hash, native L1/L2 observations, L2 input text and native output hashes
for both parents, saddles, complexes and fragment pairs (ten points total).
The outputs terminate normally. These are observations from the previously
exported archive, not new quantum calculations. With the old identity code,
the real worker writer, summary, database lookup and PES assembly reproduce
the reported exception. With the fix, both complexes and their respective PST
exits survive as separate models. The regression checks both worker orders,
specified/racemic populations, single/multiple conformer modes, omission of
complexes, final network splitting and repeated generation. Structural tests
also check atom-order and orientation independence and mirror/state identity.

The refresh regression starts from artifacts with the ambiguous old metadata,
preserves native files and live sampling allocations, then runs the owning
writers and final assembly twice. The repaired input files and scientific data
remain stable. This is native-backed local replay with a simulated scheduler;
remote completion must still be confirmed after pull and controller refresh.
