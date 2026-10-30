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
