# CH3S orbital reconnaissance and CBH/ANL–VRC–MESS

The repaired CASSCF reference panel is complete. Continue with the
[finite correlated validation batch](../../../docs/SULFUR_CORRELATED_VALIDATION.md):
verified natural-orbital partitions, CASPT2-F12 comparisons, a shared 20-job
cap, and explicit steps through production acceptance. The dense scan and
ten-job settings described below are historical, not the current launch.

Use the [current chemical orbital-review workflow](../../../docs/ORBITAL_REVIEW_WORKFLOW.md)
for fresh AVAS/CASSCF comparisons, verified checkpoint reuse, explicit
inner/outer selections and the contained Blodgett restart/export commands.
The legacy dense CASCI scan described below remains evidence for comparison;
automatic refinement is now opt-in.

The requested dissociation is CH3S → CH3 + S. The reconnaissance uses
aug-cc-pVTZ-F12 (`avtz-f12`), with Molpro's installed
basis definition and auxiliary defaults. KinBot supplies no revision suffix
or replacement basis data.

`orbital_reconnaissance.json` defines the C–S bond by atom indices. The driver
builds a dense near-bond grid, includes the reference geometry, extends to an
initial 30 Å, and requests extra intervals/range when energies or occupations
have not settled. There are initially 31 geometries (19 radial and 12 orientation probes)
and 31 basis/geometry jobs, plus one native pilot job. Adaptive refinement can add jobs. Each job uses one node; the supplied
controller limit is ten jobs. See the [audit](../../../docs/ANL_workflow_audit.md)
for the release gate and remaining chemistry checks.

Probes produce ROHF orbitals, AVAS/CASCI natural occupations, all requested
root energies, native XML exports (legacy Molden remains readable), and WFU records. Subsequent points use the
preceding saved record with START/DIAB. Hashes bind results to their inputs
and checkpoint sources. Reports list occupation-based candidate orbital
indices and electron populations at each geometry. These diagnostics do not
establish a minimal active space or certify wavefunction root continuity.

The parent doublet and declared CH3(doublet) + S(3P) products determine the
doublet/quartet coupling and three spatial roots per manifold. This follows
input evidence, without a molecule-name branch in the implementation.
These are probe/state-average counts. The number of reactive contributions
at inner and outer bottlenecks must be established separately; see the audit's
NH2O SI comparison. General multistate orchestration is implemented;
the contributing sulfur states must be selected from the evidence.

`ch3s.json` requests a **candidate** CAS(13,8). A count does not identify the
orbitals. Production requires a reviewed, basis-specific
`ch3s_state_model.json` and verified WFU records. Review must include the C–S
bonding/antibonding and sulfur valence manifold throughout dissociation,
larger-space sensitivity, and the CH3S Jahn–Teller/spin–orbit treatment.
Update the requested CAS if the evidence calls for a different space.

The [completed CH3S evidence review](../../../docs/CH3S_RECON_REVIEW_2026-10-15.md)
recommends testing the chemical CAS(5,4) baseline against CAS(9,6) and CAS(13,8).
The frozen-orbital probes do not establish the best CAS, and propagated/fresh
CASCI orbitals give different
asymptotic energies. Optimized-orbital validation is required before selecting
a production record; the candidate input remains CAS(13,8).

Once CH3S probes finish, its review can proceed while CH3SS runs. Use the
[CH3S-only evidence export](../../../docs/SULFUR_TZ_RESTART_2026-10-14.md#review-completed-ch3s-while-ch3ss-continues)
instead of the campaign-wide collector, which requires the controller lock.
Remaining coverage flags must be examined before approving the selected space.

The input supplies atomic S fine-structure levels and a ground-term centroid
correction. Its CH3S X(2E) declaration specifies spin-free spatial roots only,
not a resolved molecular partition function or molecular SO correction.

Primary references: [CH3S multistate and SO study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5995305/),
[Molpro AVAS](https://www.molpro.net/manual/doku.php?id=automated_construction_of_atomic_valence_active_spaces),
[Molpro MULTI](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi).

Use the [Blodgett commit/pull and launch runbook](../../../docs/Blodgett_workflow_runbook.md)
to prepare both systems in one campaign. Preparation records available Molpro
modules and each batch job loads its setup. Preflight and all-basis/channel
native pilots precede bulk submission. The shared limit is ten nodes total.
Full KinBot production follows review of these results. Automatic homolytic
vdW-complex discovery is deferred under the present example scope.

The [current implementation and acceptance procedure](../../../docs/WORKFLOW_READINESS_2026-10-14.md)
records recovery, verified behavior, and the chemistry review/native checks
required after reconnaissance. The earlier paused checkpoint is historical.

Molpro's Molden interface rejects h and higher angular momenta. Recon execution
now derives a hash-bound XML-export deck from the immutable prepared input;
WFU/native orbital evidence retains the complete requested basis. Existing
Molden results remain valid. Export-only failures with an evidenced converged
RHF record can recover that record through the checked retry path.
The production preset uses augmented TZ-F12 for the relaxed path, sampling,
CC high-spin anchor, both MR gaps, and direct MR leaves. There is no routine
VRC CBS extrapolation; explicit expert overrides remain supported. The CBH/ANL
thermochemistry recipe is independent and unchanged.

Finish the current TZ recon before choosing alternative production bases.
The [post-recon performance check](../../../docs/SULFUR_TZ_RESTART_2026-10-14.md#performance-check-after-reconnaissance)
will test actual correlated costs; DZ sampling or additional QZ recon is
conditional on those results. The present TZ settings remain the baseline.

See the [current safe Blodgett stop/update/restart commands and active-space guide](../../../docs/SULFUR_TZ_RESTART_2026-10-14.md).
