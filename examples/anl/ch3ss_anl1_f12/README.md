# CH3SS orbital reconnaissance and CBH/ANL–VRC–MESS

The reconnaissance treats the two dissociations separately:

| Coordinate | Products | Coupled multiplicities | Spin-free spatial roots per multiplicity |
| --- | --- | --- | --- |
| C–S | CH3 + S2(X 3Sigma) | 2, 4 | 1 |
| S–S | CH3S(X 2E) + S(3P) | 2, 4 | 6 |

Root counts follow declared fragment terms. The second count retains both
CH3S spatial components and all three S(3P) components. A model retaining a
smaller vibronic manifold requires separate evidence and explicit input;
spin multiplicity alone cannot determine that choice.
These are probe/state-average counts, not assumed reactive-flux counts.
Inner and outer contributing states require separate evidence and statistical
counting. KinBot orchestration uses independent existing ROTD_py scalar
flux calculations; it does not assume every averaged root is reactive.

`orbital_reconnaissance.json` covers aug-cc-pVDZ/TZ/QZ/5Z and their F12 variants.
Molpro supplies native basis and auxiliary definitions, with no forced
revision suffix. The two grids initially contain 75 geometries (39 radial and 36
orientation probes) and 600 basis/geometry jobs, plus 16 native pilot jobs. Each job uses one node. The controller limit is ten
jobs; diagnostics can request denser intervals or longer separations.
Neither 8 Å nor the initial 30 Å is accepted as an asymptote without a
plateau check.

Probes save ROHF and AVAS/CASCI natural orbital evidence, every requested
root energy, Molden files, and WFU records. Subsequent points stage a verified
predecessor and use START/DIAB. Refinement preserves completed checkpoint
histories. Occupation-based candidate indices and electron populations are
reported; orbital character and wavefunction root continuity need review.

The production input requests a **candidate** CAS(9,7), intended to retain
C–S/S–S bonding/antibonding and sulfur nonbonding/radical orbitals. Both
channels must support the selected subspace, or the common space must be
enlarged. `ch3ss_state_model.json` is absent until reviewed, basis-specific
WFU records exist. Missing records or changed hashes prevent MR submission.

Atomic S levels and an S2 triplet identity are supplied as example data.
S2 fine structure, CH3S spin–vibronic partition functions, molecular SO
corrections, and the long-range statistical state model remain chemistry checks.
See the [full audit](../../../docs/ANL_workflow_audit.md) before production.

Primary references: [CH3S electronic states](https://pmc.ncbi.nlm.nih.gov/articles/PMC5995305/),
[Molpro MULTI](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi).

Use the [Blodgett commit/pull and launch runbook](../../../docs/Blodgett_workflow_runbook.md)
to prepare both systems in one campaign. Preparation records available Molpro
modules and each batch job loads its setup. Preflight and all-basis/channel
native pilots precede bulk submission. The shared limit is ten nodes total.
Full KinBot production follows review of these results. Automatic homolytic
vdW-complex discovery is deferred under the present example scope.

The [current paused checkpoint](../../../docs/WORKFLOW_CHECKPOINT_2026-10-14.md)
records native pilot feedback, scratch recovery, verified behavior and all
remaining software/chemistry work. Full production is not yet released.
