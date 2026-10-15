# Sulfur workflow implementation and acceptance — 14 October 2026

The user resumed implementation after the historical pause. This document
supersedes that checkpoint's open-software list. The existing reconnaissance
can recover native exports using published commit `c14913f`. The cost-reduced
restart additionally requires the newer `exclude-bases` command. Original
inputs, completed results, excluded-basis evidence, and the directory remain intact. The resumed MESS audit
is described below. No SSH or cluster calculation has been performed here.

## Theory and cost policy

The examples request **ANL1-F12** for CBH/ANL, with augmented F12 bases.
The default VRC production preset now uses:

| Calculation | Level |
| --- | --- |
| Relaxed correction path | Selected-root CASPT2-F12b / aug-cc-pVTZ-F12, numerical gradient |
| ROTD sampling | Selected-root CASPT2-F12b / aug-cc-pVTZ-F12; one basis at each sample |
| High-spin anchor for a lower-spin target | ROHF–UCCSD(T)-F12b / aug-cc-pV{T,Q}Z-F12 |
| First spin-gap estimator | Selected-root CASPT2-F12b / aug-cc-pV{T,Q}Z-F12 |
| Second spin-gap estimator | Selected-root MRCI-F12b+Q / aug-cc-pV{T,Q}Z-F12 |
| Combination | CC high-spin energy plus the arithmetic mean of the two target-minus-high-spin MR gaps |

Each correction branch retains its configured two-point CBS extrapolation.
The inherited default exponent is 3.7 and is configurable; its use is a model
choice requiring basis sensitivity checks, not a demonstrated F12 error bound.
MRCI uses the relaxed-reference Davidson correction and the F12b energy.
An already-high-spin physical leaf uses its explicitly declared direct trusted
level rather than reversing the spin-gap construction. An excited physical
high-spin root needs a root-resolved MR potential; ground-state CC does not
automatically calculate arbitrary excited CAS roots.

Both MR methods use the declared **state-averaged CASSCF reference** and
correlate a selected state/root. "Single-state" here refers to correlation,
not separately optimized state-specific orbitals. Reviewed basis-specific WFU
records initialize the chosen subspace; CASSCF subsequently optimizes it.
The reconnaissance instead uses cheaper ROHF/AVAS/CASCI diagnostics. Counts
CAS(13,8) and CAS(9,7) remain provisional requests, not certified orbital choices.

Correction reports retain each estimator, its CBS components, and method
spreads in hartree, including the spread after each curve's asymptotic zero is
removed. These spreads diagnose method dependence; they are not statistical
uncertainties and averaging is not proof of improved accuracy.

There is no 5Z production default. The current recon examples deliberately
retain conventional and F12 D/T/Q/5Z for the user's requested comparison.
Changing the production preset does not alter a snapshotted recon campaign.
ROHF/CASCI basis comparisons alone cannot establish the accuracy gain in a
correlated capture potential. Decide whether to use 5Z production checks from
representative correlated interaction energies/gaps after orbital review.
All F12 names use Molpro's installed **default** definitions: no rev2 suffix,
basis synthesis, or auxiliary substitution. Explicit expert sampling-CBS
requests and older scalar inputs remain supported.

Following the long-running last pilot, the recommended current restart uses
native conventional and F12 **D/T/QZ**, excluding both 5Z families in place.
This gives 186 CH3S + 450 CH3SS full probes and 6 + 12 pilots: **654 initial
jobs**, versus 872 with all eight bases. Adaptive additions remain possible.
The six completed CH3S pilots and all excluded-basis artifacts are retained.
The selection is stored in manifests and survives refinement; original input
snapshots are not edited. The original eight-basis example configurations
remain available for an explicitly chosen future comparison.

## Bugs corrected during the resumed audit

1. **High-basis export failure.** Molpro completed 5Z RHF but its Molden
   exporter rejected h functions. Recon now executes a separately hashed
   XML-export deck while preserving prepared inputs. Old Molden successes
   remain readable. A scheduler-checked retry archives the failed attempt
   and can recover its converged RHF record when exact input, native output,
   and WFU hashes establish that only the export failed. The same correction
   covers 5Z-F12. Production correction decks also export XML, after CASSCF,
   instead of invoking the same unsupported Molden path.
   The first CH3SS XML pilot exposed a second export fault: its long generated
   `***` title overflowed Molpro's XML string insertion. Execution now uses a
   short deterministic ASCII title in both recon and production correction
   decks. The original recon input/hash stays immutable; no chemistry changes.
   Recovery also recognizes this specific XML failure, requiring the same
   converged-RHF and exact-input/checkpoint evidence before reusing SCF.
2. **Deferred full-PES loss of ROTD results.** Workers now save hash-bound
   per-UQ convolution contexts containing their modes, rotors, stoichiometry,
   electronic counting, and symmetry. Final PES assembly revalidates the
   contexts/native results, stages scalar or multistate Ne into global `me`,
   and renders with final network/UQ energies. Selected required channels
   cannot silently become PST. Legacy required artifacts without context
   require regeneration of worker MESS artifacts, not blind QC resubmission.
3. **Duplicate PES sampling.** Workers prepare inputs; the global driver owns
   ROTD execution. Accepted automatic homolytic selections are restored at
   aggregation instead of being lost with worker-local dictionaries.
4. **Requirement applied to unrelated exits.** Explicit selections govern
   which exits require ROTD. Unselected optional vdW exits retain their
   existing PST/well policy. With no selectors, the historical global required
   policy remains available. Physical electronic counting is independent of
   whether an optional complex is retained.
5. **Forced reference root.** Scalar leaves no longer overwrite a reviewed
   reference root with 1. Stage-specific selected-root overrides are supported.
   A ground-state-only SCF/CC method explicitly rejects an excited CAS-root
   request instead of pretending to honor it. Root-resolved MR anchors remain
   configurable.
6. **Charged AVAS initialization.** VRC correction decks explicitly set the
   requested charge and spin before an AVAS/atomic-density guess.
7. **Misleading conformer convolution.** One selected-geometry flux cannot
   stand for multiple fragment conformers. Explicit ROTD requests reject an
   ensemble or a mismatched selected record; auto mode retains the existing
   PST conformer union. An identical single selected record remains supported.
   Correct general ensemble VRC needs separate potentials/fluxes per conformer
   combination; it is future scope and is not required by these examples.
8. **Fixed separated-fragment distance.** Production scan geometry and energy
   assembly use `vrc_tst_asymptote_distance`, with legacy default 30 angstrom.
   It must exceed the scan grid and should follow the recon plateau evidence.
   The value itself is not a convergence certificate.
9. **Stale provenance.** Direct/PES writes now record separate UQ provenance
   and replace stale records, including when ROTD is switched off or a
   permitted fallback is selected. Negative Ne thresholds and the legacy
   inserted zero-flux origin remain supported; threshold energy is included
   in explicit RRHO ZeroEnergy exactly once.
10. **Insufficient scratch floor.** Recon accepts an optional positive
    `scratch_min_mb`. `KINBOT_MOLPRO_SCRATCH_MIN_MB` can raise the free-capacity
    requirement across native Molpro paths, including existing campaigns,
    and cannot lower an explicit task requirement. The selected root, free
    capacity and requirement are recorded. This is a pre-launch free-space
    check, not a quota check, reservation, or exact integral-disk forecast.
11. **Restart cost control.** `orbital_campaign exclude-bases` checks every
    manifest for live allocations/uncertain submissions before changing any
    selection. It re-audits retained results, preserves excluded job metadata
    and files, and is idempotent. Review collection includes excluded artifacts
    and failed-attempt diagnostics; large WFU files stay on the cluster.
    Refinement retains the original site setup instead of rediscovering a
    potentially different module environment midway through the campaign.
12. **Scratch health and hidden launcher errors.** The CH3SS conventional-QZ
    launcher reported `mktemp: ... Input/output error` on `/scratch`, before
    creating native output. Free capacity alone cannot detect this failure.
    The shared Molpro runtime now probes nested-directory creation, write,
    fsync, rename, read-back and deletion before native launch, and handles
    allocation/probe failures when considering another root. An explicit
    `KINBOT_MOLPRO_SCRATCH` remains strict unless the user supplies an ordered,
    path-separated `KINBOT_MOLPRO_SCRATCH_FALLBACKS` list. Each candidate must
    satisfy the same capacity floor. Host, selected root and rejected
    candidates are recorded; reconnaissance also retains launcher stdout and
    stderr and includes their tails in errors. The shared Slurm I/O retry
    policy below supersedes the earlier manual-only recovery. A probe cannot
    prevent later filesystem failures or repair a mount.
    The follow-up native QZ-F12 record independently reports `Permission
    denied: /scratch/hram`. Usable scratch capacity on one node is not evidence
    of permissions or mount health on another. Nodes for jobs `52907311` and
    `52907317` require administrator lookup in the restricted controller log
    `/var/log/slurm-llnl/slurmctld.log`, then storage/mount and directory/ACL
    checks respectively. The underlying cause of EIO remains unresolved.
13. **Shared automatic Slurm I/O recovery.** All KinBot batch generators now
    use `kinbot.slurm_retry`: ordinary QC and its custom templates, independent
    CBH/ANL tasks, orbital recon, ROTD samples, standalone Molpro single points,
    MESS/UQ and MESMER. Two batch retries are enabled by default, configurable
    through `KINBOT_SLURM_IO_RETRIES` (0 disables; maximum 10). Only explicit
    failed OS/native I/O diagnostics or a classified scratch setup error
    qualify. Syntax, basis, convergence, multiplicity and generic nonzero exits
    do not qualify. A zero-exit ordinary QC worker with an error DB row and
    I/O diagnostics is also eligible.
    Retries stay inside the same allocation and retain its node/CPU/memory/time
    limits; they never submit duplicate jobs or bypass the ten-node cap.
    Failed job-owned results are archived and exact pre-job checkpoints are
    restored, rather than consuming a WFU/CHK from a failed write. Shared DBs,
    unrelated job files and Slurm-owned output streams are left in place.
    Molpro retries prefer another explicitly allowed scratch root and create
    a fresh isolated repository. Logs, seeds, host, allocation and attempt
    outcomes remain in `.kinbot_slurm_io_retries`.
    Ordinary QC publishes its final marker only after the guard finishes.
    ROTD workers publish completed pickles atomically and repeat the original
    random sequence on retry. This prevents partial completion consumption
    and silently replacing the attempted Monte Carlo samples. Existing staged
    recon/ANL scripts acquire the guard when next submitted, after their
    existing provenance checks; running jobs are not altered.
    Failures before Slurm starts the batch script, an unreadable working
    filesystem that prevents journaling, and unknown submission outcomes
    cannot be retried by this in-allocation guard. They remain explicit
    scheduler/site failures; the guard is not a storage repair mechanism.

All state routing follows input species/fragment terms and spin coupling;
the implementation contains no CH3S/CH3SS name-dependent physics. Bond indices,
physical term data and provisional active spaces reside in the example inputs.

## NH2O SI verification

The actual ten publisher TXT SI files were checksum verified. The verifier
was rerun after the resumed changes. It passes topology, explicit-level
accounting and actual-Ne-table convolution checks. The NH2+O outer components
give degeneracies 10, 6, 2 at 0, 158.5, 226.5 inverse cm. The isolated canonical
level-counting discrepancies are below 0.05%. See `NH2O_SI_verification.json`
and `tests/verify_nh2o_si.py` for hashes and reproducibility.

This verifies the counting/staging approach, not the paper's full rates. The
original names of every SI Ne table cannot be unambiguously mapped to the
publisher-numbered files. The SI retains entrance complexes. These must not
be imposed on the sulfur examples: automatic hom_sci complex characterization
is deferred by the user. Reduced statistical/minimum two-bottleneck models
are identified as reductions, not reproductions of that retained-well network.

## Exact next procedure

1. Let the remaining current pilot allocations finish (last supplied ID
   `52907312`). The D/T/Q selection is already stored; retain it. Pull/install
   on Blodgett and select `/scratch2/$USER/kinbot-molpro`, with shared home
   scratch as the explicitly permitted fallback. Run checked `retry-failed`,
   preflight, and resume the same controller/directory with ten one-node jobs
   maximum. The command block is at the end of the runbook. Running jobs are
   not cancelled or retroactively changed by the source update; the checked
   retry refuses live campaign allocations.
2. Complete all basis/channel native pilots, then adaptive radial/orientation
   recon. Preserve successes. The user reported five completed CH3SS pilots
   after the title fix, in addition to the six completed CH3S pilots; those
   are native acceptance evidence for their particular decks, not for all
   remaining jobs. Local fixtures cannot substitute for native execution.
   Unresolved coverage after the refinement cap remains a review requirement.
3. Collect the review archive; keep WFU records on Blodgett. Review all three
   dissociations separately: CH3S C–S, CH3SS C–S, CH3SS S–S. Inspect orbital
   character, natural populations, CI character, asymptotes, orientation and
   root changes, not just a population cutoff or a few distances.
4. Select basis/channel-consistent active spaces and WFU records, check larger
   spaces, and prepare reviewed state/surface models. Choose contributing
   physical roots, offsets, statistical levels, and bottlenecks. State-average
   roots are not automatically independent reactive channels. The production
   JSONs deliberately require these files before expensive execution.
5. Establish the sulfur energy/partition-function model: atomic S levels are
   supplied; molecular CH3S Jahn–Teller/SO/vibronic and S2 fine structure, and
   molecular SO corrections, still need physical evidence. Missing molecular
   corrections are labeled spin-free approximations, never asserted exact zero.
6. Validate one licensed native correction → ROTD → Ne → MESSPF/MESS handoff
   per channel. Check reference quality, CASPT2 intruders/shift sensitivity,
   MRCI+Q size consistency, separated-fragment gaps, basis/CBS/model spread,
   numerical gradients, pivot/surface and Monte Carlo/grid convergence.
7. Validate rotor/harmonic fallbacks and collision/energy-transfer inputs;
   the example Lennard-Jones values are inputs, not fitted sulfur constants.
   Then launch the reviewed full KinBot campaigns under the shared node cap.

The software audit does not replace steps 3–7. There is no local licensed
Molpro, Slurm, or MESS binary, and no claim of native full-rate validation.
No automated minimal-active-space certificate or coupled nonadiabatic dynamics
is provided. Full conformer-specific VRC and automatic hom_sci vdW discovery
remain explicitly deferred capabilities, with incorrect silent reuse blocked.

## Acetaldehyde and compatibility

The cluster restart killed the acetaldehyde campaign and explains its stale
running labels. HIR harmonic fallbacks and the optional vdW optimization
failure were nonfatal. The earlier first-force Molpro exit 255 still requires
its own native logs for a root cause; the sulfur export/scratch errors do not
prove the same cause. The runbook's run-scoped stop/archive preserves those
attempts and does not kill unrelated jobs or delete shared scratch.

Current regression results are recorded in the runbook/release handoff.
The final suite, including shared Slurm I/O retries and atomic ROTD completion,
passes **1295 tests and 391 subtests**, all 161
non-template modules parse, and the diff and documented shell syntax checks
pass. The actual SI verifier passes again. The only suite warnings are the
existing deprecated `imp` import and Sella external-cache permissions.
The earlier unchanged master suite gave 602 passes, 41 failures and 369
subtests; fixture incompatibilities are classified in `ANL_workflow_audit.md`.
That suite is not represented as passing or proof of universal compatibility.
